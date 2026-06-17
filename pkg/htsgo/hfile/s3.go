package hfile

import (
	"fmt"
	"net/http"
	"os"
	"strings"
	"time"
)

// s3HostOverride, when non-empty, replaces the host used to reach S3. It is
// populated from the HTS_S3_HOST environment variable and may also be set
// directly by tests to point the backend at an httptest server. The value is
// a bare host (optionally host:port); the scheme defaults to https unless the
// override begins with "http://".
var s3HostOverride = os.Getenv("HTS_S3_HOST")

// nowFunc returns the current time. It is a package variable so that tests can
// pin the clock when asserting on signatures.
var nowFunc = time.Now

// s3AddressStyle is one of "auto", "virtual" or "path", mirroring htslib's
// HTS_S3_ADDRESS_STYLE env knob (and the addressing_style / host_bucket
// config-file settings). "auto" picks virtual-hosted style for DNS-safe
// buckets and path-style otherwise; "virtual" and "path" force the choice.
type s3AddressStyle int

const (
	s3StyleAuto s3AddressStyle = iota
	s3StyleVirtual
	s3StylePath
)

// resolveS3AddressStyle reads HTS_S3_ADDRESS_STYLE. Unset or unrecognised
// values mean "auto". The comparison is case-insensitive, matching htslib's
// strcasecmp.
func resolveS3AddressStyle() s3AddressStyle {
	switch strings.ToLower(os.Getenv("HTS_S3_ADDRESS_STYLE")) {
	case "virtual":
		return s3StyleVirtual
	case "path":
		return s3StylePath
	default:
		return s3StyleAuto
	}
}

// s3Handle is a read-only Handle for objects addressed by an s3:// URL. It
// delegates the actual HTTP transport to an httpHandle, installing a signer
// that adds AWS Signature Version 4 (default) or Version 2 (HTS_S3_V2) headers
// to every request.
type s3Handle struct {
	*httpHandle
}

// openS3 opens an s3://bucket/key URL for reading. By default requests are
// signed with AWS Signature Version 4; setting HTS_S3_V2 selects Signature
// Version 2 instead, as required by some older or self-hosted S3-compatible
// endpoints (and matching htslib's hfile_s3.c).
func openS3(rawurl string) (Handle, error) {
	bucket, key, err := parseS3URL(rawurl)
	if err != nil {
		return nil, err
	}

	creds := resolveAWSCredentials()
	region := creds.Region
	if region == "" {
		region = "us-east-1"
	}

	ep := s3Endpoint(bucket, key, region)

	client, err := newHTTPClient()
	if err != nil {
		return nil, err
	}
	hh := &httpHandle{
		url:     ep.fullURL,
		client:  client,
		retry:   loadRetryConfig(),
		headers: http.Header{},
	}

	if os.Getenv("HTS_S3_V2") != "" {
		hh.sign = s3SignerV2(ep, creds)
	} else {
		hh.sign = s3SignerV4(ep, creds, region)
	}

	return &s3Handle{httpHandle: hh}, nil
}

// s3SignerV4 returns a request signer that adds AWS Signature Version 4
// headers. Anonymous access (missing credentials) is sent unsigned so that
// public buckets remain readable.
func s3SignerV4(ep s3EndpointInfo, creds awsCredentials, region string) requestSigner {
	return func(req *http.Request) error {
		if creds.AccessKey == "" || creds.SecretKey == "" {
			return nil
		}
		res := signSigV4(sigV4Input{
			Method:       req.Method,
			Host:         ep.host,
			CanonicalURI: ep.canonicalURI,
			Query:        "",
			PayloadHash:  unsignedPayload,
			Region:       region,
			Service:      "s3",
			AccessKey:    creds.AccessKey,
			SecretKey:    creds.SecretKey,
			SessionToken: creds.SessionToken,
			Time:         nowFunc(),
		})
		req.Header.Set("x-amz-date", res.AmzDate)
		req.Header.Set("x-amz-content-sha256", res.ContentSHA256)
		if res.SecurityToken != "" {
			req.Header.Set("x-amz-security-token", res.SecurityToken)
		}
		req.Header.Set("Authorization", res.Authorization)
		return nil
	}
}

// s3SignerV2 returns a request signer that adds AWS Signature Version 2
// headers. The Date header it sets is part of the signed string-to-sign and so
// must be transmitted verbatim. Anonymous access is sent unsigned (but still
// carries a Date header, matching htslib).
func s3SignerV2(ep s3EndpointInfo, creds awsCredentials) requestSigner {
	return func(req *http.Request) error {
		if creds.AccessKey == "" || creds.SecretKey == "" {
			// Anonymous: still send a Date header as htslib does, but no
			// Authorization.
			req.Header.Set("Date", nowFunc().UTC().Format(sigV2DateFormat))
			return nil
		}
		res := signSigV2(sigV2Input{
			Method:       req.Method,
			Resource:     ep.canonicalResource,
			AccessKey:    creds.AccessKey,
			SecretKey:    creds.SecretKey,
			SessionToken: creds.SessionToken,
			Time:         nowFunc(),
		})
		req.Header.Set("Date", res.Date)
		if res.SecurityToken != "" {
			req.Header.Set("x-amz-security-token", res.SecurityToken)
		}
		req.Header.Set("Authorization", res.Authorization)
		return nil
	}
}

// parseS3URL splits an s3://bucket/key URL into its bucket and key components.
func parseS3URL(rawurl string) (bucket, key string, err error) {
	rest := strings.TrimPrefix(rawurl, "s3://")
	if rest == rawurl {
		return "", "", fmt.Errorf("hfile: not an s3:// URL: %q", rawurl)
	}
	slash := strings.IndexByte(rest, '/')
	if slash < 0 {
		return "", "", fmt.Errorf("hfile: s3 URL %q has no key", rawurl)
	}
	bucket = rest[:slash]
	key = rest[slash+1:]
	if bucket == "" || key == "" {
		return "", "", fmt.Errorf("hfile: s3 URL %q must be s3://bucket/key", rawurl)
	}
	return bucket, key, nil
}

// s3EndpointInfo describes how a bucket/key pair is addressed: the wire URL
// and signing host plus both the SigV4 canonical URI (the request path) and
// the SigV2 canonical resource. For virtual-hosted addressing these differ:
// the request path omits the bucket while the SigV2 resource always includes
// it ("/bucket/key").
type s3EndpointInfo struct {
	// scheme is "http" or "https".
	scheme string
	// host is the value of the Host header used for signing, e.g.
	// "bucket.s3.us-east-1.amazonaws.com" or a custom endpoint.
	host string
	// canonicalURI is the request path, already percent-encoded. Under
	// virtual-hosted addressing this is "/key"; under path-style it is
	// "/bucket/key".
	canonicalURI string
	// canonicalResource is the SigV2 CanonicalizedResource, always
	// "/bucket/key" regardless of addressing style.
	canonicalResource string
	// fullURL is the complete request URL.
	fullURL string
}

// s3Endpoint computes addressing details for a bucket/key pair. It mirrors
// hfile_s3.c: virtual-hosted style is used by default, but a bucket name
// containing dots (not DNS-safe for TLS virtual hosting) forces path-style,
// HTS_S3_ADDRESS_STYLE may force "virtual" or "path" explicitly, and
// HTS_S3_HOST overrides the host entirely (always path-style, as is
// conventional for custom S3-compatible endpoints such as MinIO, Ceph,
// Wasabi, Backblaze and Google Cloud Storage's S3 interop endpoint).
func s3Endpoint(bucket, key, region string) s3EndpointInfo {
	encodedKey := encodeS3Path(key)
	resource := "/" + bucket + "/" + encodedKey

	if s3HostOverride != "" {
		scheme := "https"
		host := s3HostOverride
		if strings.HasPrefix(host, "http://") {
			scheme = "http"
			host = strings.TrimPrefix(host, "http://")
		} else {
			host = strings.TrimPrefix(host, "https://")
		}
		host = strings.TrimSuffix(host, "/")
		// Custom endpoints use path-style: host/bucket/key.
		uri := "/" + bucket + "/" + encodedKey
		return s3EndpointInfo{
			scheme:            scheme,
			host:              host,
			canonicalURI:      uri,
			canonicalResource: resource,
			fullURL:           scheme + "://" + host + uri,
		}
	}

	style := resolveS3AddressStyle()
	virtual := style == s3StyleVirtual
	if style == s3StyleAuto {
		// Dotted bucket names break TLS virtual-hosted certs, so use
		// path-style for them; everything else is virtual-hosted.
		virtual = !strings.Contains(bucket, ".")
	}

	scheme := "https"
	var host, uri string
	if virtual {
		host = bucket + ".s3." + region + ".amazonaws.com"
		uri = "/" + encodedKey
	} else {
		host = "s3." + region + ".amazonaws.com"
		uri = "/" + bucket + "/" + encodedKey
	}
	return s3EndpointInfo{
		scheme:            scheme,
		host:              host,
		canonicalURI:      uri,
		canonicalResource: resource,
		fullURL:           scheme + "://" + host + uri,
	}
}

// encodeS3Path percent-encodes an S3 object key for use in a URI path,
// preserving '/' separators. The set of unreserved characters matches
// hfile_s3.c's escape_path (A-Z a-z 0-9 _ - ~ . and /).
func encodeS3Path(key string) string {
	var b strings.Builder
	for i := 0; i < len(key); i++ {
		c := key[i]
		if (c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') ||
			(c >= 'a' && c <= 'z') ||
			c == '_' || c == '-' || c == '~' || c == '.' || c == '/' {
			b.WriteByte(c)
		} else {
			fmt.Fprintf(&b, "%%%02X", c)
		}
	}
	return b.String()
}
