package hfile

import (
	"errors"
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

// s3Handle is a read-only Handle for objects addressed by an s3:// URL. It
// delegates the actual HTTP transport to an httpHandle, installing a signer
// that adds AWS Signature Version 4 headers to every request.
type s3Handle struct {
	*httpHandle
}

// openS3 opens an s3://bucket/key URL for reading.
func openS3(rawurl string) (Handle, error) {
	if os.Getenv("HTS_S3_V2") != "" {
		return nil, errors.New("hfile: AWS SigV2 is not supported; unset HTS_S3_V2 to use SigV4")
	}

	bucket, key, err := parseS3URL(rawurl)
	if err != nil {
		return nil, err
	}

	creds := resolveAWSCredentials()
	region := creds.Region
	if region == "" {
		region = "us-east-1"
	}

	_, host, canonicalURI, fullURL := s3Endpoint(bucket, key, region)

	client, err := newHTTPClient()
	if err != nil {
		return nil, err
	}
	hh := &httpHandle{
		url:     fullURL,
		client:  client,
		retry:   loadRetryConfig(),
		headers: http.Header{},
	}

	hh.sign = func(req *http.Request) error {
		// Anonymous access: send the request unsigned (public buckets).
		if creds.AccessKey == "" || creds.SecretKey == "" {
			return nil
		}
		res := signSigV4(sigV4Input{
			Method:       req.Method,
			Host:         host,
			CanonicalURI: canonicalURI,
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

	return &s3Handle{httpHandle: hh}, nil
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

// s3Endpoint computes the request scheme, signing host, canonical URI and full
// request URL for a bucket/key pair. It mirrors hfile_s3.c: virtual-hosted
// style is used by default, but a bucket name containing dots (not DNS-safe
// for TLS virtual hosting) forces path-style, and HTS_S3_HOST overrides the
// host entirely (always path-style, as is conventional for custom endpoints).
func s3Endpoint(bucket, key, region string) (scheme, host, canonicalURI, fullURL string) {
	encodedKey := encodeS3Path(key)

	if s3HostOverride != "" {
		scheme = "https"
		host = s3HostOverride
		if strings.HasPrefix(host, "http://") {
			scheme = "http"
			host = strings.TrimPrefix(host, "http://")
		} else {
			host = strings.TrimPrefix(host, "https://")
		}
		host = strings.TrimSuffix(host, "/")
		// Custom endpoints use path-style: host/bucket/key.
		canonicalURI = "/" + bucket + "/" + encodedKey
		fullURL = scheme + "://" + host + canonicalURI
		return scheme, host, canonicalURI, fullURL
	}

	scheme = "https"
	if strings.Contains(bucket, ".") {
		// Path-style: dotted bucket names break TLS virtual-hosted certs.
		host = "s3." + region + ".amazonaws.com"
		canonicalURI = "/" + bucket + "/" + encodedKey
	} else {
		// Virtual-hosted style.
		host = bucket + ".s3." + region + ".amazonaws.com"
		canonicalURI = "/" + encodedKey
	}
	fullURL = scheme + "://" + host + canonicalURI
	return scheme, host, canonicalURI, fullURL
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
