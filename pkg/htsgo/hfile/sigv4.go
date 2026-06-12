package hfile

import (
	"crypto/hmac"
	"crypto/sha256"
	"encoding/hex"
	"sort"
	"strings"
	"time"
)

// unsignedPayload is the payload-hash placeholder htslib uses for ranged GET
// requests. The body is empty and S3 accepts the literal string in place of a
// real SHA-256 hash, so the request need not be buffered to be signed.
const unsignedPayload = "UNSIGNED-PAYLOAD"

// emptyStringSHA256 is the hex SHA-256 of the empty string, used when the
// caller supplies an empty payload hash.
const emptyStringSHA256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"

// sigV4Input describes a single request to be signed with AWS Signature
// Version 4. Headers other than host, x-amz-date, x-amz-content-sha256 and
// x-amz-security-token are not signed, matching htslib's hfile_s3.c which
// signs exactly that fixed set.
type sigV4Input struct {
	Method       string    // HTTP method, e.g. "GET"
	Host         string    // value of the Host header, e.g. "bucket.s3.us-east-1.amazonaws.com"
	CanonicalURI string    // URI-path component, already percent-encoded, e.g. "/key"
	Query        string    // canonical (sorted) query string, may be empty
	PayloadHash  string    // x-amz-content-sha256 value, e.g. "UNSIGNED-PAYLOAD"
	Region       string    // AWS region, e.g. "us-east-1"
	Service      string    // AWS service, e.g. "s3"
	AccessKey    string    // AWS access key id
	SecretKey    string    // AWS secret access key
	SessionToken string    // optional AWS session token
	Time         time.Time // request time; used for x-amz-date and the credential scope
}

// sigV4Result holds the headers that must be added to a request for it to be
// accepted by AWS once signed with Signature Version 4.
type sigV4Result struct {
	// Authorization is the full value of the Authorization header.
	Authorization string
	// AmzDate is the value of the x-amz-date header (ISO 8601 basic format).
	AmzDate string
	// ContentSHA256 is the value of the x-amz-content-sha256 header.
	ContentSHA256 string
	// SecurityToken is the value of the x-amz-security-token header. It is
	// empty when no session token was supplied.
	SecurityToken string
}

// hmacSHA256 returns the HMAC-SHA256 of data keyed by key.
func hmacSHA256(key, data []byte) []byte {
	h := hmac.New(sha256.New, key)
	h.Write(data)
	return h.Sum(nil)
}

// hexSHA256 returns the lower-case hex SHA-256 digest of s.
func hexSHA256(s string) string {
	sum := sha256.Sum256([]byte(s))
	return hex.EncodeToString(sum[:])
}

// canonicalRequest builds the AWS SigV4 canonical request string for in and
// returns it together with the semicolon-separated list of signed header
// names. The signed header set matches htslib: host, x-amz-content-sha256,
// x-amz-date and, when a session token is present, x-amz-security-token.
func canonicalRequest(in sigV4Input, amzDate string) (canonical, signedHeaders string) {
	payload := in.PayloadHash
	if payload == "" {
		payload = emptyStringSHA256
	}

	type kv struct{ k, v string }
	headers := []kv{
		{"host", in.Host},
		{"x-amz-content-sha256", payload},
		{"x-amz-date", amzDate},
	}
	if in.SessionToken != "" {
		headers = append(headers, kv{"x-amz-security-token", in.SessionToken})
	}
	sort.Slice(headers, func(i, j int) bool { return headers[i].k < headers[j].k })

	var ch strings.Builder
	names := make([]string, 0, len(headers))
	for _, h := range headers {
		ch.WriteString(h.k)
		ch.WriteByte(':')
		ch.WriteString(h.v)
		ch.WriteByte('\n')
		names = append(names, h.k)
	}
	signedHeaders = strings.Join(names, ";")

	uri := in.CanonicalURI
	if uri == "" {
		uri = "/"
	}

	canonical = strings.Join([]string{
		in.Method,
		uri,
		in.Query,
		ch.String(),
		signedHeaders,
		payload,
	}, "\n")
	return canonical, signedHeaders
}

// signSigV4 computes the AWS Signature Version 4 headers for in. It returns the
// Authorization header value along with the x-amz-date, x-amz-content-sha256
// and (when a session token is present) x-amz-security-token header values.
func signSigV4(in sigV4Input) sigV4Result {
	t := in.Time.UTC()
	amzDate := t.Format("20060102T150405Z")
	dateShort := t.Format("20060102")

	payload := in.PayloadHash
	if payload == "" {
		payload = emptyStringSHA256
	}

	canonical, signedHeaders := canonicalRequest(in, amzDate)

	scope := strings.Join([]string{dateShort, in.Region, in.Service, "aws4_request"}, "/")
	stringToSign := strings.Join([]string{
		"AWS4-HMAC-SHA256",
		amzDate,
		scope,
		hexSHA256(canonical),
	}, "\n")

	kDate := hmacSHA256([]byte("AWS4"+in.SecretKey), []byte(dateShort))
	kRegion := hmacSHA256(kDate, []byte(in.Region))
	kService := hmacSHA256(kRegion, []byte(in.Service))
	kSigning := hmacSHA256(kService, []byte("aws4_request"))
	signature := hex.EncodeToString(hmacSHA256(kSigning, []byte(stringToSign)))

	auth := "AWS4-HMAC-SHA256 " +
		"Credential=" + in.AccessKey + "/" + scope + ", " +
		"SignedHeaders=" + signedHeaders + ", " +
		"Signature=" + signature

	return sigV4Result{
		Authorization: auth,
		AmzDate:       amzDate,
		ContentSHA256: payload,
		SecurityToken: in.SessionToken,
	}
}
