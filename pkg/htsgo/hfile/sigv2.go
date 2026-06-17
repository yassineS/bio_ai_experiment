package hfile

import (
	"crypto/hmac"
	"crypto/sha1"
	"encoding/base64"
	"sort"
	"strings"
	"time"
)

// sigV2DateFormat is the RFC 1123 GMT date used for both the Date request
// header and the Date line of the SigV2 string-to-sign. It matches the format
// htslib's hfile_s3.c emits ("%a, %d %b %Y %H:%M:%S GMT").
const sigV2DateFormat = "Mon, 02 Jan 2006 15:04:05 GMT"

// sigV2Input describes a single request to be signed with AWS Signature
// Version 2. Only the fixed header set that htslib signs is represented:
// the (empty) Content-MD5 and Content-Type, the Date, an optional
// x-amz-security-token, and the canonical resource path.
type sigV2Input struct {
	// Method is the HTTP verb, e.g. "GET".
	Method string
	// Resource is the CanonicalizedResource: "/bucket/key" with the bucket
	// always included (even under virtual-hosted addressing), already
	// percent-encoded, plus any signed subresources.
	Resource string
	// AccessKey is the AWS access key id.
	AccessKey string
	// SecretKey is the AWS secret access key.
	SecretKey string
	// SessionToken is the optional AWS session token. When non-empty it is
	// folded into the CanonicalizedAmzHeaders as x-amz-security-token and sent
	// as a request header.
	SessionToken string
	// Time is the request time; it determines the Date header and the Date
	// line of the string-to-sign, which must match.
	Time time.Time
}

// sigV2Result holds the headers that must be added to a request for it to be
// accepted once signed with Signature Version 2.
type sigV2Result struct {
	// Authorization is the full value of the Authorization header,
	// "AWS <AccessKeyId>:<base64-signature>".
	Authorization string
	// Date is the value of the Date header (RFC 1123 GMT). It MUST be sent
	// verbatim because it is part of the signed string-to-sign.
	Date string
	// SecurityToken is the value of the x-amz-security-token header, empty
	// when no session token was supplied.
	SecurityToken string
}

// hmacSHA1 returns the HMAC-SHA1 of data keyed by key.
func hmacSHA1(key, data []byte) []byte {
	h := hmac.New(sha1.New, key)
	h.Write(data)
	return h.Sum(nil)
}

// sigV2StringToSign builds the AWS Signature Version 2 string-to-sign for in
// given the already-formatted Date header value. The construction is:
//
//	HTTP-Verb + "\n" +
//	Content-MD5 + "\n" +     (empty for ranged GETs)
//	Content-Type + "\n" +    (empty for ranged GETs)
//	Date + "\n" +
//	CanonicalizedAmzHeaders +
//	CanonicalizedResource
//
// CanonicalizedAmzHeaders is the lower-cased, sorted set of x-amz-* headers,
// each rendered as "name:value\n"; for our requests this is only
// x-amz-security-token when a session token is present.
func sigV2StringToSign(in sigV2Input, date string) string {
	type kv struct{ k, v string }
	var amz []kv
	if in.SessionToken != "" {
		amz = append(amz, kv{"x-amz-security-token", in.SessionToken})
	}
	sort.Slice(amz, func(i, j int) bool { return amz[i].k < amz[j].k })

	var canonAmz strings.Builder
	for _, h := range amz {
		canonAmz.WriteString(h.k)
		canonAmz.WriteByte(':')
		canonAmz.WriteString(h.v)
		canonAmz.WriteByte('\n')
	}

	// Method\nContent-MD5\nContent-Type\nDate\nCanonicalizedAmzHeaders +
	// CanonicalizedResource. Content-MD5 and Content-Type are empty, producing
	// two blank lines. CanonicalizedAmzHeaders already carries its own trailing
	// newline (or is empty), and is immediately followed by the resource.
	return in.Method + "\n" +
		"\n" + // Content-MD5
		"\n" + // Content-Type
		date + "\n" +
		canonAmz.String() +
		in.Resource
}

// signSigV2 computes the AWS Signature Version 2 headers for in. The signature
// is base64(HMAC-SHA1(SecretKey, StringToSign)) and the Authorization header is
// "AWS <AccessKeyId>:<signature>". The returned Date must be sent as the Date
// request header because it is part of the signed string.
func signSigV2(in sigV2Input) sigV2Result {
	date := in.Time.UTC().Format(sigV2DateFormat)
	sts := sigV2StringToSign(in, date)
	sig := base64.StdEncoding.EncodeToString(hmacSHA1([]byte(in.SecretKey), []byte(sts)))
	return sigV2Result{
		Authorization: "AWS " + in.AccessKey + ":" + sig,
		Date:          date,
		SecurityToken: in.SessionToken,
	}
}
