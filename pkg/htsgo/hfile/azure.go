package hfile

import (
	"crypto/hmac"
	"crypto/sha256"
	"encoding/base64"
	"fmt"
	"net/http"
	"os"
	"sort"
	"strings"
)

// Azure Blob Storage support is an additive Go-port extension: htslib has no
// native Azure backend. Reads are served by ranged GETs over the standard Azure
// Blob REST endpoint (https://<account>.blob.core.windows.net/<container>/<blob>),
// with four authentication paths handled here:
//
//  1. SAS URL  — a fully-inline https://...?<sas> URL is just a signed HTTPS GET
//     and flows through the plain HTTP backend with no Azure-specific code. An
//     AZURE_STORAGE_SAS_TOKEN env var is appended to az:// / host URLs that lack
//     an inline SAS.
//  2. Shared Key — AZURE_STORAGE_ACCOUNT + AZURE_STORAGE_KEY, signed per-request
//     (the canonical string includes the request's Range header).
//  3. Azure AD bearer — AZURE_STORAGE_TOKEN -> Authorization: Bearer <token>.
//  4. Anonymous — unsigned, but still carrying x-ms-version.

// azureBlobEndpoint is the default Azure Blob host suffix. The wire endpoint is
// https://<account>.blob.core.windows.net.
const azureBlobHostSuffix = ".blob.core.windows.net"

// azureAPIVersion is the x-ms-version sent on every request and folded into the
// Shared Key canonicalized headers.
const azureAPIVersion = "2021-08-06"

// azureDateFormat is the RFC 1123 GMT format used for the x-ms-date header and
// the Date value embedded in the Shared Key string-to-sign.
const azureDateFormat = "Mon, 02 Jan 2006 15:04:05 GMT"

// azureEndpointOverride, when non-empty, replaces the scheme+host used to reach
// Azure Blob storage. It is populated from AZURE_STORAGE_BLOB_ENDPOINT and may
// also be set directly by tests to point the backend at an httptest server. The
// value may include a scheme ("http://host:port"); a bare host defaults to
// https.
var azureEndpointOverride = os.Getenv("AZURE_STORAGE_BLOB_ENDPOINT")

// azureHandle is a read-only Handle for objects addressed by an az:// URL (or a
// recognised *.blob.core.windows.net https URL). It delegates HTTP transport to
// an httpHandle, optionally installing a Shared Key signer or static auth
// headers.
type azureHandle struct {
	*httpHandle
}

// isAzureHTTPSURL reports whether rawurl is a direct HTTPS URL to the Azure Blob
// service (host ending in .blob.core.windows.net). Such URLs are routed through
// the Azure backend so the env-driven SAS / Shared Key / bearer auth applies.
// A URL that already carries an inline SAS query string needs no Azure auth and
// would also work through the plain HTTP backend; routing it here is harmless
// because the Azure backend leaves an existing SAS untouched.
func isAzureHTTPSURL(rawurl string) bool {
	const httpsPrefix = "https://"
	if !strings.HasPrefix(rawurl, httpsPrefix) {
		return false
	}
	rest := rawurl[len(httpsPrefix):]
	// Host is the run up to the first '/', '?' or end.
	host := rest
	if i := strings.IndexAny(host, "/?"); i >= 0 {
		host = host[:i]
	}
	host = strings.ToLower(host)
	return strings.HasSuffix(host, azureBlobHostSuffix)
}

// azureParts describes the account/container/blob decomposition of an Azure
// request, plus any inline query string carried by the source URL.
type azureParts struct {
	account   string
	container string
	blob      string // un-encoded blob name (may contain '/')
	query     string // inline query string from the source URL, without leading '?'
}

// parseAzureURL decomposes either an az://<account>/<container>/<blob> URL or a
// direct https://<account>.blob.core.windows.net/<container>/<blob> URL into its
// account, container and blob components, preserving any inline query string.
func parseAzureURL(rawurl string) (azureParts, error) {
	var p azureParts

	var rest string
	switch {
	case strings.HasPrefix(rawurl, "az://"):
		rest = rawurl[len("az://"):]
		// Split off an inline query string first.
		if q := strings.IndexByte(rest, '?'); q >= 0 {
			p.query = rest[q+1:]
			rest = rest[:q]
		}
		// rest = account/container/blob...
		slash := strings.IndexByte(rest, '/')
		if slash < 0 {
			return p, fmt.Errorf("hfile: az URL %q must be az://account/container/blob", rawurl)
		}
		p.account = rest[:slash]
		rest = rest[slash+1:]
	case isAzureHTTPSURL(rawurl):
		rest = rawurl[len("https://"):]
		// host is up to first '/'.
		slash := strings.IndexByte(rest, '/')
		if slash < 0 {
			return p, fmt.Errorf("hfile: Azure URL %q has no container/blob path", rawurl)
		}
		host := rest[:slash]
		rest = rest[slash+1:]
		// Account is the label before .blob.core.windows.net.
		p.account = strings.ToLower(host)
		if i := strings.Index(p.account, "."); i >= 0 {
			p.account = p.account[:i]
		}
		if q := strings.IndexByte(rest, '?'); q >= 0 {
			p.query = rest[q+1:]
			rest = rest[:q]
		}
	default:
		return p, fmt.Errorf("hfile: not an Azure URL: %q", rawurl)
	}

	// rest is now container/blob...
	slash := strings.IndexByte(rest, '/')
	if slash < 0 {
		return p, fmt.Errorf("hfile: Azure URL %q must include container/blob", rawurl)
	}
	p.container = rest[:slash]
	p.blob = rest[slash+1:]
	if p.account == "" || p.container == "" || p.blob == "" {
		return p, fmt.Errorf("hfile: Azure URL %q must be account/container/blob", rawurl)
	}
	return p, nil
}

// azureBaseURL returns the scheme+host for an account, honouring
// AZURE_STORAGE_BLOB_ENDPOINT (the test/endpoint override). The returned value
// has no trailing slash.
func azureBaseURL(account string) string {
	if azureEndpointOverride != "" {
		ov := azureEndpointOverride
		if !strings.HasPrefix(ov, "http://") && !strings.HasPrefix(ov, "https://") {
			ov = "https://" + ov
		}
		return strings.TrimSuffix(ov, "/")
	}
	return "https://" + account + azureBlobHostSuffix
}

// openAzure opens an az:// URL or a recognised *.blob.core.windows.net https URL
// for reading. The authentication path is chosen by environment, in priority
// order: an inline or AZURE_STORAGE_SAS_TOKEN SAS, then Shared Key
// (AZURE_STORAGE_ACCOUNT + AZURE_STORAGE_KEY), then an Azure AD bearer token
// (AZURE_STORAGE_TOKEN), then anonymous. In every case x-ms-version is sent.
func openAzure(rawurl string) (Handle, error) {
	p, err := parseAzureURL(rawurl)
	if err != nil {
		return nil, err
	}

	base := azureBaseURL(p.account)
	fullURL := base + "/" + p.container + "/" + encodeS3Path(p.blob)

	// An inline SAS, or AZURE_STORAGE_SAS_TOKEN, is appended to the query string.
	query := p.query
	if query == "" {
		if sas := os.Getenv("AZURE_STORAGE_SAS_TOKEN"); sas != "" {
			query = strings.TrimPrefix(sas, "?")
		}
	}
	if query != "" {
		fullURL += "?" + query
	}

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
	hh.headers.Set("x-ms-version", azureAPIVersion)

	switch {
	case query != "":
		// SAS path: the signature is in the URL; no extra auth header is needed.
		// x-ms-version is still sent (harmless and matches the SDK convention).
	case os.Getenv("AZURE_STORAGE_ACCOUNT") != "" && os.Getenv("AZURE_STORAGE_KEY") != "":
		acct := os.Getenv("AZURE_STORAGE_ACCOUNT")
		key := os.Getenv("AZURE_STORAGE_KEY")
		// CanonicalizedResource is built from the SAS-less account/container/blob.
		canonResource := "/" + acct + "/" + p.container + "/" + p.blob
		hh.sign = azureSharedKeySigner(acct, key, canonResource)
	case os.Getenv("AZURE_STORAGE_TOKEN") != "":
		hh.headers.Set("Authorization", "Bearer "+os.Getenv("AZURE_STORAGE_TOKEN"))
	default:
		// Anonymous: unsigned, but x-ms-version is already set.
	}

	return &azureHandle{httpHandle: hh}, nil
}

// azureSharedKeyInput describes a single request to be signed with the Azure
// Blob Shared Key scheme.
type azureSharedKeyInput struct {
	// Method is the HTTP verb, e.g. "GET".
	Method string
	// Range is the exact value of the request's Range header (e.g.
	// "bytes=0-65535"), or "" for an un-ranged request.
	Range string
	// XMSDate is the RFC 1123 GMT date placed in the x-ms-date header. The Date
	// line of the string-to-sign is left blank because x-ms-date is used.
	XMSDate string
	// XMSVersion is the x-ms-version header value.
	XMSVersion string
	// CanonicalizedResource is "/<account>/<container>/<blob>" (no query params
	// for a plain blob GET).
	CanonicalizedResource string
}

// azureSharedKeyStringToSign builds the Azure Blob Shared Key string-to-sign:
//
//	VERB\n
//	Content-Encoding\n Content-Language\n Content-Length\n Content-MD5\n
//	Content-Type\n Date\n If-Modified-Since\n If-Match\n If-None-Match\n
//	If-Unmodified-Since\n Range\n
//	CanonicalizedHeaders
//	CanonicalizedResource
//
// For a ranged GET, Content-Length is empty (not "0"), the Date line is empty
// because x-ms-date carries the timestamp, Range is the request's exact Range
// header value, and CanonicalizedHeaders is the sorted x-ms-date and
// x-ms-version lines.
func azureSharedKeyStringToSign(in azureSharedKeyInput) string {
	// CanonicalizedHeaders: all x-ms-* headers, lowercased name, sorted by name,
	// "name:value\n". We set x-ms-date and x-ms-version.
	type kv struct{ k, v string }
	hdrs := []kv{
		{"x-ms-date", in.XMSDate},
		{"x-ms-version", in.XMSVersion},
	}
	sort.Slice(hdrs, func(i, j int) bool { return hdrs[i].k < hdrs[j].k })
	var canonHeaders strings.Builder
	for _, h := range hdrs {
		canonHeaders.WriteString(h.k)
		canonHeaders.WriteByte(':')
		canonHeaders.WriteString(h.v)
		canonHeaders.WriteByte('\n')
	}

	return in.Method + "\n" +
		"\n" + // Content-Encoding
		"\n" + // Content-Language
		"\n" + // Content-Length (empty, not 0, for GET)
		"\n" + // Content-MD5
		"\n" + // Content-Type
		"\n" + // Date (empty; x-ms-date is used)
		"\n" + // If-Modified-Since
		"\n" + // If-Match
		"\n" + // If-None-Match
		"\n" + // If-Unmodified-Since
		in.Range + "\n" + // Range
		canonHeaders.String() +
		in.CanonicalizedResource
}

// azureSharedKeySign computes the Shared Key Authorization header value for in,
// using accountKey (the base64-encoded account key). The signature is
// base64(HMAC-SHA256(base64decode(accountKey), StringToSign)) and the header is
// "SharedKey <account>:<signature>". It returns an error if accountKey is not
// valid base64.
func azureSharedKeySign(account, accountKey string, in azureSharedKeyInput) (string, error) {
	rawKey, err := base64.StdEncoding.DecodeString(accountKey)
	if err != nil {
		return "", fmt.Errorf("hfile: AZURE_STORAGE_KEY is not valid base64: %w", err)
	}
	sts := azureSharedKeyStringToSign(in)
	mac := hmac.New(sha256.New, rawKey)
	mac.Write([]byte(sts))
	sig := base64.StdEncoding.EncodeToString(mac.Sum(nil))
	return "SharedKey " + account + ":" + sig, nil
}

// azureSharedKeySigner returns a requestSigner that signs each request with the
// Azure Blob Shared Key scheme. The signature is recomputed per request because
// the canonical string-to-sign includes the request's exact Range header, which
// differs for each ranged GET.
func azureSharedKeySigner(account, accountKey, canonResource string) requestSigner {
	return func(req *http.Request) error {
		date := nowFunc().UTC().Format(azureDateFormat)
		req.Header.Set("x-ms-date", date)
		// x-ms-version is already on the static headers, but ensure it is present
		// for the signature even if a caller cleared it.
		ver := req.Header.Get("x-ms-version")
		if ver == "" {
			ver = azureAPIVersion
			req.Header.Set("x-ms-version", ver)
		}
		auth, err := azureSharedKeySign(account, accountKey, azureSharedKeyInput{
			Method:                req.Method,
			Range:                 req.Header.Get("Range"),
			XMSDate:               date,
			XMSVersion:            ver,
			CanonicalizedResource: canonResource,
		})
		if err != nil {
			return err
		}
		req.Header.Set("Authorization", auth)
		return nil
	}
}
