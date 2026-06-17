package hfile

import (
	"bytes"
	"encoding/base64"
	"errors"
	"fmt"
	"io"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
	"time"
)

// payload is a deterministic test body used across the remote-backend tests.
var payload = func() []byte {
	b := make([]byte, 4096)
	for i := range b {
		b[i] = byte(i * 7)
	}
	return b
}()

// rangeServer is an httptest handler that serves payload honouring Range
// requests, recording the most recent request for assertions.
type rangeServer struct {
	body    []byte
	lastReq *http.Request
}

func (s *rangeServer) ServeHTTP(w http.ResponseWriter, r *http.Request) {
	s.lastReq = r.Clone(r.Context())

	if r.Method == http.MethodHead {
		w.Header().Set("Content-Length", strconv.Itoa(len(s.body)))
		w.WriteHeader(http.StatusOK)
		return
	}

	rng := r.Header.Get("Range")
	if rng == "" {
		w.Header().Set("Content-Length", strconv.Itoa(len(s.body)))
		w.WriteHeader(http.StatusOK)
		w.Write(s.body)
		return
	}

	var start, end int64
	if _, err := fmt.Sscanf(rng, "bytes=%d-%d", &start, &end); err != nil {
		w.WriteHeader(http.StatusBadRequest)
		return
	}
	if start >= int64(len(s.body)) {
		w.Header().Set("Content-Range", fmt.Sprintf("bytes */%d", len(s.body)))
		w.WriteHeader(http.StatusRequestedRangeNotSatisfiable)
		return
	}
	if end >= int64(len(s.body)) {
		end = int64(len(s.body)) - 1
	}
	w.Header().Set("Content-Range", fmt.Sprintf("bytes %d-%d/%d", start, end, len(s.body)))
	w.Header().Set("Content-Length", strconv.FormatInt(end-start+1, 10))
	w.WriteHeader(http.StatusPartialContent)
	w.Write(s.body[start : end+1])
}

func newRangeServer(body []byte) (*rangeServer, *httptest.Server) {
	rs := &rangeServer{body: body}
	return rs, httptest.NewServer(rs)
}

// --- ReadFile -------------------------------------------------------------

func TestReadFileRemoteAndLocal(t *testing.T) {
	// Remote: download an HTTP object in full.
	_, srv := newRangeServer(payload)
	defer srv.Close()
	got, err := ReadFile(srv.URL)
	if err != nil {
		t.Fatalf("ReadFile(remote): %v", err)
	}
	if !bytes.Equal(got, payload) {
		t.Fatalf("ReadFile(remote): got %d bytes, want %d", len(got), len(payload))
	}

	// Local: read a temp file, including the file:// form.
	dir := t.TempDir()
	path := filepath.Join(dir, "local.dat")
	want := []byte("local file contents\n")
	if err := os.WriteFile(path, want, 0o600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	for _, name := range []string{path, "file://" + path} {
		got, err := ReadFile(name)
		if err != nil {
			t.Fatalf("ReadFile(%q): %v", name, err)
		}
		if !bytes.Equal(got, want) {
			t.Fatalf("ReadFile(%q): got %q, want %q", name, got, want)
		}
	}
}

// --- HTTP backend tests -----------------------------------------------------

func TestHTTPReadAt(t *testing.T) {
	_, srv := newRangeServer(payload)
	defer srv.Close()

	h, err := Open(srv.URL)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer h.Close()

	cases := []struct{ off, n int }{
		{0, 16}, {100, 256}, {len(payload) - 10, 10}, {1000, 1},
	}
	for _, c := range cases {
		buf := make([]byte, c.n)
		got, err := h.ReadAt(buf, int64(c.off))
		if err != nil {
			t.Fatalf("ReadAt(%d,%d): %v", c.off, c.n, err)
		}
		if got != c.n {
			t.Fatalf("ReadAt(%d,%d): n=%d", c.off, c.n, got)
		}
		if !bytes.Equal(buf, payload[c.off:c.off+c.n]) {
			t.Fatalf("ReadAt(%d,%d): wrong bytes", c.off, c.n)
		}
	}
}

func TestHTTPReadAtEOF(t *testing.T) {
	_, srv := newRangeServer(payload)
	defer srv.Close()
	h, _ := Open(srv.URL)
	defer h.Close()

	// Straddle the end of the resource.
	buf := make([]byte, 100)
	n, err := h.ReadAt(buf, int64(len(payload)-10))
	if n != 10 {
		t.Fatalf("n=%d, want 10", n)
	}
	if err != io.EOF {
		t.Fatalf("err=%v, want io.EOF", err)
	}
	if !bytes.Equal(buf[:10], payload[len(payload)-10:]) {
		t.Fatal("wrong tail bytes")
	}

	// Read entirely past EOF.
	n, err = h.ReadAt(buf, int64(len(payload)))
	if n != 0 || err != io.EOF {
		t.Fatalf("past-EOF: n=%d err=%v", n, err)
	}
}

func TestHTTPSize(t *testing.T) {
	_, srv := newRangeServer(payload)
	defer srv.Close()
	h, _ := Open(srv.URL)
	defer h.Close()

	sz, err := h.Size()
	if err != nil {
		t.Fatalf("Size: %v", err)
	}
	if sz != int64(len(payload)) {
		t.Fatalf("Size=%d, want %d", sz, len(payload))
	}
}

func TestHTTPSizeViaContentRange(t *testing.T) {
	// Server that refuses HEAD, forcing the Content-Range fallback.
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.Method == http.MethodHead {
			w.WriteHeader(http.StatusMethodNotAllowed)
			return
		}
		w.Header().Set("Content-Range", fmt.Sprintf("bytes 0-0/%d", len(payload)))
		w.WriteHeader(http.StatusPartialContent)
		w.Write(payload[:1])
	}))
	defer srv.Close()

	h, _ := Open(srv.URL)
	defer h.Close()
	sz, err := h.Size()
	if err != nil {
		t.Fatalf("Size: %v", err)
	}
	if sz != int64(len(payload)) {
		t.Fatalf("Size=%d, want %d", sz, len(payload))
	}
}

func TestHTTPSequentialRead(t *testing.T) {
	_, srv := newRangeServer(payload)
	defer srv.Close()
	h, _ := Open(srv.URL)
	defer h.Close()

	var out bytes.Buffer
	buf := make([]byte, 333)
	for {
		n, err := h.Read(buf)
		out.Write(buf[:n])
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read: %v", err)
		}
	}
	if !bytes.Equal(out.Bytes(), payload) {
		t.Fatalf("sequential read mismatch: got %d bytes", out.Len())
	}
}

func TestHTTPRetry(t *testing.T) {
	t.Setenv("HTS_RETRY_MAX", "3")
	t.Setenv("HTS_RETRY_DELAY", "1")
	var attempts int
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		attempts++
		if attempts < 3 {
			w.WriteHeader(http.StatusServiceUnavailable)
			return
		}
		w.Header().Set("Content-Range", fmt.Sprintf("bytes 0-3/%d", len(payload)))
		w.WriteHeader(http.StatusPartialContent)
		w.Write(payload[:4])
	}))
	defer srv.Close()

	h, _ := Open(srv.URL)
	defer h.Close()
	buf := make([]byte, 4)
	if _, err := h.ReadAt(buf, 0); err != nil {
		t.Fatalf("ReadAt after retries: %v", err)
	}
	if attempts != 3 {
		t.Fatalf("attempts=%d, want 3", attempts)
	}
}

func TestHTTPAuthLocation(t *testing.T) {
	dir := t.TempDir()
	authFile := filepath.Join(dir, "auth")
	if err := os.WriteFile(authFile, []byte("Bearer secret-token\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	t.Setenv("HTS_AUTH_LOCATION", authFile)

	rs, srv := newRangeServer(payload)
	defer srv.Close()
	h, _ := Open(srv.URL)
	defer h.Close()

	buf := make([]byte, 8)
	if _, err := h.ReadAt(buf, 0); err != nil {
		t.Fatalf("ReadAt: %v", err)
	}
	if got := rs.lastReq.Header.Get("Authorization"); got != "Bearer secret-token" {
		t.Fatalf("Authorization=%q", got)
	}
}

// --- S3 backend tests -------------------------------------------------------

func TestS3ReadAndSign(t *testing.T) {
	rs, srv := newRangeServer(payload)
	defer srv.Close()

	host := strings.TrimPrefix(srv.URL, "http://")
	t.Setenv("HTS_S3_HOST", "http://"+host)
	t.Setenv("AWS_ACCESS_KEY_ID", "AKIDEXAMPLE")
	t.Setenv("AWS_SECRET_ACCESS_KEY", "wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY")
	t.Setenv("AWS_DEFAULT_REGION", "us-east-1")
	t.Setenv("AWS_SESSION_TOKEN", "")
	s3HostOverride = "http://" + host
	defer func() { s3HostOverride = "" }()

	h, err := Open("s3://my-bucket/path/to/object.bam")
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer h.Close()

	buf := make([]byte, 64)
	n, err := h.ReadAt(buf, 128)
	if err != nil {
		t.Fatalf("ReadAt: %v", err)
	}
	if !bytes.Equal(buf[:n], payload[128:128+64]) {
		t.Fatal("S3 bytes mismatch")
	}

	r := rs.lastReq
	auth := r.Header.Get("Authorization")
	if !strings.HasPrefix(auth, "AWS4-HMAC-SHA256 ") {
		t.Fatalf("Authorization=%q", auth)
	}
	for _, want := range []string{"Credential=AKIDEXAMPLE/", "/us-east-1/s3/aws4_request", "SignedHeaders=", "Signature="} {
		if !strings.Contains(auth, want) {
			t.Fatalf("Authorization %q missing %q", auth, want)
		}
	}
	if r.Header.Get("x-amz-date") == "" {
		t.Fatal("missing x-amz-date")
	}
	if r.Header.Get("x-amz-content-sha256") != unsignedPayload {
		t.Fatalf("x-amz-content-sha256=%q", r.Header.Get("x-amz-content-sha256"))
	}
	// Path-style URL contains the bucket.
	if !strings.Contains(r.URL.Path, "my-bucket") {
		t.Fatalf("path=%q", r.URL.Path)
	}
}

func TestS3SessionToken(t *testing.T) {
	rs, srv := newRangeServer(payload)
	defer srv.Close()
	host := strings.TrimPrefix(srv.URL, "http://")
	t.Setenv("AWS_ACCESS_KEY_ID", "AKIDEXAMPLE")
	t.Setenv("AWS_SECRET_ACCESS_KEY", "secret")
	t.Setenv("AWS_SESSION_TOKEN", "FQoGZXIvYXdzE-token")
	s3HostOverride = "http://" + host
	defer func() { s3HostOverride = "" }()

	h, _ := Open("s3://bucket/key")
	defer h.Close()
	buf := make([]byte, 4)
	if _, err := h.ReadAt(buf, 0); err != nil {
		t.Fatalf("ReadAt: %v", err)
	}
	if rs.lastReq.Header.Get("x-amz-security-token") != "FQoGZXIvYXdzE-token" {
		t.Fatalf("security token=%q", rs.lastReq.Header.Get("x-amz-security-token"))
	}
	if !strings.Contains(rs.lastReq.Header.Get("Authorization"), "x-amz-security-token") {
		t.Fatal("signed headers missing security token")
	}
}

func TestS3Endpoint(t *testing.T) {
	s3HostOverride = ""
	t.Setenv("HTS_S3_ADDRESS_STYLE", "")
	// Virtual-hosted for a simple bucket.
	ep := s3Endpoint("genomics", "a/b.bam", "eu-west-1")
	if ep.host != "genomics.s3.eu-west-1.amazonaws.com" {
		t.Fatalf("host=%q", ep.host)
	}
	if ep.canonicalURI != "/a/b.bam" {
		t.Fatalf("uri=%q", ep.canonicalURI)
	}
	// The SigV2 canonical resource always includes the bucket.
	if ep.canonicalResource != "/genomics/a/b.bam" {
		t.Fatalf("resource=%q", ep.canonicalResource)
	}
	if ep.fullURL != "https://genomics.s3.eu-west-1.amazonaws.com/a/b.bam" {
		t.Fatalf("full=%q", ep.fullURL)
	}
	// Path-style for a dotted bucket.
	ep = s3Endpoint("my.dotted.bucket", "key", "us-east-1")
	if ep.host != "s3.us-east-1.amazonaws.com" {
		t.Fatalf("dotted host=%q", ep.host)
	}
	if ep.canonicalURI != "/my.dotted.bucket/key" {
		t.Fatalf("dotted uri=%q", ep.canonicalURI)
	}
	if ep.canonicalResource != "/my.dotted.bucket/key" {
		t.Fatalf("dotted resource=%q", ep.canonicalResource)
	}
}

// TestS3AddressStyleForced verifies HTS_S3_ADDRESS_STYLE overrides the
// auto-detected addressing: "path" forces path-style for a DNS-safe bucket and
// "virtual" forces virtual-hosted style.
func TestS3AddressStyleForced(t *testing.T) {
	s3HostOverride = ""

	t.Setenv("HTS_S3_ADDRESS_STYLE", "path")
	ep := s3Endpoint("simple", "k", "us-east-1")
	if ep.host != "s3.us-east-1.amazonaws.com" || ep.canonicalURI != "/simple/k" {
		t.Fatalf("forced path: host=%q uri=%q", ep.host, ep.canonicalURI)
	}

	t.Setenv("HTS_S3_ADDRESS_STYLE", "virtual")
	ep = s3Endpoint("simple", "k", "us-east-1")
	if ep.host != "simple.s3.us-east-1.amazonaws.com" || ep.canonicalURI != "/k" {
		t.Fatalf("forced virtual: host=%q uri=%q", ep.host, ep.canonicalURI)
	}
}

// --- GCS backend tests ------------------------------------------------------

func TestGCSReadTokenless(t *testing.T) {
	rs, srv := newRangeServer(payload)
	defer srv.Close()
	gcsBaseURL = srv.URL
	defer func() { gcsBaseURL = "https://storage.googleapis.com" }()
	t.Setenv("GCS_OAUTH_TOKEN", "")
	t.Setenv("GCS_REQUESTER_PAYS_PROJECT", "")

	h, err := Open("gs://public-bucket/data/file.cram")
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer h.Close()
	buf := make([]byte, 32)
	n, err := h.ReadAt(buf, 64)
	if err != nil {
		t.Fatalf("ReadAt: %v", err)
	}
	if !bytes.Equal(buf[:n], payload[64:64+32]) {
		t.Fatal("GCS bytes mismatch")
	}
	if rs.lastReq.Header.Get("Authorization") != "" {
		t.Fatal("unexpected Authorization header for tokenless GCS")
	}
	if !strings.Contains(rs.lastReq.URL.Path, "public-bucket/data/file.cram") {
		t.Fatalf("path=%q", rs.lastReq.URL.Path)
	}
}

func TestGCSBearerAndRequesterPays(t *testing.T) {
	rs, srv := newRangeServer(payload)
	defer srv.Close()
	gcsBaseURL = srv.URL
	defer func() { gcsBaseURL = "https://storage.googleapis.com" }()
	t.Setenv("GCS_OAUTH_TOKEN", "ya29.token")
	t.Setenv("GCS_REQUESTER_PAYS_PROJECT", "proj-123")

	h, _ := Open("gs://bucket/obj")
	defer h.Close()
	buf := make([]byte, 4)
	if _, err := h.ReadAt(buf, 0); err != nil {
		t.Fatalf("ReadAt: %v", err)
	}
	if got := rs.lastReq.Header.Get("Authorization"); got != "Bearer ya29.token" {
		t.Fatalf("Authorization=%q", got)
	}
	if got := rs.lastReq.Header.Get("X-Goog-User-Project"); got != "proj-123" {
		t.Fatalf("X-Goog-User-Project=%q", got)
	}
}

// --- SigV4 AWS test-vector --------------------------------------------------

// TestSigV4AWSVector verifies the SigV4 derivation against the AWS-published
// "GET /" example from the Signature Version 4 test suite (get-vanilla). The
// canonical request, string-to-sign and final signature equal the documented
// expected strings.
func TestSigV4AWSVector(t *testing.T) {
	// Published example inputs.
	const (
		accessKey = "AKIDEXAMPLE"
		secretKey = "wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY"
		region    = "us-east-1"
		service   = "service"
		host      = "example.amazonaws.com"
		amzDate   = "20150830T123600Z"
		dateShort = "20150830"
	)

	// Canonical request for "GET /" signing host and x-amz-date only,
	// with an empty payload.
	canonical := strings.Join([]string{
		"GET",
		"/",
		"",
		"host:" + host + "\n" + "x-amz-date:" + amzDate + "\n",
		"host;x-amz-date",
		emptyStringSHA256,
	}, "\n")

	const wantCanonicalHash = "bb579772317eb040ac9ed261061d46c1f17a8133879d6129b6e1c25292927e63"
	if got := hexSHA256(canonical); got != wantCanonicalHash {
		t.Fatalf("canonical hash = %s, want %s", got, wantCanonicalHash)
	}

	scope := strings.Join([]string{dateShort, region, service, "aws4_request"}, "/")
	stringToSign := strings.Join([]string{
		"AWS4-HMAC-SHA256",
		amzDate,
		scope,
		wantCanonicalHash,
	}, "\n")

	const wantStringToSign = "AWS4-HMAC-SHA256\n20150830T123600Z\n20150830/us-east-1/service/aws4_request\nbb579772317eb040ac9ed261061d46c1f17a8133879d6129b6e1c25292927e63"
	if stringToSign != wantStringToSign {
		t.Fatalf("string-to-sign mismatch:\n got %q\nwant %q", stringToSign, wantStringToSign)
	}

	kDate := hmacSHA256([]byte("AWS4"+secretKey), []byte(dateShort))
	kRegion := hmacSHA256(kDate, []byte(region))
	kService := hmacSHA256(kRegion, []byte(service))
	kSigning := hmacSHA256(kService, []byte("aws4_request"))
	sig := hmacSHA256(kSigning, []byte(stringToSign))

	const wantSignature = "ea21d6f05e96a897f6000a1a293f0a5bf0f92a00343409e820dce329ca6365ea"
	if got := toHex(sig); got != wantSignature {
		t.Fatalf("signature = %s, want %s", got, wantSignature)
	}

	// Drive signSigV4 itself with the same inputs but using its built-in
	// header set (which additionally signs x-amz-content-sha256); confirm it
	// produces a well-formed, stable Authorization header.
	res := signSigV4(sigV4Input{
		Method:      "GET",
		Host:        host,
		Region:      region,
		Service:     service,
		AccessKey:   accessKey,
		SecretKey:   secretKey,
		PayloadHash: emptyStringSHA256,
		Time:        time.Date(2015, 8, 30, 12, 36, 0, 0, time.UTC),
	})
	if res.AmzDate != amzDate {
		t.Fatalf("AmzDate=%q", res.AmzDate)
	}
	if !strings.Contains(res.Authorization, "Credential="+accessKey+"/"+scope) {
		t.Fatalf("Authorization=%q", res.Authorization)
	}
}

// --- SigV2 AWS known-answer vector -----------------------------------------

// TestSigV2AWSVector pins signSigV2's string-to-sign and signature to the
// canonical worked example from the AWS S3 REST authentication documentation:
// GET /johnsmith/photos/puppy.jpg, Date "Tue, 27 Mar 2007 19:36:42 +0000",
// key AKIAIOSFODNN7EXAMPLE / secret wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY,
// whose documented signature is "bWq2s1WEIj+Ydj0vQ697zp+IXMU=".
func TestSigV2AWSVector(t *testing.T) {
	const (
		accessKey = "AKIAIOSFODNN7EXAMPLE"
		secretKey = "wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY"
		date      = "Tue, 27 Mar 2007 19:36:42 +0000"
		resource  = "/johnsmith/photos/puppy.jpg"
		wantSig   = "bWq2s1WEIj+Ydj0vQ697zp+IXMU="
	)

	in := sigV2Input{Method: "GET", Resource: resource, AccessKey: accessKey, SecretKey: secretKey}

	// The string-to-sign uses the example's literal Date (two blank lines for
	// the empty Content-MD5 and Content-Type, no x-amz headers).
	sts := sigV2StringToSign(in, date)
	wantSTS := "GET\n\n\n" + date + "\n" + resource
	if sts != wantSTS {
		t.Fatalf("string-to-sign mismatch:\n got %q\nwant %q", sts, wantSTS)
	}

	gotSig := base64.StdEncoding.EncodeToString(hmacSHA1([]byte(secretKey), []byte(sts)))
	if gotSig != wantSig {
		t.Fatalf("signature = %q, want %q", gotSig, wantSig)
	}

	wantAuth := "AWS " + accessKey + ":" + wantSig
	if got := "AWS " + accessKey + ":" + gotSig; got != wantAuth {
		t.Fatalf("authorization = %q, want %q", got, wantAuth)
	}

	// Drive signSigV2 end-to-end with a pinned clock matching the example
	// instant; its GMT-form Date differs from the example's +0000 spelling, so
	// the signature is recomputed against that Date (round-trip self-check).
	res := signSigV2(sigV2Input{
		Method:    "GET",
		Resource:  resource,
		AccessKey: accessKey,
		SecretKey: secretKey,
		Time:      time.Date(2007, 3, 27, 19, 36, 42, 0, time.UTC),
	})
	if res.Date != "Tue, 27 Mar 2007 19:36:42 GMT" {
		t.Fatalf("Date=%q", res.Date)
	}
	wantGMTSig := base64.StdEncoding.EncodeToString(hmacSHA1([]byte(secretKey),
		[]byte(sigV2StringToSign(in, res.Date))))
	if res.Authorization != "AWS "+accessKey+":"+wantGMTSig {
		t.Fatalf("Authorization=%q", res.Authorization)
	}
}

// TestSigV2SessionToken verifies a session token is folded into the
// CanonicalizedAmzHeaders as x-amz-security-token (and surfaced for sending as
// a request header).
func TestSigV2SessionToken(t *testing.T) {
	const (
		secretKey = "secret"
		date      = "Tue, 27 Mar 2007 19:36:42 GMT"
		resource  = "/bucket/key"
		token     = "FQoGZXIvtoken"
	)
	in := sigV2Input{Method: "GET", Resource: resource, AccessKey: "AKID", SecretKey: secretKey, SessionToken: token}
	sts := sigV2StringToSign(in, date)
	want := "GET\n\n\n" + date + "\n" + "x-amz-security-token:" + token + "\n" + resource
	if sts != want {
		t.Fatalf("string-to-sign with token:\n got %q\nwant %q", sts, want)
	}
	res := signSigV2(sigV2Input{Method: "GET", Resource: resource, AccessKey: "AKID", SecretKey: secretKey, SessionToken: token, Time: time.Unix(0, 0)})
	if res.SecurityToken != token {
		t.Fatalf("SecurityToken=%q", res.SecurityToken)
	}
}

// TestS3ReadAndSignV2 is the httptest round-trip for SigV2: it confirms the
// request carries "Authorization: AWS <key>:<sig>" and a Date header, and that
// the signature matches what signSigV2 computes for the pinned clock.
func TestS3ReadAndSignV2(t *testing.T) {
	rs, srv := newRangeServer(payload)
	defer srv.Close()
	host := strings.TrimPrefix(srv.URL, "http://")
	t.Setenv("HTS_S3_V2", "1")
	t.Setenv("AWS_ACCESS_KEY_ID", "AKIAIOSFODNN7EXAMPLE")
	t.Setenv("AWS_SECRET_ACCESS_KEY", "wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY")
	t.Setenv("AWS_SESSION_TOKEN", "")
	s3HostOverride = "http://" + host
	defer func() { s3HostOverride = "" }()

	pinned := time.Date(2007, 3, 27, 19, 36, 42, 0, time.UTC)
	nowFunc = func() time.Time { return pinned }
	defer func() { nowFunc = time.Now }()

	h, err := Open("s3://johnsmith/photos/puppy.jpg")
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer h.Close()
	if _, err := h.ReadAt(make([]byte, 16), 0); err != nil {
		t.Fatalf("ReadAt: %v", err)
	}

	r := rs.lastReq
	wantSig := base64.StdEncoding.EncodeToString(hmacSHA1(
		[]byte("wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY"),
		[]byte(sigV2StringToSign(sigV2Input{Method: "GET", Resource: "/johnsmith/photos/puppy.jpg", AccessKey: "AKIAIOSFODNN7EXAMPLE", SecretKey: "x"},
			"Tue, 27 Mar 2007 19:36:42 GMT"))))
	wantAuth := "AWS AKIAIOSFODNN7EXAMPLE:" + wantSig
	if got := r.Header.Get("Authorization"); got != wantAuth {
		t.Fatalf("Authorization=%q want %q", got, wantAuth)
	}
	if got := r.Header.Get("Date"); got != "Tue, 27 Mar 2007 19:36:42 GMT" {
		t.Fatalf("Date=%q", got)
	}
	// Path-style URL carries the bucket.
	if !strings.Contains(r.URL.Path, "johnsmith/photos/puppy.jpg") {
		t.Fatalf("path=%q", r.URL.Path)
	}
}

// --- GCS S3-interop tests ---------------------------------------------------

// TestS3AgainstGCSV4 simulates an s3://bucket/key request routed to the Google
// Cloud Storage S3-interop endpoint with HMAC keys, signed with SigV4. It
// asserts path-style addressing ("/bucket/key") and a SigV4 Authorization.
func TestS3AgainstGCSV4(t *testing.T) {
	rs, srv := newRangeServer(payload)
	defer srv.Close()
	host := strings.TrimPrefix(srv.URL, "http://")
	// Stand in for storage.googleapis.com; path-style is forced by HTS_S3_HOST.
	t.Setenv("HTS_S3_V2", "")
	t.Setenv("AWS_ACCESS_KEY_ID", "GOOG1EXAMPLEHMACKEY")
	t.Setenv("AWS_SECRET_ACCESS_KEY", "gcs-hmac-secret")
	t.Setenv("AWS_SESSION_TOKEN", "")
	t.Setenv("AWS_DEFAULT_REGION", "auto")
	s3HostOverride = "http://" + host
	defer func() { s3HostOverride = "" }()

	pinned := time.Date(2024, 1, 2, 3, 4, 5, 0, time.UTC)
	nowFunc = func() time.Time { return pinned }
	defer func() { nowFunc = time.Now }()

	h, err := Open("s3://gcs-bucket/data/file.cram")
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer h.Close()
	if _, err := h.ReadAt(make([]byte, 16), 0); err != nil {
		t.Fatalf("ReadAt: %v", err)
	}
	r := rs.lastReq
	if !strings.HasPrefix(r.Header.Get("Authorization"), "AWS4-HMAC-SHA256 ") {
		t.Fatalf("Authorization=%q", r.Header.Get("Authorization"))
	}
	if !strings.Contains(r.Header.Get("Authorization"), "Credential=GOOG1EXAMPLEHMACKEY/") {
		t.Fatalf("missing credential: %q", r.Header.Get("Authorization"))
	}
	if r.Header.Get("x-amz-date") == "" {
		t.Fatal("missing x-amz-date")
	}
	if r.URL.Path != "/gcs-bucket/data/file.cram" {
		t.Fatalf("path=%q (expected path-style /bucket/key)", r.URL.Path)
	}
}

// TestS3AgainstGCSV2 is the SigV2 (legacy HMAC) variant of the GCS S3-interop
// path: same endpoint, signed with AWS Signature Version 2.
func TestS3AgainstGCSV2(t *testing.T) {
	rs, srv := newRangeServer(payload)
	defer srv.Close()
	host := strings.TrimPrefix(srv.URL, "http://")
	t.Setenv("HTS_S3_V2", "1")
	t.Setenv("AWS_ACCESS_KEY_ID", "GOOG1EXAMPLEHMACKEY")
	t.Setenv("AWS_SECRET_ACCESS_KEY", "gcs-hmac-secret")
	t.Setenv("AWS_SESSION_TOKEN", "")
	s3HostOverride = "http://" + host
	defer func() { s3HostOverride = "" }()

	pinned := time.Date(2024, 1, 2, 3, 4, 5, 0, time.UTC)
	nowFunc = func() time.Time { return pinned }
	defer func() { nowFunc = time.Now }()

	h, err := Open("s3://gcs-bucket/data/file.cram")
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer h.Close()
	if _, err := h.ReadAt(make([]byte, 16), 0); err != nil {
		t.Fatalf("ReadAt: %v", err)
	}
	r := rs.lastReq
	wantSig := base64.StdEncoding.EncodeToString(hmacSHA1(
		[]byte("gcs-hmac-secret"),
		[]byte(sigV2StringToSign(sigV2Input{Method: "GET", Resource: "/gcs-bucket/data/file.cram"},
			pinned.UTC().Format(sigV2DateFormat)))))
	wantAuth := "AWS GOOG1EXAMPLEHMACKEY:" + wantSig
	if got := r.Header.Get("Authorization"); got != wantAuth {
		t.Fatalf("Authorization=%q want %q", got, wantAuth)
	}
	if r.Header.Get("Date") == "" {
		t.Fatal("missing Date header")
	}
	if r.URL.Path != "/gcs-bucket/data/file.cram" {
		t.Fatalf("path=%q (expected path-style /bucket/key)", r.URL.Path)
	}
}

func toHex(b []byte) string {
	const hexdigits = "0123456789abcdef"
	out := make([]byte, len(b)*2)
	for i, c := range b {
		out[i*2] = hexdigits[c>>4]
		out[i*2+1] = hexdigits[c&0xf]
	}
	return string(out)
}

// --- creds tests ------------------------------------------------------------

func TestCredentialsFileProfileSelection(t *testing.T) {
	dir := t.TempDir()
	credFile := filepath.Join(dir, "credentials")
	content := `[default]
aws_access_key_id = DEFAULTKEY
aws_secret_access_key = defaultsecret
region = us-west-2

[work]
aws_access_key_id = WORKKEY
aws_secret_access_key = worksecret
aws_session_token = worktoken
region = eu-central-1
`
	if err := os.WriteFile(credFile, []byte(content), 0o600); err != nil {
		t.Fatal(err)
	}
	t.Setenv("AWS_SHARED_CREDENTIALS_FILE", credFile)
	// Clear env so the file is consulted.
	t.Setenv("AWS_ACCESS_KEY_ID", "")
	t.Setenv("AWS_SECRET_ACCESS_KEY", "")
	t.Setenv("AWS_SESSION_TOKEN", "")
	t.Setenv("AWS_DEFAULT_REGION", "")
	t.Setenv("AWS_REGION", "")

	t.Setenv("AWS_PROFILE", "work")
	t.Setenv("AWS_DEFAULT_PROFILE", "")
	c := resolveAWSCredentials()
	if c.AccessKey != "WORKKEY" || c.SecretKey != "worksecret" || c.SessionToken != "worktoken" {
		t.Fatalf("work profile: %+v", c)
	}
	if c.Region != "eu-central-1" {
		t.Fatalf("region=%q", c.Region)
	}

	t.Setenv("AWS_PROFILE", "")
	c = resolveAWSCredentials()
	if c.AccessKey != "DEFAULTKEY" || c.Region != "us-west-2" {
		t.Fatalf("default profile: %+v", c)
	}
}

func TestCredentialsEnvPrecedence(t *testing.T) {
	dir := t.TempDir()
	credFile := filepath.Join(dir, "credentials")
	if err := os.WriteFile(credFile, []byte("[default]\naws_access_key_id = FILEKEY\naws_secret_access_key = filesecret\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	t.Setenv("AWS_SHARED_CREDENTIALS_FILE", credFile)
	t.Setenv("AWS_ACCESS_KEY_ID", "ENVKEY")
	t.Setenv("AWS_SECRET_ACCESS_KEY", "envsecret")
	t.Setenv("AWS_SESSION_TOKEN", "")
	t.Setenv("AWS_PROFILE", "")
	t.Setenv("AWS_DEFAULT_PROFILE", "")

	c := resolveAWSCredentials()
	if c.AccessKey != "ENVKEY" || c.SecretKey != "envsecret" {
		t.Fatalf("env should win: %+v", c)
	}
}

// --- local-file tests -------------------------------------------------------

func TestLocalHandleParity(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "data.bin")
	if err := os.WriteFile(p, payload, 0o600); err != nil {
		t.Fatal(err)
	}

	for _, name := range []string{p, "file://" + p} {
		h, err := Open(name)
		if err != nil {
			t.Fatalf("Open(%q): %v", name, err)
		}
		sz, err := h.Size()
		if err != nil || sz != int64(len(payload)) {
			t.Fatalf("Size=%d err=%v", sz, err)
		}
		buf := make([]byte, 50)
		ref := make([]byte, 50)
		f, _ := os.Open(p)
		for _, off := range []int64{0, 200, int64(len(payload) - 25)} {
			n1, e1 := h.ReadAt(buf, off)
			n2, e2 := f.ReadAt(ref, off)
			if n1 != n2 || !errors.Is(e1, e2) && (e1 == nil) != (e2 == nil) {
				t.Fatalf("ReadAt off=%d: hfile(n=%d,e=%v) os(n=%d,e=%v)", off, n1, e1, n2, e2)
			}
			if !bytes.Equal(buf[:n1], ref[:n2]) {
				t.Fatalf("ReadAt off=%d bytes differ", off)
			}
		}
		f.Close()
		h.Close()
	}
}

func TestSchemeHelpers(t *testing.T) {
	cases := []struct {
		name     string
		scheme   string
		isRemote bool
	}{
		{"http://x/y", "http", true},
		{"https://x/y", "https", true},
		{"s3://b/k", "s3", true},
		{"gs://b/o", "gs", true},
		{"file:///tmp/x", "file", false},
		{"/tmp/x", "", false},
		{"relative/path", "", false},
		{"C:/win/path", "", false},
	}
	for _, c := range cases {
		if got := SchemeOf(c.name); got != c.scheme {
			t.Errorf("SchemeOf(%q)=%q want %q", c.name, got, c.scheme)
		}
		if got := IsRemote(c.name); got != c.isRemote {
			t.Errorf("IsRemote(%q)=%v want %v", c.name, got, c.isRemote)
		}
	}
}

func TestOpenUnsupportedScheme(t *testing.T) {
	if _, err := Open("ftp://host/x"); err == nil {
		t.Fatal("expected error for ftp scheme")
	}
}

// --- additional coverage ----------------------------------------------------

func TestLocalSequentialRead(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "f.bin")
	if err := os.WriteFile(p, payload, 0o600); err != nil {
		t.Fatal(err)
	}
	h, err := Open(p)
	if err != nil {
		t.Fatal(err)
	}
	defer h.Close()
	got, err := io.ReadAll(h)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, payload) {
		t.Fatal("local sequential read mismatch")
	}
}

func TestReadAtZeroAndNegative(t *testing.T) {
	_, srv := newRangeServer(payload)
	defer srv.Close()
	h, _ := Open(srv.URL)
	defer h.Close()

	n, err := h.ReadAt(nil, 0)
	if n != 0 || err != nil {
		t.Fatalf("zero-length ReadAt: n=%d err=%v", n, err)
	}
	if _, err := h.ReadAt(make([]byte, 4), -1); err == nil {
		t.Fatal("expected error for negative offset")
	}
}

func TestRetryExhausted(t *testing.T) {
	t.Setenv("HTS_RETRY_MAX", "1")
	t.Setenv("HTS_RETRY_DELAY", "1")
	t.Setenv("HTS_RETRY_MAX_DELAY", "5")
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.WriteHeader(http.StatusServiceUnavailable)
	}))
	defer srv.Close()
	h, _ := Open(srv.URL)
	defer h.Close()
	if _, err := h.ReadAt(make([]byte, 4), 0); err == nil {
		t.Fatal("expected error after retries exhausted")
	}
}

func TestCurlCABundleError(t *testing.T) {
	dir := t.TempDir()
	bad := filepath.Join(dir, "bad.pem")
	if err := os.WriteFile(bad, []byte("not a certificate"), 0o600); err != nil {
		t.Fatal(err)
	}
	t.Setenv("CURL_CA_BUNDLE", bad)
	if _, err := Open("https://example.com/x"); err == nil {
		t.Fatal("expected error for invalid CA bundle")
	}
}

func TestURLParseErrors(t *testing.T) {
	t.Setenv("HTS_S3_V2", "")
	bad := []string{"s3://bucket", "s3:///key", "s3://bucket/", "gs://bucket", "gs://bucket/"}
	for _, u := range bad {
		if _, err := Open(u); err == nil {
			t.Errorf("Open(%q) should fail", u)
		}
	}
}

func TestEncodeS3Path(t *testing.T) {
	got := encodeS3Path("dir/a b+c.bam")
	if got != "dir/a%20b%2Bc.bam" {
		t.Fatalf("encodeS3Path=%q", got)
	}
}

func TestS3HostOverrideVirtualDefault(t *testing.T) {
	s3HostOverride = ""
	t.Setenv("HTS_S3_V2", "")
	t.Setenv("AWS_ACCESS_KEY_ID", "")
	t.Setenv("AWS_SECRET_ACCESS_KEY", "")
	// No override: virtual-hosted default endpoint is used (no network call).
	h, err := Open("s3://bucket/key")
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer h.Close()
	s3h := h.(*s3Handle)
	if !strings.Contains(s3h.url, "bucket.s3.us-east-1.amazonaws.com") {
		t.Fatalf("url=%q", s3h.url)
	}
}

func TestS3AnonymousUnsigned(t *testing.T) {
	rs, srv := newRangeServer(payload)
	defer srv.Close()
	host := strings.TrimPrefix(srv.URL, "http://")
	s3HostOverride = "http://" + host
	defer func() { s3HostOverride = "" }()
	t.Setenv("AWS_ACCESS_KEY_ID", "")
	t.Setenv("AWS_SECRET_ACCESS_KEY", "")
	t.Setenv("AWS_SHARED_CREDENTIALS_FILE", filepath.Join(t.TempDir(), "none"))

	h, err := Open("s3://public/obj")
	if err != nil {
		t.Fatal(err)
	}
	defer h.Close()
	if _, err := h.ReadAt(make([]byte, 4), 0); err != nil {
		t.Fatalf("ReadAt: %v", err)
	}
	if rs.lastReq.Header.Get("Authorization") != "" {
		t.Fatal("anonymous request should be unsigned")
	}
}

// TestReadFileNonRangeServer verifies ReadFile is correct even when the server
// ignores Range and returns the whole body for every request (200). The ranged
// ReadAt path would mis-stitch such responses; readWhole avoids it.
func TestReadFileNonRangeServer(t *testing.T) {
	// Large enough that io.ReadAll over ranged ReadAt would do several reads.
	body := bytes.Repeat([]byte("REFERENCE"), 100000)
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		// Deliberately ignore Range: always write the whole body with 200.
		w.WriteHeader(http.StatusOK)
		w.Write(body)
	}))
	defer srv.Close()

	got, err := ReadFile(srv.URL)
	if err != nil {
		t.Fatalf("ReadFile(non-range server): %v", err)
	}
	if !bytes.Equal(got, body) {
		t.Fatalf("ReadFile(non-range server): got %d bytes, want %d", len(got), len(body))
	}
}
