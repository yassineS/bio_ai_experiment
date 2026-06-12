package hfile

import (
	"crypto/tls"
	"crypto/x509"
	"errors"
	"fmt"
	"io"
	"net/http"
	"os"
	"strconv"
	"strings"
	"sync"
	"time"
)

// retryConfig holds the retry behaviour controlled by the HTS_RETRY_*
// environment variables, matching hfile_libcurl.c's defaults.
type retryConfig struct {
	max      int           // HTS_RETRY_MAX, default 3
	delay    time.Duration // HTS_RETRY_DELAY (ms), default 500ms
	maxDelay time.Duration // HTS_RETRY_MAX_DELAY (ms), default 60000ms
}

// loadRetryConfig reads the HTS_RETRY_* environment variables, falling back
// to htslib's defaults when they are unset or malformed.
func loadRetryConfig() retryConfig {
	rc := retryConfig{max: 3, delay: 500 * time.Millisecond, maxDelay: 60000 * time.Millisecond}
	if v := os.Getenv("HTS_RETRY_MAX"); v != "" {
		if n, err := strconv.Atoi(v); err == nil {
			rc.max = n
		}
	}
	if v := os.Getenv("HTS_RETRY_DELAY"); v != "" {
		if n, err := strconv.ParseInt(v, 10, 64); err == nil {
			rc.delay = time.Duration(n) * time.Millisecond
		}
	}
	if v := os.Getenv("HTS_RETRY_MAX_DELAY"); v != "" {
		if n, err := strconv.ParseInt(v, 10, 64); err == nil {
			rc.maxDelay = time.Duration(n) * time.Millisecond
		}
	}
	return rc
}

// requestSigner is an optional hook that lets the S3 backend add
// authentication headers to each outgoing request just before it is sent.
// It receives the *http.Request and may set headers on it.
type requestSigner func(req *http.Request) error

// httpHandle is a read-only Handle backed by ranged HTTP GET requests.
type httpHandle struct {
	url     string
	client  *http.Client
	retry   retryConfig
	sign    requestSigner // optional, used by the S3 backend
	headers http.Header   // static headers added to every request (GCS auth, etc.)

	mu       sync.Mutex // guards seqOff for sequential Read
	seqOff   int64
	sizeOnce sync.Once
	size     int64
	sizeErr  error
}

// openHTTP opens an http:// or https:// URL for reading.
func openHTTP(rawurl string) (Handle, error) {
	client, err := newHTTPClient()
	if err != nil {
		return nil, err
	}
	h := &httpHandle{
		url:     rawurl,
		client:  client,
		retry:   loadRetryConfig(),
		headers: http.Header{},
	}
	if err := applyAuthLocation(h.headers); err != nil {
		return nil, err
	}
	return h, nil
}

// newHTTPClient builds an *http.Client honouring CURL_CA_BUNDLE for a custom
// trusted CA set. The default client follows redirects, matching libcurl.
func newHTTPClient() (*http.Client, error) {
	transport := &http.Transport{
		Proxy: http.ProxyFromEnvironment,
	}
	if bundle := os.Getenv("CURL_CA_BUNDLE"); bundle != "" {
		pem, err := os.ReadFile(bundle)
		if err != nil {
			return nil, fmt.Errorf("hfile: CURL_CA_BUNDLE: %w", err)
		}
		pool := x509.NewCertPool()
		if !pool.AppendCertsFromPEM(pem) {
			return nil, fmt.Errorf("hfile: CURL_CA_BUNDLE %q contains no usable certificates", bundle)
		}
		transport.TLSClientConfig = &tls.Config{RootCAs: pool}
	}
	return &http.Client{Transport: transport}, nil
}

// applyAuthLocation honours HTS_AUTH_LOCATION: if set, the contents of the
// named file are sent verbatim as the Authorization header, matching
// hfile_libcurl.c.
func applyAuthLocation(hdr http.Header) error {
	loc := os.Getenv("HTS_AUTH_LOCATION")
	if loc == "" {
		return nil
	}
	data, err := os.ReadFile(loc)
	if err != nil {
		return fmt.Errorf("hfile: HTS_AUTH_LOCATION: %w", err)
	}
	hdr.Set("Authorization", strings.TrimSpace(string(data)))
	return nil
}

// Read implements io.Reader on top of ReadAt, maintaining a sequential
// offset. It is goroutine-safe with respect to the sequential offset, but
// concurrent sequential reads will interleave; use ReadAt for concurrent
// random access.
func (h *httpHandle) Read(p []byte) (int, error) {
	h.mu.Lock()
	off := h.seqOff
	h.mu.Unlock()

	n, err := h.ReadAt(p, off)
	if n > 0 {
		h.mu.Lock()
		h.seqOff = off + int64(n)
		h.mu.Unlock()
	}
	return n, err
}

// ReadAt issues a ranged GET for bytes [off, off+len(p)). It obeys the
// os.File.ReadAt contract and is safe for concurrent use. A short read at
// the end of the resource returns io.EOF.
func (h *httpHandle) ReadAt(p []byte, off int64) (int, error) {
	if len(p) == 0 {
		return 0, nil
	}
	if off < 0 {
		return 0, errors.New("hfile: ReadAt: negative offset")
	}

	last := off + int64(len(p)) - 1
	rangeHdr := fmt.Sprintf("bytes=%d-%d", off, last)

	resp, err := h.do(http.MethodGet, map[string]string{"Range": rangeHdr})
	if err != nil {
		return 0, err
	}
	defer resp.Body.Close()

	switch resp.StatusCode {
	case http.StatusPartialContent, http.StatusOK:
		// 206 is the expected ranged response; 200 means the server
		// ignored the Range and returned the whole body (small files).
	case http.StatusRequestedRangeNotSatisfiable:
		// Reading at or past EOF.
		return 0, io.EOF
	default:
		return 0, fmt.Errorf("hfile: GET %s: unexpected status %s", h.url, resp.Status)
	}

	n, err := io.ReadFull(resp.Body, p)
	if err == io.ErrUnexpectedEOF || err == io.EOF {
		// Fewer bytes than requested: the range ran past end of resource.
		return n, io.EOF
	}
	if err != nil {
		return n, err
	}
	return n, nil
}

// Size returns the total size of the resource. It first attempts a HEAD
// request and falls back to a "bytes=0-0" ranged GET, parsing the total
// from the Content-Range header. The result is cached.
func (h *httpHandle) Size() (int64, error) {
	h.sizeOnce.Do(func() {
		h.size, h.sizeErr = h.fetchSize()
	})
	return h.size, h.sizeErr
}

func (h *httpHandle) fetchSize() (int64, error) {
	// Try HEAD first.
	if resp, err := h.do(http.MethodHead, nil); err == nil {
		size := int64(-1)
		if resp.StatusCode == http.StatusOK && resp.ContentLength >= 0 {
			size = resp.ContentLength
		}
		io.Copy(io.Discard, resp.Body)
		resp.Body.Close()
		if size >= 0 {
			return size, nil
		}
	}

	// Fall back to a 1-byte ranged GET and read the total from Content-Range.
	resp, err := h.do(http.MethodGet, map[string]string{"Range": "bytes=0-0"})
	if err != nil {
		return 0, err
	}
	defer resp.Body.Close()
	io.Copy(io.Discard, resp.Body)

	if cr := resp.Header.Get("Content-Range"); cr != "" {
		// Format: "bytes 0-0/12345"
		if slash := strings.LastIndex(cr, "/"); slash >= 0 {
			total := strings.TrimSpace(cr[slash+1:])
			if total != "*" {
				if n, err := strconv.ParseInt(total, 10, 64); err == nil {
					return n, nil
				}
			}
		}
	}
	if resp.StatusCode == http.StatusOK && resp.ContentLength >= 0 {
		return resp.ContentLength, nil
	}
	return 0, fmt.Errorf("hfile: cannot determine size of %s", h.url)
}

// Close releases backend resources. The HTTP handle holds no persistent
// connection of its own, so this is a no-op.
func (h *httpHandle) Close() error {
	return nil
}

// do performs an HTTP request for the handle's URL with the given method
// and extra headers, retrying on transient failures with exponential
// backoff per the HTS_RETRY_* configuration.
func (h *httpHandle) do(method string, extra map[string]string) (*http.Response, error) {
	var lastErr error
	delay := h.retry.delay

	attempts := h.retry.max + 1
	if attempts < 1 {
		attempts = 1
	}

	for attempt := 0; attempt < attempts; attempt++ {
		if attempt > 0 {
			time.Sleep(delay)
			delay *= 2
			if delay > h.retry.maxDelay {
				delay = h.retry.maxDelay
			}
		}

		req, err := http.NewRequest(method, h.url, nil)
		if err != nil {
			return nil, err
		}
		for k, vals := range h.headers {
			for _, v := range vals {
				req.Header.Add(k, v)
			}
		}
		for k, v := range extra {
			req.Header.Set(k, v)
		}
		if h.sign != nil {
			if err := h.sign(req); err != nil {
				return nil, err
			}
		}

		resp, err := h.client.Do(req)
		if err != nil {
			lastErr = err
			continue // transient network error: retry
		}

		if isRetryableStatus(resp.StatusCode) {
			lastErr = fmt.Errorf("hfile: %s %s: status %s", method, h.url, resp.Status)
			resp.Body.Close()
			continue
		}
		return resp, nil
	}
	return nil, lastErr
}

// isRetryableStatus reports whether an HTTP status warrants a retry,
// matching the 429/5xx set used by hfile_libcurl.c's is_retryable.
func isRetryableStatus(code int) bool {
	switch code {
	case http.StatusTooManyRequests, // 429
		http.StatusInternalServerError, // 500
		http.StatusBadGateway,          // 502
		http.StatusServiceUnavailable,  // 503
		http.StatusGatewayTimeout:      // 504
		return true
	default:
		return false
	}
}

var _ Handle = (*httpHandle)(nil)
