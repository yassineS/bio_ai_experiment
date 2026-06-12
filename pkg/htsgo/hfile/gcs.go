package hfile

import (
	"fmt"
	"net/http"
	"os"
	"strings"
)

// gcsBaseURL is the base URL of the Google Cloud Storage XML API, to which
// "/<bucket>/<object>" is appended. The XML API supports HTTP Range requests,
// which is what the index path needs. Tests override this to point the backend
// at an httptest server.
var gcsBaseURL = "https://storage.googleapis.com"

// gcsHandle is a read-only Handle for objects addressed by a gs:// URL. It
// delegates HTTP transport to an httpHandle, optionally adding a bearer token
// and a requester-pays project header.
type gcsHandle struct {
	*httpHandle
}

// openGCS opens a gs://bucket/object URL for reading. If GCS_OAUTH_TOKEN is
// set, an "Authorization: Bearer <token>" header is added; public buckets work
// without it. If GCS_REQUESTER_PAYS_PROJECT is set, the X-Goog-User-Project
// header is added, matching hfile_gcs.c.
func openGCS(rawurl string) (Handle, error) {
	bucket, object, err := parseGCSURL(rawurl)
	if err != nil {
		return nil, err
	}

	fullURL := strings.TrimSuffix(gcsBaseURL, "/") + "/" + bucket + "/" + encodeS3Path(object)

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

	if token := os.Getenv("GCS_OAUTH_TOKEN"); token != "" {
		hh.headers.Set("Authorization", "Bearer "+token)
	}
	if project := os.Getenv("GCS_REQUESTER_PAYS_PROJECT"); project != "" {
		hh.headers.Set("X-Goog-User-Project", project)
	}

	return &gcsHandle{httpHandle: hh}, nil
}

// parseGCSURL splits a gs://bucket/object URL into its bucket and object
// components.
func parseGCSURL(rawurl string) (bucket, object string, err error) {
	rest := strings.TrimPrefix(rawurl, "gs://")
	if rest == rawurl {
		return "", "", fmt.Errorf("hfile: not a gs:// URL: %q", rawurl)
	}
	slash := strings.IndexByte(rest, '/')
	if slash < 0 {
		return "", "", fmt.Errorf("hfile: gs URL %q has no object", rawurl)
	}
	bucket = rest[:slash]
	object = rest[slash+1:]
	if bucket == "" || object == "" {
		return "", "", fmt.Errorf("hfile: gs URL %q must be gs://bucket/object", rawurl)
	}
	return bucket, object, nil
}
