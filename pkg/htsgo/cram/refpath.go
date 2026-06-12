package cram

import (
	"encoding/hex"
	"fmt"
	"os"
	"strings"
	"sync"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/hfile"
)

// ebiRefPathURL is htslib's built-in default REF_PATH endpoint: the EBI ENA
// CRAM reference registry, which serves a reference sequence addressed by the
// hex of its base-only MD5 digest. htslib uses this when REF_PATH is unset.
const ebiRefPathURL = "https://www.ebi.ac.uk/ena/cram/md5/%s"

// RefPathReference resolves a whole reference sequence by its MD5 digest from
// one or more network URL templates, mirroring htslib's REF_PATH URL-fetch
// mechanism. Each template contains a single "%s" placeholder that is filled
// with the lower-case hex MD5; the response body is the raw reference bases.
//
// Fetched sequences are upper-cased (so they hash and slice exactly like the
// FASTA and REF_CACHE sources) and memoised by digest, so a run of slices that
// share a contig triggers at most one download.
type RefPathReference struct {
	urls []string

	mu     sync.Mutex
	cached map[[16]byte][]byte
}

// NewRefPathReference returns a RefPathReference that consults the given URL
// templates in order. Each template must contain exactly one "%s". It returns
// nil when no templates are supplied.
func NewRefPathReference(urls []string) *RefPathReference {
	if len(urls) == 0 {
		return nil
	}
	return &RefPathReference{urls: urls, cached: map[[16]byte][]byte{}}
}

// RefPathFromEnv builds a RefPathReference from the REF_PATH environment
// variable. Network reference fetching is opt-in: it is enabled only when
// REF_PATH is set and contains at least one URL template (an entry with a
// "://" scheme). A bare "%s"-less URL entry has "/%s" appended so a directory-
// style endpoint still addresses by digest, matching htslib. When REF_PATH is
// unset, ok is false and the caller keeps its existing behaviour (an explicit
// FASTA / REF_CACHE, or the 'N' fallback) rather than silently reaching out to
// the network; a user opts in by setting REF_PATH (e.g. to the EBI URL
// available as DefaultRefPathURL).
func RefPathFromEnv() (*RefPathReference, bool) {
	raw := os.Getenv("REF_PATH")
	if raw == "" {
		return nil, false
	}
	var urls []string
	for _, entry := range splitRefPath(raw) {
		if !strings.Contains(entry, "://") {
			continue // a local directory entry: handled by REF_CACHE, not here.
		}
		if !strings.Contains(entry, "%s") {
			entry = strings.TrimRight(entry, "/") + "/%s"
		}
		urls = append(urls, entry)
	}
	if len(urls) == 0 {
		return nil, false
	}
	return NewRefPathReference(urls), true
}

// DefaultRefPathURL is htslib's built-in EBI REF_PATH endpoint, exposed so a
// caller (or a user via REF_PATH) can opt into network reference resolution.
const DefaultRefPathURL = ebiRefPathURL

// splitRefPath splits a REF_PATH value into its entries. htslib separates
// entries with ':' but a URL also contains colons — the "scheme://" delimiter
// and an optional ":port" in the authority — that must not be treated as
// separators. The scanner therefore suppresses ':' splitting from a
// "scheme://" through to the end of the authority (the next '/'), so
// "https://host:8080/ref/%s" stays one entry while a top-level ':' between
// entries still splits.
func splitRefPath(s string) []string {
	var entries []string
	var cur strings.Builder
	inAuthority := false // between "scheme://" and the authority-terminating '/'
	for i := 0; i < len(s); i++ {
		c := s[i]
		switch {
		case !inAuthority && c == ':' && i+2 < len(s) && s[i+1] == '/' && s[i+2] == '/':
			// The scheme colon: not a separator; enter the URL authority.
			cur.WriteString("://")
			i += 2
			inAuthority = true
		case !inAuthority && c == ':':
			entries = append(entries, cur.String())
			cur.Reset()
		case inAuthority && c == '/':
			inAuthority = false
			cur.WriteByte(c)
		default:
			cur.WriteByte(c)
		}
	}
	entries = append(entries, cur.String())
	out := entries[:0]
	for _, e := range entries {
		if e != "" {
			out = append(out, e)
		}
	}
	return out
}

// ResolveByMD5 downloads and returns the full reference sequence whose
// upper-cased, base-only MD5 digest is the given 16-byte value. Each configured
// URL template is tried in order; the first that returns a body wins. The bytes
// are upper-cased and memoised by digest. Every attempt's error is reported on
// total failure so the caller can see which endpoints were tried.
func (r *RefPathReference) ResolveByMD5(digest [16]byte) ([]byte, error) {
	r.mu.Lock()
	if bases, ok := r.cached[digest]; ok {
		r.mu.Unlock()
		return bases, nil
	}
	r.mu.Unlock()

	hexStr := hex.EncodeToString(digest[:])
	var errs []string
	for _, tmpl := range r.urls {
		url := fmt.Sprintf(tmpl, hexStr)
		body, err := hfile.ReadFile(url)
		if err != nil {
			errs = append(errs, fmt.Sprintf("%s: %v", url, err))
			continue
		}
		bases := normalizeReferenceBases(body)
		r.mu.Lock()
		r.cached[digest] = bases
		r.mu.Unlock()
		return bases, nil
	}
	return nil, errFormat("REF_PATH could not fetch reference MD5 %s: %s",
		hexStr, strings.Join(errs, "; "))
}
