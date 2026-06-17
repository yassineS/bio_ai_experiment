# hfile

`pkg/htsgo/hfile` is a pure-Go port of htslib's `hfile` virtual-filesystem
remote backends. It provides **read** access to local files and to remote
objects addressed by HTTP(S), Amazon S3 (`s3://`) and Google Cloud Storage
(`gs://`) URLs, mirroring the behaviour of htslib's `hfile_libcurl.c`,
`hfile_s3.c` and `hfile_gcs.c`.

The package depends only on the Go standard library. AWS Signature Version 4
is hand-implemented with `crypto/hmac` + `crypto/sha256`, and AWS Signature
Version 2 with `crypto/hmac` + `crypto/sha1` + `encoding/base64`; there is no
AWS SDK and no OAuth library.

## API

```go
h, err := hfile.Open(name)   // dispatches by URL scheme
if err != nil { ... }
defer h.Close()

size, _ := h.Size()          // total size in bytes

buf := make([]byte, 65536)
n, err := h.ReadAt(buf, off) // ranged random access (concurrency-safe)

io.Copy(dst, h)              // sequential read on top of ReadAt
```

`Handle` is `io.Reader` + `io.ReaderAt` + `io.Closer` + `Size() (int64, error)`.

`ReadAt` follows the `os.File.ReadAt` contract precisely: it returns a non-nil
error whenever `n < len(p)`, and returns `io.EOF` once the end of the resource
is reached. This is the contract `bgzf.NewSeekReader` relies on, and each
`ReadAt` call issues an independent ranged GET with no shared mutable offset,
so it is safe for the concurrent index-read path.

Helpers:

- `hfile.IsRemote(name) bool` — true for `http`/`https`/`s3`/`gs` URLs.
- `hfile.SchemeOf(name) string` — the lower-cased URL scheme, or `""` for a
  bare filesystem path.

## Supported schemes

| Scheme              | Backend | Endpoint                                                            |
| ------------------- | ------- | ------------------------------------------------------------------ |
| `http://` `https://`| HTTP    | the URL as given                                                   |
| `s3://bucket/key`   | S3      | `https://<bucket>.s3.<region>.amazonaws.com/<key>` (virtual-hosted)|
| `gs://bucket/object`| GCS     | `https://storage.googleapis.com/<bucket>/<object>` (XML API)       |
| `file://` / bare path | local | `*os.File` (`-` means stdin)                                      |

A bucket name containing a dot forces S3 **path-style**
(`https://s3.<region>.amazonaws.com/<bucket>/<key>`), because dotted names break
TLS virtual-hosted certificates — matching htslib. `HTS_S3_ADDRESS_STYLE` can
force `path` or `virtual` explicitly.

Write access is **not implemented**; remote handles are read-only.

## Cloud / S3-compatible provider matrix

The `s3://` backend works against AWS and any S3-compatible object store. The
default signer is **AWS Signature Version 4**; set `HTS_S3_V2` to use **AWS
Signature Version 2** (HMAC-SHA1), which some older or self-hosted gateways
require. Google Cloud Storage is reachable two ways: the native `gs://` backend
(OAuth bearer token) or the S3-interop XML API via `s3://` + HMAC keys.

| Provider                         | URL / endpoint                                              | Signing                                  | Required env vars |
| -------------------------------- | ---------------------------------------------------------- | ---------------------------------------- | ----------------- |
| **AWS S3**                       | `s3://bucket/key` (virtual-hosted by default)              | SigV4 (default)                          | `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_DEFAULT_REGION` (optional, default `us-east-1`); `AWS_SESSION_TOKEN` for STS |
| **MinIO / Ceph (RGW)**           | `s3://bucket/key` + `HTS_S3_HOST=minio.example:9000`       | SigV4 (default) or SigV2 (`HTS_S3_V2`)   | `HTS_S3_HOST`, `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`; `HTS_S3_V2=1` for legacy V2 endpoints |
| **Wasabi**                       | `s3://bucket/key` + `HTS_S3_HOST=s3.<region>.wasabisys.com`| SigV4 (default)                          | `HTS_S3_HOST`, `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY` |
| **Backblaze B2 (S3 API)**        | `s3://bucket/key` + `HTS_S3_HOST=s3.<region>.backblazeb2.com` | SigV4 (default)                       | `HTS_S3_HOST`, `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY` |
| **GCS (native)**                 | `gs://bucket/object`                                       | OAuth Bearer (or anonymous)              | `GCS_OAUTH_TOKEN` (optional for public buckets); `GCS_REQUESTER_PAYS_PROJECT` for requester-pays |
| **GCS (S3 interop)**             | `s3://bucket/key` + `HTS_S3_HOST=storage.googleapis.com`   | SigV4 (default) or SigV2 (`HTS_S3_V2`)   | `HTS_S3_HOST=storage.googleapis.com`, `AWS_ACCESS_KEY_ID`/`AWS_SECRET_ACCESS_KEY` set to a GCS **HMAC key** (`GOOG...` / secret) |

For any custom endpoint, `HTS_S3_HOST` forces **path-style** addressing
(`<host>/<bucket>/<key>`), as is conventional for S3-compatible gateways. SigV4
still places a region in the credential scope; for non-AWS endpoints any literal
(default `us-east-1`, or `auto` for GCS) is fine because the server validates
the HMAC, not the region string. SigV2 ignores the region entirely.

> **Azure Blob Storage** (`az://` / `https://<account>.blob.core.windows.net`)
> is **not yet supported** — it is a separate follow-up.

## Environment variables

### HTTP backend (`hfile_libcurl.c` parity)

| Variable             | Effect                                                              |
| -------------------- | ------------------------------------------------------------------ |
| `HTS_RETRY_MAX`      | Max retries on 5xx/429/transient network errors (default `3`).     |
| `HTS_RETRY_DELAY`    | Initial backoff in milliseconds (default `500`).                   |
| `HTS_RETRY_MAX_DELAY`| Backoff cap in milliseconds (default `60000`).                     |
| `CURL_CA_BUNDLE`     | PEM file loaded as the TLS root CA set.                            |
| `HTS_AUTH_LOCATION`  | File whose trimmed contents are sent as the `Authorization` header.|

Redirects are followed (the default `http.Client` behaviour, matching libcurl).
Backoff is exponential, doubling each attempt up to the cap.

### S3 backend (`hfile_s3.c` parity)

| Variable                      | Effect                                                        |
| ----------------------------- | ------------------------------------------------------------ |
| `AWS_ACCESS_KEY_ID`           | Access key (env takes precedence over the credentials file). |
| `AWS_SECRET_ACCESS_KEY`       | Secret key.                                                  |
| `AWS_SESSION_TOKEN`           | Optional session token (adds `x-amz-security-token`).        |
| `AWS_DEFAULT_REGION` / `AWS_REGION` | Region (default `us-east-1`).                          |
| `AWS_PROFILE` / `AWS_DEFAULT_PROFILE` | Profile in the credentials file (default `default`). |
| `AWS_SHARED_CREDENTIALS_FILE` | Path to the credentials file (default `~/.aws/credentials`). |
| `HTS_S3_HOST`                 | Overrides the host; forces path-style (custom endpoints).   |
| `HTS_S3_ADDRESS_STYLE`        | `path` or `virtual` to force addressing (default auto).     |
| `HTS_S3_V2`                   | If set, sign with AWS Signature Version 2 instead of V4.     |

Credential precedence: environment first, then the shared credentials file
(simple INI: `aws_access_key_id`, `aws_secret_access_key`, `aws_session_token`,
`region`). If no credentials resolve, requests are sent unsigned, which works
for public buckets.

**SigV4 (default).** Every ranged GET is signed afresh with AWS Signature
Version 4 over the fixed header set
`host;x-amz-content-sha256;x-amz-date[;x-amz-security-token]`, using
`UNSIGNED-PAYLOAD` as the content hash (no body buffering required).

**SigV2 (`HTS_S3_V2`).** The string-to-sign is
`Verb\n\n\n<Date>\n<CanonicalizedAmzHeaders><CanonicalizedResource>` (the two
blank lines are the empty Content-MD5 and Content-Type for a ranged GET), the
signature is `base64(HMAC-SHA1(secret, StringToSign))`, and the header is
`Authorization: AWS <AccessKeyId>:<signature>`. The `CanonicalizedResource` is
always `/<bucket>/<key>` (the bucket is included even under virtual-hosted
addressing), and `x-amz-security-token` is folded into the
`CanonicalizedAmzHeaders` when a session token is present. The `Date` header
(RFC 1123 GMT) is part of the signed string and is transmitted verbatim.

### GCS backend (`hfile_gcs.c` parity)

| Variable                     | Effect                                                  |
| ---------------------------- | ------------------------------------------------------- |
| `GCS_OAUTH_TOKEN`            | Adds `Authorization: Bearer <token>`. Public buckets work without it. |
| `GCS_REQUESTER_PAYS_PROJECT` | Adds the `X-Goog-User-Project` header.                  |

## Testing

Tests are fully self-contained: every backend is exercised against an
`net/http/httptest` mock server (S3 and GCS hosts are overridden via package
variables / `HTS_S3_HOST`), so **no real network or credentials are required**.
The SigV4 signer is verified against the AWS-published `get-vanilla` Signature
Version 4 test vector (canonical-request hash + final signature), and the SigV2
signer against the canonical AWS S3 REST worked example
(`GET /johnsmith/photos/puppy.jpg` → `bWq2s1WEIj+Ydj0vQ697zp+IXMU=`). The
SigV2, GCS S3-interop (V2 and V4) and forced-address-style paths each have
httptest round-trips under a pinned clock.

```bash
go test ./pkg/htsgo/hfile/... -count=1 -cover
```
