# hfile

`pkg/htsgo/hfile` is a pure-Go port of htslib's `hfile` virtual-filesystem
remote backends. It provides **read** access to local files and to remote
objects addressed by HTTP(S), Amazon S3 (`s3://`) and Google Cloud Storage
(`gs://`) URLs, mirroring the behaviour of htslib's `hfile_libcurl.c`,
`hfile_s3.c` and `hfile_gcs.c`.

As a **Go-port extension beyond htslib** (which has **no native Azure
backend**), it additionally reads from **Azure Blob Storage** via `az://` URLs
(or recognised `*.blob.core.windows.net` HTTPS URLs) with SAS, Shared Key,
Azure-AD bearer and anonymous authentication — see the Azure section below.

The package depends only on the Go standard library. AWS Signature Version 4
is hand-implemented with `crypto/hmac` + `crypto/sha256`, and AWS Signature
Version 2 with `crypto/hmac` + `crypto/sha1` + `encoding/base64`, and the Azure
Blob Shared Key scheme with `crypto/hmac` + `crypto/sha256` + `encoding/base64`;
there is no AWS SDK, no Azure SDK and no OAuth library.

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

- `hfile.IsRemote(name) bool` — true for `http`/`https`/`s3`/`gs`/`az` URLs.
- `hfile.SchemeOf(name) string` — the lower-cased URL scheme, or `""` for a
  bare filesystem path.

## Supported schemes

| Scheme              | Backend | Endpoint                                                            |
| ------------------- | ------- | ------------------------------------------------------------------ |
| `http://` `https://`| HTTP    | the URL as given                                                   |
| `s3://bucket/key`   | S3      | `https://<bucket>.s3.<region>.amazonaws.com/<key>` (virtual-hosted)|
| `gs://bucket/object`| GCS     | `https://storage.googleapis.com/<bucket>/<object>` (XML API)       |
| `az://account/container/blob` | Azure | `https://<account>.blob.core.windows.net/<container>/<blob>` (Go-port extension) |
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

### Azure Blob Storage (Go-port extension — htslib has no native Azure backend)

Azure Blob Storage is **not** S3-protocol: it uses Shared Key / SAS / Azure-AD
bearer auth, and htslib ships no Azure backend. This support is therefore an
**additive Go-port extension**, not an htslib parity feature. Reads are served by
ranged GETs over `https://<account>.blob.core.windows.net/<container>/<blob>`.
Address blobs as `az://<account>/<container>/<blob>`, or pass a direct
`https://<account>.blob.core.windows.net/...` URL (recognised as Azure so the
env-driven auth applies). Authentication is chosen by environment in this
priority order: inline/SAS → Shared Key → bearer → anonymous.

| Auth path                | URL / endpoint                                                              | Signing                                  | Required env vars |
| ------------------------ | -------------------------------------------------------------------------- | ---------------------------------------- | ----------------- |
| **SAS URL** (simplest)   | `https://<account>.blob.core.windows.net/<container>/<blob>?<sas>`         | none (signature is in the URL)           | *(none — the SAS travels in the URL; flows through the plain HTTPS backend)* |
| **`az://` + SAS token**  | `az://<account>/<container>/<blob>`                                        | none (SAS appended to the query)         | `AZURE_STORAGE_SAS_TOKEN` (the `?sv=...&sig=...` token) |
| **`az://` + Shared Key** | `az://<account>/<container>/<blob>`                                        | Azure Blob **SharedKey** (HMAC-SHA256), signed **per ranged GET** | `AZURE_STORAGE_ACCOUNT`, `AZURE_STORAGE_KEY` (base64 account key) |
| **`az://` + bearer (AAD)** | `az://<account>/<container>/<blob>`                                      | `Authorization: Bearer <token>`          | `AZURE_STORAGE_TOKEN` (an Azure-AD access token) |
| **Anonymous** (public)   | `az://<account>/<container>/<blob>` or the HTTPS URL                       | unsigned (still sends `x-ms-version`)    | *(none)* |

Every request carries `x-ms-version: 2021-08-06`. The endpoint can be overridden
with `AZURE_STORAGE_BLOB_ENDPOINT` (e.g. to target Azurite or an httptest
server), analogous to S3's `HTS_S3_HOST` and GCS's base-URL override.

**Shared Key string-to-sign.** Because the canonical string includes the
request's `Range` header — which differs for every ranged GET — the SharedKey
`Authorization` is computed afresh per request via the `sign` hook. The
string-to-sign (Blob service, SharedKey) is:

```
GET\n          (VERB)
\n             Content-Encoding
\n             Content-Language
\n             Content-Length   (empty, NOT 0, for a GET)
\n             Content-MD5
\n             Content-Type
\n             Date             (empty: x-ms-date carries the timestamp)
\n             If-Modified-Since
\n             If-Match
\n             If-None-Match
\n             If-Unmodified-Since
<Range>\n      Range            (the request's exact Range header value)
x-ms-date:<RFC1123 GMT>\n      } CanonicalizedHeaders: all x-ms-* headers,
x-ms-version:2021-08-06\n      }   lowercased, sorted, "name:value\n"
/<account>/<container>/<blob>  CanonicalizedResource (no query params for a plain GET)
```

The signature is `base64(HMAC-SHA256(base64decode(AccountKey), StringToSign))`
and the header is `Authorization: SharedKey <account>:<signature>` (the account
key is base64 and is decoded before the HMAC).

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

### Azure Blob backend (Go-port extension; no htslib equivalent)

| Variable                       | Effect                                                            |
| ------------------------------ | ---------------------------------------------------------------- |
| `AZURE_STORAGE_SAS_TOKEN`      | SAS token appended to the query string of `az://` / host URLs that lack an inline SAS. |
| `AZURE_STORAGE_ACCOUNT`        | Storage account name for Shared Key signing.                     |
| `AZURE_STORAGE_KEY`            | Base64 account key for Shared Key signing (decoded before HMAC). |
| `AZURE_STORAGE_TOKEN`          | Azure-AD access token; sent as `Authorization: Bearer <token>`.  |
| `AZURE_STORAGE_BLOB_ENDPOINT`  | Overrides the scheme+host (Azurite / custom endpoint / tests).   |

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

The Azure Blob Shared Key signer is verified by a **known-answer vector**: a
fixed account/key/date/`Range`/resource is pinned to an exact string-to-sign and
`Authorization` value (using the public Azurite dev key for
`devstoreaccount1`, `Range: bytes=128-191` →
`SharedKey devstoreaccount1:rFiLV/YExacX7hHP2msLtGRi89EnNIol3Xvi0s5Vn/k=`). The
SAS URL, `az://`+Shared Key (asserting two different ranges yield two different
signatures), bearer-token and anonymous paths each have httptest round-trips
under a pinned clock and `AZURE_STORAGE_BLOB_ENDPOINT`.

```bash
go test ./pkg/htsgo/hfile/... -count=1 -cover
```
