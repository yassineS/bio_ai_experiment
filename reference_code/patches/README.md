# Upstream reference patches

Small, surgical patches applied to the vendored upstream sources **at build
time** so the reference binaries run on a modern toolchain. They are NOT
committed into the submodules (the submodule pointers stay pristine); the
documented build command (`pipeline/internal/upstream` `hint`) applies the
patch before configuring/compiling.

## `vcftools-tmpfile-vla-off-by-one.patch`

vcftools allocates the `mkstemp` template into a stack VLA sized **without** the
terminating NUL:

```cpp
string new_tmp = params.temp_dir + "/vcftools.XXXXXX";
char tmpname[new_tmp.size()];          // one byte too small
strcpy(tmpname, new_tmp.c_str());      // writes size()+1 bytes
```

`strcpy` writes `size()+1` bytes (including the NUL), one byte past the buffer.
On a modern glibc built with `_FORTIFY_SOURCE` (the distro default at `-O2`),
`__strcpy_chk` detects the overrun and `abort()`s with *"buffer overflow
detected"* — which is why `--geno-r2`, `--hap-r2`, and `--012` crash before
writing any output. It is a genuine (long-standing, benign-on-older-builds)
vcftools bug, not a divergence in our port.

The patch changes every `char tmpname[new_tmp.size()]` /
`char tmpname2[new_tmp.size()]` to `... [new_tmp.size()+1]` (19 sites across
`variant_file_output.cpp` and `variant_file_format_convert.cpp`). With the fix
the binary builds with `_FORTIFY_SOURCE` on and produces correct output.
