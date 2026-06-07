---
source: src/files.c
---
`builtin_fileexistsq` calls `lstat()` on a single string-path argument and returns the symbol `True` if it succeeds (anything exists at that path), `False` otherwise. Using `lstat` rather than `stat` means a dangling symlink — itself a filesystem object — is reported as existing. The POSIX `lstat` is enabled by defining `_POSIX_C_SOURCE` before the includes for strict-C99 builds. Non-string input leaves the call unevaluated. `ATTR_PROTECTED`.
