# FileExistsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FileExistsQ["name"]
    gives True if the file with the specified name exists, and gives False otherwise.
In FileExistsQ["name"], name is interpreted relative to your current directory.
FileExistsQ does not search $Path.
FileExistsQ tests for files, directories, or any other filesystem objects.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_fileexistsq` calls `lstat()` on a single string-path argument and returns the symbol `True` if it succeeds (anything exists at that path), `False` otherwise. Using `lstat` rather than `stat` means a dangling symlink — itself a filesystem object — is reported as existing. The POSIX `lstat` is enabled by defining `_POSIX_C_SOURCE` before the includes for strict-C99 builds. Non-string input leaves the call unevaluated. `ATTR_PROTECTED`.

- `Protected`.
- `"name"` is interpreted relative to the current working directory. `$Path` is not searched.
- Implemented with `lstat()`, so dangling symlinks count as existing.
- Leaves the call unevaluated when given the wrong arity, a symbolic argument, or any non-string atom.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/files.c`](https://github.com/stblake/mathilda/blob/main/src/files.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Put[x^2 + 1, "/tmp/mathilda_demo.m"]
Out[1]= Null

In[2]:= FileExistsQ["/tmp/mathilda_demo.m"]
Out[2]= True

In[3]:= FileExistsQ["/tmp/does_not_exist.m"]
Out[3]= False
```

### Notes

`FileExistsQ["name"]` returns `True` if a file with the given name exists and `False` otherwise.
