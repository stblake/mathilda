# FileExistsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FileExistsQ["name"]`**

gives True if the file with the specified name exists, and gives False otherwise.

<details>
<summary>Notes</summary>

In FileExistsQ\["name"\], name is interpreted relative to your current directory. FileExistsQ does not search $Path. FileExistsQ tests for files, directories, or any other filesystem objects.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= Put[x^2 + 1, "/tmp/mathilda_demo.m"]
Out[1]= Null

In[2]:= FileExistsQ["/tmp/mathilda_demo.m"]
Out[2]= True

In[3]:= FileExistsQ["/tmp/does_not_exist.m"]
Out[3]= False
```

## Implementation notes

`builtin_fileexistsq` calls `lstat()` on a single string-path argument and returns the symbol `True` if it succeeds (anything exists at that path), `False` otherwise. Using `lstat` rather than `stat` means a dangling symlink — itself a filesystem object — is reported as existing. The POSIX `lstat` is enabled by defining `_POSIX_C_SOURCE` before the includes for strict-C99 builds. Non-string input leaves the call unevaluated. `ATTR_PROTECTED`.

- `Protected`.
- `"name"` is interpreted relative to the current working directory. `$Path` is not searched.
- Implemented with `lstat()`, so dangling symlinks count as existing.
- Leaves the call unevaluated when given the wrong arity, a symbolic argument, or any non-string atom.

**Attributes:** `Protected`.

## References

- Source: [`src/files.c`](https://github.com/stblake/mathilda/blob/main/src/files.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_files.c`](https://github.com/stblake/mathilda/blob/main/tests/test_files.c)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)

## Notes & additional examples

### Notes

`FileExistsQ["name"]` returns `True` if a file with the given name exists and `False` otherwise.
