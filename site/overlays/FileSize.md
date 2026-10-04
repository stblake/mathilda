### Worked examples

```mathematica
In[1]:= Put[2, 3, 5, 7, 11, "/tmp/mathilda_fs.txt"];  (* writes "2\n3\n5\n7\n11\n" = 11 bytes *)
In[2]:= FileSize["/tmp/mathilda_fs.txt"]
```

```mathematica
In[1]:= Put[2, 3, 5, 7, 11, "/tmp/mathilda_fs.txt"]; Head[FileSize["/tmp/mathilda_fs.txt"]]  (* a plain Integer, not a Quantity *)
```

```mathematica
In[1]:= Quiet[FileSize["/tmp/mathilda_no_such_file_xyz"]]  (* missing file -> $Failed *)
```

### Notes

`FileSize["name"]` gives the file's size in **bytes** as a plain `Integer` (not a
`Quantity`). The path is relative to the current working directory; `$Path` is not
searched. It is implemented with `stat`, so a symbolic link is followed and the
target file's size is reported.

A path that cannot be `stat`'d (usually because nothing is there) returns `$Failed`
with a `FileSize::nffil` message, which respects `Quiet[]` and is visible to
`Check[]`. A symbolic or non-string argument leaves the call unevaluated.
