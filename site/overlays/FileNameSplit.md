### Worked examples

```mathematica
In[1]:= FileNameSplit["/home/user/data.csv"]  (* leading "" marks an absolute path *)
```

```mathematica
In[1]:= FileNameSplit["a//b/c/"]  (* duplicate and trailing separators are dropped *)
```

```mathematica
In[1]:= FileNameJoin[FileNameSplit["/usr/local/bin"]]  (* Join inverts Split *)
```

```mathematica
In[1]:= FileNameSplit["C:\\path\\file.txt", OperatingSystem -> "Windows"]  (* a drive is an ordinary first part *)
```

### Notes

`FileNameSplit["name"]` is the structural inverse of
`FileNameJoin`: it returns the list of path components. It is a
**pure string operation** and never touches the filesystem. A leading separator
makes the path absolute and yields a leading `""` part; trailing and duplicate
separators are dropped.

The separator defaults to the host operating system's; `OperatingSystem ->
"Windows" | "MacOSX" | "Unix"` selects it. On `"Windows"` a leading UNC
`\\server\share` prefix is kept as a single part and a drive like `C:` falls out
as an ordinary first part. A non-string argument or an unknown OS leaves the call
unevaluated.
