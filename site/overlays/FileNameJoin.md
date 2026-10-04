### Worked examples

```mathematica
In[1]:= FileNameJoin[{"usr", "local", "bin"}]  (* join with the OS separator *)
```

```mathematica
In[1]:= FileNameJoin[{"a/b", "c"}]  (* components may themselves contain separators *)
```

```mathematica
In[1]:= FileNameJoin[{"", "etc", "hosts"}]  (* an empty leading part gives an absolute path *)
```

```mathematica
In[1]:= FileNameJoin[{"dir1", "dir2"}, OperatingSystem -> "Windows"]  (* backslash separator *)
```

### Notes

`FileNameJoin[{...}]` assembles a path from components, and `FileNameJoin["name"]`
canonicalizes a single name. It is a **pure string operation** — it never touches
the filesystem. Each component is split into segments and rejoined, so duplicate
and trailing separators collapse (`{"a//b", "c"}` becomes `"a/b/c"`), and an empty
(or separator-led) leading component yields an absolute path.

The separator defaults to the host operating system's; `OperatingSystem ->
"Windows" | "MacOSX" | "Unix"` selects it explicitly (`"Windows"` uses `\` and
preserves a leading UNC `\\server\share`). `FileNameJoin` inverts
`FileNameSplit`: `FileNameJoin[FileNameSplit[name]]`
reconstructs a canonicalized `name`.
