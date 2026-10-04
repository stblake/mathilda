### Worked examples

```mathematica
In[1]:= OutputStream["stdout", 1]  (* a stream object stays inert — it is a handle, not data *)
```

```mathematica
In[1]:= Head[OutputStream["file.txt", 3]]  (* its head is OutputStream *)
```

```mathematica
In[1]:= MemberQ[Attributes[OutputStream], Protected]  (* a protected stream-object head *)
```

### Notes

`OutputStream["name", n]` is the object representing an open output stream: the
string names it (a file, or `"stdout"`/`"stderr"`) and the integer `n` is an
internal handle into the process-global stream registry. A *live* object is produced
by `OpenWrite` or `OpenAppend` (which need a writable file), then passed to `Write`,
`WriteString`, or `Close`; the examples above therefore show only the inert object
itself. The integer handle is an internal registry detail — constructing an
`OutputStream[...]` by hand does not open anything.
