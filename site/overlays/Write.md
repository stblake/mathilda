### Worked examples

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_write.txt"];
In[2]:= Write[str, 1 + 1]; Write[str, a + b]; Close[str];  (* each Write evaluates, then adds a newline *)
In[3]:= ReadList["/tmp/mathilda_write.txt"]  (* 1 + 1 was evaluated to 2 before writing *)
```

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_write2.txt"]; Write[str, Hold[2 + 2]]; Close[str];  (* Hold writes the unevaluated form *)
In[2]:= ReadList["/tmp/mathilda_write2.txt", String]
```

### Notes

`Write[stream, e1, e2, ...]` writes each expression in re-readable input form,
followed by a single newline, and flushes — so the output round-trips through
`Read`/`ReadList`. Expression arguments are **evaluated**
first (`1 + 1` becomes `2`); wrap in `Hold[...]` to write the unevaluated form.

`stream` may be an `OutputStream`, a `"file"` string, or `File["file"]`; a named
file that is not already open is auto-opened (truncating) and left open. Writing
raw text with no quoting or added newline is `WriteString`'s job.
`Write` returns `Null`.
