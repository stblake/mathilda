### Worked examples

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_in.txt"]; WriteString[s, "x"]; Close[s]; in = OpenRead["/tmp/mathilda_in.txt"]; h = Head[in]; Close[in]; h  (* OpenRead returns an InputStream *)
```

```mathematica
In[1]:= FullForm[InputStream["name", 3]]  (* "name" is the source, 3 the registry handle *)
```

### Notes

`InputStream["name", n]` is the object representing an open input stream, as returned by
`OpenRead`. `"name"` is the source (a file name) and the integer `n` is an internal
handle into the stream registry. It is an inert handle — no value, no evaluation rule.

The handle is valid only while the stream is open; `Close` invalidates it. Reading from
an already-open `InputStream` continues from its current point and leaves it open, unlike
a named-file read (`Read["file", ...]`), which opens and closes the file itself. The
companion output object is `OutputStream`.
