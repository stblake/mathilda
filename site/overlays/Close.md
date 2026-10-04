### Worked examples

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_close.txt"];  (* open a stream *)
In[2]:= Close[str]  (* Close returns the file name, not Null *)
```

```mathematica
In[1]:= OpenWrite["/tmp/mathilda_close2.txt"];  (* open, discarding the object *)
In[2]:= Close["/tmp/mathilda_close2.txt"]  (* close by file name instead of by object *)
```

### Notes

`Close[stream]` closes an open `InputStream`/`OutputStream` and **returns its file
name** (a string), which is what makes `Close` convenient as the last line of a
round trip. `Close["file"]` / `Close[File["file"]]` closes the stream open for that
name — the most recently opened, if several share it.

Closing frees the registry slot and, for an output stream, flushes and closes the
underlying file. Closing a stream that is not open returns `$Failed` with a
`Close::stream` message. Streams left open are closed automatically at program
exit, so a forgotten `Close` does not leak — but an output file is only guaranteed
flushed to disk once it is closed.
