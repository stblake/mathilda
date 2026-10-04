### Worked examples

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_ow.txt"];  (* open for writing (object bound to str) *)
In[2]:= Write[str, x + y]; Write[str, x^2]; Close[str];  (* write two lines, then close *)
In[3]:= ReadList["/tmp/mathilda_ow.txt"]  (* read them back *)
```

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_ow2.txt"]; Write[str, 111]; Close[str];
In[2]:= str = OpenWrite["/tmp/mathilda_ow2.txt"]; Write[str, 222]; Close[str];  (* OpenWrite truncates: 111 is gone *)
In[3]:= ReadList["/tmp/mathilda_ow2.txt"]
```

### Notes

`OpenWrite["file"]` opens a file for writing and returns an
`OutputStream["file", n]` object; the integer `n` is an internal handle into the
stream registry. It **truncates** any existing file, so the second block above ends
with just `222`. To keep prior contents and add to the end, use
`OpenAppend` instead.

A named file that is not already open is also auto-opened (truncating) by
`Write`/`WriteString`, so an explicit `OpenWrite` is
only needed to hold the stream across several writes or to interleave with
`StreamPosition`. Finish with `Close`, which returns the file name; any
stream still open is closed at program exit, so nothing leaks.
