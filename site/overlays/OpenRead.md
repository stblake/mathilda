### Worked examples

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_or.txt"]; Write[str, 10]; Write[str, 20]; Close[str];  (* make a file *)
In[2]:= ins = OpenRead["/tmp/mathilda_or.txt"];  (* open it for reading *)
In[3]:= Read[ins, Number]  (* first object *)
In[4]:= Read[ins, Number]  (* the current point has advanced *)
In[5]:= Close[ins]
```

### Notes

`OpenRead["file"]` (or `OpenRead[File["file"]]`) opens a file for reading and
returns an `InputStream["file", n]` object. The whole file is slurped into a buffer
with a moving **current point**, so successive `Read` calls on the same
stream return successive objects, and `StreamPosition` /
`SetStreamPosition` query and move that point.

`$Failed` (with an `OpenRead::noopen` diagnostic) is returned when the file cannot
be opened. A named file handed directly to `Read`/`ReadList` is auto-opened, so an
explicit `OpenRead` is for when you want to hold the stream and advance it yourself.
Finish with `Close`, which returns the file name.
