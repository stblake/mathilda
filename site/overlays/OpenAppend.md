### Worked examples

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_oa.txt"]; Write[str, 1]; Close[str];  (* file now holds 1 *)
In[2]:= app = OpenAppend["/tmp/mathilda_oa.txt"];  (* open at the end, keeping the 1 *)
In[3]:= Write[app, 2]; Write[app, 3]; Close[app];
In[4]:= ReadList["/tmp/mathilda_oa.txt"]  (* all three lines survive *)
```

### Notes

`OpenAppend["file"]` opens a file for writing **at its end** and returns an
`OutputStream["file", n]` object. Unlike `OpenWrite`, which
truncates, it preserves the existing contents and adds after them — the round trip
above ends with `{1, 2, 3}`. If the file does not exist it is created empty, so
`OpenAppend` is safe on a first run.

`$Failed` (with an `OpenAppend::noopen` diagnostic) is returned when the file
cannot be opened. Finish with `Close`, which returns the file name.
