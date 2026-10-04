### Worked examples

```mathematica
In[1]:= Put[100, 200, 300, "/tmp/mathilda_ssp.txt"];
In[2]:= ins = OpenRead["/tmp/mathilda_ssp.txt"];
In[3]:= Read[ins, Number]  (* 100 *)
In[4]:= SetStreamPosition[ins, 0]  (* rewind to the start; returns the new position *)
In[5]:= Read[ins, Number]  (* 100 again *)
In[6]:= SetStreamPosition[ins, Infinity]  (* jump to end of file: the byte length *)
In[7]:= Close[ins]
```

### Notes

`SetStreamPosition[stream, n]` moves the current point to byte offset `n` and
returns the new position; `SetStreamPosition[stream, Infinity]` moves to the end of
the stream (so it reports the byte length). A negative `n` clamps to `0` and an `n`
past the end clamps to the length, so the position stays within the file.

Combined with `StreamPosition` this gives seek-style random
access: rewind and re-read, or skip ahead. `$Failed` is returned for a stream that
is not open.
