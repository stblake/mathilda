### Worked examples

```mathematica
In[1]:= Put[100, 200, 300, "/tmp/mathilda_sp.txt"];  (* writes "100\n200\n300\n" *)
In[2]:= ins = OpenRead["/tmp/mathilda_sp.txt"];
In[3]:= Read[ins, Number]  (* 100 *)
In[4]:= StreamPosition[ins]  (* byte offset after "100\n" is 4 *)
In[5]:= Close[ins]
```

### Notes

`StreamPosition[stream]` returns the stream's current point as an integer **byte
offset** from the start of the file. For an input stream it is the position in the
slurped buffer; for an output stream it is `ftell` of the underlying file. It pairs
with `SetStreamPosition`, which moves the point, so the two
together give random access within a file.

Reading advances the offset, which is why the value above is `4` after one
`Read` of `"100\n"`. `$Failed` is returned for a stream that is not open.
