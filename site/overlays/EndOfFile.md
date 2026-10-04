### Worked examples

```mathematica
In[1]:= Head[EndOfFile]  (* a bare sentinel symbol, not a function *)
```

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_eof.txt"]; WriteString[str, "only\n"]; Close[str];
In[2]:= ins = OpenRead["/tmp/mathilda_eof.txt"]; Read[ins, Word]; eof = Read[ins, Word]; Close[ins]; eof
```

```mathematica
In[1]:= EndOfFile === EndOfFile  (* a plain marker value, compared by identity *)
```

### Notes

`EndOfFile` is the sentinel symbol `Read` returns at the end of a stream. It is an inert
marker — no builtin, no DownValues — *produced* by the reader (`src/io/read.c`), not
consumed: the second `Read[ins, Word]` above has already consumed the only word, so it
hands back `EndOfFile`.

`ReadList` uses the same sentinel to pad the trailing slots of a type sequence that
end-of-file truncated, and stops its loop when it appears. Incremental readers test for it
by identity, as in `While[(x = Read[s, Word]) =!= EndOfFile, ...]`. It carries no numeric or
structural meaning of its own.
