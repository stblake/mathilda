### Worked examples

```mathematica
In[1]:= MatchQ[f[1, 2, 3], f[PatternSequence[1, 2], 3]]  (* the two patterns splice into the arg list *)
```

```mathematica
In[1]:= MatchQ[f[1, 2, 3, 4], f[a_, PatternSequence[b_, c_], d_]]  (* mixed with ordinary patterns *)
```

```mathematica
In[1]:= Cases[{g[1, 2], g[1, 3], g[2, 1]}, g[PatternSequence[1, _]]]  (* selects args beginning 1, _ *)
```

```mathematica
In[1]:= f[a, PatternSequence[b, c], d] /. PatternSequence -> Sequence  (* Sequence flattens in place *)
```

### Notes

`PatternSequence[p1, ..., pm]` matches a run of arguments that in turn match `p1, ..., pm`.
It is a pattern object, not a function — `Protected`, with no evaluation of its own — so it
has meaning only in a pattern position, where the matcher splices the `m` sub-patterns into
the surrounding argument list. `PatternSequence[]` matches a zero-length run.

A named `x:PatternSequence[...]` binds `x` to the matched arguments as a `Sequence`. The
last example shows the connection to `Sequence`: rewriting the head to `Sequence` splices
the arguments into `f`, giving `f[a, b, c, d]`. Like `__`/`___`, variable-length sequence
patterns backtrack over argument partitions.
