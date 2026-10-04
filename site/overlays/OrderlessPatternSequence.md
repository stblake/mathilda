### Worked examples

```mathematica
In[1]:= MatchQ[f[a, b], f[OrderlessPatternSequence[b, a]]]  (* order within the sequence is free *)
```

```mathematica
In[1]:= Cases[{f[1, 2], f[2, 1], f[3, 4]}, f[OrderlessPatternSequence[1, 2]]]  (* both orders selected *)
```

```mathematica
In[1]:= MatchQ[f[1, 2, 3], f[OrderlessPatternSequence[3, 1, 2]]]  (* any permutation of the args matches *)
```

### Notes

`OrderlessPatternSequence[p1, ..., pm]` matches a run of arguments that, in any order,
together match `p1, ..., pm`. It is a pattern object, not a function — `Protected`, with no
evaluation of its own — so it has meaning only in a pattern position.

It is the permutation analogue of `PatternSequence`, which matches its sub-patterns in
positional order: here the matcher searches over orderings and accepts the first that binds
every sub-pattern, which is why both `f[1, 2]` and `f[2, 1]` are selected above. Because it
enumerates permutations, its cost grows factorially in the number of sub-patterns, so it is
intended for a small handful.
