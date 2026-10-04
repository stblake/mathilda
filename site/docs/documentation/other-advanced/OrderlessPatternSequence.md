# OrderlessPatternSequence

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`OrderlessPatternSequence[p1, p2, ...] is a pattern object that matches a sequence of arguments, in any order, that together match p1, p2, ....`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

Order within the sequence is free

```mathematica
In[1]:= MatchQ[f[a, b], f[OrderlessPatternSequence[b, a]]]
Out[1]= True
```

Both orders selected

```mathematica
In[2]:= Cases[{f[1, 2], f[2, 1], f[3, 4]}, f[OrderlessPatternSequence[1, 2]]]
Out[2]= {f[1, 2], f[2, 1]}
```

Any permutation of the args matches

```mathematica
In[3]:= MatchQ[f[1, 2, 3], f[OrderlessPatternSequence[3, 1, 2]]]
Out[3]= True
```

## Implementation notes

**Definition.** `OrderlessPatternSequence[p1, ..., pm]` is a **pattern object** that matches
a sequence of arguments which, taken **in any order**, together match `p1, ..., pm`. It is
not a function that evaluates — it is a construct understood by the pattern matcher in
`src/match.c`. The symbol carries the `Protected` attribute and the docstring in
`src/info.c`; there is no builtin, so it is inert outside a pattern position.

**Algorithm.** When the matcher meets an `OrderlessPatternSequence[q0..q_{nq-1}]` in a
pattern argument list (optionally inside a named `x:` wrapper), it tries to match the `nq`
sub-patterns against a chosen subset of the subject's arguments **under every permutation**,
accepting the first assignment that binds all sub-patterns consistently (the
`orderless_pattern_sequence` handling in `src/match.c`). This is the permutation analogue of
`PatternSequence`, which matches the sub-patterns in positional order; it is what lets a
pattern pick out a required combination of arguments regardless of where they appear, the
way `Orderless` attributes do for whole heads.

**Data structures.** Pattern and subject are ordinary `Expr` trees; bindings live in a
`MatchEnv`. Matching enumerates permutations of the sub-patterns against candidate argument
positions, reusing `env_new`/`match`/`env_free` and the shared sequence-backtracking
machinery.

**Complexity / limits.** Because it searches over orderings, the cost grows factorially in
the number of sub-patterns `nq`, so it is meant for a handful of sub-patterns, not a long
list. Within those sub-patterns, ordinary and variable-length patterns (`_`, `__`, `___`)
are allowed and backtrack as usual.

**Attributes:** `Protected`.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`OrderlessPatternSequence[p1, ..., pm]` matches a run of arguments that, in any order,
together match `p1, ..., pm`. It is a pattern object, not a function — `Protected`, with no
evaluation of its own — so it has meaning only in a pattern position.

It is the permutation analogue of `PatternSequence`, which matches its sub-patterns in
positional order: here the matcher searches over orderings and accepts the first that binds
every sub-pattern, which is why both `f[1, 2]` and `f[2, 1]` are selected above. Because it
enumerates permutations, its cost grows factorially in the number of sub-patterns, so it is
intended for a small handful.
