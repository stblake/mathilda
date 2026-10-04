# PatternSequence

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`PatternSequence[p1, p2, ...] is a pattern object that matches a sequence of arguments, each in turn matching p1, p2, ....`**

**`PatternSequence[] matches an empty (zero-length) sequence of arguments.`**

<details>
<summary>Notes</summary>

x:PatternSequence\[...\] binds x to the matched sequence.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

The two patterns splice into the arg list

```mathematica
In[1]:= MatchQ[f[1, 2, 3], f[PatternSequence[1, 2], 3]]
Out[1]= True
```

Mixed with ordinary patterns

```mathematica
In[2]:= MatchQ[f[1, 2, 3, 4], f[a_, PatternSequence[b_, c_], d_]]
Out[2]= True
```

Selects args beginning 1, _

```mathematica
In[3]:= Cases[{g[1, 2], g[1, 3], g[2, 1]}, g[PatternSequence[1, _]]]
Out[3]= {g[1, 2], g[1, 3]}
```

Sequence flattens in place

```mathematica
In[4]:= f[a, PatternSequence[b, c], d] /. PatternSequence -> Sequence
Out[4]= f[a, b, c, d]
```

## Implementation notes

**Definition.** `PatternSequence[p1, ..., pm]` is a **pattern object** that matches a
sequence of arguments, each in turn matching `p1, ..., pm`. It is not a function that
evaluates — it is a construct understood by the pattern matcher in `src/match.c`. The
symbol carries the `Protected` attribute and the docstring in `src/info.c`; there is no
builtin, so `PatternSequence[...]` outside a pattern position stays inert.

**Algorithm.** When the matcher meets a `PatternSequence[q1..qm]` as an element of a pattern
argument list, it **splices** the `m` sub-patterns in place of the single list element, so
the enclosing head is matched against a list of consistent arity (`match_sequence` in
`src/match.c`); `PatternSequence[]` splices to nothing and matches a zero-length run. A
*named* `x:PatternSequence[q1..qm]` is handled separately: the matcher does not splice it
but enumerates the sequence partitions and binds `x` to the matched arguments as a
`Sequence[...]`, so a `PatternSequence` body that itself contains `__`/`___` backtracks
correctly.

**Data structures.** Pattern and subject are ordinary `Expr` trees; bindings live in a
`MatchEnv`. The splice is a local rewrite of the pattern argument vector, so no copy of the
subject is made. `/.` with `PatternSequence -> Sequence` turns a literal
`PatternSequence[...]` into a spliced argument sequence, since `Sequence` flattens into its
enclosing head.

**Complexity / limits.** Sequence matching backtracks over argument partitions, so a
pattern with several variable-length sequence patterns is worst-case combinatorial in the
number of partitions, as for `__`/`___`. The unnamed spliced form adds no overhead beyond
the arity change.

**Attributes:** `Protected`.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`PatternSequence[p1, ..., pm]` matches a run of arguments that in turn match `p1, ..., pm`.
It is a pattern object, not a function — `Protected`, with no evaluation of its own — so it
has meaning only in a pattern position, where the matcher splices the `m` sub-patterns into
the surrounding argument list. `PatternSequence[]` matches a zero-length run.

A named `x:PatternSequence[...]` binds `x` to the matched arguments as a `Sequence`. The
last example shows the connection to `Sequence`: rewriting the head to `Sequence` splices
the arguments into `f`, giving `f[a, b, c, d]`. Like `__`/`___`, variable-length sequence
patterns backtrack over argument partitions.
