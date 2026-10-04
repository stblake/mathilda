---
source: src/match.c
---
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
