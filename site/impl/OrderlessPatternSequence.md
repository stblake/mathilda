---
source: src/match.c
---
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
