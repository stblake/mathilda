---
source: src/patterns.c
---
**Algorithm.** `builtin_first_case` is a thin front-end over `Cases`: it builds
`Cases[expr, pattern]`, evaluates it, and returns a copy of the first element of
the resulting match list. Delegating this way means `FirstCase` inherits `Cases`'s
whole surface for free — the pattern may be a plain pattern or a transformation
rule `patt -> rhs` / `patt :> rhs`, in which case the returned value is the
transformed first match, not the matching element. A visible `NDArray` is an atom
to the matcher, so it is materialised to a list first (`patterns_delist_visible`).

If `Cases` yields no match, the three-argument form `FirstCase[expr, pattern,
default]` returns `default` and the two-argument form returns
`Missing["NotFound"]`.

**Data structures.** One synthesised `Cases[...]` call `Expr`, evaluated and
freed; the result is a copy of the match list's first argument (or of the
supplied default).

**Complexity / limits.** Because it delegates to the two-argument `Cases`, the
full match list is materialised before the first element is taken — there is no
first-match early stop (contrast `FirstPosition`, which passes `n = 1` to
`Position`). A level spec is not accepted: the third argument is the default, not
a level.
