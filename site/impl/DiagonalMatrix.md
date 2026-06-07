---
source: src/linalg/construct.c
---
`builtin_diagonalmatrix` builds a matrix placing the entries of the given `List` on the `k`-th diagonal. With one argument the entries go on the main diagonal (`k = 0`); a second integer argument `k` selects a super-/sub-diagonal (`j - i == k`), sizing the matrix to `(s + |k|) × (s + |k|)` where `s` is the list length; an optional third argument fixes the output dimensions as `n` or `{m, n}`. Off-diagonal cells are `Integer` `0`; diagonal cells are deep-copied verbatim from the input list, so symbolic, exact, and inexact entries flow through unchanged. Malformed `k`/dimension specs return the call unevaluated. The result is a `List` of `List`s.
