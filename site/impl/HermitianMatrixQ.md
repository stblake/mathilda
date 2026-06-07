---
source: src/list.c
---
`builtin_hermitian_matrix_q` tests whether a matrix equals its conjugate transpose, i.e. `m[i,j] == Conjugate[m[j,i]]`. After validating that the argument is a non-empty square `List` of `List`s with no deeper nesting (returning `False` otherwise), it walks the upper triangle including the diagonal (the pair test is symmetric under `(i,j)↔(j,i)`) and checks each pair with one of three predicates: the default structural test (`hermitian_pair_structural`, exact for symbolic/exact-numeric entries), a user `SameTest -> f`, or `Tolerance -> t` (accepting pairs with `Abs[a - Conjugate[b]] <= t`). `SameTest`/`Tolerance` of `Automatic` fall through to the structural test; any unrecognised option leaves the call unevaluated. Returns `True`/`False`.
