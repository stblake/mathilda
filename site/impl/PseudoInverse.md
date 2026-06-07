---
references:
  - "A. Ben-Israel and T. N. E. Greville, *Generalized Inverses: Theory and Applications*, 2nd ed. (Springer, 2003)."
source: src/linalg/inv.c
---
**Algorithm.** `builtin_pseudoinverse` computes the Moore-Penrose pseudoinverse via an **exact full-rank decomposition**, not via SVD. `pseudoinverse_exact` row-reduces `A` (`mat_rref`), uses `find_pivots` to recover the rank `r` and the pivot columns, and forms a rank factorisation `A = B·C`: `B` is the `m × r` matrix of `A`'s pivot columns (`extract_columns`), `C` is the `r × n` matrix of the non-zero RREF rows (`extract_rows`). The pseudoinverse is then assembled from the closed form
`A⁺ = Cᴴ (C Cᴴ)⁻¹ (Bᴴ B)⁻¹ Bᴴ`,
using `hermitian_transpose`, `mat_mult`, and `mat_invert` on the small `r × r` Gram matrices, with the product finally `expand_matrix`-ed. The zero matrix (`r == 0`) returns the `n × m` zero matrix; an invertible square `A` collapses to the ordinary inverse.

**Data structures.** Standard `Expr*` `List`-of-`List` matrices throughout; the `r × r` Gram matrices `C Cᴴ` and `Bᴴ B` are full-rank by construction so `mat_invert` always succeeds. Inexact (`Real`/`MPFR`) input goes through the `common_rationalize_input` → exact-pipeline → `common_numericalize_result` round-trip at the input precision, giving a well-defined rank during row reduction.

**Limits.** A `Tolerance` option is parsed but currently has no effect on the exact pipeline. Non-rectangular input emits `PseudoInverse::matrix`.
