---
source: src/linalg/norm.c
---
**Algorithm.** `builtin_norm` uses `get_tensor_dims` to classify the argument and then builds a symbolic norm expression that it evaluates. A rank-0 scalar reduces to `Abs[x]` (the `p` argument is rejected for scalars). For a rank-1 vector (or any tensor when `p` is the string `"Frobenius"`) it flattens to `Expr**` and assembles: `Norm[v, Infinity]` → `Max[Abs[v_i]]`; otherwise (default `p = 2`, or `"Frobenius"` → 2, or a user `p`) → `Power[Sum[Abs[v_i]^p], 1/p]`. Every `Abs`, `Power`, `Plus`, `Max` is created and run through the evaluator (`eval_and_free`), so exact/symbolic/complex entries produce exact symbolic norms.

**Limits.** Genuine matrix p-norms (e.g. the spectral 2-norm via the largest singular value) are not implemented — a rank-≥2 argument with a non-`"Frobenius"` `p` falls through to NULL and the call stays symbolic. Jagged arrays (`rank < 0`) also return NULL.
