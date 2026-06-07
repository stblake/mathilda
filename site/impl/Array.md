---
source: src/list.c
---
**Algorithm.** `builtin_array` (in `src/list.c`) builds an N-dimensional list of `f[i1, ..., iN]` applications. It accepts `Array[f, n]`, `Array[f, {n1,...,nN}]` (dimensions), and an optional third argument giving the index origin/range per dimension. The handler normalizes the dimension spec into a `dim_count`-length `n_array` of integer counts and a parallel `r_array` of per-dimension range specs (either a shared scalar, a per-dimension list, or `NULL` for the default origin of 1). All counts must be non-negative integers or the builtin returns `NULL` (unevaluated).

The work is done by the recursive `array_helper`, which descends one dimension per call accumulating index values in `current_args`. At each level it computes the index for slot `i`: for a `{a, b}` range spec it interpolates `a + i*(b-a)/(n-1)` (evaluated symbolically through `Plus`/`Times`/`Divide` so exact rationals survive), otherwise it produces the arithmetic sequence `r_base + i` from the origin. At the deepest level it builds `f[i1, ..., iN]` and calls `evaluate` on it. Each level wraps its children in a `List[...]`.

**Data structures.** Plain `Expr**` working arrays (`n_array`, `r_array`, `current_args`); results assembled bottom-up as nested `List` expressions. Unlike `Table`, `Array` does not bind any iteration symbol — indices are passed positionally as arguments to `f`.
