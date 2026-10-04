---
source: src/special_functions/productlog.c
references:
  - "DLMF §4.13 — the Lambert W-function."
  - "R. M. Corless et al., On the Lambert W function, Adv. Comput. Math. 5 (1996) 329-359."
---
**Algorithm.** `builtin_productlog` evaluates the Lambert W function:
`ProductLog[z]` is the principal branch `W_0`, `ProductLog[k, z]` the `k`-th
branch (`k` an explicit machine integer). Exact special values:
`ProductLog[0] = 0`, `ProductLog[E] = 1`, `ProductLog[-1/E] = -1`,
`ProductLog[-Pi/2] = I Pi/2`, `ProductLog[±Infinity] = Infinity`, and
`ProductLog[k, 0] = -Infinity` for `k != 0`. Numeric uses a unified complex-MPFR
core (`ncpx`): a region-chosen seed — the branch-point series in
`p = sqrt(2(e z + 1))` near `z = -1/e` (branches 0 and -1), the Maclaurin seed
`z(1 - z + 3/2 z^2)` for the principal branch near 0, else the asymptotic
`L1 - L2 + L2/L1` with `L1 = log z + 2 pi i k`, `L2 = log L1` — refined by
Halley's cubically-convergent iteration (Corless 1996). A real seed keeps the
iteration exactly real, so real-valued branches return a real leaf.
`D[ProductLog[z], z] = ProductLog[z]/(z(1 + ProductLog[z]))` lives in
`calculus/deriv.c`; Series at 0, at `-1/E`, and at Infinity in
`calculus/series.c`.

**Data structures.** `Expr`; the shared `ncpx` (`mpfr_t` re/im) complex toolkit
for the Halley core. ND: real-only unary kernel `NDKU_ProductLog = { NULL,
ndk_ProductLog_r, ... }` (the principal branch, via `sf_machine_productlog`),
registered `REG_U`, so `packed_aware`. Attributes: `Listable`, `NumericFunction`,
`Protected`.

**Complexity / limits.** Halley converges cubically (bounded at 100 iterations
plus one polishing step). A real leaf is returned only on the real-valued
domains (`W_0` for `x >= -1/e`, `W_{-1}` for `-1/e <= x < 0`); otherwise a
`Complex[..]`. Only numeric inputs evaluate — exact non-special arguments stay
symbolic. `Compile[]` lowers at both scalar and rank-1 array shapes
(`Compiled -> True`).
