---
source: src/special_functions/polylog.c
references:
  - "DLMF §25.12 — the polylogarithm and Jonquiere's function (§25.12.11/12 for the zeta expansion)."
---
**Algorithm.** `builtin_polylog` evaluates `PolyLog[n, z] = Li_n(z)`. Exact
closed forms: `Li_n(0) = 0`; `Li_1(z) = -Log[1-z]`; `Li_0(z) = z/(1-z)`;
`Li_{-m}(z)` for `m >= 1` the Eulerian-number rational function (cap
`POLYLOG_NEGINT_CAP = 400`); for integer `n >= 2`, `Li_n(1) = Zeta[n]`,
`Li_n(-1) = (2^(1-n)-1) Zeta[n]`, and the special values `Li_2(1/2)`,
`Li_3(1/2)`. Numeric (at least one inexact operand, all numeric): a real order
with real `-1 < z < 1` takes a direct real-MPFR power series fast path;
`|z| <= 1/2` the direct complex series; `1/2 < |z|` with `|ln z| < 2 pi` the
Jonquiere/zeta expansion (`zeta(s-k)` and `Gamma(1-s)`, with `zeta` reflected
through the functional equation in the left half-plane); otherwise symbolic. The
branch cut `[1, Infinity)` is taken continuous from below (a negative-zero
imaginary part). `PolyLog[n, p, z]` (Nielsen) is accepted but left symbolic.

**Data structures.** `Expr`; a local `pcx` (`mpfr_t` re/im) toolkit; real paths
use `mpfr_zeta`/`mpfr_gamma`, complex-order paths reuse the `Zeta`/`Gamma`
builtins. ND: binary kernel `NDK_BIN2(PolyLog, sf_machine_polylog)` (order and
argument), registered `REG_B`, so `packed_aware`. Attributes: `Listable`,
`NumericFunction`, `Protected`.

**Complexity / limits.** Each numeric regime is gated by `|z|`/`|ln z|`; outside
them the call stays symbolic. Integer orders take the exact/closed-form or
fast-path routes. `Compile[]` lowers at both scalar and rank-1 array shapes
(`Compiled -> True`).
