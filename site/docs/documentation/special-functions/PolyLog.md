# PolyLog

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PolyLog[n, z]`**

gives the polylogarithm Li\_n(z) = Sum\_{k\>=1} z^k/k^n.

**`PolyLog[n, p, z]`**

gives the Nielsen generalized polylogarithm S\_{n,p}(z) (accepted but

**`PolyLog[1, z] = -Log[1-z], PolyLog[0, z] = z/(1-z), negative integer`**

**`PolyLog[n, -1] = (2^(1-n)-1) Zeta[n] for integer n >= 2, with exact forms`**

<details>
<summary>Notes</summary>

left unevaluated). Special arguments reduce in closed form: PolyLog\[n, 0\] = 0, orders give rational functions, PolyLog\[n, 1\] = Zeta\[n\] and for PolyLog\[2, 1/2\] and PolyLog\[3, 1/2\]. Inexact real or complex arguments evaluate numerically at machine or arbitrary (MPFR) precision via a power series or the Jonquiere/zeta expansion. There is a branch cut from 1 to Infinity. Listable.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= PolyLog[3, 1/2]
Out[1]= 1/6 Log[2]^3 - 1/12 Log[2] Pi^2 + 7/8 Zeta[3]

In[2]:= PolyLog[2, 0.9]
Out[2]= 1.29971
```

### Applications (8)

```mathematica
In[3]:= PolyLog[2, 1]
Out[3]= 1/6 Pi^2

In[4]:= PolyLog[1, z]
Out[4]= -Log[1 - z]

In[5]:= PolyLog[2, 1/2]
Out[5]= -1/2 Log[2]^2 + 1/12 Pi^2

In[6]:= PolyLog[3, 1/2]
Out[6]= 1/6 Log[2]^3 - 1/12 Log[2] Pi^2 + 7/8 Zeta[3]

In[7]:= PolyLog[0, z]
Out[7]= z/(1 - z)

In[8]:= PolyLog[-2, z]
Out[8]= (z + z^2)/(1 - z)^3

In[9]:= N[PolyLog[2, 1/2], 40]
Out[9]= 0.58224052646501250590265632015968010874422

In[10]:= N[PolyLog[3, 1/2 + I/2], 30]
Out[10]= 0.48615953708556007896672148708 + 0.5700774070887689781956097575898*I
```

## Algorithm

Mathilda -- the polylogarithm PolyLog.

```text
  PolyLog[n, z]     Li_n(z) = Sum_{k>=1} z^k / k^n      (|z| < 1; analytic
                    continuation elsewhere, branch cut [1, Infinity))
  PolyLog[n, p, z]  Nielsen generalized polylogarithm S_{n,p}(z).  Accepted
                    for surface compatibility but left symbolic -- there is
                    no closed-form / numeric engine for it here.
```

Evaluation is layered so each kind of argument takes the cheapest exact or fastest numeric route, mirroring src/gamma.c and src/zeta.c:

```text
  exact special values (any z) ->  closed forms
      PolyLog[n, 0]  = 0
      PolyLog[1, z]  = -Log[1 - z]
      PolyLog[0, z]  = z/(1 - z)
      PolyLog[-m, z] = Eulerian-number rational function   (m >= 1)
      PolyLog[n, 1]  = Zeta[n]                              (integer n >= 2)
      PolyLog[n, -1] = (2^(1-n) - 1) Zeta[n]                (integer n >= 2)
      PolyLog[2, 1/2] = Pi^2/12 - Log[2]^2/2
      PolyLog[3, 1/2] = Log[2]^3/6 - Pi^2 Log[2]/12 + 7 Zeta[3]/8
  numeric (>= 1 inexact operand, all numeric):
      real s, real -1 < z < 1   -> direct real MPFR power series (fast path)
      |z| <= 1/2                -> direct power series (complex MPFR)
      1/2 < |z|, |ln z| < 2 Pi  -> Jonquiere / zeta expansion (DLMF 25.12.11/12)
      otherwise                 -> stays symbolic
  everything else -> stays symbolic (return NULL)

The zeta expansion needs zeta(s - k) and Gamma(1 - s).  For real order these
```

are evaluated directly with MPFR (mpfr_zeta / mpfr_gamma); for the rare complex-order case they are obtained by evaluating the Zeta / Gamma builtins (which carry their own arbitrary-precision complex kernels).

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Gamma](../../special-functions/Gamma/), [Zeta](../../special-functions/Zeta/)

- DLMF §25.12 — the polylogarithm and Jonquiere's function (§25.12.11/12 for the zeta expansion).
- Source: [`src/special_functions/polylog.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/polylog.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_cherry_dilog.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog.c)
- Tests: [`tests/test_cherry_dilog_exp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog_exp.c)
- Tests: [`tests/test_cherry_polylog_exp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_polylog_exp.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)

## Notes & additional examples

### Notes

`PolyLog[n, z]` is the polylogarithm Li_n(z) = Σ_{k≥1} z^k / k^n. Special
arguments reduce in closed form: `PolyLog[n, 0] = 0`, `PolyLog[1, z] = -Log[1-z]`,
`PolyLog[0, z] = z/(1-z)`, negative integer orders give rational functions, and
`PolyLog[n, 1] = Zeta[n]`, `PolyLog[n, -1] = (2^(1-n)-1) Zeta[n]` for integer
`n ≥ 2`, with exact forms for `PolyLog[2, 1/2]` and `PolyLog[3, 1/2]`. Inexact
real or complex arguments evaluate numerically at machine or MPFR precision via
a power series or the Jonquière/zeta expansion. There is a branch cut from `1`
to `Infinity`. Listable. The Nielsen generalized form `PolyLog[n, p, z]` is
accepted but left unevaluated.
