# Series

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Series[f, {x, x0, n}]`**

generates a power-series expansion of f about x = x0 up to order (x - x0)^n.

**`Series[f, x -> x0]`**

generates the leading term of a power-series expansion of f about x = x0.

**`Series[f, {x, x0, nx}, {y, y0, ny}, ...]`**

iteratively expands f, first in x, then in y, etc.

**`Series[f, {x, Infinity, n}] expands around x = Infinity by substituting x -> 1/u.`**

<details>
<summary>Notes</summary>

Series handles Taylor, Laurent (negative powers), and Puiseux (fractional powers) expansions, as well as logarithmic and symbolic-exponent cases such as x^x and (1+x)^n. The Assumptions -\> assm option (also read from an ambient Assuming\[...\] scope or $Assumptions) uses the sign/reality/domain of parameters and the expansion variable to simplify coefficients (Sqrt\[a^2\] -\> a, Abs\[a\] -\> a, Log\[a^p\] -\> p Log\[a\] for a \> 0), pick the Log branch of the integral family at x = 0, and expand non-analytic heads (Abs\[x\], Sign\[x\], UnitStep\[x\], Conjugate\[x\]). The result of Series is a SeriesData object; use Normal to convert it back to an ordinary expression by dropping the O-term. Series is Protected and HoldAll so the expansion variable is not evaluated.

</details>

## Examples (20)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= Series[Exp[x], {x, 0, 10}]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + 1/120 x^5 + 1/720 x^6 + 1/5040 x^7 + 1/40320 x^8 + 1/362880 x^9 + 1/3628800 x^10 + O[x]^11

In[2]:= Series[f[x], {x, a, 3}]
Out[2]= f[a] + Derivative[1][f][a] (x - a) + 1/2 Derivative[2][f][a] (x - a)^2 + 1/6 Derivative[3][f][a] (x - a)^3 + O[x - a]^4

In[3]:= Series[Cos[x]/x, {x, 0, 10}]
Out[3]= 1/x - 1/2 x + 1/24 x^3 - 1/720 x^5 + 1/40320 x^7 - 1/3628800 x^9 + O[x]^11

In[4]:= Series[Sqrt[Sin[x]], {x, 0, 10}]
Out[4]= Sqrt[x] - 1/12 x^(5/2) + 1/1440 x^(9/2) - 1/24192 x^(13/2) - 67/29030400 x^(17/2) + O[x]^(21/2)

In[5]:= Series[x^x, {x, 0, 4}]
Out[5]= 1 + Log[x] x + 1/2 Log[x]^2 x^2 + 1/6 Log[x]^3 x^3 + 1/24 Log[x]^4 x^4 + O[x]^5

In[6]:= Series[(1 + x)^n, {x, 0, 4}]
Out[6]= 1 + n x + 1/2 n (-1 + n) x^2 + 1/6 n (-2 + n) (-1 + n) x^3 + 1/24 n (-3 + n) (-2 + n) (-1 + n) x^4 + O[x]^5

In[7]:= Series[Sin[1/x], {x, Infinity, 10}]
Out[7]= 1/x - 1/6 (1/x)^3 + 1/120 (1/x)^5 - 1/5040 (1/x)^7 + 1/362880 (1/x)^9 + O[1/x]^11

In[8]:= Series[Sin[x + y], {x, 0, 3}, {y, 0, 3}]
Out[8]= y - 1/6 y^3 + O[y]^4 + (1 - 1/2 y^2 + O[y]^4) x + (-1/2 y + 1/12 y^3 + O[y]^4) x^2 + (-1/6 + 1/12 y^2 + O[y]^4) x^3 + O[x]^4
```

### Worked examples (4)

```mathematica
In[9]:= Series[Erf[x], {x, Infinity, 2}]
Out[9]= 1 + E^(-x^2) ((-1/Sqrt[Pi])/x + O[1/x]^3)

In[10]:= Series[LogGamma[x], {x, Infinity, 2}]
Out[10]= 1/2 Log[2 Pi] - x + Log[x] (-1/2 + x) + 1/12/x + O[1/x]^3

In[11]:= Series[PolyGamma[0, x], {x, Infinity, 4}]
Out[11]= Log[x] + -1/2/x - 1/12 (1/x)^2 + 1/120 (1/x)^4 + O[1/x]^5

In[12]:= Series[PolyGamma[1, x], {x, Infinity, 5}]
Out[12]= 1/x + 1/2 (1/x)^2 + 1/6 (1/x)^3 - 1/30 (1/x)^5 + O[1/x]^6
```

### Applications (8)

```mathematica
In[13]:= Series[Sin[x], {x, 0, 5}]
Out[13]= x - 1/6 x^3 + 1/120 x^5 + O[x]^6

In[14]:= Series[1/(1 - x), {x, 0, 4}]
Out[14]= 1 + x + x^2 + x^3 + x^4 + O[x]^5

In[15]:= Series[Log[1 + x], {x, 0, 4}]
Out[15]= x - 1/2 x^2 + 1/3 x^3 - 1/4 x^4 + O[x]^5

In[16]:= Normal[Series[Exp[x], {x, 0, 3}]]
Out[16]= 1 + x + 1/2 x^2 + 1/6 x^3

In[17]:= Series[Tan[x], {x, 0, 7}]
Out[17]= x + 1/3 x^3 + 2/15 x^5 + 17/315 x^7 + O[x]^8

In[18]:= Series[1/(Exp[x] - 1), {x, 0, 4}]
Out[18]= 1/x - 1/2 + 1/12 x - 1/720 x^3 + O[x]^5

In[19]:= Series[x^x, {x, 0, 3}]
Out[19]= 1 + Log[x] x + 1/2 Log[x]^2 x^2 + 1/6 Log[x]^3 x^3 + O[x]^4

In[20]:= Series[(1 + 1/x)^x, {x, Infinity, 2}]
Out[20]= E + (-1/2 E)/x + 11/24 E (1/x)^2 + O[1/x]^3
```

## Options & behaviour

**Coefficient arithmetic** automatically promotes to BigInt-backed `Rational` when 64-bit numerators or denominators would overflow, so previously-overflowing Laurent/Puiseux cases like `Series[1/Sin[x]^10, {x, 0, 2}]` and `Series[Sqrt[Log[1 + x]], {x, 0, 12}]` now produce exact coefficients (at the cost of slower evaluation for large orders).

Inverse trigonometric and inverse hyperbolic heads (`ArcSin`, `ArcCos`, `ArcTan`, `ArcCot`, `ArcSinh`, `ArcCosh`, `ArcTanh`, `ArcCoth`) are handled by direct series kernels at `u = 0` rather than by naive repeated differentiation, which would blow up expression size exponentially for higher orders. `ArcCosh` uses the principal-branch identity `ArcCosh[u] = I*ArcCos[u]`, so its expansion at `x = 0` has the expected `I*Pi/2` constant term and imaginary coefficients.

Forward reciprocal heads (`Sec`, `Csc`, `Cot`, `Sech`, `Csch`, `Coth`) are rewritten as `1/Cos[x]`, `1/Sin[x]`, `Cos[x]/Sin[x]`, etc., before expansion. Inverse reciprocal heads (`ArcSec`, `ArcCsc`, `ArcSech`, `ArcCsch`) are rewritten via the identities `ArcSec[z] = ArcCos[1/z]`, `ArcCsc[z] = ArcSin[1/z]`, `ArcSech[z] = ArcCosh[1/z]`, `ArcCsch[z] = ArcSinh[1/z]`, so that a blowing-up inner series (e.g. `z = 1/x`) collapses to a convergent kernel case rather than triggering spurious `Power::infy` warnings.

Expansions where the inner series diverges at the expansion point (e.g. `Series[f[1/x], {x, 0, n}]`) are handled via dedicated at-infinity identities:
- `ArcCoth[1/u] = ArcTanh[u]`, `ArcCot[1/u] = ArcTan[u]` (handled at the series level via inner-series inversion).
- `ArcTan[u] = Pi/2 - ArcTan[1/u]` when the argument blows up (the `+Infinity` real direction). Unlike `ArcCot`/`ArcCsc`/`ArcSec`/`ArcCoth` — whose reciprocal-argument value at the singular point is finite, so the generic `x -> 1/u` Taylor path already handles them — `ArcTan[1/u]` probes `ArcTan[ComplexInfinity]` (Indeterminate) at `u = 0`, so `ArcTan` needs this explicit kernel identity to compose inside `Plus`/`Times` (e.g. `Series[x (Pi/2 - ArcTan[x]), {x, Infinity, n}]`). The bare `Series[ArcTan[x], {x, Infinity, n}]` = `Pi/2 - 1/x + 1/(3 x^3) - 1/(5 x^5) + O[1/x]^(n+1)` is also emitted directly.
- `ArcTanh[1/u] = I*Pi/2 + ArcTanh[u]` (principal branch).
- `ArcSinh[1/v] = -Log[v] + Log[1 + Sqrt[1 + v^2]]` and `ArcCosh[1/v] = -Log[v] + Log[1 + Sqrt[1 - v^2]]` (handled by rewriting at the expression level; the symbolic `-Log[x]` term rides the existing `Log[x]` symbolic-coefficient path).

**Asymptotic expansions at Infinity for special functions with essential singularities**: Some functions have no Laurent series at `Infinity` because their leading behaviour is an essential singularity (a factor like `E^x`). For these the generic `x -> 1/u` substitution would hand a pole to naive Taylor, so they are emitted from dedicated asymptotic identities with the essential factor kept symbolic:
- `Series[ExpIntegralEi[x], {x, Infinity, n}]` returns `E^x (1/x + 1/x^2 + 2/x^3 + ... + O[1/x]^(n+1))`, i.e. `Times[Power[E, x], SeriesData[Power[x, -1], 0, {0!, 1!, ..., (n-1)!}, 1, n+1, 1]]` (DLMF 6.12.2: `Ei(x) ~ E^x Σ_{k≥0} k!/x^(k+1)`). The `E^x` factor rides the expression-level `Times` exactly as a symbolic `x^alpha` prefactor does.
- `Series[Erf[x], {x, Infinity, n}]`, `Series[Erfc[x], ...]`, `Series[Erfi[x], ...]` return the error-function asymptotic expansions (DLMF 7.12.1), each an `Exp[±x^2]` essential-singularity prefactor times a Laurent series in `1/x` with only odd powers populated.  (the leading `1` is the limit); `Erfc = 1 - Erf` uses the negated multiplier without the constant; `Erfi` has all-positive coefficients and a growing `Exp[+x^2]` prefactor. The general coefficient is `a_k = (2k-1)!!/(2^k Sqrt[Pi])` at `x^-(2k+1)`, with alternating signs for `Erf`/`Erfc`.
- `Series[LogGamma[x], {x, Infinity, n}]` returns the Stirling expansion (DLMF 5.11.1): an *additive* growth head `(x - 1/2) Log[x] - x + Log[2 Pi]/2` (kept symbolic, as the `E^x` prefactor is for `ExpIntegralEi`) plus a Bernoulli Laurent tail `Σ_{k≥1} B_{2k}/(2k(2k-1)) x^-(2k-1)` in `1/x`, i.e. `1/(12x) - 1/(360 x^3) + 1/(1260 x^5) - ...`. . (`Gamma[x] = Exp[LogGamma[x]]` still has no series at infinity — it diverges — so `Series[Gamma[x], {x, Infinity, n}]` stays unevaluated.)
- `Series[PolyGamma[m, x], {x, Infinity, n}]` returns the polygamma Stirling expansion (DLMF 5.11.2). The digamma `m = 0` mirrors `LogGamma`: an additive `Log[x]` growth head plus a Laurent tail `-1/(2x) - Σ_{k≥1} B_{2k}/(2k) x^-(2k)`, so . For `m ≥ 1` it is a pure Laurent series in `1/x` with leading power `x^-m` (decays to 0): `PolyGamma[m, x] ~ (-1)^(m-1) [ (m-1)!/x^m + m!/(2 x^(m+1)) + Σ_{k≥1} B_{2k}(2k+m-1)!/(2k)! x^-(2k+m) ]`. E.g. . `PolyGamma[x]` (unindexed) normalises to `PolyGamma[0, x]`. A requested order below the leading power `m` yields a pure `O[1/x]^(n+1)`.
- `Series[Zeta[x], {x, Infinity, n}]` returns the truncated Dirichlet head `1 + 2^-x + 3^-x + ... + (n+1)^-x`. Because `Zeta[x] = Σ_{k≥1} k^-x` and each `k^-x = Exp[-x Log[k]]` is exponentially smaller than the previous as `x → +∞`, the natural asymptotic scale here is *exponential* (`2^-x`), not a power of `1/x` — so this is a plain `Plus` of powers, not a `SeriesData` with an `O`-term. `n` counts the correction terms kept (bases `2..n+1`); the leading `1` is always present. This is exactly the exp-log form the Gruntz mrv engine consumes to resolve limits like `(Zeta[x] - 1) 2^x → 1` and `Log[Zeta[x] - 1]/x → -Log[2]`.

**Logarithmic series at `x = 0` for special functions**: Some functions have a finite-radius series at the origin built around a single `Log[x]` branch term, with the log baked into the `x^0` coefficient. Naive Taylor-via-`D` cannot reach these (`f(0)` is infinite or the derivatives have poles), so they are emitted from dedicated identities:
- `Series[ExpIntegralEi[x], {x, 0, n}]` returns `EulerGamma + Log[x] + x + x^2/4 + x^3/18 + ... + O[x]^(n+1)`, i.e. `SeriesData[x, 0, {EulerGamma + Log[x], 1, 1/4, ..., 1/(n*n!)}, 0, n+1, 1]` (DLMF 6.6.2: `Ei(x) = EulerGamma + Log[x] + Σ_{k≥1} x^k/(k k!)`). The `EulerGamma + Log[x]` branch term occupies the `x^0` slot. With `Assumptions -> x < 0` the branch term becomes `EulerGamma + Log[-x]`.

**Generalized series in `Log[x]` at `x = 0`**: `LogIntegral[x] = Ei(Log[x])` has *no* ordinary Taylor (or Laurent/Puiseux) series at `x = 0`: as `x -> 0+`, `L = Log[x] -> -Infinity` drives it into Ei's asymptotic regime, so the result is a generalized series whose `x^1` coefficient is itself a series in `1/Log[x]`.
- `Series[LogIntegral[x], {x, 0, n}]` returns `x Σ_{k=0}^{2n+1} k!/Log[x]^(k+1) + O[x]^(n+1)`. For `n = 2` this is `((120 + 24 Log[x] + 6 Log[x]^2 + 2 Log[x]^3 + Log[x]^4 + Log[x]^5) x)/Log[x]^6 + O[x]^3`, matching Mathematica's `Assumptions -> x > 0` output. Since the only `x`-dependence is the prefactor `E^L = x`, every term carries exactly `x^1`; the coefficient slot is emitted as `Together[Σ k!/Log[x]^(k+1)]`. Mathematica's no-assumptions form additionally wraps this in a `Floor[Arg[...]]` branch discriminator (to track `Log[1/x]` vs `Log[x]` across the cut); Mathilda emits the principal `x > 0` form by default, or — with `Assumptions -> x < 0` — the `x < 0` form, where every `Log[x]` becomes `Log[-x]` and an additive `I Pi` (the `x^0` term, from `Log[x] = Log[-x] + I Pi`) leads the series: `I Pi + ((120 + 24 Log[-x] + ... + Log[-x]^5) x)/Log[-x]^6 + O[x]^3`.

### Internal padding for symbolic expansion points

The engine computes series at a padded internal order (user order + 12 by default) so that intermediate Laurent/Puiseux operations don't lose accuracy. When the expansion point `x0` is not a literal number, padding is capped at 2 — at a symbolic point the series coefficients are themselves symbolic expressions (e.g. `Cosh[a]`, `Sinh[a]`), and the `O(N^2)` convolution inside `so_inv`/`so_div` would otherwise spin indefinitely on exponentially growing expression trees. This makes cases like `Series[Coth[x], {x, a, 1}]`, `Series[Tanh[x], {x, a, 1}]`, `Series[Sec[x], {x, a, 1}]`, and `Series[1/Cosh[x], {x, a, 1}]` terminate in milliseconds.

### Constant inputs

If `f` is free of the expansion variable (e.g. `Series[0, {x, 0, 4}]`, `Series[Sin[y], {x, 0, 4}]`, `Series[a + b^2, {x, 0, 3}]`), `Series` returns `f` verbatim instead of wrapping it in a trivial `SeriesData`.

### Symbolic prefactors

A factor of `x^alpha` with `alpha` symbolic (non-integer, non-rational) is pulled outside the expansion so the remaining body is expanded as an ordinary power series. For example, `Series[x^a Exp[x], {x, 0, 5}]` returns `x^a (1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + O[x]^6)` — a `Times[Power[x, a], SeriesData[...]]` at the expression level, so the `SeriesData` pretty-printer still renders the body and the outer `Times` decorates it with the symbolic prefactor.

**Expansion at regular points of `Arc*` heads**: When `so_apply_kernel_at_zero` can't apply (because the inner series constant `c` is not `0`), the engine falls back to naive Taylor via `D`. This makes `Series[ArcSin[x], {x, 1/2, 3}]`, `Series[ArcTan[x], {x, 2, 2}]`, `Series[ArcSinh[x], {x, 1, 2}]`, etc. work without special-casing each non-zero expansion point.

### Maxima-style algebraic fast paths

- **Monomial binomial** `(a + b x^m)^alpha` with `alpha` non-integer. The generic path factors out `a`, forms `u = (b/a) x^m`, and feeds `(1+u)^alpha` through Horner composition. When `u` is a single-term series (exactly one non-zero coefficient in `SeriesObj` terms) we skip the `O(N^2)` convolution and emit `binomial(alpha, k) * (b/a)^k` directly at exponent `k*m`. This covers `Sqrt[1+x]`, `(1 - x^2)^(1/2)`, `(1+x)^(-1/2)`, `(2 + 3x)^(1/3)`, `(1+x)^n` with symbolic `n`, and Puiseux bases like `(1+x^(1/2))^(1/2)`.
- **split-two-term probe** (`series_split_two_term` in `series.h`). Structural decomposition of `e` into `a + b*x^(p/q)` without running the full series-expansion pipeline. Feeds the Log fast path and the Apart gate; exposed for unit testing.
- **Log fast path**: when `arg` matches `a + b x^(p/q)` with `a, b` both free of `x` and `a != 1`, rewrite `Log[a + b x^c]` as `Log[a] + Log[1 + (b/a) x^c]` and let the `Log1p` kernel compose with a pure monomial. Maxima's `sp2log` uses the same identity.
- **Apart preprocessing**: if the input contains `Power[p(x), -n]` for `p` a polynomial in `x`, run `Apart[f, x]` to decompose into partial fractions before expanding. Composite denominators like `1/((1-x)(1-2x)(1-3x))` then break up into geometric-series pieces that hit the monomial fast path. Gated by a polynomial check so non-rational denominators (e.g. `1/(Exp[x] - 1 - x)`) bypass Apart and fall through to the generic `so_inv` path.

**Branch-point expansion for inverse trig / hyperbolic heads** (MMA-faithful). All eight `Arc*` heads expand at their branch points with the same wrapped output shape:

The `(-1)^Floor[...]` factor is the MMA branch discriminator — it is `1` on the principal sheet near `x0` and flips sign across the branch cut. Two mathematical families:

- **Family A — square-root branches** (derivative ~ `1/Sqrt[(x - x0) * linear]`): `ArcSin` / `ArcCos` at `x = ±1`; `ArcSinh` at `x = ±I`; `ArcCosh` at `x = ±1`. Output is a Puiseux series with `den = 2` and no `Log` term. Derived from the identities `ArcCos[1 - s] = Sqrt[2s] · Σ b_k s^k / (2k+1)`, `b_k = (2k)! / (8^k (k!)^2)`, `ArcSinh[σI + u] = σ·I·π/2 + 2 ArcSinh[Sqrt[u / (2σI)]]`, and the principal-branch `ArcCosh[1 + u] = 2 ArcSinh[Sqrt[u/2]]`.
- **Family B — logarithmic branches** (derivative has a simple pole at `x0`): `ArcTan` / `ArcCot` at `x = ±I`; `ArcTanh` / `ArcCoth` at `x = ±1`. Output contains an explicit `Log[x - x0]` term with its own coefficient (e.g. `-1/2` for `ArcTanh@1`, `-I/2` for `ArcTan@I`) plus a regular power series with `den = 1`. Derived from the identities `ArcTanh[x] = (1/2) Log[(1 + x)/(1 - x)]` and `ArcTan[z] = (1/(2 I)) Log[(1 + I z)/(1 - I z)]`.

The handler fires when the inner series at `x0` is exactly `c + q (x - x0)` (constant plus linear, with `q ≠ 0`) and `c` matches the branch-point value. For nested cases (e.g. `Sin[ArcSinh[x]]` near `x = I`), composition is preserved by emitting a constant-inside `SeriesObj` instead of the wrapper.

Examples:

`Normal[Series[ArcTanh[x], {x, 1, 3}]]` preserves the `Log[x - 1]` term and the branch discriminator — `Normal` collapses the inner `SeriesData` but the surrounding `Plus`/`Times` pass through unchanged (matching MMA).

### Naive-Taylor fallback

For expansion points that are not branch points (e.g. `Series[ArcSinh[x], {x, 1 + I, 3}]`, `Series[ArcSin[x], {x, 1/2, 3}]`), `series_expand` falls back to naive Taylor via repeated `D`. The fallback caps iterations and bails out on `Infinity` / `Indeterminate` derivatives so unknown heads cannot spin the engine. The singularity probe **evaluates** the substituted `f^(k)(x0)` with arithmetic warnings muted before testing for infinities, so an *unevaluated* pole produced by the substitution (e.g. `f[1/x]` at `x = 0` becoming `f[1/0]`) collapses to `ComplexInfinity` / `Indeterminate` and is detected cleanly — it no longer spills a spurious `Power::infy: 1/0` to stderr before the guard fires. So heads with no recognised expansion at the requested point (e.g. `Series[Gamma[x], {x, Infinity, 2}]`) return unevaluated silently rather than with a stray warning.

## Algorithm

============================================================================ series.c - Series and SeriesData ============================================================================

This module implements the power-series machinery for Mathilda.

SeriesData[x, x0, {a0, ..., a_{k-1}}, nmin, nmax, den] is the data head that represents a truncated power series. The i-th coefficient multiplies (x - x0)^((nmin + i)/den) and an O[x - x0]^(nmax/den) term captures the dropped higher-order terms.

```text
Series[f, {x, x0, n}]  expands f as a power series in (x - x0) up to
```

order n. Series also accepts the leading-term form Series[f, x -> x0] and the iterated multivariate form Series[f, {x, x0, nx}, {y, y0, ny}, ...]. The algorithm is a recursive "series algebra": primitive subexpressions become SeriesObj's, algebraic heads (Plus, Times, Power) combine them, and elementary heads (Exp, Log, Sin, Cos, Sinh, Cosh, Tan, Tanh) apply their known series kernels. Unknown heads fall back to naive Taylor via D[...]. Expansion about Infinity is handled by substituting x -> 1/u internally and presenting the result with Power[x, -1] as the series variable.

Normal[s] drops the O-term from a SeriesData and returns an ordinary sum.

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Limit x (Log[x+1]-Log[x]) at Infinity | 1.24 s | 2.09 s | -- |
| Limit Sin[x]/x at 0 | 0.963 s | 1.37 s | -- |
| Limit (Exp[x]-1-x)/x^2 at 0 | 0.218 s | 1.63 s | -- |
| Series Exp[Sin[x]] to order 20 | -- | -- | -- |
| Series 1/(1-x-x^2) to order 60 | -- | -- | -- |
| Series Log[1+Sin[x]] to order 24 | -- | -- | -- |

## Implementation notes

**Algorithm.** Series computes a truncated power series via recursive *series
algebra*. `builtin_series` (`ATTR_HOLDALL`) parses each spec — full form
`{x, x0, n}` or leading-term form `x -> x0` (via `parse_series_spec`) — threads a
`List` first argument, and handles the multivariate form by expanding left-to-
right (each coefficient of the outer series is recursively expanded in the next
variable). Inexact inputs are rationalised then numericalised back
(`internal_rationalize_then_numericalize`).

`do_series_single` does the work. It evaluates `f`, returns it verbatim if free of
`x`, optionally Apart-decomposes a rational function in `x` (each partial-fraction
term hits the cheap monomial binomial path instead of a Newton inversion), pulls
out a symbolic `x^alpha` prefactor when expanding at 0 or Infinity, and handles
expansion at Infinity by substituting `x -> 1/u`. It expands to a padded internal
order (`order + pad`, pad = 12 for numeric `x0`, 2 for symbolic `x0` to keep
symbolic convolution from exploding) and then truncates back to the user O-term;
the leading-term form widens the O-term to the first non-zero coefficient.

The recursion is `series_expand(e, ctx)`: a subexpression free of `x` becomes a
constant series; `x` itself becomes the identity series; `Plus`/`Times`/`Power`
combine child series with `so_add`/`so_mul`/`so_inv`/`so_pow_int`; the elementary
heads `Exp, Log, Sin, Cos, Tan, Sinh, Cosh, Tanh` and the inverse trig/hyperbolic
family compose their known Taylor kernels (`kernel_coefs` →
`so_compose_scalar_kernel`) with the inner argument's series, with reciprocal
heads (`Sec`/`Csc`/`Cot`/...) rewritten via `rewrite_reciprocal_head` and a large
set of branch-point / at-infinity identities for the inverse functions. Any
unrecognised head falls back to naive Taylor `series_taylor_via_D`: coefficients
`a_k = (D^k f at x0)/k!`, capped at `MAX_NAIVE_ORDER` = 20 and bailing out
(`has_infinity`) at singularities.

**Data structures.** The internal `SeriesObj` struct holds the expansion variable
`x`, expansion point `x0`, an owned array of coefficient `Expr*`, the leading
exponent numerator `nmin`, the (exclusive) O-term exponent numerator `order`, and
a common exponent denominator `den` (>= 1) — so coefficient `i` multiplies
`(x - x0)^((nmin+i)/den)`, supporting Laurent and Puiseux (fractional-exponent)
series. `so_rescale`/`so_align_den` reconcile denominators before arithmetic.
`SeriesObj` is converted to the user-facing `SeriesData[x, x0, {coefs}, nmin,
nmax, den]` head by `so_to_expr`.

**Complexity / limits.** Kernel-path heads are exact and fast; the naive D-path is
capped at order 20 and fails at true branch points where derivatives blow up
(Puiseux at such points is out of scope except for the explicit inverse-function
branch-point handlers). Symbolic expansion points cap the internal pad tightly to
avoid `O(N^2)` symbolic coefficient blow-up.

- `HoldAll` and `Protected` (so the expansion variable is not evaluated before `Series` has a chance to shield it).
- Threaded over lists: `Series[{f1, f2, ...}, spec]` becomes `{Series[f1, spec], Series[f2, spec], ...}`.
- Handles Taylor expansions for smooth functions, Laurent expansions where the function has a pole at `x0`, Puiseux expansions for fractional-power cases such as `Sqrt[Sin[x]]`, and logarithmic expansions for cases like `x^x` where `Log[x]` survives as a symbolic coefficient.
- Symbolic parameters in exponents are supported: `Series[(1 + x)^n, {x, 0, 4}]` returns the binomial expansion with `n` kept unexpanded.
- Approximate numeric coefficients flow through series arithmetic unchanged.
- For unknown heads (e.g. `f[x]` where `f` has no rules), the engine falls back to naive Taylor via `D` at the expansion point; the coefficients appear as `Derivative[k][f][x0]`.
- **Assumptions** — `Series` honours assumptions through the same three channels as `Limit` and `PossibleZeroQ`: the `Assumptions -> assm` **option** (matched by its LHS symbol `Assumptions`, so it is not confused with a leading-term spec `x -> x0`), an ambient **`Assuming[assm, ...]`** scope, and a direct **`$Assumptions`** assignment (e.g. via `Block[{$Assumptions = assm}, ...]`). The option overrides the ambient value; a `True`/`Automatic`/inconsistent/uninformative assumption reproduces the assumption-free result byte-for-byte. Assumptions reach the expansion in three places:
  - *Sign of the expansion variable* — selects the branch of the logarithmic expansions at `x = 0` (`ExpIntegralEi`, `LogIntegral`, `CosIntegral`, `CoshIntegral`). Under `x < 0` (proved by the assumption context) `Log[x]` is emitted as `Log[-x]`; otherwise the principal `x > 0` form is used.
  - *Non-analytic heads of the expansion variable* — `Abs[x]`, `Sign[x]`, `UnitStep[x]`, `Sqrt[x^2]`, `Conjugate[x]` have no Taylor series, but collapse to an analytic form once the assumptions pin the argument's sign or reality: `Series[Abs[x], {x, 0, 3}, Assumptions -> x > 0]` = `x + O[x]^4` (and `-x + O[x]^4` for `x < 0`); `Sign[x] -> ±1`, `UnitStep[x] -> 1`/`0`, `Conjugate[x] -> x` for real `x`. Without an assumption these fall back to the naive `Derivative[Abs][0]` Taylor row, as before.
  - *Coefficient cleanup* — the final coefficients are simplified under the assumptions via the shared `apply_assumption_rules` rewriter (no `Simplify` search): `Sqrt[a^2] -> a`, `Abs[a] -> a` (for `a > 0`; `-a`/`Abs[a]` for `a < 0`/real), `(a^2)^r -> a^(2r)`, `Log[a^p] -> p Log[a]`, `Sign[a] -> ±1`, and integer/parity facts such as `Cos[n Pi] -> (-1)^n`. So `Series[Sqrt[a^2 + x], {x, 0, 2}, Assumptions -> a > 0]` = `a + x/(2a) - x^2/(8 a^3) + O[x]^3`, whereas the no-assumption form keeps `Sqrt[a^2]`. Only subterms keyed on a proven symbol are touched (a coefficient free of such symbols — the entire no-assumption case — is returned unchanged), and `Element[a, Complexes]` correctly leaves `Sqrt[a^2]` alone.

  The assumption is forwarded into each inner variable for multivariate expansions.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [SeriesData](../../power-series/SeriesData/), [HoldAll](../../expression-information/HoldAll/), [D](../../calculus/D/), [Limit](../../calculus/Limit/), [PossibleZeroQ](../../expression-information/PossibleZeroQ/), [$Assumptions](../../simplification/$Assumptions/), [ExpIntegralEi](../../special-functions/ExpIntegralEi/), [LogIntegral](../../special-functions/LogIntegral/)

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 3.
- Joel S. Cohen, *Computer Algebra and Symbolic Computation: Mathematical Methods* (A K Peters, 2003).
- Source: [`src/calculus/series.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/series.c)
- Specification: [`docs/spec/builtins/power-series.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/power-series.md)
- Tests: [`tests/test_airyai.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airyai.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_besseli.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besseli.c)
- Tests: [`tests/test_besselj.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besselj.c)

## Notes & additional examples

### Notes

`Series[f, {x, x0, n}]` builds a power-series expansion to order `(x - x0)^n`, returning a `SeriesData` object that prints with a trailing `O[x]^(n+1)` term. It handles Taylor, Laurent (negative powers), and Puiseux (fractional powers) cases, as well as expansion around `Infinity` via the `x -> 1/u` substitution. Apply `Normal` to drop the order term and recover an ordinary polynomial, as in the `Exp` example above. `Series` is `HoldAll`, so the expansion variable is held unevaluated while the expansion point and order are read off.
