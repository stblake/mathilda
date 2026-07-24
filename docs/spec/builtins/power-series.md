# Power Series

## SeriesData
Represents an explicit truncated power series.
- `SeriesData[x, x0, {a0, a1, ..., a_{k-1}}, nmin, nmax, den]`

The `ai` are the coefficients of the series about the expansion point `x0`.
The power of `(x - x0)` attached to `ai` is `(nmin + i)/den`, and a trailing
`O[x - x0]^(nmax/den)` term indicates the order at which higher terms have
been dropped.

**Features**:
- `Protected`.
- `SeriesData` is a pure data head; it has no evaluator and is normally
  produced by `Series`.
- Standard printing renders the series as an ordinary mathematical sum:
  `a0 + a1 (x - x0) + a2 (x - x0)^2 + ... + O[x - x0]^p`. Zero
  coefficients are suppressed, and `x0 == 0` is displayed as simply `x`
  without the subtraction.
- `InputForm[...]` switches to the literal `SeriesData[x, x0, {...},
  nmin, nmax, den]` form, which round-trips through the parser.
- `FullForm[...]` shows the raw tree structure.

```mathematica
In[1]:= SeriesData[x, 0, {1, 1, 1/2, 1/6, 1/24, 1/120}, 0, 6, 1]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + 1/120 x^5 + O[x]^6

In[2]:= InputForm[%]
Out[2]= SeriesData[x, 0, {1, 1, 1/2, 1/6, 1/24, 1/120}, 0, 6, 1]

In[3]:= SeriesData[x, 0, Table[i^2, {i, 10}], 0, 10, 1]
Out[3]= 1 + 4 x + 9 x^2 + 16 x^3 + 25 x^4 + 36 x^5 + 49 x^6 + 64 x^7 + 81 x^8 + 100 x^9 + O[x]^10

In[4]:= SeriesData[x, 2, {a, b, c}, 0, 3, 1]
Out[4]= a + b (x - 2) + c (x - 2)^2 + O[x - 2]^3

In[5]:= SeriesData[x, 0, {1, 2, 3}, 1, 7, 2]
Out[5]= Sqrt[x] + 2 x + 3 x^(3/2) + O[x]^(7/2)
```

### Calculus on SeriesData

`D` and `Integrate` operate on a `SeriesData` term-by-term, returning a new
`SeriesData` (matching Mathematica).

- **Differentiation w.r.t. the series variable** applies the power rule to each
  term: the coefficient of `(x - x0)^((nmin+i)/den)` is multiplied by
  `(nmin+i)/den` and the exponent drops by one. Both `nmin` and `nmax` decrease
  by `den`; the differentiated constant term becomes a leading zero and is
  trimmed.
- **Integration w.r.t. the series variable** raises each term: the coefficient
  is multiplied by `den/(nmin+i+den)` and the exponent rises by one. Both `nmin`
  and `nmax` increase by `den`. The constant of integration is taken to be `0`.
  A genuine `(x - x0)^-1` term (nonzero residue) integrates to a `Log` that
  `SeriesData` cannot represent, so `Integrate` is left **unevaluated** in that
  case.
- **With respect to a different variable** (when the expansion point `x0` is
  free of that variable), `D`/`Integrate` thread into the coefficients and keep
  the powers of `(x - x0)` unchanged.

```mathematica
In[1]:= Integrate[Series[Exp[x], {x, 0, 8}], x]
Out[1]= x + 1/2 x^2 + 1/6 x^3 + ... + 1/362880 x^9 + O[x]^10

In[2]:= D[Series[Exp[x], {x, 0, 8}], x]
Out[2]= 1 + x + 1/2 x^2 + ... + 1/5040 x^7 + O[x]^8

In[3]:= Integrate[Series[1/x^2 + 1, {x, 0, 3}], x]
Out[3]= -1/x + x + O[x]^5
```

### Arithmetic on SeriesData

`Plus`, `Times`, and `Power` combine `SeriesData` objects (and `Divide` /
`Subtract`, which reduce to `Times[a, Power[b,-1]]` and `Plus[a, Times[-1,b]]`).
All operands are converted to series about the common `(x, x0)` and folded with
the internal series algebra; the result is truncated to the **minimum** O-term
order of the operands.

- **Operands that are not series** are expanded about the same point and order:
  a scalar (free of `x`) folds into the constant (`a0`) coefficient; the bare
  variable, polynomials, and transcendental functions are series-expanded with
  the full engine. Mixed exact/approximate coefficients (integer, bigint,
  rational, machine real, MPFR) combine per-coefficient, so adding a real folds
  it into `a0` while leaving the other coefficients exact, as in Mathematica.
- **Incompatible operands** — series about a different variable or expansion
  point — leave the `Plus`/`Times`/`Power` **unevaluated** (a symbolic sum or
  product of the `SeriesData` objects).
- **Power**: integer exponents (including negative, via series inversion) and
  exponents free of the series variable are handled directly. An exponent that
  depends on the series variable, or an ordinary base raised to a series
  exponent (e.g. `2^Series[...]`), is computed as `Exp[exp*Log[base]]`.

```mathematica
In[1]:= Series[Exp[x], {x, 0, 2}] + 1
Out[1]= 2 + x + x^2/2 + O[x]^3

In[2]:= x Series[Sin[x], {x, 0, 5}]
Out[2]= x^2 - x^4/6 + x^6/120 + O[x]^7

In[3]:= Series[Exp[x], {x, 0, 2}] Series[Exp[x], {x, 0, 3}]
Out[3]= 1 + 2 x + 2 x^2 + O[x]^3

In[4]:= Series[Exp[x], {x, 0, 2}]^3
Out[4]= 1 + 3 x + 9 x^2/2 + O[x]^3

In[5]:= 2^Series[x, {x, 0, 3}]
Out[5]= 1 + Log[2] x + Log[2]^2/2 x^2 + Log[2]^3/6 x^3 + O[x]^4
```

## Series
Produces the power-series expansion of an expression about a point.
- `Series[f, {x, x0, n}]` — Taylor/Laurent/Puiseux expansion up to order `(x - x0)^n`.
- `Series[f, x -> x0]` — leading-term form. The engine scans the internal expansion for the first non-zero coefficient at exponent `e1` and the next non-zero at `e2 > e1`; the reported `O`-term lands at exponent `e2` (or `e1 + 1` when no further non-zero term exists). So `Series[Sin[x] - x, x -> 0]` returns `-x^3/6 + O[x]^5`, `Series[f[x], x -> 0]` returns `f[0] + O[x]^1`, and analytic-at-x0 inputs collapse to their constant plus `O[x - x0]^1`.
- `Series[f, {x, x0, nx}, {y, y0, ny}, ...]` — iterated multivariate expansion. Each inner coefficient is itself a `SeriesData` in the next variable.
- `Series[f, {x, Infinity, n}]` — expansion at infinity, substituting `x -> 1/u` internally. The emitted `SeriesData` uses `Power[x, -1]` as its variable, so the series prints with `1/x` as the base and an `O[1/x]^(n+1)` term.

**Features**:
- `HoldAll` and `Protected` (so the expansion variable is not evaluated before `Series` has a chance to shield it).
- Threaded over lists: `Series[{f1, f2, ...}, spec]` becomes `{Series[f1, spec], Series[f2, spec], ...}`.
- Handles Taylor expansions for smooth functions, Laurent expansions where the function has a pole at `x0`, Puiseux expansions for fractional-power cases such as `Sqrt[Sin[x]]`, and logarithmic expansions for cases like `x^x` where `Log[x]` survives as a symbolic coefficient.
- Symbolic parameters in exponents are supported: `Series[(1 + x)^n, {x, 0, 4}]` returns the binomial expansion with `n` kept unexpanded.
- Approximate numeric coefficients flow through series arithmetic unchanged.
- For unknown heads (e.g. `f[x]` where `f` has no rules), the engine falls back to naive Taylor via `D` at the expansion point; the coefficients appear as `Derivative[k][f][x0]`.
- The `Assumptions -> assm` option selects the branch of the logarithmic expansions at `x = 0` (`ExpIntegralEi`, `LogIntegral`). When `assm` forces the expansion variable negative (e.g. `x < 0`, `x < -2`, `0 > x`, or an `And[...]` containing such a relation), `Log[x]` is emitted as `Log[-x]`; otherwise the principal `x > 0` form is used. The option is matched by its LHS symbol `Assumptions` (so it is not confused with a leading-term spec `x -> x0`) and is forwarded into each inner variable for multivariate expansions.

**Coefficient arithmetic** automatically promotes to BigInt-backed `Rational` when 64-bit numerators or denominators would overflow, so previously-overflowing Laurent/Puiseux cases like `Series[1/Sin[x]^10, {x, 0, 2}]` and `Series[Sqrt[Log[1 + x]], {x, 0, 12}]` now produce exact coefficients (at the cost of slower evaluation for large orders).

Inverse trigonometric and inverse hyperbolic heads (`ArcSin`, `ArcCos`, `ArcTan`, `ArcCot`, `ArcSinh`, `ArcCosh`, `ArcTanh`, `ArcCoth`) are handled by direct series kernels at `u = 0` rather than by naive repeated differentiation, which would blow up expression size exponentially for higher orders. `ArcCosh` uses the principal-branch identity `ArcCosh[u] = I*ArcCos[u]`, so its expansion at `x = 0` has the expected `I*Pi/2` constant term and imaginary coefficients.

Forward reciprocal heads (`Sec`, `Csc`, `Cot`, `Sech`, `Csch`, `Coth`) are rewritten as `1/Cos[x]`, `1/Sin[x]`, `Cos[x]/Sin[x]`, etc., before expansion. Inverse reciprocal heads (`ArcSec`, `ArcCsc`, `ArcSech`, `ArcCsch`) are rewritten via the identities `ArcSec[z] = ArcCos[1/z]`, `ArcCsc[z] = ArcSin[1/z]`, `ArcSech[z] = ArcCosh[1/z]`, `ArcCsch[z] = ArcSinh[1/z]`, so that a blowing-up inner series (e.g. `z = 1/x`) collapses to a convergent kernel case rather than triggering spurious `Power::infy` warnings.

Expansions where the inner series diverges at the expansion point (e.g. `Series[f[1/x], {x, 0, n}]`) are handled via dedicated at-infinity identities:
- `ArcCoth[1/u] = ArcTanh[u]`, `ArcCot[1/u] = ArcTan[u]` (handled at the series level via inner-series inversion).
- `ArcTanh[1/u] = I*Pi/2 + ArcTanh[u]` (principal branch).
- `ArcSinh[1/v] = -Log[v] + Log[1 + Sqrt[1 + v^2]]` and `ArcCosh[1/v] = -Log[v] + Log[1 + Sqrt[1 - v^2]]` (handled by rewriting at the expression level; the symbolic `-Log[x]` term rides the existing `Log[x]` symbolic-coefficient path).

**Asymptotic expansions at Infinity for special functions with essential singularities**: Some functions have no Laurent series at `Infinity` because their leading behaviour is an essential singularity (a factor like `E^x`). For these the generic `x -> 1/u` substitution would hand a pole to naive Taylor, so they are emitted from dedicated asymptotic identities with the essential factor kept symbolic:
- `Series[ExpIntegralEi[x], {x, Infinity, n}]` returns `E^x (1/x + 1/x^2 + 2/x^3 + ... + O[1/x]^(n+1))`, i.e. `Times[Power[E, x], SeriesData[Power[x, -1], 0, {0!, 1!, ..., (n-1)!}, 1, n+1, 1]]` (DLMF 6.12.2: `Ei(x) ~ E^x Σ_{k≥0} k!/x^(k+1)`). The `E^x` factor rides the expression-level `Times` exactly as a symbolic `x^alpha` prefactor does.
- `Series[Erf[x], {x, Infinity, n}]`, `Series[Erfc[x], ...]`, `Series[Erfi[x], ...]` return the error-function asymptotic expansions (DLMF 7.12.1), each an `Exp[±x^2]` essential-singularity prefactor times a Laurent series in `1/x` with only odd powers populated. `Series[Erf[x], {x, Infinity, 2}]` → `1 + E^(-x^2) (-1/(Sqrt[Pi] x) + O[1/x]^3)` (the leading `1` is the limit); `Erfc = 1 - Erf` uses the negated multiplier without the constant; `Erfi` has all-positive coefficients and a growing `Exp[+x^2]` prefactor. The general coefficient is `a_k = (2k-1)!!/(2^k Sqrt[Pi])` at `x^-(2k+1)`, with alternating signs for `Erf`/`Erfc`.
- `Series[LogGamma[x], {x, Infinity, n}]` returns the Stirling expansion (DLMF 5.11.1): an *additive* growth head `(x - 1/2) Log[x] - x + Log[2 Pi]/2` (kept symbolic, as the `E^x` prefactor is for `ExpIntegralEi`) plus a Bernoulli Laurent tail `Σ_{k≥1} B_{2k}/(2k(2k-1)) x^-(2k-1)` in `1/x`, i.e. `1/(12x) - 1/(360 x^3) + 1/(1260 x^5) - ...`. `Series[LogGamma[x], {x, Infinity, 2}]` → `1/2 Log[2 Pi] - x + Log[x] (-1/2 + x) + 1/(12 x) + O[1/x]^3`. (`Gamma[x] = Exp[LogGamma[x]]` still has no series at infinity — it diverges — so `Series[Gamma[x], {x, Infinity, n}]` stays unevaluated.)

**Logarithmic series at `x = 0` for special functions**: Some functions have a finite-radius series at the origin built around a single `Log[x]` branch term, with the log baked into the `x^0` coefficient. Naive Taylor-via-`D` cannot reach these (`f(0)` is infinite or the derivatives have poles), so they are emitted from dedicated identities:
- `Series[ExpIntegralEi[x], {x, 0, n}]` returns `EulerGamma + Log[x] + x + x^2/4 + x^3/18 + ... + O[x]^(n+1)`, i.e. `SeriesData[x, 0, {EulerGamma + Log[x], 1, 1/4, ..., 1/(n*n!)}, 0, n+1, 1]` (DLMF 6.6.2: `Ei(x) = EulerGamma + Log[x] + Σ_{k≥1} x^k/(k k!)`). The `EulerGamma + Log[x]` branch term occupies the `x^0` slot. With `Assumptions -> x < 0` the branch term becomes `EulerGamma + Log[-x]`.

**Generalized series in `Log[x]` at `x = 0`**: `LogIntegral[x] = Ei(Log[x])` has *no* ordinary Taylor (or Laurent/Puiseux) series at `x = 0`: as `x -> 0+`, `L = Log[x] -> -Infinity` drives it into Ei's asymptotic regime, so the result is a generalized series whose `x^1` coefficient is itself a series in `1/Log[x]`.
- `Series[LogIntegral[x], {x, 0, n}]` returns `x Σ_{k=0}^{2n+1} k!/Log[x]^(k+1) + O[x]^(n+1)`. For `n = 2` this is `((120 + 24 Log[x] + 6 Log[x]^2 + 2 Log[x]^3 + Log[x]^4 + Log[x]^5) x)/Log[x]^6 + O[x]^3`, matching Mathematica's `Assumptions -> x > 0` output. Since the only `x`-dependence is the prefactor `E^L = x`, every term carries exactly `x^1`; the coefficient slot is emitted as `Together[Σ k!/Log[x]^(k+1)]`. Mathematica's no-assumptions form additionally wraps this in a `Floor[Arg[...]]` branch discriminator (to track `Log[1/x]` vs `Log[x]` across the cut); Mathilda emits the principal `x > 0` form by default, or — with `Assumptions -> x < 0` — the `x < 0` form, where every `Log[x]` becomes `Log[-x]` and an additive `I Pi` (the `x^0` term, from `Log[x] = Log[-x] + I Pi`) leads the series: `I Pi + ((120 + 24 Log[-x] + ... + Log[-x]^5) x)/Log[-x]^6 + O[x]^3`.

**Internal padding for symbolic expansion points**: The engine computes series at a padded internal order (user order + 12 by default) so that intermediate Laurent/Puiseux operations don't lose accuracy. When the expansion point `x0` is not a literal number, padding is capped at 2 — at a symbolic point the series coefficients are themselves symbolic expressions (e.g. `Cosh[a]`, `Sinh[a]`), and the `O(N^2)` convolution inside `so_inv`/`so_div` would otherwise spin indefinitely on exponentially growing expression trees. This makes cases like `Series[Coth[x], {x, a, 1}]`, `Series[Tanh[x], {x, a, 1}]`, `Series[Sec[x], {x, a, 1}]`, and `Series[1/Cosh[x], {x, a, 1}]` terminate in milliseconds.

**Constant inputs**: If `f` is free of the expansion variable (e.g. `Series[0, {x, 0, 4}]`, `Series[Sin[y], {x, 0, 4}]`, `Series[a + b^2, {x, 0, 3}]`), `Series` returns `f` verbatim instead of wrapping it in a trivial `SeriesData`.

**Symbolic prefactors**: A factor of `x^alpha` with `alpha` symbolic (non-integer, non-rational) is pulled outside the expansion so the remaining body is expanded as an ordinary power series. For example, `Series[x^a Exp[x], {x, 0, 5}]` returns `x^a (1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + O[x]^6)` — a `Times[Power[x, a], SeriesData[...]]` at the expression level, so the `SeriesData` pretty-printer still renders the body and the outer `Times` decorates it with the symbolic prefactor.

**Expansion at regular points of `Arc*` heads**: When `so_apply_kernel_at_zero` can't apply (because the inner series constant `c` is not `0`), the engine falls back to naive Taylor via `D`. This makes `Series[ArcSin[x], {x, 1/2, 3}]`, `Series[ArcTan[x], {x, 2, 2}]`, `Series[ArcSinh[x], {x, 1, 2}]`, etc. work without special-casing each non-zero expansion point.

**Maxima-style algebraic fast paths**:
- **Monomial binomial** `(a + b x^m)^alpha` with `alpha` non-integer. The generic path factors out `a`, forms `u = (b/a) x^m`, and feeds `(1+u)^alpha` through Horner composition. When `u` is a single-term series (exactly one non-zero coefficient in `SeriesObj` terms) we skip the `O(N^2)` convolution and emit `binomial(alpha, k) * (b/a)^k` directly at exponent `k*m`. This covers `Sqrt[1+x]`, `(1 - x^2)^(1/2)`, `(1+x)^(-1/2)`, `(2 + 3x)^(1/3)`, `(1+x)^n` with symbolic `n`, and Puiseux bases like `(1+x^(1/2))^(1/2)`.
- **split-two-term probe** (`series_split_two_term` in `series.h`). Structural decomposition of `e` into `a + b*x^(p/q)` without running the full series-expansion pipeline. Feeds the Log fast path and the Apart gate; exposed for unit testing.
- **Log fast path**: when `arg` matches `a + b x^(p/q)` with `a, b` both free of `x` and `a != 1`, rewrite `Log[a + b x^c]` as `Log[a] + Log[1 + (b/a) x^c]` and let the `Log1p` kernel compose with a pure monomial. Maxima's `sp2log` uses the same identity.
- **Apart preprocessing**: if the input contains `Power[p(x), -n]` for `p` a polynomial in `x`, run `Apart[f, x]` to decompose into partial fractions before expanding. Composite denominators like `1/((1-x)(1-2x)(1-3x))` then break up into geometric-series pieces that hit the monomial fast path. Gated by a polynomial check so non-rational denominators (e.g. `1/(Exp[x] - 1 - x)`) bypass Apart and fall through to the generic `so_inv` path.

**Branch-point expansion for inverse trig / hyperbolic heads** (MMA-faithful). All eight `Arc*` heads expand at their branch points with the same wrapped output shape:

```
Plus[ f(x0),
      Times[ log_coef, Log[x - x0] ],         (Family B only)
      Times[ (-1)^Floor[(Pi/2 - Arg[x - x0])/(2 Pi)],
             SeriesData[...] ] ]
```

The `(-1)^Floor[...]` factor is the MMA branch discriminator — it is `1` on the principal sheet near `x0` and flips sign across the branch cut. Two mathematical families:

- **Family A — square-root branches** (derivative ~ `1/Sqrt[(x - x0) * linear]`): `ArcSin` / `ArcCos` at `x = ±1`; `ArcSinh` at `x = ±I`; `ArcCosh` at `x = ±1`. Output is a Puiseux series with `den = 2` and no `Log` term. Derived from the identities `ArcCos[1 - s] = Sqrt[2s] · Σ b_k s^k / (2k+1)`, `b_k = (2k)! / (8^k (k!)^2)`, `ArcSinh[σI + u] = σ·I·π/2 + 2 ArcSinh[Sqrt[u / (2σI)]]`, and the principal-branch `ArcCosh[1 + u] = 2 ArcSinh[Sqrt[u/2]]`.
- **Family B — logarithmic branches** (derivative has a simple pole at `x0`): `ArcTan` / `ArcCot` at `x = ±I`; `ArcTanh` / `ArcCoth` at `x = ±1`. Output contains an explicit `Log[x - x0]` term with its own coefficient (e.g. `-1/2` for `ArcTanh@1`, `-I/2` for `ArcTan@I`) plus a regular power series with `den = 1`. Derived from the identities `ArcTanh[x] = (1/2) Log[(1 + x)/(1 - x)]` and `ArcTan[z] = (1/(2 I)) Log[(1 + I z)/(1 - I z)]`.

The handler fires when the inner series at `x0` is exactly `c + q (x - x0)` (constant plus linear, with `q ≠ 0`) and `c` matches the branch-point value. For nested cases (e.g. `Sin[ArcSinh[x]]` near `x = I`), composition is preserved by emitting a constant-inside `SeriesObj` instead of the wrapper.

Examples:

```mathematica
In[]:= Series[ArcSinh[x], {x, I, 3}]
Out[]= I Pi/2 + (-1)^Floor[(Pi/2 - Arg[x - I])/(2 Pi)] (
         (1 - I) Sqrt[x - I] + (1 + I)/12 (x - I)^(3/2)
         + (-3 + 3 I)/160 (x - I)^(5/2) + O[x - I]^(7/2) )

In[]:= Series[ArcTanh[x], {x, 1, 3}]
Out[]= Log[2]/2 + I Pi/2 - Log[x - 1]/2 + (-1)^Floor[(Pi/2 - Arg[x-1])/(2 Pi)] (
         (x - 1)/4 - (x - 1)^2/16 + (x - 1)^3/48 + O[x - 1]^4 )

In[]:= Series[ArcCosh[x + 1], {x, 0, 5}]
Out[]= (-1)^Floor[(Pi/2 - Arg[x])/(2 Pi)] (
         Sqrt[2] Sqrt[x] - Sqrt[2]/12 x^(3/2) + 3 Sqrt[2]/160 x^(5/2)
         - 5 Sqrt[2]/896 x^(7/2) + 35 Sqrt[2]/18432 x^(9/2) + O[x]^(11/2) )
```

`Normal[Series[ArcTanh[x], {x, 1, 3}]]` preserves the `Log[x - 1]` term and the branch discriminator — `Normal` collapses the inner `SeriesData` but the surrounding `Plus`/`Times` pass through unchanged (matching MMA).

**Naive-Taylor fallback**: For expansion points that are not branch points (e.g. `Series[ArcSinh[x], {x, 1 + I, 3}]`, `Series[ArcSin[x], {x, 1/2, 3}]`), `series_expand` falls back to naive Taylor via repeated `D`. The fallback caps iterations and bails out on `Infinity` / `Indeterminate` derivatives so unknown heads cannot spin the engine. The singularity probe **evaluates** the substituted `f^(k)(x0)` with arithmetic warnings muted before testing for infinities, so an *unevaluated* pole produced by the substitution (e.g. `f[1/x]` at `x = 0` becoming `f[1/0]`) collapses to `ComplexInfinity` / `Indeterminate` and is detected cleanly — it no longer spills a spurious `Power::infy: 1/0` to stderr before the guard fires. So heads with no recognised expansion at the requested point (e.g. `Series[Gamma[x], {x, Infinity, 2}]`) return unevaluated silently rather than with a stray warning.

```mathematica
In[1]:= Series[Exp[x], {x, 0, 10}]
Out[1]= 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + x^6/720 + x^7/5040 + x^8/40320 + x^9/362880 + x^10/3628800 + O[x]^11

In[2]:= Series[f[x], {x, a, 3}]
Out[2]= f[a] + Derivative[1][f][a] (x - a) + 1/2 Derivative[2][f][a] (x - a)^2 + 1/6 Derivative[3][f][a] (x - a)^3 + O[x - a]^4

In[3]:= Series[Cos[x]/x, {x, 0, 10}]
Out[3]= 1/x - x/2 + x^3/24 - x^5/720 + x^7/40320 - x^9/3628800 + O[x]^11

In[4]:= Series[Sqrt[Sin[x]], {x, 0, 10}]
Out[4]= Sqrt[x] - x^(5/2)/12 + x^(9/2)/1440 - x^(13/2)/24192 - 67 x^(17/2)/29030400 + O[x]^(21/2)

In[5]:= Series[x^x, {x, 0, 4}]
Out[5]= 1 + Log[x] x + Log[x]^2/2 x^2 + Log[x]^3/6 x^3 + Log[x]^4/24 x^4 + O[x]^5

In[6]:= Series[(1 + x)^n, {x, 0, 4}]
Out[6]= 1 + n x + 1/2 n (-1 + n) x^2 + 1/6 n (-2 + n) (-1 + n) x^3 + 1/24 n (-3 + n) (-2 + n) (-1 + n) x^4 + O[x]^5

In[7]:= Series[Sin[1/x], {x, Infinity, 10}]
Out[7]= 1/x - 1/6 (1/x)^3 + 1/120 (1/x)^5 - 1/5040 (1/x)^7 + 1/362880 (1/x)^9 + O[1/x]^11

In[8]:= Series[Sin[x + y], {x, 0, 3}, {y, 0, 3}]
Out[8]= (y - y^3/6 + O[y]^4) + (1 - y^2/2 + O[y]^4) x + (-y/2 + y^3/12 + O[y]^4) x^2 + (-1/6 + y^2/12 + O[y]^4) x^3 + O[x]^4

In[9]:= Series[{Sin[x], Cos[x], Tan[x]}, {x, 0, 5}]
Out[9]= {x - x^3/6 + x^5/120 + O[x]^6, 1 - x^2/2 + x^4/24 + O[x]^6, x + x^3/3 + 2 x^5/15 + O[x]^6}
```

## Normal
Converts a `SeriesData` back into an ordinary expression by dropping its O-term.
- `Normal[expr]`

**Features**:
- `Protected`.
- Returns the Plus of the coefficient-times-power terms (zero coefficients skipped). For non-`SeriesData` input, `Normal` is the identity.
- Recurses through the whole expression, dropping the O-term of **every** `SeriesData` at any depth. This matters for expansions around `+-Infinity`, whose `SeriesData` is wrapped inside `Plus`/`Times` (e.g. the trig- or exponential-prefactored asymptotic forms of `BesselJ`, `BesselY`, `BesselK`, `BesselI`, `AiryAi`, `AiryBiPrime`); the surrounding factors are preserved and recombined by the evaluator.

```mathematica
In[1]:= Normal[Series[Exp[x], {x, 0, 5}]]
Out[1]= 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120

In[2]:= Normal[a + b]
Out[2]= a + b

In[3]:= Normal[Series[BesselJ[0, x], {x, Infinity, 2}]]
Out[3]= Sqrt[2/Pi] Sqrt[1/x] Cos[1/4 Pi - x] - 1/8 Sqrt[2/Pi] (1/x)^(3/2) Sin[1/4 Pi - x]
```

## SeriesCoefficient
Gives the coefficient of `(x - x0)^k` in the power-series expansion of `f`.
- `SeriesCoefficient[f, {x, x0, k}]`

**Features**:
- `HoldAll`, `Protected`.
- Computed by expanding with `Series` and extracting the `k`-th coefficient from
  the resulting `SeriesData`; general for any head `Series` can expand, with a
  concrete integer index `k` and a finite expansion point.
- Composite results (a prefactor times a `SeriesData`, e.g. asymptotic
  expansions at Infinity) and non-integer indices are left unevaluated; the
  symbolic-index general term (a Piecewise) is not produced.

```mathematica
In[1]:= SeriesCoefficient[BesselJ[0, x], {x, 0, 4}]
Out[1]= 1/64

In[2]:= SeriesCoefficient[Exp[x], {x, 0, 5}]
Out[2]= 1/120
```

