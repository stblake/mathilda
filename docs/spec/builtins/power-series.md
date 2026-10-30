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
`SeriesData` (matching the reference).

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
  the powers of `(x - x0)` unchanged. This branch introduces no exponent shift,
  so a coefficient that becomes zero (e.g. the `a`-derivative of an
  `a`-independent term) is a **genuine** zero coefficient and is retained —
  `nmin` is not raised — exactly as `Series` keeps genuine leading zeros. This
  makes `D[Series[f, {x, x0, n}], a]` equal `Series[D[f, a], {x, x0, n}]`.
  (Only the same-variable power rule above trims a leading zero, because its
  exponent shift opens a phantom boundary slot.)

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
  it into `a0` while leaving the other coefficients exact, as in the reference.
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
- `ArcTan`/`ArcCot` of an argument that blows up at the expansion point take the constant `Sign[k] Pi/2`, where `k` is the argument's leading coefficient. The sign comes from the N value of a numeric `k` (so `1/(1 - Sqrt[3])` gives `-Pi/2`), from the assumptions for a symbolic `k`, and otherwise from the reference's `Sqrt[k^2]/k` factor (v0.458). `Limit` inherits this, so `Limit[2 ArcTan[x/(1 - Sqrt[3])], x -> Infinity]` is `-Pi`.
- `Series[f, {x, Infinity, n}]` — expansion at infinity, substituting `x -> 1/u` internally. The emitted `SeriesData` uses `Power[x, -1]` as its variable, so the series prints in powers of `1/x` with an `O[1/x]^(n+1)` term. An integer power `(1/x)^k` is displayed as the equivalent power of `x`, as the reference does: `Series[x + 1/x, {x, Infinity, 2}]` prints `x + 1/x + O[1/x]^3` and `Series[x^2 + x, {x, Infinity, 1}]` prints `x^2 + x + O[1/x]^2` (these previously printed `1/1/x + ...` and `1/(1/x)^2 + 1/1/x + ...`, the first of which re-reads as `1/x`). A fractional power stays `(1/x)^(p/q)` (`Series[Sqrt[x + 1], {x, Infinity, 1}]` = `1/Sqrt[1/x] + 1/2 Sqrt[1/x] + O[1/x]^(3/2)`), because `(1/x)^(1/2)` is not `x^(-1/2)` on the negative reals a `-Infinity` expansion lives on. When the series is a summand of a larger sum (`Series[PolyGamma[0, x], {x, Infinity, 4}]` is `Log[x] + SeriesData[...]`) its terms are printed inline, so the signs read `Log[x] - 1/(2 x) - ...`, not `Log[x] + -1/(2 x) ...`. Every printed form re-parses, `O`-terms dropped, to the series' `Normal`.

**Features**:
- `Protected` only — **`Series` does not hold its arguments**, matching the reference (and `Limit`, `D` and `NSeries` here). So a *named* argument behaves exactly like the literal it stands for: `s = {x, 0, 3}; Series[Exp[x], s]` expands, `Series[Exp[x], Sequence @@ {{x, 0, 3}}]` splices, `l = {Sin[x], Cos[x]}; Series[l, {x, 0, 2}]` threads into a list of two series just as the literal list does, and a call that cannot be expanded echoes the evaluated argument rather than the name. (Through v0.235 these all froze: `Series` carried a stray `HoldAll`.)
- **A variable naming a symbol is resolved**, because that is what ordinary argument evaluation does: `Module[{w}, w = Unique["w"]; Normal[Series[Sqrt[1 + w^4], {w, 0, 4}]]]` expands in the generated symbol and gives `1 + w11^4/2`. This matters for any `.m` routine that expands in a freshly generated symbol — without the resolution the symbol does not occur in the expression at all, the expression is *constant* in it, and the "series" is the input itself, returned silently with no message. `SeriesCoefficient` resolves the same way.
- **The expansion variable must be a symbol.** Anything else — a variable carrying a numeric or compound value, or a literal non-symbol — declines with `Series::ivar` / `SeriesCoefficient::ivar` (`5 is not a valid variable.`) and leaves the call unevaluated, as in the reference and as `NSeries` already did. Nothing downstream type-checks the variable, so this is a correctness gate rather than tidiness: with `x = 5` in scope the spec is `{5, 0, 3}`, the "free of x" early-out cannot fire for a literal that occurs in the body, and the Taylor engine would otherwise answer `Series[Sin[5], {5, 0, 3}]` = `5 - 1/6 5^3 + O[5]^4`. A second argument that is no spec at all (`Series[f, x]`), and an unrecognised option name, are left unevaluated *silently*, as before. The message is capped at three consecutive appearances per head and re-armed by the next well-formed call, because the iterator heads are `HoldAll` and bind a fresh value per step — `Plot[Normal[Series[Sin[x], {x, 0, 5}]], {x, 0, 1}]` would otherwise print one line per sample point. That spelling needs `Plot[Evaluate[...], ...]`, exactly as in the reference.
- Threaded over lists: `Series[{f1, f2, ...}, spec]` becomes `{Series[f1, spec], Series[f2, spec], ...}`.
- Handles Taylor expansions for smooth functions, Laurent expansions where the function has a pole at `x0`, Puiseux expansions for fractional-power cases such as `Sqrt[Sin[x]]`, and logarithmic expansions for cases like `x^x` where `Log[x]` survives as a symbolic coefficient.
- `Gamma`, `PolyGamma` and `Beta` expand at their **poles** (an argument whose
  value at `x0` is a non-positive integer `−n`). The argument is shifted off the
  pole by the functional equations `Gamma[u] = Gamma[u+n+1]/(u (u+1) ... (u+n))`
  and `PolyGamma[m,u] = PolyGamma[m,u+n+1] − (−1)^m m! Σ_{k=0}^{n} (u+k)^−(m+1)`,
  and `Beta` is written as its Gamma ratio. So `Series[Gamma[s], {s,0,1}]` is
  `1/s − EulerGamma + (EulerGamma² + Pi²/6) s/2 + O[s]^2` and
  `Series[PolyGamma[s], {s,0,1}]` is `−1/s − EulerGamma + Pi² s/6 + O[s]^2`.
  Before, these were left unevaluated, and so `Limit` could not resolve
  `Limit[Gamma[s] Sin[s], s → 0]`.
- Symbolic parameters in exponents are supported: `Series[(1 + x)^n, {x, 0, 4}]` returns the binomial expansion with `n` kept unexpanded.
- Approximate numeric coefficients flow through series arithmetic unchanged.
- For unknown heads (e.g. `f[x]` where `f` has no rules), the engine falls back to naive Taylor via `D` at the expansion point; the coefficients appear as `Derivative[k][f][x0]`.
- **Assumptions** — `Series` honours assumptions through the same three channels as `Limit` and `PossibleZeroQ`: the `Assumptions -> assm` **option** (matched by its LHS symbol `Assumptions`, so it is not confused with a leading-term spec `x -> x0`), an ambient **`Assuming[assm, ...]`** scope, and a direct **`$Assumptions`** assignment (e.g. via `Block[{$Assumptions = assm}, ...]`). The option overrides the ambient value; a `True`/`Automatic`/inconsistent/uninformative assumption reproduces the assumption-free result byte-for-byte. Assumptions reach the expansion in three places:
  - *Sign of the expansion variable* — selects the branch of the logarithmic expansions at `x = 0` (`ExpIntegralEi`, `LogIntegral`, `CosIntegral`, `CoshIntegral`). Under `x < 0` (proved by the assumption context) `Log[x]` is emitted as `Log[-x]`; otherwise the principal `x > 0` form is used.
  - *Non-analytic heads, resolved by the sign at the expansion point* — `Abs[g]`, `Sign[g]`, `UnitStep[g]`, `Sqrt[g^2]`, `Conjugate[g]` have no Taylor series, but collapse to an analytic form once the assumptions pin the argument's sign or reality: `Series[Abs[x], {x, 0, 3}, Assumptions -> x > 0]` = `x + O[x]^4` (and `-x + O[x]^4` for `x < 0`); `Sign[x] -> ±1`, `UnitStep[x] -> 1`/`0`, `Conjugate[x] -> x` for real `x`. When the argument is not the bare variable, the branch is fixed by the sign of `g` **at the expansion point** `g(x0)`, since `Abs`/`Sign`/`UnitStep` are locally smooth/constant wherever the argument is nonzero: `Series[Abs[x - a], {x, 0, 2}, Assumptions -> a > 0]` = `a - x` (because `(x-a)|_{x=0} = -a < 0`), and `x - a` for `a < 0`. The sign query at the point goes through the Reduce/CAD bridge, so a compound value like `-a` or `b - c` (under `b > c`) is decided, not just a bare signed symbol. Without an assumption these fall back to the naive `Derivative[Abs][0]` Taylor row, as before.
  - *Even-order roots under `x < 0`* — when the radicand vanishes to an even order `2m` at the origin, the series machinery factors the leading monomial and takes the principal `+` root, giving `x^m`; for real `x < 0` the true value is `|x|^m`, so an **odd** extracted exponent flips the whole series sign. Under `x < 0`: `Series[Sqrt[1 - Cos[x]], {x, 0, 2}]` = `-x/Sqrt[2]` (vs `x/Sqrt[2]` for `x > 0`), and `Series[Sqrt[x^4 + a x^2], {x, 0, 3}, Assumptions -> a > 0 && x < 0]` = `-Sqrt[a] x - x^3/(2 Sqrt[a])`. This is the "hidden square" generalization of the explicit `Sqrt[x^2] -> ±x` handler (the radicand need not be literally `g^2`).
  - *Coefficient cleanup* — the final coefficients are simplified under the assumptions via the shared `apply_assumption_rules` rewriter (no `Simplify` search): `Sqrt[a^2] -> a`, `Abs[a] -> a` (for `a > 0`; `-a`/`Abs[a]` for `a < 0`/real), `(a^2)^r -> a^(2r)`, `Log[a^p] -> p Log[a]`, `Sign[a] -> ±1`, the principal imaginary branch of a half-integer power of a negative real `a^(p/2) -> I^p (-a)^(p/2)` (so `Sqrt[a] -> I Sqrt[-a]`, `1/Sqrt[a] -> -I/Sqrt[-a]` for `a < 0`), and integer/parity facts such as `Cos[n Pi] -> (-1)^n`. So `Series[Sqrt[a^2 + x], {x, 0, 2}, Assumptions -> a > 0]` = `a + x/(2a) - x^2/(8 a^3) + O[x]^3`, whereas the no-assumption form keeps `Sqrt[a^2]`; and `Series[Sqrt[a - x], {x, 0, 1}, Assumptions -> a < 0]` = `I Sqrt[-a] + (I/2) x/Sqrt[-a]` (vs the real branch `Sqrt[a] - x/(2 Sqrt[a])` for `a > 0`). Only subterms keyed on a proven symbol are touched (a coefficient free of such symbols — the entire no-assumption case — is returned unchanged), and `Element[a, Complexes]` correctly leaves `Sqrt[a^2]` alone.

  The assumption is forwarded into each inner variable for multivariate expansions.

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
- `Series[Erf[x], {x, Infinity, n}]`, `Series[Erfc[x], ...]`, `Series[Erfi[x], ...]` return the error-function asymptotic expansions (DLMF 7.12.1), each an `Exp[±x^2]` essential-singularity prefactor times a Laurent series in `1/x` with only odd powers populated. `Series[Erf[x], {x, Infinity, 2}]` → `1 + E^(-x^2) (-1/(Sqrt[Pi] x) + O[1/x]^3)` (the leading `1` is the limit); `Erfc = 1 - Erf` uses the negated multiplier without the constant; `Erfi` has all-positive coefficients and a growing `Exp[+x^2]` prefactor. The general coefficient is `a_k = (2k-1)!!/(2^k Sqrt[Pi])` at `x^-(2k+1)`, with alternating signs for `Erf`/`Erfc`.
- `Series[LogGamma[x], {x, Infinity, n}]` returns the Stirling expansion (DLMF 5.11.1): an *additive* growth head `(x - 1/2) Log[x] - x + Log[2 Pi]/2` (kept symbolic, as the `E^x` prefactor is for `ExpIntegralEi`) plus a Bernoulli Laurent tail `Σ_{k≥1} B_{2k}/(2k(2k-1)) x^-(2k-1)` in `1/x`, i.e. `1/(12x) - 1/(360 x^3) + 1/(1260 x^5) - ...`. `Series[LogGamma[x], {x, Infinity, 2}]` → `1/2 Log[2 Pi] - x + Log[x] (-1/2 + x) + 1/(12 x) + O[1/x]^3`. (`Gamma[x] = Exp[LogGamma[x]]` still has no series at infinity — it diverges — so `Series[Gamma[x], {x, Infinity, n}]` stays unevaluated.)
- `Series[PolyGamma[m, x], {x, Infinity, n}]` returns the polygamma Stirling expansion (DLMF 5.11.2). The digamma `m = 0` mirrors `LogGamma`: an additive `Log[x]` growth head plus a Laurent tail `-1/(2x) - Σ_{k≥1} B_{2k}/(2k) x^-(2k)`, so `Series[PolyGamma[0, x], {x, Infinity, 4}]` → `Log[x] - 1/(2x) - 1/(12 x^2) + 1/(120 x^4) + O[1/x]^5`. For `m ≥ 1` it is a pure Laurent series in `1/x` with leading power `x^-m` (decays to 0): `PolyGamma[m, x] ~ (-1)^(m-1) [ (m-1)!/x^m + m!/(2 x^(m+1)) + Σ_{k≥1} B_{2k}(2k+m-1)!/(2k)! x^-(2k+m) ]`. E.g. `Series[PolyGamma[1, x], {x, Infinity, 5}]` → `1/x + 1/(2 x^2) + 1/(6 x^3) - 1/(30 x^5) + O[1/x]^6`. `PolyGamma[x]` (unindexed) normalises to `PolyGamma[0, x]`. A requested order below the leading power `m` yields a pure `O[1/x]^(n+1)`.
- **Expansion at `-Infinity`** is handled by reflection: `Series[f, {x, -Infinity, n}]` is computed as `Series[f[x -> -x], {x, +Infinity, n}]` and then reflected back `x -> -x`, so the entire `+Infinity` machinery (the special-function asymptotics above and the generic `x -> 1/u` path) is reused verbatim. The result is reported in the canonical `1/x` variable: `Series[Sqrt[x^2 + a x], {x, -Infinity, 1}, Assumptions -> Element[a, Reals]]` = `-x - a/2 + a^2/(8 x) + O[1/x]^2` (the `+Infinity` twin with the odd-exponent terms sign-flipped), and `Series[Sqrt[a^2 x^2 + b x], {x, -Infinity, 1}, Assumptions -> a > 0 && Element[b, Reals]]` = `-a x - b/(2 a) + b^2/(8 a^3 x) + O[1/x]^2`. Because `-Infinity` never reaches the finite-`x0` arithmetic, this path emits none of the spurious `Infinity::indet` warnings the pre-reflection code did.
- **Parametric (scaled) arguments at infinity.** For a head applied to a homogeneous linear argument `c*x` (with `c` free of `x` and its sign known under the assumptions), the sign of `c` selects the branch/direction:
  - *`ArcTan[c x]`* emits `Sign[c] Pi/2 - 1/(c x) + ...`: `Series[ArcTan[a x], {x, Infinity, 1}, Assumptions -> a > 0]` = `Pi/2 - 1/(a x)`, and `-Pi/2 - 1/(a x)` for `a < 0` (the argument runs to `-Infinity`). Composes with the `-Infinity` reflection, so `{x, -Infinity, 1}` flips accordingly.
  - *Saturating heads* `Erf`/`Erfc`/`Tanh`/`FresnelS`/`FresnelC` limit to a sign-dependent constant: `Erf[c x] -> Sign[c]`, `Tanh[c x] -> Sign[c]`, `Erfc[c x] -> 0`/`2`, `FresnelS[c x]`/`FresnelC[c x] -> ±1/2` (at order 0). Emitted directly for the scaled form (the bare `H[x]` emitter, which keeps the exponential/oscillatory tail, is unchanged).
  - *`ArcCosh[c x]`/`ArcSinh[c x]`* take the correct principal branch for `c < 0` via the reflections `ArcCosh[c x] = I Pi + ArcCosh[-c x]` and `ArcSinh[c x] = -ArcSinh[-c x]`: `Series[ArcCosh[a x], {x, Infinity, 1}, Assumptions -> a < 0]` = `I Pi + Log[-2 a x] + ...` (previously the wrong `-I Pi` branch).
  - *Growth/decay heads* `ExpIntegralEi[c x]` and `ExpIntegralE[1, c x]` (`c > 0`) carry the scaled asymptotic: `Series[ExpIntegralEi[a x], {x, Infinity, 1}, Assumptions -> a > 0]` = `E^(a x)/(a x) + O[...]`, and `ExpIntegralE[1, a x]` = `E^(-a x)/(a x) + O[...]`.
- `Series[ExpIntegralE[1, x], {x, Infinity, n}]` returns the first-kind exponential-integral asymptotic (DLMF 6.12.1): `E_1(x) ~ E^(-x) Σ_{k≥0} (-1)^k k!/x^(k+1)` = `E^(-x)(1/x - 1/x^2 + 2/x^3 - ...)`, the decaying `E^(-x)` kept as a symbolic prefactor (mirroring `ExpIntegralEi`).
- `Series[Zeta[x], {x, Infinity, n}]` returns the truncated Dirichlet head `1 + 2^-x + 3^-x + ... + (n+1)^-x`. Because `Zeta[x] = Σ_{k≥1} k^-x` and each `k^-x = Exp[-x Log[k]]` is exponentially smaller than the previous as `x → +∞`, the natural asymptotic scale here is *exponential* (`2^-x`), not a power of `1/x` — so this is a plain `Plus` of powers, not a `SeriesData` with an `O`-term. `n` counts the correction terms kept (bases `2..n+1`); the leading `1` is always present. This is exactly the exp-log form the Gruntz mrv engine consumes to resolve limits like `(Zeta[x] - 1) 2^x → 1` and `Log[Zeta[x] - 1]/x → -Log[2]`.

**Logarithmic series at `x = 0` for special functions**: Some functions have a finite-radius series at the origin built around a single `Log[x]` branch term, with the log baked into the `x^0` coefficient. Naive Taylor-via-`D` cannot reach these (`f(0)` is infinite or the derivatives have poles), so they are emitted from dedicated identities:
- `Series[ExpIntegralEi[x], {x, 0, n}]` returns `EulerGamma + Log[x] + x + x^2/4 + x^3/18 + ... + O[x]^(n+1)`, i.e. `SeriesData[x, 0, {EulerGamma + Log[x], 1, 1/4, ..., 1/(n*n!)}, 0, n+1, 1]` (DLMF 6.6.2: `Ei(x) = EulerGamma + Log[x] + Σ_{k≥1} x^k/(k k!)`). The `EulerGamma + Log[x]` branch term occupies the `x^0` slot. With `Assumptions -> x < 0` the branch term becomes `EulerGamma + Log[-x]`.

**Generalized series in `Log[x]` at `x = 0`**: `LogIntegral[x] = Ei(Log[x])` has *no* ordinary Taylor (or Laurent/Puiseux) series at `x = 0`: as `x -> 0+`, `L = Log[x] -> -Infinity` drives it into Ei's asymptotic regime, so the result is a generalized series whose `x^1` coefficient is itself a series in `1/Log[x]`.
- `Series[LogIntegral[x], {x, 0, n}]` returns `x Σ_{k=0}^{2n+1} k!/Log[x]^(k+1) + O[x]^(n+1)`. For `n = 2` this is `((120 + 24 Log[x] + 6 Log[x]^2 + 2 Log[x]^3 + Log[x]^4 + Log[x]^5) x)/Log[x]^6 + O[x]^3`, matching the reference's `Assumptions -> x > 0` output. Since the only `x`-dependence is the prefactor `E^L = x`, every term carries exactly `x^1`; the coefficient slot is emitted as `Together[Σ k!/Log[x]^(k+1)]`. The reference's no-assumptions form additionally wraps this in a `Floor[Arg[...]]` branch discriminator (to track `Log[1/x]` vs `Log[x]` across the cut); Mathilda emits the principal `x > 0` form by default, or — with `Assumptions -> x < 0` — the `x < 0` form, where every `Log[x]` becomes `Log[-x]` and an additive `I Pi` (the `x^0` term, from `Log[x] = Log[-x] + I Pi`) leads the series: `I Pi + ((120 + 24 Log[-x] + ... + Log[-x]^5) x)/Log[-x]^6 + O[x]^3`.

**Internal padding for symbolic expansion points**: The engine computes series at a padded internal order (user order + 12 by default) so that intermediate Laurent/Puiseux operations don't lose accuracy. When the expansion point `x0` is not a literal number, padding is capped at 2 — at a symbolic point the series coefficients are themselves symbolic expressions (e.g. `Cosh[a]`, `Sinh[a]`), and the `O(N^2)` convolution inside `so_inv`/`so_div` would otherwise spin indefinitely on exponentially growing expression trees. This makes cases like `Series[Coth[x], {x, a, 1}]`, `Series[Tanh[x], {x, a, 1}]`, `Series[Sec[x], {x, a, 1}]`, and `Series[1/Cosh[x], {x, a, 1}]` terminate in milliseconds.

**Requested order is guaranteed at a numeric expansion point** (v0.6088).
- **Why it was needed.** The pad of 12 was a guess. Every order of pole that flows through `so_inv`/`so_mul` costs
  one order of relative precision: an `Apart`-merged denominator such as `e^7 + e^8`, or `1/e^k` squared up by
  `so_pow_int`. So a pole of order ≥ 7 used up the pad, and the expansion came back short of the requested O-term.
- **What went wrong.** `Series[Sqrt[(1+e)^2+e^6]/(e^7 (1+e)), {e,0,0}]` gave `1/e^7 + O[e]^-1`, and
  `SeriesCoefficient` read the missing `e^-1` coefficient as 0.
- **The fix.** `do_series_single` now measures any shortfall and re-expands with the pad widened by exactly that
  much, for at most 4 rounds. The symbolic-centre pad stays tight on purpose; `residue.c` compensates there.
- **A related fix.** The single-term branch of `so_inv` no longer overstates the precision of a truncated monomial:
  its O-term moves to `order - 2 nmin`.

**Constant inputs**: If `f` is free of the expansion variable (e.g. `Series[0, {x, 0, 4}]`, `Series[Sin[y], {x, 0, 4}]`, `Series[a + b^2, {x, 0, 3}]`), `Series` returns `f` verbatim instead of wrapping it in a trivial `SeriesData`.

**Symbolic prefactors**: A factor of `x^alpha` with `alpha` symbolic (non-integer, non-rational) is pulled outside the expansion so the remaining body is expanded as an ordinary power series. For example, `Series[x^a Exp[x], {x, 0, 5}]` returns `x^a (1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + O[x]^6)` — a `Times[Power[x, a], SeriesData[...]]` at the expression level, so the `SeriesData` pretty-printer still renders the body and the outer `Times` decorates it with the symbolic prefactor.

**Expansion at regular points of `Arc*` heads**: When `so_apply_kernel_at_zero` can't apply (because the inner series constant `c` is not `0`), the engine falls back to naive Taylor via `D`. This makes `Series[ArcSin[x], {x, 1/2, 3}]`, `Series[ArcTan[x], {x, 2, 2}]`, `Series[ArcSinh[x], {x, 1, 2}]`, etc. work without special-casing each non-zero expansion point.

**Maxima-style algebraic fast paths**:
- **Monomial binomial** `(a + b x^m)^alpha` with `alpha` non-integer. The generic path factors out `a`, forms `u = (b/a) x^m`, and feeds `(1+u)^alpha` through Horner composition. When `u` is a single-term series (exactly one non-zero coefficient in `SeriesObj` terms) we skip the `O(N^2)` convolution and emit `binomial(alpha, k) * (b/a)^k` directly at exponent `k*m`. This covers `Sqrt[1+x]`, `(1 - x^2)^(1/2)`, `(1+x)^(-1/2)`, `(2 + 3x)^(1/3)`, `(1+x)^n` with symbolic `n`, and Puiseux bases like `(1+x^(1/2))^(1/2)`.
- **split-two-term probe** (`series_split_two_term` in `series.h`). Structural decomposition of `e` into `a + b*x^(p/q)` without running the full series-expansion pipeline. Feeds the Log fast path and the Apart gate; exposed for unit testing.
- **Log fast path**: when `arg` matches `a + b x^(p/q)` with `a, b` both free of `x` and `a != 1`, rewrite `Log[a + b x^c]` as `Log[a] + Log[1 + (b/a) x^c]` and let the `Log1p` kernel compose with a pure monomial. Maxima's `sp2log` uses the same identity.
- **Apart preprocessing**: if the input contains `Power[p(x), -n]` for `p` a polynomial in `x`, run `Apart[f, x]` to decompose into partial fractions before expanding. Composite denominators like `1/((1-x)(1-2x)(1-3x))` then break up into geometric-series pieces that hit the monomial fast path. Gated by a polynomial check so non-rational denominators (e.g. `1/(Exp[x] - 1 - x)`) bypass Apart and fall through to the generic `so_inv` path.

**Branch-point expansion for inverse trig / hyperbolic heads** (reference-faithful). All eight `Arc*` heads expand at their branch points with the same wrapped output shape:

```
Plus[ f(x0),
      Times[ log_coef, Log[x - x0] ],         (Family B only)
      Times[ (-1)^Floor[(Pi/2 - Arg[x - x0])/(2 Pi)],
             SeriesData[...] ] ]
```

The `(-1)^Floor[...]` factor is the reference branch discriminator — it is `1` on the principal sheet near `x0` and flips sign across the branch cut. Two mathematical families:

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

`Normal[Series[ArcTanh[x], {x, 1, 3}]]` preserves the `Log[x - 1]` term and the branch discriminator — `Normal` collapses the inner `SeriesData` but the surrounding `Plus`/`Times` pass through unchanged (matching the reference).

**Naive-Taylor fallback**: For expansion points that are not branch points (e.g. `Series[ArcSinh[x], {x, 1 + I, 3}]`, `Series[ArcSin[x], {x, 1/2, 3}]`), `series_expand` falls back to naive Taylor via repeated `D`. The fallback caps iterations and bails out on `Infinity` / `Indeterminate` derivatives so unknown heads cannot spin the engine. The singularity probe **evaluates** the substituted `f^(k)(x0)` with arithmetic warnings muted before testing for infinities, so an *unevaluated* pole produced by the substitution (e.g. `f[1/x]` at `x = 0` becoming `f[1/0]`) collapses to `ComplexInfinity` / `Indeterminate` and is detected cleanly — it no longer spills a spurious `Power::infy: 1/0` to stderr before the guard fires. So heads with no recognised expansion at the requested point (e.g. `Series[Gamma[x], {x, Infinity, 2}]`) return unevaluated silently rather than with a stray warning.

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
Out[7]= 1/x - 1/(6 x^3) + 1/(120 x^5) - 1/(5040 x^7) + 1/(362880 x^9) + O[1/x]^11

In[8]:= Series[Sin[x + y], {x, 0, 3}, {y, 0, 3}]
Out[8]= y - 1/6 y^3 + O[y]^4 + (1 - 1/2 y^2 + O[y]^4) x + (-1/2 y + 1/12 y^3 + O[y]^4) x^2 + (-1/6 + 1/12 y^2 + O[y]^4) x^3 + O[x]^4

In[9]:= Series[{Sin[x], Cos[x], Tan[x]}, {x, 0, 5}]
Out[9]= {x - 1/6 x^3 + 1/120 x^5 + O[x]^6, 1 - 1/2 x^2 + 1/24 x^4 + O[x]^6, x + 1/3 x^3 + 2/15 x^5 + O[x]^6}
```

## Normal
Converts a `SeriesData` back into an ordinary expression by dropping its O-term.
- `Normal[expr]`

**Features**:
- `Protected`.
- Returns the Plus of the coefficient-times-power terms (zero coefficients skipped). For non-`SeriesData` input, `Normal` is the identity.
- `Normal[SparseArray[...]]` gives the dense nested List the specification
  denotes (see [SparseArray](data-structures.md#sparsearray)), at any depth.
- Recurses through the whole expression, dropping the O-term of **every** `SeriesData` at any depth. This matters for expansions around `+-Infinity`, whose `SeriesData` is wrapped inside `Plus`/`Times` (e.g. the trig- or exponential-prefactored asymptotic forms of `BesselJ`, `BesselY`, `BesselK`, `BesselI`, `AiryAi`, `AiryBiPrime`); the surrounding factors are preserved and recombined by the evaluator.

```mathematica
In[1]:= Normal[Series[Exp[x], {x, 0, 5}]]
Out[1]= 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120

In[2]:= Normal[a + b]
Out[2]= a + b

In[3]:= Normal[Series[BesselJ[0, x], {x, Infinity, 2}]]
Out[3]= Sqrt[2/Pi] Sqrt[1/x] Cos[1/4 Pi - x] - 1/8 Sqrt[2/Pi] (1/x)^(3/2) Sin[1/4 Pi - x]

In[4]:= Normal[SparseArray[{1 -> a, 3 -> b}, 4]]
Out[4]= {a, 0, b, 0}
```

## SeriesCoefficient
Gives the coefficient of `(x - x0)^k` in the power-series expansion of `f`.
- `SeriesCoefficient[f, {x, x0, k}]`

**Features**:
- `Protected` only; the arguments are evaluated and the expansion variable must be
  a symbol (see `Series` above — the same `ivar` decline applies, under the head
  `SeriesCoefficient`).
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

