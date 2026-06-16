# Series

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Series[f, {x, x0, n}]
    generates a power-series expansion of f about x = x0 up to order (x - x0)^n.
Series[f, x -> x0]
    generates the leading term of a power-series expansion of f about x = x0.
Series[f, {x, x0, nx}, {y, y0, ny}, ...]
    iteratively expands f, first in x, then in y, etc.
Series handles Taylor, Laurent (negative powers), and Puiseux (fractional powers)
    expansions, as well as logarithmic and symbolic-exponent cases such as x^x
    and (1+x)^n.
Series[f, {x, Infinity, n}] expands around x = Infinity by substituting x -> 1/u.
The result of Series is a SeriesData object; use Normal to convert it back to
    an ordinary expression by dropping the O-term.
Series is Protected and HoldAll so the expansion variable is not evaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

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
- The `Assumptions -> assm` option selects the branch of the logarithmic expansions at `x = 0` (`ExpIntegralEi`, `LogIntegral`). When `assm` forces the expansion variable negative (e.g. `x < 0`, `x < -2`, `0 > x`, or an `And[...]` containing such a relation), `Log[x]` is emitted as `Log[-x]`; otherwise the principal `x > 0` form is used. The option is matched by its LHS symbol `Assumptions` (so it is not confused with a leading-term spec `x -> x0`) and is forwarded into each inner variable for multivariate expansions.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 3.
- Joel S. Cohen, *Computer Algebra and Symbolic Computation: Mathematical Methods* (A K Peters, 2003).
- Source: [`src/calculus/series.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/series.c)
- Specification: [`docs/spec/builtins/power-series.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/power-series.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Series[Sin[x], {x, 0, 5}]
Out[1]= x - 1/6 x^3 + 1/120 x^5 + O[x]^6
```

```mathematica
In[1]:= Series[1/(1 - x), {x, 0, 4}]
Out[1]= 1 + x + x^2 + x^3 + x^4 + O[x]^5
```

```mathematica
In[1]:= Series[Log[1 + x], {x, 0, 4}]
Out[1]= x - 1/2 x^2 + 1/3 x^3 - 1/4 x^4 + O[x]^5
```

```mathematica
In[1]:= Normal[Series[Exp[x], {x, 0, 3}]]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3
```

### Notes

`Series[f, {x, x0, n}]` builds a power-series expansion to order `(x - x0)^n`, returning a `SeriesData` object that prints with a trailing `O[x]^(n+1)` term. It handles Taylor, Laurent (negative powers), and Puiseux (fractional powers) cases, as well as expansion around `Infinity` via the `x -> 1/u` substitution. Apply `Normal` to drop the order term and recover an ordinary polynomial, as in the `Exp` example above. `Series` is `HoldAll`, so the expansion variable is held unevaluated while the expansion point and order are read off.
