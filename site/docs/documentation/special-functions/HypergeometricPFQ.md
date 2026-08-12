# HypergeometricPFQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergeometricPFQ[{a1, ...}, {b1, ...}, z]`**

is the generalized hypergeometric function pFq(a;b;z), the series

<details>
<summary>Notes</summary>

Sum (prod\_i Pochhammer\[a\_i, k\] / prod\_j Pochhammer\[b\_j, k\]) z^k / k!. Common upper/lower parameters cancel; a non-positive integer upper parameter terminates the series to a polynomial. Evaluates to machine, arbitrary-precision (MPFR), and complex numbers by direct summation in the convergent regime (p\<=q for all z; p==q+1 for |z|\<1).

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= HypergeometricPFQ[{a1, a2, a3}, {b1, b2, b3}, 0]
Out[1]= 1

In[2]:= HypergeometricPFQ[{}, {}, z]
Out[2]= E^z

In[3]:= HypergeometricPFQ[{a, b, c}, {a, d, e}, z]
Out[3]= HypergeometricPFQ[{b, c}, {d, e}, z]

In[4]:= HypergeometricPFQ[{1, 1}, {3, 3, 3}, 2.]
Out[4]= 1.07893

In[5]:= HypergeometricPFQ[{1, 2, 3, 4}, {5, 6, 7}, {0.1, 0.3, 0.5}]
Out[5]= {1.01164, 1.03627, 1.06296}

In[6]:= N[HypergeometricPFQ[{1, 1, 1}, {3/2, 3/2, 3/2}, 10], 50]
Out[6]= 530.19188827362590438855961685444087792733053398358

In[7]:= D[HypergeometricPFQ[{a1, a2}, {b1, b2, b3}, x], x]
Out[7]= (a1 a2 HypergeometricPFQ[{1 + a1, 1 + a2}, {1 + b1, 1 + b2, 1 + b3}, x])/(b1 b2 b3)
```

### Applications (3)

```mathematica
In[1]:= HypergeometricPFQ[{1, 1}, {2}, z]
Out[1]= -Log[1 - z]/z
```

Common upper and lower parameters cancel, and any of the specialised forms
(`0F1`, `1F1`, `2F1`) is just a particular shape of `HypergeometricPFQ`. A
higher `3F2` evaluates numerically by direct summation:

```mathematica
In[1]:= N[HypergeometricPFQ[{1, 2, 3}, {4, 5}, 1/2], 40]
Out[1]= 1.1898747542564229318256831180919799547257
```

A non-positive integer upper parameter terminates the series to a polynomial:

```mathematica
In[1]:= HypergeometricPFQ[{-3, 1}, {1}, z]
Out[1]= 1 - 3 z + 3 z^2 - z^3
```

## Options & behaviour

> **Packed arrays.** With a real buffer as `z` and every parameter inexact,
> the series is summed elementwise over the buffer. `Hypergeometric0F1`,
> `1F1` and `2F1` all rewrite to `HypergeometricPFQ`, so they take the same
> path. **Exact** parameters do not: the parameter-cancellation, terminating
> and closed-form reductions answer with exact or symbolic results
> (`Hypergeometric2F1[1, 1, 2, z]` is `-Log[1 - z]/z`) and must keep running.

## Implementation notes

- Attributes `NumericFunction`, `Protected`.
- `HypergeometricPFQ[a, b, 0]` is `1`; threads over a `List` third argument.
- Parameters common to the upper and lower lists cancel; a non-positive
  integer upper parameter terminates the series to a polynomial (valid for
  symbolic `z`).
- Reduces to elementary functions for simple parameters: `0F0 -> E^z`,
  `1F0(a) -> (1-z)^(-a)`, `0F1(1/2) -> Cosh[2 Sqrt[z]]`,
  `0F1(3/2) -> Sinh[2 Sqrt[z]]/(2 Sqrt[z])`, `1F1(1;2) -> (E^z-1)/z`,
  `2F1(1,1;2) -> -Log[1-z]/z`. The central-binomial / arcsin family also closes:
  `2F1(1,1;1/2;z) -> 1/(1-z) + Sqrt[z] ArcSin[Sqrt[z]]/(1-z)^(3/2)`,
  `2F1(1,1;3/2;z) -> ArcSin[Sqrt[z]]/(Sqrt[z] Sqrt[1-z])`,
  `2F1(2,1;3/2;z) -> 1/(2(1-z)) + ArcSin[Sqrt[z]]/(2 Sqrt[z] (1-z)^(3/2))`.
  A table-backed very-well-poised reduction closes the classic
  central-binomial-cubed Ramanujan `1/Pi` series:
  `4F3({1/2,1/2,1/2,5/4}, {1/4,1,1}, -1) -> 2/Pi` (general `1/Pi` summation is an
  open problem, so this is a recognizer for the known class, not universal).
- Numeric evaluation at machine, arbitrary (MPFR), and complex precision by
  direct series summation, with output precision tracking the input. The
  series is summed only where it converges — `p <= q` (entire) and `p == q+1`
  with `|z| < 1`; outside that regime (and `p > q+1`) the call stays
  unevaluated (analytic continuation beyond the unit disk is not yet
  implemented).
- **Cancellation is handled, not ignored.** For negative real `z` the series
  alternates and its largest term exceeds the sum by a factor of order `e^|z|`,
  so a doubles-only summation silently loses that many bits. The machine path
  measures the loss (`max|term| / |sum|`) and re-sums through the MPFR path at
  enough working precision to absorb it, then rounds back. The answer stays
  machine precision — it is simply the correctly rounded one. `1F1(1;2;-40)`
  agrees with its closed form `(E^z-1)/z` to one ulp; before this it was wrong
  in the second decimal place.
- Has a machine kernel, so bodies containing it compile
  (`Compile`, `Plot`, `NIntegrate`, `NDSolve`, …) rather than falling back to
  the interpreter. `Hypergeometric0F1`/`1F1`/`2F1` share that one kernel.
- `D[HypergeometricPFQ[{a},{b},x], x]
   = (prod a_i / prod b_j) HypergeometricPFQ[{a_i+1},{b_j+1},x]`.

**Attributes:** `NumericFunction`, `Protected`.

## See also

[List](../../other-advanced/List/), [Compile](../../control-flow/Compile/), [Plot](../../graphics/Plot/), [NIntegrate](../../numerical-calculus/NIntegrate/), [NDSolve](../../numerical-calculus/NDSolve/), [Hypergeometric0F1](../../special-functions/Hypergeometric0F1/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergeopfq.c)

## Notes & additional examples

### Notes

`HypergeometricPFQ[{a1, ...}, {b1, ...}, z]` is the generalized hypergeometric function `pFq`, the series `Sum[(Product[Pochhammer[a_i, k]] / Product[Pochhammer[b_j, k]]) z^k / k!]`. It converges for all `z` when `p <= q`, and for `|z| < 1` when `p == q + 1`; a non-positive integer upper parameter truncates it to a polynomial. The specialised heads `Hypergeometric0F1`, `Hypergeometric1F1`, and `Hypergeometric2F1` are convenience wrappers.
