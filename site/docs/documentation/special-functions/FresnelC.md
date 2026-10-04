# FresnelC

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FresnelC[z]`**

gives the Fresnel integral C(z) = Integral\_0^z Cos\[Pi t^2/2\] dt.

**`FresnelC[+-Infinity] = +-1/2, FresnelC[+-I Infinity] = +-I/2.`**

<details>
<summary>Notes</summary>

An entire, odd function with no branch cuts. FresnelC\[0\] = 0, Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[FresnelC\[z\], z\] = Cos\[Pi z^2/2\]. Listable.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

Odd and entire, so the value at zero is zero

```mathematica
In[1]:= FresnelC[0]
Out[1]= 0
```

The limit value is 1/2

```mathematica
In[2]:= FresnelC[Infinity]
Out[2]= 1/2
```

Arbitrary precision via the MPFR series

```mathematica
In[3]:= N[FresnelC[1], 20]
Out[3]= 0.779893400376822829476
```

The integrand Cos[Pi x^2/2]

```mathematica
In[4]:= D[FresnelC[x], x]
Out[4]= Cos[1/2 Pi x^2]
```

Threads element-wise over the packed real list

```mathematica
In[5]:= FresnelC[{0.5, 1.0, 1.5}]
Out[5]= {0.492344, 0.779893, 0.445261}
```

## Algorithm

Mathilda -- the Fresnel integrals (Pi/2-normalized).

```text
  FresnelC[z] = Int_0^z Cos[Pi t^2 / 2] dt
  FresnelS[z] = Int_0^z Sin[Pi t^2 / 2] dt
```

Both are entire and odd, with no branch cuts. FresnelC and FresnelS share one numeric kernel: the pair (C, S) is computed together and each builtin returns its component. Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values     ->  0, +-1/2, +-I/2, Indeterminate
  machine / arbitrary real ->  MPFR series (small |x|) or asymptotic (large)
  complex (any precision)  ->  the paired A/B ncpx series with guard bits
  everything else          ->  stays symbolic (return NULL)
```

Convergent Maclaurin series (valid for all z):

```text
  C(z) = Sum_{m>=0} (-1)^m (Pi/2)^(2m)   z^(4m+1) / ((2m)!   (4m+1)),
  S(z) = Sum_{m>=0} (-1)^m (Pi/2)^(2m+1) z^(4m+3) / ((2m+1)! (4m+3)).
```

Equivalently A(z) = C(z) + i S(z) = Sum_{k>=0} (i Pi/2)^k z^(2k+1)/(k!(2k+1)) and B(z) = C(z) - i S(z) is the same with i -> -i; then C = (A+B)/2 and S = (A-B)/(2i). The real path sums C and S as two real series directly; the complex path sums A and B together in one ncpx loop.

The partial sums reach magnitude ~e^((Pi/2)|z|^2) before the O(1)-sized answer emerges, so the MPFR paths add ~(Pi/2)|z|^2/ln2 guard bits to absorb that cancellation exactly. For large real |x| the convergent series is infeasible; there we use the asymptotic expansion (DLMF 7.12)

```text
  C(x) = 1/2 + f(x) sin(Pi x^2/2) - g(x) cos(Pi x^2/2),
  S(x) = 1/2 - f(x) cos(Pi x^2/2) - g(x) sin(Pi x^2/2),
    f(x) ~ (1/(Pi x))   Sum_j (-1)^j (4j-1)!! / (Pi x^2)^(2j),
    g(x) ~ (1/(Pi^2 x^3)) Sum_j (-1)^j (4j+1)!! / (Pi x^2)^(2j),
```

summed to the smallest term (optimal truncation). The asymptotic constant 1/2 is the value only in a sector around the real axis (Stokes phenomenon), so it is used for real inputs only; complex inputs always use the convergent series (correct everywhere, merely costlier for large |z|).

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_fresnelc` handles `FresnelC[z] = Int_0^z cos(π t^2/2) dt` (the π/2-normalized / Wolfram convention), entire and odd. It shares one numeric kernel with `FresnelS`: the pair `(C, S)` is computed together and this builtin returns the `C` component. Exact special values first: `FresnelC[0] = 0`, `FresnelC[±Infinity] = ±1/2`, `FresnelC[±I Infinity] = ±I/2`, `ComplexInfinity`/`Indeterminate -> Indeterminate`. A **numeric real** argument (machine or arbitrary precision) uses the convergent Maclaurin series for small/moderate `|x|` (with `~(π/2)|x|^2/ln2` guard bits to absorb the `~e^{(π/2)|z|^2}` partial-sum cancellation), or the **asymptotic expansion** DLMF 7.12 (`C(x) = 1/2 + f(x) sin(π x^2/2) - g(x) cos(π x^2/2)`, summed to optimal truncation) for large `|x|` — the asymptotic constant `1/2` holds only in a sector around the real axis (Stokes), so it is used for real inputs only. A **complex** argument always uses the convergent paired `A/B` series (`A = C + iS`, `B = C - iS`) in the shared `ncpx` toolkit (correct everywhere). A symbolic negative-leading argument folds by oddness. A `USE_MPFR=0` build uses a double-complex `A/B` series.

**Data structures.** `Expr`; the shared complex-MPFR toolkit `ncpx` (`numeric_complex.h`), folding by oddness to `Re z >= 0`. The ND kernel is a real `REG_U` registration (`NDKU_FresnelC`, `ndk_FresnelC_r` → `sf_machine_fresnel_c` in `src/special_functions/sf_machine.c`): element-wise over a packed or visible real `NDArray`. `Compile[]` lowers `FresnelC` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes.

**Complexity / limits.** `O(1)` per element at machine precision; MPFR term count and guard bits scale with `|z|^2` and precision. Entire function (no branch cuts). `D[FresnelC[x], x] = Cos[π x^2/2]`. Symbolic arguments stay symbolic. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [FresnelS](../../special-functions/FresnelS/), [Piecewise](../../control-flow/Piecewise/)

- DLMF §7.2(iii) — the Fresnel integral C(z) = Int_0^z cos(π t^2/2) dt.
- DLMF §7.12 — the asymptotic expansion of the Fresnel integrals.
- Source: [`src/special_functions/fresnel.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/fresnel.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_fresnelc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fresnelc.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

Mathilda uses the Pi/2-normalized (Wolfram) convention
`FresnelC[z] = Int_0^z Cos[Pi t^2/2] dt`. The function is entire and odd, with
`FresnelC[±Infinity] = ±1/2` and `FresnelC[±I Infinity] = ±I/2`.

`FresnelC` and `FresnelS` share one numeric kernel: the pair `(C, S)` is computed
together and each builtin returns its component. The real path uses a convergent
Maclaurin series for small/moderate arguments and an asymptotic expansion (DLMF
7.12) for large ones; complex arguments always use the convergent paired `A/B`
series, correct across the whole plane. `FresnelC` carries a real `NDArray`
kernel and lowers under `Compile[]` at both scalar and rank-1 array shapes.
