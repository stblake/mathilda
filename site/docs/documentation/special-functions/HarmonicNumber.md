# HarmonicNumber

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HarmonicNumber[n]`**

gives the n-th harmonic number H\_n = Sum\_{i=1}^n 1/i.

**`HarmonicNumber[n, r]`**

gives the order-r harmonic number H\_n^(r) = Sum\_{i=1}^n 1/i^r.

**`Zeta[r]; a non-positive integer order r gives the Faulhaber polynomial in n.`**

<details>
<summary>Notes</summary>

Non-negative integer n expands to the exact finite sum (a rational for integer r, an explicit sum for symbolic r); HarmonicNumber\[Infinity, r\] is Inexact arguments evaluate numerically at machine or arbitrary (MPFR) precision, including complex order, via Zeta\[r\] - Zeta\[r, n+1\] (and the digamma form for r = 1). Listable.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Worked examples (1)

```mathematica
In[1]:= HarmonicNumber[Infinity]
Out[1]= ComplexInfinity
```

### Applications (6)

The exact fourth harmonic number 25/12

```mathematica
In[2]:= HarmonicNumber[4]
Out[2]= 25/12
```

Generalized order-2 harmonic number, an exact rational

```mathematica
In[3]:= HarmonicNumber[10, 2]
Out[3]= 1968329/1270080
```

The limit is Zeta[2] = Pi^2/6

```mathematica
In[4]:= HarmonicNumber[Infinity, 2]
Out[4]= 1/6 Pi^2
```

Numericalised via EulerGamma + PolyGamma[0, n+1]

```mathematica
In[5]:= N[HarmonicNumber[10], 20]
Out[5]= 2.92896825396825396826
```

Non-positive order gives the Faulhaber polynomial in n

```mathematica
In[6]:= HarmonicNumber[n, -2]
Out[6]= 1/6 n + 1/2 n^2 + 1/3 n^3
```

Threads over the list

```mathematica
In[7]:= HarmonicNumber[Range[4]]
Out[7]= {1, 3/2, 11/6, 25/12}
```

## Algorithm

Mathilda -- HarmonicNumber: generalized (order-r) harmonic numbers.

```text
  HarmonicNumber[n]     H_n     = Sum_{i=1}^n 1/i
  HarmonicNumber[n, r]  H_n^(r) = Sum_{i=1}^n 1/i^r
```

Rather than carry bespoke numeric kernels, HarmonicNumber reduces to the primitives the system already provides and lets the evaluator finish the job:

```text
  n a non-negative integer  ->  explicit finite sum  Sum_{i=1}^n i^-r
                                (combines to an exact rational for integer r,
                                 stays an explicit Plus for symbolic/complex r)
  n -> Infinity             ->  Zeta[r]
  r a non-positive integer  ->  Faulhaber polynomial (built from BernoulliB)
  inexact argument          ->  N[ Zeta[r] - Zeta[r, n+1] ]   (r != 1)
                                N[ EulerGamma + PolyGamma[0, n+1] ]  (r == 1)
  otherwise                 ->  symbolic (return NULL)

The analytic identity  H_n^(r) = Zeta[r] - Zeta[r, n+1]  (and its r == 1
```

digamma special case) carries arbitrary precision and complex arguments straight through Zeta / PolyGamma.

Memory: builtin_harmonicnumber takes ownership of `res` but must not free it

```text
(the evaluator does).  Every Expr* built here is owned and either handed to
```

expr_new_function (which adopts it) or released via eval_and_free.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_harmonicnumber` serves `HarmonicNumber[n] = Sum_{i=1}^n 1/i` and the generalized `HarmonicNumber[n, r] = Sum_{i=1}^n 1/i^r`. Rather than carry bespoke numeric kernels it reduces to existing primitives and lets the evaluator finish: an exact **non-negative integer `n`** (within `HN_EXPAND_CAP = 100000`) expands to the explicit finite sum `Sum_{i=1}^n i^{-r}` (an exact rational for integer `r`, an explicit `Plus` for symbolic/complex `r`); **`n -> Infinity`** gives `Zeta[r]`; a **non-positive integer `r = -m`** gives the Faulhaber polynomial in `n` built from `BernoulliB`; and an **inexact / numericizable** argument uses the analytic identity `H_n^{(r)} = Zeta[r] - Zeta[r, n+1]` (or, for `r == 1`, `EulerGamma + PolyGamma[0, n+1]`) wrapped in `N[…]` at the input precision — so arbitrary precision and complex arguments pass straight through `Zeta`/`PolyGamma`. A `numericizable` guard keeps `HarmonicNumber[x, 2.5]` symbolic in a free symbol `x`. Everything else stays symbolic.

**Data structures.** `Expr` trees driven through `eval_and_free`; GMP for the integer/exact paths. The ND kernel is a real `REG_U` registration (`NDKU_HarmonicNumber`, `ndk_HarmonicNumber_r` → `sf_machine_harmonic`, which evaluates `γ + ψ(x+1)` via the machine digamma): element-wise over a packed or visible real `NDArray`. `Compile[]` lowers `HarmonicNumber` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes.

**Complexity / limits.** The finite-sum expansion is `O(n)` and capped at `n <= 100000` (an exact integer argument is never numerically contaminated, so beyond the cap it simply stays symbolic — there is no cheaper fallback). Numeric reduction cost is that of `Zeta`/`PolyGamma`. Diagnostics route through `mth_message` (`argt`). Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [BernoulliB](../../special-functions/BernoulliB/)

- DLMF §25.11 — the Hurwitz zeta function and the relation H_n^{(r)} = ζ(r) - ζ(r, n+1).
- DLMF §5.15 — the digamma function (H_n = γ + ψ(n+1)).
- Source: [`src/special_functions/harmonicnumber.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/harmonicnumber.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_harmonicnumber.c`](https://github.com/stblake/mathilda/blob/main/tests/test_harmonicnumber.c)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

`HarmonicNumber[n] = Sum_{i=1}^n 1/i` and `HarmonicNumber[n, r] = Sum_{i=1}^n 1/i^r`.
An exact non-negative integer `n` expands to the explicit finite sum (an exact
rational for integer `r`); `n -> Infinity` gives `Zeta[r]`; a non-positive integer
order `r = -m` gives a Faulhaber polynomial in `n` built from `BernoulliB`.

Inexact or numericizable arguments reduce through the analytic identity
`H_n^{(r)} = Zeta[r] - Zeta[r, n+1]` (and, for `r == 1`, `EulerGamma + PolyGamma[0, n+1]`),
so arbitrary precision and complex arguments pass straight through `Zeta` and
`PolyGamma`. `HarmonicNumber` carries a real `NDArray` kernel and lowers under
`Compile[]` at both scalar and rank-1 array shapes.
