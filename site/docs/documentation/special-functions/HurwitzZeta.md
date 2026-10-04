# HurwitzZeta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HurwitzZeta[s, a]`**

is the Hurwitz zeta function zeta(s, a) = Sum\_{k\>=0} (k + a)^-s.

<details>
<summary>Notes</summary>

Identical to Zeta\[s, a\] for Re(a) \> 0, but built on the principal-branch power (k + a)^-s, so it differs from Zeta for non-positive real a and has poles at a = 0, -1, -2, ... . HurwitzZeta\[s, 1\] is Zeta\[s\], HurwitzZeta\[s, 1/2\] is (2^s - 1) Zeta\[s\], and a positive integer a reduces to Zeta\[s\] minus a finite power sum. A non-positive integer a gives ComplexInfinity for positive integer s and the Bernoulli-polynomial value for non-positive integer s. Real, complex, machine and arbitrary-precision (MPFR) arguments evaluate numerically via an Euler-Maclaurin kernel. Listable.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= HurwitzZeta[s, 1/2]
Out[1]= (-1 + 2^s) Zeta[s]

In[2]:= HurwitzZeta[3, -3.5]
Out[2]= 0.0307784
```

### Applications (6)

Reduces to Zeta[2] when the second argument is 1

```mathematica
In[3]:= HurwitzZeta[2, 1]
Out[3]= 1/6 Pi^2
```

The half-integer shift: (2^s - 1) Zeta[s]

```mathematica
In[4]:= HurwitzZeta[s, 1/2]
Out[4]= (-1 + 2^s) Zeta[s]
```

Read off the Bernoulli polynomial, -BernoulliB[1, 0]

```mathematica
In[5]:= HurwitzZeta[0, 0]
Out[5]= 1/2
```

Positive integer second argument: Zeta[4] minus a finite power sum

```mathematica
In[6]:= HurwitzZeta[4, 5]
Out[6]= -22369/20736 + 1/90 Pi^4
```

Arbitrary precision via Euler-Maclaurin

```mathematica
In[7]:= N[HurwitzZeta[3, 1/4], 20]
Out[7]= 64.6638699687684601666
```

Inexact argument triggers the numeric kernel

```mathematica
In[8]:= HurwitzZeta[2, 1.5]
Out[8]= 0.934802
```

## Algorithm

Mathilda -- the Hurwitz zeta function.

```text
  HurwitzZeta[s, a]   zeta(s,a) = Sum_{k>=0} (k+a)^-s        (Re s > 1)
```

defined elsewhere by analytic continuation. HurwitzZeta agrees with the two-argument Zeta for Re(a) > 0, but unlike Zeta it sums the *principal branch* powers (k+a)^-s rather than ((k+a)^2)^(-s/2). The consequences:

```text
  - the two functions disagree for non-positive real a, and
  - HurwitzZeta retains the singular summands that Zeta discards, so it has
    poles at a = 0, -1, -2, ... .
```

The evaluator routes each kind of argument to the cheapest exact or fastest numeric path:

```text
  s == 1 (exact)             ->  ComplexInfinity (pole, for any a)
  a == 1                     ->  Zeta[s]          (Riemann closed forms)
  a == 1/2                   ->  (2^s - 1) Zeta[s]
  a positive integer m >= 2  ->  Zeta[s] - Sum_{k=1}^{m-1} k^-s
  a non-positive integer:
      s positive integer     ->  ComplexInfinity (pole)
      s non-positive integer ->  -BernoulliB[1-s, a]/(1-s)   (polynomial)
  any inexact operand        ->  Euler-Maclaurin complex-MPFR kernel
  everything else            ->  stays symbolic (return NULL)
```

MPFR has no Hurwitz zeta, so the numeric kernel is implemented here from the Euler-Maclaurin summation formula (DLMF 25.11.5):

```text
  zeta(s,a) = Sum_{k=0}^{N-1} (a+k)^-s
            + (a+N)^(1-s)/(s-1)
            + 1/2 (a+N)^-s
            + Sum_{j>=1} B_{2j}/(2j)! (s)_{2j-1} (a+N)^(-s-2j+1)
```

with (s)_{2j-1} the rising factorial. N grows with the working precision and

```text
|s|; the correction series is truncated at its optimal (smallest) term. The
```

kernel uses the principal branch for every (a+k)^-s, which is exactly the HurwitzZeta convention. (The structure mirrors src/special_functions/zeta.c; the self-contained Bernoulli cache and complex-MPFR toolkit are replicated here so the two files stay independent.)

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_hurwitzzeta` handles `HurwitzZeta[s, a] = Sum_{k>=0} (k+a)^{-s}`. Unlike the two-argument `Zeta`, it sums the *principal-branch* powers `(k+a)^{-s}`, so it has poles at `a = 0, -1, -2, …`. Exact reductions (when neither operand is inexact): `s == 1` gives `ComplexInfinity` (a pole for any `a`); `a == 1` gives `Zeta[s]`; `a == 1/2` gives `(2^s - 1) Zeta[s]`; a **positive integer `a = m`** gives `Zeta[s] - Sum_{k=1}^{m-1} k^{-s}` (within `HZ_HURWITZ_A_CAP = 100000`); a **non-positive integer `a`** with non-positive integer `s = -m` gives the entire value `-BernoulliB[m+1, a]/(m+1)` (and `ComplexInfinity` if `s` is a positive integer there). Any **inexact operand** routes to the complex-MPFR kernel `hz_hurwitz_cx`, which implements the Euler–Maclaurin summation formula DLMF 25.11.5 (`N` head terms, the `(a+N)^{1-s}/(s-1)` and `1/2 (a+N)^{-s}` corrections, and the Bernoulli correction series `Sum_j B_{2j}/(2j)! (s)_{2j-1} (a+N)^{-s-2j+1}` truncated at its optimal term), using the principal branch for every `(a+k)^{-s}`. A head term `k + a = 0` with `Re(s) > 0` signals a pole and returns `ComplexInfinity`. Everything else stays symbolic.

**Data structures.** `Expr`; a self-contained Bernoulli-number cache (GMP `mpq_t`) and a file-local complex-MPFR toolkit `hcx` (pairs of `mpfr_t`, alias-safe) replicating the structure of `zeta.c`. The ND kernel is a binary `REG_B` registration (`NDKB_HurwitzZeta`, `ndk_HurwitzZeta_c` → `sf_machine_hurwitz_zeta`): element-wise over real arrays, declining complex operands to the `List` path. `Compile[]` lowers `HurwitzZeta[s, a]` (e.g. `HurwitzZeta[2, x]`) at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes.

**Complexity / limits.** The Euler–Maclaurin head-term count `N` and the correction depth grow with the working precision and `|s|`; each `(a+k)^{-s}` is a complex power. Poles at non-positive integer `a` give `ComplexInfinity`. Precision follows Mathematica contagion (the minimum over inexact leaves, floored at machine). Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Zeta](../../special-functions/Zeta/)

- DLMF §25.11 — the Hurwitz zeta function ζ(s,a) = Sum_{k>=0} (k+a)^{-s} and the Euler–Maclaurin summation formula 25.11.5.
- Source: [`src/special_functions/hurwitzzeta.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/hurwitzzeta.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flint_bridge.c)
- Tests: [`tests/test_hurwitzzeta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hurwitzzeta.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

`HurwitzZeta[s, a] = Sum_{k>=0} (k+a)^{-s}` (continued to all `s != 1`). Unlike the
two-argument `Zeta`, it sums the *principal-branch* powers `(k+a)^{-s}`, so it
retains poles at `a = 0, -1, -2, …` and disagrees with `Zeta[s, a]` for
non-positive real `a`.

Exact reductions cover `s == 1` (`ComplexInfinity`), `a == 1` (`Zeta[s]`),
`a == 1/2` (`(2^s - 1) Zeta[s]`), positive integer `a`
(`Zeta[s] - Sum_{k=1}^{a-1} k^{-s}`), and non-positive integer `s` at integer `a`
(the entire value `-BernoulliB[1-s, a]/(1-s)`). Any inexact operand routes to an
Euler–Maclaurin complex-MPFR kernel (MPFR has no native Hurwitz zeta).
`HurwitzZeta` carries a real `NDArray` kernel and lowers under `Compile[]` at both
scalar and rank-1 array shapes.
