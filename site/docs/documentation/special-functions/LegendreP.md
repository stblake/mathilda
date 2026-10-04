# LegendreP

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LegendreP[n, x]`**

gives the Legendre polynomial P\_n(x).

**`LegendreP[n, m, x] gives the associated Legendre function P_n^m(x).`**

**`LegendreP[n, m, a, x] gives the Legendre function of type a (a in`**

<details>
<summary>Notes</summary>

{1, 2, 3}, default 1). Integer n yields the explicit polynomial; a non-integer order with an inexact argument evaluates numerically at machine or arbitrary (MPFR) precision, real or complex. Listable.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= LegendreP[3, x]
Out[1]= -3/2 x + 5/2 x^3

In[2]:= LegendreP[10, 2, x]
Out[2]= (1 - x^2) (3465/128 - 45045/32 x^2 + 675675/64 x^4 - 765765/32 x^6 + 2078505/128 x^8)
```

### Applications (8)

```mathematica
In[3]:= LegendreP[2, x]
Out[3]= -1/2 + 3/2 x^2

In[4]:= Table[LegendreP[n, x], {n, 0, 4}]
Out[4]= {1, x, -1/2 + 3/2 x^2, -3/2 x + 5/2 x^3, 3/8 - 15/4 x^2 + 35/8 x^4}
```

An even-degree value at the origin

```mathematica
In[5]:= LegendreP[4, 0]
Out[5]= 3/8
```

Every P_n equals 1 at the endpoint

```mathematica
In[6]:= LegendreP[n, 1]
Out[6]= 1
```

The associated function P_2^1

```mathematica
In[7]:= LegendreP[2, 1, x]
Out[7]= -3 x Sqrt[1 - x^2]
```

Orthogonality on the interval

```mathematica
In[8]:= Integrate[LegendreP[2, x] LegendreP[3, x], {x, -1, 1}]
Out[8]= 0
```

Non-integer order via the Gauss 2F1 series

```mathematica
In[9]:= N[LegendreP[1.5, 0.3]]
Out[9]= -0.0897873
```

```mathematica
In[10]:= D[LegendreP[3, x], x]
Out[10]= -3/2 + 15/2 x^2
```

## Algorithm

Mathilda -- Legendre polynomials and associated Legendre functions.

```text
  LegendreP[n, x]        Legendre polynomial / function P_n(x).
  LegendreP[n, m, x]     associated Legendre function P_n^m(x) (type 1).
  LegendreP[n, m, a, x]  Legendre function of type a (a in {1, 2, 3}).
```

Evaluation is layered so each argument shape takes the cheapest exact or numeric route:

```text
  LegendreP[n, x]
    exact integer n            ->  the explicit degree-|n'| polynomial in x
                                   (n' = n, or -1-n for n < 0, since
                                   P_{-1-n} = P_n) with exact rational
                                   coefficients, built from the three-term
                                   recurrence; an inexact x then evaluates
                                   the monomials numerically.
    x == 1                     ->  1 (for any order n).
    non-integer n, some arg    ->  numeric Gauss series
       inexact                      P_n(x) = 2F1(-n, n+1; 1; (1-x)/2),
                                   real or complex, machine or MPFR
                                   precision (requires |(1-x)/2| < 1).
    everything else            ->  stays symbolic (return NULL).

  LegendreP[n, m, x] / [n, m, a, x]   (integer n, integer m >= 0)
    type 1 (default, a == 1)   ->  (-1)^m (1-x^2)^(m/2) d^m/dx^m P_n(x)
                                   (the Rodrigues derivative form; 0 when
                                   m > |n'|).
    types 2, 3                 ->  C(x) * R_a(x), where
                                   C(x) = 2F1Reg(-n, n+1, 1-m, (1-x)/2)
                                   is the (terminating, exact) regularized
                                   Gauss polynomial and the prefactor is
                                     R_2 = (1+x)^(m/2) (1-x)^(-m/2),
                                     R_3 = (1+x)^(m/2) (-1+x)^(-m/2).
    non-integer / negative m   ->  stays symbolic (return NULL).
```

Attributes: Listable, NumericFunction, Protected.

Deferred (left symbolic): symbolic Series / SeriesCoefficient, D[] rules, the non-integer associated and Legendre-function forms, and analytic continuation of the numeric series for |(1-x)/2| >= 1.

## Implementation notes

**Algorithm.** `builtin_legendre_p` handles `LegendreP[n, x]`,
`LegendreP[n, m, x]` and `LegendreP[n, m, a, x]`. For `P_n(x)` an exact integer
order builds the explicit degree-`|n'|` polynomial from the three-term
recurrence `k P_k = (2k-1) x P_{k-1} - (k-1) P_{k-2}` with exact `mpq`
coefficients (using `P_{-1-n} = P_n`, cap `LEG_POLY_CAP = 2000`); `x == 1 -> 1`
for any order. A non-integer order with an inexact argument is evaluated by the
Gauss series `P_n(x) = 2F1(-n, n+1; 1; (1-x)/2)` summed in the `ncpx` MPFR-complex
toolkit (real/complex, machine/arbitrary precision; requires `|(1-x)/2| < 1`).
The associated forms (integer `n`, integer `m >= 0`): type 1 is the Rodrigues
derivative `(-1)^m (1-x^2)^(m/2) d^m/dx^m P_n(x)` (0 when `m > |n'|`); types 2
and 3 are the regularized Gauss polynomial `2F1Reg(-n, n+1, 1-m, (1-x)/2)` times
a `(1±x)^(±m/2)` prefactor. Non-integer/negative `m`, and `|(1-x)/2| >= 1`, stay
symbolic.

**Data structures.** `Expr`; GMP `mpq_t` coefficient arrays for the polynomial;
`ncpx` (`mpfr_t` re/im) for the numeric series. ND: binary kernel
`NDKB_LegendreP = { ndk_LegendreP_c, ... }` (complex, threads `n` and `x`
element-wise), registered `REG_B`, so `packed_aware`. Attributes: `Listable`,
`NumericFunction`, `Protected`.

**Complexity / limits.** The recurrence is `O(n^2)` in growing big rationals,
capped at `n = 2000`; the numeric series costs `O(wp)` terms inside the
convergence disk. Deferred (left symbolic): symbolic `Series`/`D` rules, the
non-integer associated forms, and continuation for `|(1-x)/2| >= 1`. `Compile[]`
lowers at both scalar and rank-1 array shapes (`Compiled -> True`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/), [Series](../../power-series/Series/), [SeriesCoefficient](../../power-series/SeriesCoefficient/)

- DLMF §14 — Legendre and related functions (§14.3 series, §14.7 integer degree).
- Source: [`src/special_functions/legendre.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/legendre.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve_m14_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m14_stress.c)
- Tests: [`tests/test_legendre.c`](https://github.com/stblake/mathilda/blob/main/tests/test_legendre.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

For an exact integer order `LegendreP[n, x]` is the explicit Legendre polynomial
`P_n(x)`, built from the three-term recurrence with exact rational coefficients
(orders up to `n = 2000`); `LegendreP[n, 1]` is `1` for every order. The
Legendre polynomials are orthogonal on `[-1, 1]`, so the integral of a product
of two of different degree vanishes.

A non-integer order is evaluated numerically (when an argument is inexact)
through the Gauss hypergeometric series `P_n(x) = 2F1(-n, n+1; 1; (1-x)/2)`,
valid for `|(1-x)/2| < 1`.

`LegendreP[n, m, x]` is the associated Legendre function, and
`LegendreP[n, m, a, x]` selects the Legendre function of type `a` (one of 1, 2,
3). The derivatives of the integer-order polynomials follow from the explicit
form.
