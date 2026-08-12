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

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= LegendreP[3, x]
Out[1]= -3/2 x + 5/2 x^3

In[2]:= LegendreP[10, 2, x]
Out[2]= (1 - x^2) (3465/128 - 45045/32 x^2 + 675675/64 x^4 - 765765/32 x^6 + 2078505/128 x^8)
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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[N](../../arithmetic/N/), [Series](../../power-series/Series/), [SeriesCoefficient](../../power-series/SeriesCoefficient/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_legendre.c`](https://github.com/stblake/mathilda/blob/main/tests/test_legendre.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_possiblezeroq_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_possiblezeroq_stress.c)
