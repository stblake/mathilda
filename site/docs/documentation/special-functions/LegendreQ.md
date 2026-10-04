# LegendreQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LegendreQ[n, x]`**

gives the Legendre function of the second kind Q\_n(x).

**`LegendreQ[n, m, x] gives the associated Legendre function Q_n^m(x).`**

**`LegendreQ[n, m, a, x] gives the Legendre function of type a (a in`**

<details>
<summary>Notes</summary>

{1, 2, 3}, default 1). For integer n the explicit closed form P\_n(x) (Log\[1+x\] - Log\[1-x\])/2 + v\_n(x) is generated; a non-integer order with an inexact argument on the cut (|x| \< 1) evaluates numerically at machine or arbitrary (MPFR) precision, real or complex. D\[LegendreQ\[n,x\],x\] and the origin Series are supported. Listable.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= LegendreQ[2, x]
Out[1]= -3/2 x + (1/2 Log[1 + x] - 1/2 Log[1 - x]) (-1/2 + 3/2 x^2)

In[2]:= LegendreQ[2, 0.5]
Out[2]= -0.818663
```

### Applications (6)

The second-kind solution carries a logarithm

```mathematica
In[3]:= LegendreQ[0, x]
Out[3]= 1/2 Log[1 + x] - 1/2 Log[1 - x]
```

```mathematica
In[4]:= LegendreQ[1, x]
Out[4]= -1 + x (1/2 Log[1 + x] - 1/2 Log[1 - x])

In[5]:= LegendreQ[2, x]
Out[5]= -3/2 x + (1/2 Log[1 + x] - 1/2 Log[1 - x]) (-1/2 + 3/2 x^2)
```

An exact value at the origin for non-integer order

```mathematica
In[6]:= LegendreQ[1/2, 0]
Out[6]= -(1/2 Gamma[3/4] Sqrt[Pi])/(Sqrt[2] Gamma[5/4])
```

```mathematica
In[7]:= N[LegendreQ[0, 1/2], 20]
Out[7]= 0.549306144334054845701
```

Non-integer order on the cut

```mathematica
In[8]:= N[LegendreQ[2.5, 0.4]]
Out[8]= -0.337214
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

**Algorithm.** `builtin_legendre_q` handles `LegendreQ[n, x]`,
`LegendreQ[n, m, x]` and `LegendreQ[n, m, a, x]`. Unlike `P_n`, `Q_n` carries a
logarithm even at integer order: for integer `n >= 0` (cap `LEG_POLY_CAP = 2000`,
`n < 0` singular → symbolic) it builds `Q_n(x) = P_n(x) L(x) + v_n(x)` with
`L(x) = (1/2)(Log[1+x] - Log[1-x])` and the polynomial `v_n` from the same
three-term recurrence as `P_n` but seeded `v_0 = 0`, `v_1 = -1`. A non-integer
order emits the exact special value `Q_v(0) = -(Sqrt[Pi]/2) Sin[v Pi/2]
Gamma[(v+1)/2] / Gamma[v/2+1]` for an exact zero argument (which also lets the
origin `Series` fall out of Taylor-via-`D`), and for an inexact argument on the
cut `|x| < 1` evaluates the two Frobenius `2F1` series
`Q_v(0) 2F1(-v/2, (v+1)/2; 1/2; x^2) + Q_v'(0) x 2F1((1-v)/2, (v+2)/2; 3/2; x^2)`
built as an `Expr` over `Hypergeometric2F1`/`Gamma`/`Sin`/`Cos`, so it inherits
their machine, MPFR and complex numerics. Associated `Q_n^m` (all three types,
integer `n, m >= 0`) differentiate `Q_n` in a fresh dummy variable then
substitute.

**Data structures.** `Expr`; GMP `mpq_t` coefficient arrays; the numeric path is
an `Expr` over `2F1`/`Gamma`. ND: binary kernel
`NDKB_LegendreQ = { ndk_LegendreQ_c, ... }`, registered `REG_B`, so
`packed_aware`. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Recurrence `O(n^2)` big rationals, cap `n = 2000`; the
numeric series converges only inside `|x| < 1` (the result is accepted only if
the series actually collapsed to a number, else symbolic). `Compile[]` lowers at
both scalar and rank-1 array shapes (`Compiled -> True`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/), [Log](../../elementary-functions/Log/), [LegendreP](../../special-functions/LegendreP/)

- DLMF §14 — Legendre functions of the second kind (§14.3, §14.7).
- Source: [`src/special_functions/legendre.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/legendre.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_dsolve_m14_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m14_stress.c)
- Tests: [`tests/test_legendreq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_legendreq.c)

## Notes & additional examples

### Notes

`LegendreQ[n, x]` is the Legendre function of the second kind `Q_n(x)`. Unlike
`P_n`, it carries a logarithm even at integer order: for integer `n` the result
is `P_n(x) L(x) + v_n(x)` with `L(x) = (1/2)(Log[1+x] - Log[1-x]) = ArcTanh[x]`
and a polynomial correction `v_n` (orders up to `n = 2000`).

A non-integer order is handled through two Frobenius `2F1` series on the cut
`|x| < 1`; the exact special value `Q_v(0)` is returned for an exact zero
argument. The associated functions `LegendreQ[n, m, x]` and
`LegendreQ[n, m, a, x]` (type `a` in 1, 2, 3) are obtained by differentiating
`Q_n`.
