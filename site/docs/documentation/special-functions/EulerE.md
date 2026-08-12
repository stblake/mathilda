# EulerE

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EulerE[n]`**

gives the Euler number E\_n.

**`EulerE[n, x]`**

gives the Euler polynomial E\_n(x).

<details>
<summary>Notes</summary>

Non-negative integer n gives the exact integer E\_n (odd n give 0, E\_0 = 1, E\_2 = -1, E\_4 = 5); an inexact integer-valued n evaluates it at machine or arbitrary (MPFR) precision. EulerE\[n, x\] expands the degree-n polynomial with exact rational coefficients, staying symbolic in x or evaluating numerically when x is inexact; EulerE\[n, 1/2\] folds to 2^-n EulerE\[n\]. Listable.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Table[EulerE[k], {k, 0, 10}]
Out[1]= {1, 0, -1, 0, 5, 0, -61, 0, 1385, 0, -50521}

In[2]:= EulerE[3, z]
Out[2]= 1/4 - 3/2 z^2 + z^3
```

### Applications (5)

```mathematica
In[1]:= EulerE[4]
Out[1]= 5
```

```mathematica
In[1]:= Table[EulerE[2 n], {n, 0, 6}]
Out[1]= {1, -1, 5, -61, 1385, -50521, 2702765}
```

```mathematica
In[1]:= EulerE[6, x]
Out[1]= -3 x + 5 x^3 - 3 x^5 + x^6
```

```mathematica
In[1]:= Sum[EulerE[2 k] / (2 k)! Pi^(2 k), {k, 0, 5}]
Out[1]= 1 - 1/2 Pi^2 + 5/24 Pi^4 - 61/720 Pi^6 + 277/8064 Pi^8 - 50521/3628800 Pi^10
```

```mathematica
In[1]:= N[EulerE[5, 1/3], 40]
Out[1]= -0.24897119341563786008230452674897119341565
```

## Algorithm

Mathilda -- Euler numbers and polynomials.

```text
  EulerE[n]      Euler number      E_n
  EulerE[n, x]   Euler polynomial  E_n(x)
```

The Euler polynomials are the coefficients of the generating function

```text
  2 e^(x t) / (e^t + 1) = Sum_{n>=0} E_n(x) t^n / n!,
```

and the Euler numbers are E_n = 2^n E_n(1/2). For odd n the numbers vanish; E_0 = 1, E_2 = -1, E_4 = 5, E_6 = -61, ...

Evaluation is layered so each argument shape takes the cheapest exact or numeric route:

```text
  exact non-negative integer n   ->  exact integer E_n (cached recurrence)
  inexact integer-valued n       ->  the integer, numericalised (Real/MPFR)
  EulerE[n, x]                   ->  the degree-n polynomial in monomial
                                     form with exact rational coefficients,
                                     then evaluated (exact x stays exact,
                                     inexact x evaluates numerically)
  EulerE[n, 1/2], symbolic n     ->  2^-n EulerE[n]
  everything else                ->  stays symbolic (return NULL)
```

Euler numbers E_n are computed by the recurrence

```text
  E_0 = 1,   E_{2m} = -Sum_{k=0}^{m-1} C(2m, 2k) E_{2k}   (m >= 1),
```

with odd-index numbers identically zero, using exact GMP integers in a lazily-grown, process-lifetime cache.

The Euler polynomial coefficients are obtained from the Taylor expansion

```text
about x = 1/2.  Writing E_n(x) = Sum_i c_i x^i, the algebra collapses (the
```

powers of two in C(n,j) E_{n-j}/2^{n-j} (x-1/2)^j combine to a single 2^{n-i}) to the all-integer inner sum

```text
  S_i = Sum_{j=i}^{n} (-1)^{j-i} C(n,j) C(j,i) E_{n-j},   c_i = S_i / 2^{n-i}.
```

Attributes: Listable, Protected.

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## See also

[N](../../arithmetic/N/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_eulere.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eulere.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

`EulerE[n]` is the integer Euler number `E_n` (odd `n` vanish, `E_0 = 1`),
and `EulerE[n, x]` is the degree-`n` Euler polynomial with exact rational
coefficients. The truncated `secant`-style series above is the partial sum
of `sec(Pi/2)`'s generating expansion; the polynomial form stays symbolic in
`x` or, given an inexact argument, evaluates to arbitrary precision.
