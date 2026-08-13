# BernoulliB

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BernoulliB[n]`**

gives the Bernoulli number B\_n.

**`BernoulliB[n, x]`**

gives the Bernoulli polynomial B\_n(x).

<details>
<summary>Notes</summary>

Non-negative integer n gives the exact rational B\_n (odd n \> 1 give 0, B\_0 = 1, B\_1 = -1/2); an inexact integer-valued n evaluates it at machine or arbitrary (MPFR) precision. BernoulliB\[n, x\] expands the degree-n polynomial with exact rational coefficients, staying symbolic in x or evaluating numerically when x is inexact. Listable.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Table[BernoulliB[k], {k, 0, 10}]
Out[1]= {1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0, 5/66}

In[2]:= BernoulliB[3, z]
Out[2]= 1/2 z - 3/2 z^2 + z^3
```

### Applications (6)

```mathematica
In[3]:= BernoulliB[10]
Out[3]= 5/66

In[4]:= BernoulliB[20]
Out[4]= -174611/330

In[5]:= Table[BernoulliB[2 k], {k, 0, 6}]
Out[5]= {1, 1/6, -1/30, 1/42, -1/30, 5/66, -691/2730}

In[6]:= BernoulliB[50]
Out[6]= 495057205241079648212477525/66

In[7]:= BernoulliB[4, x]
Out[7]= -1/30 + x^2 - 2 x^3 + x^4

In[8]:= Sum[k^2, {k, 1, n}]
Out[8]= 1/6 n (1 + n) (1 + 2 n)
```

## Algorithm

Mathilda -- Bernoulli numbers and polynomials.

```text
  BernoulliB[n]      Bernoulli number      B_n
  BernoulliB[n, x]   Bernoulli polynomial  B_n(x)
```

The Bernoulli polynomials are the coefficients of the generating function

```text
  t e^(x t) / (e^t - 1) = Sum_{n>=0} B_n(x) t^n / n!,
```

and the Bernoulli numbers are B_n = B_n(0). For odd n the numbers vanish except B_1 = -1/2.

Evaluation is layered so each argument shape takes the cheapest exact or numeric route:

```text
  exact non-negative integer n   ->  exact rational B_n (cached recurrence)
  inexact integer-valued n       ->  the rational, numericalised (Real/MPFR)
  BernoulliB[n, x]               ->  Sum_j C(n,j) B_{n-j} x^j, built with
                                     exact rational coefficients then
                                     evaluated (exact x stays exact,
                                     inexact x evaluates numerically)
  everything else                ->  stays symbolic (return NULL)
```

Bernoulli numbers B_n are computed by the recurrence

```text
  B_0 = 1,   B_m = -1/(m+1) Sum_{k=0}^{m-1} C(m+1, k) B_k   (m >= 1),
```

using exact GMP rationals in a lazily-grown, process-lifetime cache. (The same construction lives in src/zeta.c and src/polygamma.c; it is replicated here to keep this module self-contained.)

Attributes: Listable, Protected.

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_bernoullib.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bernoullib.c)
- Tests: [`tests/test_boolean.c`](https://github.com/stblake/mathilda/blob/main/tests/test_boolean.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_possiblezeroq_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_possiblezeroq_stress.c)

## Notes & additional examples

### Notes

`BernoulliB[n]` gives the Bernoulli number `B_n`; `BernoulliB[n, x]` gives the Bernoulli polynomial. Non-negative integer `n` returns the exact rational (`B_0 = 1`, `B_1 = -1/2`, odd `n > 1` give 0); inexact integer-valued `n` evaluates numerically at machine or MPFR precision. `BernoulliB` is Listable.
