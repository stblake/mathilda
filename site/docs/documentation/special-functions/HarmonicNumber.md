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

## Examples

_No verified examples yet for this function._

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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[BernoulliB](../../special-functions/BernoulliB/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_harmonicnumber.c`](https://github.com/stblake/mathilda/blob/main/tests/test_harmonicnumber.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_possiblezeroq_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_possiblezeroq_stress.c)
