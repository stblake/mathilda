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

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= HurwitzZeta[s, 1/2]
Out[1]= (-1 + 2^s) Zeta[s]

In[2]:= HurwitzZeta[3, -3.5]
Out[2]= 0.0307784
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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Zeta](../../special-functions/Zeta/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flint_bridge.c)
- Tests: [`tests/test_hurwitzzeta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hurwitzzeta.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
