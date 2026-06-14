# NSum

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NSum[f, {i, imin, imax}]
    gives a numerical approximation to the sum of f for i from imin to imax.

NSum[f, {i, imin, imax, di}] uses step di. imax may be Infinity. NSum[f, {i, ...}, {j, ...}, ...] evaluates a multidimensional sum (an inner bound may depend on an outer index). The index is localised (HoldAll). Method -> Automatic picks Euler-Maclaurin for monotone series, the Cohen-Villegas-Zagier method for alternating series, and Wynn's epsilon (partial-sum acceleration) otherwise; large finite sums use the difference of two infinite tails. Machine or arbitrary precision via WorkingPrecision.

Options: Method (Automatic | EulerMaclaurin | AlternatingSigns | WynnEpsilon), WorkingPrecision (default MachinePrecision), NSumTerms (head terms summed explicitly, default 15), NSumExtraTerms, WynnDegree, VerifyConvergence (default True; a divergent sum gives ComplexInfinity), AccuracyGoal, PrecisionGoal.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NSum[(-5)^i/i!, {i, 0, Infinity}, NSumTerms -> 25] - Exp[-5]
Out[1]= 1.44069e-15

In[2]:= NSum[1/i^2, {i, 1, Infinity}] - Pi^2/6 // N
Out[2]= 2.22045e-16

In[3]:= NSum[1/n^(11/10), {n, 1, Infinity}, WorkingPrecision -> 40] - Zeta[11/10]
Out[3]= 1.4693679385278593849609206715278070972733e-39

In[4]:= NSum[(-1)^x/(1 + (x - 12)^2), {x, 0, Infinity}, Method -> "AlternatingSigns", WorkingPrecision -> 30]
Out[4]= 0.2751938594139530395689715615907

In[5]:= NSum[1/2^i, {i, 0, Infinity, 2}]
Out[5]= 1.33333

In[6]:= NSum[Log[x]/x^(2 + 2 I), {x, 1, Infinity}]
Out[6]= -0.182175 - 0.136618*I

In[7]:= NSum[1/i^2, {i, 100, 10^6}]
Out[7]= 0.0100492

In[8]:= NSum[(-1)^n (2/n)^k/k^2, {n, 2, Infinity}, {k, 1, n}]
Out[8]= 0.770188
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
