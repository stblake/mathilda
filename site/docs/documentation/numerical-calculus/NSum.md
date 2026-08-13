# NSum

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NSum[f, {i, imin, imax}]`**

gives a numerical approximation to the sum of f for i from imin to imax.

**`NSum[f, {i, imin, imax, di}] uses step di. imax may be Infinity. NSum[f, {i, ...}, {j, ...}, ...] evaluates a multidimensional sum (an inner bound may depend on an outer index). The index is localised (HoldAll). Method -> Automatic picks Euler-Maclaurin for monotone series, the Cohen-Villegas-Zagier method for alternating series, and Wynn's epsilon (partial-sum acceleration) otherwise, with Levin's u-transform as a last resort; large finite sums use the difference of two infinite tails. Method -> "Levin" forces Levin's transformation ("LevinU" | "LevinT" | "LevinV" select the u/t/v variant). Machine or arbitrary precision via WorkingPrecision.`**

<details>
<summary>Notes</summary>

Options: Method (Automatic | EulerMaclaurin | AlternatingSigns | WynnEpsilon | "Levin"), WorkingPrecision (default MachinePrecision), NSumTerms (head terms summed explicitly, default 15), NSumExtraTerms, WynnDegree, VerifyConvergence (default True; a divergent sum gives ComplexInfinity), AccuracyGoal, PrecisionGoal.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= NSum[1/i^2, {i, 1, Infinity}] - Pi^2/6 // N
Out[1]= 2.22045e-16

In[2]:= NSum[1/2^i, {i, 0, Infinity, 2}]
Out[2]= 1.33333

In[3]:= NSum[Log[x]/x^(2 + 2 I), {x, 1, Infinity}]
Out[3]= -0.182175 - 0.136618*I

In[4]:= NSum[1/i^2, {i, 100, 10^6}]
Out[4]= 0.0100492

In[5]:= NSum[(-1)^n (2/n)^k/k^2, {n, 2, Infinity}, {k, 1, n}]
Out[5]= 0.770188
```

### Options (3)

```mathematica
In[6]:= NSum[(-5)^i/i!, {i, 0, Infinity}, NSumTerms -> 25] - Exp[-5]
Out[6]= -2.4182e-15

In[7]:= NSum[1/n^(11/10), {n, 1, Infinity}, WorkingPrecision -> 40] - Zeta[11/10]
Out[7]= 1.4693679385278593849609206715278070972733e-39

In[8]:= NSum[(-1)^x/(1 + (x - 12)^2), {x, 0, Infinity}, Method -> "AlternatingSigns", WorkingPrecision -> 30]
Out[8]= 0.2751938594139530395689715615907
```

### Applications (5)

```mathematica
In[9]:= NSum[1/n^2, {n, 1, Infinity}]
Out[9]= 1.64493

In[10]:= NSum[(-1)^(n+1)/n, {n, 1, Infinity}]
Out[10]= 0.693147

In[11]:= NSum[1/n^2, {n, 1, Infinity}, WorkingPrecision -> 30]
Out[11]= 1.644934066848226436472415166649

In[12]:= NSum[1/n^4, {n, 1, Infinity}, WorkingPrecision -> 30]
Out[12]= 1.082323233711138191516003696543

In[13]:= NSum[1/n^2, {n, 1, Infinity}, Method -> "Levin"]
Out[13]= 1.64493
```

## Algorithm

```text
nsum.c — NSum[f, {i, imin, imax (, di)}, opts]   (see nsum.h)
```

Strategy -------- NSum holds its arguments, evaluates the iterator bounds, then Block-localises

```text
the index and evaluates the summand once per term.  Terms are reindexed to
```

k = 0, 1, 2, … with the actual index value x_k = imin + k·di, so a step di is handled uniformly and multidimensional sums fall out by making the summand of the outer sum an inner NSum[...] (HoldAll + localisation lets a dependent inner bound such as {k,1,n} see the bound outer index).

Methods are layered: this file currently provides Direct (small finite sums) and WynnEpsilon (partial-sum extrapolation, shared seqaccel kernels), machine

```text
and MPFR, real and complex.  Euler–Maclaurin and Cohen–Villegas–Zagier are
```

added on top of the same term machinery.

Memory: receives `res` owned by the evaluator; returns a fresh Expr* on

```text
success or NULL (unevaluated).  Never frees `res`.  Every temporary index
```

binding is removed on all return paths.

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Block](../../scoping-constructs/Block/), [Chop](../../elementary-functions/Chop/), [Integrate](../../calculus/Integrate/), [D](../../calculus/D/), [BernoulliB](../../special-functions/BernoulliB/), [NLimit](../../numerical-calculus/NLimit/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_nprod.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nprod.c)
- Tests: [`tests/test_nsum.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nsum.c)

## Notes & additional examples

### Notes

`NSum[f, {i, imin, imax}]` numerically sums a series, with `imax` allowed to be
`Infinity`. The first two cases recover the Basel sum `Pi^2/6 = 1.64493...` and
the alternating harmonic sum `Log[2] = 0.693147...`. With `WorkingPrecision -> 30`
the Basel sum is computed to 30 digits, and `Sum[1/n^4]` returns
`Pi^4/90 = 1.082323233711...`. `Method -> Automatic` chooses Euler–Maclaurin for
monotone series, the Cohen–Villegas–Zagier method for alternating series, and
Wynn's epsilon otherwise, with Levin's u-transform as a last resort. Any
accelerator can be forced: `Method -> "Levin"` (`"LevinU"`/`"LevinT"`/`"LevinV"`)
selects Levin's transformation, which reaches full `WorkingPrecision` on smooth
series. With `VerifyConvergence -> True` (default) a divergent sum gives
`ComplexInfinity`.
