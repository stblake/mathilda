# AccuracyGoal

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`AccuracyGoal`**

is an option for numerical operations (NLimit, NSum, NProduct, NIntegrate, NDSolve, NResidue, ND, NSeries, FindRoot, NRoots, NSolve) specifying how many digits of ABSOLUTE accuracy to seek.

<details>
<summary>Notes</summary>

With AccuracyGoal -\> a and PrecisionGoal -\> p, Mathilda seeks a result whose numerical error in a value of size x is below 10^-a + |x| 10^-p: a is the absolute term, p the relative term, of the combined tolerance. AccuracyGoal effectively bounds the absolute error. AccuracyGoal -\> MachinePrecision (the default) seeks ~$MachinePrecision (about 15.95) digits. AccuracyGoal -\> Automatic seeks near-full working precision (two digits below WorkingPrecision). AccuracyGoal -\> Infinity disables the absolute criterion, leaving PrecisionGoal to govern termination. Each operation refines adaptively -- growing its term count, sampling, or recursion up to a resource cap -- until the goal is met. If it cannot be met at the cap, a Head::accgl message is issued and the best approximation obtained is returned. Set WorkingPrecision at least as large as AccuracyGoal; otherwise the result may fall well short of it.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A bare option symbol, not a function

```mathematica
In[1]:= Head[AccuracyGoal]
Out[1]= Symbol
```

It only sits on the left of a rule

```mathematica
In[2]:= FullForm[AccuracyGoal -> 8]
Out[2]= Rule[AccuracyGoal, 8]
```

Consumed by the numerical op

```mathematica
In[3]:= NIntegrate[Exp[-x^2], {x, 0, 1}, AccuracyGoal -> 8]
Out[3]= 0.746824
```

## Implementation notes

**Definition.** `AccuracyGoal` is an **option name** shared by the numerical operations
(`NIntegrate`, `NDSolve`, `NSum`, `NProduct`, `NLimit`, `ND`, `NSeries`, `NResidue`,
`FindRoot`, `NRoots`, `NSolve`) specifying how many digits of **absolute** accuracy to
seek. It is a bare option symbol, not a function: no builtin, no DownValues, not even the
`Protected` attribute — its only registration is the docstring in `info_init`
(`src/info.c`). Each numerical subsystem (`src/numerical_calculus/`, `src/numerical_roots/`)
reads it out of its option list and turns it into a tolerance.

**Representation.** `AccuracyGoal` stays an inert `EXPR_SYMBOL`, appearing only on the left
of a rule, so `FullForm[AccuracyGoal -> 8]` is `Rule[AccuracyGoal, 8]`. The value is a
digit count `a`, or `MachinePrecision` (the default, ~15.95 digits), `Automatic` (two
digits below `WorkingPrecision`), or `Infinity` (disable the absolute criterion).

**Usage & limits.** With `AccuracyGoal -> a` and `PrecisionGoal -> p`, an operation seeks a
result whose error in a value of size `x` is below `10^-a + |x| 10^-p`: `a` is the absolute
term, `p` the relative term of one combined tolerance, so `AccuracyGoal` effectively bounds
the absolute error. The operation refines adaptively — more terms, finer sampling, deeper
recursion — up to a resource cap; if the goal is not met it issues a `Head::accgl` message
(which routes through the `mth_message` funnel, so `Quiet`/`Check` see it) and returns its
best approximation. Set `WorkingPrecision` at least as large as `AccuracyGoal`, or the
result may fall well short.

**Attributes:** none registered.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`AccuracyGoal` is the option shared by the numerical operations (`NIntegrate`, `NDSolve`,
`NSum`, `FindRoot`, `NSolve`, …) that says how many digits of **absolute** accuracy to
seek. It is an inert option-name symbol — no builtin, no DownValues — so it does nothing on
its own; each numerical subsystem reads it and turns it into a tolerance.

With `AccuracyGoal -> a` and `PrecisionGoal -> p` the operation seeks an error below
`10^-a + |x| 10^-p` in a value of size `x`, so `AccuracyGoal` bounds the absolute term.
The value may be a digit count, `MachinePrecision` (the default), `Automatic`, or
`Infinity` (disable the absolute criterion). If the goal cannot be met at the resource cap
the operation emits a `Head::accgl` message and returns its best approximation; keep
`WorkingPrecision` at least as large as the goal.
