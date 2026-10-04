# PrecisionGoal

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`PrecisionGoal`**

is an option for numerical operations (NIntegrate, NDSolve, NLimit, NSum, and the others accepting AccuracyGoal) specifying how many digits of RELATIVE precision to seek.

<details>
<summary>Notes</summary>

With PrecisionGoal -\> p and AccuracyGoal -\> a, Mathilda seeks a result whose numerical error in a value of size x is below 10^-a + |x| 10^-p; p is the relative term of that combined tolerance and effectively bounds the relative error. PrecisionGoal -\> Automatic (the default) seeks near-full working precision (two digits below WorkingPrecision). PrecisionGoal -\> Infinity disables the relative criterion, leaving AccuracyGoal to govern termination. When the goal cannot be met the operation issues a Head::accgl message and returns its best approximation. Set WorkingPrecision at least as large as PrecisionGoal.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A bare option symbol, not a function

```mathematica
In[1]:= Head[PrecisionGoal]
Out[1]= Symbol
```

It only sits on the left of a rule

```mathematica
In[2]:= FullForm[PrecisionGoal -> 10]
Out[2]= Rule[PrecisionGoal, 10]
```

Consumed by the numerical op

```mathematica
In[3]:= NIntegrate[1/x, {x, 1, 2}, PrecisionGoal -> 10]
Out[3]= 0.693147
```

## Implementation notes

**Definition.** `PrecisionGoal` is an **option name** shared by the numerical operations
(`NIntegrate`, `NDSolve`, `NLimit`, `NSum`, and the others that accept `AccuracyGoal`)
specifying how many digits of **relative** precision to seek. It is a bare option symbol,
not a function: no builtin, no DownValues, not even the `Protected` attribute — its only
registration is the docstring in `info_init` (`src/info.c`). Each numerical subsystem
(`src/numerical_calculus/`, `src/numerical_roots/`) reads it from its option list and turns
it into a tolerance.

**Representation.** `PrecisionGoal` stays an inert `EXPR_SYMBOL`, appearing only on the
left of a rule, so `FullForm[PrecisionGoal -> 10]` is `Rule[PrecisionGoal, 10]`. The value
is a digit count `p`, or `Automatic` (the default: two digits below `WorkingPrecision`), or
`Infinity` (disable the relative criterion, leaving `AccuracyGoal` to govern termination).

**Usage & limits.** With `PrecisionGoal -> p` and `AccuracyGoal -> a`, an operation seeks a
result whose error in a value of size `x` is below `10^-a + |x| 10^-p`; `p` is the relative
term and effectively bounds the relative error. The operation refines adaptively up to a
resource cap, and if the goal is not met it issues a `Head::accgl` message (routed through
the `mth_message` funnel, so `Quiet`/`Check` see it) and returns its best approximation.
Set `WorkingPrecision` at least as large as `PrecisionGoal`.

**Attributes:** none registered.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`PrecisionGoal` is the option shared by the numerical operations (`NIntegrate`, `NDSolve`,
`NSum`, …) that says how many digits of **relative** precision to seek. It is an inert
option-name symbol — no builtin, no DownValues — so it does nothing on its own; each
numerical subsystem reads it and turns it into a tolerance.

With `PrecisionGoal -> p` and `AccuracyGoal -> a` the operation seeks an error below
`10^-a + |x| 10^-p` in a value of size `x`, so `PrecisionGoal` bounds the relative term.
The value may be a digit count, `Automatic` (the default, two digits below
`WorkingPrecision`), or `Infinity` (disable the relative criterion). If the goal cannot be
met at the resource cap, a `Head::accgl` message is emitted and the best approximation
returned; keep `WorkingPrecision` at least as large as the goal.
