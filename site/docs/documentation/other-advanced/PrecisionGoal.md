# PrecisionGoal

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
PrecisionGoal
    is an option for numerical operations (NIntegrate, NDSolve, NLimit,
    NSum, and the others accepting AccuracyGoal) specifying how many
    digits of RELATIVE precision to seek.

With PrecisionGoal -> p and AccuracyGoal -> a, Mathilda seeks a result
whose numerical error in a value of size x is below 10^-a + |x| 10^-p;
p is the relative term of that combined tolerance and effectively
bounds the relative error.

PrecisionGoal -> Automatic (the default) seeks near-full working
precision (two digits below WorkingPrecision). PrecisionGoal -> Infinity
disables the relative criterion, leaving AccuracyGoal to govern
termination.

When the goal cannot be met the operation issues a Head::accgl message
and returns its best approximation. Set WorkingPrecision at least as
large as PrecisionGoal.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
