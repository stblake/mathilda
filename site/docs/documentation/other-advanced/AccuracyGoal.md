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

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
