---
source: src/info.c
---
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
