---
source: src/info.c
---
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
