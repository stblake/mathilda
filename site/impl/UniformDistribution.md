---
source: src/ml/dist.c
---
**Definition.** `UniformDistribution[{lo, hi}]` is a **distribution object** representing the
continuous uniform distribution on `[lo, hi]`; `UniformDistribution[]` is uniform on
`{0, 1}`. It is an inert, `Protected` head, not a function: `UniformDistribution[{0, 1}]`
evaluates to itself, and its meaning is supplied by the consumers in `src/ml/dist.c` (`PDF`,
`RandomVariate`), which parse it into an internal `MlDist` record.

**Representation.** The object stays an `EXPR_FUNCTION` with head `UniformDistribution`, so
`Head[UniformDistribution[{0, 1}]]` is `UniformDistribution`. When a consumer reads it, the
parser (`ml_parse_distribution` in `dist.c`) sets the internal kind to `ML_D_UNIFORM` with
endpoints `a = lo`, `b = hi`; the no-argument form defaults to `a = 0, b = 1`. A spec whose
upper endpoint does not exceed the lower (`b > a` fails) is rejected, so the object never
yields a degenerate density.

**Usage & limits.** `PDF[UniformDistribution[{lo, hi}], x]` returns `1/(hi - lo)` for
`lo <= x <= hi` (endpoints included, Wolfram's convention for a continuous uniform) and `0`
strictly outside; `RandomVariate[UniformDistribution[{lo, hi}]]` draws a sample. These
numeric consumers evaluate at machine precision, so `PDF[...]` comes back as a machine real.
Symbolic moment functions such as `Mean`, `Variance` and `CDF` are **not** yet wired to this
object and are left unevaluated. It sits beside `NormalDistribution` in the same module.
