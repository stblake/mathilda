---
source: src/expand_power.c
---
**Algorithm.** `builtin_powerexpand` distributes powers over products and collapses nested
powers and logarithms, applying the rewrites `(a b)^c -> a^c b^c`, `(a^b)^c -> a^(b c)`,
`Log[a b] -> Log[a] + Log[b]`, `Log[a^b] -> b Log[a]`, and `Arg[a b] -> Arg[a] + Arg[b]`. Since
`Sqrt[x]` is stored as `Power[x, 1/2]` and `Log[1/z]` as `Log[Power[z, -1]]`, these cover
`Sqrt[x y]`, `Sqrt[z^2] -> z`, `Log[1/z] -> -Log[z]` with no special-casing. The transform
(`pe_rec`) is applied top-down to a fixed point, so an outermost rule fires first and the result
is reprocessed. Three modes are selected by the `Assumptions` option: **Automatic** (default,
the textbook transforms, valid for positive-real bases / integer exponents), **`-> True`** (emit
universally-correct formulas with a branch-correction term built from `Floor`/`Arg`/`Im`/`E`/
`I`/`Pi`), and **`-> assum`** (emit the True-mode formula then refine the correction terms under
the assumptions via `pe_refine`, degrading gracefully to the symbolic form where the reasoning
runs out). It threads over `List`, equations, inequalities, and logic heads, and supports the
variable-restricted form `PowerExpand[expr, {x1, …}]`.
