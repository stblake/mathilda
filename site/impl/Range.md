---
source: src/list.c
---
`builtin_range` (in `src/list.c`) generates the arithmetic sequence for `Range[imax]` (origin 1, step 1), `Range[imin, imax]`, and `Range[imin, imax, di]`. Bounds may be integers, reals, or rationals (parsed via `is_rational`); a `double` view of each is used only for the loop-termination test (`val <= max_val + 1e-14`, or the reversed test for negative step, with a 1,000,000-element cap). The element values themselves are built exactly: a running `curr_e` starts at `imin` and is advanced each step by `evaluate(Plus[curr_e, di])`, so integer and rational ranges stay exact while any real bound promotes the elements to `EXPR_REAL`. A zero step, or an empty oriented range, yields `{}`; the result is wrapped as `List[...]`.
