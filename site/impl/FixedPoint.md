---
source: src/funcprog.c
---
**Algorithm.** `builtin_fixedpoint` calls `fixedpoint_impl(res, as_list=false)`. It seeds an `ExprBuf` history with a copy of the start value and drives the generic iterator `iter_run` with the `fixedpoint_step` callback. Each step applies `f` to the most recent value (`apply_unary`) and compares the new value to the previous one: by default with structural equality (`expr_eq`), or with a user `SameTest -> g` option (applying `g` and testing for `True`). When they agree the iteration halts (`ITER_STEP_HALT_ADD`); otherwise it continues. A `MaxIterations` count (integer or `Infinity`) is parsed by `parse_fp_opts`; unbounded runs are capped at `ITER_SAFETY_CAP` (exceeding it returns NULL). Throw/Abort/Quit/Return markers produced by `f` are propagated out early. `ebuf_finalize(..., as_list=false)` returns just the final value; the sibling `FixedPointList` returns the whole history.

**Data structures.** `ExprBuf` (a growable `Expr*` vector) is the iteration history; `FixedPointCtx` carries `f`, the optional same-test, and the throw-propagation flag.
