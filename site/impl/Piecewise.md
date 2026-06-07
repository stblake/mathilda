---
source: src/cond.c
---
**Algorithm.** `builtin_piecewise` (in `src/cond.c`, not a separate `piecewise.c`) is `ATTR_HOLDALL`. It takes `Piecewise[{{val,cond},...}]` or `Piecewise[{...}, default]`; the first argument must be a `List` of two-element `{value, condition}` `List`s (validated by `piecewise_is_pair`), else `NULL`. It walks the clauses, calling `evaluate` only on each condition (values stay held): `{v, False}` clauses are dropped, `{v, True}` is kept and stops the scan (all later clauses and the default become unreachable), and any other condition is kept in its evaluated form.

Surviving clauses are then compacted: a run of consecutive clauses with structurally equal values (`expr_eq`) is merged into one clause whose condition is `Or[c1, c2, ...]`, reduced through `eval_and_free` so `True`/`False` simplifications fire.

**Result selection.** Zero survivors yield the default (copied, or `0` if omitted). A single survivor whose condition simplified to `True` returns its (held) value directly. Otherwise it rebuilds a symbolic `Piecewise[{{v,c},...}, default]`, dropping the default iff the final clause is unconditionally `True` (unreachable). If the rebuilt expression is `expr_eq` to the input, it returns `NULL` to signal "no change" for efficient fixed-pointing.

**Data structures.** Two parallel growable `Expr**` arrays (`out_vals`, `out_conds`) accumulate survivors before the merge/select phase.
