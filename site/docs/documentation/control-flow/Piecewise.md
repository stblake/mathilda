# Piecewise

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Piecewise[{{val_1, cond_1}, {val_2, cond_2}, ...}]
    represents a piecewise function with values val_i in the regions
    defined by the conditions cond_i.
Piecewise[{{val_1, cond_1}, ...}, val]
    uses the default value val if none of the cond_i apply. The
    default for val is 0.

The cond_i are evaluated in turn until one yields True. If all
preceding cond_i yield False, the corresponding val_i of the
first True cond_i is returned. If any preceding cond_i does not
literally yield False, the Piecewise expression is returned in
symbolic form.

Only those val_i explicitly included in the returned form are
evaluated (Piecewise has attribute HoldAll). Pairs of the form
{val_i, False} are dropped, and all clauses after the first
{val_i, True} are dropped together with the default value.
Consecutive clauses with structurally identical values are
merged: their conditions are combined with Or.

Piecewise[conds] automatically evaluates to Piecewise[conds, 0].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Piecewise[{{Sin[x]/x, x < 0}, {1, x == 0}}, -x^2/100 + 1]
Out[1]= Piecewise[{{Sin[x]/x, x < 0}, {1, x == 0}}, (-x^2)/100 + 1]

In[2]:= Piecewise[{{e1, True}, {e2, d2}, {e3, d3}}]
Out[2]= e1

In[3]:= Piecewise[{{a, d1}, {b, d2}, {c, False}, {d, d4}}, ef]
Out[3]= Piecewise[{{a, d1}, {b, d2}, {d, d4}}, ef]

In[4]:= Piecewise[{{a, d1}, {b, d2}, {b, d3}, {c, d4}}, ef]
Out[4]= Piecewise[{{a, d1}, {b, d2 || d3}, {c, d4}}, ef]

In[5]:= Piecewise[{{Sin[x]/x, x < 0}, {1, x == 0}}, -x^2/100 + 1] /. x -> 5
Out[5]= 3/4
```

## Implementation notes

**Algorithm.** `builtin_piecewise` (in `src/cond.c`, not a separate `piecewise.c`) is `ATTR_HOLDALL`. It takes `Piecewise[{{val,cond},...}]` or `Piecewise[{...}, default]`; the first argument must be a `List` of two-element `{value, condition}` `List`s (validated by `piecewise_is_pair`), else `NULL`. It walks the clauses, calling `evaluate` only on each condition (values stay held): `{v, False}` clauses are dropped, `{v, True}` is kept and stops the scan (all later clauses and the default become unreachable), and any other condition is kept in its evaluated form.

Surviving clauses are then compacted: a run of consecutive clauses with structurally equal values (`expr_eq`) is merged into one clause whose condition is `Or[c1, c2, ...]`, reduced through `eval_and_free` so `True`/`False` simplifications fire.

**Result selection.** Zero survivors yield the default (copied, or `0` if omitted). A single survivor whose condition simplified to `True` returns its (held) value directly. Otherwise it rebuilds a symbolic `Piecewise[{{v,c},...}, default]`, dropping the default iff the final clause is unconditionally `True` (unreachable). If the rebuilt expression is `expr_eq` to the input, it returns `NULL` to signal "no change" for efficient fixed-pointing.

**Data structures.** Two parallel growable `Expr**` arrays (`out_vals`, `out_conds`) accumulate survivors before the merge/select phase.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/cond.c`](https://github.com/stblake/mathilda/blob/main/src/cond.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
