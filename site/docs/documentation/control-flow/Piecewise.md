# Piecewise

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Piecewise[{{val_1, cond_1}, {val_2, cond_2}, ...}]`**

represents a piecewise function with values val\_i in the regions defined by the conditions cond\_i.

**`Piecewise[{{val_1, cond_1}, ...}, val]`**

uses the default value val if none of the cond\_i apply. The default for val is 0.

**`Piecewise[conds] automatically evaluates to Piecewise[conds, 0].`**

<details>
<summary>Notes</summary>

The cond\_i are evaluated in turn until one yields True. If all preceding cond\_i yield False, the corresponding val\_i of the first True cond\_i is returned. If any preceding cond\_i does not literally yield False, the Piecewise expression is returned in symbolic form. Only those val\_i explicitly included in the returned form are evaluated (Piecewise has attribute HoldAll). Pairs of the form {val\_i, False} are dropped, and all clauses after the first {val\_i, True} are dropped together with the default value. Consecutive clauses with structurally identical values are merged: their conditions are combined with Or.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

### Applications (5)

```mathematica
In[6]:= Piecewise[{{-x, x < 0}, {x, x >= 0}}]
Out[6]= Piecewise[{{-x, x < 0}, {x, x >= 0}}, 0]

In[7]:= Piecewise[{{1, x > 0}, {-1, x < 0}}, 0] /. x -> 5
Out[7]= 1

In[8]:= Piecewise[{{1, x > 0}}] /. x -> -2
Out[8]= 0

In[9]:= Piecewise[{{a, x == 1}, {a, x == 2}, {b, x == 3}}]
Out[9]= Piecewise[{{a, x == 1 || x == 2}, {b, x == 3}}, 0]

In[10]:= D[Piecewise[{{x^2, x < 0}, {x^3, x >= 0}}], x]
Out[10]= Piecewise[{{2 x, x < 0}, {3 x^2, x >= 0}}, 0]
```

## Implementation notes

**Algorithm.** `builtin_piecewise` (in `src/cond.c`, not a separate `piecewise.c`) is `ATTR_HOLDALL`. It takes `Piecewise[{{val,cond},...}]` or `Piecewise[{...}, default]`; the first argument must be a `List` of two-element `{value, condition}` `List`s (validated by `piecewise_is_pair`), else `NULL`. It walks the clauses, calling `evaluate` only on each condition (values stay held): `{v, False}` clauses are dropped, `{v, True}` is kept and stops the scan (all later clauses and the default become unreachable), and any other condition is kept in its evaluated form.

Surviving clauses are then compacted: a run of consecutive clauses with structurally equal values (`expr_eq`) is merged into one clause whose condition is `Or[c1, c2, ...]`, reduced through `eval_and_free` so `True`/`False` simplifications fire.

**Result selection.** Zero survivors yield the default (copied, or `0` if omitted). A single survivor whose condition simplified to `True` returns its (held) value directly. Otherwise it rebuilds a symbolic `Piecewise[{{v,c},...}, default]`, dropping the default iff the final clause is unconditionally `True` (unreachable). If the rebuilt expression is `expr_eq` to the input, it returns `NULL` to signal "no change" for efficient fixed-pointing.

**Data structures.** Two parallel growable `Expr**` arrays (`out_vals`, `out_conds`) accumulate survivors before the merge/select phase.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/cond.c`](https://github.com/stblake/mathilda/blob/main/src/cond.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_cond.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cond.c)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_fresnelc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fresnelc.c)
- Tests: [`tests/test_fresnels.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fresnels.c)

## Notes & additional examples

### Notes

`Piecewise[{{v1, c1}, {v2, c2}, ...}, def]` is the symbolic conditional: the `ci`
are tested in order and the value of the first `True` condition is returned, or
`def` (default `0`) if all are `False`. If any earlier condition is not literally
`False` the whole expression is held symbolic. It has attribute `HoldAll`, so only
the `vi` actually returned are evaluated. On construction it canonicalises:
`{vi, False}` clauses are dropped, everything after the first `{vi, True}` is
discarded, and consecutive clauses with identical values are merged via `Or`.
Differentiation and the elementary symbolic machinery thread through the branches.
