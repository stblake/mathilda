# If

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`If[cond, t]`**

gives t if cond evaluates to True; gives Null otherwise.

**`If[cond, t, f]`**

gives t if cond is True, f if False, and is left unevaluated if cond is neither.

**`If[cond, t, f, u]`**

also supplies u as the result when cond is neither True nor False.

<details>
<summary>Notes</summary>

If has attribute HoldRest: only the branch chosen by cond is evaluated.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= If[True, x, y]
Out[1]= x

In[2]:= If[a < b, 1, 0, Indeterminate]
Out[2]= Indeterminate
```

### Applications (4)

```mathematica
In[3]:= If[3 > 2, yes, no]
Out[3]= yes

In[4]:= If[x > 0, pos, neg]
Out[4]= If[x > 0, pos, neg]

In[5]:= If[1 == 2, a, b, c]
Out[5]= b

In[6]:= If[PrimeQ[7], prime, composite]
Out[6]= prime
```

## Implementation notes

**Algorithm.** `builtin_if` accepts 2–4 arguments and is registered with `ATTR_HOLDREST`, so only the condition (arg 0) is evaluated by the standard evaluator before the builtin runs; the branches stay held. The handler inspects the evaluated condition: if it is the interned symbol `True` it returns a copy of arg 1; if it is `False` it returns a copy of arg 2 (the else branch), or `Null` when no else branch was supplied. The returned branch is still "held" data that the outer fixed-point evaluator then reduces. If the condition is neither `True` nor `False`, a 4-argument call returns a copy of the fourth argument (the "neither" / indeterminate branch), and otherwise the call returns `NULL`, leaving `If[...]` unevaluated so a symbolic condition flows through unchanged. Truth testing is pointer equality against `SYM_True`/`SYM_False`, never a string compare.

- `HoldRest`, evaluating only the chosen branch.
- Remains unevaluated if the condition is undetermined and `u` is not provided.
- `If[condition, t]` returns `Null` if `condition` evaluates to `False`.

**Attributes:** `HoldRest`, `Protected`.

## References

**See also:** [HoldRest](../../other-advanced/HoldRest/)

- Source: [`src/cond.c`](https://github.com/stblake/mathilda/blob/main/src/cond.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_boolean.c`](https://github.com/stblake/mathilda/blob/main/tests/test_boolean.c)
- Tests: [`tests/test_catch_throw.c`](https://github.com/stblake/mathilda/blob/main/tests/test_catch_throw.c)
- Tests: [`tests/test_cherry_ei.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_ei.c)

## Notes & additional examples

### Notes

`If[cond, t, f]` evaluates `t` when `cond` is `True` and `f` when it is `False`.
`If` has the `HoldRest` attribute, so only the selected branch is evaluated — the
other is never run. When the condition is neither `True` nor `False` (e.g. the
symbolic `x > 0`), `If` returns unevaluated rather than guessing. The optional
fourth argument `u` is returned for that indeterminate case in the four-argument
form `If[cond, t, f, u]`.
