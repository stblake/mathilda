# Switch

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Switch[expr, form_1, value_1, form_2, value_2, ...]
    evaluates expr, then compares it with each of the form_i in turn,
    evaluating and returning the value_i corresponding to the first
    match found.
Only the value_i corresponding to the first form_i that matches
expr is evaluated. Each form_i is itself evaluated only when the
match is tried.
If the last form_i is the pattern _, then the corresponding
value_i is always returned if this case is reached -- it acts as
a default branch.
If none of the form_i match expr, the Switch is returned
unevaluated.
Switch has attribute HoldRest, so the form/value pairs are not
evaluated until Switch examines them.
Break, Return, and Throw inside the chosen value behave as they
do in any other held context.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Switch[42, _Integer, "int", _Real, "real", _, "other"]
Out[1]= "int"

In[2]:= t[e_] := Switch[e, _Plus, Together, _Times, Apart, _, Identity]; t[(1+x)/(1-x) + x/(1+x)]
Out[2]= Together

In[3]:= Switch[#, 1, one, 2, two, _, other] & /@ {1, 2, 3}
Out[3]= {one, two, other}
```

## Implementation notes

- Attribute `HoldRest`; the form/value pairs are held until `Switch` examines them.
- Each `form_i` is evaluated immediately before its match is tried; only the chosen `value_i` is evaluated.
- A trailing form of `_` (Blank) acts as a catch-all default clause.
- If no `form_i` matches `expr`, the call is returned unevaluated.
- Wrong arity (no form/value pair, or an odd number of arguments after `expr`) is a usage error; the expression is returned unevaluated.
- Pattern variables bound by `form_i` (e.g. `{x_, y_}`) are *not* substituted into `value_i`; the form acts purely as a discriminator.
- `Break`, `Return`, and `Throw` inside the chosen value propagate as they do in any other held context.

**Attributes:** `HoldRest`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
