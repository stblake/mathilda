# Switch

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Switch[expr, form_1, value_1, form_2, value_2, ...]`**

evaluates expr, then compares it with each of the form\_i in turn, evaluating and returning the value\_i corresponding to the first match found.

<details>
<summary>Notes</summary>

Only the value\_i corresponding to the first form\_i that matches expr is evaluated. Each form\_i is itself evaluated only when the match is tried. If the last form\_i is the pattern \_, then the corresponding value\_i is always returned if this case is reached -- it acts as a default branch. If none of the form\_i match expr, the Switch is returned unevaluated. Switch has attribute HoldRest, so the form/value pairs are not evaluated until Switch examines them. Break, Return, and Throw inside the chosen value behave as they do in any other held context.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Switch[42, _Integer, "int", _Real, "real", _, "other"]
Out[1]= "int"

In[2]:= t[e_] := Switch[e, _Plus, Together, _Times, Apart, _, Identity]; t[(1+x)/(1-x) + x/(1+x)]
Out[2]= Together

In[3]:= Switch[#, 1, one, 2, two, _, other] & /@ {1, 2, 3}
Out[3]= {one, two, other}
```

### Applications (5)

```mathematica
In[4]:= Switch[5, _Integer, "int", _Real, "real", _, "other"]
Out[4]= "int"

In[5]:= Switch[2.5, _Integer, "int", _Real, "real", _, "other"]
Out[5]= "real"

In[6]:= Switch[7, _?PrimeQ, "prime", _, "composite"]
Out[6]= "prime"

In[7]:= Table[Switch[Mod[n, 15], 0, "FizzBuzz", 3 | 6 | 9 | 12, "Fizz", 5 | 10, "Buzz", _, n], {n, 1, 15}]
Out[7]= {1, 2, "Fizz", 4, "Buzz", "Fizz", 7, 8, "Fizz", "Buzz", 11, "Fizz", 13, 14, "FizzBuzz"}

In[8]:= Switch[x, 1, "a", 2, "b"]
Out[8]= Switch[x, 1, "a", 2, "b"]
```

## Implementation notes

**Algorithm.** `builtin_switch` is `ATTR_HOLDREST`: the discriminant (arg 0) is evaluated by the standard evaluator, while the `form_i, value_i` pairs (args 1..) arrive held. It requires `argc >= 3` and odd (at least one pair, an even number of held args), otherwise returns `NULL`. It scans the pairs in order, calling `evaluate` on each `form_i` just before trying it — later forms and all values are never touched. Matching is structural via a fresh `MatchEnv` and `match(expr, form_eval, env)`; on the first match it returns a copy of the still-held `value_i`, which the outer evaluator reduces. Pattern bindings captured in the form (e.g. `x_`) are deliberately discarded with the local `MatchEnv` and never substituted into the value (the form is a pure discriminator). A trailing `_, default` pair works as a default solely because `Blank[]` matches anything. If no form matches, the call is left unevaluated (`NULL`).

- Attribute `HoldRest`; the form/value pairs are held until `Switch` examines them.
- Each `form_i` is evaluated immediately before its match is tried; only the chosen `value_i` is evaluated.
- A trailing form of `_` (Blank) acts as a catch-all default clause.
- If no `form_i` matches `expr`, the call is returned unevaluated.
- Wrong arity (no form/value pair, or an odd number of arguments after `expr`) is a usage error; the expression is returned unevaluated.
- Pattern variables bound by `form_i` (e.g. `{x_, y_}`) are *not* substituted into `value_i`; the form acts purely as a discriminator.
- `Break`, `Return`, and `Throw` inside the chosen value propagate as they do in any other held context.

**Attributes:** `HoldRest`, `Protected`.

## References

**See also:** [HoldRest](../../other-advanced/HoldRest/), [Break](../../control-flow/Break/), [Return](../../control-flow/Return/), [Throw](../../control-flow/Throw/)

- Source: [`src/cond.c`](https://github.com/stblake/mathilda/blob/main/src/cond.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_cond.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cond.c)

## Notes & additional examples

### Notes

`Switch[expr, form_1, value_1, ...]` evaluates `expr`, compares it against each
`form_i` in turn, and returns the corresponding `value_i` for the first match.
Because of the `HoldRest` attribute, the form/value pairs are held: each form is
evaluated only when its match is tried, and only the chosen value is evaluated.
A trailing `_` form acts as a catch-all default; with no matching form and no
default, the `Switch` stays unevaluated. `Break`, `Return`, and `Throw` inside
the selected value behave as in any other held context.
