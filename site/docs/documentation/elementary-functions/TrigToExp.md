# TrigToExp

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
TrigToExp[expr]
    rewrites circular and hyperbolic trigonometric functions (and their
    inverses) in expr in terms of exponentials and logarithms.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= TrigToExp[Cos[x]]
Out[1]= 1/2 E^(-I x) + 1/2 E^(I x)

In[2]:= TrigToExp[ArcSin[x]]
Out[2]= -I Log[I x + Sqrt[1 - x^2]]

In[3]:= TrigToExp[{Sin[x], Cos[x], Tan[x]}]
Out[3]= {(-1/2*I) E^(I x) + (1/2*I) E^(-I x), 1/2 E^(-I x) + 1/2 E^(I x), -I E^(I x)/(E^(-I x) + E^(I x)) + I E^(-I x)/(E^(-I x) + E^(I x))}
```

## Implementation notes

**Algorithm.** `builtin_trigtoexp` rewrites circular, hyperbolic, and their
inverse functions into complex-exponential / logarithmic form. It applies the
`trig_to_exp_rules` rule list (a `ReplaceAll`) — e.g. `Sin[x] :> I/2 E^(-I x) -
I/2 E^(I x)`, `Cos[x] :> (E^(-I x) + E^(I x))/2`, `Tan`/`Cot`/`Sec`/`Csc`,
the `Sinh`/`Cosh`/`Tanh`/... hyperbolic analogues, and the `ArcSin :> -I Log[...]`
family of inverse identities — then runs `Expand` to distribute the result into
a flat sum. The `Times`/`Power`-level trig canonicalizer is suppressed for the
duration (`trig_canon_suppress_inc`/`dec`) so the intermediate `Sin/Cos` forms
are not prematurely re-collapsed back to `Tan` etc. before the rule fires.

**Data structures.** The rules are a single `parse_expression`'d `List` of
`RuleDelayed`s, built once in `trigsimp_init` and stored in the static
`trig_to_exp_rules`. Results are memoized through the active `FactorMemo`
(keyed under a `TrigToExp[arg]` head) so repeated Simplify candidate-set calls
on identical subexpressions hit the cache instead of re-running `ReplaceAll +
Expand`, which costs 30–130 ms on Tan-rich inputs.

- `Listable`, `Protected`.
- Operates on both circular and hyperbolic functions, and their inverses.
- Rewrites recursively, including trig heads nested inside the *argument* of an

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/simp/trigsimp.c`](https://github.com/stblake/mathilda/blob/main/src/simp/trigsimp.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= TrigToExp[Sin[x]]
Out[1]= (-1/2*I) E^(I x) + (1/2*I) E^(-I x)
```

```mathematica
In[1]:= TrigToExp[Cos[x]]
Out[1]= 1/2 E^(-I x) + 1/2 E^(I x)
```

Hyperbolic functions become real exponentials:

```mathematica
In[1]:= TrigToExp[Cosh[x]]
Out[1]= 1/2 E^x + 1/2 E^(-x)
```

The inverse functions are rewritten via their logarithmic forms — here the
classic `arctan` identity:

```mathematica
In[1]:= TrigToExp[ArcTan[x]]
Out[1]= (-1/2*I) Log[1 + I x] + (1/2*I) Log[1 - I x]
```

```mathematica
In[1]:= TrigToExp[ArcSin[x]]
Out[1]= -I Log[I x + Sqrt[1 - x^2]]
```
