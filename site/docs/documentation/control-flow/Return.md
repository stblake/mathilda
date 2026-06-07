# Return

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Return[expr]
    returns the value expr from a function.
Return[]
    returns the value Null.
Return[expr, h]
    returns expr from the nearest enclosing construct whose head is h.

Return[expr] exits control structures within the definition of a function,
and gives the value expr for the whole function.
Return takes effect as soon as it is evaluated, even if it appears inside
other functions.

The recognised boundary heads are Function, Module, Block, With, Do, For,
and While. CompoundExpression and other Hold-free heads pass Return through
unchanged so that it can reach the boundary.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Function[x, If[x > 0, Return[positive], Return[negative]]][3]
Out[1]= positive

In[2]:= Module[{}, Do[Return[5], {3}]]
Out[2]= 5

In[3]:= Module[{}, Do[Return[5, Module], {3}]]
Out[3]= 5

In[4]:= Module[{}, Do[Return[5, Block], {3}]]
Out[4]= Return[5, Block]

In[5]:= f[n_] := Module[{s = 0}, Do[s = s + i; If[s > 10, Return[i]], {i, 1, n}]]; f[10]
Out[5]= 5
```

## Implementation notes

**Algorithm.** `Return` has no C builtin handler; it is a control-flow marker (just `ATTR_PROTECTED` in `attr.c`) recognized by `eval_classify_return` in `src/eval.c`. A `Return[...]` expression simply evaluates to itself and bubbles up through `CompoundExpression` and enclosing constructs unchanged until it reaches a boundary that inspects it. `eval_classify_return(e, boundary_head, &out)` uses pointer equality on the interned `SYM_Return` head and returns one of three actions: `CONSUME` with an out-value (the loop yields it), `PROPAGATE` (keep bubbling), or `NONE`. `Return[]` consumes with a fresh `Null`; `Return[expr]` consumes at the nearest boundary regardless of head, copying `expr`; `Return[expr, h]` consumes only at a boundary whose head symbol equals `h` (else PROPAGATE). The classifier is deliberately side-effect free — the single allocation on CONSUME is the copied payload, needed because the caller frees the original marker. The iteration loops (`Do`/`For`/`While`) call it via `iter_flow_classify`, and the body of user functions / `CompoundExpression` is the other principal boundary. Extra arguments beyond the second are ignored.

- `Protected`. No Hold attributes — arguments are evaluated before the marker takes effect.
- `Return` takes effect as soon as it is evaluated, even when it appears inside other functions (Plus, Times, etc. in argument positions still see the substituted value, but at the top level the marker propagates immediately through `CompoundExpression`).
- 1-arg `Return[expr]` is consumed by the *innermost* boundary on the call stack. The 2-arg `Return[expr, h]` form lets the user skip past intervening boundaries to a specific enclosing construct.
- If no enclosing boundary matches, `Return[expr]` (or `Return[expr, h]`) survives at the top level as a literal expression.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
