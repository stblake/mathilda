# Label

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Label[tag]`**

Marks a point in a CompoundExpression to which control can be transferred with Goto\[tag\]. As a statement it evaluates to Null.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Module[{i = 0, s = 0}, Label[top]; i = i + 1; s = s + i; If[i < 5, Goto[top]]; s]
Out[1]= 15

In[2]:= f[a_] := Module[{x = 1., xp}, Label[begin]; If[Abs[xp - x] < 10^-8, Goto[end]]; xp = x; x = (x + a/x)/2; Goto[begin]; Label[end]; x]; f[2]
Out[2]= 1.41421
```

## Implementation notes

- Both are `Protected`. `tag` is evaluated (conventionally a literal symbol or
  integer) and compared structurally to each `Label`'s tag.
- Like `Catch`/`Throw`, `Goto` is implemented by sentinel propagation through the
  evaluator's normal return paths (no `setjmp`/`longjmp`), so a `Goto` fired
  inside a nested call (e.g. an `If` branch) still reaches the enclosing
  `CompoundExpression`. Leak-free.
- A `Goto` loop is a genuine loop with no artificial iteration cap; termination
  is the program's responsibility (as with `While`).
- A `Goto[tag]` that reaches the top level with no matching `Label` anywhere
  emits a `Goto::nolabel` message (stderr) and returns the inert `Goto[tag]`
  node. The message fires only when truly unmatched — a `Goto` that legitimately
  propagates from an inner to an outer `CompoundExpression` mid-evaluation is
  silent.

**Attributes:** `Protected`.

## See also

[Goto](../../control-flow/Goto/), [CompoundExpression](../../assignment-and-rules/CompoundExpression/), [Catch](../../control-flow/Catch/), [Throw](../../control-flow/Throw/), [If](../../control-flow/If/), [While](../../control-flow/While/)

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_goto_label.c`](https://github.com/stblake/mathilda/blob/main/tests/test_goto_label.c)
