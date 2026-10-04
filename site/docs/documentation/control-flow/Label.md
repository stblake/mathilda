# Label

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Label[tag]`**

Marks a point in a CompoundExpression to which control can be transferred with Goto\[tag\]. As a statement it evaluates to Null.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Module[{i = 0, s = 0}, Label[top]; i = i + 1; s = s + i; If[i < 5, Goto[top]]; s]
Out[1]= 15

In[2]:= f[a_] := Module[{x = 1., xp}, Label[begin]; If[Abs[xp - x] < 10^-8, Goto[end]]; xp = x; x = (x + a/x)/2; Goto[begin]; Label[end]; x]; f[2]
Out[2]= 1.41421
```

### Applications (2)

Label marks the target that Goto returns to

```mathematica
In[3]:= Module[{k = 1}, Label[a]; k = 2 k; If[k < 16, Goto[a]]; k]
Out[3]= 16
```

As a bare statement a Label evaluates to Null

```mathematica
In[4]:= Label[done]
```

## Implementation notes

**Algorithm.** `Label[tag]` marks a jump target inside a `CompoundExpression`.
`Label` is `Protected`. `builtin_label` validates arity (exactly one argument) and
returns `expr_new_symbol(SYM_Null)`, so evaluated as an ordinary statement a
`Label` is a no-op worth `Null`. Its real role is passive: the *raw held*
`Label[tag]` node, as it appears literally in the enclosing `CompoundExpression`,
is what `builtin_compoundexpression` scans for when it consumes a `Goto[tag]`
sentinel (see [Goto](Goto.md)). Tags are compared structurally (conventionally a
literal symbol or integer).

**Data structures.** None of its own; it is a plain two-node `EXPR_FUNCTION` that
`CompoundExpression` reads by position among its statements.

**Complexity / limits.** A `Label` is meaningful only as an explicit element of a
`CompoundExpression`; the matching `Goto` resolves it by a linear scan of that
compound expression's statements, then of enclosing ones. A `Label` reached in
normal top-to-bottom flow simply evaluates to `Null` and execution continues.

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

## References

**See also:** [Goto](../../control-flow/Goto/), [CompoundExpression](../../assignment-and-rules/CompoundExpression/), [Catch](../../control-flow/Catch/), [Throw](../../control-flow/Throw/), [If](../../control-flow/If/), [While](../../control-flow/While/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_goto_label.c`](https://github.com/stblake/mathilda/blob/main/tests/test_goto_label.c)

## Notes & additional examples

### Notes

`Label[tag]` marks a point that `Goto[tag]` can jump to. It must appear as an
explicit element of a `CompoundExpression` — it is the literal `Label[tag]` node
that the compound expression scans for when it consumes a `Goto` sentinel.

Evaluated in ordinary top-to-bottom flow a `Label` is a no-op worth `Null`, so
reaching one by falling through simply continues to the next statement; only a
`Goto` gives it an effect.
