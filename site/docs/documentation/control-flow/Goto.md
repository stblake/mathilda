# Goto

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Goto[tag]`**

Scans the CompoundExpression it appears in directly for Label\[tag\], then enclosing ones, and transfers control to that point.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Module[{i = 0, s = 0}, Label[top]; i = i + 1; s = s + i; If[i < 5, Goto[top]]; s]
Out[1]= 15

In[2]:= f[a_] := Module[{x = 1., xp}, Label[begin]; If[Abs[xp - x] < 10^-8, Goto[end]]; xp = x; x = (x + a/x)/2; Goto[begin]; Label[end]; x]; f[2]
Out[2]= 1.41421
```

### Applications (3)

A backward jump forms a loop

```mathematica
In[3]:= Module[{i = 0, s = 0}, Label[top]; i++; s += i; If[i < 5, Goto[top]]; s]
Out[3]= 15
```

Doubling until the guard fails

```mathematica
In[4]:= Module[{k = 1}, Label[a]; k = 2 k; If[k < 16, Goto[a]]; k]
Out[4]= 16
```

With no matching Label it is left inert

```mathematica
In[5]:= Goto[nowhere]
Out[5]= Goto[nowhere]
```

## Implementation notes

**Algorithm.** `Goto` is `Protected` with no `Hold` attributes, so `tag` is
evaluated by the argument loop. `builtin_goto` validates arity (exactly one
argument) and returns `NULL`, because — exactly as with `Throw` — the plain
`Goto[tag]` node is itself the in-flight sentinel. `eval_is_inflight_goto`
recognises it by head, and `evaluate_step`'s argument loop short-circuits on it,
so a `Goto` fired inside a nested call (an `If` branch, say) bubbles up to the
enclosing `CompoundExpression`. That `CompoundExpression` (`builtin_compound`
expression) is the construct that actually consumes the sentinel: it scans its own
statements for a raw held `Label[tag]` whose tag is structurally equal and resumes
evaluation there — a forward jump skipping intervening statements, or a backward
jump forming a loop. If the current compound expression has no matching `Label`,
the sentinel propagates to the enclosing one.

**Data structures.** The sentinel is the `Goto[tag]` node itself; targets are raw
held `Label[tag]` nodes scanned linearly. Propagation is by the normal return path
(no `setjmp`/`longjmp`), so it is leak-free.

**Complexity / limits.** A `Goto` loop is a genuine loop with no artificial
iteration cap; termination is the program's responsibility, as with `While`. A
`Goto[tag]` that reaches top level with no matching `Label` anywhere emits
`Goto::nolabel` (stderr) and returns the inert `Goto[tag]` node; the message fires
only when truly unmatched, so a `Goto` that legitimately propagates from an inner
to an outer `CompoundExpression` mid-evaluation stays silent.

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

**See also:** [Label](../../control-flow/Label/), [CompoundExpression](../../assignment-and-rules/CompoundExpression/), [Catch](../../control-flow/Catch/), [Throw](../../control-flow/Throw/), [If](../../control-flow/If/), [While](../../control-flow/While/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_goto_label.c`](https://github.com/stblake/mathilda/blob/main/tests/test_goto_label.c)

## Notes & additional examples

### Notes

`Goto[tag]` transfers control to the `Label[tag]` in the `CompoundExpression` the
`Goto` appears in directly, then in enclosing ones — a forward jump skips the
statements in between, a backward jump forms a loop. Like `Throw`, it propagates
by a sentinel through the evaluator's normal return path, so a `Goto` fired inside
a nested call (an `If` branch, as above) still reaches the enclosing
`CompoundExpression`.

A `Goto` loop has no artificial iteration cap; termination is the program's
responsibility, exactly as for `While`. A `Goto[tag]` that reaches top level with
no matching `Label` anywhere prints `Goto::nolabel` and returns the inert
`Goto[tag]` node.
