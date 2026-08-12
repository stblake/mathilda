# Throw

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Throw[value]`**

Stops evaluation and returns value to the nearest enclosing Catch. Throw\[value, tag\] is caught only by Catch\[expr, form\] whose form matches tag. Throw\[value, tag, f\] returns f\[value, tag\] if uncaught.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Catch[a; b; Throw[c]; d; e]
Out[1]= c

In[2]:= f[x_] := If[x > 10, Throw[overflow], x!]; Catch[f[2] + f[11]]
Out[2]= overflow

In[3]:= Catch[Do[If[i! > 10^10, Throw[i]], {i, 100}]]
Out[3]= 14

In[4]:= Catch[Throw[a, u], u]
Out[4]= a

In[5]:= Catch[Throw[v, tg], tg, {#1, #2} &]
Out[5]= {v, tg}
```

## Implementation notes

- `Throw` is `Protected`; `Catch` is `HoldFirst, Protected` (it drives evaluation of its body itself, so it can intercept a throw; `form` and `f` evaluate normally).
- Implemented by sentinel propagation through the evaluator's normal return paths (no `setjmp`/`longjmp`), so every frame runs its own cleanup — leak-free.
- The first `Throw` evaluated wins; a tagless `Throw[value]` is not caught by a form-`Catch`.
- An uncaught `Throw[value]`/`Throw[value, tag]` returns `Hold[Throw[...]]` with a `Throw::nocatch` message; an uncaught `Throw[value, tag, f]` returns `f[value, tag]`.

**Attributes:** `Protected`.

## See also

[Catch](../../control-flow/Catch/), [Return](../../control-flow/Return/), [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Map](../../data-structures/Map/), [Sum](../../calculus/Sum/), [Table](../../lists-and-iteration/Table/)

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_catch_throw.c`](https://github.com/stblake/mathilda/blob/main/tests/test_catch_throw.c)
- Tests: [`tests/test_core_algebra.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core_algebra.c)
- Tests: [`tests/test_fixedpoint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fixedpoint.c)
- Tests: [`tests/test_scan.c`](https://github.com/stblake/mathilda/blob/main/tests/test_scan.c)
