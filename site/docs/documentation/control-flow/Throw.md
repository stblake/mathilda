# Throw

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Throw[value]`**

Stops evaluation and returns value to the nearest enclosing Catch. Throw\[value, tag\] is caught only by Catch\[expr, form\] whose form matches tag. Throw\[value, tag, f\] returns f\[value, tag\] if uncaught.

## Examples (8)

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

### Applications (3)

The thrown value is handed to the nearest Catch

```mathematica
In[6]:= Catch[Throw[5]]
Out[6]= 5
```

The throw unwinds out of the Plus

```mathematica
In[7]:= g[x_] := If[x > 10, Throw[overflow], x!]; Catch[g[2] + g[11]]
Out[7]= overflow
```

Uncaught at top level: reported and wrapped in Hold

```mathematica
In[8]:= Throw[orphan]
Out[8]= Hold[Throw[orphan]]
```

## Implementation notes

**Algorithm.** `Throw` is `Protected` with no `Hold` attributes, so `value`,
`tag` and `f` are evaluated by the ordinary argument loop before the throw
propagates. `builtin_throw` does almost nothing: it validates arity (1–3) and
returns `NULL`, because the plain `Throw[...]` node **is** the in-flight sentinel.
`eval_is_inflight_throw` recognises it by head (`SYM_Throw`, arity 1–3), and
`evaluate_step`'s argument-evaluation loop short-circuits when an evaluated
argument is such a sentinel — freeing the sibling arguments and the head and
returning the sentinel up the normal return path — so a `Throw` anywhere inside a
surrounding expression unwinds to the nearest enclosing `Catch`.

**Data structures.** No separate payload type: the sentinel is the user-visible
`Throw[value]` / `Throw[value, tag]` / `Throw[value, tag, f]` tree. Propagation is
by return value, never `longjmp`, so per-frame cleanup always runs (leak-free).

**Complexity / limits.** Unlike `Return`, which only escapes a scope boundary, a
`Throw` passes through *any* enclosing head. A tagged throw is caught only by a
`Catch[expr, form]` whose `form` matches the (re-evaluated) tag. If it reaches top
level uncaught, `eval_report_uncaught_throw` (`evaluate()`) emits `Throw::nocatch`
and returns `Hold[Throw[...]]`, except that an uncaught `Throw[value, tag, f]`
returns `f[value, tag]`.

- `Throw` is `Protected`; `Catch` is `HoldFirst, Protected` (it drives evaluation of its body itself, so it can intercept a throw; `form` and `f` evaluate normally).
- Implemented by sentinel propagation through the evaluator's normal return paths (no `setjmp`/`longjmp`), so every frame runs its own cleanup — leak-free.
- The first `Throw` evaluated wins; a tagless `Throw[value]` is not caught by a form-`Catch`.
- An uncaught `Throw[value]`/`Throw[value, tag]` returns `Hold[Throw[...]]` with a `Throw::nocatch` message; an uncaught `Throw[value, tag, f]` returns `f[value, tag]`.

**Attributes:** `Protected`.

## References

**See also:** [Catch](../../control-flow/Catch/), [Return](../../control-flow/Return/), [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Map](../../data-structures/Map/), [Sum](../../calculus/Sum/), [Table](../../lists-and-iteration/Table/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_catch_throw.c`](https://github.com/stblake/mathilda/blob/main/tests/test_catch_throw.c)
- Tests: [`tests/test_core_algebra.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core_algebra.c)
- Tests: [`tests/test_fixedpoint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fixedpoint.c)
- Tests: [`tests/test_scan.c`](https://github.com/stblake/mathilda/blob/main/tests/test_scan.c)

## Notes & additional examples

### Notes

`Throw[value]` stops evaluation and returns `value` to the nearest enclosing
`Catch`. Its arguments are evaluated first (`Throw` is not held), and the
`Throw[...]` node then propagates up through every intervening expression — the
second example shows it escaping a partially-evaluated `Plus`.

`Throw[value, tag]` is caught only by a `Catch[expr, form]` whose `form` matches
`tag`. An uncaught `Throw[value]` or `Throw[value, tag]` prints `Throw::nocatch`
and comes back as `Hold[Throw[...]]`; an uncaught `Throw[value, tag, f]` instead
returns `f[value, tag]`.
