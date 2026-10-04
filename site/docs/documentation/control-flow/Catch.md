# Catch

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Catch[expr]`**

Returns the argument of the first Throw generated while evaluating expr, or expr if none. Catch\[expr, form\] catches only a Throw\[value, tag\] whose tag matches form (tag is re-evaluated per comparison); Catch\[expr, form, f\] returns f\[value, tag\].

## Examples (11)

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

### Applications (6)

The first i whose square exceeds 20

```mathematica
In[6]:= Catch[Do[If[i^2 > 20, Throw[i]], {i, 1, 10}]]
Out[6]= 5
```

Evaluation stops at the Throw, so d is never reached

```mathematica
In[7]:= Catch[a; b; Throw[c]; d]
Out[7]= c
```

A tagged throw caught by a matching form

```mathematica
In[8]:= Catch[Throw[x, u], u]
Out[8]= x
```

The three-argument form returns f[value, tag]

```mathematica
In[9]:= Catch[Throw[v, tg], tg, f]
Out[9]= f[v, tg]
```

A non-matching tag is left uncaught

```mathematica
In[10]:= Catch[Throw[val, mytag], othertag]
Out[10]= Hold[Throw[val, mytag]]
```

With no throw, Catch is just the value of its body

```mathematica
In[11]:= Catch[1 + 1]
Out[11]= 2
```

## Implementation notes

**Algorithm.** `Catch` is `HoldFirst, Protected`, so it drives evaluation of its
own body and can intercept a throw. `builtin_catch` calls `evaluate` on the held
first argument; if the result is not an in-flight `Throw` sentinel
(`eval_is_inflight_throw`) it is returned verbatim, so `Catch[expr]` with no throw
yields `expr`'s value. When a sentinel comes back, the 1-argument form returns a
copy of the thrown value. The 2- and 3-argument forms are eligible only for a
*tagged* throw: the tag is re-evaluated (Wolfram semantics) and matched against
`form` through the pattern matcher (`match`); a non-match re-returns the same
sentinel so it propagates to an outer `Catch`, and a tagless `Throw[value]` is
never caught here. The 3-argument form returns the evaluated `f[value, tag]`.

**Data structures.** The in-flight marker *is* an ordinary `Throw[...]` node
(arity 1–3), carried up the evaluator's normal return path rather than by
`setjmp`/`longjmp` — so every intervening frame runs its own cleanup and the
construct is leak-free (proven byte-identical to a no-throw control loop under
valgrind). `env_new`/`match`/`env_free` back the tag comparison.

**Complexity / limits.** The first `Throw` evaluated wins. `evaluate_step`'s
argument loop short-circuits on a sentinel, so a throw deep inside `Plus`,
`Times`, `Map`, `Sum`, `Table` or a function application still reaches the nearest
`Catch`; a few consuming sites (`Which`/`Switch`, `Scan`, `SelectFirst`, the
`iter_run` family) carry an explicit propagation check the arg loop cannot
backstop. An uncaught throw is handled by `eval_report_uncaught_throw`.

- `Throw` is `Protected`; `Catch` is `HoldFirst, Protected` (it drives evaluation of its body itself, so it can intercept a throw; `form` and `f` evaluate normally).
- Implemented by sentinel propagation through the evaluator's normal return paths (no `setjmp`/`longjmp`), so every frame runs its own cleanup — leak-free.
- The first `Throw` evaluated wins; a tagless `Throw[value]` is not caught by a form-`Catch`.
- An uncaught `Throw[value]`/`Throw[value, tag]` returns `Hold[Throw[...]]` with a `Throw::nocatch` message; an uncaught `Throw[value, tag, f]` returns `f[value, tag]`.

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [Throw](../../control-flow/Throw/), [Return](../../control-flow/Return/), [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Map](../../data-structures/Map/), [Sum](../../calculus/Sum/), [Table](../../lists-and-iteration/Table/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_catch_throw.c`](https://github.com/stblake/mathilda/blob/main/tests/test_catch_throw.c)
- Tests: [`tests/test_core_algebra.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core_algebra.c)
- Tests: [`tests/test_scan.c`](https://github.com/stblake/mathilda/blob/main/tests/test_scan.c)
- Tests: [`tests/test_sow_reap.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sow_reap.c)

## Notes & additional examples

### Notes

`Catch[expr]` catches the first `Throw` generated anywhere while `expr` is
evaluated and returns its value, or returns `expr`'s own value if nothing is
thrown. Because the throw propagates through *any* enclosing head — `Plus`,
`Map`, `Table`, a function call — `Catch` is the tool for a non-local exit, where
`Return` only escapes a scope boundary.

`Catch[expr, form]` catches only a `Throw[value, tag]` whose `tag` matches `form`
(the tag is re-evaluated before each comparison); a throw with a non-matching tag,
and any tagless `Throw[value]`, propagates to an outer `Catch`. `Catch[expr, form,
f]` returns `f[value, tag]`. The construct uses sentinel propagation rather than
`setjmp`/`longjmp`, so it is leak-free.
