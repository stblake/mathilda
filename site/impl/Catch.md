---
source: src/funcprog.c
---
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
