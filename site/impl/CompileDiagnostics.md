---
source: src/compile/compiled_function.c
---
**Algorithm.** `CompileDiagnostics[argspec, expr]` is `HoldAll, Protected` and
takes exactly what `Compile` takes. `builtin_compile_diagnostics` strips the same
trailing `WorkingPrecision`/`"BigIntegers"` options (so it can report whether the
*arbitrary-precision* subset lowers), parses the argspec, and calls the same
`compile_expr_prec` the compiler uses. On success it returns a rule list carrying
`"Compiled" -> True`, the `"ResultType"` (from `ct_name`: `Real`, `Integer`,
`Complex`, `Boolean`, or `MPFRReal`/`MPFRComplex`/`BigInteger` under the managed
options), the `"Instructions"` count, the `"CommonSubexpressions"` the optimiser
hoisted, and `"InstructionsUnoptimized"` — obtained by recompiling with
`COMPILE_NO_OPT`, because the optimiser rewrites in place and the pre-pass count
is otherwise gone. On failure it returns `"Compiled" -> False`, a `"Reason"`
(`compiled_bail_reason`) and the `"Subexpression"` (`compiled_bail_expr`) — the
*innermost* node that could not be lowered, not the construct that contains it.

**Data structures.** A small fixed `Expr*` array of `Rule[key, value]` pairs built
by the `diag_rule` helper and wrapped in a `List`; no compiled program survives
the call (both the optimised and unoptimised programs are freed after their counts
are read).

**Complexity / limits.** It exists because a bail is otherwise invisible: the
caller quietly interprets and the answer stays correct but an order of magnitude
slower, and the compilable subset is a *cliff* — one unsupported head costs the
whole body. A malformed argspec leaves the call unevaluated, exactly as `Compile`
does. `MATHILDA_COMPILE_DIAG=1` prints the same report to stderr whenever an
auto-compiled numeric builtin falls back.
