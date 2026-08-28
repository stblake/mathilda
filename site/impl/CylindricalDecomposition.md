---
references:
  - "Collins & Hong, \"Partial Cylindrical Algebraic Decomposition for Quantifier Elimination\", J. Symbolic Computation 12 (1991)."
source: src/solve/reduce_companions.c
---
**Algorithm.** `builtin_cylindrical_decomposition` is a thin, Reals-only companion to
`Reduce`: it does no decomposition of its own but forwards the problem to `Reduce`'s
real-domain engine, where the cylindrical algebraic decomposition actually lives. This
keeps a single implementation of the CAD/real-quantifier-elimination machinery and
guarantees the two heads agree on every real problem.

*Argument handling.* It first peels any trailing option `Rule`s (a symbol LHS, e.g.
`Modulus -> p`) off the end, leaving the positional arguments. It accepts `[expr, vars]`
or `[expr, vars, Reals]` -- since the domain is always the reals, an explicit `Reals` is
redundant but allowed, and any other third positional declines soundly (returns `NULL`,
leaving the input unevaluated). The variable specification is validated by
`cad_valid_vars`.

*Delegation.* It builds `Reduce[expr, vars, Reals, <peeled options...>]`, forwarding the
option rules verbatim so that all of `Reduce`'s options (`Modulus`, `Cubics`, `Quartics`,
`WorkingPrecision`, ...) are reused with no per-option logic here, and evaluates it under
`mth_msg_suppress_push` / `pop` so that any diagnostics emitted while `Reduce` probes
internally are not misattributed to `CylindricalDecomposition`.

*Sound decline.* If `Reduce` declines -- returning an unevaluated `Reduce[...]` -- the
companion returns `NULL` rather than echoing an inner `Reduce[...]` under its own head,
so `CylindricalDecomposition[...]` is simply left unevaluated. Otherwise `Reduce`'s
result (a quantifier-free `And`/`Or` formula, or a bare `True` / `False`) is returned
directly.

**Complexity / limits.** Entirely those of the delegate `Reduce` over the reals: CAD is
doubly exponential in the number of variables in the worst case, so the head is practical
for few variables and low-degree polynomials, and returns unevaluated where the real
decomposition cannot be produced exactly.
