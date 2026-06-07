---
source: src/calculus/deriv.c
---
**Algorithm.** `Derivative[n1, ..., nm]` is primarily a *tag* head: the actual
differentiation work is done inside the `D` dispatch (src/calculus/deriv.c).
`builtin_derivative` itself returns NULL, leaving `Derivative[n]` in canonical
unevaluated form; the builtin exists chiefly so attributes can be registered on
the symbol. Two pieces of real logic apply the tag:

- `derivative_of_pure_function(deriv_head, pure_fn)` differentiates
  `Derivative[n1,...,nm][Function[{t1,...,tm}, body]]` by partial-differentiating
  the body `ni` times in each slot `ti` via the shared `compute_deriv` core.
- `derivative_of_symbol(deriv_head, fsym)` reduces `Derivative[...][f]` when the
  symbol `f` carries DownValues: it mints fresh temporary slot symbols, builds
  and evaluates `f[t1,...,tm]` (triggering the DownValue rewrite), wraps the
  substituted body in a synthetic `Function`, and delegates to
  `derivative_of_pure_function`. If the call did not rewrite (no matching
  DownValue) it aborts to NULL. All `ni` must be nonnegative integers.

For unknown functions, `compute_deriv`'s chain rule emits `Derivative[...][f]`
factors, so the tag composes naturally through the rest of differentiation.

**Data structures.** `Expr*` trees; uses a static counter to generate
collision-free temporary slot-variable names (`Derivative$<id>$<k>`) without
registering them in the symbol table.

**Complexity / limits.** Only nonnegative integer derivative orders are
reduced; symbolic or negative orders stay unevaluated.
