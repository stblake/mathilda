---
source: src/poly/eliminate.c
references:
  - "T. Becker, V. Weispfenning, *Gröbner Bases* (Springer, 1993)."
  - "D. Cox, J. Little, D. O'Shea, *Ideals, Varieties, and Algorithms* (Springer)."
---
**Algorithm.** `builtin_eliminate` (in `src/poly/eliminate.c`) removes a set of variables from a system of equations by **Gröbner elimination**. It accepts `Eliminate[eqns, vars]` where `eqns` is a `List`/`And` of `lhs == rhs` equations. A pre-pass tries to handle simple invertible transcendental equations of the shape `f[poly] == const` (e.g. `Exp`, `Log`, trig) by a one-layer principal-branch rewrite, emitting an `Eliminate::ifun` diagnostic to warn that branches may be lost.

The algebraic path moves each equation to `lhs − rhs` form, collects all variables, and orders them so the variables to be eliminated form the leading block (the same elimination-block layout used by `GroebnerBasis`'s 3-arg form). It converts the polynomials to `GBPoly` via `gb_from_expr` and runs `gb_buchberger` under the lex/elimination order. By the elimination theorem, the basis polynomials free of the eliminated variables generate the elimination ideal; those survivors are re-presented as balanced `Equal[posPart, negPart]` equations and joined with `And` when there are several.

**Data structures.** Reuses the Gröbner subsystem's `GBPoly` (GMP `mpq_t` coefficients + row-major exponent matrix) and the Buchberger/Gebauer–Möller core in `groebner.c`. The memory contract is the standard builtin one — every early-return path frees temporaries and never frees `res`.

**Complexity / limits.** Inherits Buchberger's worst-case doubly-exponential cost (lex/elimination orders are the expensive case). `gb_from_expr` cannot atomise Power-headed main variables, so genuinely transcendental systems outside the simple-invertible pre-pass fall back to leaving the input unevaluated.
