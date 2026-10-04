---
references:
  - "G. E. Collins and H. Hong, *Partial Cylindrical Algebraic Decomposition for Quantifier Elimination*, J. Symbolic Computation **12** (1991) 299-328."
  - "S. McCallum, *An Improved Projection Operation for Cylindrical Algebraic Decomposition*, in *Quantifier Elimination and Cylindrical Algebraic Decomposition* (Springer, 1998) 242-268."
  - "S. Basu, R. Pollack and M.-F. Roy, *Algorithms in Real Algebraic Geometry*, 2nd ed. (Springer, 2006), Ch. 5 & 11 (sign determination, CAD)."
  - "A. Schrijver, *Theory of Linear and Integer Programming* (Wiley, 1986), §12.2 (Fourier-Motzkin elimination)."
source: src/solve/reduce.c
---
**Algorithm.** `builtin_reduce` wraps `reduce_impl` in an internal-function
message mute (`mth_msg_ifun_suppress`), so the speculative probes the engine
runs internally do not leak diagnostics. `reduce_impl` peels trailing option
`Rule`s, validates the variable spec, and rewrites a `List` in the `expr` slot
into an `And`. Two heads short-circuit the pipeline: a top-level `Exists`/`ForAll`
routes to `reduce_qe_dispatch` (quantifier elimination), and a `Modulus -> p`
option routes to `reduce_modular` (residue enumeration over `Z/pZ`). Otherwise the
statement is normalised to a DNF `RForm` by `reduce_form_from_expr` and simplified
(`rform_simplify`); the domain is chosen (default `Complexes`, but `Reals` when an
ordering `<`/`<=` is present with no explicit domain, or when a univariate
real-function / multivariate piecewise/radical statement is force-routed and
pre-rewritten by `reduce_realfn_preprocess` / `reduce_piecewise_preprocess`), and
the per-domain dispatch runs.

The dispatch is: **Integers/Rationals** → `reduce_integers`, which reuses the
`Solve` Diophantine engine and reformats its answer as an `Or` of `And`s carrying
`Element[C[k], dom]` for each free parameter. **Complexes** (equations only) → a
single univariate polynomial equation goes to `reduce_eq_univariate` (which carries
the full leading-coefficient case tree) or, when the residual is transcendental, to
`reduce_eq_transcendental` via `Solve`; a linear system to `reduce_eq_system`; and
anything else (nonlinear, or equations plus `!=`) to the zero-dimensional engine
`reduce_zerodim`. **Reals, one variable** → the exact real-algebraic sign-diagram
engines `reduce_univar` / `reduce_univar_general`, with trig/hyperbolic/exponential
equations pre-empted through `Solve` (the sign diagram is unsound on them) and
periodic-equation-plus-bounds shapes handled by `reduce_periodic_region` /
`reduce_trig_ineq_region`. **Reals, several variables** → `reduce_fm`
(Fourier-Motzkin for linear systems), falling through to `reduce_cad` (Cylindrical
Algebraic Decomposition, McCallum projection) and finally `reduce_zerodim` for a
zero-dimensional system with irrational fibres.

**Data structures.** The working form is an `RForm`: a DNF of `RConj`
conjunctions of `RAtom`s, each a polynomial (an `Expr`) paired with a relation code
(`R_EQ`, `R_LT`, `R_LE`, `R_ELEM`). The CAD (`reduce_cad.c`) keeps per-level
`PolySet` projection stacks and a `CADRegion` cell tree, and delegates every
numeric step — `Discriminant`, `Resultant`, `FactorList`, `Solve[..,Reals]`
real-root isolation, and the exact `qqbar` sign oracle — to the evaluator and
`reduce_real_util`.

**Complexity / limits.** CAD cost is doubly exponential in the number of
variables (McCallum projection then lifting); the 2-variable case is fully
realised and the n-variable path is an iterated projection (`cad_project_out` for
each level). The governing invariant is **soundness over completeness**: an
undecidable sign or ordering (`qqbar` returning `-2`), an unwired case, or an
irrational section all return `NULL`, leaving the input unevaluated rather than
emitting a wrong or partial formula. Consequently a returned formula describes the
entire solution set exactly, and a literal `False` is a genuine proof of emptiness.
