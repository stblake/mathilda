---
references:
  - "M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005)."
  - "M. Bronstein, \"poor man's integrator\" (pmint), 2004 — parallel Risch / Risch–Norman heuristic."
  - "B. M. Trager, \"Algebraic Factoring and Rational Function Integration\", ACM SYMSAC 1976."
  - "K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992)."
  - "K. Roach, \"Symbolic Integration\", 1992 (undefined-function integrands, §1.7)."
  - "CRC Standard Mathematical Tables, 31st ed. — reference integral tables."
source: src/calculus/integrate.c
---
**Algorithm.** `Integrate[f, x]` computes an *indefinite, single-variable*
antiderivative only; there is no constant of integration and no definite-bound
support, and results are not guaranteed to be simplified. The dispatcher
`builtin_integrate` (src/calculus/integrate.c) routes through a fixed cascade of
context-qualified sub-integrators, each of which returns the antiderivative on
success or comes back as its own unevaluated head (e.g.
`Integrate\`BronsteinRational[...]`) to signal "fall through". The Automatic
cascade order is:

1. **Undefined-function integrator** (`try_undefined` → src/calculus/integrate_unknown.c),
   Roach 1992 §1.7: integrands rational in unknown `u[x]` and their derivatives
   (e.g. `x f'[x] + f[x] -> x f[x]`). Cheaply gated.
2. **Bronstein rational integration** (`try_rational` → src/calculus/intrat.c),
   gated by `PolynomialQ`/rational-function tests. This is the mathematically
   substantial stage: polynomial part by term-wise power rule, then **Mack's
   linear Hermite reduction** (`HermiteReduce`) to split off the rational
   (algebraic) part of the antiderivative, then the **Lazard–Rioboo–Trager
   log part** (`IntRationalLogPart`) built on a subresultant PRS, with Rioboo's
   recursive `LogToAtan`/`LogToReal` conversions producing real-elementary
   `Log`/`ArcTan` forms (a `RootSum`-based `NaiveLogPart` is the universal
   fallback). A pre-Hermite derivative-recognition fast path catches `c·D'/D^k`.
3. **Radical substitutions**: linear-radical (`try_linrad`), quadratic-radical
   Euler substitution (`try_quadrad`), and linear-ratio Möbius substitution
   (`try_linratiorad`).
4. **Derivative-divides** substitution (`try_derivdivides`): recognises
   `c·h(u(x))·u'(x)` and reduces to `Integrate[h[u], u]`.
5. **Risch–Norman heuristic** (`try_risch` → src/calculus/intrischnorman.c),
   Bronstein's *pmint*: the parallel-Risch ansatz for transcendental
   integrands (exp/log/trig via `tan` rewriting). Builds a candidate monomial
   set, a vector field with `splitFactor`/deflation, a linear system solved by
   `RowReduce`, plus a log-candidate sum; retries over `K = I`.
6. **CRC integral table** (`try_crctable`): a large rule set lazily
   `Get`-loaded from src/internal/CRCMathTablesIntegrals.m on first use only.
   Recursion-capped at `MAX_CRC_DEPTH` (256); a leaked internal
   `IntegrateTable[...]` head anywhere in the result counts as a miss.

An explicit `Method -> "<name>"` option (parsed by `parse_method_option`)
bypasses the cascade and dispatches to one stage strictly (no fallback).
Inexact integrands are force-rationalised by `common_rationalize_input`
(min-precision-aware, ½-ulp fallback), integrated exactly, then numericalised
back via `common_numericalize_result` for inexact-in/inexact-out semantics.
Applied 1-D `InterpolatingFunction` objects integrate to a fresh antiderivative
`InterpolatingFunction` before any rationalisation (`integrate_interp`).

**Data structures.** Everything is a Mathilda `Expr*` tree; the rational stage
leans on Mathilda's own polynomial builtins (`Together`, `Numerator`,
`Denominator`, `SubresultantPolynomialRemainders`, `PolynomialGCD`). Stages
communicate by evaluating context-qualified heads (`call_stage`) and
distinguishing success from passthrough by inspecting the result head
(`result_is_unresolved`, `result_contains_head`).

**Complexity / limits.** Indefinite, single-variable only. The rational stage
is complete for rational integrands (Hermite + LRT); transcendental coverage is
heuristic (pmint may give up) backed by a finite table. No definite integrals,
no multivariate integration, no constant of integration.
