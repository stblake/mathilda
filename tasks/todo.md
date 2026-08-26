# Reduce Phase 8 — FindInstance (in progress)

Scope (user-approved): FindInstance over Complexes / Reals / Integers / Rationals /
Booleans. CylindricalDecomposition deferred.

- [x] SYM_FindInstance (3 sites) in sym_names.{h,c}
- [x] reduce_companions.c: builtin_find_instance + helpers (Reduce-walker + Solve-fallback
      + Booleans SAT via static to_dnf), verify-gated
- [x] reduce_companions.h: declare builtin_find_instance; update header comment
- [x] reduce_companions_init: register FindInstance + Protected + docstring
- [x] Build clean (-Wall -Wextra, gcc-16) — also cleared the file's 4 pre-existing
      warnings (3x discarded-qualifiers via xcopy wrapper, 1x alloc-size guard);
      make check-c99 clean
- [x] tests/test_reduce.c: test_find_instance group (21 asserts, FullForm-pinned)
- [x] Script check: instances verified across C/R/Z/Q/Booleans; {} for unsat;
      unevaluated for decline (fi1/fi2/fi3 scratch scripts)
- [x] Docs: solutions-of-equations.md ## FindInstance + changelog 2026-08-24.md
- [x] version.h 0.103 -> 0.104
- [x] reduce_tests green; reduce_corpus 158/158; valgrind byte-identical to macOS baseline
- [x] REDUCE_PLAN.md status updated (Phase 8: FindInstance done; CylindricalDecomposition left)

## Review

FindInstance shipped over Complexes / Reals / Integers / Rationals / Booleans.
Key architectural call: rather than exposing the three `static` CAD paths
(nu==1/2/>=3) to reach cell samples, witnesses are extracted from the PUBLIC
already-cylindrical outputs of Reduce and Solve and every candidate is VERIFIED
against the original expr (expr /. point === True). That single gate makes all the
heuristic point-picking sound: a bad pick fails verification and is skipped; {} is
emitted only when Reduce proves emptiness; unwitnessed-and-unrefuted declines.
Result: many cases match Mathematica exactly, all are true solutions, and it even
finds instances where Reduce declines a full reduction (Solve fallback). No new
low-level machinery — reuses Reduce, Solve, rru_rational_between, and the module's
own to_dnf (Booleans). Leak-clean, warning-clean, 158/158 corpus.

Follow-ups (out of scope, sound declines): Rationals *inequalities* decline
(Reduce[Rationals] handles only equations); `Element[x,Reals]` inside the statement
declines (Reduce doesn't parse membership in the statement); bare `x != 0` declines
(Reduce itself declines). CylindricalDecomposition is the last Phase-8 companion.

## Design notes

- Witness extraction from PUBLIC Reduce/Solve outputs (already cylindrical), not CAD
  internals (3 static paths). Every candidate verified: expr /. point === True.
- clause_point: params(Element)->0, free vars (not Equal-isolated) sampled first via
  rru_rational_between, then pin equations resolved by Solve[atom, s, dom].
- Solve fallback instantiates free listed vars -> 0 (grid for extra instances).
- Booleans: to_dnf clause -> partial assignment, unconstrained -> False, verify
  (Xor/Implies evaluate on Booleans since v0.102).
