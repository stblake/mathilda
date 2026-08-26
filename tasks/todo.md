# Reduce Phase 8 — LogicalExpand  (DONE)

- [x] SYM_LogicalExpand + SYM_NotElement (3 sites each) in sym_names.{h,c}
- [x] New src/solve/reduce_companions.{c,h}: standalone Expr-level DNF engine + builtin_logical_expand + builtin_not_element + reduce_companions_init
- [x] Wire reduce_companions_init() into reduce_init() (reduce.c)
- [x] Attributes set inline (Protected) in reduce_companions_init (Reduce family convention; no attr.c change)
- [x] Add reduce_companions.c to tests/CMakeLists.txt COMMON_SRC
- [x] test_logical_expand group in tests/test_reduce.c (23 assertions, FullForm-pinned)
- [x] Build clean (-Wall -Wextra, gcc-16); make check-c99 clean
- [x] Reproduce every spec example; all match (bar the one documented consensus-minimization case)
- [x] reduce_tests green (all pass); reduce_corpus_tests 158/158; valgrind byte-identical to baseline (no LE-attributable leak)
- [x] Docs: solutions-of-equations.md (## LogicalExpand + ## NotElement) + changelog 2026-08-24.md; version -> 0.101
- [~] Corpus le-* rows: SKIPPED — the reduce corpus oracle is Reduce-specific (calls Reduce, samples numeric grids); it does not fit LogicalExpand's opaque-Boolean inputs. Verification is via C-side FullForm assertions instead.
- [~] Audit EXEMPT: NOT NEEDED — Reduce/LogicalExpand/NotElement are symbolic heads with no numeric fast path; check-packed-aware passes and neither head is flagged.

## Review

Implemented `LogicalExpand` as a standalone Expr-level DNF distributor
(src/solve/reduce_companions.c) — NOT reusing RForm/RAtom, because those force
polynomial-relational leaves while LogicalExpand must treat symbols / x==a /
Element as OPAQUE Boolean atoms (no domain reasoning), exactly like Mathematica.
A mutually-recursive to_dnf / to_dnf_neg handles And/Or/Not/Implies/Xor and
Element-over-Alternatives; negation folds into the complementary relation head
(== <-> !=, < <-> >=, Element <-> NotElement, a <-> !a). True/False collapse is
sound AND complete with no truth-table enumeration: a DNF over independent atoms
is UNSAT iff it distributes to zero surviving clauses, so phi empty => False and
DNF(!phi) empty => True.

Notes / follow-ups discovered (out of scope, flagged for the user):
- Xor and Implies are INERT on literal Booleans in Mathilda (no eval builtin), so
  the user's `Table[e1==e2,...]` equivalence check does not reduce. A small,
  in-family follow-up (make Xor/Implies evaluate on True/False) would close this.
- Output is a sound, equivalent DNF but not always minimal: Mathematica's extra
  consensus contractions can produce a shorter, logically identical cover (only
  the e1=Implies[Xor[a,b,c],(a||b)&&c] example differs). Matches the rest of the
  Reduce family's "correct, not necessarily minimal" convention.
- Pre-existing, unrelated: `make check-compile-coverage` fails on Image*/
  Interpolation* heads (not touched by this change; LogicalExpand/NotElement are
  not flagged).
- Remaining Phase 8: FindInstance, CylindricalDecomposition (need a public seam
  over the static cad_build/CADRegion tree + reduce_zerodim_solve witnesses).

## Follow-up (done): Xor / Implies evaluate on Booleans  (v0.102)

- [x] builtin_xor + builtin_implies in src/boolean.c; registered in boolean_init with docstrings
- [x] attr.c: Xor = Flat|Orderless|OneIdentity|Protected; Implies = Protected
- [x] Xor folds Booleans, cancels duplicate args (a Xor a = False), Xor[]=False, Xor[e]=e,
      Xor[True,a]=Not[a]; Implies[False,_]=True, Implies[True,q]=q, Implies[p,False]=Not[p], Implies[p,p]=True
- [x] test_xor_evaluation + test_implies_evaluation in tests/test_boolean.c
- [x] Docs: control-flow.md (## Xor, ## Implies) + changelog; version -> 0.102
- [x] Regression: reduce_tests, corpus 158/158, boolean/comparisons/compile/core/patterns/fullsimplify all green; check-c99 clean
- [x] The user's Table[e1 == LogicalExpand[e1], ...] equivalence check now reduces to {True} for all 8 assignments
