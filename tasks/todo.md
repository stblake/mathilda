# P1 — Structure-theorem decision + logarithmic-derivative decision

Roadmap item P1 from `RISCH_BRONSTEIN_GAPS.md`. Bronstein §9.3 (Cor. 9.3.1/9.3.2),
§5.12 (In-Field Integration), §7.3 (Parametric Log Derivative). Substrate
(`risch_structure.c` complex case + `risch_rational_span`) already landed.
User-approved scope: **A + B + C-full (full oracle replacement)**. Extensive unit +
stress tests for each Bronstein routine. (Book page N ≈ PDF page N+17.)

## Phase A — Real structure case (Cor. 9.3.2 iii/iv)  [risch_structure.c]  ✅ DONE
- [x] Add `ArcTan` monomial kind to `rs_decode_tower`; expose per-monomial kind.
- [x] Partition generators by the four DISJOINT index sets E/L/T/A (E∪L for
      log/exp, T∪A for tan/arctan). Map coefficients back to full monomial order.
- [x] `Risch`TanReducible` (eq. 9.15) + `Risch`ArcTanReducible` (eq. 9.14).
- [x] Docstrings + attributes + changelog.
- [x] `tests/test_risch_structure_real.c` — all green (real reductions, genuine-new
      declines, disjoint-index, deeper real towers, back-compat, robustness).

## Phase B — LogarithmicDerivativeOfRadical (§5.12 / eq. 7.44)  [risch_structure.c]  ✅ DONE
- [x] `Risch`LogarithmicDerivativeOfRadical[f, x, mons]` -> `{n, u}` or `False`,
      structure-theorem test (Cor. 9.3.1/9.3.2 ii) via E∪L span + witness
      u = prod base^{n r}, D[u]/u == n f. (Full §5.12 recursion w/ log-monomial
      radicals = documented completeness item, pinned as sound decline.)
- [x] Docstrings/attrs/changelog.
- [x] `tests/test_risch_logderiv_radical.c` — all green (structural witnesses,
      end-to-end D[u]/u−n f==0 identity, declines, §5.12 scope boundary, stress).

BASELINE (pre-Phase-C, must not regress): integrate_risch_transcendental_tests &
integrals_tests CLEAN; intrat_tests 2 FAIL + intrischnorman_tests 6 FAIL are
PRE-EXISTING (files untouched).

## Phase C — Full oracle replacement in the live integrator  [integrate_risch_transcendental.c]
- [ ] Replace `rt_class_primitive` (exp commensurability) + `rt_expand_logs` (log
      syntactic split) as the PRIMARY independence engine in `rt_tower_build`:
      incremental structure-theorem tower construction — for each candidate
      log/exp kernel, decide independence vs the existing generators via the P1
      oracle and record the Q-relation; new => new generator.
- [ ] Catch additive Q-relations the syntactic path misses (Log[x^2+x]).
- [ ] Lift/relax RT_MAXK where independence is proven.
- [ ] Full transcendental + rational regression suite; diff vs baseline binary;
      wrong answer impossible (diff-back / tower-identity gate => decline not wrong).
- [ ] valgrind; docs; changelog.

## Verification (each phase)
- Clean `-std=c99 -Wall -Wextra`; scoped `*_tests` green (grep FAIL:); REPL pins;
  valgrind vs Sin[1.0] baseline. Commit per phase.

## Review
(to be filled in)
