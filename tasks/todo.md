# Task: General log/exp transcendental solving in Solve & Reduce

Plan: /Users/user/.claude/plans/both-solve-and-reduce-gleaming-journal.md

## Steps

- [x] 1. Add `simp_log_fuse_all` (unconditional log fusion) to `src/simp/simp_log.c` + `.h`
- [x] 2. Splice fusion into `solveinv_drive` (`src/solve/solveinv.c:~1167`)
- [x] 3. Build + verify: all 20 log-family Solve cases solve
- [x] 4. Add `reduce_eq_transcendental` to `src/solve/reduce_eq.c` + `.h`
- [x] 5. Gate + call at `src/solve/reduce.c:325` (+ includes)
- [x] 6. Build + verify: log & exp families solve under Reduce
- [~] 7. Radical case #14 — NOT fixed: pre-existing `solverad` limitation on a
        rational-in-`Sqrt[t]` equation with a `Power`-valued constant RHS (some
        such inputs even *hang* solverad). Case #14 declines fast/safely.
        Documented as a known limitation.
- [x] 8. Tests: 31 back-sub cases in `tests/solve_corpus.m`; exact+count
        assertions in `test_solve.c`/`test_reduce.c`
- [x] 9. Regressions: solve/reduce/solve_corpus(130)/reduce_corpus(170)/dsolve/
        logexp_simplify/fullsimplify all green
- [x] 10. valgrind: no leak/error in new funcs (182 def-lost = macOS objc/dyld
         baseline); `make check-c99` clean
- [x] 11. Docs: weekly changelog + `solutions-of-equations.md` (Solve+Reduce);
         `tasks/lessons.md`; memory note. (Docstrings unchanged — terse, no
         capability enumeration to update.)

## Review

**Delivered.** Solve now inverts the log family (17/20 of the supplied cases;
the other 3 are documented edges) and the exp family (already worked, 19/20).
Reduce inverts both families via a new route. Core change is small and general:
one unconditional log-fuser (`simp_log_fuse_all`) spliced into `solveinv_drive`,
mirroring the exp gathering the evaluator already does; Reduce re-enters Solve
(`reduce_eq_transcendental`) and renders the rule-list as a formula.

**Solve status (40 supplied eqns):** log 1–13,16,17,19,20 solve; #14 radical
declines (solverad limit), #15 symbolic-exponent left unevaluated on purpose
(branch-preserving), #18 degenerate identity → unevaluated. exp 21–39 solve; #40
degenerate `1==C[1]`. **Reduce:** all log+exp solve; forward-trig and polynomials
unchanged.

**Files:** `src/simp/simp_log.{c,h}`, `src/solve/solveinv.c`,
`src/solve/reduce_eq.{c,h}`, `src/solve/reduce.c`; tests + docs + memory.

**Not done (deliberate):** solverad rational-radical fix for #14 (fragile/hangs;
out of scope, follow-up); no PowerExpand pre-pass for #15/#40 (would lose
solutions / never fire).
