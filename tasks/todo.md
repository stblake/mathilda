# Issue #69 — Get line continuation + Reduce/Solve nonlinear systems

## Problem
1. **Line continuation `\` in `.m` files / `Get`** — a backslash immediately
   before a newline should join lines (Mathematica syntax). The lexer rejected
   it with "Unexpected character: '\'". **REGRESSION-class bug (real).**
2. **Reduce/Solve on a zero-dimensional nonlinear polynomial system with
   inequality constraints** returns unevaluated, e.g. the 3-circle system
   `u^2+v^2==9 && u^2+(a+v)^2==36 && (a+u)^2+v^2==25 && u>0 && v>0 && a>0`.
   Confirmed via a 0.089 rebuild that this is **NOT a regression** — it was never
   handled. The user wants a *general algorithmic* fix, for both Solve and Reduce.

## Root cause (part 2)
- `reduce.c` complexes branch: nonlinear equation systems fall to
  `reduce_eq_system`, which is LINEAR-only → declines.
- `reduce.c` reals branch: nonlinear systems fall to `reduce_cad`, which
  DECLINES on irrational fibre samples (`is_rational_number` gate). A
  zero-dimensional system with irrational algebraic solutions is exactly that.
- Solve's front-end does not handle a mix of equations + inequalities for
  nonlinear systems.

## Algorithm (general, not a hack)
Zero-dimensional decomposition + exact algebraic filtering:
1. Split a conjunct into equations E and side relations O (`<,<=,!=`).
2. Solve E over Complexes (reuse the existing polynomial-system solver) → a
   FINITE, fully-determined solution set (else decline: positive-dimensional).
3. For each solution branch, decide EXACTLY (FLINT qqbar oracle):
   - Reals domain: each coordinate is real (`flint_qqbar_is_real`), else drop.
   - Each side relation holds at the branch (`rru_sign_of` / `flint_qqbar_equal`),
     else drop; undecidable → decline (sound-over-complete).
4. Emit surviving branches: Reduce → Or of And(var==val); Solve → List of rules.

This is complete for zero-dimensional systems (finite solution set, enumerated
and filtered exactly) and sound (declines on anything undecidable).

## Tasks
- [x] Fix line continuation in `skip_whitespace` (src/parse.c).
- [x] Confirm part 2 is not a regression (0.089 worktree build).
- [ ] Add `flint_qqbar_is_real()` realness oracle (src/poly/flint_qqbar.{c,h}).
- [ ] New engine `src/solve/reduce_zerodim.{c,h}`: core branch solver + Reduce
      formatter + Solve list formatter.
- [ ] Wire `reduce_zerodim` into `reduce.c` (complexes + reals fallbacks).
- [ ] Wire Solve delegation for the mixed eqn+ineq nonlinear case.
- [ ] Tests: unit tests + the issue's exact system; regression sweep vs 0.089.
- [ ] Docs: docstring, docs/spec, changelog; book Index if needed.

## Verification
- Line continuation: `Get`/`-file` on the issue's reduce.m.
- `Reduce[3-circle system] → u==r && v==r && a==r` (branch 3, all positive).
- `Solve[3-circle system with ineqs] → {{u->..., v->..., a->...}}`.
- No regression on the reduce sweep (identical to 0.089 on prior-passing cases).

## Review (done)
All tasks complete. Summary of changes:
- `src/parse.c` `skip_whitespace`: `\`+newline (LF/CR/CRLF) is a line
  continuation (no statement break); a stray `\` still errors. Fixes the issue.
- `src/poly/flint_qqbar.{c,h}`: new `flint_qqbar_is_real()` realness oracle
  (1 real / 0 non-real / -1 undecidable).
- `src/solve/reduce_zerodim.{c,h}`: NEW shared engine. Solves zero-dimensional
  polynomial equation systems exactly and filters branches by side relations +
  realness using the qqbar oracle. `reduce_zerodim` (Reduce → Or/And) and
  `reduce_zerodim_solve` (Solve → rule-lists).
- `src/solve/reduce.c`: wired into Complexes (after linear) and Reals (after
  Fourier-Motzkin + CAD) branches; docstring refreshed.
- `src/solve/solve.c`: equations-with-constraints pre-pass (mirrors the
  Integers pre-pass); only fires when a side constraint is present.
- Tests: reduce_corpus 158/158, solve_corpus 99/99, reduce_tests/solve_tests
  green; new corpus + unit cases (incl. issue's system, complex-rejection→False).
  Two previously-"decline" corpus/unit cases promoted to solved (correct
  improvements). CMake COMMON_SRC updated.
- Docs: docs/spec/builtins/solutions-of-equations.md (Reduce + Solve),
  changelog 2026-08-24.md (both fixes).
- Verified: no regression vs 0.089 (byte-identical on prior-passing cases),
  `make check-c99` clean, valgrind no new Mathilda-frame leaks.

## Not in scope (pre-existing, not regressions)
- Positive-dimensional nonlinear systems over Reals with irrational fibres (e.g.
  `x^2+y^2+z^2==1 && x>0 && y>0 && z>0`) still decline — the n-var CAD's
  rational-sample limitation, unchanged by this work.
