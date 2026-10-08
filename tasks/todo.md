# Fix: Minimize on positive-dimensional (flat-valley) minima — v0.312

Target bug: `Minimize[(x y - 3)^2 + 1, {x, y}]` returns unevaluated; should be
`{1, {x -> 3, y -> 1}}`. Root cause: positive-dimensional stationary variety →
`Solve` returns a parametric branch → `mz_parse_points` rejects the whole
result → zero candidates → silent `NULL`.

## Steps
- [x] 0. Build; premise check → **plan premise was WRONG**: `Solve` returns `Solve::nsdim` *unevaluated*, not a parametric branch. Pivoted to the QE route (anticipated in the plan).
- [x] 1-3. (Superseded) Implemented `mz_qe_infimum` fallback instead of the salvage: `Reduce[ForAll[{vars}, f>=b],{b},Reals]` for the infimum + `FindInstance` witness, verified exactly. Wired in `mz_run` after `mz_exact_poly` declines (unconstrained).
- [x] 4. Target case → `{1, {x->-1, y->-3}}` (f@pt=1, on x y=3); Maximize mirror `{-1,...}`; line valley `(x+y-2)^2`→`{0,...}`; the renamed `_solves` case → `{-5/4,...}`.
- [x] 5. Saddle `x^2-y^2` → `{-Infinity,...}` (QE answers it); `x y` and 3-var decline (Reduce unevaluated — safe).
- [x] 6. Full `tests/test_minimize.c` green; renamed `_declines`→`_solves`, added `test_flat_valley_hyperbola`, `test_positive_dim_unbounded`.
- [x] 7. valgrind: QE path identical to 2+2 baseline (420/60/16 blocks) — zero new leaks.
- [x] 8. `src/version.h`→v0.312 ($VersionNumber=0.312); calculus.md + changelog 2026-10-05.md + memory updated. check-messages green. Commit+tag: pending user go-ahead.

## Follow-up (same v0.312): algebraic-infimum flat valleys
- Reported: `Minimize[(x y - 3)^4 - x y + 1, {x, y}]` still unevaluated. Infimum is
  an algebraic `Root` ≈ -2.4725. QE (`Reduce[ForAll...]`) gives it fine, but the
  v0.312 witness `FindInstance[f == Root, {x,y}, Reals]` returns unevaluated, and
  `Solve[f_slice == Root]` yields a NESTED Root (Root-with-Root-coeffs) Mathilda
  can't `N`/`Simplify`/`Element[_,Reals]` — so verification failed → decline.
- Fix: `mz_qe_witness` second strategy = pin all-but-one var to a trial const and
  recursively `mz_univar_poly` the univariate SLICE → clean single-Root minimiser
  (`{x -> Root[4#^3-36#^2+108#-109,1], y -> 1}`), verified `rru_sign_compare==0`.
  `mz_point_attains`/`mz_point_from_findinstance` helpers factored out.
- Verified: user query → `{Root[256#^3+1536#^2+3072#+2075,1], {x->Root[...], y->1}}`
  in 0.15s; sextic rational valley `{3,...}`; Maximize mirror; all prior cases;
  full `tests/test_minimize.c` green (+`test_flat_valley_algebraic`); valgrind
  identical to baseline (no new leaks); docs/changelog/memory updated.

## Review
Root cause: positive-dimensional (non-isolated) unconstrained minimum. `Solve[∇f==0]`
is positive-dimensional → `Solve::nsdim`, so the critical-point method collects no
candidate and the unconstrained path returned a silent NULL. Fix: a QE fallback
(`mz_qe_infimum`) reads the global infimum off `Reduce[ForAll[...]]` and realises a
verified witness with `FindInstance`. Additive, single-file (`minimize.c`) + the
`mz_run` wire-in; soundness preserved by the existing Reduce oracle + exact witness
verification. Bonus: unbounded-below cases (saddle) now report `-Infinity` instead of
declining. Scope: unconstrained, 2-var within the CAD's regime; constrained /
non-attained / 3-var remain safe declines.
