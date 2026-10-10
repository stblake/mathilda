# DSolve M65 — §2.2.38 corpus (Problems 3701–3800)

Goode & Annin 4th ed., 2nd-order-linear dominated (67 _linear, 18
_with_linear_symmetries, 7 _missing_x, 6 _missing_y, 5 _Emden, 4 _exact,
2 _homogeneous, 1 _Gegenbauer). Full-milestone wave (like M64).

## Tasks
- [x] Fetch upstream §2.1.38 HTML (browser UA) → scratchpad
- [x] Convert → `DSolve_test_status/DE_examples_2238.m` (100 scalar, 9 IVP, 0 sys)
- [x] `make check-corpus-indvar` green
- [x] Spot-check records (all clean; standard Goode & Annin equations)
- [x] Build `dsolve_corpus_tests`; baseline = 98/100, 0 FAIL (3746, 3764) — 2 identical runs
- [x] Register `dsolve_corpus_2_2_38_tests` gate (baseline 1)
- [x] Root-cause fix: Integrate linear-argument stage → 3746 fixed (98→99/100)
- [x] Integrate regression suite clean (intrep test updated — a correct improvement)
- [~] Re-run full dsolve corpus set — RUNNING (background)
- [x] Regenerate `reports/2.2.38.md` + `.tsv`; update `STATUS.md` (+ wave-history)
- [x] Gate baseline = 1 (final non-PASS = 3764, latency residue)
- [x] M65 block in `DSOLVE_PLAN.md`; changelog `2026-10-05.md`; calculus.md doc
- [x] `make check-c99` + `check-messages` green
- [~] Clean `make -j` (after regression; version.h touched)
- [x] Bump `src/version.h` → v0.341
- [ ] test_dsolve unit test
- [ ] commit + tag v0.341 (await user ok; on main branch)

## Review

**Result: §2.2.38 (Problems 3701–3800, Goode & Annin 4th ed.) — 98 → 99/100,
0 FAIL, 0 crash, 0 timeout.** One general `Integrate` fix.

**The fix (`src/calculus/integrate_linarg.c`):** a linear-argument substitution
stage run just before Weierstrass in the indefinite cascade. Weierstrass always
substitutes `Tan[x/2]`, so a scaled trig argument like `Sec[3x]^2` was multiple-
angle expanded into a degree-12 rational in `Tan[x/2]` (and DSolve spun on it).
The new stage pulls a shared non-trivial linear argument `a·x+b` (a a non-zero
number; kernel in a denominator) out via `u=a·x+b`, closes the bare-argument
integral through the recursive cascade, and scales by 1/a. It declines on bare
arguments and polynomial trig, so existing outputs are untouched.

**Diagnosis notes (for lessons):**
- The symptom was "DSolve 3746 times out"; the cause was two cascade slots deep in
  `Integrate` (Weierstrass grabbing a scaled trig integrand before the clean rules).
  Pinned `DSolve\`VariationOfParameters` + direct `Integrate[...]` probes localised it.
- Placing the new stage *before* Weierstrass (rather than rewriting Weierstrass)
  kept the blast radius to exactly one behaviour change across the whole
  integrate/risch/trig/simp/CRC suite — and that one was an *improvement*
  (`∫₀^∞ e^{−cx}J₀(ax)dx` now closes to a ConditionalExpression; `test_integrate_intrep`
  updated).
- `tests/CMakeLists.txt` uses an explicit source list for the common test lib, so a
  NEW src file (`integrate_linarg.c`) must be added there too, else the test binaries
  link a stale object (symptom: ld error but a prior binary still "PASS"es).

**Residue 1, honest:** 3764 solves correctly but cold DSolve ~12.6 s > 8 s wall
(latency, like §2.2.32's 3161/3164/3165).
