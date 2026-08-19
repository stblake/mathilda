# NSum on integer-only summands (Prime, etc.)

## Symptom
- `NSum[1/Prime[n], {n,1,Infinity}]` spews `Prime::intpp: ...Prime[16.1579]...`
  (sampling the summand *outside* the integers) and returns 1.67108.
- `NSum[1./Prime[n]^2, {n,1,Infinity}]` returns 0.448548, a poor approx of the
  true prime-zeta P(2)=0.4522474200...

## Root cause
`ns_choose_method` sees monotone-decreasing terms → picks **EulerMaclaurin**.
EM samples the summand at *continuous* real x (exp-sinh tail integral) and on a
*complex* contour (circle-DFT derivatives). `Prime[16.1579]` has no value:
- emits `Prime::intpp` per sample (the spew), and
- fails to numericalize, so EM bails and falls back to WynnEpsilon with only the
  default 15 head terms → weak extrapolation (0.4485 / 1.671).

The summand has **no real-analytic continuation** — EM is not just wasteful, it
is *invalid* here (its integral & derivative model does not exist).

## Fix (3 parts, all in src/numerical_calculus/nsum.c)
1. **Mute per-sample arithmetic messages.** Wrap `ns_eval_expr_at` (the single
   funnel every summand sample passes through) in
   `arith_warnings_mute_push/pop`. NSum probes at many trial points; those are
   internal and must not surface `Prime::intpp`/overflow messages (Wolfram's
   NSum doesn't). Kills the spew regardless of method.
2. **Detect non-continuable summands, keep EM away.** In `ns_build_profile`
   (non-black-box path only, so multidim adds no evals) probe the summand at two
   non-integer index points (imin+0.5, imin+0.3). If neither numericalizes →
   `prof.continuous=false`. `ns_choose_method` then never returns EulerMaclaurin
   for a non-continuous summand — routes monotone tails to WynnEpsilon.
3. **Give the discrete path enough terms.** For a non-continuous, not-user-
   pinned sum, raise the head-term floor (~100) and sequence length (~40) so
   Wynn extrapolates a genuinely small tail. P(2): 0.4485 → ~0.4522.

## Divergence detection — deliberately NOT a heuristic
Verified with the condensation ladder: divergent `1/(n ln n)` (ratio→0.941) is
numerically indistinguishable from convergent `1/n^1.1` (ratio→0.933) over any
finite sample. A divergence flag would false-positive on legitimate slow-
convergent series — which is why Mathematica also returns a finite value for
`NSum[1/Prime[n]]`. The sound-but-narrow test (flag only when condensation
terms do NOT →0, i.e. harmonic-like) does NOT catch P(1), so it is not added by
default. Instead: honesty — non-converged extrapolation already warns `ncvg`.

## Verify
- No `Prime::intpp` on either input.
- P(2) ≈ 0.4522 (≥3 digits); 1/n^2 still 1.64493; Log[1+1/n^2] EM path intact.
- tests/test_nsum.c still green; valgrind clean.

## Review (done 2026-08-20)
Shipped all three parts in `src/numerical_calculus/nsum.c`:
1. `ns_eval_expr_at` wrapped in `arith_warnings_mute_push/pop` → 0 `Prime::intpp`
   (was ~42).
2. `NsProfile.continuous` set by `ns_summand_is_continuous` (2 non-integer probes
   at imin+0.5, imin+0.3), gates the two EM returns in `ns_choose_method`. Runs
   only on the non-black-box path → multidim adds 0 evals.
3. Non-continuous, not-user-pinned → `nsum_terms` floor 100; machine-only
   `extra_terms` floor 24 (MPFR length left to its bit-scaled auto).
Results: P(2) 0.4485 → 0.452173; P(4) full accuracy; P(1) 2.30467 with honest
`NSum::accgl` (err ~1.5e-3). All continuous EM cases unchanged.
Divergence flag NOT added — condensation ratio 0.941 (div `1/(n ln n)`) vs 0.933
(conv `1/n^1.1`) proves the boundary is undecidable from finite samples.
Tests: `test_integer_only_summand` + Prime case in `test_memory_loop`; nsum(17)
& nprod(13) green; check-c99 clean. Docstring, spec doc, changelog updated.
Lesson → memory `project_nsum_integer_only_summand`.
