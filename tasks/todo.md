# Simplify log-exp-under-Assumptions campaign (v0.331+)

Goal: close the assumption-gated log/exp/trig gaps in the 50-expression stress corpus.
Baseline 40/50. Tier-1 milestone 48/50, Tier-2 stretch 50/50. Fully general + sound;
no `PossibleZeroQ` (0s come from exact transforms). Plan:
`~/.claude/plans/similar-to-previous-campaigns-replicated-unicorn.md`.

## Done — Tier 1 (40/50 → 48/50), one coherent commit, v0.331
- [x] T3a: `Log[x_ + I Sqrt[1-x_^2]] :> I ArcCos[x]` in `exp_to_trig_rules`
      (`src/simp/trigsimp.c`). Unconditional, sound. Closes #32, #39.
- [x] T1: `sr_proves_sle` helper + `ArcSin[Sin]`, `ArcCos[Cos]`, `ArcCot[Cot]`
      range-gated collapse (`simp_assume_rewrite.c`). Closes #27, #28.
- [x] T2: `ArcCosh[u] -> Log[u + Sqrt[u^2-1]]` gated u≥1 (`simp_assume_rewrite.c`).
      Closes #35. (Combined radical unsound off [1,∞) → gate essential.)
- [x] T4: Exp peel `E^(2 I Pi n + rest) :> E^rest` (+ odd `I Pi n` analogue),
      integer bucket. Closes #45.
- [x] T5: `Power[-1, Plus[k_Integer, n]] :> Power[-1,k] Power[-1,n]`, integer
      bucket. Closes #50.
- [x] T6: `Power[E^w, r] -> E^(r w)` gated `Im[w] ∈ (−π,π]` via Reduce oracle
      (`simp_assume_rewrite.c`). Closes #25.

## Tier 2 — investigated, dropped to documented known-gap (no clean sound fix)
- [~] T7 (#37) `Log[(1+Ix)/(1-Ix)] = 2I ArcTan[x]`: a search-adoption problem, not
      a missing identity — `TrigReduce` already produces the split `Log[1+Ix] -
      Log[1-Ix]` but the search won't adopt the higher-complexity split to let the
      ArcTan rule fold it. Fix needs a seed/ordering change (broad blast radius).
- [~] T8 (#30) `Log[1+E^(Ix)]` half-angle: a pattern-specific trusted rewrite,
      not a general mechanism. Both documented in the v0.331 changelog.

## Verification — all green
- [x] corpus `scratchpad/logexp_stress.m` → 48/50
- [x] 17 new unit + soundness-control tests in `tests/test_logexp_simplify.c` pass
- [x] 33 Simplify-adjacent suites pass (simp/fullsimplify/radical/trig/normalize/
      series/refine/assuming/powerexpand/complexexpand/possiblezeroq/… ) — no regressions
- [x] `make check-messages` green (BASELINE empty)
- [x] valgrind: definitely/indirectly-lost identical to startup baseline → zero new leaks
- [x] build clean under gcc-16 `-std=c99 -Wall -Wextra -Werror=...`

## Review
Eight exact, assumption-aware rewrites closed the tractable 8 of 10 corpus gaps,
all in two files (`trigsimp.c` ExpToTrig table + `simp_assume_rewrite.c`). Design
hinges confirmed by probing: Simplify's 0s come from exact transforms only (no
PossibleZeroQ); every gated rule is paired with a soundness control (inert without
its assumption, refuses the wrong region). The two unconditional rules (ArcCos log
form, (−1)^(k+n) split) are principal-value general and verified numerically. The
two remaining cases are a search-adoption issue (#37) and a pattern-specific
identity (#30), deliberately left as documented known gaps rather than shipped as
fragile/risky rewrites. v0.330 → v0.331.
