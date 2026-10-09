# Task: POSSIBLE_ZEROQ_IMPROVEMENTS.md — trim + implement #5

## Part 1 — Trim the backlog file  ✅ DONE
- [x] Retire resolved sections (#1, #3, #4, #2's PZQ subset) into a short
      "Resolved & removed" ledger; keep #2's residual Simplify note and #5.

## Part 2 — Implement #5 (deep-cancellation screen false-negative)

### Diagnosis (confirmed empirically against HEAD)
`PossibleZeroQ[Cosh[14]^30 - (1 - Tanh[14]^2)^(-15)]` is identically 0 but HEAD
returns **False** (wrong). Root cause: at a huge operand scale the machine
`(1 - near-1)` subtraction inside Tanh loses >52 bits, so the rung-0 residual
(5.3e170) is a large fraction of scale (2.36e173, ratio ~2.3e-3 > the 2^-12
`ZT_OBVIOUS_NONZERO_BITS` gate). `decide_numeric` (and `screen_point` in the SZ
sampler) **trust that machine FALSE directly** and never climb — a wrong answer.

Key empirical fact: a cancellation-zero and a genuine small non-zero (`...+1`)
have IDENTICAL residuals through 500 bits (the +1 is masked by the artifact);
they only diverge at the top (1000-bit) rung. So the fix must climb the FULL
ladder and trust the top rung's floor — NOT early-exit. A genuine huge SUM (no
cancellation) plateaus at rung-1 and stays False cheaply.

### Plan
- [ ] Add `#define ZT_HUGE_SCALE_BITS 40` (near `ZT_OBVIOUS_NONZERO_BITS`).
- [ ] Add `decide_numeric_huge_scale(e, mag0, scale)`: climb 200→500→1000 bits;
      return the rung's verdict on the first plateau (`m >= prev*0.5` — genuine
      value resolved); else trust the top rung's verdict after shrinking to the
      floor; lenient TRUE if MPFR is unavailable.
- [ ] In `decide_numeric`: when scale is huge AND rung-0 residual trips the
      obvious gate, divert to `decide_numeric_huge_scale` instead of the
      immediate FALSE. Normal-scale path byte-for-byte unchanged.
- [ ] In `screen_point`: at huge scale + obvious residual, return UNKNOWN (defer
      to `decide_numeric` via sz_trial's existing escalation) rather than FALSE.
- [ ] Regression test `tests/test_zero_test.c` Group 20: the pure cancellation
      zeros → True; `...+1` and the SUM → False; Weierstrass round-trip → True;
      simple doc forms unchanged; a stable-verdict check.
- [ ] Verify: rebuild; confirm reproducers flip True and controls stay False;
      run full test_zero_test suite; spot-check DSolve corpus / Weierstrass
      timing (perf of the extra climb); run `make check-messages`.
- [ ] Docs: mark #5 RESOLVED in POSSIBLE_ZEROQ_IMPROVEMENTS.md; changelog note
      in docs/spec/changelog/2026-10-05.md; bump version.h (+0.001); git tag.

## Review

Done (v0.324). All items complete.

- **Reproducer (new).** Found a *deterministic* repro the prior session lacked:
  symbol-free `PossibleZeroQ[Cosh[14]^30 - (1 - Tanh[14]^2)^(-15)]` (≡ 0) returned
  `False` on HEAD (residual 5.3e170, scale 2.36e173, ratio 2.3e-3 > the 2^-12
  gate). High powers of `(1-Tanh^2)^-k` widen the bad region so a fixed point
  trips it without hitting a pole.
- **Key insight.** A cancellation-zero and a genuine small non-zero (`…+1`) are
  bit-identical through 500 bits; they diverge only at 1000 bits. So the fix
  climbs the full ladder and trusts the trend/top-rung floor — never an
  early-exit. Genuine huge sums (no cancellation) plateau at rung 1 → `False`
  cheaply.
- **Fix.** `ZT_HUGE_SCALE_BITS` (2^40) gate; `decide_numeric_huge_scale` helper;
  diversion in `decide_numeric`; defer-to-ladder in `screen_point`. Normal-scale
  path byte-for-byte unchanged.
- **Verified.** Reproducers flip `False→True`; controls (`…+1`, huge sum) stay
  `False`; Weierstrass round-trip `True` (now robust, not luck); all Group 1–20
  zero_test tests pass; `make check-messages` clean; huge-scale climb sub-10 ms;
  §2.2.5 DSolve solve unaffected.
- **Docs.** #5 moved to the backlog's Resolved ledger (full section retired);
  v0.324 changelog; version bump + tag; two memory updates.
- **Outstanding in the backlog file:** only #2 (a `Simplify`-search deficiency,
  not a `PossibleZeroQ` one). The file is now a one-item backlog.
