# Four more infinite-Sum families + ipow fix

## Goal
Close four held infinite-`Sum` classes, each a general algorithm for its type
(honest bounds where full generality is impossible):
- `Sum[Sin[k]/k] -> (Pi-1)/2`  (trigonometric / Fourier)
- `Sum[HarmonicNumber[k,2]/k^3] -> 3 Zeta[2] Zeta[3] - 9/2 Zeta[5]`  (non-diagonal Euler)
- `Sum[HarmonicNumber[k]^2/k^2] -> 17/4 Zeta[4]`  (quadratic Euler)
- `Sum[(4k+1) Binomial[2k,k]^3/(-64)^k] -> 2/Pi`  (Ramanujan 1/Pi, table-backed)

## Work
- [x] Phase 0 — `ipow` overflow (`src/power.c`): `64^22` returned 0 (`__int128`
      wrap for perfect-power bases past 2^127). Keep `|b| <= INT64_MAX` each
      iteration; route overflow to `bigint_pow`.
- [x] Phase 1 — new `Sum`Trigonometric` (`src/sum/sum_trigonometric.c`):
      TrigReduce+Expand linearise -> `Im/Re PolyLog[s, e^{i a}]`; s=1 collapses
      to `(Pi-a)/2` / `-Log[2 Sin[a/2]]`. Wired in `src/sum/sum.c`.
- [x] Phase 2 — non-diagonal linear Euler, odd weight (`src/sum/sum_euler.c`):
      Borwein double-zeta `euler_Z_odd`, parity branch + reflection.
- [x] Phase 3 — quadratic Euler `H_k^2/k^q`, q=2..5 (`src/sum/sum_euler.c`):
      `is_harmonic_sq` detector + per-weight table; q>=6 (weight 8) NULL.
- [x] Phase 4 — Ramanujan 1/Pi (`src/sum/sum_hypergeometric.c` gate +
      `src/special_functions/hypergeopfq.c` VWP `4F3 -> 2/Pi`). Boundary
      passthrough gated on genuine convergence (z != 1 AND sigma < 0) and a real
      reduction, so `Sum[k]`/`Sum[(k+1)(-1)^k]` stay held.

## Review
- All four targets produce the exact closed forms; numeric cross-checks match
  (`NSum` for the Ramanujan series -> `0.63661977236758...` = `N[2/Pi,20]`).
- ipow: `64^22`, `64^30`, `(-64)^25` correct; `64^22 == 2^132`; near-limit
  `7^22`/`12^13` still exact (no over-eager bailout).
- Tests green: sum (64), sum_trigonometric (14, new), hypergeopfq,
  sum_alternating (12), power_corpus, sum_euler. No regressions.
- Fixed one regression introduced mid-work (`Sum[k] -> ComplexInfinity` from the
  relaxed pFq gate) via the sigma<0 & z!=1 convergence test.
- Valgrind (`--leak-check=full`) on the three new-path sums: no leak stack
  references any new function; residue is the documented pre-existing
  `builtin_times <- dispatch_def` ownership-leak class + macOS dyld baseline.
- Docs: `docs/spec/builtins/calculus.md` (Sum`Trigonometric + extended Euler),
  `special-functions.md` (VWP pFq), changelog `2026-06-29.md`.
