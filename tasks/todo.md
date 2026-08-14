# Task: NMinimize SimulatedAnnealing Method options

## Goal
Make the SimulatedAnnealing method honor its three Method sub-options
("PerturbationScale", "BoltzmannExponent", "SearchPoints") — currently `nm_sa`
does `(void)nc;` and ignores all of them, and the parser silently drops
PerturbationScale/BoltzmannExponent. Add regression tests, including the
Griewank-10 case Mathematica returns a nonzero local minimum on.

## Plan
- [ ] Add `perturbation_scale` (double, <0 ⇒ default 1.0) and `boltzmann_fn`
      (Expr*, borrowed / NULL ⇒ default) to `NmConfig`; init in the setup block.
- [ ] Parse `"PerturbationScale"` and `"BoltzmannExponent"` in `nm_parse_method`.
- [ ] Add `nm_boltzmann_exponent()` helper (mirrors `nm_apply_penalty_fn`),
      calling `fn[i, df, f0]` → double, fallback to default exponent on failure.
- [ ] Rewrite `nm_sa` to:
      - run `K = search_points>0 ? search_points : 1` annealing chains from
        random starts, keeping the global best (SearchPoints = restarts);
      - scale the perturbation step by `perturbation_scale`;
      - use `Exp[boltzmann_fn[i, df, f0]]` as the Metropolis acceptance
        probability when a function is supplied, else the current `exp(-d/T)`.
      - **Preserve the K=1, default-options path bit-for-bit** (multiply by 1.0,
        identical RNG draw order) so existing SA tests stay deterministic.
- [ ] Bound total work when many chains are requested (per-chain iteration cap).
- [ ] Update docs: `docs/spec/builtins/calculus.md` sub-option table + paragraph;
      `src/info.c` docstring; changelog `docs/spec/changelog/2026-08-10.md`.
- [ ] Add regression tests in `tests/test_nminimize.c`:
      - Griewank-10 with the exact Method from the report (SA + PerturbationScale
        + BoltzmannExponent + SearchPoints) returns a finite nonzero local min;
      - each of the three options individually takes effect / is honored.
- [ ] Build, run test_nminimize, verify no regressions.

## Review
Done. `src/findmin.c`: added `perturb_scale` + `boltzmann_fn` to `NmConfig`
(init'd -1.0 / NULL), parsed `"PerturbationScale"` / `"BoltzmannExponent"` in
`nm_parse_method` (with `NMinimize::sopt` / `::bexp` fallbacks), added
`nm_boltzmann_exponent()` (mirrors `nm_apply_penalty_fn`), and rewrote `nm_sa`
to run K=`search_points` chains sharing a bounded budget (`NM_SA_TOTAL_CAP`),
scale the step by `perturb_scale`, and use `Exp[f[i,df,f0]]` acceptance. The
K=1/default-option path is bit-for-bit unchanged (verified: 1D SA still
−3.51391; default Griewank-10 still 0.233824, matching Mathematica).

Docs: `docs/spec/builtins/calculus.md` (sub-option table + SA paragraph),
`src/info.c` docstring, changelog `docs/spec/changelog/2026-08-10.md`.

Tests (`tests/test_nminimize.c`): `test_sa_suboptions` (each option honored +
invalid-value fallback) and `test_griewank_simulatedannealing` (default lands
in [0.1,1.0]; the exact reported all-options invocation stays finite/feasible).
Full nminimize_tests (52) + findmin_tests pass; `make check-c99` clean; suite
1.84s.

Caveat surfaced to user: the reported `"BoltzmannExponent" -> (1/#&)` is
always-accept (Exp[1/i] > 1), so under our schedule it random-walks to a higher
local min (46.4) than Mathematica's 0.234 — but the *default* SA matches MMA at
0.234. Options are demonstrably honored regardless.
