# Gruntz `Limit` — current state

Status snapshot as of 2026-07-25. The 2026-07-24 mrv engine + Phase 2/3 work is
committed (a309b74, d2c75a8). The 2026-07-25 additions below — ArcTan-at-oo,
and Gruntz coverage of `PolyGamma`, `Zeta`, modified Bessel `BesselK`/`BesselI`,
and `Max`/`Min` — are **uncommitted** in the working tree. Everything is green:
`gruntz_tests` pass; `limit`/`limit_assumptions`/`nlimit`/`nseries`/`residue`/
`series` show no new failures; valgrind == `1+1` baseline.

## 2026-07-25 session — new at-infinity coverage

- **ArcTan** at oo: `Series[ArcTan[x],{x,Infinity,n}]` + kernel composition
  (`Limit[x(Pi/2-ArcTan[x])]->1`).
- **PolyGamma[m, x]** (DLMF 5.11.2): digamma `Log[x]` head + Laurent; m>=1 pure
  `x^-m` Laurent. Series hook + 2-arg isolation.
- **Zeta[x]**: exp-log Dirichlet head `1+2^-x+3^-x+...`; isolation (`+oo` only).
- **BesselK/BesselI** (2-arg, any/symbolic order): monotonic `Exp[∓z]` envelopes.
  BesselJ/BesselY excluded (oscillatory).
- **Max/Min**: `resolve_maxmin` dominance pre-pass (leading-term sign of `a-b`);
  cascade `layer2_series`/`layer5_lhospital` now bail on Max/Min.

Details in the "Known gaps" section below and `docs/spec/changelog/2026-07-20.md`.

## What it is

A faithful C port of Dominik Gruntz's most-rapidly-varying (mrv) limit algorithm
(1996 ETH thesis; structure follows SymPy's `gruntz.py`), in
`src/calculus/gruntz.c` (~1589 lines) + `gruntz.h`. Single exported entry:

```c
Expr* gruntz_limit(Expr* f, Expr* x, Expr* point, int dir, int depth);
```

Reached two ways from `src/calculus/limit.c`:
- explicit `Limit[f, x->x0, Method -> "Gruntz"]` (string or bare-symbol method), and
- as a last-resort layer in the `Automatic` cascade (`LIMIT_M_GRUNTZ`,
  `layer_gruntz`, `TRY` after `layer2_series`).

Inputs borrowed; returns owned `Expr*` or `NULL` (leave unevaluated). Engine
works at `x -> +oo`; the driver reduces finite / `-oo` / one-sided by
substitution. `gruntz.c` is in the tests' `COMMON_SRC`.

## Engine (Phase 1 — exp-log core)

- `mrv`/`mrv_max` — most-rapidly-varying comparability set, with SubsSet dummy
  bookkeeping for nested mrv elements.
- `compare` — class order via `sign(limitinf(ln a / ln b))`.
- `rewrite_omega` — mrv elements as `A·w^c`, `w -> 0+`.
- `mrv_leadterm` — leading `(c0, e0)` of the series in `w`; §3.3.4 move-up
  (`x -> e^x`, `Log[E^u]->u`) when `x ∈ Ω`.
- `limitinf`, `gz_sign` oracle.
- `series_leadterm` — escalating-order `Series` in `w` with `logx` semantics
  (freeze `Log[w]`, `expand_logs` surfaces inner `-Log[w]`); a `leadterm_lo`
  fallback handles irrational/symbolic `w`-exponents `SeriesData` can't hold
  (e.g. `Log[3]/Log[5]` in `(3^x+5^x)^(1/x)`), but cannot do leading-order
  cancellation.
- Reuses the evaluator: `Series`, `PowerExpand`, `Together`/`Cancel`,
  `zero_test_decide`, `numericalize` — no new low-level machinery.

**Robustness invariants (hold under a 67-case home-made stress battery):**
memoization (limitinf/sign/leadterm) keyed on subexpr; a work budget
(`GRUNTZ_MAX_WORK`) + Series-input size/shape guards ⇒ **never hangs**; a
`gz_result_ok` guard rejects buried divergences (`Sin[ComplexInfinity]`); an mrv
generic-function whitelist ⇒ unsupported heads are left **unevaluated, never
wrong**.

## Phase 2 — essential-singularity isolation (thesis §5.2)

`gruntz_semitractable_limit` / `isolate_semitractable` / `asymptotic_expansion`,
run from the driver when `contains_semitractable(e0)`. For each `F[g]` with
`g -> ±oo`, replace `F[g]` by `Normal[Series[F[dummy],{dummy,Infinity,N}]]` with
`dummy -> g` — the singular part becomes an explicit `Exp`/log head the Phase-1
engine handles. Supported heads (`is_semitractable_head`):
**`Erf`, `Erfc`, `ExpIntegralEi`, `LogGamma`.**

- `Erf`/`Erfc` at `-oo`: reflect (`Erf[g] = -Erf[-g]`, `Erfc[g] = 2 - Erfc[-g]`).
- `Ei` asymptotic valid both directions; `LogGamma` `+oo` only.
- Order escalation `{3,4,5}`, accept once **any two orders concur** (a lone
  order could be a truncation artefact); each order runs `Expand` first so the
  constant parts (e.g. the leading `1` in `Erf -> 1 + …`) cancel across a
  difference before the engine sees them. Two gotchas that cost time: reset
  `g_work=0` per order; "any two agree" (not "two consecutive") because higher
  orders put long truncated series in denominators the engine can't expand.

Supporting Series work landed this session:
`Series[LogGamma[x], {x, Infinity, n}]` = Stirling (DLMF 5.11.1), additive head
`(x-1/2)Log[x] - x + Log[2Pi]/2` + Bernoulli `1/x` tail
(`try_series_loggamma_at_infinity` in `src/calculus/series.c`, wired into the
at-Infinity dispatch). Also a `Normal` array-leak fix in `seriesdata_to_normal`.

## Phase 3 — deep log-tower cancellation (thesis 8.19) — DONE (2026-07-24)

The pure log-tower `(Log[Log[x]+Log[Log[x]]] - Log[Log[x]]) /
Log[Log[x]+Log[Log[Log[x]]]] · Log[x] -> 1` now resolves (`test_log_tower`). It
needs a leading-order cancellation of two nearly-equal logs; three fixes in the
exp-log core made the mrv `Series` path see it:

1. **Factor the w-pole out of each `Log` before expanding** (`expand_logs`).
   `Series[Log[1/w - a], {w,0,n}]` splits the pole into the branch-contaminated
   constants `Log[-a] + Log[-1/a]` (= `2 Pi I`, not `0`, because the frozen
   `lw_sym` has lost its sign). Instead we write `Log[G] = k Log[w] + Log[H]`
   with `H = Expand[G w^{-k}] -> ` a nonzero constant (`k` = the w-valuation from
   a cheap `Series[G]`), so `Series[Log[H]]` is a clean pole-free Taylor series
   and the singular part stays the explicit `k Log[w]` the freeze step turns
   into `k lw_sym`. Only taken when `H(0)` is a clean constant free of `w`/`lw_sym`.
2. **Run the outer `Series` in the positive log-scale `P = -Log[w]`**
   (`series_leadterm`). `Log[w] -> -oo`, so `-lw_sym > 0` and `Log[-lw_sym]`
   (= `Log[x]`) is real — but `Series` mis-branches it to `Log[lw_sym] + I Pi`.
   Substituting `lw_sym -> -P` before `Series` (so `Log[-lw_sym] -> Log[P]`,
   clean) and restoring `P -> -logw` afterwards keeps the whole expansion real.
   This also let the `contains_log_of_sym` guard stop bailing on a `Log[lw_sym]`.
3. **Accept a w-free level as the sub-scale coefficient.** When the expanded
   `f3` is free of `w` the level collapsed to a constant in the sub-scale; it is
   the sub-scale limit, which `limitinf` resolves recursively — no `Series`-in-w
   to corrupt, so accept `(c0 = f3, e0 = 0)` rather than abstaining.

Leak-clean (valgrind == macOS baseline), no new frames; `gruntz_tests` +
`limit`/`limit_assumptions`/`nlimit`/`nseries` green.

### Never-hang guard (complex-contamination reject)

Plain `Limit[8.19]` (Automatic) and `Method->"Series"/"Asymptotic"` now resolve
too. They had **hung**: the cascade recurses through the Gruntz layer on deep
sub-limits (a `x -> 1/x` substitution then nested-log recursion reaches depth
~23), and there `series_leadterm` fed `Series` a **complex-contaminated** input
and looped. Root cause: a 3+-level tower re-introduces `Log[-Log[w]] ->
Log[Log[w]] + I Pi` at a depth the single positive-scale (`P = -Log[w]`)
substitution cannot fix — the mrv engine is real-valued, so a surviving
`Complex[]` is always a branch artefact. Fix (`contains_complex_head`): reject a
frozen-scale `Series` input carrying `I` — `series_leadterm` abstains, and
`expand_logs` leaves such a log unexpanded. The engine then abstains FAST on the
sub-limit; the asymptotic/series cascade layers fail; and the **top-level**
Gruntz layer closes 8.19 to 1. Restores the "never hang, never wrong"
invariant. (A still-deeper 3-level tower like
`Log[Log[Log[x]+Log[Log[x]]]] - Log[Log[Log[x]]]` remains an honest abstention —
bounded ~15s, never wrong — resolving it needs multi-level sign tracking.)

Regression tests: `test_log_tower_no_hang` (plain Automatic 8.19 -> 1;
`Method->"Asymptotic"`/`"Series"` terminate). This also un-hangs the
`test_series_infinity_no_inv_var_leak` case in `tests/test_series.c`, which
evaluates plain `Limit[8.19]` and previously hung the whole `series_tests` binary
(exposing 3 **pre-existing** `D`/`Integrate`-of-`SeriesData` normalisation
failures that were hidden behind that hang — unrelated to this work).

## Generalized stress battery + tutorial (2026-07-25)

- **`tests/test_gruntz_stress.c`** (`gruntz_stress_tests`): 141 hand-verified
  *generalizations* of the thesis examples across 11 families (dominant base,
  cancellation `→-c`, exp-tower ratios `→E^a`, trig-at-vanishing-arg, nested-log
  ratios, Hardy sub-polynomial, conjugate radicals, finite-point power series,
  `Erfc`/`Ei`/`Zeta`/`PolyGamma`/`LogGamma` singularities, `Max`/`Min`) + 4 pinned
  honest abstentions. Each case pins ONLY where an independent hand derivation and
  the engine's output agree; abstentions/wrong answers dropped. Complements (does
  not duplicate) `test_gruntz.c`, which pins the exact thesis expressions.
- **`site/docs/tutorials/computing-limits-gruntz.md`**: advanced tutorial over the
  same families; all 53 transcripts verified exact against the binary.

## Coverage (all verified, all in `tests/test_gruntz.c`)

- **Thesis exp-log (Table 8.1):** 8.1, 8.5–8.9, 8.11–8.13, 8.17, 8.20–8.22.
- **Deep log-tower (8.19):** `test_log_tower` — the building-block
  cancellations (`Log[x+Log[x]]-Log[x] -> 0`, `x(...) -> Infinity`,
  `(...)x/Log[x] -> 1`) plus the full 8.19 ratio `-> 1`.
- **Trig at a vanishing arg:** 8.21/8.22 (ex. 5.1). **Worked:** 3.13, 3.15, 5.4.
- **Home-made stress:** `test_stress_elementary` (20), `test_stress_nested_exp`
  (26), `test_stress_reductions` (finite/one-sided/`-oo`).
- **Phase-2 special:** `test_thesis_special` — 8.23
  `(Erf[x-E^-x]-Erf[x])E^x E^(x^2) -> -2/Sqrt[Pi]`, the `Erfc` spelling, `Erf[x]->1`,
  `Erfc[x]->0`, `Erf[Sqrt[x]]->1`, `x Erfc[x]E^(x^2)->1/Sqrt[Pi]`, `Erf[-oo]->-1`,
  `Ei[x]E^-x x->1`, `Ei[Log[x]]Log[x]/x->1`, `Erfc[x+1/x]/Erfc[x]->E^-2`.
- **Phase-2 LogGamma:** `test_thesis_loggamma` — `LogGamma[x]/(x Log[x])->1`,
  `LogGamma[2x]/(x Log[x])->2`, `LogGamma[x]-(x-1/2)Log[x]+x->Log[2Pi]/2`,
  `(LogGamma[x]-x Log[x]+x)/Log[x]->-1/2`, `x(LogGamma[x+1]-LogGamma[x]-Log[x])->0`.
- **Automatic fallback:** hard exp-log + Phase-2 cases resolve without a Method.
- **Honest abstentions** (`test_honest_abstentions`, pinned unevaluated): 8.31
  Stirling-difference, 8.19 log-tower cancellation, symbolic-sign, bare `Sin[x]@oo`.
- **Phase-2 Gamma** (`test_thesis_gamma`, now RESOLVING): Stirling ratio
  `Gamma[x]/(Sqrt[2Pi] x^(x-1/2) E^-x)->1`, thesis 5.5 `Log[Gamma[Gamma[x]]]/E^x->Infinity`.

## Known gaps (left UNEVALUATED — never wrong, never hang)

- **`Gamma` isolation — RE-ENABLED (FLINT zero test landed, 2026-07-24).**
  `Gamma` is back in `is_semitractable_head`; the `cls==1` branch builds
  `Gamma[g] = Exp[LogGamma[g]]` from LogGamma's additive Stirling series. The
  Stirling ratio and thesis 5.5 now resolve. The thesis-8.31 Stirling-difference
  `Gamma[x+1]/Sqrt[2Pi] - E^-x(x^(x+1/2)+x^(x-1/2)/12)` still abstains — the
  `x^x`-scale tower needs a deeper `Series` cancellation than the machinery
  reaches — but is no longer a hang: its abstention cost fell **173.9s → 14.9s
  (~11.7×)** with the FLINT-backed `is_zero_poly` (the direct confirmation that
  the zero test, not `Series`, was the wall the earlier hang hit). Remaining 8.31
  residual is `Series`-depth, a separate follow-up. See [[project_flint_zero_test]].
- `PolyGamma` at infinity — **DONE (2026-07-25)**. Added
  `try_series_polygamma_at_infinity` (DLMF 5.11.2: digamma `Log[x]` head +
  Laurent tail; `m≥1` pure `x^-m` Laurent) in `series.c`, and 2-arg isolation in
  `gruntz.c` (`is_semitractable_head` recognises `PolyGamma[m,z]`,
  `asymptotic_expansion_pg`, slot-1 argument). Resolves `PolyGamma[x]/Log[x]->1`,
  `x PolyGamma[1,x]->1`, `x(PolyGamma[x]-Log[x])->-1/2`, nested `psi(psi(x))`,
  etc. Symbolic order & `-oo` (poles) stay honest abstentions. See
  `test_polygamma_at_infinity` / `test_series_polygamma_at_infinity`.
- `Zeta` at infinity — **DONE (2026-07-25)**. `try_series_zeta_at_infinity`
  emits the exp-log Dirichlet head `1 + 2^-x + 3^-x + ...` (exponential scale,
  not `1/x`); `Zeta` admitted to `is_semitractable_head` (`+oo` only). Resolves
  `(Zeta[x]-1) 2^x->1`, `Log[Zeta[x]-1]/x->-Log[2]`, `x(Zeta[x]-1)->0`, etc.
  Same-mrv-class ratios `(Zeta[x]-1)/(Zeta[x+1]-1)->2` now resolve too
  (`gz_split_const_exp` normalises `2^-(x+1) -> (1/2)2^-x`). See
  `test_zeta_at_infinity` / `test_series_zeta_at_infinity`.
- Modified Bessel `BesselK`/`BesselI` at infinity — **DONE (2026-07-25)**.
  2-arg isolation (shared `asymptotic_expansion_2arg` with PolyGamma); monotonic
  `Exp[∓z]` envelopes, any/symbolic order. Resolves
  `Exp[x]Sqrt[x]BesselK[0,x]->Sqrt[Pi/2]`, `BesselI[0,x]Exp[-x]Sqrt[x]->1/Sqrt[2Pi]`,
  etc. `BesselJ`/`BesselY` (oscillatory `Cos[x-π/4]`) stay honest abstentions —
  the monotonic mrv engine can't expand bare oscillation. See
  `test_bessel_at_infinity`.
- `Max`/`Min` at infinity — **DONE (2026-07-25)**. `resolve_maxmin` pre-pass in
  the driver rewrites `Max[a,b]` to the eventually-dominant arg via the
  leading-term sign of `a-b` (`gz_sign`, not the limit — so `Max[1/x,2/x]=2/x`),
  recursing through nesting/factors. Also fixed the cascade so
  `layer2_series`/`layer5_lhospital` bail on `Max`/`Min` (were emitting
  `Derivative[Max]` garbage) → Automatic falls through to Gruntz. Resolves
  `x Max[1/x,2/x]->2`, `Min[x,Log[x]]->Infinity`, nested/n-ary. Bounded-osc
  comparisons (`Max[Sin[x],2]`) abstain. See `test_maxmin_at_infinity`.
- `BesselJ`/`BesselY` at infinity — **decay-to-0 DONE (2026-07-25, Automatic)**.
  `magnitude_upper_bound` gives them the `Sqrt[2/(Pi x)]` envelope so the squeeze
  layer resolves `BesselJ[0,x]->0`, `BesselJ[0,x]/x->0`; `x BesselJ[0,x]` stays
  unevaluated (no limit). Exclusive `Method->"Gruntz"` still abstains (oscillatory).
- `Max`/`Min` of bounded oscillation vs a dominating constant/∞ — **DONE
  (2026-07-25, Automatic)**. `layer_maxmin_bounded`: `Max[Sin[x],2]->2`,
  `Max[Sin[x],x]->Infinity`, `Min[Cos[x],-x]->-Infinity`. Ambiguous cases
  (`Max[Cos[x],1/2]`) correctly stay unevaluated.
- `Max`/`Min` of semi-tractable functions — **DONE (2026-07-25)**.
  `resolve_maxmin` isolates the difference before the sign test:
  `Max[PolyGamma[x],Log[x]]-Log[x]->0`, `Max[Gamma[x],x^10]->Infinity`.
- Deep log-tower cancellation (8.19) — **DONE**, both `Method->"Gruntz"` and
  plain `Limit` (Automatic) resolve it to 1 (Phase 3 above).
- `ArcTan[x]@oo` — **DONE (2026-07-25)**. Added `Series[ArcTan[x],
  {x,Infinity,n}]` (`try_series_arctan_at_infinity` in `series.c`) plus a
  kernel-level at-infinity branch `ArcTan[u] = Pi/2 - ArcTan[1/u]` (mirrors the
  `ArcCot` branch) so it composes inside `Plus`/`Times`.
  `Limit[ArcTan[x],x->oo,Method->"Gruntz"] -> Pi/2` and
  `Limit[x(Pi/2-ArcTan[x]),x->oo] -> 1` (both Automatic and Gruntz). See
  `test_arctan_at_infinity` / `test_series_arctan_at_infinity`.
- Genuinely undecidable / no-limit (correct honest abstentions, never wrong):
  bare oscillation (`Sin[x]@oo` has no limit — `Indeterminate` under Automatic,
  abstains under exclusive Gruntz); symbolic-sign `Log[x]/x^s` and symbolic-order
  `PolyGamma[n,x]` (value depends on `Sign[s]` / `n`, need assumptions). NOTE
  (2026-07-25): `Log[x]/x^s` previously returned the **malformed** `Infinity x^-s`
  (an x-dependent "value") — `layer2_series` now discards any leading-term value
  that still mentions the expansion variable, so it abstains cleanly. See
  `test_symbolic_exponent_no_garbage`.
- thesis-8.31 `Gamma` Stirling difference: bounded ~15s abstention (needs a
  deeper `x^x`-tower `Series` cancellation), pinned in `test_honest_abstentions`.

## Gamma re-enabling — DONE (2026-07-24)

Implemented in `src/calculus/gruntz.c`:
1. `"Gamma"` added to `is_semitractable_head`.
2. `isolate_semitractable` `cls==1` branch:
   `Gamma[g] -> simp(mk_exp(asymptotic_expansion("LogGamma", g2, nterms)))`
   (Gamma stays OUT of the `-oo` reflection branch — poles).
Stirling ratio + 5.5 moved to `test_thesis_gamma` (passing); 8.31 kept in
`test_honest_abstentions` (still abstains, now 14.9s not a hang). Residual 8.31
gap is `Series`-machinery depth on the `x^x` tower, a separate follow-up.

## Key files

- `src/calculus/gruntz.{c,h}` — engine.
- `src/calculus/limit.c` — `LIMIT_M_GRUNTZ`, `layer_gruntz`, method parsing.
- `src/calculus/series.c` — `try_series_loggamma_at_infinity` (+ at-Infinity
  dispatch); `seriesdata_to_normal` leak fix.
- `tests/test_gruntz.c`; `tests/test_series.c` (`test_series_loggamma_at_infinity`).
- Docs: `docs/spec/builtins/calculus.md` (Gruntz method), `power-series.md`
  (LogGamma Stirling hook), weekly changelog `docs/spec/changelog/2026-07-20.md`.
- Memory note: `project_gruntz_limit_algorithm`.
- Thesis PDF cached at
  `…/41228d4d-…/tool-results/webfetch-1784834672324-hemqie.pdf`
  (`pdftotext -layout` → §5.2 transforms 5.6–5.16).
