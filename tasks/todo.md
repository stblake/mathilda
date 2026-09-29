# Route every user-facing message through the Quiet/Check funnel

Plan: `/Users/user/.claude/plans/we-need-to-do-misty-boot.md`

## Commit 0 — tooling (no version bump)  ✅ DONE (22bbd3da)
- [x] `tools/check_message_routing.py` — multi-line-aware detection, EXEMPT + BASELINE, ratchet
- [x] Seed BASELINE with the full current backlog (294)
- [x] `make check-messages` target + `.PHONY`
- [x] Wire into `.github/workflows/build.yml`
- [x] Verify: gate green on unchanged tree

## Commit 1 — funnel + highest leverage (bump+tag v0.222)  ✅ DONE
- [x] `mth_message` / `mth_message_gated` / `mth_message_v` / `mth_message_cont` in message.h/.c (guarded printf attr)
- [x] Migrate `common.c:builtin_arg_error` (4 fprintf) → dropped common.c from BASELINE (290 left)
- [x] Build clean, check-messages ratchets, check-c99, behavior verified (Check[Fourier[],CAUGHT]→CAUGHT, Quiet suppresses)

## Commit 2 — confirmed correctness bugs (bump+tag v0.223)  ✅ DONE
- [x] power.c (6), plus.c (2), times.c (1) → mth_message_gated(g_arith_warnings_muted,...)
- [x] linalg/matpow.c (3) → expr_to_string + mth_message
- [x] solve/solvenlsys.c warn_nsdim, solve/solveinv.c emit_ifun (note; ifun-suppress stays moot)
- [x] findmin_common.c fm_warn → mth_message_v(g_fm_quiet,...) (adds note)
- [x] Verified: Check[MatrixPower[..,1/2]]→CAUGHT, Check[Power[0,-2]]→CAUGHT, Quiet suppresses, Quiet[Check]→FAILED. BASELINE 290→279.

## Commit 3 — subsystem helpers → wrappers (bump+tag)
- [ ] fit_warn, fm_warn, DRY re-body fs_msg/dt_msg/ops_msg/root_warn/inv_warn/purefunc

## Commits 4..N — inline sites by module (bump+tag each)
- [ ] solve/
- [ ] linalg/
- [ ] poly/
- [ ] calculus/ (incl. Integrate::nonelem + RischTranscendental routed)
- [ ] numerical_calculus/
- [ ] numerical_roots/
- [ ] strings/ + strings/regex/
- [ ] product/
- [ ] top-level A: int.c, real.c, special_functions/, numbertheory/
- [ ] top-level B: eval/core/context/refine/rootreduce/numberform/funcprog/names/interp/complex_expand/options_builtin/vectors/vectoranal/partitions/bitwise/list/precision/numeric/random/nc_accuracy/expand*

## Final commit — assert-empty + docs (bump+tag)
- [ ] BASELINE == {}; flip gate to assert-empty
- [ ] tests: test_messages.c (or extend test_parallelmixedtower); update deliberate-raw test
- [ ] docs: CLAUDE.md, SPEC.md §9, docs/design/message_routing.md, changelog
- [ ] valgrind spot-check; rebuild code-review graph

## Review
(to be filled in)

## `Integrate` — `ParallelMixedTower`: the real form of the logarithmic part (logrewrite.m, 2026-09-28)

Port of the research campaign of 2026-09-28 (logrewrite.{py,wl,mac}: Rioboo's collapse of conjugate
pairs of logarithms to real logarithms, arctangents and hyperbolic arctangents) to the Mathilda
package, with pre/post measurements. `.m` files, one C test and docs only; no C change.

### Plan
- [x] 1. Probe every kernel construct of logrewrite.wl in Mathilda first (three probe scripts): FreeQ on
  complex atoms, ComplexExpand on nested radicals, CountRoots, $InputFileName, nested TimeConstrained,
  Do-iterator capture, PolynomialExtendedGCD / PolynomialRemainder over algebraic constants, D and N of
  ArcTan/ArcTanh forms for the verify gate
- [x] 2. Baseline (pre) on the unchanged package: Charlwood 50 (charlwood_wl.py mathilda), review corpus
  (review_wl.py, SYSTEM=mathilda), tests/build/parallelmixedtower_tests
- [x] 3. `src/internal/mixed/logrewrite.m`: one function per function of logrewrite.wl, loaded with
  `LoadModule["mixed/logrewrite.m"]` inside the private context; host adaptations (`_Complex`,
  `RectPow`, `SturmCount`, `RootRadicals`, `CanRaw`, lr-prefixed iterators); default `LogToReal` in
  ParallelMixed.m when the file is missing; hook at the end of iPIM (`lrTerms` -> `LogToReal`)
- [x] 4. Develop against a copy of src/internal selected with MATHILDA_HOME; compare 42 integrands with
  Mathematica (trace lines and derivative checks at points including x = 2, 3, 5/2)
- [x] 5. Fix found in every port (a defect of the campaign, not of the port): the arctangent argument over
  the radical is a polynomial in y (YQuot / _y_quot / pm_y_quot) -- a radical in a denominator was
  read on the principal branch by the surface layer (Charlwood P5 wrong for Tan[x] < 0 in Mathematica
  and Mathilda); re-verified in the four ports
- [x] 6. Install, C test `test_method_real_form`, MATHILDA_DIVERGENCES.md A19-A24 / B7 / D, changelog
- [x] 7. Post measurements on the installed package (Charlwood, review corpus, C tests) and the review
  section below

### Review
- Both measurements on build 0.222 on a quiet machine, the baseline through `MATHILDA_HOME` pointing at
  a copy of the previous module tree (the binary was rebuilt from 0.221 to 0.222 during the session):
  Charlwood 48/50 -> 49/50 (A19 passes the gate in its real form), kernel 20.4 s -> 23.6 s, median
  0.209 -> 0.204 s, results with I 16 -> 0; review 303/68/0 -> 303/68/0, 693 s -> 602 s, 44 -> 0
  returned integrals with I; `parallelmixedtower_tests` 11 tests pass incl. `test_method_real_form`.
- The port exposed a defect of the campaign itself (the arctangent argument with the radical in a
  denominator, P5 wrong for Tan[x] < 0 in Mathematica too), fixed in all four ports and re-verified.
- Six new host divergences (A19-A24) and one behavioural one (B7), each with a one-line repro and a
  workaround in `logrewrite.m`; A24 (`Can`'s field detour returning `Dot[{}, Inverse[{}], {}]`) has
  no standalone repro yet -- it appeared on A40 after P4 in one kernel and is worth a core look.
- Not done: the `E^x` / `E^(2 x)` independent-generator defect of BuildTower (both ports certify
  `1/(1 + E^x + E^(2 x))` as not elementary) was noticed and left alone; the rewrite's 30 s budget is
  not interruptible by the package budget (B7). Nothing committed; the user was committing their own
  work in this repository at the time.


## `Integrate` — `ParallelMixedTower`: Charlwood 50 at least as fast as the faster of Maxima and Mathematica on every row (plan, 2026-09-28)

Goal: on each of Charlwood's fifty integrals, Mathilda's kernel time under the protocol of
`charlwood_wl.py mathilda` (one process per integral, lazy load, two warm-ups, `Integrate[f, x,
Method -> "ParallelMixedTower"]`) is at most min(Maxima, Mathematica) measured the same way with
the same package state (real-form rewrite active in all three).  Nothing here changes what the
package computes: every item is a cost fix, verified by the same corpus runs as before.

### Where Mathilda stands (build 0.228, HEAD f74307de, 2026-09-28)

- 49/50 verified, **26.0 s** of kernel time (A39 the same honest `failed` as every system).
  Stored: `mixed/charlwood_mathilda.json` (0.190 moved to `mixed/stress/charlwood_mathilda_0190.json`,
  the fresh run also in `mixed/stress/charlwood_mathilda_0228.json`).
- Like-for-like targets re-measured today WITH the rewrite (the research-directory result files
  `charlwood_wl.json` of Sep 27 and `maxima/charlwood_maxima.json` of Sep 23 predate
  logrewrite): Maxima **9.34 s** (`mixed/stress/charlwood_maxima_20260928.json`; A1 3.46 -> 0.15 s,
  A13 0.43 -> 0.09 since the review changes), Mathematica **17.18 s**
  (`mixed/stress/charlwood_wl_20260928.json`).  Sum of the per-row bests **9.20 s**; Maxima is the
  faster of the two on 45 rows.
- Mathilda is slower than the best on **50/50** rows (only 2 within 1.25x: P4, P8); median ratio 2.4;
  worst A32 43x, A28 13.9x, A27 13.6x, A40 7.0x.
- Run-to-run noise is ~1-2 % (five passes of eight integrals in `mixed/stress/stage_timing`); the
  differences that matter are systematic.

### Where the time goes

Stage timing from an instrumented copy of the package (`mixed/stress/stage_timing/instrument.py`,
`run_stages.py`, result `stages_0228.json`; PMStage marks summed per stage over one timed call):

| stage | ms over the suite | share | largest rows (ms) |
|---|---:|---:|---|
| residue analysis of the prime loop (RealisePoints, RealiseClass, residues) | 4845 | 20.6% | A28 3106, A27 1110, A29 227, A34 40, A19 35 |
| ansatz assembly: parts + entries + dense fill (AnsatzSystem) | 3373 | 14.3% | A3 687, A2 632, P4 614, A16 250, A39 138 |
| real form of the logarithmic part (logrewrite.m) | 3151 | 13.4% | P8 1237, A40 973, P4 374, A3 103, A2 101 |
| unit / S'-unit logand columns, recomputed on every rung | 2507 | 10.7% | P4 427, P5 214, A3 199, A2 192, A35 179 |
| degree bounds (PlaceBounds), per rung | 2300 | 9.8% | A3 145, A2 141, P4 132, A16 104, P5 102 |
| numeric verify-or-decline gate (33 points, 30 digits) | 1625 | 6.9% | P3 139, A40 125, A12 87, A13 76, A27 75 |
| residual f - sum tau D u/u (Can with algebraic constants) | 1364 | 5.8% | A32 1012, A28 160, A27 131, A29 42, A34 13 |
| hypertangent residue + tower specials (FactorList of the derivations) | 853 | 3.6% | A39 72, A3 51, A2 49, P5 41, P4 39 |
| entry of each rung (TransConsts, AlgPolesQ, TPad, Can /@ f) | 754 | 3.2% | A2 58, A3 51, P4 50, A16 47, A39 36 |
| S'-units over the specials (norm search, Miller functions) | 742 | 3.2% | P8 181, P5 77, A37 72, A39 34, A19 34 |
| number field of the ansatz (ToNumberField + basis), per rung | 623 | 2.6% | P4 196, A40 120, A35 76, A3 47, A2 46 |
| tower construction | 413 | 1.8% | A39 20, A35 18, P3 17, A11 16, A6 15 |
| post-solve residual check and back-mapping | 390 | 1.7% | A35 133, A40 51, P4 49, A3 22, A2 21 |
| exact row reduction | 226 | 1.0% | P4 52, A3 48, A2 46, A16 13 |
| **all stages (instrumented run, rewrite active)** | **23517** | | |

Per-row times in ms (Mathilda 0.228; after P1; the two targets re-measured today; the target is the smaller):

| id | integrand | Mathilda 0.228 | after fix P1 | Maxima | Mathematica | target | ratio |
|---|---|---:|---:|---:|---:|---:|---:|
| A28 | `(x**2 + 1)/((1 - x**2)*sqrt(x**4 …` | 4321 | 3147 | 227 | 576 | 227 (Max) | 13.9 |
| A40 | `atan(x*sqrt(1 - x**2))` | 1736 | 1583 | 225 | 2314 | 225 (Max) | 7.0 |
| A27 | `sqrt(sin(x))/(sin(x)**2 + 1)` | 1525 | 1402 | 103 | 526 | 103 (Max) | 13.6 |
| A2 | `atan(x + sqrt(1 - x**2))` | 1462 | 1408 | 704 | 1148 | 704 (Max) | 2.0 |
| A3 | `x*atan(x + sqrt(1 - x**2))/sqrt(1…` | 1501 | 1392 | 699 | 1153 | 699 (Max) | 2.0 |
| P5 | `cos(x)**2/sqrt(cos(x)**4 + cos(x)…` | 604 | 576 | 133 | 389 | 133 (Max) | 4.3 |
| A35 | `sqrt(-sqrt(sec(x) - 1) + sqrt(sec…` | 663 | 738 | 384 | 309 | 309 (Mat) | 2.4 |
| A16 | `x*log(x + sqrt(1 - x**2))/sqrt(1 …` | 643 | 665 | 345 | 471 | 345 (Max) | 1.9 |
| A19 | `atan(x)/(x**2*sqrt(1 - x**2))` | 480 | 441 | 136 | 496 | 136 (Max) | 3.2 |
| P8 | `sqrt(tan(x)**2 + 2*tan(x) + 2)` | 1748 | 1618 | 1356 | 1560 | 1356 (Max) | 1.2 |
| A39 | `asin(x*sqrt(1 - x**2))` | 556 | 511 | 249 | 921 | 249 (Max) | 2.1 |
| P3 | `-asin(sqrt(x) - sqrt(x + 1))` | 389 | 341 | 102 | 180 | 102 (Max) | 3.3 |
| A12 | `x**3*asec(x)/sqrt(x**4 - 1)` | 398 | 379 | 147 | 271 | 147 (Max) | 2.6 |
| A1 | `x*log(x + sqrt(x**2 + 1))*log(x**…` | 484 | 379 | 153 | 193 | 153 (Max) | 2.5 |
| A11 | `x**3*asin(x)/sqrt(1 - x**4)` | 397 | 372 | 149 | 243 | 149 (Max) | 2.5 |
| P4 | `log(x*sqrt(x**2 + 1) + 1)` | 2174 | 1922 | 1734 | 2270 | 1734 (Max) | 1.1 |
| A13 | `x*log(x + sqrt(x**2 + 1))*atan(x)…` | 264 | 275 | 89 | 130 | 89 (Max) | 3.1 |
| A5 | `log(x + sqrt(x**2 + 1))/(1 - x**2…` | 352 | 319 | 134 | 356 | 134 (Max) | 2.4 |
| A37 | `atan(x*sqrt(x**2 + 1))` | 351 | 401 | 216 | 360 | 216 (Max) | 1.9 |
| A6 | `asin(x)/(x**2 + 1)**(3/2)` | 389 | 340 | 156 | 266 | 156 (Max) | 2.2 |
| A22 | `asin(x)/(x**2*sqrt(1 - x**2))` | 204 | 198 | 36 | 130 | 36 (Max) | 5.5 |
| A29 | `(1 - x**2)/((x**2 + 1)*sqrt(x**4 …` | 367 | 353 | 202 | 313 | 202 (Max) | 1.7 |
| A31 | `sqrt(sin(x) + 1)*log(sin(x))` | 246 | 231 | 96 | 106 | 96 (Max) | 2.4 |
| A7 | `log(x + sqrt(x**2 - 1))/(x**2 + 1…` | 235 | 237 | 116 | 201 | 116 (Max) | 2.0 |
| P10 | `x**3*exp(asin(x))/sqrt(1 - x**2)` | 300 | 207 | 90 | 94 | 90 (Max) | 2.3 |
| A4 | `asin(x)/(sqrt(1 - x**2) + 1)` | 155 | 158 | 44 | 89 | 44 (Max) | 3.6 |
| A34 | `sin(x)/sqrt(1 - sin(x)**6)` | 154 | 158 | 48 | 159 | 48 (Max) | 3.3 |
| A21 | `atan(x)/(x**2*sqrt(x**2 + 1))` | 155 | 135 | 34 | 96 | 34 (Max) | 4.0 |
| A24 | `log(x)/(x**2*sqrt(x**2 + 1))` | 137 | 125 | 29 | 83 | 29 (Max) | 4.3 |
| P1 | `log(x)*asin(x)` | 198 | 193 | 134 | 101 | 101 (Mat) | 1.9 |
| A38 | `-atan(sqrt(x) - sqrt(x + 1))` | 143 | 167 | 79 | 85 | 79 (Max) | 2.1 |
| A30 | `log(sin(x))/(sin(x) + 1)` | 124 | 127 | 44 | 69 | 44 (Max) | 2.9 |
| A26 | `x*log(x)/sqrt(x**2 + 1)` | 113 | 112 | 30 | 81 | 30 (Max) | 3.7 |
| A36 | `x*log(x**2 + 1)*atan(x)**2` | 109 | 127 | 59 | 46 | 46 (Mat) | 2.8 |
| A15 | `x*log(x + sqrt(x**2 + 1))/sqrt(x*…` | 114 | 108 | 28 | 71 | 28 (Max) | 3.9 |
| A25 | `x*asec(x)/sqrt(x**2 - 1)` | 155 | 135 | 59 | 104 | 59 (Max) | 2.3 |
| A20 | `x*atan(x)/sqrt(1 - x**2)` | 152 | 163 | 90 | 214 | 90 (Max) | 1.8 |
| P6 | `sqrt(tan(x)**4 + 1)*tan(x)` | 162 | 131 | 62 | 167 | 62 (Max) | 2.1 |
| A33 | `tan(x)/sqrt(tan(x)**4 + 1)` | 165 | 158 | 93 | 138 | 93 (Max) | 1.7 |
| A18 | `x*atan(x)/sqrt(x**2 + 1)` | 129 | 117 | 55 | 86 | 55 (Max) | 2.1 |
| A32 | `sec(x)/sqrt(sec(x)**4 - 1)` | 1252 | 99 | 42 | 82 | 42 (Max) | 2.4 |
| A8 | `log(x)/(x**2*sqrt(x**2 - 1))` | 111 | 110 | 57 | 88 | 57 (Max) | 1.9 |
| A14 | `x*log(sqrt(1 - x**2) + 1)/sqrt(1 …` | 119 | 119 | 66 | 89 | 66 (Max) | 1.8 |
| A23 | `x*log(x)/sqrt(x**2 - 1)` | 125 | 108 | 60 | 87 | 60 (Max) | 1.8 |
| P2 | `x*asin(x)/sqrt(1 - x**2)` | 77 | 75 | 29 | 42 | 29 (Max) | 2.6 |
| A10 | `x*log(x + sqrt(x**2 - 1))/sqrt(x*…` | 83 | 77 | 31 | 39 | 31 (Max) | 2.5 |
| P9 | `sin(x)*atan(sqrt(sec(x) - 1))` | 108 | 88 | 46 | 53 | 46 (Max) | 1.9 |
| P7 | `tan(x)/sqrt(sec(x)**3 + 1)` | 54 | 49 | 23 | 36 | 23 (Max) | 2.1 |
| A17 | `log(x)/(x**2*sqrt(1 - x**2))` | 123 | 112 | 94 | 86 | 86 (Mat) | 1.3 |
| A9 | `sqrt(x**3 + 1)/x` | 28 | 29 | 19 | 9 | 9 (Mat) | 3.2 |
| **total** | | **26034** | **22415** | **9336** | **17176** | **9197** | |

C-level sampling (`sample`, 1 ms, three calls per integral, `mixed/stress/profiles_20260928/`):
A28 spends 68 % under `builtin_plus` -> `expr_compare` -> `collect_symbols_in` (sort.c:271);
P8 35 % in GMP integer gcds under `PolynomialGCD`/`Cancel`; the retry-heavy rows (A2, A3, P4, P5,
A16, A19) 10-20 % in `rb_collect_symbol_names` (match.c:1544) + its `strcmp`, called from
`replace_bindings` on every down-value application; 10-16 % everywhere in malloc/free churn.

### Root causes (each verified in isolation)

- **F1 `Cancel[e, Extension -> Automatic]` on AlgebraicNumber input.**  rat.c
  `builtin_cancel_compute` runs `cancel_auto_gcd_quotient` (PolynomialGCD + PolynomialQuotient,
  each re-detecting the splitting field) before the native `flint_field_cancel`: 315 ms on a
  LeafCount-73 fraction over Q(Sqrt 2) versus 0.2 ms for the 1-argument `Cancel[e]`, which reduces
  identically (checked on a fraction with a common factor over Q(Sqrt 2, I)).  `Can`'s field
  detour hits this on every call (A32: three residual `Can`s = 0.7 s; the detour is a net loss on
  small fields and a net win on P8/P4 -- the suite is 29.0 s with `$CanFieldEnabled = False`);
  `CanL`/`KGcd` of logrewrite.m hit it on nested radicals.  **Verified fix**: dropping the option
  in `Can`'s detour gives 26.0 -> 22.4 s, 49/50 unchanged (A32 1.25 -> 0.10, A28 -1.2, P4 -0.25,
  P8 -0.13, A40 -0.15, A27 -0.12, A3 -0.11; `mixed/stress/charlwood_mathilda_0228_canfix.json`).
- **F2 `VanishOrder`** (`SeriesCoefficient[ub, {e, 0, k}]`, k < 8, on Sqrt[q(rho + e)] with Sqrt[2]
  constants): 3.33 s of A28's 3.4 s, 0.13 s of A29 -- the kernel's Series is the Plus-ordering hot
  spot of the A28 profile.  `DeepResidues` (line 491), `RealiseClass` (1190), `FundamentalUnit`
  (1273), `ResidueFree` (1512/1516), `SpecialData` (1617) and `InfinityData` (1841-1847) still use
  Series too.
- **F3 `RealiseClass`**: `NullSpace[rows, ZeroTest -> (RootReduce[Together[#]] === 0 &)]` over
  Gaussian rationals: 1.08 s of A27's 1.4 s.
- **F4 the rewrite** (`LogToReal`, 3.15 s): P8 1.24 s and A40 0.97 s sit in the conjugate-pair
  loop (`LogToAtan` -> `RiobooAtan`/`RealZeroFreeQ`/`SturmCount`/`YQuot`, all through `CanL`) and
  in the final `CanL[Expand[arg]]` pass -- F1 on nested radicals and Root objects.
- **F5 the retry ladder recomputes rung-invariant work** (P4, A2, A3, A16 run 5 rungs; P5, A19,
  A22, A4 3): logand columns 2.5 s, PlaceBounds 2.3 s, FieldData 0.6 s, rung entry 0.75 s.
- **F6 the assembly** (3.4 s): `AppendTo` in the entries/rhs/parts loops is O(n) per call in
  Mathilda (5,000 appends 0.26 s, 10,000 1.5 s, 20,000 6.0 s; `Sow`/`Reap` 5,000 in 4 ms), the
  dense fill by `aug[[i, j]] = v` copies a row per assignment (150 x 70 in 31 ms, 300 x 140 in
  470 ms), and `CoefficientRules` on the unexpanded product runs once per column per coordinate.
- **F7 the verify gate** (1.6 s, up to 140 ms on P3): 33 points, the integrand evaluated twice
  per point and D[surf] once, all at 30 digits.
- **F8 kernel overheads** (the floor of the small rows: A9 28 ms against 7, P7 50 against 21,
  A15 105 against 29): (a) `expr_compare`'s polynomial-degree tier outside an Orderless sort
  walks both subtrees on every comparison (`g_symmemo == NULL` path of sort.c:654 does not consult
  the persistent `symset_cache`; Plus's already-sorted pre-check at plus.c:877 and every
  `Times`/`Plus` re-evaluation pay it); (b) `replace_bindings` collects the symbol names of every
  binding value with an O(k^2) `strcmp` dedup on every rule application (match.c:1540-1580, added
  with the A11 capture fix in v0.212 -- the small rows doubled between 0.190 and 0.228: A15
  37 -> 105 ms, P2 46 -> 77, A9 17 -> 28); (c) `AppendTo` copies the whole list (core.c:2140);
  (d) `Series` of a square root of a polynomial with algebraic constants (F2).

### Plan

Phase A -- the package (`src/internal/mixed/ParallelMixed.m`, `logrewrite.m`).  After every item:
`python3.11 charlwood_wl.py mathilda` from the research directory (49/50, 0 wrong, time per row),
`SYSTEM=mathilda python3.11 review_wl.py` (304/67/0, no regression), `tests/build/parallelmixedtower_tests`.
Expected effect in parentheses (suite total after the item, from today's stage data).

- [ ] **P1 `Can` field detour without `Extension -> Automatic`** (verified today: 26.0 -> 22.4 s).
      Also guard `FieldData` against the A24 garbage (`piv === {}` or a singular `B` -> `$Failed`,
      so `Can` falls back to `CanRaw`), and reuse `$CanFieldMemo` for `AnsatzSystem`'s `FieldData`.
- [ ] **P2 `VanishOrder` by exact arithmetic** (A28 -2.9 s, A29 -0.13 -> ~19.4 s): ord_P(a + b y)
      at an unramified constant place is the multiplicity of rho in N(u) = a^2 - q b^2 when the
      conjugate a - b y does not vanish at P; when both vanish, divide (g - rho) out of a and b and
      recurse.  Same-answer check against the Series version on the residue corpus.  Then audit
      the remaining `Series` uses (lines 491, 1190, 1273, 1512, 1516, 1617, 1841-1847): the ones on
      Sqrt[polynomial] with algebraic constants go through the package's own exact truncated
      arithmetic (`LaurentPolyTimes` / `PMTruncate`, the A13 route), the rest stay.
- [ ] **P3 `RealiseClass` / `FindElement` null spaces over one field** (A27 -1.0 s -> ~18.4 s):
      map the rows with `FieldData` and take `NullSpace` without `ZeroTest` (AlgebraicNumber
      arithmetic is canonical, zero is syntactic; plain `NullSpace` for Gaussian rationals), as
      `AnsatzSystem` already does; the `RR`/`RRad` calls on the rows go with it.
- [ ] **P4 the rewrite through the fast canonical form** (P8 -1.0, A40 -0.8, P4 -0.3 -> ~16.3 s):
      `CanL` = the memoised field detour of P1 (with the guard) where the constants map into one
      number field, `CanRaw` otherwise; `KGcd` = `PolynomialGCD` on the AlgebraicNumber form
      without the option (probe first that the 2-argument form reduces over the field); memoise
      `RootRadicals` per Root object and `RealZeroFreeQ` per polynomial (it is asked for A and B of
      every pair).  Re-measure `lr-Ipairs` and `lr-finalCanL` with the stage tool.
- [ ] **P5 rung-invariant work cached in `$analyses[key]`** (-3.0 s -> ~13.3 s): the logand columns
      of `AnsatzSystem` per (unit, coordinate) (they do not depend on the bounds), `FieldData` by
      atom set (the existing memo), the `VP`/`VInf` data of `PlaceBounds` (SpecialData/InfData are
      memoised; the loop around them is not), `TPad`/`Can /@ f`/`TransConsts`/`AlgPolesQ` at rung
      entry.
- [ ] **P6 the assembly without quadratic primitives** (3.4 -> ~1.2 s, -2.2 s -> ~11.1 s):
      `Sow`/`Reap` (or `Table`) for `parts`, `entries`, `rhs`; the matrix built once from the
      rules (a `Table` over the row/column association, or `Normal[SparseArray[...]]` once its
      cost is probed) instead of 10^4 `Part` assignments; `CoefficientRules` once per column.
- [ ] **P7 the verify gate** (-0.5 s with the integrand evaluated once per point; -1.2 s if the
      grid is cut to 8-12 points with early exit -- the grid was set to 33 deliberately after an
      off-grid wrong surface, so the point count is the user's call; the exact post-solve residual
      check of `AnsatzSystem` already rules out an inconsistent solve, the gate only pins the branch).
- [ ] **P8 the floor of the small rows** (~-1.5 s over the 30 rows under 150 ms): `PlaceBounds`
      costs 20-45 ms on towers with no curve (`InfData`'s `VInf` per generator), the tower-specials
      pass 8-12 ms (`FactorList` of every derivation, twice per rung with the ConicToLine retry),
      rung entry 8 ms, the gate 5-25 ms; each has a trivial-case shortcut.
- [ ] Checkpoint: four-way table (`mixed/charlwood_four.py`), per-row ratio; expected ~10-11 s
      against the 9.2 s sum of bests, with the sub-100 ms rows still 1.5-3x over their targets.

Phase B -- the kernel (each with a differential toggle in the style of `MATHILDA_NO_SYMSET_CACHE`,
the message/valgrind gates, and a bump + tag).

- [ ] **C1 `Cancel`/`PolynomialGCD` with `Extension -> Automatic` on AlgebraicNumber-coefficient
      input** take the native field path (`flint_field_cancel` / the FLINT number-field gcd)
      before `cancel_auto_gcd_quotient`; memoise `extension_autodetect` by the set of algebraic
      atoms (nested radicals re-derive the compositum on every call).  Makes P1/P4 host-independent.
- [ ] **C2 `expr_compare` tier 3 outside a sort** consults `expr_symset_cache_get`/`put` in the
      `g_symmemo == NULL` path, and the already-sorted pre-checks of `Plus` (plus.c:877) and
      `Times` run inside a `symmemo_begin`/`symmemo_free` scope; consider caching the degree vector
      per node the same way (`expr_poly_degree` is the next walk).  Expected 10-20 % on A27, A32,
      A40, P4; the A28 pathology disappears with P2 regardless.
- [ ] **C3 `replace_bindings`' danger set**: names are interned, so the `strcmp` fallback in
      `rb_collect_symbol_names` goes (pointer compare only); compute the set lazily, only when
      `rb_rec` meets a scoping construct, from the per-node `symset_cache` of the binding values
      rather than a fresh walk; or precompute per definition the locals of its scoping constructs
      and test the values for those names only.  Expected 10-20 % on the retry-heavy rows and most
      of the 0.190 -> 0.228 regression of the small rows.
- [ ] **C4 `AppendTo` amortised O(1)**: grow in place with spare capacity when the value is
      uniquely owned (refcount 1); `Sow`/`Reap` in P6 is the package-side fix that does not wait
      for this.
- [ ] **C5 `NullSpace`/`RowReduce` with a `ZeroTest`** on exact algebraic entries: a native exact
      path over Q(i) and AlgebraicNumber (after P3 the package no longer needs it).
- [ ] **C6 `Series` of Sqrt[polynomial] with algebraic constants**: 400 ms per coefficient today
      (after P2 the package no longer needs it; a kernel-level look at why Plus ordering dominates).
- [ ] Re-measure the fifty and the review corpus, four-way table, `MATHILDA_DIVERGENCES.md`
      (new section on the cost divergences B8-B12: F1, F2, F3, F6, F8b), changelog, lessons.

### Measurement protocol and pitfalls found today

- Measure through the lazy load only.  A `Get` of the package after the lazy load re-defines the
  trivial default `LogToReal` and `LoadModule["mixed/logrewrite.m"]` is then a no-op: the
  re-read package runs WITHOUT the real-form rewrite (P8 0.40 s instead of 1.65 s, suite 19.7 s
  instead of 26.0 s).  `review_wl.py`'s `MATHILDA_PKG` mode has this property; `MATHILDA_HOME`
  pointing at a module tree does not.  The stage tool Gets its own instrumented `logrewrite.m` for
  this reason.
- The Maxima and Mathematica columns must be re-measured whenever the packages change; today's
  files carry the date in their names.  Table 1 of the Maxima paper is NOT touched by this work.
- `$`-prefixed package globals (`$CanFieldEnabled`, `$CanFieldMemo`, `$analyses`) live outside
  the private context (bare names reach them); `$CanFieldEnabled = False` is the A/B switch of P1.
- Mathilda has no evaluated-flag: a declined `Integrate` re-runs on every reference (the runner's
  `/. Integrate -> PMDeclined`).

### Review  (2026-09-28, v0.229)

**Result.** Charlwood's fifty, one process per integral under the protocol of `charlwood_wl.py
mathilda`: kernel time **23.53 s -> 13.10 s** (1.80x), 49/50 verified throughout, 0 wrong. Targets
re-measured the same day: Maxima 9.34 s, Mathematica 17.18 s, the sum of the per-row bests 9.20 s.
Mathilda is now faster than Mathematica over the suite as a whole, and on four rows it is faster
than BOTH: **P8 0.64 s against 1.36 / 1.56, P4 1.24 against 1.73 / 2.27, A33 0.08 against 0.09 /
0.14, A17 0.07 against 0.09 / 0.09**; A28 (1.02x) and A38 (1.04x) sit on the target. Biggest
single-row wins: A28 3.45 s -> 0.23 (14.9x), P8 1.70 -> 0.64, P4 2.01 -> 1.24, A40 1.48 -> 1.13,
A27 1.43 -> 0.96. Worst remaining ratios: A27 9.3x and A40 5.0x of Maxima. Files:
`mixed/stress/session_base_0228.json` (before) and `charlwood_mathilda_0229.json` (after);
`mixed/stress/vs_targets.py RUN.json [PREV.json]` prints the per-row table against the targets.

**The review corpus got BETTER, not just faster**: 304 correct / 67 gap -> **329 correct / 42 gap,
0 wrong answers**, 582 s -> 249 s of kernel time. That is the `Series` fix below, not the cost work.
(The one row the runner marks WRONG, R272, previously TIMED OUT at 120 s and now returns an honest
`needs torsion realisation` decline at 113 s; `KNOWN_GAPS['mathilda']` records it as `timeout` and
the tolerance rule does not cover that reason. No antiderivative is returned, so no wrong answer.)

**The largest single find was a correctness bug, not a cost one.** `Series`/`SeriesCoefficient` are
`HoldAll` and did not resolve a series variable carrying a symbol-valued OwnValue, which Mathematica
does. A `.m` routine that expands in a `Unique[...]` symbol therefore got its INPUT back, silently:
`RealiseClass` built the rows of its exact linear algebra out of a non-series, so every unbalanced
configuration at infinity was solved against a matrix carrying `Sqrt[1 + w^4]` where a rational
belonged, and realised nothing. 25 more integrals of the review corpus now close.
(MATHILDA_DIVERGENCES.md A25.)

**What landed** (each measured against the suite and the corpus; the plan's numbering):

| item | effect | note |
|---|---|---|
| P1 `Can` field detour without the option + `FieldData` span guard + shared memo | -1.3 s | A32 1.15 s -> 0.08 |
| P2 `VanishOrder` by the norm's multiplicity | -3.0 s | A28 3.24 -> 0.22, at Maxima's 0.23 |
| Series held-variable fix (kernel) | +1.4 s | a correctness fix that does MORE work; pays for itself in the corpus |
| P5 rung-invariant logand columns + P6 assembly without quadratic primitives | -5.4 s | with C2/C3 |
| C2 `expr_compare` degree tier uses the persistent symset cache | (in the -5.4) | |
| C3 `replace_bindings` danger set built once, lazily | (in the -5.4) | most of the v0.212 small-row regression |
| P3b `RealiseClass` cheapest configuration first | -0.2 s | A27 1.19 -> 1.05 |
| P4 `KGcd` over the number field + `RootRadicals`/`RealZeroFreeQ` memos | -0.2 s | |
| residual `f - Sum tau Du/u`: one `Can` per logand, not two | -0.2 s | A27 1.05 -> 0.95 |
| **`CanRaw` asks for `Extension -> Automatic` only when a Gaussian constant is present** | **-1.4 s** | P8 1.62 -> 0.59, A40 1.40 -> 1.11 |
| P6b the 1-part and y-part columns share the derivative contraction | small | strict work removal |
| P8a `Can` once in the tower-specials scan, `ClassifyPrime` memoised | small | strict work removal |

**What was tried and REJECTED** (each measured, each left out):

- **P3 as planned** -- `NullSpace` over one number field instead of `ZeroTest -> RootReduce`. A27
  1.38 s -> 2.99 s: the entries are large radical expressions and mapping them into the field costs
  more than the zero tests. A cheap numeric screen in the `ZeroTest` was also neutral. The honest
  answer for that null space is a native exact path (C5), not a package rewrite.
- **P4's `CanL` through `Can`'s field detour**. The suite 14.9 s -> 23.5 s (A40 1.44 -> 8.07): the
  rewrite's constants change from pair to pair, so nearly every call builds a fresh number field.
- **`$CanFieldEnabled = False`** (now that `CanRaw` is cheap) and **a larger size gate on the
  detour**. A27 likes it (0.98 -> 0.61) and P8 does not at all (0.68 -> 3.39): P8's compositum
  fractions are exactly what the detour is for. Left at the original gate.
- **Pre-contracting `FieldData`'s `back` to one dot product.** `(v . Binv) . basisRad` multiplies
  each radical monomial by ONE rational; contracting first makes every coordinate carry a sum of
  radicals. Suite +1.7 s.
- **The verify gate at 30-digit sample points instead of exact rationals.** Isolated it looked like
  85 ms -> 31-65 ms for the 33 evaluations; over the suite it was +0.95 s.
- **P7's single integrand evaluation** is kept (it is strictly less work) but is worth ~0 ms: the
  gate's 12 % is `N[D[surf] /. x -> pt, 30]`, not the integrand.

**Kernel bugs found on the way** (MATHILDA_DIVERGENCES.md A25-A27):

- A25 `Series` held-variable resolution -- FIXED.
- A26 multivariate `PolynomialGCD` with `AlgebraicNumber` coefficients returns a non-divisor.
  Half-fixed: `collect_variables` no longer enrols an `AlgebraicNumber` as a polynomial variable
  (`Variables[a + x]` is `{x}` now) and `Extension -> Automatic` no longer answers `1` on the
  univariate case. The classical content computation over `K` still returns the second operand when
  the main-variable degrees are equal. Needs a native `flint_field_gcd`; every component exists.
- A27 `SparseArray` is not implemented, and `Normal` of the unevaluated head returns the
  `SparseArray[...]` expression -- a silent non-matrix. The first attempt at P6 assembled the
  augmented matrix this way and every integral failed; the fill is now row by row (158 ms -> 15 ms
  on 300 x 140, which is also faster than `SparseArray` would have been).

**Where the remaining 3.9 s sits** (stage timing, `mixed/stress/stage_timing`, after the changes):
the verify gate 12 %, the ansatz `parts` + `entries` 22 %, the logand columns 10 %, `PlaceBounds`
9 %, the conjugate-pair rewrite 6 %. Per row the gap is A40 (-0.88 s) and A27 (-0.84) followed by
A3/A2 (-0.35 each) and thirty rows each 20-80 ms over. The C-level profile of both A40 and A27 is
now generic evaluator overhead -- 15-17 % malloc/free, 5 % symbol interning, 9 % ordering -- not a
package hot spot, so the next real step is C4 (`AppendTo`/args-array churn), C5 (a native exact
null space over `Q(i)`/`AlgebraicNumber`) and the `flint_field_gcd` of A26, not more `.m` work.

**Not done from the plan**: C1 (only its two correctness halves), C4, C5, C6, and the four-way
table checkpoint.

**The full 510-binary C test suite: 8 failures, ALL verified PRE-EXISTING** (the five changed
`src/*.c` files and `src/internal` reverted to `HEAD`, rebuilt, identical failure in every case).
Zero regressions. Two new regression tests were added for the fixes: the held-variable `Series`
case (`tests/test_series.c`) and `AlgebraicNumber` as a constant rather than a polynomial variable
(`tests/test_algebraicnumber.c`).

- `crc_corpus_tests`: 6 diff-nonzero against a baseline of 3. The three new ones are the
  `Sqrt[(a + b x)/(c + d x)]` family through `Integrate[.., Method -> "CRCTable"]`, which returns an
  `Abs`-carrying form whose formal derivative does not close. The test's own comment claiming that
  family "now closes cleanly after the number-field Cancel improvements" is stale for the
  `Method -> "CRCTable"` route it exercises -- it does close through the full `Integrate`.
- `dsolve_stress_tests`: `DSolve`UndeterminedCoefficients[y'' - 2 y' + y == Cos[2 x], y, x]`
  returns unevaluated where the test asserts a `List`. The full `DSolve` still solves it
  (`C[1] E^x + C[2] x E^x + Sin[ArcTan[-4, -3] + 2 x]/5`) by another method, so only that
  method's direct entry point regressed -- and before this session.
- `dsolve_tests`: exits 142 = SIGALRM, the 120 s watchdog every test binary gets from
  `test_utils.h`. It reaches ~48 of its tests in 120 s. This is the pre-existing in-suite
  failure recorded in memory (`project_dsolve_tests_m19_insuite_abort`); the suite is simply
  longer than the watchdog now. NOT a slowdown from this work: 60 representative `DSolve`
  calls lifted out of `test_dsolve.c` run in 1.62 s on a pristine HEAD kernel and 1.66 s on
  this one (2 %, at the noise floor of a loaded machine).
- `intrischnorman_tests`: `Integrate[1/Log[x], x]` gives `ExpIntegralEi[Log[x]]` where the test
  wants `LogIntegral[x]` (the same function). Reproduced identically on a pristine HEAD kernel
  AND a pristine `src/internal`, so PRE-EXISTING.
- `moebiusmu_tests` / `primenu_tests`: `MoebiusMu[10^50 + 1]` and
  `PrimeNu[2491230487120948712093481230948273409812734091238]` come back with the wrong parity /
  count inside the test binary (`FactorInteger::nofac: ... composite but no factor was found
  within the search bounds` -- the ECM budget), while a plain `-file` run answers `-1` and `8`
  correctly. Both test binaries rebuilt from a HEAD kernel fail identically: PRE-EXISTING.
- `risch_rde_tower_tests`: asserts `Integrate[E^(Log[x]^2), x]` stays unintegrated; it now returns
  an `Erf` form. Rebuilt from a HEAD kernel with a HEAD `src/internal`: the same `Erf` form, so
  PRE-EXISTING (the assertion is stale, not a regression).
- `dsolve_corpus_tests` needs more than the 20 min cap this run gave it (1204 cases, ~500 in
  17 min); re-run it with `timeout 3600` for a verdict.

---

# `PolynomialGCD` — a native multivariate GCD over a number field (plan, 2026-09-29)

Plan: `/Users/user/.claude/plans/cheeky-seeking-flute.md`.  Target: the `flint_field_gcd`
named as the fix in `MATHILDA_DIVERGENCES.md` A26.

- [x] Spike the FLINT call sequence before building on it (`fq_nmod_mpoly_gcd` in each residue
      field of `M mod p`, CRT back) — verified on the A26 pair for an INERT prime (p = 5, 11)
      and a SPLIT one (p = 7, 17); both reconstruct `x + sqrt2 y` exactly
- [x] `flint_field_gcd` — modular (Encarnación) multivariate GCD over `K = Q(theta)`, with the
      `{G, M}` Gröbner certificate and a `MATHILDA_NO_FIELD_GCD=1` A/B gate
- [x] Radical / `Root` input normalised to one common field and rendered back through the
      PRODUCT BASIS of the caller's own atoms (radicals in, radicals out)
- [x] Two hooks: `poly_gcd_internal`'s inner FLINT fast path (covers `Factor`, `SquareFreeQ`,
      `FactorTerms`, `PolynomialLCM`, `Cancel`/`Together`, Risch/DSolve) and
      `builtin_polynomialgcd`'s FLINT block, ahead of the Phase C/D tower paths
- [x] Safety net: the multivariate classical path CHECKS its answer and returns `1` when it
      does not divide both operands (a pre-emptive refusal was too blunt — see below)
- [x] `KGcd` (`logrewrite.m`) drops its multivariate A26 workaround
- [x] Tests, docs (A26 rewritten, `algebra.md`, changelog, docstring), version `0.230`
- [ ] `flint_field_reduce_core`'s univariate-only restriction — **dropped, not needed**: see below

## Review

**What was actually wrong.** Three defects, of which A26 recorded only the first.

1. *A wrong answer.* `PolynomialGCD[f, g]` returned `g`; `PolynomialGCD[g, f]` returned `f` —
   always the second operand, never a common divisor. `PolynomialGCD[x + a y, (x+a y)(x+2)]`
   returned something of higher degree than its own first operand. Mechanism, traced to source:
   `poly_content` bottoms out in `my_number_gcd` → `get_int_content` (`poly.c:1251-1291`), which
   is INTEGER content and answers 1 for any `K`-coefficient, so both operands enter the PRS
   non-primitive over `K` and the first `pseudo_rem` (`lc(B)*A - lc(A)*B`) vanishes identically.
2. *A missed factor.* `Extension -> Automatic` on RADICAL input reached the Phase D tower path,
   which computes the `Q[gamma,x,y]`-GCD, not the `Q(gamma)[x,y]`-GCD (`qafactor.c:2798` says so
   itself). Correct only while the cofactors are free of algebraic constants; returns 1 otherwise.
   Not recorded in A26 — found by probing, not by reading.
3. *An unreachable engine.* Every explicit `Extension -> <value>` form (a value, a list, `All`,
   `None`) went to the classical path. Only `Extension -> Automatic` reached a field engine.

**Verified against Mathematica 13.2** rather than assumed. That mattered: Mathematica's default
`PolynomialGCD` on `AlgebraicNumber` coefficients returns **1**, not the gcd — so "what the right
answer is" was a question with a surprising answer, and the safety net was designed to match it.

**Why modular/`fq_nmod`.** FLINT has no multivariate gcd over a number field — `gr_mpoly.h`
declares no arithmetic at all (a skeleton: no add, no mul, no gcd) and there is no
`nf_elem_mpoly`. `fq_nmod_mpoly_gcd` is the only multivariate gcd with characteristic-p
coefficients, so Encarnación is both the textbook answer and the *least* code: FLINT does the
multivariate work, this writes reduction, CRT, reconstruction and certification.
`flint_bridge.h:358` already advertised "multivariate via a modular fq_nmod GCD" — that comment
had been aspirational since it was written; it is now true.

**Splitting the residue ring is mandatory, not an optimisation.** For a non-cyclic Galois group —
`Q(sqrt2, sqrt3)`, group `(Z/2)^2` — *no* prime keeps `M` irreducible, since the group has no
element of order 4. Any design that requires an inert prime fails on exactly the compositum
fields this work needs. The spike deliberately exercised both an inert and a split prime before
a line of the engine was written.

**The certificate is one call.** `G` is monic ⇒ its leading monomial is tau-free; `M`'s is
`tau^n`; coprime leading monomials ⇒ `{G, M}` is a Gröbner basis by Buchberger's first criterion
⇒ `fmpq_mpoly_divrem_ideal` *decides* divisibility over `K`. No cofactor reconstruction. So an
unlucky prime or a premature reconstruction costs an iteration and can never produce a wrong
answer.

**Results** (all now correct; previously the first two were non-divisors and the third missed the
factor):

| input | before | after |
|---|---|---|
| `PolynomialGCD[(x+a y)(x+1), (x+a y)(x+2)]` | `(x+a y)(x+2)` | `x + a y` |
| `PolynomialGCD[g, f]` (order swapped) | `f` | `x + a y` |
| `d(x^2+√2y+3), d(x^2+2√2y−1)`, `Ext → Automatic` | `1` | `d` |
| explicit `Extension -> Sqrt[2]`, multivariate | second operand | `x + √2 y` |

`Cancel`/`Together` now match Mathematica in all five option forms — which is why the planned
lift of `flint_field_reduce_core`'s `gens.count != 1` restriction was **dropped**: the
`poly_gcd_internal` hook already reaches `Cancel`, and Mathematica also leaves the default
`AlgebraicNumber` case uncancelled, so there was nothing left to fix. Less code than planned.

**Two bugs found in my own work, both worth recording.**

- *Use-after-free.* On the radical path `theta` is borrowed from the normalised operand, which I
  freed before rendering — the result came back as `x`, a clean-looking non-divisor that the
  certificate had already passed on the *correct* value. Fixed by a single cleanup at the bottom.
- *A safety net that was too blunt.* The first version refused the classical PRS outright
  whenever an algebraic constant was present with more than one variable. That is sound but
  over-broad: it broke `dsolve_m12_stress`, whose answer runs through `Q(Sqrt[17])` and whose
  gcds the PRS was handling correctly. Replaced by a post-check — run the PRS, then verify the
  answer divides both operands and fall back to 1 only when it does not. Precise, and it only
  ever replaces a genuinely wrong answer. Caught by the full suite, not by the targeted tests.
- *Normalising in the wrong basis.* Clearing denominators in `theta`'s power basis and then
  rendering in the caller's radicals injected a junk constant: `x + I y` came back as
  `3 x + 3 I y`, because `Expand` folds `2*I` into `Complex[0,2]`, so the atom set is `{I, 2I}`
  and the primitive element is `theta = 3I`. The compositum did the same via
  `sqrt2 = (theta^3 − 9 theta)/2` → a factor 2. Both are correct *associates*, which is exactly
  what makes them easy to ship. Fix: stay monic (canonical and basis-independent, and already
  what the modular reconstruction produces). Recorded as a memory.

**Performance — the interesting part.** The engine is called from `poly_gcd_internal` on every
gcd and declines on almost all of them, so the decline path is a hot path. Three rounds:

1. First cut: Charlwood 13.1 → **18.0 s** (A40 alone 1.13 → 3.06). Fixed by two cheap gates — a
   structural `fg_has_algebraic` scan before anything expensive, and declining below two
   polynomial variables (univariate is already correct and faster elsewhere; this engine exists
   for the multivariate case nothing else can do). → **13.4 s**.
2. Adding nested-radical support put A40 back to 2.13 s. The obvious suspect was the qqbar
   minpoly lookup in atom collection, so I memoised it: 496 constructions → **19**, and the time
   **did not move at all**. Worth recording — the cheap-looking thing was not the cost.
3. Instrumenting properly settled it: `flint_qqbar_to_number_field_common` was **113 calls
   costing 1.065 s**, while all 43 modular GCDs those calls enabled cost **0.003 s** together.
   Building the field, not computing in it. Caching it by atom set (the C analogue of the .m
   layer's `FieldDataMemo`) → 0.180 s, A40 → **1.247 s**.

Final: **12.5 s**, 49/50, 0 wrong — *faster than the 13.10 s baseline* while fixing the wrong
answers. P8 0.64 → 0.557, P4 1.24 → 1.18, A27 0.96 → 0.888, A28 0.23 → 0.212, A40 1.13 → 1.247.
Both memos were kept: the atom memo does not show on this benchmark but removes an obvious
repeated cost, and the field memo is the whole of the win.

**Also found, left alone** (recorded in the plan, not fixed — out of scope):
- `qafactor.c:4770` builds a three-argument `PolynomialGCD[num, den, S]` intending `S` as the
  variable. `PolynomialGCD` is variadic over *polynomials*, so it computes `gcd(num, den, S)` —
  almost always 1.
- `Modulus` is advertised in `Options[PolynomialGCD]` (`options_builtin.c:607`) but
  `builtin_polynomialgcd` never parses it, so `PolynomialGCD[a, b, Modulus -> 5]` returns
  unevaluated.

# `flint_field_gcd` — stress test and harden to library standard (plan, 2026-09-29)

Plan: `/Users/user/.claude/plans/cheeky-seeking-flute.md`

- [x] Instrument first (`MATHILDA_FIELD_GCD_STATS=1`): per-stage timers + prime/certificate counters
- [x] Raise the coefficient ceiling (62-bit primes, bit budget not prime count)
- [x] Fix the degree-6 decline (upstream, qqbar compositum primitive-element choice)
- [x] Fix the compositum render-back (greedy spanning atom selection)
- [x] Fix the `fg_field_images` allocation-failure leak
- [x] `tests/bench_field_gcd.c` — 22 self-certifying cases + 2 gates, in ctest
- [x] `test_field_gcd_stress_regressions` in `tests/test_algebraicnumber.c`
- [x] Optimise only what the profile condemned (prime choice, prime pool, push_term)
- [x] Adversarial pass (zero/constant operands, deep towers, huge exponents, mixed spellings)
- [x] Docs: A26a, `algebra.md`, changelog, `flint_bridge.h` contract; v0.232 + tag

## Review

**The bar was "would a FLINT developer accept this", so the question asked was not "does it answer
the A26 repros" but "what is its cost curve, where does it silently give up, and can it be made to
crash". Three of the four things found were silent — they all ANSWERED.**

A decline is the failure mode that matters here and it is not a slow answer: the engine returns
NULL, `poly_gcd_internal`'s post-check answers 1, that is a valid common divisor, nothing
downstream complains, and the real gcd is gone. Two whole classes were doing this.

| | v0.230 | v0.232 |
|---|---|---|
| coefficient ceiling | declines past **~831 bits** | no decline at 13,288 bits |
| `[K:Q] = 6` via radicals | declines for **every** generator tried | works |
| compositum `{√2, √3}` | correct but a degree-4 `Root` per coefficient | `1 + x^3 + Sqrt[2] x y + y^2` |
| 628-term operands | 63 ms | 30 ms |
| residue-field work | — | 2.4–2.6× less |
| Charlwood 50 | 12.5 s, 49/50 | 12.5 s, 49/50 |

**The ceiling was arithmetic.** 64 primes × 29 bits = 1856 bits of modulus, and rational
reconstruction needs about twice the coefficient size; measured, 10^250 passed and 10^300 did not.
62-bit primes plus a cap demoted to a grind-backstop fixed it. This is only safe because the
certificate rather than the budget is what makes the answer correct — extra primes cost time, never
correctness.

**The degree-6 decline was upstream and beautifully specific.** `qqbar_express_in_field_esc` tries
64 bits and then refuses to escalate when the generator has degree `<= 6`. Degrees 2–5 resolve
inside 64 bits; 7+ escalate; **6 alone** needs more than 64 bits and is denied. So the compositum's
trial-membership search rejected every candidate multiplier and the whole field failed. Fixed by
choosing the primitive element from DEGREES — `alpha + c*b` always lies in `Q(alpha, b)`, so it
generates the compositum exactly when the degrees agree — which removes the membership test rather
than tuning it. The shared `in_field` was left alone: its degree gate exists because ungated
escalation had regressed the DSolve callers.

**Two suspects from reading the code that measurement cleared.** Recording these because the
reading was persuasive and wrong both times, which is the same lesson this engine taught last round:
- the *certificate* looked like the expensive step — an exact multivariate ideal division over Q run
  after every prime, up to 64 times. It is ~5 µs against ~150 µs for one prime's residue gcds.
  Gating it on a stabilised reconstruction therefore **cost an extra prime of the dominant work**
  and made the small cases 2× slower; reverted to certifying each distinct candidate once,
  immediately, so an easy input finishes on the first prime.
- hoisting the per-component input reduction (`fg_lift` runs once per irreducible factor) looked
  like an `r`× saving. Lift is 8% of the time. Dropped, deliberately, unimplemented.

Only the `O(L²)` keyed-insert image construction was worth replacing (`push_term` + `sort_terms`).

**The win came from somewhere the plan had not listed at all**: `fg_gcd_mod_p` runs one multivariate
gcd per irreducible factor of `M mod p`, so the per-prime cost *is* how far the prime splits.
Choosing least-split primes cut the residue work 2.4–2.6×. The direction was A/B'd rather than
assumed, because it is genuinely not obvious — a split prime has single-word coefficient arithmetic
where an inert one does degree-`n` polynomial arithmetic, so this only wins if FLINT's per-call cost
dominates. It does, at these sizes: choosing the *most*-split prime was 1.6× slower at `n = 2` and
2.6× slower at `n = 8`. Noted in the code that the balance is worth re-measuring on much larger
operands.

**One self-inflicted trap, caught by the harness abort-trapping on its own 903-digit row:**
`is_associate` interpolated each operand twice into a buffer sized for one copy each. My bug, in the
test code, and a useful reminder that a fixed `char buf[2048]` in `check_gcd` would have silently
truncated the large-coefficient assertions into something that no longer tested what it named — both
are now sized from the strings.

**Plan item 1 was wrong and is dropped.** I expected a reachable hard abort: nothing on the `fg_`
path guards `_get_term_exp_ui`, and FLINT's failure there is `flint_throw`, which is `FLINT_NORETURN`.
It is unreachable — `to_mpoly` gates exponents on `EXPR_INTEGER`, i.e. `int64`, and every
non-negative `int64` fits a `ulong`; `x^(2^63)` becomes a bigint and is refused before FLINT sees
it. Verified against four adversarial exponent shapes, all of which decline cleanly. The invariant
is now documented rather than guarded with dead code.

**Verification.** 22/22 harness cases certify with the scaling ratio at 2.33 for 1.94× the terms
(limit 3.0); `algebraicnumber_tests`, `numberfield_tests`, `flint_bridge_tests` pass; Charlwood
49/50 at 12.5 s, unchanged; `make check-c99` and `make check-messages` green.

**Also found, pre-existing, recorded not fixed** (identical at v0.230 and with
`MATHILDA_NO_FIELD_GCD=1`, so outside this engine):
- `PolynomialGCD[0, f]` with algebraic coefficients answers 1; Mathematica answers `f`.
- A coefficient mixing an inexact real with an algebraic constant returns a garbage near-zero float
  instead of declining: `PolynomialGCD[Expand[(x + 1.5 Sqrt[2] y)(x+1)], ...]` → `3.71618e-16`. The
  pure-float case is correct, so it is specifically the mixture.
- Mixed spellings of one field, and two distinct `AlgebraicNumber` generators, still answer 1
  (`field_scan` reports a conflict rather than building the compositum).
