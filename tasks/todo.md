# DSolve M23 — §2.2.3 corpus (Problems 201–300) + coverage

## Wave 0 — Corpus infrastructure (DONE)
- [x] Fetch §2.2.3 HTML (indexsubsection12.htm = "Problems 201 to 300", Table 2.19)
- [x] Generate `DE_examples_223.m` (100 scalar, 35 IVPs, 0 systems)
- [x] Regen §2.2.1/§2.2.2 byte-identical (converter unchanged → no regression)
- [x] Verify corpus parses through binary (harness ran clean)
- [x] Measure baseline: **98/100 PASS, 0 FAIL, 2 UNEVAL** (204, 232)
- [x] Register `dsolve_corpus_2_2_3_tests` in `tests/CMakeLists.txt` (baseline 1)
- [x] reports/2.2.3.{tsv,md}; STATUS.md §2.2.3 block + README row

## Gap analysis (2 non-PASS)
- **2.2.3-204** `9√x y^(4/3) − 12 x^(1/5) y^(3/2) + (8 x^(3/2) y^(1/3) − 15 x^(6/5)√y) y' == 0`
  — EXACT (M_y==N_x), potential `F = 6 x^(3/2) y^(4/3) − 10 x^(6/5) y^(3/2)`.
  Explicit `Solve[F==C, Y]` HANGS on the mixed fractional powers (4/3, 3/2). The
  implicit entry returns `F(x,y[x])==C[1]` (as Maple/Mma do) — just needs a gate.
- **2.2.3-232** `y y'' == 6 x^4` — Emden–Fowler `_with_linear_symmetries`; declines
  cleanly (bounded UNEVAL). Has particular soln y=±x^3. Investigate M12 route.

## Wave 1 — Exact radical potential → implicit (fixes 204) — DONE
- [x] Gate explicit `ds_solve` in `dsolve_exact_try` on `ds_is_rational_in(Fpot, Yn)`
- [x] Verify 204 solves (implicit `F(x,y[x])==C[1]`); explicit path preserved (Rule)
- [x] Anti-overfit unit `t_m23_exact_radical` (impl-fn-rule verify + explicit guard)

## Wave 2 — 2.2.3-232 (Emden–Fowler nonlinear 2nd-order) — RESIDUE
- [x] Investigated: scaling symmetry `X=x∂ₓ+3y∂_y` reduces to autonomous
      `r r''+5r r'+6r²==6` → Abel 2nd-kind first-order (`r p p'==6−5rp−6r²`),
      non-elementary (deferred-M13 AIR). Documented as bounded UNEVAL, 0 wrong answer.

## Wrap-up
- [x] DSolve unit + all stress suites + §2.2.1/§2.2.2/§2.2.3 gates green; check-c99 green
- [~] §2.1.2 gate running (regression check for the Exact routing change)
- [x] DSOLVE_PLAN.md M23 entry + docs/spec/changelog + STATUS.md + calculus.md spec

## Review
**Result: §2.2.3 98→99/100 (0 FAIL, +1). Baseline was already high (elementary
Table 2.19: const-coeff any-order + Euler + basic first-order — all covered).
One surgical solver fix (Exact radical potential → implicit first integral, the
radical twin of M22's transcendental-exact wave), reusing the verified implicit
substrate. Sole residue 232 (Emden–Fowler → Abel 2nd kind). Converter unchanged
→ §2.2.1/§2.2.2 byte-identical (no regression). Delivery: uncommitted (user reviews).**
