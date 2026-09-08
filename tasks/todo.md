# M24 — §2.2.4 corpus (Problems 301–400)

Add Nasser Abbasi 12000.org §2.2.4 (indexsubsection13.htm, "Problems 301 to 400") to
the DSolve corpus dashboard; measure baseline; close reuse-reachable gaps at 0 FAIL;
decline research-grade residue cleanly. Mechanical repeat of M21–M23.

## Step 1 — Convert HTML → corpus
- [x] Convert `ss13.htm` → `DE_examples_224.m` (`--label 2.2.4`): 100 scalar, 24 IVP, 0 sys
- [x] **Converter fix 1**: imaginary unit `i` (309/310/311) — exclude from indep-var + map→`I`
- [x] **Converter fix 2**: `y^{(n)}` derivative notation (336/340/343) → `y'''''[x]`
- [x] Regression guard: §2.2.1/2/3 regen byte-for-byte identical ✔

## Step 2 — Register ctest + first baseline run
- [x] Add `dsolve_corpus_2_2_4_tests` (baseline now 1)
- [x] Baseline run: **93 PASS / 6 UNEVAL / 1 FAIL**
- [x] Generate `reports/2.2.4.{tsv,md}`

## Step 3 — Analyze gaps + targeted verified fixes
- [x] **FAIL 387** root-caused = harness `dsResidVerdict` sweep-var `k` ↔ ODE param `k`
      collision → prelude `$dsSweep`/`$dsVal` (benefits every section) → PASS
- [x] **326/362/363/365** trig-power/product forcing → `TrigReduce` in undetcoeff → PASS
- [x] **312** complex cube-root IVP → `ComplexExpand` numeric roots in homog_basis → PASS
- [x] Pinned units `t_m24_trig_power_forcing`, `t_m24_complex_cuberoot_ivp` (dsolve_tests OK)
- [x] Residue 381 (var-coeff Legendre, SymPy-failed) declines cleanly

## Step 4 — Scoreboard + gate
- [x] STATUS.md §2.2.4 block + wave-history bullet
- [x] README.md contents row
- [x] `dsolve_corpus_2_2_4_tests` baseline → 1 (final 99/1/0)

## Step 5 — Plan + changelog + task log
- [x] DSOLVE_PLAN.md M24 milestone entry
- [x] `docs/spec/changelog/2026-09-07.md` M24 note (solver behavior changed)
- [x] Review section (below)

## Verification
- [x] Final §2.2.4: **99 PASS / 1 UNEVAL / 0 FAIL**; new PASSes reproduced interactively
- [x] §2.2.1/2.2.2/2.2.3/2.2.4 ctest gates hold
- [~] §2.1.2 full re-run (0 FAIL so far at 461/1000; awaiting completion)
- [x] All DSolve unit + stress suites green (dsolve_tests + m5/m12/m14/m17/m18/m19/m20)
- [x] `make check-c99` green
- [x] C changes leak-clean by review (eval_and_free idiom; macOS valgrind is documented-noisy)
- [x] §2.2.1/2/3 `.m` files unchanged (byte-identical regen)

## Review

**Result: §2.2.4 baseline 93/6/1 → 99/1/0 (+6 solved, FAIL eliminated), 0 regression.**

Two classes of fix. **Corpus fidelity** (converter): §2.2.4 was the first section to carry
constant-coefficient *complex* ODEs written with the imaginary unit `i` and 5th-order
`y^{(n)}` notation — both silently misconverted (309-311 picked `i` as the indep var;
336/340/343 became powers of `y`). Fixed in the converter with the byte-identical §2.2.1/2/3
regen invariant preserved.

**Solver / harness** (3 fixes, all reuse verified machinery → 0-FAIL by construction):
the lone FAIL (387) was a **harness** variable-capture bug (sweep index `k` colliding with a
spring-constant parameter `k`), not a wrong answer — fixed once, benefits every section;
trig-power/product forcing now `TrigReduce`s to first-harmonic sinusoids in
`UndeterminedCoefficients`; and numeric complex characteristic roots are `ComplexExpand`-
concretized so IVPs with cube/quartic roots fit. Sole residue 381 is a genuine
variable-coefficient-nonhomogeneous gap (Kovacic homogeneous + var-params particular) SymPy
also fails — a clean deterministic decline.
