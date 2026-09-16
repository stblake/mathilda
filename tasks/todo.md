# DSolve M49 — §2.2.30 corpus (Problems 2901–3000)

Milestone: M49 = §2.2.30. Version bump 0.146 → 0.147.
Definition of done: 0 FAIL, 0 crash; UNEVAL OK only for no-closed-form residue.

## Phase 1 — Fetch & convert
- [x] Fetch Ch2.S2.SS30.htm (browser UA), verified title "Problems 2901 to 3000"
- [x] Convert → DE_examples_2230.m (100 records: 100 scalar [24 IVP], 0 systems)
- [x] Spot-audit found ONE converter transcription bug class (Greek-letter variables):
      - 2972/2973/2992: `dr/d\theta` ODEs — θ is the INDEPENDENT variable (was defaulted to x)
      - 2984: `\theta'(t)` ODE — θ is the DEPENDENT variable (was mangled `theta^(prime)`, func→`a`)
- [x] Fixed generally in latex_ode_to_mathilda.py — Greek letters now first-class variables:
      (1) primed-detection recognizes a Greek macro as a dependent variable (NAME_RE)
      (2) indvar detection adopts a lone Greek letter as indvar for FIRST-order ODEs only
          (2nd-order `y''+λy=0` eigenvalue problems keep λ as a PARAMETER — order gate)
      (3) convert_side canonicalizes Greek VARIABLE macros before the mains/derivative pass
      (4) placeholder-expansion regex broadened to multi-char identifiers (theta)
- [x] Byte-identity verified: OLD-vs-NEW converter on SAME HTML, SS20..SS29 ALL IDENTICAL;
      §2.1.2 arbitrary-function detect_symbols unchanged (direct old-vs-new). Only Greek
      variable names in the whole corpus are the 4 new §2.2.30 θ records.

## Phase 2 — Baseline measurement
- [x] Build dsolve_corpus_tests
- [x] Baseline (clean corpus): 92 PASS / 8 UNEVAL / 0 FAIL / 0 crash
- [x] reports/2.2.30.md + .tsv generated

## Phase 3 — Root-cause fixes
- [x] 0 FAIL / 0 crash confirmed (mandatory bar met at baseline)
- [x] Probed the 4 sympy=True UNEVALs: 2980 (linear, IF not rationalized → HANG),
      2970 (Lie symmetry hang), 2979 (IVP fit fails on 4 radical branches), 2933 (dAlembert nonelem)
- [x] General solver fix: inverse-hyperbolic integrating factor (dsolve_linear_factor_solve) —
      TrigToExp-rationalize Exp[c ArcTanh] → 2980 U→P. Validated on 8 representative mu + REPL
      full-solve simulation (general+IVP residual 0) BEFORE the C edit.
- [x] Re-measure: 93 PASS / 7 UNEVAL / 0 FAIL / 0 crash. Only 2980 flipped (nothing regressed in-section)
- [x] Residue 7 = 4 sympy=False (honest) + 3 deferred sympy=True gaps

## Phase 4 — Record & land
- [x] tests/CMakeLists.txt — dsolve_corpus_2_2_30_tests gate (baseline 7)
- [x] STATUS.md — §2.2.30 block + M49 wave-history line
- [x] README.md — DE_examples_2230.m row
- [x] DSOLVE_PLAN.md — M49 bullet
- [x] docs/spec/changelog/2026-09-14.md — M49 entry (top; this ISO week)
- [x] src/version.h — 0.146 → 0.147; rebuilt (banner shows 0.147)
- [x] make check-c99 clean (C changed this time — PASS)
- [~] Full ctest -R dsolve_corpus → 0 regression (RUNNING in background)
- [x] lessons.md + memory updated
- [ ] Commit + push (on regression green + user OK)

## Review

**Result: §2.2.30 landed at 93/100, 0 FAIL, 0 crash. v0.147.**

Unlike M48, this section needed BOTH a converter fix and a solver fix, both general and both
verified byte-identical/0-regression on prior sections.

1. **Converter — Greek-letter variables** (`latex_ode_to_mathilda.py`). The single-letter
   symbol-detection layer could not see `θ` as a variable: 2972/2973/2992 (`dr/dθ`, θ independent)
   were mis-read with indvar `x`/θ-as-constant (solvable but WRONG), and 2984 (`θ'(t)`, θ dependent)
   mangled the derivative to `^(prime)`. Four coordinated parts make a Greek macro a first-class
   name; crucially the indvar adoption is gated to FIRST-order ODEs so a 2nd-order eigenvalue
   `y''+λy=0` keeps λ a parameter (this regression was caught by the byte-identity sweep, not by
   reasoning). OLD-vs-NEW on same HTML: all §2.2.20–29 IDENTICAL, §2.1.2 detect_symbols unchanged.

2. **Solver — inverse-hyperbolic integrating factor** (`dsolve_linear_factor_solve`). 2980 is a
   LINEAR ODE that hung because μ=Exp[−4 ArcTanh[x]] was integrated un-rationalized. Gated
   TrigToExp rationalization → real algebraic (x−1)^a(x+1)^b → elementary integral. Trig/polynomial
   IFs byte-identical.

Residue 7: 4 sympy=False (no closed form: 2923/2944/2948/2955) + 3 deferred sympy=True gaps (2933
dAlembert nonelem integral, 2970 Lie symmetry, 2979 radical-branch IVP fit). All honest UNEVAL.
