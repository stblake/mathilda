# DSolve M48 — §2.2.29 corpus (Problems 2801–2900)

Milestone: M48 = §2.2.29. Version bump 0.145 → 0.146.
Definition of done: 0 FAIL, 0 crash; UNEVAL OK only for no-closed-form residue.

## Phase 1 — Fetch & convert
- [x] Fetch Ch2.S2.SS29.htm (browser UA), verified title "Problems 2801 to 2900"
- [x] Convert → DE_examples_2229.m (100 records: 75 scalar [21 IVP], 25 systems)
- [x] Spot-audit found THREE converter transcription bugs (fixed generally):
      - 2806/2807: 4-var systems x,y,z,h — `h` dropped from function list (ARBFUN)
      - 2824: `\textit{x\_}1` italic-glued subscript → `x_1^(prime)` (mangled)
      - 2890: `x-k\sqrt{...}` → bogus `kSqrt[...]` glued head
- [x] Byte-identity verified: §2.2.27/§2.2.28 old-vs-new converter (same HTML) IDENTICAL;
      §2.1.2 arbitrary-function detect_symbols unchanged (direct comparison)

## Phase 2 — Baseline measurement
- [x] Build dsolve_corpus_tests
- [x] Baseline (clean corpus): 88/100 PASS, 0 FAIL, 0 crash, 12 UNEVAL
- [x] reports/2.2.29.md + .tsv generated

## Phase 3 — Root-cause fixes (converter-only; no solver change needed)
- [x] `\sqrt` letter-juxtaposition glue fix (_implicit_mult) — general; also fixed
      latent priors 2.2.24-2360 + 2.2.26-2536 (both P→P, faithful eq; gates unchanged)
- [x] `\textit{x\_}N` italic-glued subscript fix (normalize_subscripts) → 2824
- [x] System-variable-vs-arbitrary-function fix (detect_symbols, start-anchored) → 2806/2807
- [x] Re-run §2.2.29: 88/100, 0 FAIL. Residue 12 all no-closed-form (11 sympy=False + 2819 elliptic)
- [x] Regenerate corrected priors 2224/2226 (diff = ONLY the tSqrt→t Sqrt fix)

## Phase 4 — Record & land
- [x] tests/CMakeLists.txt — dsolve_corpus_2_2_29_tests gate (baseline 12); 2224/2226 comments
- [x] STATUS.md — §2.2.29 block + M48 wave-history line
- [x] README.md — DE_examples_2229.m row
- [x] DSOLVE_PLAN.md — M48 bullet
- [x] docs/spec/changelog/2026-09-14.md — M48 entry (top)
- [x] src/version.h — 0.145 → 0.146; rebuilt (banner shows 0.146)
- [ ] make check-c99 clean (no C changed this milestone, but run anyway)
- [ ] Full ctest -R dsolve_corpus → 0 regression (running)
- [x] Review section (below)

## Residue 12 (all no elementary closed form)
- 7 nonlinear systems 2811/2813–2818 (sympy=False)
- 4 autonomous 2nd-order 2820/2821/2822/2823 (Duffing/hyperelliptic, sympy=False)
- 2819 (z''+z³==0, sympy=True; elliptic-integral implicit solution, DSolve times out — deferred)

## Review

**Result: §2.2.29 landed at 88/100, 0 FAIL, 0 crash (scalars 70/75). v0.146.**

Unlike prior milestones, **no DSolve solver change was needed** — the existing first-order +
constant-coefficient linear-system stack solves the whole section out of the box. The entire
milestone was **converter transcription fidelity**: a spot-audit of the generated corpus caught
three distinct bugs, each fixed generally in `tools/latex_ode_to_mathilda.py`:

1. **`\sqrt` letter-juxtaposition glue** — a variable letter juxtaposed with `\sqrt` glued to the
   `Sqrt` head as a bogus single symbol (`kSqrt[...]`, `tSqrt[...]`). This is a GENERAL fix: it
   also corrected two latent wrong-equation records in prior sections (2.2.24-2360, 2.2.26-2536).
   Restricted to a LETTER lookbehind so digit-glued (`2\sqrt{x}`, already multiplication) stays
   byte-identical.
2. **`\textit{x\_}N` italic-glued subscript** — a source-side rendering inconsistency (2824) that
   defeated subscript+prime handling; folded to the clean `x_{N}` early.
3. **System variable `h` vs arbitrary function** — the subtlest: `h∈ARBFUN={f,g,h}` was always
   excluded from the function list, but 2806/2807 use `h` as a genuine 4th system variable. The
   discriminator is start-anchored: a symbol that HEADS its own derivative row is a dependent
   variable, while `f'(x)` embedded in a scalar ODE's body (§2.1.2) stays arbitrary. This required
   care to not break §2.1.2's arbitrary-function Abel/Riccati records — verified via direct
   old-vs-new `detect_symbols` comparison.

**Verification of byte-identity** (the discipline for shared-path converter changes): rather than
re-fetch every prior section (the live site drifts and §2.1.2's old tex4ht URL is now a redirect
stub), I ran the **old (HEAD) converter and the new converter on the SAME fetched HTML** and
diffed — isolating my change from site drift. §2.2.27/§2.2.28 came back IDENTICAL; §2.1.2's
arbitrary-function classification was confirmed unchanged by a direct function-level comparison.
The two latent-prior corrections (2360/2536) were the ONLY diffs in those sections and both stayed
PASS, so their gates are unchanged (9 / 11).

**Residue 12** is honest: 11 `sympy=False` (7 nonlinear systems + 4 autonomous 2nd-order with no
elementary closed form), plus 2819 (`z''+z³==0`) whose closed form is elliptic/implicit — DSolve
times out at the 8 s harness bound, and an implicit/elliptic solution is not back-substitution-
verifiable anyway. Deep elliptic-ODE work, deferred per the Standard-milestone scope.

**Verification:** `./Mathilda -v` → 0.146; full `ctest -R dsolve_corpus` (pending confirmation);
`make check-c99` (no C changed, run as backstop). Committed to main on user request.
