# Bose–Einstein / Fermi–Dirac / PolyLog half-line Mellin family

**Goal.** Close the statistical-mechanics half-line integrals
`∫₀^∞ x^(s−1)/(e^(cx) ∓ 1) dx` (Planck / Debye / Fermi–Dirac / Sommerfeld),
plus the general fugacity kernel `1/(e^(cx)+γ)` landing on `PolyLog`.

**Design decision (elegance).** These are NOT a new method / file / cascade
stage. The kernel `1/(e^(cx)∓1)` is a *missing recognizer* in the existing
Mellin-transform table in `src/calculus/integrate_ramanujan.c`. Adding one
recognizer to `try_recognizers` makes it work through both the Automatic cascade
and `Method -> "RamanujanMasterTheorem"`, reusing all existing machinery: term
splitting, monomial `x^k` substitution, strip gating, `ConditionalExpression`
output, term-by-term summation. Zero new plumbing.

## The mathematics (correct-by-construction)

Kernel `1/(e^(cx) + γ)`, c>0, γ real with −1 ≤ γ ≤ 1. With `u = e^(−cx)`:
```
1/(e^(cx)+γ) = u/(1+γu) = (−1/γ) Σ_{j≥1} (−γ)^j u^j    (|γ|e^(−cx) < 1)
∫₀^∞ x^(s−1) u^j dx = Γ(s) (jc)^(−s) = Γ(s) c^(−s) j^(−s)
⇒ M(s) = Γ(s) c^(−s) (−1/γ) Σ_{j≥1} (−γ)^j/j^s
        = Γ(s) c^(−s) (−1/γ) PolyLog(s, −γ)
```
- **Bose** γ = −1: `M = Γ(s) c^(−s) ζ(s)`         (via `PolyLog(s,1)`)
- **Fermi** γ = +1: `M = Γ(s) c^(−s) (−PolyLog(s,−1)) = Γ(s) c^(−s) η(s)`

**Convergence strip.** x→∞: e^(−cx) decay ⇒ any s. x→0: denominator `1+γ` ≠ 0
for −1<γ≤1 ⇒ `Re s > 0`; for γ = −1 the denominator ~ c·x ⇒ `Re s > 1`.
No interior pole on (0,∞) when |γ| ≤ 1 (the `e^(cx) = −γ` root is at
x = ln(−γ)/c ≤ 0). Gate exactly these.

**Verified in REPL (2026-07-09):**
`Γ[4]·PolyLog[4,1] = π⁴/15` (Debye), `Γ[4]·(−PolyLog[4,−1]) = 7π⁴/120` (Fermi),
`PolyLog[1,−1] = −Log[2]` ⇒ `∫₀^∞ 1/(e^x+1) = Log 2` at s=1 with NO 0·∞.

## Tasks

### 0. FIRST — fix the never-wrong violation in the residue family
- [ ] `Integrate[Sqrt[x]/(1+x)^2, {x,0,Infinity}]` returns **`0`** (correct `π/2`).
      Attributed: `Method->"Residue"` → `0`; `Method->"RamanujanMasterTheorem"` → `π/2`.
      Bug is in `residue_family_mellin` (`src/calculus/integrate_residue.c`),
      mis-summing a **double pole** at x=−1 for non-integer p=1/2. Runs before RMT
      in the cascade, so it shadows the correct answer.
- [ ] Fix: correctly handle order≥2 poles in the keyhole residue sum **or**, if
      the family cannot, detect order>1 poles and return `NULL` (fall through to
      RMT, which is already correct). Prefer correctness; NULL-gate is the floor.
- [ ] Regression test in `tests/test_integrate_residue.c` (and confirm the
      Automatic cascade now returns `π/2`).

### 1. Recognizer `rec_expgeom` in `integrate_ramanujan.c`
- [ ] Add `static bool rec_expgeom(K, x, sv, M, P)` next to `rec_polylog`:
      - K must be `Power[base, −1]`, base a `Plus`.
      - Split base into x-free constant `γ` and exactly one x-term
        `A·Exp[arg]` (reuse `exp_exponent`; extract `c1 = D[arg,x]` x-free,
        `c0 = arg − c1 x` folded via `A_eff = A·e^(c0)`).
      - Require `c1 > 0` and `A_eff > 0` (prove via Simplify).
      - `γ' = γ / A_eff`; require provably real, `−1 ≤ γ' ≤ 1`.
      - `M = (1/A_eff) · Γ(sv) · c1^(−sv) · (−1/γ') · PolyLog(sv, −γ')`.
      - `P = And[c1>0, sv>0]`; if `is_zero_now(γ'+1)` also `And[..., sv>1]`.
- [ ] Insert `if (rec_expgeom(K,x,sv,M,P)) return true;` into `try_recognizers`
      (after `rec_polylog`). The monomial `x^k` wrapper in `dispatch_kernel`
      then gives `1/(e^(√x)−1)` etc. for free.
- [ ] Canonicalize symbolic-s Bose output: post-rewrite `PolyLog[s,1] :> Zeta[s]`
      in the result (Simplify does NOT do this for symbolic s). `−γ'=−1` (Fermi)
      stays `−PolyLog[s,−1]`, which evaluates cleanly at integer s incl. s=1.

### 2. Docs / attributes
- [ ] Extend the `Integrate\`RamanujanMasterTheorem` docstring kernel list with
      "exponential-geometric 1/(e^(cx)+γ) → Γ·PolyLog (Bose–Einstein/Fermi–Dirac)".
      Keep terse (no examples in docstrings).
- [ ] `docs/spec/builtins/calculus.md`: add the kernel + the two headline results.
- [ ] `docs/spec/changelog/2026-07-06.md` (Monday of this ISO week): changelog note.

### 3. Tests — `tests/test_integrate_ramanujan.c`
- [ ] `∫₀^∞ x^3/(e^x−1) = π⁴/15`     (Debye)
- [ ] `∫₀^∞ x^3/(e^x+1) = 7π⁴/120`   (Fermi)
- [ ] `∫₀^∞ 1/(e^x+1) = Log 2`        (Fermi, s=1, the 0·∞-avoidance case)
- [ ] `∫₀^∞ x/(e^x−1) = π²/6`         (Bose, s=2)
- [ ] `∫₀^∞ x^(s−1)/(e^x−1), Assumptions->s>1  = Γ[s] Zeta[s]`  (symbolic, canonical)
- [ ] scaled rate: `∫₀^∞ x^3/(e^(2x)−1) = π⁴/240`  (c=2 ⇒ c^(−s) factor)
- [ ] monomial: `∫₀^∞ 1/(e^(√x)−1)` closes via the x^k wrapper
- [ ] out-of-scope stays unevaluated (no false positives): e.g. `1/(e^x−2)`
      (|γ'|=2 > 1 → NULL), divergent `1/(e^x−1)` at s=1 (Bose, strip Re s>1).
- [ ] add to `tests/CMakeLists.txt` COMMON_SRC only if a new .c file is created
      (none expected — edits are in-file).

### 4. Verify
- [ ] Automatic cascade: all §3 cases return the closed form (residue no longer
      shadows them; if it ever matches these it must agree, not return 0).
- [ ] Regression: existing `test_integrate_ramanujan` + `test_integrate_residue`
      pass; spot-check prior working integrals unaffected.
- [ ] valgrind clean on `∫₀^∞ x^3/(e^x−1)` (diff vs `Sin[1.0]` baseline noise).
- [ ] "Would a staff engineer approve?" — the strip gate makes each result
      unconditionally correct; no NIntegrate crosscheck (project rule).

### 5. Frullani (IN SCOPE) — `frullani_try` in `integrate_ramanujan.c`
Runs as a whole-integrand pre-pass at the top of `integrate_ramanujan_try`
(BEFORE Expand — each split term is individually divergent, so it must match
whole). `∫₀^∞ (f(ax)−f(bx))/x dx = (f(0⁺)−f(∞)) Log(b/a)`, a,b>0.
- [ ] `G = Simplify[x*F]`; require `Plus` of two terms `t1, t2`.
- [ ] Find x-free `ρ>0, ρ≠1` with `Simplify[(t1 /. x->ρ x) + t2] == 0`
      (numeric-anchor at x=1 via Solve to get candidate ρ, then symbolic verify
      for all x — robust even for non-injective f like Cos).
- [ ] `f0 = Limit[t1, x->0]`, `finf = Limit[t1, x->Infinity]`; both must be
      finite (else NULL — relies on Limit, which handles the exp cases).
- [ ] Return `(f0 − finf) * Log[ρ]`. Test: `(e^(−2x)−e^(−5x))/x = Log[5/2]`.

### 6. Log-weighted Mellin (IN SCOPE) — extend `split_term` / `mellin_term`
`∫₀^∞ x^(s−1) Log[x]^k R(x) dx = ∂^k/∂s^k M_R(s)` (Log^k is dominated by x^±ε,
so R's OPEN convergence strip transfers unchanged).
- [ ] In `split_term`, detect factors `Log[x]` / `Power[Log[x], k]` (arg exactly
      the variable x, k a positive integer); accumulate `kw` (total log power),
      do NOT add them to `kernels`. Must NOT match `Log[1+x]` (that's rec_log).
- [ ] In `mellin_term`, if `kw>0`: build `M_R(sv)` via `dispatch_term` in a fresh
      spectral symbol, then `M = (D[M_R, {sv,kw}] /. sv->s) * C`; strip = M_R's P.
- [ ] Test: `∫₀^∞ Log[x]/(1+x^2) = 0` (symmetry ⇒ M_R'(1)=0);
      `∫₀^∞ Log[x]/(1+x)^2 = 0`; a nonzero case e.g. `∫₀^∞ x Log[x] e^(−x) = 1−γ`.

## Review (2026-07-09) — DONE

All four deliverables landed and verified.

- **§0 Residue double-pole fix** (`src/calculus/integrate_residue.c`): keyhole
  sum is now over `Res[z^(s-1) R(z), z_k]` (full integrand), via the shift
  `w=z−z_k` + analytic `(1+w/z_k)^(s-1)` through the existing `residue_compute`.
  Simple poles keep the fast `z_k^(s-1) Res(R)` path (gated by `Q'(z_k)≠0`) — this
  was needed because the general path is much heavier and every pre-existing case
  is a simple pole. `√x/(1+x)² → π/2`, `√x/(1+x)³ → π/8`, `x^{1/3}/(1+x)² →
  2π/(3√3)`. Regression tests in `test_mellin`.
- **§1 Bose/Fermi/PolyLog** (`rec_expgeom` in `integrate_ramanujan.c`): all cases
  pass — numeric, symbolic (`Γ[s]Zeta[s]`, tight `s>1` strip), fugacity, monomial
  `1/(e^√x−1)`, and out-of-scope/divergent decline cleanly.
- **§5 Frullani** (`frullani_try`): concrete + symbolic scales. Symbolic needed a
  same-sign gate (`k1,k2` proved separately in `>0` orientation) — Simplify can't
  discharge `b/a>0` or a normalised `-a<0`.
- **§6 Log-weighted Mellin**: `Log[x]^k` weight = `∂^k_s` of the base transform;
  threaded through `split_term`/`mellin_term`. `Log[x]/(1+x²)→0`,
  `x Log[x] e^{-x}→1-γ`. (`Log²/(1+x²)` is exact but shows `¼π PolyGamma[1,1/2]` —
  pre-existing PolyGamma reduction gap; value verified numerically.)

**Verification**: `integrate_ramanujan_tests` 0 FAILs; `integrate_residue_tests`
0 FAILs (with the pre-existing `test_rectangular` hang excluded — see below).
Automatic cascade returns all headline results. Valgrind: my files appear in
**zero** leak stacks (13,872B/429 blocks vs. the 13,440B/420 trivial `Sin[1.0]`
baseline; delta traces to pre-existing cascade machinery).

**Pre-existing issues found, NOT introduced here:**
- ~~`integrate_residue_tests` times out on `test_rectangular`'s
  `Integrate[Exp[a x]/(Exp[x]-1), {x,-∞,∞}, Assumptions->0<a<1]`~~ **FIXED
  2026-07-09** (see below). Suite now completes in ~4.6s, 0 FAILs.
- `PolyGamma[1,1/2]` does not reduce to `π²/2` (Simplify/FunctionExpand gap).

## Follow-up (2026-07-09): residue-suite timeout FIXED

The `test_rectangular` 130s+ hang was a *divergent* whole-line integral (real
axis pole at x=0) that DiffUnderInt pursued into a non-terminating
`Integrate[x^k f]` escalation + `Exp→Cosh/Sinh` Risch-Norman expression swell.
- `whole_line_divergent_pole` gate in `integrate_diffunderint.c`: declines
  two-sided-infinite integrands with a real denominator pole (numeric
  sign-change scan; `Solve[…,Reals]` is unreliable on `Exp[x]±1`). Half-line /
  finite untouched. Regression test in `test_declines_cleanly`.
- pmint symbolic-RowReduce guards in `intrischnorman.c`: size cap (>180 rows /
  >8000 leaves; legit max ~132×37 / 5455) + per-call wall-clock budget.
- Result: 130s+ → <1s (still unevaluated, correct). Residue suite green (4.6s).
  Valgrind clean (new code in zero leak stacks). No regressions in
  diffunderint / ramanujan / residue suites. `intrischnorman_tests` line-41 fail
  confirmed pre-existing (identical at HEAD).

## Deferred (genuinely later)
- Decaying-exp variant `1/(1−e^(−x))` (c1<0) — normalize to the growing form.
- `Log[x^m]` weight (currently only `Log[x]`); interaction with the x^k monomial
  substitution wrapper.
