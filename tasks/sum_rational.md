# Stage 5 — `Sum`Rational`: infinite rational-function summation

## Goal

Close `Sum[p(i)/q(i), {i, imin, Infinity}]` for rational summands in the index,
producing the digamma/zeta closed forms Mathematica gives. Target outputs
(verbatim from WL, our acceptance set):

| Input | Expected |
|---|---|
| `Sum[1/i^2, {i,1,Infinity}]` | `Pi^2/6` |
| `Sum[1/i^3, {i,1,Infinity}]` | `Zeta[3]` |
| `Sum[1/i^4, {i,1,Infinity}]` | `Pi^4/90` |
| `Sum[1/(i (i^2+1)), {i,1,Infinity}]` | `1/2 (2 EulerGamma + PolyGamma[0,1-I] + PolyGamma[0,1+I])` |
| `Sum[1/(i^2 (i^2+1)), {i,1,Infinity}]` | `1/6 (3 + Pi^2 - 3 Pi Coth[Pi])` |
| `Sum[1/((i^2+3 i+1)(i^2+1)), {i,1,Infinity}]` | four `PolyGamma[0, …]` at radical/complex args |
| `Sum[1/((i^2+3 i+2)(a i^2+b)), {i,1,Infinity}]` | EulerGamma terms + `PolyGamma[0, …]` at symbolic poles |

## The algorithm (one master identity)

Convergence requires `deg q >= deg p + 2`. Decompose into **linear** partial
fractions over the denominator's splitting field:

```
p(i)/q(i)  =  Σ_j Σ_{k=1}^{m_j}  c_{j,k} / (i − ρ_j)^k
```

Then sum each term from `i = imin` to `∞`:

- **k ≥ 2:**  `Σ 1/(i − ρ)^k  =  ((−1)^k/(k−1)!) · PolyGamma[k−1, imin − ρ]`
  ( = `Zeta[k, imin − ρ]` ).
- **k = 1:**  individually divergent; convergence ⇒ `Σ_j c_{j,1} = 0`. Closed
  form: `−Σ_j c_{j,1} · (PolyGamma[0, imin − ρ_j] + EulerGamma)`. Keep the
  `EulerGamma` per term (matches WL; it cancels formally because Σc=0).

Verified end-to-end against `Sum[1/(i(i^2+1))]`: master identity gives
`0.671866…`, partial sum to 3000 gives `0.6718659…`. ✓

## Prerequisite audit (done 2026-06-19)

- ✅ `PolyGamma[n,z]` symbolic + numeric — `src/special_functions/polygamma.c`
- ✅ `Zeta[s]`, `Zeta[s,a]` (Hurwitz), even-`s` → π powers — `zeta.c`
- ✅ `EulerGamma`, `BernoulliB` present
- ✅ `Apart[r, Extension -> {roots}, i]` returns **linear** poles `1/(i−ρ)` for
  numeric/radical/complex roots (probed: `Extension -> I`, `-> {I, Sqrt[5]}`)
- ✅ `HarmonicNumber[n]`, `HarmonicNumber[n, r]` — implemented 2026-06-19
  (`src/special_functions/harmonicnumber.c`). Phase E (finite rational sums) is
  now unblocked.
- ⚠️ `Apart[…, Extension -> {Sqrt[b]}, i]` with **symbolic** coefficients falls
  through — the fully-symbolic case (Out[193]) is unreachable via `Apart`;
  needs the direct-residue route (Phase D).

## Naming / wiring (locked)

- Stage: **`Sum`Rational`**, `Method -> "Rational"`.
- File: `src/sum/sum_rational.c`, init `sum_rational_init()`.
- Dispatch: insert in `dispatch_def` (`sum.c`) **before** `Sum`Hypergeometric`
  (rational `p/q` should be claimed by the digamma path, not handed to PFQ).
- Scope: **infinite only** (`imax == Infinity`, `imin` a concrete integer).
  Finite/indefinite rational sums are a separate path (Phase E, deferred).

---

## Plan

### Phase A — scaffold + recognizer  (rational, rational poles)  ✅ DONE 2026-06-19
- [x] Create `src/sum/sum_rational.c`; register `Sum`Rational` builtin in
      `sum_rational_init()`; add `void sum_rational_init(void); sum_rational_init();`
      to `sum_init()` in `sum.c`; add the file to `tests/CMakeLists.txt`
      COMMON_SRC (memory: a new `src/*.c` builtin must be listed or `*_tests`
      fail to link).
- [x] Add `"Rational"` to `parse_method_option` + `SUM_METHOD_RATIONAL` enum,
      and one `try_def("Sum`Rational", …)` line in the Automatic cascade
      **before** `Sum`Hypergeometric`.
- [x] Use `sum_stage_args` to unpack `[f, var, imin, imax]`. Bail (`NULL`) unless
      `definite`, `imax` is Infinity (`is_infinity_sym`), and `imin` is
      `EXPR_INTEGER`.
- [x] Recognize rational summand: `Together` then split `Numerator`/
      `Denominator`; require both polynomial in `var`; reject if
      `deg q < deg p + 2` (divergent — leave held).
- [x] **Rational-pole sub-case first** (no extension): `Apart[f, var]`. For each
      term `c/(i−ρ)^k` with `ρ ∈ Q`, apply the master identity. Integer/rational
      `ρ` with `k≥2` → `Zeta[k, imin−ρ]`; this already reduces `1/i^s` →
      `Zeta[s]` and even-`s` → π powers via existing `zeta.c`. Terms with a
      non-linear (irreducible quadratic) base bail `NULL` → left for Phase B.
- [x] **Acceptance:** `1/i^2 → Pi^2/6`, `1/i^3 → Zeta[3]`, `1/i^4 → Pi^4/90`.
      Also verified: `1/(i(i+2)) → 3/4` (multi-term k=1 telescope),
      shifted bounds (`{i,2,Inf}`), divergence gate (`i`, `1/i` held),
      complex pole held (`1/(i^2+1)`). 52/52 existing sum tests still pass.

### Phase B — splitting-field decomposition (complex / radical roots)  ✅ DONE 2026-06-19
- [x] Compute the pole set: `Solve[Denominator==0, var]`, collect the distinct
      algebraic generators (`sr_walk_gens`/`sr_generators`: `I` + each surd
      `Power[c,p/q]`; bails on `Root[]`).
- [x] Call `Apart[f, var, Extension -> {generators}]` → linear poles. `sr_scan`/
      `sr_find_pole` descends nested `Times` (extension output isn't flattened),
      `sr_residue` recovers each residue via `Cancel[term*base^k]`.
- [x] Emit `Σ` of master-identity terms; post-pass `Simplify` (not `Together` —
      it cancels the spurious overall `I` and groups EulerGamma).
- [x] **Acceptance:** `1/(i(i^2+1))` → `1/2(2 EulerGamma + PolyGamma[0,1-I] +
      PolyGamma[0,1+I])`; `1/((i^2+3i+1)(i^2+1))` → four-`PolyGamma` form. Both
      numerically verified vs partial sums (0.12777931731762 vs 0.1277793173).

### Phase C — conjugate-pair → `Coth`/`Csc`  ✅ DONE 2026-06-19
Implemented at the over-Q residue level (cleaner than an output post-pass): an
irreducible quadratic `(A i+B)/(i^2+p i+q)` with disc<0 splits about `i=α=-p/2`
into an **odd** part (`A(i-α)` → conjugate digamma sum, kept as PolyGamma) and an
**even** part (constant → `Coth`). `sr_quadratic_contribution` (k=1):
`-(Ã/2)(ψ(t-βI)+ψ(t+βI)+2γ) + S((πβ Coth[πβ]-1)/(2β²) - Σ_{j=1}^{t-1}1/(j²+β²))`,
`t=imin-α` (gated to a positive integer; disc≥0 / higher mult / non-integer t →
extension fallback). This also yields Phase B's `1/(i(i^2+1))` over Q without an
extension; the extension path remains the fallback for disc≥0 radical factors.
- [x] **Acceptance:** `1/(i^2 (i^2+1))` → `1/6 (3 + Pi^2 - 3 Pi Coth[Pi])`.
      Also: `1/(i^2+1)` → `1/2(-1+Pi Coth[Pi])`; β=2 `1/(i^2(i^2+4))` →
      `Coth[2Pi]`; shifted t=2 `1/(i^2+1),{i,2,∞}` → `1/2(-2+Pi Coth[Pi])`
      (finite correction). All numerically verified.

### Phase D — symbolic-coefficient roots (stretch)
- [ ] `Apart` can't decompose with symbolic `a,b` (probed). Use the **direct
      residue** route instead: roots `ρ_j` from `Solve[q==0, var]` (radical /
      `Root`), residues at simple poles `c = p(ρ)/q'(ρ)` (via `D` + `sum_subst`),
      higher order via derivative formula. Emit master-identity terms directly,
      skipping `Apart`.
- [ ] **Acceptance:** `1/((i^2+3i+2)(a i^2+b))` → WL Out[193] shape (EulerGamma
      terms at `i=-1,-2` + PolyGamma at `±Sqrt[b/a] I`). Numeric spot-check at
      sample `a,b`.

### Phase E — finite / indefinite rational (deferred, separate work)
- [ ] Out of scope here. Needs `HarmonicNumber[n]`/`HarmonicNumber[n,r]` (and
      `PolyGamma` differences). Track as its own task; do **not** route finite
      rational sums into `Sum`Rational`.

### Phase D — symbolic-coefficient roots (stretch)  — NOT DONE (deferred)
Out of scope for this pass. `Sum[1/((i^2+3i+2)(a i^2+b))]` stays held: needs the
direct-residue route (`Solve[q==0]` radical roots, `c = p(ρ)/q'(ρ)` via `D`),
since `Apart` can't decompose with symbolic `a,b`. Tracked as a follow-up.

### Phase F — docs, attributes, changelog, valgrind  ✅ DONE 2026-06-19
- [x] Docstring for `Sum`Rational` via `symtab_set_docstring` (terse). Stage
      symbol `Protected | ReadProtected` like its siblings.
- [x] `docs/spec/builtins/calculus.md` Sum`Rational` section + Method list +
      this week's `docs/spec/changelog/2026-06-15.md` feature entry.
- [x] `SUM_DEVEL_PLAN.md`: Stage 5 infinite path marked shipped (finite/Phase E
      and symbolic/Phase D noted as still open).
- [x] `tests/test_sum_rational.c` (17 checks, Phases A–C) wired into
      `tests/CMakeLists.txt`; `sum_rational_tests` + `sum_tests` (52) green.
- [x] `valgrind` clean: 0 frames in `sum_rational.c`; total 12,992B/404 blocks
      vs the 12,800B/400 macOS dyld baseline (+192B/4 in the pre-existing
      Solve/Apart-extension path). Fixed one real leak found en route (the
      `gens` array in `sr_generators` — `expr_new_function` copies, doesn't adopt).

## Risks / notes
- Memory contract: `Sum`Rational` takes ownership of `res`, must not free it;
  free every owned `Expr*` on all paths (mirror `sum_hypergeometric.c`).
- Watch the `evaluate(X)`-doesn't-free-X leak class (use `eval_and_free`-style
  discipline / NULL-out-before-free).
- `Apart` Extension generator discovery is the fiddly part of Phase B — if
  `Solve` returns `Root[…]` objects, may need to extract radical generators or
  fall back to Phase D's direct-residue path uniformly.
- Convergence gate must run *before* any expensive `Solve`/`Apart`.

## Review (2026-06-19)

Shipped Phases A–C of `Sum`Rational` (`src/sum/sum_rational.c`). Key design
decisions and deviations from the original plan:

- **Phase C done at the residue level, not as an output post-pass.** The plan
  proposed detecting conjugate `PolyGamma` pairs in the result and rewriting via
  the reflection identity. That turned out fragile (the `I→-I` symmetric/
  antisymmetric split gives `anti=0` because the real result is already manifestly
  conjugation-symmetric). The clean route: handle an irreducible quadratic
  `(A i+B)/(i^2+p i+q)` (disc<0) directly from the **over-Q** `Apart` — its
  constant-numerator (even) part collapses to `Coth`, its linear-numerator (odd)
  part is the conjugate digamma sum. This naturally produces WL's exact forms and
  also resolves Phase B's `1/(i(i^2+1))` over Q (no extension).
- **Extension route is the fallback** for disc≥0 irreducible factors (real
  radical roots, e.g. the four-pole's `i^2+3i+1`). Generators (`I` + surds) are
  auto-discovered from `Solve[q==0]`; `Apart[f, i, Extension -> {...}]` fully
  splits to linear poles.
- **Post-pass is `Simplify`, not `Together`** — `Together` left a spurious overall
  `I` and an ugly complex denominator; `Simplify` rationalises the residue
  coefficients into WL's clean form.
- Generalises beyond the acceptance set: β≠1 (`Coth[2 Pi]`), shifted lower bounds
  (`t = imin-α` with the finite correction `Σ_{j=1}^{t-1} 1/(j²+β²)`).
- **Gotcha:** `Apart`'s extension output is not flattened — a pole term can be a
  nested `Times[c, Times[..., Power[quad,-1]]]`; `sr_find_pole` recurses.
- **Gotcha:** `Apart[expr, var, Extension -> α]` needs `var` *before* the option;
  with the option first it returns unevaluated. (The prereq audit's
  `Extension -> I` probe was mis-ordered.)
- **Memory:** `expr_new_function` *copies* the passed arg array; the caller must
  free it. The `gens` array in `sr_generators` leaked until fixed.

Deferred: Phase D (symbolic coefficients, direct-residue route) and Phase E
(finite/indefinite rational → HarmonicNumber).
