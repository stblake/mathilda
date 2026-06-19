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

### Phase A — scaffold + recognizer  (rational, rational poles)
- [ ] Create `src/sum/sum_rational.c`; register `Sum`Rational` builtin in
      `sum_rational_init()`; add `void sum_rational_init(void); sum_rational_init();`
      to `sum_init()` in `sum.c`; add the file to `tests/CMakeLists.txt`
      COMMON_SRC (memory: a new `src/*.c` builtin must be listed or `*_tests`
      fail to link).
- [ ] Add `"Rational"` to `parse_method_option` + `SUM_METHOD_RATIONAL` enum,
      and one `try_def("Sum`Rational", …)` line in the Automatic cascade
      **before** `Sum`Hypergeometric`.
- [ ] Use `sum_stage_args` to unpack `[f, var, imin, imax]`. Bail (`NULL`) unless
      `definite`, `imax` is Infinity (`is_infinity_sym`), and `imin` is
      `EXPR_INTEGER`.
- [ ] Recognize rational summand: `Together` then split `Numerator`/
      `Denominator`; require both polynomial in `var`; reject if
      `deg q < deg p + 2` (divergent — leave held).
- [ ] **Rational-pole sub-case first** (no extension): `Apart[f, var]`. For each
      term `c/(i−ρ)^k` with `ρ ∈ Q`, apply the master identity. Integer/rational
      `ρ` with `k≥2` → `Zeta[k, imin−ρ]`; this already reduces `1/i^s` →
      `Zeta[s]` and even-`s` → π powers via existing `zeta.c`.
- [ ] **Acceptance:** `1/i^2 → Pi^2/6`, `1/i^3 → Zeta[3]`, `1/i^4 → Pi^4/90`.

### Phase B — splitting-field decomposition (complex / radical roots)
- [ ] Compute the pole set: `Solve[Denominator==0, var]`, collect the distinct
      algebraic generators (e.g. `I`, `Sqrt[5]`) for the `Extension` list.
- [ ] Call `Apart[f, Extension -> {generators}, var]` → linear poles. Parse each
      `c_{j,k}/(var − ρ_j)^k` term: read multiplicity `k` and residue `c`.
- [ ] Emit `Σ` of master-identity terms; wrap in `Together`/`Simplify` for a
      clean form. Sum simple-pole residues to confirm `Σc = 0` (sanity gate).
- [ ] **Acceptance:** `1/(i(i^2+1))` → `1/2(2 EulerGamma + PolyGamma[0,1-I] +
      PolyGamma[0,1+I])`; `1/((i^2+3i+1)(i^2+1))` → four-`PolyGamma` form.
      Verify each numerically vs a 3000-term partial sum *inside the language*
      (memory: don't parse printed reals).

### Phase C — conjugate-pair → `Coth`/`Csc` post-pass
- [ ] Detect conjugate pole pairs `ρ = α ± βI` with equal multiplicity/residue
      magnitude and rewrite `PolyGamma[k, a+βI] + PolyGamma[k, a−βI]` via the
      reflection identity `Σ_{n≥1} 1/(n^2+β^2) = (πβ coth(πβ) − 1)/(2β^2)`
      (and the `k=1` analogue → `Csc^2`).
- [ ] Apply only when it strictly simplifies (WL keeps the raw PolyGamma form in
      Out[186] but collapses to `Coth` in Out[191]) — gate on "no leftover
      explicit `I`".
- [ ] **Acceptance:** `1/(i^2 (i^2+1))` → `1/6 (3 + Pi^2 - 3 Pi Coth[Pi])`.

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

### Phase F — docs, attributes, changelog, valgrind
- [ ] Docstring for `Sum`Rational` via `symtab_set_docstring` (terse; no
      examples — memory). Stage symbol Protected like its siblings.
- [ ] Update `docs/spec/builtins/` Sum entry + this week's
      `docs/spec/changelog/2026-06-15.md` (Monday of the ISO week).
- [ ] Update `SUM_DEVEL_PLAN.md`: mark Stage 5 shipped.
- [ ] `tests/sum_rational_tests.c` covering the acceptance set (Phases A–C
      mandatory, D if landed). Run only that binary (memory: scope tests).
- [ ] `valgrind` clean vs the macOS dyld baseline; build foreground.

## Risks / notes
- Memory contract: `Sum`Rational` takes ownership of `res`, must not free it;
  free every owned `Expr*` on all paths (mirror `sum_hypergeometric.c`).
- Watch the `evaluate(X)`-doesn't-free-X leak class (use `eval_and_free`-style
  discipline / NULL-out-before-free).
- `Apart` Extension generator discovery is the fiddly part of Phase B — if
  `Solve` returns `Root[…]` objects, may need to extract radical generators or
  fall back to Phase D's direct-residue path uniformly.
- Convergence gate must run *before* any expensive `Solve`/`Apart`.

## Review
_(to be filled after implementation)_
