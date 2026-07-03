# Plan: fast, rigorous polynomial arithmetic over algebraic-extension coefficients (FLINT-backed)

**Status:** in progress (2026-06-30). Revised to assume FLINT as the
arithmetic engine. **M0** done (FLINT 3.6.0 installed; toolchain probed — both
`fmpq_mpoly_gcd` and `gr_ctx_init_nf`+`gr_poly_gcd` over `Q(√2)` verified).
**M1 build integration** done (`USE_FLINT` autodetect + graceful degrade).
**M1 bridge** foundation done: `src/poly/flint_bridge.{c,h}` with the rational
multivariate case `Q[x_1..x_n]` (`Expr ⇄ fmpq_mpoly`), exercised via the
scaffolding builtin `` Flint`GCD `` and valgrind-clean. **M2 (number-field GCD)**
first increment done: `flint_numberfield_gcd` — rigorous univariate GCD over
`Q(√d)` (single radical) via `gr_ctx_init_nf`+`gr_poly_gcd`, e.g.
`gcd(x²−2, x−√2)=x−√2`; tested + valgrind-clean. Generalised to **any simple
algebraic generator** (`absfield_gcd_core` + `GenSpec`, reduce mod minpoly via
`fmpq_poly_rem`) and added **Q(ζ_n)** (`flint_cyclotomic_gcd`, `Φ_N` via
`fmpz_poly_cyclotomic`); `` Flint`GCD `` dispatches rational → Q(√d) → Q(ζ_n).
The **tower collapse** is done (`flint_tower_gcd` — Q(√d₁,…,√d_r) via
linear-algebra primitive-element collapse + product-basis readback).

✅ **Wired into real consumers (M5, GCD + division).** `PolynomialGCD`
(`builtin_polynomialgcd`) and `Cancel`/`Together` (`flint_cancel_fraction` +
`cancel_exact_div_wrapper` fallback) now route the extension cases through the
bridge: `flint_extension_gcd` and a new `flint_extension_divexact` (exact division
via `gr_poly_divrem`, shared `absfield_op_core`). `Cancel[(x²−2)/(x−√2)]` now
gives `x+√2`; Q(ζ_n) and radical-tower fractions reduce too; plain Q[x] is
unchanged. Two stale test expectations updated to the canonical FLINT form
(`#ifdef USE_FLINT`-guarded). Regression suites green.

✅ **M3 (parametric radical `Q(t_1..t_p)(√k)`, k a symbol) done** (2026-07-01).
Key simplification vs the original design: because `√k` is *transcendental* over
`Q(t)` with `k = (√k)²`, the field `Q(a,b,k)(√k)` **is** the rational function
field `Q(a,b,√k)` — there is no genuine algebraic relation to reduce. So this
whole `p ≥ 1, r = 1, deg = 2` regime collapses to ordinary multivariate GCD /
exact division over `Q` after the substitution `√k → S, k → S²`
(`flint_parametric_sqrt_gcd` / `flint_multivariate_divexact`), then readback
`S → √k`. **No `fq_nmod` residue field, CRT, or outer loop are needed here** —
those remain only for genuine constant-coefficient minimal polynomials that
define a `GF(p^d)` (a later milestone). Wired into `flint_extension_gcd` /
`flint_extension_divexact`, so `PolynomialGCD` / `Cancel` / `Together` reduce over
`Q(a,b,k)(√k)` sub-second (`Cancel[(x²−k)/(x−√k)] = √k+x`; the `(x−√k)`-cofactor
Cancel verified via `PossibleZeroQ`). Fixed a `flint_cancel_fraction` sign-unit
artifact (numeric divided-out denominator now distributed). Tests: `test_parametric`
in `test_flint_bridge.c`; regressions green; valgrind-clean.

✅ **Genuine-algebraic parametric tower — zero-reduction increment done**
(2026-07-04). `flint_algebraic_field_normalize` (`src/poly/flint_bridge.c`) handles
the higher-degree / polynomial-radicand parametric regime `Q(params)(α_1..α_r)`
(cube roots, polynomial/symbol radicands, roots of unity — e.g. the Goursat tower
`Q(x,k)(∛(x(1-x)(1-kx)), ∛k, (-1)^(1/3))`) by modelling `K = Q[params,gens]/I`
with the monic minimal polynomials as a Gröbner basis (generators = leading LEX
vars → pairwise-coprime leading monomials) and reducing via
`fmpz_mpoly_divrem_ideal`. Wired into `flint_cancel_fraction` (Cancel/Together,
hence Simplify): it returns `0` exactly when the input is identically zero in `K`
(numerator reduces to 0 mod `I`, denominator not) — a **rigorous** field zero test
with **no numeric oracle** and no primes/CRT, else `NULL` (classical fallback). So
`Simplify[D[Integrate[f],x] - f] = 0` for the cube-root Goursat family, while
`Sqrt[x^2]-x`, `(...)^(1/3)+k^(1/3)`, and perturbed variants stay put. This is a
zero-reduction / Together increment, not full Cancel-to-lowest-terms.

**Remaining:** full Cancel over the genuine-algebraic tower (rationalise the
denominator: field inversion via tower conjugates / xgcd, or the `fq_nmod` outer
loop with CRT — needed for the non-zero reduced form, not the zero test);
number-field **factoring** over a tower (FLINT's `gr` layer cannot —
`gr_poly_roots` → `GR_UNABLE`; Trager stays in `qafactor.c`); then re-measure the
full Goursat `Cancel`/`Integrate` end-to-end (M7 acceptance).
**Motivation:** the two classic Goursat named integrals (and a broad class of
algebraic integrands, `Together`/`Cancel`/`Apart`/`Factor` over splitting
fields) are blocked by a single root cause — Mathilda has no efficient *and*
rigorous polynomial arithmetic over coefficient rings that contain algebraic
generators (`Sqrt[k]`, `ζ_n`, …). This document specifies the architecture to
fix that by **building on FLINT** (≥ 3.0) rather than hand-rolling a modular
multivariate-GCD / factoring engine, replacing the current pseudo-variable
workarounds.

---

## 1. Problem statement

A great deal of CAS work happens in a ring

```
R[x_1, …, x_n]      where   R = Q(t_1, …, t_p)(α_1, …, α_r)
```

i.e. polynomials in *variables* `x_i` whose *coefficients* live in a field built
from rational functions in transcendental **parameters** `t_j` (free symbols
like `a, b, k`) extended by **algebraic generators** `α_l` (each with a minimal
polynomial `m_l` over the field below it: `Sqrt[k]` ↦ `y² − k`, `ζ_n` ↦ `Φ_n`,
nested radicals ↦ towers).

The core operations — `+ − ×`, exact division, **GCD**, **extended GCD**,
content/primitive part, square-free decomposition, **factorisation** — must be
both:

* **rigorous**: compute in the *actual* ring `R[x]`, using the relations
  `m_l(α_l) = 0`. A GCD over `R` may be strictly larger than the GCD computed
  with `α_l` treated as a free indeterminate, because `α_l`'s relation enables
  cancellations the free-variable view misses (e.g. `y² − k` is irreducible in
  `Q[k, y]` but splits as `(y − √k)(y + √k)` in `Q(k)(√k)`); and
* **efficient**: polynomial in input size and parameter count, *not* exponential
  in the number of parameters/generators.

### 1.1 Why today's code fails

Mathilda currently has two partial mechanisms, each insufficient:

1. **The QA substrate** (`src/poly/qa.{c,h}`, `qaupoly.{c,h}`,
   `qafactor.{c,h}`). `QAExt` defines `Q(α)` via a minimal polynomial with
   **`mpq_t` (rational) coefficients** — so `α` must be algebraic over **Q**, not
   over a function field. `QAUPoly` is **univariate** (a single `x`) with `QANum`
   coefficients. There is no room for parametric (`a,b,k`-dependent) coefficients
   and no second polynomial variable. `qaupoly_gcd` is a classical Euclidean GCD
   over `Q(α)` — correct but with the usual coefficient blow-up. Net: the
   substrate cannot even *represent* `Q(a,b,k)(√k)[x]`.

2. **The pseudo-variable path** (`extension_autodetect`,
   `qa_cancel_with_poly_radical`, the recursive `poly_gcd_internal`). When the QA
   substrate doesn't apply, the code falls back to treating the algebraic
   generator as an ordinary polynomial variable `S` and doing multivariate
   arithmetic in `Q[t_1,…,t_p, S]`. This is **non-rigorous** (it computes in the
   polynomial ring, not the quotient by `m(S)`, so it *under-reduces*) and
   **non-performant** (`poly_gcd_internal`'s recursive content / pseudo-remainder
   algorithm is exponential in `p`; this is the measured hang in
   `CANCEL_IMPROVEMENT_PLAN.md`).

### 1.2 Evidence (this session)

Profiling the Goursat Example-1 descent surfaced a *chain* of bottlenecks, every
one of them the same missing capability:

| Site | Operation | Symptom |
|------|-----------|---------|
| `goursat_v4`→`canonic`→`Cancel` Phase E (`qafactor.c:4168`) | `PolynomialExtendedGCD` over `Q[a,b,k,S]` | hang (deep `together_recursive`+GCD) |
| `together_recursive`→`cancel_recursive` (`rat.c:605`) | multivariate `PolynomialGCD` (`poly_gcd_internal`) | exponential coefficient growth |
| reduced `∫G/√Q` → `IntegrateRational` (`intrat.c:3431`) | `PolynomialQuotientRemainder`→`exact_poly_div` over `Q[√k,u]` | hang (multivariate exact division) |

Interim prototypes (a heuristic evaluation GCD "GCDHEU", and a q=2 conjugate
rationalisation in Phase E) were implemented and **reverted**: GCDHEU treats
`√k` as a free variable (under-reduces, and its non-canonical multivariate sign
broke `FactorSquareFree`/`Factor` structural tests), and gating it on the leaky
`budget_blown` heuristic merely displaced the failure. The lesson — confirmed
with the user — is that band-aids on the pseudo-variable path cannot be made
rigorous. The fix must be the real ring arithmetic below.

### 1.3 Decision: buy the engine (FLINT), build only the bridge

Reimplementing rigorous *and* competitively fast modular multivariate GCD,
number-field arithmetic, and Trager factoring from scratch is "the largest CAS
subsystem added to date" — and matching Magma/Pari/Maple on constant factors is
20+ years of tuning that a hand-rolled C99 engine will not reach. **FLINT** is
the mature library that already implements exactly this machinery, tuned to that
level, and its hard dependencies (**GMP**, **MPFR**) are *already linked* by
Mathilda. The pico-CAS "hand-rolled" preferences for the surface language do not
require hand-rolling the number-theory core any more than they require
hand-rolling GMP; ECM and LAPACK are existing precedents for vendoring/linking a
heavyweight numeric dependency behind an autodetect-and-degrade switch.

Therefore this plan **uses FLINT as the engine** and limits Mathilda's own new
code to:

* the **`Expr ⇄ FLINT` bridge** (recognise membership in `R[x]`, convert in and
  out) — the genuinely Mathilda-specific work; and
* a **thin outer loop** (prime selection + CRT + rational reconstruction +
  trial-division verify) for the one regime FLINT does not expose as a single
  call: GCD/division over a *parametric* algebraic extension
  `Q(t_1,…,t_p)(α)[x]`. Even there, FLINT supplies the expensive inner
  arithmetic (multivariate GCD over a finite field, with its own
  Brown/Zippel parameter interpolation), so this loop is small.

What FLINT removes from the original plan: the word-size `Fp`/`Fp[α]` inner
arithmetic, packed-exponent sparse representation, Brown/Zippel parameter
interpolation, multivariate Hensel lifting (Wang/EEZ) for the factor norm, and
all of their tuning and valgrind exposure. These were the bulk of the risk.

---

## 2. Mathematical foundation (implemented by FLINT)

The efficient, rigorous route is **modular / evaluation–interpolation
arithmetic** over algebraic extensions. FLINT implements every piece below; this
section records *which FLINT capability* covers each so the mapping in §3.4 is
auditable, not so we re-implement them.

* **Multivariate GCD over Q/Z** (Brown dense, Zippel sparse, Hensel): FLINT
  `fmpq_mpoly_gcd` / `fmpz_mpoly_gcd`. Replaces `poly_gcd_internal`.
* **Multivariate factoring over Q/Z** (Wang/EEZ multivariate Hensel, the hard
  part the old plan glossed): FLINT `fmpq_mpoly_factor` / `fmpz_mpoly_factor`.
* **Number-field arithmetic** `Q(α)`: FLINT's bundled **ANTIC** (`nf_t`,
  `nf_elem_t`), with the minimal polynomial as an `fmpq_poly`. Number-field
  *polynomial* GCD/factoring is reached either through ANTIC helpers or FLINT
  3's generic-rings layer (`gr_ctx_init_nf` → `gr_poly`/`gr_mpoly`); the exact
  entry point is pinned in M1 (see §6) and isolated behind the bridge so the
  choice is invisible to the rest of Mathilda.
* **Finite-field multivariate polynomials** `GF(p^d)[params][x]`: FLINT
  `fq_nmod_mpoly` + `fq_nmod_mpoly_gcd`. This is the workhorse for the
  *parametric* extension regime (§3.4): it already performs Brown/Zippel
  interpolation in the parameters `t_j` internally, so our outer loop only
  handles the lift from `GF(p^d)` back to `Q(α)`.
* **Rational reconstruction** and **CRT**: FLINT `fmpq_reconstruct_fmpz`,
  `fmpz_CRT`. The standard Encarnación recovery of a `Q(α)` result from its
  images mod good primes.
* **Termination / bounds**: FLINT's GCD/factor routines are internally bounded
  and self-verifying; for the hand-written parametric outer loop we add an
  explicit prime budget from Landau–Mignotte (degree in each `t_j`, integer
  height, `deg α`) and a final trial-division check, so it is finite and
  provably complete, not heuristic.

Every result the outer loop produces is **trial-division verified** against the
inputs (FLINT's own routines are already rigorous), so correctness never depends
on the probabilistic choices (lucky primes); those only affect *speed*.

---

## 3. Architecture

### 3.1 Dependency & build integration

* **Optional, autodetected, graceful-degrade** — mirrors `USE_GRAPHICS` /
  LAPACK. A `USE_FLINT` switch; when FLINT is absent (or too old), the build
  proceeds with the existing classical path and a `$(warning …)`.

  ```make
  USE_FLINT ?= 1
  ifeq ($(USE_FLINT), 1)
    ifneq ($(shell $(PKG_CONFIG) --exists 'flint >= 3.0' && echo y),)
      CFLAGS  += -DUSE_FLINT $(shell $(PKG_CONFIG) --cflags flint)
      LDFLAGS += $(shell $(PKG_CONFIG) --libs flint)
    else
      $(warning FLINT >= 3.0 not detected; USE_FLINT=0 (algebraic-extension \
                GCD/Factor use the classical fallback))
      override USE_FLINT := 0
    endif
  endif
  ```

* **Version floor: FLINT ≥ 3.0** — ANTIC (number fields) was merged into FLINT
  at 3.0. Older `apt` packages (2.x) lack it; the autodetect checks the version,
  not mere presence.
* **C99 isolation** — all FLINT calls live in one translation unit
  (`src/poly/flint_bridge.{c,h}`). FLINT headers may warn under
  `-std=c99 -Wall -Wextra`; if so, relax flags for that single object, and never
  let FLINT headers leak into the rest of the tree.
* **License** — FLINT is LGPL-3, compatible with Mathilda's GPLv3. System
  dependency (like readline), not a `src/external/` vendored tree; **no
  `src/external/` edits**.
* **The fallback must stay rigorous** — with `USE_FLINT=0`, the classical path
  must remain *correct* (even if slow); only the fast path is gated. "This is now
  fast" claims are conditional on FLINT being present, and tests assert that.

### 3.2 Coefficient-ring descriptor — `AExtField`

A descriptor for `R = Q(t_1,…,t_p)(α_1,…,α_r)` that the bridge builds from
`extension_autodetect`, and which maps to the appropriate FLINT context:

```
typedef struct AExtField {
    int      n_params;          /* transcendentals t_j (free symbols)        */
    char**   params;            /* interned names, canonical order           */
    int      n_gens;            /* algebraic generators α_l, tower order      */
    AMinPoly* minpoly;          /* α_l's minimal poly over Q(t)(α_1..α_{l-1}) */
} AExtField;
```

`AMinPoly` carries each `m_l` as a **multivariate** polynomial in
`(t_1..t_p, α_1..α_l)` (so `y² − k`, `Φ_n(ζ)`, nested-radical relations are all
expressible). Dispatch on `(p, r)`:

| Regime | `(p, r)` | FLINT mapping |
|--------|----------|---------------|
| rational functions only | `r = 0` | `fmpq_mpoly` over `{t_j, x_i}` |
| pure number field | `p = 0, r = 1` | ANTIC `nf_t` from `m`'s `fmpq_poly` |
| number-field tower | `p = 0, r ≥ 2` | collapse to absolute primitive element θ → single `nf` (§3.2.1) |
| **parametric extension** | `p ≥ 1, r ≥ 1` | outer loop over `fq_nmod_mpoly` (§3.4) |

#### 3.2.1 Towers collapse to a single absolute field

FLINT/ANTIC has **no native nested/relative extension type**: `nf_init` and
`gr_ctx_init_nf` take a single `fmpq_poly` minimal polynomial, i.e. an
*absolute* field `Q(θ)` over `Q`. (`gr_ctx_init_gr_poly` can build
`Q(α_1)[y]`, but there is no turnkey constructor to quotient it by a relative
minimal polynomial into the field `Q(α_1)(α_2)`.) A tower
`Q(α_1,…,α_r)` — including nested radicals like `Sqrt[1+Sqrt[2]]` — is therefore
handled the standard CAS way: **collapse it to a single primitive element** `θ`
(primitive element theorem; always exists in char 0), compute `θ`'s absolute
minimal polynomial via **resultants** (`fmpq_poly_resultant`, which FLINT
provides — `θ = α_1 + c·α_2`, `m_θ = Res_y(m_1(y), m_2(θ − c·y))`, retry on small
integer `c` until square-free of the right degree), and record each `α_l` as a
polynomial in `θ` for round-tripping. The collapse is the bridge's job; FLINT
supplies the resultant and the resulting absolute-field arithmetic. Caveat: the
absolute degree is the **product** of the relative degrees (`Q(√2,√3,√5)` → 8),
so deep towers grow `deg θ` multiplicatively — fine for the Goursat targets
(single `√k`, shallow nested radicals), costly for tall towers. For the
parametric case the collapse happens *after* evaluating the `t_j` at integers
(§3.4), so each prime/point sees an absolute number field.

`QAExt`/`QANum` remain the `p = 0, r = 1` legacy representation, but their
*operations* are re-pointed at ANTIC (M2), so the slow `qaupoly_gcd` Euclid goes
away while the surface type can stay until fully retired.

### 3.3 The bridge — `Expr ⇄ FLINT`

The one substantial piece of new Mathilda code (analogous to today's
`gb_from_expr` / `rebuild_from_bp`, extension-aware):

* **In**: recognise an `Expr` as a member of `R[x]` (via `extension_autodetect`,
  extended to record **symbolic-radicand** generators `√k` as parametric
  `AMinPoly`s instead of rejecting them), then emit the matching FLINT object
  (`fmpq_mpoly`, `nf_elem` polynomial, or the per-prime `fq_nmod_mpoly`).
* **Out**: render a FLINT result back to a canonical `Expr`. Honour
  `feedback_together_no_expand`: carry `x`-degree without forcing expansion of
  parameter sub-expressions; do not expand `Power[Plus[...], n]`.
* **Canonical normalisation**: unit-normalise the FLINT output deterministically
  (FLINT's monic/`fmpq` normalisation gives this for free), so the surface form
  is run-independent and a true drop-in for `Factor`/`FactorSquareFree`
  consumers — removing the sign non-determinism that sank the reverted GCDHEU.

### 3.4 The engine — operation → FLINT

```
amod_gcd(A, B, R)        amod_divexact(A, B)      amod_extgcd(A, B)
amod_content / amod_primpart / amod_squarefree    amod_factor(A, R)
```

are thin wrappers in `flint_bridge.c`:

* **`r = 0` / `p = 0` cases** are *direct* FLINT calls: `fmpq_mpoly_gcd`,
  `fmpq_mpoly_divides`, `fmpq_mpoly_factor`; ANTIC for `Q(α)` univariate
  GCD/factor. No outer loop.
* **Parametric case `p ≥ 1, r ≥ 1`** (the actual Goursat ring, e.g.
  `Q(a,b,k)(√k)[x]`) — Encarnación-style outer loop, but with FLINT doing the
  expensive interior:
  1. Pick a prime `p` for which `m mod p` stays **irreducible** ⇒ residue field
     `GF(p^{deg m})`, a true field (one prime ideal of full degree over `p`).
     Reject primes where `m mod p` factors or `lc` vanishes (*unlucky* —
     correctness-safe to skip).
  2. Build `fq_nmod` = `GF(p^{deg m})`; map the parameters `t_j` to the
     variables of an `fq_nmod_mpoly`. Call `fq_nmod_mpoly_gcd` — **FLINT performs
     the Brown/Zippel interpolation over the `t_j` itself**.
  3. **CRT** the per-prime images and **rational-reconstruct**
     (`fmpq_reconstruct_fmpz`) to `Q(t)(α)`.
  4. **Trial-divide** both inputs to verify; on failure add a prime (bounded by
     the Landau–Mignotte estimate) and retry.
  5. Output is canonically unit-normalised (§3.3).
  `amod_divexact` reuses steps 1–4 with exact division as the per-prime
  operation; `amod_factor` for the parametric case uses **Trager** (norm via
  `amod` resultant → `fmpq_mpoly_factor` over `Q(t)` → GCD back up), reusing the
  same loop.

The only genuinely new arithmetic code is steps 1, 3, 4 (prime selection, CRT +
reconstruction, verify) — all FLINT-assisted and a few hundred lines, versus the
multi-kLoC modular engine the pre-FLINT plan required.

---

## 4. Integration points (what gets re-pointed)

Replace the pseudo-variable / classical paths, one consumer at a time, behind
the same public builtins, all routed through `flint_bridge.c`:

* `PolynomialGCD` / `poly_gcd_internal` (`poly.c`): dispatch to `amod_gcd` when
  the inputs lie in `R[x]` with `r ≥ 1` *or* `p ≥ 2` parameters (the regime
  where pseudo-remainder blows up). Keep the current code for trivial
  univariate-over-Q; the size threshold is **measured in M0**, not guessed.
* `Cancel` / `Together` Phase E (`qafactor.c` `qa_cancel_with_poly_radical`,
  `rat.c` `together_recursive`): combine + reduce via `amod_*` over the detected
  `R`, deleting the substitute-`S`/`PolynomialExtendedGCD`/multivariate-`Together`
  pipeline that hangs today.
* `Apart` (`parfrac.c`) and `IntegrateRational` (`intrat.c`): their
  `PolynomialQuotientRemainder` / `exact_poly_div` over `Q(√k)` route through
  `amod_divexact` — fixing the third bottleneck (`∫G/√(quadratic over Q(√k))`).
* `Factor` (`mvfactor*.c`): `amod_factor` for inputs over `R` (direct FLINT for
  `Q`/`Q(α)`; Trager-over-FLINT for the parametric case).
* `extension_autodetect`: keep as the *detector* that produces an `AExtField`
  (it already recognises integer-base radicals, roots of unity, nested radicals);
  extend it to record **symbolic-radicand** generators (`√k`) as parametric
  `AMinPoly`s instead of rejecting them, and hand the field to the bridge.

---

## 5. Rigor & correctness obligations

* **Soundness**: FLINT's GCD/factor routines are rigorous; the parametric outer
  loop additionally trial-division verifies every result, so a bad prime can only
  cost time, never correctness.
* **Completeness**: FLINT bounds its own loops; the outer loop's prime count is
  bounded below by Landau–Mignotte / degree bounds, so it provably finds the
  *full* GCD (no silent under-reduction — the defect of the pseudo-variable
  hack).
* **Zero-divisor handling**: choose primes where `m mod p` is irreducible (true
  residue field); detect and skip the unlucky ones. Never divide by a zero
  divisor.
* **Determinism**: canonical unit normalisation (§3.3) ⇒ identical surface form
  across runs and a true drop-in for structural consumers (`Factor`,
  `FactorSquareFree`).
* **Fallback rigor**: with `USE_FLINT=0` the classical path must stay correct
  (slow is acceptable, wrong is not). CI builds at least one `USE_FLINT=0`
  configuration.
* **No edits under `src/external/`**; C99 (FLINT confined to one TU);
  valgrind-clean — FLINT objects (`fmpq_mpoly_t`, `nf_elem_t`, `fq_nmod_*`) are
  init/clear-paired exactly like `mpz_t`, and the bridge owns every conversion
  buffer.

---

## 6. Milestones (incremental, each independently testable)

0. **M0 — profile the real operands.** Instrument the three §1.2 sites
   (`qafactor.c:4168`, `rat.c:605`, `intrat.c:3431`) to dump the actual
   GCD/division operands (degree in `x`, parameter count, term support, integer
   height, `deg α`). Run the *same* polynomials through `gp` (Pari) and, if
   available, Magma to capture reference wall-clocks. Output: concrete size
   tuples + the benchmark targets that make §6/§7's "fast" measurable, and the
   data to set the M4 dispatch threshold. *No code change.*
1. **M1 — build integration + bridge + ring sanity.** `USE_FLINT` makefile
   block; `flint_bridge.{c,h}`; `Expr ⇄ {fmpq_mpoly, nf_elem, fq_nmod_mpoly}`;
   `AExtField` from `extension_autodetect`. **Pin the ANTIC vs `gr` entry point
   for number-field polynomial GCD/factor here.** Unit tests: round-trip and
   `+ − ×`/reduce over `Q(k)(√k)`, `Q(ζ_6)`, a nested radical. *No behaviour
   change to builtins yet.*
2. **M2 — `amod_gcd`/`amod_factor`, number field (`p = 0`).** Direct FLINT/ANTIC
   for `Q(α)`, `Q(ζ_n)`. Re-point `qaupoly_gcd`'s callers. Acceptance: results
   **associate-equal** to the old QA output (verified by exact ratio-is-a-unit /
   `PossibleZeroQ[a − b]`), **in canonical form from day one** — not
   byte-identical (the new normalisation deliberately differs; the old QA corpus
   strings that encoded the Euclid normalisation are updated here, with the
   equality noted).
3. **M3 — parametric outer loop (`p ≥ 1, r ≥ 1`).** Prime selection (irreducible
   `m mod p`), `fq_nmod_mpoly_gcd` per prime, CRT + rational reconstruction,
   trial-division verify. Benchmark against `poly_gcd_internal` and the M0 Pari
   targets on the Goursat cofactors (target: the `Q[a,b,k,S]` GCD that hangs →
   within a small factor of Pari, certainly sub-second).
4. **M4 — division / extended GCD / content / squarefree.** `amod_divexact` etc.
   on the engine; set the `poly.c` dispatch threshold from M0 data. Fixes the
   `IntegrateRational` `exact_poly_div` bottleneck.
5. **M5 — `Cancel`/`Together`/`Apart` re-point.** Delete the Phase-E
   substitute-`S` pipeline. Acceptance: `Cancel[…, Extension→Automatic]` over
   `Q(a,b,k,√k)` is fast and **fully reduced** (rigorous); number-field and
   cyclotomic regressions unchanged.
6. **M6 — `Factor` re-point.** Direct FLINT for `Q`/`Q(α)`; Trager-over-FLINT
   (norm → `fmpq_mpoly_factor` → GCD back) for the parametric case. Acceptance:
   `FactorSquareFree`/`Factor` corpus passes with canonical output.
7. **M7 — Goursat acceptance.** Both classic integrals close end-to-end with the
   ArcTan / ArcTanh answers (`CANCEL_IMPROVEMENT_PLAN.md` §5), under the normal
   `TimeConstrained` budget, with `D[result,x] − integrand ≡ 0` verified.

---

## 7. Testing strategy

* New `tests/test_flint_bridge.c` (round-trip + ring ops), `test_amod_gcd.c`
  (GCD vs a brute-force reference over small fields; randomised + adversarial
  unlucky primes), `test_amod_factor.c`.
* **Both build configs in CI**: `USE_FLINT=1` (fast path) and `USE_FLINT=0`
  (fallback stays rigorous, even if slow). Tests that assert *speed* are gated on
  `USE_FLINT`.
* Regression guards: full `poly`/`mpoly`/`facpoly`/`rat`/`intrat`/`radical`
  suites stay green (canonical normalisation in M2/M6 keeps `Factor` output
  stable — the failure mode the reverted GCDHEU exposed).
* Corpus: re-run `intrat_corpus`, `crc_corpus`, and the Goursat suite.
* Performance asserts (wall-clock ceilings, derived from M0/Pari) on the
  `Q[a,b,k,√k]` GCD, the Phase-E `Cancel`, and the two integrals.
* valgrind on a representative `Cancel`/`Factor`/`Integrate` over `Q(a,b,k,√k)`
  — FLINT objects init/clear-paired.

---

## 8. Risks & scope

* **New heavyweight dependency.** FLINT (LGPL-3) becomes a build/runtime dep for
  the fast path. Mitigated by autodetect + graceful degrade + a rigorous
  `USE_FLINT=0` fallback, and by precedent (ECM, LAPACK). Distribution docs note
  the `flint >= 3.0` requirement.
* **Number-field polynomial API uncertainty.** Whether `Q(α)` polynomial
  GCD/factor is cleanest via ANTIC helpers or FLINT 3's generic `gr` rings is
  resolved in **M1** and hidden behind the bridge, so it cannot leak risk into
  consumers.
* **Parametric outer loop is the only hand-written algorithm.** Bad-prime /
  bad-point handling for `α` over `Q(t)`. Mitigated by trial-division
  verification (soundness never at risk) and bound-driven retry; the expensive
  interior is FLINT's, not ours.
* **C99 / `-Wextra` header friction** from FLINT — confined to one TU with
  relaxed flags if needed.
* **Determinism vs. legacy test strings** — a small, enumerated set of
  `Factor`/`FactorSquareFree` (and QA) expectations is updated to canonical form
  in M2/M6 (documented, with the mathematical equality shown).
* **Scope is now bridge + thin loop**, not "largest subsystem ever" — the
  modular engine, sparse representation, multivariate Hensel, and their tuning
  move into FLINT. M0–M2 are self-contained and verifiable before any builtin
  re-points (M5+).

---

## 9. Relationship to prior notes

* `CANCEL_IMPROVEMENT_PLAN.md` — the *symptom* writeup (the Cancel hang). This
  document is the *rigorous cure*: M5 deletes the Phase-E pipeline that plan
  describes.
* `project_cyclotomic_extension_support`, `project_cancel_autodetect_symbolic_radicand`
  (memory) — the cyclotomic and symbolic-radicand cases become two `AExtField`
  instances routed through the same FLINT bridge, not separate special cases.
* `feedback_together_no_expand` (memory) — preserved: the bridge must not expand
  `Power[Plus[...], n]`; it carries `x`-degree without forcing expansion of
  parameter sub-expressions.
* Earlier (pre-FLINT) revision of this plan proposed a hand-rolled modular engine
  (word-size `Fp`, packed-exponent sparse `APoly`, Brown/Zippel, multivariate
  Hensel). That design is superseded: FLINT provides all of it, tuned. The
  `AExtField` descriptor and the `Expr` bridge survive from it; the engine does
  not.
</content>
</invoke>
