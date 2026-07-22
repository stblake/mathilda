# Together / Cancel / Simplify — unified rational-normalization rewrite

**Status:** planning. Phases are scoped to be executed one-per-session from a
fresh context. Read §0 (primer) first every time, then the one phase you are
executing. Do not skip the per-phase acceptance gate.

Related memory: `[[project_ratcanon_unified_classifier]]`. Prior work
(the classifier + shadow harness now on branch
`feature/ratcanon-unified-classifier`) is the *starting point*, not the target
— see §0.3.

---

## 0. Context primer (read every session)

### 0.1 The problem

`Together`/`Cancel` (and `Simplify`, which drives them) resolve through ~8
terminal routing paths that dispatch among a **zoo of partial FLINT engines**
(`flint_cancel_fraction`, `flint_gaussian_*`, `flint_parametric_field_*`,
`flint_algebraic_field_{normalize,canonical,together}`, the QA tower path,
Phase E `qa_cancel_with_poly_radical`, and the classical cascade). Each engine
builds its **own** generator set, has its **own** output-form convention, covers
a **different subset**, and falls through to the **weak classical multivariate
GCD** (`poly_gcd_internal`) on its gaps — which is what actually hangs/SIGSEGVs.

The accumulated "fixes" are per-case hacks that reconcile the zoo's
inconsistencies: sign-normalizers (two conventions), decline-heuristics (dodge
one engine's bad form or another's hang), dependent-generator bail guards, and
arbitrary caps (degree-32, leaf-200, iter-50, `TRIGRAT_RADICAL_MAX`). This is
whack-a-mole; it will never converge.

### 0.2 The general algorithm (the target)

A `Together`/`Cancel` is a **normal-form computation in a rational-function
field over one differential/algebraic tower**. There is exactly one such normal
form and FLINT computes it directly. Build *that*, once, and route everything
through it — then delete the zoo.

1. **One tower (structure theorem — a decision, not a guess).** From `E` build
   the minimal ordered generator set:
   - constant field `K = Q(θ)` — the number field of the algebraic constants
     (`I`, `Sqrt[d]`, `(-1)^(p/q)`), one primitive element θ, minimal poly `m_θ`;
   - monomials `t_1 < … < t_n` over `K(x)`, each either **transcendental**
     (Log/Exp/Tan/inverse-trig — independence *verified* by
     `risch_rational_span`, Bronstein §9.3, after `Log[ab]→Log a+Log b` and
     reducing exp exponents to a common fundamental) or an **algebraic function**
     (`Sqrt[k]`, `k^(1/q)`) with relation `t_i^q - g_i = 0`.
   Independence is *proven*, so the generators handed to the GCD are genuinely
   independent — the weak-GCD blow-up mode cannot arise.

2. **One representation.** Rewrite `E` as a single `num/den ∈ K[x,t_1..t_n]/I`,
   where `I` is the ideal of the algebraic relations (`m_θ`, `t_i^q - g_i`),
   the algebraic generators ordered as the leading LEX variables so `I` is a
   Gröbner basis. Transcendental `t_i` are free variables.

3. **One reduction (FLINT's complete primitives; the weak GCD is never called):**
   combine to a single fraction (`fmpz_mpoly_q`); reduce `num`,`den` mod `I`
   (`fmpz_mpoly_divrem_ideal`); cancel `gcd(num,den)` (`fmpz_mpoly_gcd`, over the
   quotient ring). `Together` combines; `Cancel` maps the reduction over the
   top-level additive terms.

4. **One output convention** (§0.4), applied once at map-back.

**Why every hack disappears:** one convention ⇒ no sign-normalizer; recursive
kernel-argument normalization during tower build ⇒ no "decline if kernel arg has
a denominator"; radicals are first-class generators with explicit relations ⇒
no `$pc_rad` placeholder hack; independence proven ⇒ weak GCD never runs ⇒ no
hang guards; one path ⇒ no cascade.

### 0.3 What already exists (reuse — do NOT rebuild)

- **Tower builders:** `rt_tower_build_min(f,x,T,min_n)` → `RtTower`
  (`src/calculus/risch_tower.{c,h}`, transcendental Log/Exp/Tan/prim tower with
  triangularity check); `qa_resolve_extension_tower(alphas,n)` → `QATower`
  (`src/poly/qafactor.h:206`, algebraic-constant compositum, primitive element).
- **Independence decision:** `risch_rational_span(theta,gens,m,vars)`
  (`src/calculus/risch_structure.h:41`, Bronstein §9.3 Q-span).
- **Reduction primitives (proof the engine works — generalize, don't reinvent):**
  `flint_algebraic_field_normalize/_canonical/_together` (`src/poly/flint_bridge.c`
  :3241/:3416/:3664) already use `fmpz_mpoly_q_canonicalise` +
  `fmpz_mpoly_divrem_ideal` + `fmpz_mpoly_gcd` + `gr_ctx_init_nf` over a tower
  with a relation ideal. They are narrowly scoped; the rewrite generalizes ONE
  of them into THE engine.
- **The shadow parity harness** (`MATHILDA_RATCANON=shadow`, in `src/rat.c`):
  compares any new engine against the current cascade via an exact
  rational-function equality oracle (`rc_math_equal`: numerator of
  `Together[a-b]==0`, `RootReduce` fallback) plus a surface-form counter. This is
  the one genuinely-reusable piece already built — it is the **cross-phase
  validation oracle**. Keep it working; every phase must keep DIVERGED=0.
- **Expr↔FLINT bridge helpers:** `collect_all_symbols`, `expr_to_mpolyq`,
  `fmpz_mpoly_to_expr` (`src/poly/flint_bridge.c`); `extract_num_den`,
  `is_superficially_negative`, `negate_expr` (static in `src/rat.c`).

### 0.4 The output-convention rule (the thing that replaces the sign hack)

State it once; apply it in one function; test it as a spec. WL-faithful
`Together`/`Cancel`:

- **WL-faithful, NOT rationalized.** `Together[1/(x-Sqrt[2])]` stays
  `1/(x-Sqrt[2])`. Do NOT clear radicals from the denominator via norm/conjugate
  (that is `RootReduce`/`flint_algebraic_field_canonical`, out of scope). Reduce
  mod relations and gcd-cancel only.
- **Top-level `num`/`den` are expanded** polynomials in the generators
  (WL: `Together[1/(x-Sqrt2)+1/(x+Sqrt2)] = 2x/(x^2-2)`, expanded).
- **Kernel arguments are normalized but NOT force-expanded**:
  `Sqrt[1/2(x-y)]` stays factored inside; `Sqrt[1+(-1+u^2)/(1+u^2)]` has its
  radicand Together'd to `2u^2/(1+u^2)`. Kernels are atoms whose arguments were
  normalized during tower build.
- **Denominator sign:** one convention, chosen once (recommend: leading
  coefficient positive; flip only an all-negative denominator, term-wise, so
  `-(N)/(-D)` never surfaces and mixed-sign `x-1`, `a^2-ab` are untouched).
- **`Cancel` ≠ `Together`:** `Cancel[a/b+c/d]` leaves the sum uncombined
  (per top-level additive term); `Together` combines to one fraction.

The Phase-1 SPEC turns this into ~40 concrete `input → exact expected form`
assertions that every later phase is tested against.

### 0.5 Cross-phase invariants (never violate)

- **No arbitrary caps in the decision procedure** (`[[feedback_no_arbitrary_caps_decision_procedures]]`).
  Bounds are derived (structure theorem, degree bounds), never magic constants.
- **Builtin ownership contract** (`[[feedback_builtin_res_ownership]]`): never
  `expr_free(res)`; NULL-out reused args.
- **`Together`/`Cancel` never `Expand[Power[Plus,n]]`** as a pre-pass
  (`[[feedback_together_no_expand]]`).
- **Shadow harness stays green** (DIVERGED=0) at every phase boundary.
- **C99 strict**, guarded `M_*` constants, valgrind-clean (macOS dyld baseline
  noise only, `[[project_macos_valgrind_baseline_noise]]`).
- **Test conventions** (`[[project_tests_common_src_list]]`,
  `[[feedback_mathilda_test_cmake_conventions]]`): new `src/*.c` → `tests/CMakeLists.txt`
  `COMMON_SRC`; new `test_*.c` → its own `add_executable`; run the `*_tests`
  binary directly; scope runs to the change (`[[feedback_scope_tests_to_change]]`).

### 0.6 File placement decision

Put the engine in a new module `src/poly/ratcanon.{c,h}` so it gets its own unit
test binary. It needs a few helpers currently static in `rat.c`
(`extract_num_den`, `is_superficially_negative`, `negate_expr`) — expose them via
a small `src/rat_internal.h` (Phase 2, step 1) rather than duplicating. The
reduction (Phase 3) lives in `flint_bridge.c` (needs the FLINT headers) exposed
through `flint_bridge.h`; `ratcanon.c` orchestrates build → reduce → map-back.

---

## Phase 1 — Prototype + SPEC (de-risk; no production wiring)

**Goal.** Prove the one-tower-one-reduction covers the representative regimes
*before* committing to the rewrite, and produce the executable form-SPEC that all
later phases test against.

**Deliverables.**
1. A scratch builtin `RatCanon`Prototype`[expr]` (guarded, undocumented,
   in `src/poly/ratcanon.c`) that, for the four inputs below only, builds the
   tower by hand-wiring the existing builders and runs one FLINT reduction,
   returning the reduced form. Throwaway — deleted at end of Phase 3.
2. `tests/test_ratcanon_spec.c` — ≥40 `input → FullForm expected` assertions
   encoding §0.4, spanning: plain-Q; Log-tower; Exp commensurate
   (`(E^(2x)-1)/(E^x-1)→1+E^x`); Tan; Q(i) (`1/(x-I)+1/(x+I)→2x/(1+x^2)`);
   Q(√d) (`(x^2-2)/(x-Sqrt2)→x+Sqrt2`; `1/(x-Sqrt2)` stays); Q(ζn) cyclotomic;
   symbolic `Sqrt[k]` (`1/(b-a Sqrt k)+1/(b+a Sqrt k)`); cube root
   (`(y-1)/(y^(1/3)-1)→1+y^(1/3)+y^(2/3)`); nested kernel arg
   (`Sqrt[1+(-1+u^2)/(1+u^2)]`); `Cancel` vs `Together` (`a/b+c/d`); WL-faithful
   (radical kept in denom); a cyclotomic disguised-zero. Each row is the
   canonical WL form. **This file is the contract for Phases 3–4.**

**Prototype inputs to cover (proof of coverage):** `(1+x)/(1-x^2)`,
`1/(1+Log[x])+1/Log[x]`, `1/(x-I)+1/(x+I)`, `1/(b-a Sqrt[k])+1/(b+a Sqrt[k])`.

**Tests.** `test_ratcanon_spec` compiles and the prototype passes its four rows;
the remaining rows are `EXPECT`-marked but may be `KNOWN_TODO` until Phase 3.

**Acceptance gate.** Prototype returns the correct reduced form for all four
representative inputs via a single fmpz_mpoly_q(+ideal) reduction (not by calling
the old engines). If any regime cannot be reduced this way, STOP and record the
gap here before proceeding — that is a design finding, not a detail.

### Phase 1 — LANDED (2026-07-22)

- `src/poly/ratcanon.{c,h}` + builtin `RatCanonPrototype`; registered in
  `core.c`, added to `tests/CMakeLists.txt` COMMON_SRC.
- `tests/test_ratcanon_spec.c` — 44 hard assertions (plain-Q, transcendental,
  Gaussian, number field incl. hard cancels + WL-faithful, cyclotomic, symbolic
  radical, cube root, Cancel-vs-Together, `Q(√d)` d=5/7 no-crash) + 6
  prototype-correctness rows. Target `ratcanon_spec_tests`. All green.
- Prototype pipeline (substitute every kernel/constant/radical → fresh free
  symbol → `flint_rational_together` → map back → `eval` applies relations)
  reduces all four representative regimes correctly. Verified math-equal via
  `Together[RatCanonPrototype[e] - e] == 0`. Leak-clean (valgrind == baseline).

**KEY DESIGN FINDING (input to Phase 3).** The prototype's simple pipeline
(reduce with all generators FREE, apply relations only via `eval` at map-back)
covers plain-Q, transcendental, Q(i), and the *sum-of-conjugates* algebraic
cases (where combining clears the radical). It does NOT do the hard cancellation
of a *pre-formed* algebraic fraction — e.g. `Cancel[(x^2-2)/(x-Sqrt2)] -> x+Sqrt2`
— because that requires GCD in the quotient field `K[mainvar]`, not free-ring
gcd + relation-eval. The real builtins already do this via the number-field /
tower GCD engines (`flint_extension_gcd` family). Therefore **Phase 3's one
reduction = fmpz_mpoly_q for the free/transcendental part + field-GCD over
`K[mainvar]` for the algebraic part** (generalizing `flint_algebraic_field_*` +
the number-field GCD into ONE tower-field GCD). This is the crux and the main
technical risk; `test_ratcanon_spec.c`'s number-field rows are its acceptance
gate. Note also: `flint_algebraic_field_together` reduces with all gens free
(relations via eval), which is why it expands radicands and misses cube-root
cancellation — do NOT reuse it as the algebraic reducer.

---

## Phase 2 — The tower IR + builder (`rat_canon_build`)

**Goal.** One data structure and one builder that turns any supported `E` into
the canonical tower, reusing the existing builders. No reduction, no builtin
change.

**Deliverables.**
1. `src/rat_internal.h` exposing `extract_num_den`, `is_superficially_negative`,
   `negate_expr` from `rat.c` (make them non-static; declare here).
2. `src/poly/ratcanon.h`:
   ```c
   typedef enum { RCG_TRANSCENDENTAL, RCG_ALGEBRAIC } RcGenKind;
   typedef struct {
       Expr*  kernel;    /* Log[u]/E^w/Tan[u]/Sqrt[k]/... (surface form)     */
       Expr*  var;       /* fresh indeterminate t_i                          */
       RcGenKind kind;
       Expr*  relation;  /* NULL (transcendental) or t_i^q - g_i (algebraic) */
   } RcGen;
   typedef struct {
       QATower* K;       /* constant field Q(θ); NULL ⇒ K = Q               */
       RcGen*  gens; size_t n;   /* ordered: algebraic (leading) then transc. */
       Expr*   num; Expr* den;   /* E as num/den over the generators         */
       Expr*   var;              /* main variable x, or NULL                 */
   } RatCanonForm;
   RatCanonForm* rat_canon_build(const Expr* e);  /* NULL ⇒ unsupported     */
   void          rat_canon_free(RatCanonForm* f);
   ```
3. `rat_canon_build`: pre-normalize (`rt_powers_to_exp`, `rt_expand_logs`,
   `rt_expand_exp_sums`); collect algebraic constants → `K` via
   `extension_autodetect`/`qa_resolve_extension_tower`; collect transcendental +
   algebraic-function generators, **verify independence** via
   `risch_rational_span`, order them (algebraic leading); **recursively
   normalize each kernel's argument** (`rat_canon_build`+reduce on the argument —
   forward-declared, stubbed until Phase 3, so kernels are atoms with normalized
   args); substitute kernels → `t_i`; emit `num`,`den`.

**Tests — `tests/test_ratcanon_build.c` (extensive):**
- **Structure:** for ~30 inputs, assert `n`, each gen's `kind`, `relation`, and
  that `K` matches the expected number field (degree, minimal poly).
- **Independence:** `Log[x]`+`Log[x^2]` collapses to one gen (via expand);
  `E^x`+`E^(2x)` → one gen + `t^2`; `Log[2]`,`Log[3]`,`Log[6]` handled per the
  structure theorem; `Sqrt[2]`,`Sqrt[3]`,`Sqrt[6]` compositum has the right
  degree.
- **Round-trip:** substitute each `t_i`→`kernel`, θ→algebraic-number back into
  `num/den`; assert `expr_eq`(reconstructed, evaluated input).
- **Decline:** unsupported heads (e.g. `Gamma[x]`) → `rat_canon_build`==NULL.
- **Ownership:** valgrind-clean build/free over the corpus.

**Acceptance gate.** All structure + round-trip tests pass; valgrind clean;
builtins unchanged (shadow harness still DIVERGED=0 since nothing is wired).

### Phase 2 — LANDED (2026-07-22)

- `src/rat_internal.h` exposes `extract_num_den` / `is_superficially_negative` /
  `negate_expr` (de-static'd in `rat.c`).
- `src/poly/ratcanon.{c,h}`: `RatCanonForm` + `RcGen` + `rat_canon_build` /
  `rat_canon_free` / `rat_canon_roundtrip` / `rat_canon_subst_back`. One
  generator-profiling build: pre-normalize (`rt_expand_logs`/`rt_expand_exp_sums`
  — they BORROW their arg, free intermediates), substitute every kernel to a
  fresh `$rcgN$` symbol, `extract_num_den`, order algebraic-leading, pick a
  heuristic `var`. Commensurate exponentials collapse to a common fundamental
  (rational-gcd of the base coefficients; `E^(c·u) -> fund^(c/g)`); `Log[ab]`
  expanded; algebraic gens (`I`, `Sqrt`, `r^(p/q)`, roots of unity) carry the
  explicit relation `sym^q - radicand` / `sym^2 + 1`.
- `tests/test_ratcanon_build.c` (`ratcanon_build_tests`): structure (counts +
  kind + algebraic-leading order), independence (commensurate exps → 1 gen,
  `Log[x^2]` → `Log[x]`), relation-validity (every algebraic relation vanishes at
  its kernel; transcendental gens have none), round-trip
  (`Together[roundtrip - e] == 0`). All green.
- rat_tests + ratcanon_spec green; shadow DIVERGED=0; valgrind == baseline.

**DESIGN SIMPLIFICATION (vs the struct sketch above).** `RatCanonForm` drops the
`QATower* K` field: algebraic constants are NOT collapsed to one primitive
element. Each algebraic generator (constant or function) carries its own minimal
polynomial as a per-gen `relation`, and Phase 3 reduces mod the *set* of
relations (an ideal), ordering the algebraic gens as leading LEX vars so they
form a Gröbner basis — exactly `flint_algebraic_field_normalize`'s model. No
`qa_resolve_extension_tower` / `risch_rational_span` call is needed in the
builder; independence of the transcendental gens is achieved structurally
(log-expansion + exp-commensurability), and the round-trip test is the
correctness gate. Phase 3 may still build a primitive element internally if the
field-GCD needs `gr_ctx_init_nf`, but that is a reduction detail, not IR state.
Forward trig (`Sin/Cos/Sec/Csc/Sinh/Cosh`) is left un-substituted (algebraically
dependent — out of tower scope); only `Tan/Cot/Tanh/Coth` and Log/inverse-trig
are transcendental generators. Phase 3/4 declines forms carrying un-substituted
kernels to the classical path.

---

## Phase 3 — The reduction engine (`rat_canon_reduce`)

**Goal.** One reduction from `RatCanonForm` to the WL-faithful normal form,
generalizing `flint_algebraic_field_*` into THE engine, with the §0.4 convention
in one place.

**Deliverables.**
1. In `flint_bridge.{c,h}`: `Expr* flint_tower_reduce(num, den, gens, ngens,
   mode)` — build one `fmpz_mpoly` context over `{algebraic gens (leading), x,
   transcendental gens, free params}`; convert `num`,`den` (`expr_to_mpolyq`);
   `fmpz_mpoly_q_canonicalise`; reduce mod the relation ideal
   (`fmpz_mpoly_divrem_ideal` with the `t_i^q-g_i`/`m_θ` set — the code at
   `flint_bridge.c:3294/3479` is the template); `fmpz_mpoly_gcd`-cancel;
   read back (`fmpz_mpoly_to_expr`). WL-faithful: keep algebraic gens in the
   denominator, do NOT norm-rationalize.
2. `src/poly/ratcanon.c`: `Expr* rat_canon_reduce(RatCanonForm* f, RatCanonMode
   mode)` — map generators to `flint_tower_reduce`, map back (`t_i`→kernel,
   θ→algebraic number), apply `rc_output_normalize` (§0.4: sign convention +
   `Cancel` per-term vs `Together` combine). Complete the recursive
   kernel-argument normalization stubbed in Phase 2.
3. `Expr* rat_canon_normalize(const Expr* e, RatCanonMode mode)` =
   `build` → `reduce` → `free`; NULL if build declines.

**Tests — `tests/test_ratcanon_reduce.c` (extensive) + `test_ratcanon_spec`:**
- **Form SPEC:** every Phase-1 `test_ratcanon_spec` row now asserts exactly
  (flip `KNOWN_TODO` → `EXPECT`). This is the primary gate.
- **WL-faithful:** `1/(x-Sqrt2)` stays; `(x^2-2)/(x-Sqrt2)→x+Sqrt2`;
  `Cancel[a/b+c/d]` uncombined.
- **Regime matrix:** plain-Q, Log, Exp-commensurate, Tan, Q(i), Q(√d) for
  `d=2,3,5,6,7,8` (incl. the `TOGETHER_ALGEBRAIC_OVERFLOW` shape — must not
  crash), Q(ζn), `Sqrt[k]`, cube root, nested kernel arg, mixed
  (Q(i,√3)+Log), disguised zeros → `0`.
- **Open bugs close:** `Cancel[(a+b Sqrt k-a^2 k)/(b-a Sqrt k+b^2 k-a^2 k^2),
  Extension->Automatic]` <1s; the two Goursat named integrals
  (`CANCEL_IMPROVEMENT_PLAN.md §4`) integrate.
- **No hang:** each reduction wrapped in a wall-clock assertion (a test-side
  timeout, not a production cap).
- **Shadow parity:** run the full suite with `MATHILDA_RATCANON=shadow` pointed
  at `rat_canon_normalize`; DIVERGED=0 and user-facing form-diff=0.
- **Ownership:** valgrind on the reduce corpus + open-bug reproducers.

**Acceptance gate.** `test_ratcanon_spec` 100% green; DIVERGED=0 and form-diff=0
across the heavy suites in shadow; the two open bugs closed; valgrind clean.
Builtins still unchanged.

### Phase 3 — LANDED (2026-07-22)

- `flint_tower_reduce(frac, alg_syms, relations, n_alg)` in `flint_bridge.{c,h}`:
  one `fmpz_mpoly_q` combine + `fmpz_mpoly_divrem_ideal` reduce mod the relation
  set (algebraic gens leading) + re-canonicalise. Modeled on
  `flint_algebraic_field_normalize`.
- `rat_canon_reduce` / `rat_canon_normalize` + `rco_sign_normalize` (§0.4 sign)
  in `ratcanon.c`; test builtins `RatCanonNormalize` / `RatCanonCancel`.
  `rat_canon_normalize` splits a Cancel over top-level additive terms.
- `tests/test_ratcanon_reduce.c` (`ratcanon_reduce_tests`): accepted-form,
  decline, and parity-vs-builtin rows. All green. rat_tests/spec/build green;
  valgrind == build-baseline (a `builtin_of` test-helper borrow-leak was fixed).

**SCOPE DECISION (important for Phase 3b / Phase 4).** `rat_canon_reduce`
**accepts only a fully radical-FREE reduced result** and otherwise **declines to
NULL** (→ classical path in Phase 4). Covered by the one FLINT reduction:
plain-Q, transcendental (Log/Exp/Tan/inverse-trig), Gaussian, and the
sum-of-conjugates radical cases where the combine + ideal reduction eliminates
every generator (`1/(x-Sqrt2)+1/(x+Sqrt2) -> 2x/(x^2-2)`; the Sqrt[k] open bug).
Declined (residual radical remains): WL-faithful-kept radicals
(`1/(x-Sqrt2)`), coprime multi-radical sums, and **pre-formed algebraic-fraction
cancellations** (`(x^2-2)/(x-Sqrt2) -> x+Sqrt2`, cube roots). Rationale: the
field GCD in `K[x]` that those need is genuinely multi-case (number field vs
symbolic-radicand vs cube root — the existing `flint_extension_gcd` family), and
wiring it into `rat_canon_reduce` reintroduced the engines' inconsistent output
forms (`-(-Sqrt2-x)` etc.). Keeping the reduction **radical-free-only** preserves
ONE clean output convention; the field-GCD completion is deferred to **Phase 3b**
(a proper unified tower-field GCD) — until then Phase 4's classical fallback
covers the declined cases exactly as today. Parity: where `rat_canon_normalize`
does not decline, it is math-equal to the classical builtin (verified on ~35
inputs; pure numeric constants also decline — nothing to normalize).

Also note: `rat_canon` treats `Tan/Tanh/Cot/Coth` as transcendental generators,
so it reduces them via FLINT and emits `-2/(-1+Tan^2)` where the classifier's
classical Tan path emits `2/(1-Tan^2)` — math-equal, one consistent FLINT sign
convention. Phase 4 must update the SPEC's Tan row to the FLINT form (a
consistency improvement, reviewed via the baseline diff).

---

## Phase 4 — Switch the builtins and DELETE the zoo

**Goal.** `Together`/`Cancel` route entirely through `rat_canon_normalize`; the
classical cascade shrinks to a thin fallback for `build`-declines only; the
redundant engines and their hacks are removed so the code actually shrinks.

**Deliverables.**
1. `builtin_cancel_compute`/`builtin_together_compute` become:
   ```c
   Expr* r = rat_canon_normalize(arg, mode);   /* handles options structurally */
   if (r) return r;
   return mode==RC_CANCEL ? cancel_recursive(arg) : together_recursive(arg); /* sole fallback */
   ```
   Remove the Phase-2 (branch) prepended dispatcher, `rat_canon_dispatch`,
   `rat_canon_classify`, `rc_sign_normalize`, and the `MATHILDA_RATCANON`
   enable/shadow scaffolding from `rat.c` — superseded by the real engine.
   (Keep `rc_output_normalize` inside `ratcanon.c`; it is the *defined
   convention*, not a reconciliation hack.)
2. **Delete, one commit each, shadow/spec-verified after each:** the redundant
   FLINT fast-path blocks in `rat.c` (now subsumed), `flint_gaussian_*`,
   `flint_parametric_field_normalize` as a Together/Cancel entry,
   `flint_algebraic_field_together` as a Together entry, Phase E
   `qa_cancel_with_poly_radical`, `poly_find_radical_gen`/`poly_subst_radical_*`
   pass, `pick_best_tower_generator`, the `rat_has_dependent_*` bail predicates,
   and the degree-32/leaf-200/iter-50/`TRIGRAT_RADICAL_MAX` caps that only
   existed to survive the weak GCD.

**Tests.**
- Full suite green in normal mode (`rat`, `poly`, `simp`, `trigrat`,
  `radical_simplify`, `parfrac`, all `integrate_*`, `mpoly`).
- `test_ratcanon_spec` remains 100%.
- Regression vs the pre-rewrite baseline: keep a checked-in
  `tests/ratcanon_baseline.txt` of `input → form` for ~200 corpus inputs
  captured from `main` before Phase 4; assert no unintended form change (each
  intended change is reviewed and the baseline updated).
- valgrind on the reproducers; every deletion commit re-runs the affected
  `*_tests`.

**Acceptance gate.** Net LOC in `rat.c`+`flint_bridge.c`+`qafactor.c` DECREASES;
zero `poly_gcd_internal` calls reachable from `Together`/`Cancel` on the covered
regimes; full suite green; baseline diff reviewed.

---

## Phase 5 — Simplify de-dup + residual fixes

**Goal.** Remove the last duplication and close the remaining known defects.

**Deliverables.**
1. Delete `builtin_simplify`'s top-level autodetect+Together fast-path cascade
   (`src/simp/simp_builtins.c` ~780–928) — now redundant; Simplify seeds
   `Together`/`Cancel` as ordinary transform candidates (`simp_search.c`) and the
   memoized per-node descent calls the one engine.
2. Fix the `flint_algebraic_field_normalize` unbounded recursion (`Q(√d) d≥5`
   SIGSEGV, `TOGETHER_ALGEBRAIC_OVERFLOW.md`) if the Phase-3 engine did not
   already subsume that path; add its exact partial-fraction trigger as a test.
3. Reconcile the residual internal Gaussian coefficient-clearing forms
   (`(1/2 I)/(I+Log)`) if any downstream consumer depends on them (it did not in
   Phase-3 shadow; confirm and document).

**Tests.** `simplify_tests` (soft-fail count unchanged or reduced); a new
`test_together_algebraic_overflow.c` with the `d=5..9` shapes (no SIGSEGV);
full-suite green; valgrind.

**Acceptance gate.** Simplify has no private rational-normalization routing; the
SIGSEGV reproducer passes; suite green.

---

## Sequencing & fresh-context notes

- Order is strict: **1 → 2 → 3 → 4 → 5**. Each ends at a green, committed state.
- Every session: read §0, run `git log --oneline` to see which phase last landed,
  read that phase's acceptance gate to confirm it held, then execute the next.
- The shadow harness (`MATHILDA_RATCANON=shadow`) and `test_ratcanon_spec` are the
  two oracles that carry state across sessions — trust them over intuition.
- If a phase uncovers that the single engine cannot cover a regime WL-faithfully
  (e.g. a multivariate algebraic GCD FLINT lacks), STOP and record it in the
  relevant phase §; do not reintroduce a special-case engine to paper over it —
  that is the failure mode this rewrite exists to end.
