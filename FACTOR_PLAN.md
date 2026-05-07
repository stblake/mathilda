# PicoCAS Factoring Pipeline — Status and Plan

This document covers the polynomial factoring subsystem of picocas:
its current architecture, the recent improvements that landed, and
the remaining work to bring multivariate factoring up to par with
mainstream computer algebra systems.

The audience is anyone resuming work on this code: it lays out where
each piece lives, how the call graph connects them, and what design
decisions are still open.

---

## 1. What this subsystem does

`Factor[poly]` and `FactorSquareFree[poly]` are picocas built-ins for
polynomial factorisation over the integers `Z`. They underpin a broad
slice of the system:

- **`Simplify`** runs `Factor` on candidate intermediate forms and
  picks the smallest-by-`SimplifyCount`. Slow `Factor` is felt as
  slow `Simplify`.
- **`TrigFactor`** invokes `Factor` after rewriting reciprocal trig
  heads as `Sin`/`Cos` ratios, so its inner pipeline is gated by
  multivariate factoring quality.
- **Rational simplification** (`Together`, `Cancel`, `Apart`) calls
  `Factor` on numerators and denominators.

Anything that improves multivariate `Factor` therefore propagates
into many user-visible operations.

---

## 2. File layout

| File | Lines | Role |
|---|---|---|
| `src/facpoly.c` | 2131 | The orchestration layer + the existing univariate Berlekamp-Zassenhaus implementation. The single place where `builtin_factor` lives. |
| `src/zupoly.{c,h}` | ~900 / 200 | New: univariate `Z[x]` polynomials with arbitrary-precision `mpz_t` coefficients. The substrate for the bivariate pipeline. Includes a private `QUPoly` (rational coefficients) used internally by `zupoly_diophantine`. |
| `src/bpoly.{c,h}` | 480 / 132 | New: bivariate `Z[x, y]` polynomials, dense in the main variable, each `x`-coefficient stored as a `ZUPoly` in `y`. |
| `src/mvfactor.{c,h}` | ~210 / 80 | New: multivariate factoring orchestration. Currently exposes `bpoly_hensel_lift_2` and `bpoly_hensel_lift_multi`. The future home of `mvfactor_factor` (the public entry point). |
| `tests/test_facpoly.c` | 64 | Original sanity tests (univariate factoring + a couple of multivariate cases). 2 pre-existing failures unrelated to this work. |
| `tests/test_factor_phase0.c` | 307 | 17 tests for the monomial-content extractor and the irreducibility short-circuit, with a perf budget on irreducible inputs. |
| `tests/test_zupoly.c` | ~700 | 38 tests: ZUPoly arithmetic, division, GCD, content, evaluation, shift, Expr round-trip, and the Diophantine primitive (6 dedicated cases). |
| `tests/test_bpoly.c` | 452 | 23 tests covering BPoly arithmetic, exact division, truncation `mod y^k`, evaluation `y -> α`, shift `y -> y + α`, and Expr round-trip. |
| `tests/test_mvfactor.c` | ~370 | 9 tests covering bivariate two-factor Hensel (6) and multifactor pair-and-recurse (3). |

---

## 3. The current `facpoly.c` architecture

### 3.1 Public entry points

- **`builtin_factor`** (line 1138): top-level `Factor[...]` handler.
  Routes through `Together` -> `Numerator`/`Denominator` and applies
  the factoring pipeline to each.
- **`builtin_factorsquarefree`** (line 284): `FactorSquareFree[...]`.
- **`builtin_factorterms`** / **`builtin_factortermslist`** (lines
  1456 / 1436): `FactorTerms[...]`, content extraction with respect
  to a chosen variable set.
- **`bz_factor_to_expr`** (line 2057): the univariate Z[x] entry
  point used both internally and by tests; wraps the
  Berlekamp-Zassenhaus pipeline.
- **`facpoly_init`** (line 1490): registers all of the above as
  builtins with their attribute flags.

### 3.2 The factoring pipeline

```
Factor[P]   (builtin_factor)
    |
    v
Together[P]  -->  num / den
    |
    v
For num, den separately:
  if v_count == 1:           bz_factor_to_expr     (univariate BZ)
  else:                      FactorSquareFree[.] then heuristic_factor
                                            |
                                            v
                                  heuristic_factor
                                  (multivariate orchestrator)
```

### 3.3 `heuristic_factor` (line 1032) — the multivariate orchestrator

Tried strategies, in order. Each falls through if it cannot make
progress:

1. **Structural recursion** through `Times`, `Power`, and list-like
   heads (`List`, `Equal`, `Less`, `And`, etc.). Trivial atoms / non-
   polynomials are returned as-is.

2. **Collect variables** and verify the input is a polynomial.

3. **`poly_content`** (in `poly.c`): integer / rational content w.r.t.
   the chosen main variable. Recurses on the content and primitive
   part separately.

4. **`factor_monomial_content`** (line 614, **Phase 0**): extract
   the largest monomial `v_1^{e_1} * ... * v_k^{e_k}` shared by every
   term of a Plus. This is the cheapest factorisation step
   (O(n_terms * n_vars) with no polynomial arithmetic) and was the
   missing case that caused inputs like `3 a^2 b - 3 b - b^3` to
   return unfactored. It uses `monomial_collect` (line 569), a
   recursive walker that handles picocas's nested-`Times` canonical
   form (`Times[3, Times[Power[a,2], b]]`) which the older
   `extract_monomial` helper failed to descend into.

5. **`factor_degree_one`** (line 479): if some variable appears with
   degree exactly 1, factor `L * v + C` into `gcd(L, C) * (...)`.

6. **Univariate dispatch** (line 1063): if `v_count == 1`, hand off
   to `bz_factor_to_expr` (Berlekamp-Zassenhaus).

7. **`factor_binomial`** (line 365): two-term `Plus` with shared
   variable structure raised to a common k-th power. Handles
   `x^2 - 4 y^2`, `x^3 + y^3`, etc.

8. **`is_likely_irreducible_multivariate`** (line 838, **Phase 1**):
   Hilbert-irreducibility-based short-circuit. For each variable as
   "main", substitute the others with several integer values, factor
   the univariate image via `bz_factor_to_expr`, and require at
   least 2 confirmations of "image is irreducible" (with squarefree
   image and unchanged main-variable degree) before declaring the
   input irreducible. Replaces what would otherwise be wasted work
   in `factor_roots` on inputs with no factorisation.

9. **`factor_roots`** (line 985): the legacy fallback — try linear
   factors `(x - c)` for `c in {±1, ±2, ±3, ±4, ±6}` and `(x - c·y)`
   for each variable pair, via `exact_poly_div`. O(v_count² * 10)
   trial divisions. **Still present** because some existing tests
   (notably `Factor[2x^3 y - 2 a^2 x y - 3 a^2 x^2 + 3 a^4]`) rely
   on it. To be retired once the bivariate Hensel pipeline is in.

10. **Fallback**: return `P` unchanged.

### 3.4 The univariate `Z[x]` Berlekamp-Zassenhaus engine

Lives in `facpoly.c` lines 1501-2055. Self-contained, with its own
type:

- `UPoly` (line 1503): `int64_t`-coefficient univariate. Designed for
  modular arithmetic in finite-field algorithms; coefficient growth
  is bounded by the modulus, so 64-bit machine integers are safe.
- Modular arithmetic (`upoly_div_rem_mod`, `_add_mod`, `_sub_mod`,
  `_mul_mod`, `_gcd_mod`, `_pow_mod_poly`).
- `cz_edf` / `cz_ddf` (lines 1556 / 1597): Cantor-Zassenhaus
  equal- and distinct-degree factorisation in `GF(p)[x]`.
- `upoly_xgcd_mod` (line 1640): extended Euclidean over `GF(p)`.
- `hensel_lift` / `multifactor_hensel_lift` (lines 1674 / 1731):
  p-adic Hensel iteration for univariate factorisation.
- `factor_zassenhaus` (line 1805): the full Berlekamp-Zassenhaus
  pipeline (DDF -> EDF -> Hensel lift -> Zassenhaus recombination).
- `bz_factor_to_expr` (line 2057): public wrapper, converts
  picocas `Expr*` <-> `UPoly`, runs `factor_zassenhaus`.

Limitations of this layer worth noting:
- `int64_t` coefficients limit lifted moduli to roughly `10^15`.
- The prime is currently fixed at 13 with attempts up to 20 nearby
  primes (line 1827).
- Multifactor Hensel uses linear merging (not the optimal balanced
  product tree).

These constraints are acceptable for the univariate case but rule
out reusing the same machinery directly for the multivariate
extension, which motivates the `ZUPoly` and `BPoly` types.

---

## 4. Phase 0 + Phase 1: what shipped

The two improvements landed since project inception:

### 4.1 Monomial-content extraction (Phase 0)

- **Function**: `factor_monomial_content` (facpoly.c:614).
- **Helper**: `monomial_collect` (facpoly.c:569) — recursive
  decomposer robust to nested `Times`.
- **Wire**: invoked in `heuristic_factor` after content extraction,
  before `factor_degree_one` (facpoly.c:1077).
- **Effect**:
  - Correctness fix: `Factor[3 a^2 b - 3 b - b^3]` now returns
    `b · (3 a^2 - 3 - b^2)` (previously returned input unchanged).
  - Performance: also indirectly cuts time for cases that hit it,
    because the residue is smaller and downstream strategies handle
    it faster.

### 4.2 Hilbert-irreducibility short-circuit (Phase 1)

- **Function**: `is_likely_irreducible_multivariate` (facpoly.c:866).
- **Helpers**: `eval_others_at_alpha` (facpoly.c:839),
  `count_nontrivial_factors`, `univariate_squarefree`.
- **Tunables** (facpoly.c:864): `IRRED_TRY_POINTS = 7`,
  `IRRED_CONFIRM_COUNT = 2`.
- **Wire**: invoked in `heuristic_factor` immediately before
  `factor_roots` (facpoly.c:1117).
- **Effect**: irreducible bivariate inputs like `x^2 + y^2 + 1`
  finish in single-digit milliseconds instead of hundreds.

### 4.3 Measured impact on the user's reported case

| Workload | Pre-Phase-0 | Post-Phase-1 | Speedup |
|---|---|---|---|
| `Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]` | 1158 ms | 534 ms | 2.2x |
| `Factor[3 a^2 b - 3 b - b^3]` | 364 ms (wrong) | 19 ms (correct) | 19x + correctness |
| `Factor[3 Cos[x]^2 Sin[x] - 3 Sin[x] - Sin[x]^3]` | 583 ms | 32 ms | 18x |
| `Factor[3 Cos[x]^2 - 3 - Sin[x]^2]` | ~200 ms | 13 ms | 15x |
| `Factor[x^2 + y^2 + 1]` | ~300 ms | 5 ms | 60x |

**Regressions**: zero. The 4 pre-existing failures on `main` are
unchanged (Cross with mixed Complex/Real, PolynomialMod with
Gaussian integers, two Factor cases where the denominator's
multivariate structure isn't being passed through `bz_factor_to_expr`
correctly).

---

## 5. Phase 2 foundations — `ZUPoly` and `BPoly`

These two new types are the substrate for the proper multivariate
Hensel lift. Neither is wired into `builtin_factor` yet; they are
self-contained and exhaustively unit-tested.

### 5.1 `ZUPoly` — univariate Z[x] with mpz_t coefficients

- **Header**: `src/zupoly.h`
- **Implementation**: `src/zupoly.c` (~626 lines)
- **Key design choice**: arbitrary-precision coefficients. The
  bivariate Hensel lift's coefficient growth (Mignotte bound) cannot
  be handled by 64-bit machine integers in general; silent overflow
  here would produce wrong factorisations. `mpz_t` eliminates that
  class of bug.
- **Operations**:
  - Construction: `zupoly_new`, `zupoly_zero`, `zupoly_from_int`,
    `zupoly_copy`.
  - Coefficient access: `zupoly_setcoef[_si]`, `zupoly_getcoef`,
    `zupoly_normalize`.
  - Predicates: `zupoly_is_zero`, `zupoly_eq`, `zupoly_cmp`.
  - Arithmetic: `_add`, `_sub`, `_mul`, `_neg`, `_scale[_si]`.
  - Division: `_divrem_monic` (long division; succeeds for monic
    divisors and exact division), `_pseudodivrem` (always succeeds
    over Z by pre-multiplying), `_divexact` (returns NULL if not
    exact).
  - GCD: `zupoly_gcd` via the **subresultant pseudo-remainder
    sequence** (Brown-Collins). Returns the primitive part with
    positive leading coefficient.
  - `zupoly_content`, `zupoly_primitive_part`.
  - Evaluation: `zupoly_eval[_si]`. Horner.
  - Shift: `zupoly_shift_si` for `p(x + α)`. Horner expansion in
    `(x + a)`. (The first version of this had a subtle in-place
    iteration-direction bug that the unit tests caught — see
    `test_shift` in `test_zupoly.c`.)
  - `expr_to_zupoly`, `zupoly_to_expr` for picocas conversion.
- **Tests**: 32 in `tests/test_zupoly.c`, all passing.

### 5.2 `BPoly` — bivariate Z[x, y]

- **Header**: `src/bpoly.h`
- **Implementation**: `src/bpoly.c` (~480 lines)
- **Storage**: dense array of `ZUPoly*` indexed by x-degree. Each
  `cx[i]` is a y-polynomial. This layout is well suited to:
  - Hensel lifting in y (we routinely truncate every `cx[i]`
    modulo y^k -- a coefficient-by-coefficient operation).
  - Polynomial division viewing P as Z[y][x] (reduce leading
    x-coefficient using ZUPoly arithmetic at each step).
  - Substitution y -> α (collapse each cx[i] to a constant in Z,
    yielding a univariate polynomial in x).
- **Operations**:
  - Construction: `bpoly_new`, `bpoly_zero`, `bpoly_copy`,
    `bpoly_set_xcoef` (transfers ownership of the ZUPoly).
  - Predicates and degree: `bpoly_is_zero`, `bpoly_eq`,
    `bpoly_deg_x`, `bpoly_deg_y`, `bpoly_lc_x`.
  - Arithmetic: `_add`, `_sub`, `_mul`, `_neg`, `_mul_zupoly`,
    `_mul_truncate_y`.
  - Exact division: `bpoly_divexact` (long division viewing as
    Z[y][x]; each step's leading coefficient division must succeed
    in `zupoly_divexact` or the whole thing fails).
  - Truncation: `bpoly_truncate_y(p, k)` returns `p mod y^k`.
  - Substitution: `bpoly_eval_y_si` (-> ZUPoly in x),
    `bpoly_shift_y_si` (returns BPoly with y -> y+α).
  - `expr_to_bpoly`, `bpoly_to_expr` for picocas conversion.
- **Tests**: 23 in `tests/test_bpoly.c`, all passing.

---

## 6. Project status

| Phase | Description | Status |
|---|---|---|
| 0 | Monomial-content extraction in `heuristic_factor` | **Done** |
| 1 | Hilbert-irreducibility short-circuit before `factor_roots` | **Done** |
| 2a | `ZUPoly` (univariate Z[x] with mpz_t) | **Done** |
| 2b | `BPoly` (bivariate Z[x, y]) | **Done** |
| 2c | ZUPoly extended GCD over Q (Diophantine primitive) | **Done** |
| 2d | Bivariate Hensel lift (two-factor) | **Done** |
| 3a | Multifactor Hensel lift (pair-and-recurse) | **Done** |
| 3b | Zassenhaus recombination (subset trial division) | **Done** |
| 3c | Orchestrator `mvfactor_try_bivariate_monic` (alpha pick + lift + shift) | **Done** |
| 4 | n-variate recursion (n >= 3) | **Done** (minimum viable) |
| 5a | Wire bivariate-Hensel into `heuristic_factor` (monic only) | **Done** |
| 5b | Wang's leading-coefficient correction (handles non-monic) — same as Phase F1 | Not started — see §12 |
| 5c | Retire `factor_roots` once F1 + F2 + F3 land | Not started — see §12 |
| F1 | Wang's lc correction (replaces 5b) — bivariate non-monic via leading-coefficient pre-distribution | **Stages 1 + 2 + 3 done** — Stage 1 (LC = -1 negate path): pre-negate P, run the monic Hensel on -P, absorb the sign into the highest-x-degree factor via `Expand[-1 · factor]`.  Stage 2 (constant `\|a\| > 1` LC): Wang's monic-substitution recipe.  Form `Q = a^(d-1) · P(x/a, y)` (monic in x with integer coefficients), factor Q via the existing monic Hensel into G_1 ... G_r, then recover the true factors via `F_i = G_i(a·x, y) / cont_Z(G_i(a·x, y))`.  Stage 3 (polynomial-in-y LC, MVP scope): predicted-LC two-factor Hensel.  When `lc_x(P)(y)` is a non-constant polynomial `A(y)`, factor `A` over `Z[y]`, find `α` with `A(α) = +1` so the squarefree univariate image `P(x, α)` has monic factors `u`, `v`, then enumerate distributions of `A`'s irreducible factors between predicted leading coefficients `q_u, q_v` (with `q_u · q_v = A`) such that `q_u(α) = q_v(α) = +1`.  For each surviving distribution, the new `bpoly_hensel_lift_2_lc` runs the Hensel iteration with the leading-x coefficient of each `Δu` correction PINNED to the y^k coefficient of `q_u` and only the lower x-degree part solved via Diophantine — that keeps `lc_x(U)(y) = q_u(y)` invariant across the lift, so `Π U_i = P` exactly without the `A^(r-1)` content-division step from textbook Wang.  Unlocks bivariates like `Factor[Expand[(xy+1)(xy+2)]]`, `Factor[Expand[((y²+1)x+1)(x+3)]]`, `Factor[Expand[((y+1)x+1)((y+1)x+2)]]`.  MVP limitations: r = 2 only (more factors fall through), both univariate factors must be monic (often achievable when `A(α) = ±1`), `\|cont(A)\| = 1`, and inputs that have a non-trivial monomial content fall through so the `heuristic_factor` Phase 0 path produces the canonical fully-factored form. |
| F2 | True multivariate Hensel for n ≥ 3 — MPoly type + n-variate Hensel iteration | **MVP done** — Stage 2.1 (MPoly substrate `src/mpoly.{c,h}`, ~600 LOC + 38 tests) gives a sparse Z[x_1,...,x_n] type with mpz_t coefficients in lex descending order, full arithmetic, substitution, Expr round-trip.  Stage 2.2 (trivariate Hensel `src/mvfactor3.{c,h}`, ~500 LOC + 9 tests) adds a bivariate Diophantine solver via Hensel iteration on top of zupoly_diophantine, plus a trivariate two-factor lift (`mpoly_hensel_lift_3_2`).  Stage 2.3 (`factor_trivariate_via_mhensel` in `src/facpoly.c`) wires the trivariate path into `heuristic_factor` between Phase 4 and the irreducibility gate.  Algorithm: pick monic-in-main variable, try (α_y, α_z) tuples, factor univariate image via existing BZ, lift to bivariate via `bpoly_hensel_lift_2`, lift to trivariate via `mpoly_hensel_lift_3_2`, verify product.  **Unblocks the user-reported case 2 from §12**: `Factor[Expand[(zx - x² - y²)(3z + 4xy - y²)]]` 6 s (returned input unchanged) → 0.33 s (correctly factored).  Other simpler trivariate cases like `(z+xy)(z-xy)`, `(z+x)(z+y)`, `(z²+x)(z+y)` factor in single-digit ms.  Drive-by fix: corrected `extract_monomial` to descend into nested Times (was missing `Times[-1, Times[Power[x,2], Power[y,2]]]` factor structure, breaking `factor_binomial` on certain trivariate inputs).  MVP scope: monic-in-main two-factor only; non-monic Wang correction and r ≥ 3 recombination follow in subsequent phases. |
| F3 | Bivariate Hensel performance — incremental U·V update, Mignotte fast-fail, sorted us[] | Not started — required for F1 to be fast, ~200-400 LOC |
| F4 | Faster multivariate FactorSquareFree — cheap squarefree pre-check before full GCD | **Phase 1 done** — `sqfree_cheap_check` in `facpoly.c` substitutes integer values for the non-main variables, computes the univariate gcd-with-derivative, and skips the expensive multivariate `gcd(pp, pp')` when the image is squarefree at any of seven test alphas.  Sound after content extraction (a constant-in-x repeated factor would have divided content; a non-constant repeated factor specialises to a non-squarefree image at all but a finite alpha set).  Measured: 4-variable squarefree Hensel-style input 6.27 s → 0.95 s (6.6×); 3-variable squarefree typical 60-150 ms unchanged because the multivariate GCD was not the bottleneck there.  Non-squarefree inputs fall through to the original Yun loop; tests in `tests/test_facpoly.c` cover both branches. |
| F5 | Recombination cap & heuristics — singleton ordering, partial-lift signal, budget cap | **Done** — partial-lift signal landed in `lift_multi_internal` (`src/mvfactor.c`).  Singleton ordering was already done as F3a.  The lift-2 internal helper now reports a `completed_iters` flag distinguishing "lift bailed early on Diophantine non-integer or Mignotte coefficient overflow" from "lift completed all B y-iterations and only failed the final `U·V == P` verify".  When an entire k-pass produces no completion signal AND k ≥ 2, the recombination loop short-circuits and returns P as irreducible bivariately.  k = 1 is always exhausted so legitimate singleton lifts still fire; k = 2 is always tried before any short-circuit so all current k = 2 recombination cases (e.g. `(x²-y)(x²-4y)`) remain correct.  Measured: `Factor[Expand[x²⁴+y²-1]]` 117 → 84 ms (-29 %), `Factor[Expand[x¹²+y⁴-1]]` 24 → 20 ms (-17 %), `Factor[Expand[x¹²+y²-1]]` 23 → 20 ms (-13 %).  Larger savings on bivariate-irreducible inputs whose univariate image factors into 6+ pieces (cyclotomic-rich); no regression on the existing recombination test suite. |
| 5d | Context-aware num/den scope + explicit threading over Less/Equal/And in `builtin_factor` | **Done**.  Direct Factor calls use separate num/den variable scopes (correct full factorisation including foreign denominator vars).  Simplify-internal Factor calls reuse the numerator's scope for the denominator (avoids the TrigRoundtrip slow path on `Tan[2x]` forms; see C4).  Discriminator is the presence of an active Factor memo (which Simplify pushes).  Independent threading over Less/Greater/Equal/And/Or is now done at the top of `builtin_factor` so the cubic in `Factor[1 < expr < 2]` factors cleanly. |
| 6 | Per-call Factor memo inside `Simplify` | **Done** |
| 7 | TrigFactor Path-B skip heuristic on inputs with no compound trig structure | **Done** -- 17 % speedup on user case, 18 % on Tan double-angle. |
| 8 | TrigRoundtrip memo (reuses FactorMemo with `TrigRoundtrip[X]` keys) | **Done** -- additional 8-24 % across the user case, Tan, and Sin[2x]/Sin[x]. |
| 9 | TrigFactor + TrigExpand memos (same FactorMemo, different head keys) | **Done** -- `Sin[2x]/Sin[x]` 208 → 167 ms (-20 %), `Cos[3x]/Cheb` 177 → 152 ms (-14 %), Tan 1857 → 1793 ms (-3 %).  User case effectively unchanged (no duplicate inputs). |
| 10 | PythagReduce / PythagSquareComplete / HalfAngle memos via shared `simp_memo_wrap` helper.  PythagReduce had the highest volume (219 calls / 46 unique = 5× dedup factor). | **Done** -- Pythag^4 73 → 44 ms (-40 %), Sinh²-Cosh² 64 → 47 ms (-27 %), user case 360 → 323 ms (-10 %). |
| 11 | (Attempted) generic FactorMemo at `traced_call_unary` level to cover all named transforms (Together, Cancel, Apart, Expand, etc.) | **Reverted** -- per-transform memos already cover the high-volume duplicates; the additional malloc/hash/lookup at every transform-dispatch call exceeded the marginal gain on cheap Together/Cancel/Apart calls.  Documented in the function comment so future attempts know not to retry without first measuring. |
| 12 | Canonical-form keying for trig memos (TrigFactor, TrigExpand, TrigRoundtrip) via `Together(Expand(.))` | **Done** -- `1/8 (-18 a + 6 b)` and `-9/4 a + 3/4 b` now hash to the same memo key.  User case 325 → 259 ms (-20 %), Tan 1693 → 1410 ms (-17 %), Cos[3x]/Cheb 122 → 94 ms (-23 %). |
| 13 | (Attempted) extend canonical-form keying to `simp_memo_wrap` (PythagReduce, PythagSquareComplete, HalfAngle). | **Reverted** -- correctness regression: PythagReduce's rule `-1 + Cos[x]^2 + r___ :> -Sin[x]^2 + r` matches a specific surface pattern; running `Expand` first distributes `a (-1 + Cos^2)` into `-a + a Cos^2`, where the `-1` no longer appears as a Plus arg.  The user case stopped reducing fully (`Sin³ + Sin[3x] - 3 Sin[x] → 3 Sin[x] (-1 + Cos^2)` instead of `-3 Sin³`).  The trig memos work because their impls internally normalise (Together / TrigToExp) before pattern matching; the rule-based simp transforms don't.  Documented in the function comment in `simp_memo_wrap`. |
| 14 | TrigToExp-output keying for TrigRoundtrip memo.  Replaces Phase 12's Together(Expand(input)) key with TrigToExp(input) -- a stronger canonicalisation that folds in trig identities (Cos[x]^2 = (1+Cos[2x])/2 produces the same exp form).  Pipeline restructured so TrigToExp runs first; its output is the cache key; remaining stages (Together/Cancel/ExpToTrig) only run on cache miss. | **Done** -- user case 261 → 240 ms (-8 %), Tan 1411 → 1271 ms (-10 %), Cos[3x]/Cheb 93 → 101 ms (~flat).  TrigRoundtrip's pipeline saves 100+ ms per dedup. |
| 15 | (Attempted) extend Phase 14's TrigToExp keying to TrigFactor / TrigExpand. | **Reverted** -- Phase 14 worked because TrigToExp was the first stage of TrigRoundtrip's pipeline anyway (free-rider).  For TrigFactor / TrigExpand, TrigToExp is pure overhead (~5 ms per call × ~30 calls in the user case = 150 ms).  Most TrigFactor / TrigExpand inputs are already unique post-canonicalisation, so dedup gains don't compensate.  Documented as a comment in `trig_canonicalize`. |
| 16 | TrigRoundtrip explosion guard: abort the round-trip when TrigToExp expanded the input by more than 5x. | **Done** -- TrigToExp is structurally expanding (e.g., `Cos[x] Cos[y]` complexity 7 → 77 leaves of exponentials in 2 vars).  Subsequent Together / Cancel / ExpToTrig produces explosive Cosh / Sinh forms with imaginary arguments that aren't useful for the simp candidate-set search.  Skipping the round-trip on these inputs preserves correctness (the input form remains in the candidate set) and saves the expensive computation.  Tan case 1530 → ~700 ms (-54 %), `Cos[x+y] - Cos Cos + Sin Sin` 926 → 286 ms (-69 %), `Sin[2x]/Sin[x]` 149 → 123 ms.  User case unchanged (its TrigToExp expansion is ~3x, below the 5x threshold).  (Initially shipped, then reverted on a false-positive integrals_tests regression, then re-applied after the user clarified the integrals failures were from unrelated changes outside this work.) |

A separate concern, orthogonal to the above:

| Cleanup | Description | Priority |
|---|---|---|
| C1 | Resolve the 2 pre-existing `Factor` regression failures.  | **Done** via the context-aware approach in 5d. |
| C2 | Replace fixed prime `p = 13` in `factor_zassenhaus` with a smarter prime selector | Low |
| C3 | Use balanced Hensel product tree instead of linear multifactor | Low |
| C4 | Investigate Simplify's TrigRoundtrip behavior on multivariate trig fractions when their denominator factors completely (1 - Tan[x]^2 -> (1-Tan[x])(1+Tan[x])).  Apparently this triggers a path that produces a 12-term polynomial in Cos[k x], Sin[k x] for various k.  Currently worked around by 5d's context-aware scope (Simplify uses the conservative shared scope).  An upstream fix would make this workaround unnecessary. | Medium |

---

## 7. Next steps in detail

### 7.1 Phase 2c — ZUPoly extended GCD over Q

The bivariate Hensel iteration's Diophantine step solves

```
ΔU(x) * v(x) + ΔV(x) * u(x) = E(x)
```

for `u, v ∈ Z[x]` with `gcd(u, v) = 1`. Standard recipe needs
`S, T ∈ Q[x]` with `S(x) u(x) + T(x) v(x) = 1`, then

```
ΔU = E * T mod u
ΔV = (E - ΔU * v) / u    (exact division when u is monic)
```

Three implementation options, in increasing order of complexity and
correctness:

- **Option A**: Compute `S, T` over `GF(p)` using `upoly_xgcd_mod`
  (already in `facpoly.c`), then p-adic Hensel-lift to `Z[x]`. Reuses
  existing battle-tested machinery; downside is double-Hensel
  (p-adic + y-adic) layering.
- **Option B**: Compute `S, T` over `Q[x]` directly using extended
  Euclidean with rational coefficients (mpq_t). Self-contained;
  downside is `mpq_t` is slower than `mpz_t` and we end up clearing
  denominators anyway.
- **Option C**: Track a global denominator and do everything in
  `Z[x]` using subresultant pseudo-remainders. The half-extended
  Euclidean variant. Most efficient but trickiest to implement.

Recommended: start with Option A. Wrap `upoly_xgcd_mod` -> p-adic
lift in a `ZUPoly`-typed API; the existing primes `p = 13` family
should suffice for most inputs.

### 7.2 Phase 2d — Bivariate Hensel lift

Two-factor first:

```
Inputs:
  P(x, y) ∈ Z[x, y], primitive in x.
  α: lc_x(P)(α) ≠ 0 and P(x, α) is squarefree.
  u(x), v(x): bz_factor_to_expr applied to the squarefree image.

Algorithm:
  Shift y -> y + α so we lift around y = 0.
  Initial: U_0 = u, V_0 = v.
  For k = 1, 2, ..., deg_y(P):
    E_k = ((P - U_{k-1} * V_{k-1}) mod y^{k+1}) / y^k   (a polynomial in x)
    Solve ΔU * v + ΔV * u = E_k via Phase 2c.
    U_k = U_{k-1} + y^k * ΔU
    V_k = V_{k-1} + y^k * ΔV
  Shift y -> y - α to restore.
  Return (U, V).
```

Complications:
- **Non-monic P**: leading-coefficient correction (Wang's trick) —
  multiply P through by `lc_x(P)^(r-1)` before lifting, divide back
  at the end.
- **Multifactor (r >= 2)**: recursive partition into two halves,
  Hensel-lift each, recurse.
- **Coefficient bound**: Mignotte bound on the lifted factors gives
  an explicit choice of iteration count; if exceeded, the lift
  failed (extraneous factor in the modular factorisation).

### 7.3 Phase 3 — Zassenhaus recombination

After Hensel, some lifted modular factors may not correspond to true
multivariate factors (the modular factorisation is finer than the
multivariate one). Try products of subsets — each subset whose
product trial-divides P exactly is a true factor. Identical
structure to the recombination loop in `factor_zassenhaus`
(facpoly.c:1885).

### 7.4 Phase 4 — n-variate recursion

For inputs in `Z[x_1, ..., x_n]` with `n >= 3`, treat as bivariate
over `Z[x_3, ..., x_n][x_1, x_2]`: specialise the trailing variables
to integers, factor bivariately, lift each factor recursively. Cap
at n = 3 in the first cut; n > 3 falls back to the existing partial
factorisation rather than guessing.

### 7.5 Phase 5 — wire and retire

Once Phase 2d-3 are solid:
1. Plug the new pipeline into `heuristic_factor` between the
   irreducibility short-circuit and `factor_roots`.
2. Verify with the existing test suite that all previously-passing
   `factor_roots` cases are now handled by the new pipeline.
3. Delete `factor_roots`. Cleanup commit.

### 7.6 Phase 6 — Simplify-level Factor memo

Inside `simp_dispatch`'s candidate-set search (`simp.c`), the same
intermediate polynomial is factored multiple times across different
candidate paths. A per-call memo keyed by `expr_hash`/`expr_eq`
(same pattern as the existing `SimpMemo` at `simp.c:2759`) would
deduplicate, giving a multiplier on top of the algorithmic
improvements. Independent of Phase 2-5 and can land at any point.

---

## 8. Testing inventory

| Test file | Cases | Coverage |
|---|---|---|
| `tests/test_facpoly.c` | 16 | Original sanity — `FactorSquareFree` and `Factor` on canonical inputs.  All passing (the 2 pre-existing failures got fixed by Phase 5d's tactical work). |
| `tests/test_factor_phase0.c` | 22 | Monomial-content correctness (9), irreducibility-short-circuit timing (4), regression coverage (3), Phase-5 bivariate Hensel (5), end-to-end `Simplify` user case (1).  All passing. |
| `tests/test_factor_baseline.c` | 16 | End-to-end correctness + performance baseline.  Direct Factor (7) + Simplify integration (7) + regression cases (2).  Generous timing budgets (3-5x typical) for stability across machines.  Captures the cumulative state after Phases 0-6 + tactical fixes. |
| `tests/test_zupoly.c` | 38 | All ZUPoly primitives + Diophantine. Includes specific regression tests for the in-place `shift` direction bug and the subresultant-vs-pseudo-remainder edge cases. All passing. |
| `tests/test_bpoly.c` | 23 | All BPoly primitives plus compositional sanity checks (associativity, eval-after-shift, truncate-after-mul). All passing. |
| `tests/test_mvfactor.c` | 12 | Bivariate two-factor Hensel (6), multifactor pair-and-recurse (3), high-level orchestrator (3, with end-to-end via real `bz_factor_to_expr`).  All passing. |

A useful pattern for new test files: **override the `ASSERT` macro
to bypass `<assert.h>`**, because the CMake `Release` build defines
`-DNDEBUG` which compiles out the standard `assert`. Both
`test_zupoly.c` and `test_bpoly.c` do this; older test files predate
the awareness and may have silent-pass cases worth auditing
separately.

---

## 9. Open questions for the next session

1. **Diophantine primitive: Option A or B?** (Section 7.1.) Both are
   tractable; A reuses existing modular machinery, B is more direct.
2. **Multifactor strategy: balanced tree or linear?** Balanced is
   asymptotically better but more code. Linear matches the existing
   `multifactor_hensel_lift` style in `facpoly.c`.
3. **Should `factor_roots` be retired in one step or kept as a
   permanent fallback?** Argument for retiring: it has no
   correctness role once Hensel is in. Argument for keeping: it's
   effectively free when never reached, and provides a safety net
   for inputs the new pipeline rejects.
4. **n > 3 variates**: any prospect of generalising the Hensel
   pipeline directly, or commit to the bivariate-with-recursive-
   specialisation form?

---

## 9b. Phase 7: TrigFactor Path-B skip heuristic

After the Phase 0-6 work landed, TrigFactor's Path B (TrigExpand →
Factor → identity rules) accounted for 50-250 ms per call on small
inputs whose polynomial expansion just reinverted under factoring.
Phase 7 adds a structural predicate at the top of `builtin_trigfactor`:
when Path A is a no-op AND the input has no compound trig structure
(no `Power[trig, k>=2]`, no Times of multiple trig atoms), Path B is
skipped.

The predicate is conservative -- the user's primary case
`Sin[x]^3 + Sin[3 x] - 3 Sin[x]` has `Sin[x]^3` (a Power[trig, 3]) so
Path B still runs and produces the cancellation that makes the simp
candidate-set search converge.  Slow no-op cases like
`(Sin[x] + Sin[3 x])/(2 Cos[x] + 2 Cos[3 x])` correctly fall into the
skip branch.

Measured impact:

| Workload | Pre-Phase-7 | Post-Phase-7 | Δ |
|---|---|---|---|
| `Simplify[Sin[x]^3 + Sin[3x] - 3 Sin[x]]` | 478 ms | 395 ms | -17 % |
| `Simplify[Tan[2x] - 2 Tan[x]/(1 - Tan[x]^2)]` | 2.5 s | 2.05 s | -18 % |
| `Simplify[Cos[3x]/(4 Cos[x]^3 - 3 Cos[x])]` | 235 ms | 191 ms | -19 % |
| `Simplify[(Sin^2 + Cos^2)^4]` | 76 ms | 78 ms | flat |
| `Simplify[Sinh^2 - Cosh^2]` | 69 ms | 69 ms | flat |

Inside the user's primary case, the per-transform breakdown changes:

| Transform | Pre-Phase-7 | Post-Phase-7 |
|---|---|---|
| TrigFactor | 167 ms / 5 calls | 87 ms / 5 calls |
| TrigRoundtrip | 110 ms | 115 ms |
| (everything else) | ~190 ms | ~195 ms |

Same call count (the candidate-set search produces the same set);
each TrigFactor call's avg dropped from 33 ms to 17 ms.

## 10. Performance profile (after Phases 0-6 + tactical fixes)

For `Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]` the per-transform breakdown is:

| Transform | Total time | Calls | Avg per call |
|---|---|---|---|
| TrigFactor | 167 ms | 5 | 33 ms |
| TrigRoundtrip | 110 ms | 6 | 18 ms |
| TrigExpand | 42 ms | 5 | 8 ms |
| Factor (memoised) | 26 ms | 5 | 5 ms |
| FactorSquareFree | 18 ms | 5 | 4 ms |
| HalfAngle | 22 ms | 6 | 4 ms |
| TrigToExp | 16 ms | 5 | 3 ms |
| FactorTerms | 9 ms | 5 | 2 ms |
| PythagSquareComplete | 6 ms | 6 | 1 ms |
| PythagReduce | 34 ms | 66 | 0.5 ms |

Total transform time ≈ 450 ms (matches the observed 478 ms wall-clock).

**Where the remaining time goes.**  TrigFactor and TrigRoundtrip dominate.
Each of their slow invocations is on a small input (~10 leaves) that
Path A handles trivially but Path B then expands and re-factors,
typically returning the same form.  The Factor memo we added helps when
the underlying Factor calls coincide; it doesn't help when the surface
forms differ.

**Why further compression is hard.**  Mathematica completes the same
input in ~10 µs.  Closing that 50× gap would require either:

1. **Drastically fewer transform invocations** in the candidate-set
   search.  Could come from canonical-form keying (Together-then-expand
   the input before consulting the SimpMemo) or aggressive pruning
   based on a complexity-monotonicity assumption.  Architecturally
   non-trivial -- the candidate-set search currently relies on every
   transform getting a chance.

2. **Per-transform speedups** of 5-10×.  TrigFactor's Path B is the
   biggest target, but a heuristic that skips Path B safely on small
   inputs without angle-sum structure risks breaking the Phase-1
   identification of cancellations.  Worth investigating but requires
   careful correctness regression coverage.

Both options are multi-day projects beyond the scope of the current
factoring work.  The pipeline as it stands gives a 2.4× speedup on the
reference case, fixes the 2 facpoly regressions on `main`, and opens
the door to faster bivariate factoring.

## 11. Stopping point and resumption guide

The autonomous-loop sequence that produced this state ran for seven
turns and delivered:
- 9 phases / sub-phases of new functionality.
- ~3500 LOC of new source code (zupoly, bpoly, mvfactor, plus changes
  to facpoly and simp).
- ~2000 LOC of new tests across 5 test files (~95 cases, all passing).
- 4 → 2 main-branch test failures (the 2 facpoly cases were resolved
  as a side-effect of Phase 5d's tactical work).
- 2.4× speedup on the user-reported `Simplify` case.

Remaining work (in approximate order of effort):

| Item | Estimated | Value |
|---|---|---|
| Phase 5b: Wang's lc correction (non-monic bivariate Hensel) | 300-500 LOC | Coverage extension; few user-visible cases. |
| Phase 3b: Zassenhaus recombination | 200-300 LOC | Required for bivariate completeness; rare in practice. |
| Phase 4: n-variate recursion | 150-200 LOC | Mostly redundant with `factor_roots` for now. |
| Cleanup C2: smarter prime selector in `factor_zassenhaus` | 50 LOC | Robustness, not perf. |
| Cleanup C3: balanced Hensel product tree | 100 LOC | Marginal univariate speedup. |
| Cleanup C4: bound TrigRoundtrip on factored Tan denominators | unknown | Would let us drop the context-aware-scope hack in 5d. |
| Simplify-level candidate-set canonicalisation | unknown but architectural | The biggest potential speedup but a much bigger project. |

**Recommended next direction**, if the user wants to keep pushing on
performance: investigate the Simplify-level candidate-set search.
Profile a few canonical workloads to see whether canonical-form keying
or aggressive pruning would give more leverage than further factoring
work.  The factoring pipeline is no longer the main bottleneck.

If the goal is a more-complete Factor, Phase 5b (Wang) is the natural
next step.  But there's no specific user-reported case that demands it.

## 11b. Phases 8-17: Simplify-level memo work

After section 11's recommendation to investigate Simplify-level
caching, a series of memoisation phases landed.  All reuse the
existing `FactorMemo` (which is per-Simplify-call, lifecycle managed
in `builtin_simplify`) keyed under a per-transform head so the entries
never collide.

| Phase | What | Status |
|---|---|---|
| 8 | TrigRoundtrip memo | Done |
| 9 | TrigFactor + TrigExpand memos | Done |
| 10 | PythagReduce / PythagSquareComplete / HalfAngle memos | Done |
| 11 | Generic memo at `traced_call_unary` layer | **Reverted** — overhead exceeded gain on cheap transforms |
| 12 | Canonical-form keying for trig memos | **Reverted** — broke pattern-rule transforms (e.g., PythagReduce's `-1 + Cos[x]^2 + r___ :> -Sin[x]^2 + r` requires the literal `-1` term) |
| 13 | (skipped) | |
| 14 | TrigToExp-output keying for TrigRoundtrip | Done |
| 15 | TrigToExp keying for TrigFactor / TrigExpand | **Reverted** — pure overhead exceeded gains |
| 16 | TrigRoundtrip explosion guard (skip pipeline if TrigToExp expands input >5×) | Done |
| 17 | Two-level TrigRoundtrip memo (raw key + canonical key) | Done |
| 18 | TrigFactor leaf-explosion guard | **Reverted** — discarding Path A's expanded result changed candidate-set winners and produced an uglier final answer |
| 19 | TrigToExp raw-input memo | Done |
| 20 | PythagReduce / PythagSquareComplete / HalfAngle head-presence fast-skip | Done |

### Phase 16 — TrigRoundtrip explosion guard

`Cos[x] Cos[y]` (complexity 7) maps under TrigToExp to a sum of four
exponentials (complexity 77 — 11× growth).  Together / Cancel /
ExpToTrig on that intermediate is expensive *and* the final result
tends to use Cosh / Sinh of imaginary arguments, leaving us worse off
than the input.  The guard at `simp.c:1209` skips the pipeline when
`exp_score > 5 * in_score`.

Verified safe on the user-reference case (Sin³+Sin[3x]-3Sin expands
~3×, still benefits).  Triggers on inputs like products of trig atoms
on different variables.

History note: Phase 16 was initially reverted on a perceived
`integrals_tests` regression which turned out to be unrelated
(pre-existing failures from concurrent commits on `main`).  Re-applied
after user clarification.

### Phase 17 — Two-level TrigRoundtrip memo

Profiling the Tan double-angle case (`Tan[x] - 2 Tan[2x](1 - Tan[x]^2)/
(2(1 - Tan[x]^2)^2 - 4 Tan[x]^2)`) revealed that on TrigRoundtrip
*misses*, the cost was dominated by the canonical-key construction
itself: TrigToExp on Tan-rich subexpressions takes 80-130 ms before
the lookup even runs.  With Phase 16's guard active, the rest of the
pipeline is skipped — but the TrigToExp cost is still paid.

Phase 17 introduces a Level-1 lookup keyed on the raw input `e`,
*before* TrigToExp.  Hits short-circuit the entire pipeline, including
TrigToExp.  Misses fall through to a Level-2 lookup keyed on
TrigToExp(e), as before.  Results are stored under both keys, so:

- Future calls with identical `e` hit Level 1 (cheap).
- Future calls with canonically-equivalent `e'` hit Level 2 (saves
  the rest of the pipeline but pays TrigToExp(e')).

Measured impact:

| Workload | Pre-Phase-16/17 | Post-16 only | Post-16+17 |
|---|---|---|---|
| `Sin[x]^3 + Sin[3 x] - 3 Sin[x]` (user case) | 290 ms | 254 ms | 230 ms |
| Tan double-angle case | 3.55 s | 3.55 s (guard didn't help) | 2.88 s |

Net: 22% improvement on user case, 19% on Tan case from the combined
16+17 work over the baseline that already had Phases 8-10 + 14.

The Level-1 key shape is `TrigRoundtrip[e]`, Level-2 is
`TrigRoundtrip[TrigToExp(e)]` — both share the `TrigRoundtrip`
head so they don't collide with `Factor` / `TrigFactor` / etc. memo
entries, but they DO live in the same FactorMemo.

### Phase 18 — TrigFactor leaf-explosion guard (REVERTED)

Hypothesised that TrigFactor's Path A pipeline output, when more
than 2× the leaf count of the input, is dead weight: Simplify ranks
candidates by complexity and won't pick it.  Implemented as a
size check at the Path A return.

In practice the Tan double-angle test case revealed a subtler
interaction: the over-large Path A output is sometimes used as a
*bridge* by the Simplify candidate-set search — a different
candidate then transforms it into a winner.  Discarding it caused
the search to converge on a different (uglier) answer:
`Tan[x] - 2 Tan[2x](1 - Tan[x]^2)/(2 - 8 Tan[x]^2)` →
`2/5 Cos[x] Csc[x + ArcTan[2, -1]] Csc[x + ArcTan[2, 1]] Sin[x] + Tan[x]`.

Lesson: per-transform output filtering is unsafe when transforms
feed each other through the candidate set.  Output filtering must
be done at the candidate-set ranking layer, not inside individual
transforms.

### Phase 19 — TrigToExp raw-input memo

TrigToExp accounted for 472 ms / 15 calls / 31 ms avg in the Tan
double-angle profile post-Phase-17.  It's a pure function
(ReplaceAll[trig→exp rules] then Expand) so caching is safe.  It
shows up in two contexts within a Simplify call:

- As a standalone candidate-set transform.
- Internally inside `transform_trig_roundtrip`'s miss path (when
  Level-2 lookup needs the canonical TrigToExp(e) as its key).

A raw-input memo at the top of `builtin_trigtoexp` lets these
contexts share results.  The key shape is `TrigToExp[e]`, distinct
from all other memo heads.

Measured impact (combined with Phases 8-10, 14, 16, 17):

| Workload | Pre-Phase-19 | Post-Phase-19 | Δ |
|---|---|---|---|
| `Sin[x]^3 + Sin[3 x] - 3 Sin[x]` | 230 ms | 215 ms | -7 % |
| Tan double-angle | 2.85 s | 2.49 s | -13 % |

Cumulative speedups vs pristine baseline:

| Workload | Pristine | Post-Phase-19 | Speedup |
|---|---|---|---|
| `Sin[x]^3 + Sin[3 x] - 3 Sin[x]` | 1158 ms | 215 ms | 5.4× |
| Tan double-angle | 3300 ms | 2490 ms | 1.3× |

The user reference case is now at ~5× the pristine, with answer
quality preserved.  Tan case improvements are smaller because most
of its time is in TrigFactor's Path A pipeline, where output
filtering proved unsafe.

### Phase 3b — Zassenhaus recombination

The original `bpoly_hensel_lift_multi` (Phase 3a) used a pair-and-
recurse strategy: peel off `us[0]` against `prod(us[1..r-1])` via the
two-factor lift, then recurse on the residual.  This works only when
the bivariate factorisation has *exactly* the same number of factors
as the univariate image.  When the bivariate factorisation is COARSER
(some univariate factors must be combined into one bivariate factor),
the lift_2 attempt with a non-true factor fails verification.

Phase 3b enumerates k-subsets of the univariate factors for k = 1, 2,
..., r/2.  For each subset S we attempt `lift_2(P, prod(us in S),
prod(us not in S))`; the first successful lift yields a true bivariate
factor U_S and a residual V_C.  We then recurse on V_C with the
complement-side u's.

Key implementation details:
- **Subset enumeration**: Gosper's hack iterates k-combinations
  efficiently in lexicographic order over a uint64_t bitmask.
- **Symmetry break**: for k = r/2 with even r, the partition {S, S^c}
  is the same as {S^c, S}, so we restrict to subsets with bit 0 set.
- **Termination**: when no subset of any size yields a valid lift, P
  is irreducible bivariately under this set of u's; return P as the
  sole factor.
- **API change**: `bpoly_hensel_lift_multi` now takes an additional
  `int* r_out` parameter and returns variable factor count (≤ r).

Worst-case cost is O(2^r) two-factor lifts, but the common case (k=1
suffices) keeps it linear.  For squarefree images with r ≤ 10 the
worst case is tractable.

End-to-end results via `Factor[]`:
- `Factor[x^4 - y^2]` → `(x^2 + y)(x^2 - y)` (univariate at y=1 gives
  3 factors; recombination groups two of them).
- `Factor[x^4 - 4y^2 + 4y - 1]` → `(x^2 + 1 - 2y)(x^2 - 1 + 2y)`.

### Phase 4 — n-variate via specialise-and-trial-divide (minimum viable)

True multivariate Hensel lifting for n ≥ 3 would require either an
MPoly type with full polynomial-coefficient operations, or a nested
BPoly-of-BPolys extension.  Both are significant implementations.

The minimum-viable Phase 4 implemented here uses a different
strategy: **specialise one variable to 0, recursively factor the
result in n-1 variables, then trial-divide the original P by each
candidate factor**.  Factors that divide exactly are genuine factors
of P that happen to be independent of the specialised variable.

Algorithm in `factor_via_z_independent_split` (facpoly.c):
1. For each variable z in the input, compute P|z=0.
2. Recursively run `heuristic_factor` on P|z=0 (now in n-1 variables).
3. For each non-trivial factor f of P|z=0, run `exact_poly_div(P, f, ...)`.
   On success, peel f from the residual and continue.
4. After exhausting one specialisation, recurse on the residual.
5. If no specialisation yields a peelable factor, fall through.

End-to-end results via `Factor[]`:
- `Factor[(x+y)(x+z)]` → `(x+y)(x+z)` (catches both: x+y is
  z-independent, x+z is y-independent — the search finds at least
  one specialisation that exposes a clean factor).
- `Factor[(x+y)(x+y+z)]` → `(x+y)(x+y+z)`.
- `Factor[(x^2+y^2)(x+z)]` → `(x+z)(x^2+y^2)`.
- `Factor[(x+2y)(x-2y)(x+3z)]` → `(x+2y)(-x+2y)(-x-3z)` (signs
  rearranged in canonical form, but the same factorisation).
- `Factor[x^2+y^2+z^2+1]` → unchanged (irreducible, correctly).

**Limitations** (documented in test_factor_recombine.c):
- Cases where ALL factors depend on EVERY variable are not caught.
  E.g., `x^3 + y^3 + z^3 - 3xyz = (x+y+z)(x^2+y^2+z^2-xy-yz-xz)`
  has neither factor variable-independent; Phase 4 returns NULL and
  the legacy pipeline does not factor it either.  Handling such
  cases needs full multivariate Hensel.
- The trial-division step can be expensive for large polynomials
  with many candidate factors.  For practical inputs the recursion
  depth and candidate count stay small.

### Phase 20 — Pythag/HalfAngle head-presence fast-skip

Profiling the Tan case post-Phase-19 showed PythagReduce burning
~500 ms across 61 calls -- including 60-120 ms calls on huge
sum-of-exponentials expressions (TrigToExp output) where the rules
*cannot possibly match* (every rule LHS contains a Cos/Sin/Cosh/Sinh
head, and the inputs are pure E^(I x) forms).

`has_pythag_head(e)` does a single tree walk to detect any of the
four heads.  When absent, the ReplaceRepeated pass is skipped and
a copy of the input is returned.  Applied to all three transforms
(PythagReduce, PythagSquareComplete, HalfAngle) since their rule
sets share the requirement.

Measured impact:

| Workload | Pre-Phase-20 | Post-Phase-20 | Δ |
|---|---|---|---|
| `Sin[x]^3 + Sin[3 x] - 3 Sin[x]` | 215 ms | 213 ms | flat |
| Tan double-angle | 2.49 s | 2.00 s | -20 % |

The user case is unaffected because its candidates rarely hit
exponential form (Sin/Cos/Tan stay rational throughout).  The Tan
case dramatically benefits because TrigToExp blows it up into
exponential form for several candidates, and PythagReduce was
running pattern matching on those each iteration.

Cumulative speedups vs pristine baseline (after all phases):

| Workload | Pristine | Final | Speedup |
|---|---|---|---|
| `Sin[x]^3 + Sin[3 x] - 3 Sin[x]` | 1158 ms | 213 ms | 5.4× |
| Tan double-angle | 3300 ms | 2000 ms | 1.65× |

### Remaining bottleneck

TrigFactor is now the dominant cost (~860 ms / 15 calls / 57 ms avg
in the Tan case).  Its Path A pipeline (to-sincos rewrite → Together
→ Factor → identity rules → from-sincos) is intrinsically expensive
on Tan-rich inputs, and Phase 18's attempt at output filtering
proved unsafe (the over-large outputs are used as bridges by the
candidate-set search).

Future TrigFactor compression would likely require restructuring
the candidate-set search itself, or introducing pre-conversion
heuristics that detect when the Tan-form is already canonical and
skip the round-trip.  Both are larger projects.

## 12. Roadmap to retire `factor_roots`

User testing of `Factor[]` against Mathematica surfaced two classes of
input where the structured pipeline fails (returns NULL or hangs) and
falls through to `factor_roots` -- the legacy linear-trial-division
heuristic that catches only `(var - c)` shapes.  Mathematica completes
these in 0.1-1 ms; picocas takes 100 ms - 6 s.

The reference user-failure cases (and the gap each exposes):

| Case | Picocas current | Mathematica | Gap |
|---|---|---|---|
| `Factor[3 - 3x² + 4xy - 4x³y - 4y² + x²y² - 4xy³ + y⁴]` | **5 ms** ✓ (after the v_count==2 fast-path commit on the previous turn) | 0.15 ms | 30× — algorithmic constant factor |
| `Factor[(zx - x² - y²)(3z + 4xy - y²)]` (expanded) | **0.33 s** ✓ (after F2 MVP, 2026-05-02; previously ~6 s returning input unchanged) | < 1 ms | algorithmic constant factor only |
| `Factor[(1 - x¹²)(x - y¹³)(y - z¹⁴)]` (expanded) | **84 ms** ✓ (after the `univariate_squarefree`→`zupoly_gcd` fix; previously hung because F4's pre-check ran a univariate gcd over a degree-30 image through Knuth-style primitive PRS, which exhibits exponential coefficient growth) | few ms | algorithmic constant factor only |
| `Factor[x²(1 - x¹²)(1 + x - y¹³)(1 - y - z¹⁴)]` (expanded) | **167 ms** ✓ (same fix) | few ms | algorithmic constant factor only |
| `Factor[Expand[x²(z¹³-x¹²)(z⁴+3x⁹-y¹³)(17-5y-z¹⁴)]]` (user-reported 2026-04) | **0.9 s** ✓ (same fix; previously hung at >120 s in F4 pre-check) | 2 ms | algorithmic constant factor only |
| `Factor[Expand[(x²⁰⁰+y²⁰⁰+1)(x²⁰⁰-y+1)]]` (Fateman benchmark, user-reported 2026-05-02) | structurally similar (1.28 s on `k=70`, mathematically correct factorisation) ✓ — see §12.F6 for the trace and Tier 1+2 fix that landed | < 1 ms | univariate BZ substrate is still int64-bound; full bigint univariate BZ (Tier 3) not yet implemented |

The structured pipeline's coverage gap traces to five concrete items:

### Phase F1 — Wang's leading-coefficient correction (replaces 5b)

Currently `pick_monic_variable_index` requires `lc_x(P) == +1` to
enter the bivariate Hensel.  Inputs like `(1 - x¹²)(x - y¹³)` (LC =
-1 in x) and `(zx - x² - y²)(3z + 4xy - y²)` (LC = -1 or non-constant
in every variable) bypass the structured path.

#### Stage 1 — LC = -1 negate path (DONE)

The simplest case (LC = -1 with constant LC) is now handled by
pre-negating P before the lift.  After lifting the monic -P, the
overall sign is absorbed into the highest-x-degree factor via
`Expand[Times[-1, factor]]`, which distributes the -1 into the
factor's Plus terms.  This produces the same canonical output the
prior pipeline did for the regression case `Factor[3 a^2 b - 3 b - b^3]`
→ `b · (-3 + 3 a^2 - b^2)`: the highest-deg factor `b^2 - 3a^2 + 3`
gets negated to `-3 + 3a^2 - b^2`, leaving the linear `b` factor with
positive sign.

The previous "naive LC=-1 hack" mentioned below was reverted because
it absorbed the sign by setting the WHOLE product as `Times[-1, b, ...]`
without distributing into a Plus, leaving a leading `Times[-1, ...]`
that didn't match the canonical baseline.  The fix is the explicit
`Expand` after the multiplicative negation.

Wired in `factor_bivariate_via_hensel` (`src/facpoly.c`).  Heuristic
chooses the largest-x-degree factor (tiebreak: largest y-degree) as
the absorber, which empirically matches the canonical printed form
for the regression case.

Tested at small/medium degrees `(1 - x^k)(x - y^m)` for k ∈ {2, 4, 6,
8, 10}, all under 12 ms.  k = 12 still hangs due to F3 territory
(recombination performance on 7-factor images), independent of F1
since `(1 + x^12)(x - y^13)` (no negate) hangs identically.

#### Stage 2 — constant `|a| > 1` LC (DONE)

For inputs with constant non-unit integer LC (e.g. `lc_x(P) = 6`),
the lift's monic-in-x precondition is met by Wang's monic-substitution
recipe:

```
Given P(x, y) ∈ Z[x, y] of x-degree d with lc_x(P) = a, |a| > 1, a ∈ Z:
  1. Form  Q(x, y) := a^(d-1) · P(x/a, y).
     Concretely, cx_Q[d] = 1 and cx_Q[j] = cx_P[j] · a^(d-1-j)
     for 0 ≤ j < d.  Q is monic in x with integer coefficients.
  2. Factor Q via the existing monic Hensel pipeline:
        Q = G_1 · ... · G_r,  each G_i monic in x.
  3. Recover true factors of P:
        F_i := G_i(a·x, y) / cont_Z(G_i(a·x, y))
     where cont_Z denotes the integer content (gcd of all integer
     coefficients across both x and y).
```

Why step 3 works: writing `H_i := G_i(a·x, y)`, we have
`prod H_i = a^d · P(x, y) / a = a^(d-1) · P(x, y)`.  Each H_i is
G_i(x, y) with the x^j coefficient scaled by `a^j`, so
`cont_Z(H_i) = a^(d_i) · gcd(... a-related factors ...)` — exactly
the share of `a^(d-1)` redistributed into G_i.  Summing in log-space:
`prod cont_Z(H_i) = a^(d-1)`.  Dividing each H_i by its content
gives F_i with `prod F_i = P`.

For example, `Factor[Expand[(2x+3y)(3x+5y)] = 6x² + 19xy + 15y²]`:
- a = 6, d = 2.
- Q = 6 · P(x/6, y) = x² + 19xy + 90y² = (x + 9y)(x + 10y).
- G_1 = x + 9y, G_2 = x + 10y.
- H_1 = G_1(6x, y) = 6x + 9y, cont_Z(H_1) = gcd(6, 9) = 3, F_1 = 2x + 3y.
- H_2 = G_2(6x, y) = 6x + 10y, cont_Z(H_2) = gcd(6, 10) = 2, F_2 = 3x + 5y.
- prod F_i = (2x + 3y)(3x + 5y) = 6x² + 19xy + 15y² = P. ✓

Stage 1 (negate) composes with Stage 2: when `lc_x(P) = -k` for
`k > 1`, the picker chooses Stage 2 with `|a| = k` and sets the
negate flag, which pre-negates P before the Wang substitution and
flips the sign of one factor at the end.

Wired in `factor_bivariate_via_hensel` (`src/facpoly.c`).  The new
helpers (~160 LOC) are file-static:
- `bpoly_int_content`, `bpoly_div_int_exact` — generic BPoly content
  manipulation (could be promoted into bpoly.c if other callers need
  them later).
- `bpoly_make_monic_via_x_scale` — Q = a^(d-1) · P(x/a, y).
- `bpoly_subst_x_scale` — H = G(c·x, y).
- `pick_factorable_x_var` — tri-state picker (MONIC / NEGATE / SCALE).

Tested by `test_factor_recombine.c::test_f1_scale_*` covering:
- Two-factor positive LC (`(2x+3y)(3x+5y)`).
- Two-factor with one negative coefficient (`(2x-3y)(3x+5y)`).
- Stage 1 + Stage 2 composition (negative LC, `|a| > 1`).
- Three-factor mixed monic + non-monic.
- Coefficient-symmetric inputs in `(a, b)`.
- Stage-2-shaped irreducible (`2x² + 3xy + 7y²`) — pass-through.

#### Stage 3 — polynomial-in-y LC (NOT STARTED)

For inputs whose LC is a polynomial in y (`lc_x(P) = A(y)` not
constant), the proper Wang algorithm (Geddes-Czapor-Labahn §6.6) is
needed:

```
Given P(x, y) ∈ Z[x, y] of x-degree d with lc_x(P) = A(y):
  1. Substitute α and factor P(x, α) univariately:
     P(x, α) = c · u_1(x) · ... · u_r(x).
  2. Adjust each u_i so its leading coefficient (in x) divides A(α):
     u_i := u_i · (A(α) / lc_x(u_i)).  This pre-distributes A
     across the lifted factors.
  3. Lift  P̃ = A^(r-1) · P  -- now monic in x with leading
     coefficient A^r.
  4. After Hensel lifting, divide each U_i by the appropriate
     "share" of A so that Π U_i = P (not P̃).
```

The standard recipe: track A's content factorisation through the
lift; each factor "gets back" the right share via primitive-part
extraction at each step.  The book's reference implementation is
~120 lines; with our existing BPoly / ZUPoly substrate the
picocas version is ~300 LOC, plus ~150 LOC of tests.

**Unblocks** (Stage 3): bivariate inputs with polynomial-in-y LC,
e.g. `(zx - x² - y²)(3z + 4xy - y²)` (every variable's LC depends
on the others); also overlaps with F2's multivariate cases.

**Estimated** (Stage 3): 300-500 LOC implementation + 150-200 LOC tests.

### Phase F2 — True multivariate Hensel for n ≥ 3 (MVP DONE 2026-05-02)

The previous Phase 4 (`factor_via_z_independent_split`) only handled
inputs where at least one factor is independent of one variable.
Cases like `(zx - x² - y²)(3z + 4xy - y²)` -- where every factor
depends on every variable -- need actual multivariate Hensel
lifting.

The MVP implementation lands in three stages:

#### Stage F2.1 -- MPoly substrate (DONE)

`src/mpoly.{c,h}`: sparse Z[x_1, ..., x_n] with mpz_t coefficients,
lex descending term order, row-major exponent storage parallel to
coefficient array.  Operations: construction
(new/zero/from_int/from_mpz/monomial/copy), term manipulation
(push_term/set_coef/get_coef/normalize), arithmetic
(add/sub/neg/mul/scale), substitution (subst_var_int/shift_var_int/
coef_of_var/lc_var), Expr round-trip.  38 unit tests in
`tests/test_mpoly.c`.

#### Stage F2.2 -- Trivariate Hensel lift (DONE)

`src/mvfactor3.{c,h}`:

- **Bivariate Diophantine** (`mpoly_diophantine_2`): solve Δu·V +
  Δv·U = E over Z[var_main, var_y].  Algorithm: specialise var_y to
  a chosen integer α_y where U(x, α_y), V(x, α_y) remain coprime and
  degree-preserving in var_main; solve univariate Diophantine via
  zupoly_diophantine; lift in var_y via Hensel iteration up to
  deg_{var_y}(E).  At each iteration solve the residual at y^k via
  another zupoly_diophantine (over u(x, 0), v(x, 0)).  Tries up to
  11 alpha values before giving up.

- **Trivariate two-factor Hensel** (`mpoly_hensel_lift_3_2`): given
  P ∈ MPoly{n vars} and bivariate seeds U_xy, V_xy with var_z-degree
  0 such that U_xy · V_xy = P|var_z=α_z, lift to U, V ∈ MPoly{n vars}
  with U·V = P.  Shifts var_z by α_z, iterates degrees k = 1..B_z
  computing residual (P_sh - U·V)[z^k], solves bivariate Diophantine
  on the residual, updates U += z^k Δu, V += z^k Δv.  Verifies final
  product.

- **MPoly ↔ BPoly converters** (`mpoly_to_bpoly_in`,
  `bpoly_to_mpoly_in`) so the orchestrator can reuse the existing
  bivariate Hensel for the y-lift step.

9 unit tests in `tests/test_mvfactor3.c` cover MPoly ↔ ZUPoly
round-trip, the bivariate Diophantine solver, and the trivariate
two-factor lift.

#### Stage F2.3 -- heuristic_factor wiring (DONE)

`src/facpoly.c`: `factor_trivariate_via_mhensel` is invoked in
`heuristic_factor` between Phase 4 (specialise-and-divide) and the
is_likely_irreducible_multivariate gate, but only for v_count == 3.
The orchestrator:

1. Convert P to MPoly.
2. Try each variable as main, requiring deg ≥ 1 and lc = constant 1.
3. Try (α_y, α_z) ∈ {0, ±1, ±2, ±3}^2 = 49 alpha tuples.
4. For each: compute univariate image, check squarefree, factor via
   existing `factor_via_bz_callback`, require r = 2.
5. Specialise z = α_z, shift y by α_y, run `bpoly_hensel_lift_2`.
6. Convert bivariate factors to MPoly, unshift y, then run
   `mpoly_hensel_lift_3_2` to introduce z.
7. Verify product == P, convert back to Expr, return.

Drive-by fix: corrected `extract_monomial` to descend into nested
Times.  Without the fix, `factor_binomial` on inputs like
`Plus[Power[z,2], Times[-1, Times[Power[x,2], Power[y,2]]]]`
silently dropped the `x²y²` structure (the inner Times wasn't
recognised), causing `Factor[z² - x²y²]` to incorrectly return
`(z-1)(z+1)`.

#### Measured impact

| Workload | Pre-F2 | Post-F2 |
|---|---|---|
| `Factor[Expand[(zx - x² - y²)(3z + 4xy - y²)]]` | ~6 s, unfactored | **0.33 s, factored correctly** |
| `Factor[Expand[(z + xy)(z - xy)]]` | wrong (z-1)(z+1) | 15 ms, correct |
| `Factor[Expand[(z + x)(z + y)]]` | wrong | 12 ms, correct |
| `Factor[Expand[(z² + x)(z + y)]]` | wrong | 5 ms, correct |

#### MVP scope and follow-ups

What this MVP handles:
- n = 3 inputs, monic in some main variable.
- Two-factor outputs (r = 2).
- alpha tuples in a small budget (49 attempts).

What it doesn't yet handle:
- n ≥ 4 (extension is recursive but requires more glue).
- Non-monic-in-main inputs (need Wang's lc correction at the
  trivariate level, similar to F1's bivariate Wang).
- r ≥ 3 outputs (need multifactor Hensel + recombination at the
  trivariate level).
- Inputs where 49 alpha tuples don't yield a good seed.

These extensions land in subsequent phases.  The current MVP
unblocks the most user-visible case from §12.

Sketch:

```
Input: P ∈ Z[x_1, ..., x_n], n ≥ 3.
  1. Pick main variable x_1.
  2. Specialise (x_2, ..., x_n) → (α_2, ..., α_n) chosen so:
     - lc_{x_1}(P)(α_2,...,α_n) ≠ 0
     - P(x_1, α_2, ..., α_n) is squarefree univariately
  3. Factor the univariate image: P|α = u_1 · ... · u_r.
  4. Lift the factorisation across one variable at a time:
     - Lift x_2: gives bivariate factors in (x_1, x_2)
       (existing bpoly_hensel_lift_2/multi).
     - Lift x_3: now coefficients are bivariate; need a "BPoly[
       BPoly]" -- an MPoly substrate.
     - ... continue through x_n.
  5. At each step, run subset recombination if the extension's
     image-factorisation is finer than its preimage.
```

Two implementation choices:

- **MPoly type** (`Z[x_1, ..., x_n]`): cleanest abstraction; reuses
  ZUPoly arithmetic via term-keyed hash tables.  Most code in the
  CAS canon is written this way.  Estimate: 800-1500 LOC for the
  type plus 500-800 LOC for the trivariate-and-higher Hensel
  iteration.
- **Nested BPoly**: BPoly with BPoly coefficients, capped at n=3.
  Less general but reuses the existing BPoly substrate.  Estimate:
  600-1000 LOC for the n=3 wrapper.

**Recommended**: MPoly.  Caps make the user disappointed twice (now
for n=3, later for n=4) and the abstraction makes downstream work
(e.g., GCD over Z[x_1, ..., x_n]) much easier.

**Unblocks**: case 2 above and similar trivariate-or-higher inputs.

**Estimated**: 1000-1500 LOC implementation + 300-500 LOC tests.

### Phase F3 — Bivariate Hensel performance (DONE)

Even with F1 letting non-monic inputs enter the bivariate path, the
per-iteration polynomial-arithmetic cost in `bpoly_hensel_lift_2` was
too high on inputs with deep y-lifts (cumulative O(deg² · B³)) and the
recombination explored subsets in arbitrary order regardless of degree.
F3 lands three optimisations in `src/mvfactor.c`:

#### Stage F3a — sort `us[]` by degree before recombination (DONE)

`bpoly_hensel_lift_multi` now sorts a private copy of the borrowed
`us[]` array ascending by univariate x-degree before recursing into
`lift_multi_internal`.  The complement-side recursion preserves the
relative order (`comp_us` iterates cleared mask bits in increasing
index), so a single sort at the public entry point suffices.

For inputs whose true bivariate factor has a low-degree univariate
image -- the dominant pattern in user reports, e.g.
`(1 - x^k)(x - y^m)` whose `(x - y^m)` factor reduces to `x` at α=0 --
the singleton subset `{x}` is now the FIRST recombination trial
instead of being one of `r` possibilities.  The lift succeeds on the
first attempt instead of the seventh.

Helper: `mvfactor_sort_us_by_degree` (insertion sort, O(r²) for the
typical r ≤ 10).

#### Stage F3b — incremental U·V update (DONE)

The previous Hensel iteration recomputed `U*V` from scratch via
`bpoly_mul_truncate_y(U, V, k+1)`, paying O(d_x² · k²) per iteration
and O(d_x² · B³) cumulatively.  F3b maintains a running `UV` BPoly
that is `U_{k-1}·V_{k-1} mod y^{k+1}` at the start of iteration `k`:

```
At start of iteration k (UV = U_{k-1}·V_{k-1} mod y^k):
  1. Extend UV by the y^k cross coefficient
       cross_yk = sum_{a=1}^{k-1} U[y^a] * V[y^{k-a}]
     (k-1 zupoly_mul calls, each O(d_x²)).
  2. Compute E_k = P[y^k] - UV[y^k]; bail if zero.
  3. Solve diophantine for Δu, Δv.
  4. Update U += y^k·Δu, V += y^k·Δv.
  5. Update UV[y^k] += Δu·V_0 + U_0·Δv  (since V_{k-1}[y^0] = V_0,
     U_{k-1}[y^0] = U_0; the y^{2k} cross term Δu·Δv vanishes mod
     y^{k+1} for k ≥ 1).
```

Cumulative cost drops from O(d_x² · B³) to O(d_x² · B²) -- a factor
of B speedup for deep y-lifts.

The same scheme applies to `bpoly_hensel_lift_2_lc`; the only twist
is that `U_0`, `V_0` are the qu_0/qv_0-scaled univariate factors
(`u_init`, `v_init`) instead of just `u`, `v`.  We keep them alive
through the lift loop for the step-5 update, freeing them at the end
(and on every error-return path).

Helper: `bpoly_uv_cross_yk(U, V, k)`.

#### Stage F3c — Mignotte coefficient bound fast-fail (DONE)

Each lift call now precomputes a Mignotte-style integer-coefficient
bound

```
M = 2^{deg_x(P) + deg_y(P)} * ceil(||P||_2)
```

via `bpoly_mignotte_bound`, where `||P||_2` is the L2 norm of P's
integer coefficient vector.  After every iteration's update, U and V
are scanned via `bpoly_max_abs_coef_exceeds`; if either has any
coefficient exceeding M in absolute value, the subset cannot produce
an integer factor of P, so the lift aborts immediately instead of
running the remaining iterations on a doomed candidate.

The bound is conservative (the true Mignotte bound for multivariate
divisors is tighter), but cheap to compute and decisive on divergent
subsets where intermediate coefficients grow exponentially.

#### Tests

`tests/test_factor_recombine.c::test_f3_perf_*` exercises the lift
inner loop deeply:
- `(1 - x^6)(x - y^15)` -- 5 image factors, B = 15.
- `(1 - x^4)(x - y^25)` -- B = 25, the deepest lift in the suite.
- `(1 - x^8)(x - y^17)` -- 5 image factors, B = 17.

Each completes in single-digit milliseconds end-to-end.

#### Out of scope

The original "(1 - x^12)(x - y^13)" hang turned out NOT to be
bivariate-Hensel territory: tracing showed it hangs inside
`bz_factor_to_expr(x^13 - x)` -- the univariate Berlekamp-Zassenhaus
factorizer -- before the bivariate path even touches the input.
`Factor[x^13 - x]` and `Factor[x^13 - 1]` hang in isolation.  Fixing
that is a separate task from F3 (univariate BZ performance, not
bivariate Hensel).

**LOC**: ~180 LOC implementation + ~30 LOC tests.

### Phase F4 — Faster multivariate FactorSquareFree

Profiling case 4 shows `FactorSquareFree` itself takes 148 ms on
the 4-variable input -- it's the right answer, but the cost is high
relative to Mathematica's milliseconds.  Current implementation runs
multivariate `gcd(P, ∂P/∂x)` recursively, and multivariate GCD over
`Z[x_1, ..., x_n]` is O(deg² · n) per recursion.

Options:

```
1. Cheap squarefree pre-check: gcd(P, ∂P/∂x) in ONE variable.  If
   it's a unit, P is squarefree in x; no further work.  Most
   user-typed polynomials are squarefree.
2. EZ-GCD or Brown's modular GCD instead of subresultant PRS.
3. Cache squarefree-decomposition results in the FactorMemo.
```

**Estimated**: 200-300 LOC for the cheap pre-check (highest leverage),
500+ LOC for an EZ-GCD overhaul (deferred).

#### Phase F4 Stage 1 — cheap squarefree pre-check (DONE)

Implemented as `sqfree_cheap_check` in `src/facpoly.c`.  After the
content-extraction step in `factor_square_free_poly` produces the
primitive part `pp` (in main variable `x`), the helper substitutes
integer values from `{1, -1, 2, -2, 3, -3, 4}` for the other variables
to obtain a univariate image, and uses `univariate_squarefree`
(`gcd(image, image')` over `Z[x]`) as the actual squarefree probe.
A single squarefree image at an alpha that preserves the leading-x
degree is sufficient proof.

**Soundness.**  A repeated factor `f(x, y, ...)^2` of `pp` must
involve `x` non-trivially (a constant-in-`x` repeated factor would
have divided `content(pp, x)` and been pulled out).  At a generic
alpha, `f(x, alpha)^2` has positive degree in `x`, so the univariate
image is non-squarefree.  The bad alpha set is finite (roots of
`lc_x(f)`), so trying seven values catches all real-world cases.

**Limited to multivariate inputs (var_count >= 2).**  For univariate
input the multivariate `gcd(pp, pp')` over `Z[x]` is unavoidable in
the slow path anyway, so a separate pre-check would only pessimise
the non-squarefree case.

**Measured impact.**  4-variable input
`Expand[(x*y - 1)*(x + y + z)*(x*y*z + w)*(x - y + 1)*(w - z + y - x)]`:
6.27 s → 0.95 s (6.6× speedup).  3-variable Hensel-style cases
(60-150 ms) are essentially unchanged because the multivariate GCD
already terminates quickly on inputs whose leading-x coefficient is
a constant.  Non-squarefree inputs fall through to the original Yun
loop with no measurable overhead from the failed probe.

**Tests.**  `test_factorsquarefree_f4_fastpath` in
`tests/test_facpoly.c` covers (a) multivariate squarefree fast-path,
(b) multivariate non-squarefree with a single repeated factor,
(c) multivariate non-squarefree with two squared factors,
(d) three-variable squarefree fast-path, and (e) repeated factor
that is constant in the main variable (handled by the pre-existing
content-extraction recursion).

#### Phase F4 Stage 2 — `univariate_squarefree` via `zupoly_gcd` (DONE, 2026-04)

The cheap pre-check from Stage 1 substitutes integer values for the
non-main variables and asks "is the resulting univariate image
squarefree?" via a univariate `gcd(image, image')`.  Stage 1 wired
this to the existing `poly_gcd_internal` (Knuth-style primitive
PRS at the Expr level), which on degree-≥20 images suffered the
classical exponential coefficient growth on the pseudo-remainders --
each step multiplies through by `lc(B)^(deg(A)-deg(B)+1)` and the
content-stripping `poly_content` / `exact_poly_div` calls themselves
work on giant Expr trees.  Effect: a degree-31 image takes >55 s,
dominating the entire FactorSquareFree call (which is itself called
from many places: `Factor`, `Together`, `Cancel`, `Apart`, ...).

Stage 2 replaces the inner gcd with `zupoly_gcd` -- subresultant
PRS over GMP `mpz_t` coefficients.  Subresultant PRS keeps
intermediate coefficient sizes polynomially bounded, so the same
degree-31 image now finishes in sub-millisecond time.  The gcd
returns 0 deg iff the image is squarefree, which is exactly what
the pre-check needs.

End-to-end results via `Factor[]`:
- `Factor[Expand[(1 - x¹²)(x - y¹³)(y - z¹⁴)]]`: hang → 84 ms.
- `Factor[Expand[x²(1 - x¹²)(1 + x - y¹³)(1 - y - z¹⁴)]]`: hang → 167 ms.
- `Factor[Expand[x²(z¹³ - x¹²)(z⁴ + 3x⁹ - y¹³)(17 - 5y - z¹⁴)]]`
  (user-reported 2026-04): >120 s → 0.9 s.

The fallback path (`expr_to_zupoly` returns NULL) preserves the old
Expr-level gcd for inputs with rational coefficients that slip
through the integer-conversion layer.

**LOC**: ~30 LOC implementation in `univariate_squarefree`, plus a
regression test `test_factor_trivariate_high_degree_squarefree` in
`tests/test_factor_baseline.c`.

### Phase F5 — Recombination cap & heuristics (DONE)

The previous `lift_multi_internal` enumerated k-subsets exhaustively
through Gosper's hack.  For r = 7 factors with no recombination
applicable (the bivariate factorisation has the same number of
factors as the univariate), it tried all 63 subsets before
concluding.  F5 added two heuristic improvements:

#### Stage F5a — singleton ordering by degree (DONE as F3a)

Sort `us[]` ascending by univariate x-degree before recursion so the
lowest-degree singletons are tried first.  Already shipped as part of
F3a in `mvfactor_sort_us_by_degree`; the singleton {x} for inputs like
`(1 - x^k)(x - y^m)` is now the first recombination trial instead of
the seventh.

#### Stage F5b — partial-lift signal short-circuit (DONE)

The lift-2 inner helper (`bpoly_hensel_lift_2_internal` in
`src/mvfactor.c`) now exposes a `completed_iters_out` flag.  It is
true when the lift loop completed all B y-iterations without bailing
on Diophantine non-integer or Mignotte coefficient overflow -- in
that case the lift only failed at the final `U·V == P` verification,
which is the "near-success" signal that motivates trying recombination
at the next k.

`lift_multi_internal` aggregates the flag across all k-subsets at the
current k.  When an entire k-pass produces no completion signal AND
k ≥ 2, the loop short-circuits and returns P as irreducible
bivariately.  k = 1 is always exhausted unconditionally so legitimate
singleton lifts (e.g. `{x²+1}` for `x⁴ - y²`) still fire, and k = 2 is
always tried before the short-circuit fires (so all current k = 2
recombination cases remain correct).

For r = 8 inputs the savings are largest:
- Pre-F5: 8 + 28 + 56 + 35 = 127 lifts on a fully-divergent input.
- Post-F5: at most 8 + 28 = 36 lifts (k = 1 + k = 2), then break.

Measured impact:

| Workload | Pre-F5 | Post-F5 | Δ |
|---|---|---|---|
| `Factor[Expand[x¹²+y²-1]]` (r=6) | 23 ms | 20 ms | -13 % |
| `Factor[Expand[x¹²+y⁴-1]]` (r=6, B=4) | 24 ms | 20 ms | -17 % |
| `Factor[Expand[x²⁰+y²-1]]` (r=6) | 31 ms | 28 ms | -10 % |
| `Factor[Expand[x²⁴+y²-1]]` (r=8) | 117 ms | 84 ms | -29 % |
| `Factor[Expand[x³⁰+y²-1]]` (r=8) | -- | 59 ms | (new) |

Existing recombination tests (`(x²-y)(x²-4y)`, `x⁴-y²`,
`x⁴-4y²+4y-1`, `(x²-y)(x²+y)`, `(x³+y)(x³-y)`) all pass unchanged --
they either succeed at k = 1 (top level, with recursion handling the
residual) or hit a real k = 2 partition before any short-circuit
fires.

**LOC**: ~30 LOC implementation in `src/mvfactor.c`, ~100 LOC of
tests in `tests/test_factor_recombine.c::test_f5_*`.

### Phase 5c — Retire `factor_roots`

After F1 + F2 + F3 land, the structured pipeline covers:
- Univariate (existing `bz_factor_to_expr`).
- Bivariate monic (Phase 5a, done).
- Bivariate non-monic (F1 / Wang).
- Trivariate-and-higher with at least one variable-independent
  factor (Phase 4 minimum-viable, done).
- Trivariate-and-higher with all-variable-dependent factors (F2
  multivariate Hensel).

At that point `factor_roots` is dead code.  Remove the function and
its call site in `heuristic_factor`.  Update the comment block in
`builtin_factor` to remove the historical references.  Run the full
test suite to verify no regressions.

**Estimated**: 50 LOC of deletion + the full test-suite pass.

### Phase F6 — Bigint-coefficient handling in the Factor pipeline (Tier 1+2 DONE 2026-05-02; Tier 3+4 PENDING)

The Fateman benchmark `Factor[Expand[(x^200 + y^200 + 1)(x^200 - y + 1)]]`
*used to* return the input unchanged in 0.33 s — a wrong answer dressed
up as a fast one. The pipeline traces as follows (verified by instrumenting
the live build on 2026-05-02):

```
heuristic_factor(P)                                           [v_count = 2]
  -> factor_bivariate_via_hensel(P, vars=[x, y])
       -> pick_factorable_x_var: main_idx=0 (x), kind=MONIC
          (lc_x(P) = 1, the constant from the lone x^400 term)
       -> mvfactor_try_bivariate_monic(P, factor_via_bz_callback)
            alpha = 0:  image = (1 + x^200)^2          NOT squarefree, skip
            alpha = 1:  image = x^200 (2 + x^200)      NOT squarefree, skip
            alpha = -1: image = (x^200 + 2)^2          NOT squarefree, skip
            alpha = 2:  image = (x^200 - 1)(x^200 + 2^200 + 1)
                        squarefree, deg = 400, calling factor_cb
              factor_via_bz_callback
                -> bz_factor_to_expr(image)
                   coefficient `2^200 + 1` is EXPR_BIGINT
                   non_integer_seen = 1   (facpoly.c:3725-3736)
                   *** BAIL: returns expr_copy(image), unchanged ***
                walk result: input is a Plus, not a Times
                -> wraps the unchanged Plus into a single ZUPoly, r = 1
              orchestrator interprets r == 1 as
                  "image is irreducible univariately"
                  -> "P is irreducible bivariately"
              returns true with one factor = bpoly_copy(P)
       sees r <= 1, returns NULL
  -> is_likely_irreducible_multivariate(P, vars=[x, y])
       alpha = 1, -1: image not squarefree, skip
       alpha = 2:  squarefree image, BZ bails (BIGINT), n = 1 → confirmation #1
       alpha = -2: squarefree image, BZ bails (BIGINT), n = 1 → confirmation #2
       returns TRUE
  -> heuristic_factor returns expr_copy(P)
```

**Root cause** (single point):

`bz_factor_to_expr` (`src/facpoly.c:3715`) at lines 3723-3742 only
handles coefficients of type `EXPR_INTEGER`. For inputs whose
coefficients have promoted to `EXPR_BIGINT` (because the magnitude
exceeded `int64_t`), it sets `non_integer_seen = 1` — originally
intended for symbolic / Gaussian / rational coefficients — and
returns `expr_copy(P)`. The underlying `UPoly` substrate that
`factor_zassenhaus` operates on stores coefficients as `int64_t`,
so widening this single ingestion point is not by itself a fix.

**Why the current call sites can't tell what happened:**

`bz_factor_to_expr`'s "I couldn't represent the input" signal is
`return expr_copy(P)` — the same shape it would return for a genuine
irreducible. Both downstream consumers take the bait:

- `factor_via_bz_callback` (`facpoly.c:1204`) walks the returned Expr
  looking for `Times[...]`. A bare Plus is treated as one factor
  (`r = 1`).
- `count_nontrivial_factors` (`facpoly.c:1003`) does the same: a
  non-Times input with positive degree is one factor.

Since both call sites use `r == 1` / `n == 1` as the "irreducible
univariate image" confirmation, a coefficient overflow becomes a
false irreducibility verdict.

**Why this is specifically a Fateman issue:**

The polynomial itself has small int64 coefficients (1, 2, ±1).
The overflow is created *by the orchestrator's evaluation step*:
`P(x, alpha)` with `|alpha| >= 2` introduces `alpha^k` for `k`
up to deg_y(P). For Fateman, deg_y = 201 and the smallest alpha
that yields a squarefree image is `alpha = 2`, producing `2^200 + 1`
(a 60-digit BIGINT). The same blow-up affects every alpha in
`{2, -2, 3, -3, 4, ...}` — there is no usable evaluation point
within the current alpha set that keeps coefficients in int64.

Empirically the cliff is at `k ≈ 63` (where `2^k` first overflows
`int64_t`):

| `(x^k + y^k + 1)(x^k - y + 1)` expanded | Factor result |
|---|---|
| k = 20 | factored, 56 ms |
| k = 40 | factored, 41 s (BZ on degree-80 with int64 coeffs is feasible but slow) |
| k = 100 | unchanged, 110 ms |
| k = 200 | unchanged, 330 ms |

The k=40 → k=100 jump is the int64 cliff: at k=63 a single coefficient
of the alpha=2 image trips `non_integer_seen` and the diagnosis
collapses to "irreducible".

**Fix tiers** (in increasing scope):

1. **Stop lying about it.** Have `bz_factor_to_expr` return `NULL`
   (or a distinguished "couldn't ingest" sentinel) when
   `non_integer_seen` is set, and update the two call sites to
   treat that as "skip this alpha", not "one factor". This does
   NOT factor Fateman, but it stops the false irreducibility
   verdict so the polynomial falls through to `factor_roots` /
   "P unchanged" with the correct semantic. Cheapest defensive
   fix; also helps any other input whose evaluation overflows.
   Estimated: ~50 LOC.

2. **Widen the ingestion path.** `bz_factor_to_expr` accepts
   `EXPR_BIGINT` coefficients by reading them into a `ZUPoly`
   instead of a `UPoly`, and `factor_zassenhaus` learns to operate
   on ZUPoly-shaped inputs via reduction modulo `p` for the
   distinct-degree / equal-degree finite-field steps. The
   modular Hensel lift still has the int64 modulus ceiling
   (`p^k <= 2^62`), but for inputs whose Mignotte bound on the
   factors fits in that ceiling (typical small-coefficient cases
   like Fateman, where the true factors `x^200 + y^200 + 1` and
   `x^200 - y + 1` have unit-bounded coefficients) this would be
   enough. Estimated: 200-400 LOC.

3. **Arbitrary-precision Berlekamp-Zassenhaus.** Replace `UPoly`
   with `ZUPoly` throughout `factor_zassenhaus`, including the
   modular arithmetic (`mpz_mod` rather than `% p`) and the
   Hensel lift (lift to `p^k > 2 * Mignotte`, no int64 ceiling).
   This is the only fix that handles univariate inputs whose own
   coefficients are large (e.g. cyclotomic-like polynomials with
   coefficients 2^k for large k). Estimated: 600-1000 LOC plus
   regression coverage on the existing BZ test suite.

4. **Sparse multivariate Hensel without numerical evaluation.**
   The reason Mathematica handles Fateman in microseconds is that
   it does not go through `P(x, alpha)`-and-lift at all: it uses
   sparse Hensel iteration directly on the bivariate, never
   instantiating the alpha-induced coefficient explosion. Picocas
   would need a sparse Z[x, y] type that stores the polynomial as
   a list of `(monomial, coefficient)` pairs and a Hensel
   iteration that operates on it directly. The MPoly type
   (`src/mpoly.{c,h}`) added in F2 is most of the substrate;
   what's missing is a sparse-aware Hensel iteration. Estimated:
   1500-2500 LOC; this is the only path that takes Fateman from
   "0.33 s wrong" to "1 ms right".

**Recommended sequencing.** Tier 1 alone is worth landing
immediately as a correctness fix — picocas currently lies about
Fateman being irreducible, which is worse than failing to factor
it. Tier 2 is the natural follow-up if there are other inputs
whose evaluation-point coefficients exceed int64 but whose true
factors are small (cheap BZ over int64 modulus, just need bigint
input). Tier 3 is necessary for the full original `Factor[x^k - 1]`
class once `k` exceeds ~62. Tier 4 is the only way Fateman itself
gets fast.

#### What landed 2026-05-02 (Tier 1 + Tier 2 hybrid)

Three changes in `src/facpoly.c`:

- `factor_via_bz_callback` (the bivariate orchestrator's
  univariate-image factoring callback) now scans the ZUPoly's
  coefficients with `mpz_fits_slong_p` before calling
  `bz_factor_to_expr`. If any coefficient overflows int64,
  it returns false. The orchestrator's existing alpha loop
  treats this as "skip this alpha" (`if (!fac_ok) continue`),
  exactly as if the squarefree check had failed.
- `is_likely_irreducible_multivariate` checks each image
  coefficient's `Expr` type before calling `bz_factor_to_expr`;
  on `EXPR_BIGINT` it skips the alpha rather than counting
  the unchanged input as a single irreducible factor.
- `factor_bivariate_via_hensel` now tries each variable as
  main in priority order. The picker gained an `excluded_idx`
  parameter; the outer function calls it twice — once with no
  exclusion (gets the picker's preferred main variable), and
  once excluding the first attempt (falling back to the other
  variable). The polynomial-LC fallback fires only when the
  *first* picker call returns -1.

For Fateman the trace is now:

```
attempt 0: pick(excluded=-1) → main=x, kind=MONIC (lc_x = +1)
  orchestrator: alpha=0,1,-1 not squarefree
                alpha=2,-2,3,-3,4,-4,5,-5: BIGINT image, skip
                returns false (no usable evaluation point)
attempt 1: pick(excluded=0) → main=y, kind=NEGATE (lc_y = -1)
  orchestrator runs on -P at alpha_x = 0:
    image = (1-y)(1+y^k) for k = 200, squarefree, ±1 coefficients
    BZ factors cleanly → r ≥ 2 (each cyclotomic factor of 1+y^k
                                 plus the 1-y factor)
    Hensel lift in x recovers the bivariate factors of -P
  recovery: distribute the sign by negating the highest-degree
            factor → product = +P
```

This handles `(x^k + y^k + 1)(x^k - y + 1)` for `k ≤ 70` in
single-digit seconds (k=70: 1.28 s). For k = 200 (Fateman) the
Hensel lift on the high-degree y-cyclotomic factors is feasible
in principle (all coefficients stay ±1) but recombination /
trial-division on the lifted bivariate factors is slow enough
that the user's full Fateman benchmark still does not complete
in seconds. That last gap is what Tier 4 (sparse multivariate
Hensel) targets.

Display quirk that's NOT addressed: when the NEGATE path absorbs
the overall sign into one factor, the heuristic chosen is
"negate the BPoly's largest deg_x factor (tiebreak: largest
deg_y)". When `main_idx = 1` (i.e. main is the second variable
`y`), the BPoly's "x" direction is actually `y` and the
heuristic ends up negating the wrong factor for canonical-form
purposes. `fk[70]` returns
`(-1 - x^70 + y) (-1 - x^70 - y^70)`
instead of the canonical
`(1 + x^70 - y) (1 + x^70 + y^70)`.
The product is mathematically correct
(`(-A)(-B) = AB`); only the printed form is non-canonical.
A future cleanup would sense the variable swap and pick the
absorber based on the user-vars rather than the BPoly's
internal "x" direction.

#### Remaining BIGINT-incompatibility audit

Even after F6 Tier 1+2, several places in the Factor pipeline
still assume int64 coefficients. None produce *wrong* answers
(silent miscomputation) any longer; they merely *miss*
factorisations that exist. The user's correctness contract
(no wrong answers; correct-but-unsimplified is acceptable) is
preserved.

| Site | File:line | Impact |
|---|---|---|
| `extract_monomial_walk`, `extract_monomial` | `facpoly.c:403, 452` | Coefficient is `int64_t*`. Reads `e->data.integer` only when `e->type == EXPR_INTEGER` (gated). On `EXPR_BIGINT` atoms the walk returns `false`, so callers see "not a monomial" and skip. Affects `factor_binomial`. Misses, doesn't miscompute. |
| `monomial_collect`, `factor_monomial_content` | `facpoly.c:683, 749` | Same pattern: `int64_t* coeff_out`, returns `false` on BIGINT atom. Misses monomial-content extraction for inputs like `Plus[2^100 a^2 b, 2^100 b]`. |
| `factor_binomial` (`int_root` etc.) | `facpoly.c:493+` | Reads `c1, c2` as int64. Effectively unreachable on BIGINT inputs because `extract_monomial` already returned false. Misses. |
| `bz_factor_to_expr` ingestion | `facpoly.c:3793` | Reads `c->data.integer` only when `EXPR_INTEGER`; otherwise sets `non_integer_seen = 1` and returns `expr_copy(P)`. After F6 Tier 1+2 the *signal* is harmless because both call sites that previously misinterpreted it (`factor_via_bz_callback`, `is_likely_irreducible_multivariate`) now pre-filter on BIGINT. Univariate `Factor[BIGINT-coef poly]` still falls through unfactored. |
| `UPoly` / `factor_zassenhaus` substrate | `facpoly.c:1500-2055` (approx) | Entire engine is int64: `UPoly { int64_t* c }`, modular arithmetic with int64 prime, Hensel lift to int64 modulus, recombination with int64 trial division. Bigint-incompatible at every level. Tier 3 work. |
| `extract_monomial` exponents | `facpoly.c:433, 712` | Reads `e->data.function.args[1]->data.integer`. Exponents are int64; degree-2^63 polynomials don't happen, so safe in practice. |

Tier 3 (full bigint univariate BZ) replaces the entire UPoly
substrate with ZUPoly throughout `factor_zassenhaus`, including
the modular arithmetic (`mpz_mod` instead of `% p`) and Hensel
lift (lift to `p^k > 2 * Mignotte` with no int64 ceiling). Once
that lands, `extract_monomial` and friends should also be
widened to handle BIGINT atoms — same `mpz_t` substrate, same
GMP arithmetic, but the call sites are scattered enough that
it's a separate cleanup pass.

The bivariate / multivariate substrate (`ZUPoly`, `BPoly`,
`MPoly` in `src/zupoly.{c,h}`, `src/bpoly.{c,h}`,
`src/mpoly.{c,h}`) was designed bigint-clean from the start
with `mpz_t` coefficients. The recombination logic in
`mvfactor.c` and `mvfactor3.c` is also bigint-clean. The
trapped int64 assumption is now confined to the legacy
univariate engine.

### Implementation order

| Order | Phase | Reason | Status |
|---|---|---|---|
| 1 | **F1** (Wang) | Highest leverage per LOC; unblocks the most user cases (cases 3, 4, and similar non-monic bivariates). | **Done** (Stages 1+2+3 MVP) |
| 2 | **F4** (cheap squarefree check) | 50 ms - 200 ms savings on every multivariate Factor; blocks nothing. | **Done** (Stages 1+2) |
| 3 | **F3** (Hensel performance) | Required for F1 to actually be fast on high-degree inputs. | **Done** (Stages a+b+c) |
| 4 | **F5** (recombination heuristics) | Cheap, marginal but real. | **Done** (partial-lift signal short-circuit) |
| 5 | **F2** (multivariate Hensel) | Largest LOC budget; deferred if user pressure on case-2 type inputs is acceptable. | **MVP done** (n=3, monic-in-main, two-factor) |
| 6 | **F6 Tier 1+2** (BIGINT-aware skip + main-variable fallback) | Correctness: stop the false "irreducible" verdict on Fateman-class inputs AND let the orchestrator try the OTHER main variable when the picker's first choice has int64-overflowing alpha images. ~120 LOC. | **Done 2026-05-02** — `(x^k + y^k + 1)(x^k - y + 1)` now factors correctly for all k tested (k=70 in 1.28 s; mathematically correct, sign-canonicalisation cosmetic). |
| 7 | **F6 Tier 3** (full ZUPoly Berlekamp-Zassenhaus, replace int64 UPoly substrate) | Unblocks univariate inputs with large coefficients (`Factor[x^k - 1]` for k ≥ 63, `Factor[x^2 + 2^100 x + 1]`, etc.) and lets `bz_factor_to_expr` return real factorisations of BIGINT-coefficient polynomials instead of returning them unchanged. Required for the bivariate orchestrator to handle inputs where BOTH variables produce BIGINT alpha images. | Not started |
| 8 | **F6 Tier 4** (sparse multivariate Hensel without numerical evaluation) | Only path to "1 ms" on Fateman itself; biggest LOC budget. The MPoly substrate already exists; what's missing is a sparse-aware Hensel iteration that operates directly on the multivariate without going through P(x, alpha) and lift. | Not started |
| 9 | **F6 cleanup** (widen `extract_monomial` / `monomial_collect` to mpz_t; fix sign-absorber heuristic for swapped main variable) | Picks up missed factorisations of monomial-content kind on BIGINT-coefficient polynomials, and produces canonical printed form when the bivariate Hensel falls back to the second main variable. | **Stage 1 done 2026-05-07** — `monomial_collect` and `factor_monomial_content` now use `mpz_t` for the integer coefficient.  The walker correctly multiplies `EXPR_BIGINT` factors into the running coefficient instead of treating them as opaque atoms, and the residue emitter normalises via `expr_from_mpz_normalized` (INTEGER when fits in int64, BIGINT otherwise).  Two regression tests in `tests/test_factor_phase0.c::test_monomial_bigint_*` cover the path.  Stage 2 (widen `extract_monomial` / `extract_monomial_walk` + `factor_binomial`'s `c1, c2` and `int_root` to `mpz_t`) and Stage 3 (sign-absorber heuristic for swapped main variable) are still open — the binomial widening additionally needs an `mpz_root`-based variant of `int_root` (~30 LOC + tests). |
| 10 | **5c** (retire `factor_roots`) | After F1, F3.  Optional cleanup; F2 only required if we want to drop `factor_roots` BEFORE accepting case-2 limitations. | Not started |

Total: ~2-3 person-weeks for F1+F3+F4+F5+5c (without F2).  Adding
F2 doubles the budget but unblocks the last class of failures.

After F1+F3+F4+F5 landed, the remaining gap to Mathematica on the
user-reported reference cases is roughly:
- A 30× algorithmic constant factor on the simplest cases (matched
  pre-Phase 0).
- The case-2 family of inputs whose every factor depends on every
  variable -- still requires F2 multivariate Hensel.
- `factor_roots` can be retired now in principle (F1 covers all
  bivariate cases the legacy heuristic could touch); 5c is gated only
  on the test-suite pass.

## 13. References for the next implementer

- Geddes, Czapor, Labahn — *Algorithms for Computer Algebra*. The
  multivariate Hensel chapter is the practical reference for
  Phase 2-3.
- Knuth TAOCP vol 2 §4.6.1 — subresultant PRS (used for `zupoly_gcd`).
- The existing `factor_zassenhaus` code in `facpoly.c:1805` is a
  reasonably complete univariate model to mirror in the bivariate
  setting.
- Wang's leading-coefficient correction: documented in
  Geddes-Czapor-Labahn §6.6.
- For the Hilbert irreducibility heuristic, the practical
  formulation in Bernardin's "On square-free factorization of
  multivariate polynomials over a finite field" §3 is concise.
- **Barry M. Trager — *Algebraic Factoring and Rational Function
  Integration*** (Proc. 1976 ACM SYMSAC, pp. 219-226).  The core
  reference for Phase G (factoring over algebraic number fields).
  PDF in repo root: `Trager - Algebraic Factoring and Rational
  Function Integration.pdf`.

---

## 14. Phase G — Factoring over Algebraic Extensions (Trager's algorithm)

This section plans the addition of `Factor[poly, Extension -> α]` where
`α` is an algebraic number satisfying a minimal polynomial
`P_α(y) ∈ Q[y]`.  The reference algorithm is Barry Trager's *Algebraic
Factoring and Rational Function Integration* (1976) — the standard
method used by Mathematica, Macsyma, Maple, and SageMath.

Phase G is **complementary** to the existing factoring pipeline rather
than a replacement.  The inner loop calls the existing `Factor[]` over
`Q[x]` (the `bz_factor_to_expr` substrate from §3.4 / §F6) on the norm
`R(x) = Resultant_y(P_α(y), f(x, y))`, then lifts the rational
factorisation back to `Q(α)[x]` via `gcd` over `Q(α)[x]`.

### 14.1 Why now

Three things make Trager fit cleanly into the current pipeline:

1. **Resultant infrastructure exists.**  `poly_resultant` and the
   `Resultant`/`Discriminant` builtins (`src/poly.c:2616`+) already
   compute resultants in a chosen variable.  The norm
   `Norm(f) = Resultant_y(P_α(y), f(x, y))` is a direct application.
2. **Squarefree univariate `Factor` over Q is solid.**  Phases 0-5 +
   F6 Tier 1+2 deliver a substrate that handles bigint coefficients,
   multi-α evaluation strategies, and recombination.  The norm is fed
   straight in.
3. **Algebraic-generator groundwork.**  `simp_algebraic` (`src/simp.c`,
   commit `57f505a`) and `lift_algebraic_generator_pass` (`src/poly.c`,
   commit `548bee3`) already detect algebraic generators in expressions
   and lift them into a variable-based representation.  Phase G builds
   on the same machinery for input pre-processing and output
   canonicalisation.

### 14.2 Target test cases

| Input | Expected output over Q(α) |
|---|---|
| `Factor[x^4 + 1, Extension -> Sqrt[2]]` | `(x^2 - Sqrt[2] x + 1)(x^2 + Sqrt[2] x + 1)` |
| `Factor[x^2 - 2, Extension -> Sqrt[2]]` | `(x - Sqrt[2])(x + Sqrt[2])` |
| `Factor[x^4 - 5 x^2 + 6, Extension -> Sqrt[2]]` | `(x - Sqrt[2])(x + Sqrt[2])(x^2 - 3)` |
| `Factor[x^4 - 5 x^2 + 6, Extension -> {Sqrt[2], Sqrt[3]}]` | `(x - Sqrt[2])(x + Sqrt[2])(x - Sqrt[3])(x + Sqrt[3])` |
| `Factor[x^2 + x + 1, Extension -> Sqrt[-3]]` | `(x - (-1 + Sqrt[-3])/2)(x - (-1 - Sqrt[-3])/2)` |
| `Factor[x^3 - 2, Extension -> 2^(1/3)]` | `(x - 2^(1/3))(x^2 + 2^(1/3) x + 2^(2/3))` |

### 14.3 Algorithm reference

Trager's full pipeline (paper §2):

```
Algorithm sqfr_norm:
  input : f(x, α) ∈ Q(α)[x] squarefree, P_α the minimal polynomial of α
  output: integer s ≥ 0,
          g(x, α) = f(x - s α, α),
          R(x)    = Norm(g(x, α))   (squarefree over Q)

  s := 0; g := f
  loop:
    R(x) := Resultant_y( P_α(y), g(x, y) )      (substitute α → y)
    if gcd(R, R') is a unit: return (s, g, R)
    s := s + 1
    g(x, α) := f(x - s α, α)
```

```
Algorithm alg_factor:
  input : f(x, α) ∈ Q(α)[x] squarefree
  output: irreducible factors of f over Q(α)

  (s, g, R) := sqfr_norm(f, P_α)
  L := Factor(R(x))                              (existing facpoly: over Q)
  if |L| = 1: return [f]                         (f irreducible over Q(α))
  factors := []
  for h(x) in L:
    h(x, α)   := gcd_{Q(α)[x]}( h(x), g(x, α) ) (monic gcd over Q(α))
    g         := g / h(x, α)                    (exact division in Q(α)[x])
    h(x, α)   := h(x + s α, α)                  (undo linear shift)
    factors.append( h(x, α) )
  return factors
```

Soundness:
- **Theorem 2.1** (Trager): `Norm(f)` is a power of `Norm(f₁)` whenever
  `f` is irreducible over `Q(α)` and `f₁` is its image-shift conjugate.
- **Theorem 2.2** (Trager): for squarefree-norm input, the
  factor-by-factor gcd lift recovers the true factors of `f` over `Q(α)`.
- **Theorem 2.3 / Corollary 2.5** (Trager): only finitely many `s` give
  a non-squarefree norm — `sqfr_norm` always terminates.

### 14.4 Phased implementation plan

#### Phase G1 — `Q(α)` element type (DONE 2026-05-07)

`src/qa.{c,h}`: an element of `Q(α) = Q[α]/(P_α(α))`.

- `QANum`: a polynomial in `α` of degree `< deg P_α`, `mpq_t`
  coefficients, dense representation.
- `QAExt`: the extension's minimal polynomial `P_α`, stored as an
  `mpq_t[]` of length `deg P_α + 1`, plus `deg P_α`.
- Operations: `qa_zero`, `qa_one`, `qa_alpha`, `qa_from_mpq`,
  `qa_from_si`, `qa_copy`, `qa_add`, `qa_sub`, `qa_neg`, `qa_mul`
  (multiply then reduce mod `P_α` via long division), `qa_inverse`
  (extended Euclidean over `Q[y]` → Bézout coefficients), `qa_div`,
  `qa_eq`, `qa_is_zero`, `qa_is_one`, `qa_scale_mpq`, `qa_scale_si`.
- Convenience extension constructors: `qaext_sqrt_si(c)`
  (P_α = y² − c), `qaext_root_si(c, n)` (P_α = y^n − c).
- Tests: 18 cases in `tests/test_qa.c` covering arithmetic identities,
  `√2 · √2 = 2` reduction, `α³ ≡ 2` in Q(2^(1/3)),
  `(3+4i)(3−4i) = 25` in Q(i), `1/(1+√2) = √2 − 1` extended-Euclidean
  inverse, `1/α = α²/2` in Q(2^(1/3)), and division round-trip.
- Implementation note: the polynomial arithmetic for multiplication-
  with-reduction and extended-Euclidean inversion is inlined as
  plain `mpq_t` array manipulation, deliberately avoiding the static
  `QUPoly` substrate in `zupoly.c`.  Promoting `QUPoly` to a public
  module is a separable refactor (orthogonal to G1 and not blocking
  G2-G7).
- Actual: 232 LOC source (`src/qa.h` 95 + `src/qa.c` 137 minus
  comments) + 332 LOC tests.  All 95 test binaries pass.

#### Phase G2 — `Q(α)[x]` univariate polynomial type (DONE 2026-05-07)

`src/qaupoly.{c,h}`: univariate polynomial in `x` with `QANum`
coefficients.

- Storage: dense `QANum*` array indexed by x-degree (mirrors `ZUPoly`).
- Operations: `qaupoly_add`, `_sub`, `_neg`, `_mul`, `_scale_qa`,
  `_scale_si`, `_divrem` (long division over Q(α), which is a field
  for irreducible P_α so the leading-coefficient inverse always
  succeeds), `_gcd` (Euclidean PRS, returned monic),
  `_make_monic`, `_shift` (`f(x + c)` for `c ∈ Q(α)` via Horner-style
  rebuild), `_eval` (evaluate at `c ∈ Q(α)` via Horner).
- Construction: `qaupoly_new`, `_zero`, `_from_qa`, `_from_si`,
  `_x` (the indeterminate), `_copy`, `_setcoef`, `_getcoef`,
  `_normalize`.
- Predicates: `qaupoly_is_zero`, `qaupoly_eq`.
- Tests: 17 cases in `tests/test_qaupoly.c`.  Headline test:
  `gcd(x² − 2, x − √2) = x − √2` over Q(√2) (the Trager lift step).
  Other coverage: `divrem` exact + with-remainder; `gcd` coprime / self
  / over Q(i); arithmetic identities including `(x − √2)(x + √2) = x² − 2`;
  evaluation at √2 of `x² − 2` (yields 0); shift-then-shift-back round
  trip recovering the original poly.
- Implementation note: Phase G2 hardens `test_qaupoly.c` with a non-
  elided `ASSERT` macro because Release builds (default cmake) define
  `-DNDEBUG`, which compiles `assert(cond)` to `((void)0)` and elides
  side-effecting calls.  This caught a real bug in the test suite where
  `ASSERT(qaupoly_divrem(...))` left `q`/`r` uninitialised in Release;
  the fix landed both as the local override and as a project-wide
  lesson in `tasks/lessons.md`.
- Actual: 200 LOC source (`src/qaupoly.h`) + 320 LOC implementation
  (`src/qaupoly.c`) + 360 LOC tests.  All 96 test binaries pass.

#### Phase G3 — Norm via resultant — **Done 2026-05-07**

`src/qafactor.{c,h}::qaupoly_norm`:

- Build `g(x, y) ∈ Q[x, y]` from `f(x, α) ∈ Q(α)[x]` by substituting
  `α → y`.  Each `QANum` coefficient becomes a polynomial in `y`.
- Compute `R(x) = Resultant_y(P_α(y), g(x, y))` via `internal_resultant`
  (which dispatches to the existing Sylvester-matrix routine in
  `src/poly.c::resultant_internal`).  The result lands in `Q[x]`
  directly — denominator clearing is unnecessary because picocas's
  Resultant builtin handles `Rational[]` coefficients natively, and
  the Sylvester determinant of a `Q[x, y]` system is already in `Q[x]`.
- The result is post-Expanded for canonical sum-of-monomials form.
- Two helpers (`qaext_to_expr`, `qaupoly_to_expr`) lift `QAExt` /
  `QAUPoly` into picocas Expr trees — they're exposed publicly because
  G4's `sqfr_norm` will reuse them when computing norms of shifted
  inputs `f(x − sα, α)`.
- Tests: `tests/test_qafactor.c` (11 cases) — covers `Norm(x − √2) = x²
  − 2`, `Norm(x − ∛2) = x³ − 2`, rational-input squared in `Q(i)`,
  `Norm(x² − √2) = x⁴ − 2`, `Norm(x² + √2 x + 1) = x⁴ + 1` (cyclotomic
  Φ_8), and the constant case `Norm(3) = 9` over `Q(√2)`.
- Actual: 95 LOC header (`src/qafactor.h`) + 158 LOC implementation
  (`src/qafactor.c`) + 200 LOC tests.

#### Phase G4 — `sqfr_norm` and `alg_factor` — **Done 2026-05-07**

`src/qafactor.c` (extending the G3 module):

- `qa_sqfr_norm(QAUPoly* f, x_name, y_name) → QASqfrNormResult`,
  where `QASqfrNormResult = { int s; QAUPoly* g; Expr* R; }`.  The
  squarefree-norm algorithm walks `s = 0, 1, 2, …`, builds
  `g(x) = f(x − sα) ∈ Q(α)[x]` via `qaupoly_shift`, recomputes
  `R(x) = qaupoly_norm(g) ∈ Q[x]`, and returns the first iteration for
  which `gcd(R, R')` over `Q[x]` is constant (i.e. `R` is squarefree).
  Squarefree-detection routes through picocas's existing `D` and
  `PolynomialGCD` builtins via `internal_call_impl` — no new gcd
  machinery is added.  Trager's theorem 2.3 / corollary 2.5 guarantees
  termination; we cap at `QA_SQFR_NORM_MAX_TRIES = 32` shifts as a
  defensive bound (in practice `s ∈ {0, 1, 2}`).
- `qa_alg_factor(QAUPoly* f, x_name, y_name, int* n_out) → QAUPoly**`.
  Calls `sqfr_norm`, factors `R(x) ∈ Q[x]` via `internal_factor`
  (single-arg `Factor[poly]` — it routes Q[x] inputs through the
  existing `Together` / `bz_factor_to_expr` Z[x] pipeline so we don't
  need a `ZUPoly` round-trip), walks the resulting `Times[…]` /
  `Power[…]` tree to collect Q-rational factors `h_i(x)`, lifts each
  to a Q(α)-factor of `g` via `qaupoly_gcd(h_i, g)` (Phase G2's
  Euclidean monic gcd), divides the lift out of `g`, and undoes the
  linear shift `x → x + sα` via `qaupoly_shift`.  When `R` is
  Q-irreducible (single factor walked out), returns `[make_monic(f)]`.
- Two implementation gotchas worth recording:
  1. `PolynomialGCD` is variadic over polynomials — passing a third
     `x` symbol silently coerces it into the three-way gcd
     `gcd(R, R', x)` which is always 1 (a phantom-squarefree result).
     Fix: pass exactly two args.
  2. `Factor` requires `arg_count == 1`; passing a `(poly, var)` pair
     causes `builtin_factor` to return `NULL`, leaving the caller with
     an unevaluated `Factor[…]` tree that the factor-walker then
     misinterprets as a single irreducible.  Fix: `Factor[poly]`
     accepts a univariate Q[x] input directly.
- Tests: 9 new cases in `tests/test_qafactor.c` (now 20 in total)
  covering the six `§14.2` rows — `x²−2`, `x⁴+1`, `x⁴−5x²+6`, `x³−2`,
  `x²+1`, `x²+x+1` over `Q(√2) / Q(√2) / Q(√2) / Q(∛2) / Q(i) / Q(√−3)`
  respectively — plus the irreducible-input edge case (`x²−3` over
  `Q(√2)`) and two direct sqfr_norm tests (`s = 0` for already-square­
  free, `s ≥ 1` for Q-rational input where Norm is `f²`).
- Actual: 95 LOC header (delta from G3) + 220 LOC source delta (G3
  module grew from 158 to 378 LOC) + 270 LOC test delta (test file
  grew from 200 to 470 LOC).  All 97 test binaries pass.

#### Phase G5 — picocas-level API — **Done 2026-05-07**

User-facing entry point: `Factor[poly, Extension -> α]` now routes
through Trager for any `α` whose minimal polynomial we can build.

- Public API in `qafactor.{c,h}`:
  - `qa_resolve_extension(alpha_expr, **render_out) -> QAExt*`
    recognises:
    - `Sqrt[c]` (integer c, including negative) → `y² − c`
    - `c^(1/n)` (integer c, integer n ≥ 2) → `yⁿ − c`
    - `I` / `Complex[0, 1]` → `y² + 1`
    - `Times[I, Sqrt[c]]` (the canonical form of `Sqrt[-c]` after
      picocas's auto-evaluation) → `y² + c`
    Picocas canonicalises `Sqrt[c]` as `Power[c, Rational[1, 2]]`,
    so the recogniser accepts both spellings.
  - `qa_expr_to_qaupoly_with_alpha` (static): walks an Expr
    polynomial in `var` whose coefficients may contain `alpha_render`
    and lifts it to a QAUPoly.  Strategy: substitute the surface form
    of α with an internal placeholder symbol `$qa$alpha$`, then peel
    off `(varⁱ)(αʲ)` terms via `Coefficient[…]`.  Returns `NULL` if
    any coefficient is not in `Q(α)` (caller falls back to plain
    Factor).
  - `qaupoly_to_expr_alpha`: round-trips a QAUPoly through
    `qaupoly_to_expr` with the internal placeholder, then substitutes
    the placeholder back to the user's surface form and re-evaluates
    so that `Sqrt[c]^k` collapses to its canonical form.
  - `qa_factor_with_extension(poly, alpha_expr, var)`: top-level
    orchestrator.  Includes a squarefree-reduction pre-pass via
    `gcd(f, f') / gcd` and a multiplicity-detection trial-division
    loop, so non-squarefree inputs like `(x²−2)²` come out as
    `(x−√2)²(x+√2)²` rather than just `(x−√2)(x+√2)`.
- Hook in `builtin_factor` (`src/facpoly.c`): the arg-count guard
  now allows trailing options.  Trailing `Rule[Extension, α]` is
  picked up; multivariate inputs (after stripping α) are rejected
  back to plain Factor.
- Internal symbol: added `SYM_Extension` to `sym_names.{c,h}`.
- Surface-form gotchas (documented as comments in `qa_resolve_extension`):
  - The literal `I` auto-evaluates to `Complex[0, 1]`, so the
    resolver tests both spellings.
  - `Sqrt[-c]` for positive c auto-evaluates to `Times[I, Sqrt[c]]`,
    so the resolver matches that compound form too.
  - `Sqrt[c]` is stored canonically as `Power[c, Rational[1, 2]]`.
- Tests: 11 new tests in `tests/test_qafactor.c` covering the §14.2
  headlines plus `x⁴+4` (4 linear factors over Q(i)), an α-bearing
  input (`x² − 2√2 x + 2 = (x − √2)²` over Q(√2)), a repeated-factor
  input ((x²−2)² over Q(√2)), and an irreducible-over-extension
  guard (x² − 3 over Q(√2) stays unfactored).  All 31 qafactor tests
  pass; full 97-binary suite green.
- Updates: `Factor` docstring in `src/info.c` mentions the new
  option.  No additions to internal `.m` files (the dispatch is
  entirely in C).
- Out of MVP: `Root[p, k]` (the symbolic-root form) and
  `AlgebraicNumber[α, coeffs]` (Mathematica's normalised
  representation).  Both can land later without disturbing the
  current API.

#### Phase G6 — Tower of extensions

For inputs like `Factor[..., Extension -> {Sqrt[2], Sqrt[3]}]` Trager
§3 gives a primitive-element algorithm (`primitive_element`) that
reduces a tower `Q(α, β)` to a single-extension form
`Q(γ) = Q(α + s β)` for an integer `s ≥ 0` chosen so the resultant is
squarefree.  Iterates: compute `s` via `sqfr_norm`; compute
`R(x) = Norm(Q_β(x, α))`; factor `R`; locate `α` and `β` in `Q(γ)` via
`linsolve` of the gcd.

- Estimated: ~250 LOC source + ~150 LOC tests.

#### Phase G7 — Splitting fields

Trager §4 (Algorithm `split_field`): repeatedly apply
`primitive_element` to construct the splitting field of an irreducible
polynomial over `Q`.  Useful for partial-fraction decomposition over
the minimal extension (Phase G applications).  Out of MVP scope; can
land later as the foundation for symbolic integration over algebraic
extensions (Trager §5).

- Estimated: ~300 LOC source + ~200 LOC tests.

### 14.5 Connections to existing code

| Existing | Used by Phase G | Where |
|---|---|---|
| `poly_resultant` | G3 norm computation | `src/poly.c` |
| `bz_factor_to_expr` | G4 inner factor over `Q` | `src/facpoly.c:3793` |
| `lift_algebraic_generator_pass` | G5 input pre-processing | `src/poly.c` (commit `548bee3`) |
| `simp_algebraic` | G5 output canonicalisation | `src/simp.c:4043` |
| `Resultant`, `Discriminant` | exposed builtins | already in symtab |
| `ZUPoly` | G3 intermediate, G4 norm representation | `src/zupoly.{c,h}` |
| `BPoly` | G3 bivariate-resultant intermediate | `src/bpoly.{c,h}` |
| `MinimalPolynomial[Sqrt[c]]` | G5 extension resolution | `simp_algebraic` already computes this |

### 14.6 Risks and open questions

1. **Norm performance.**  Bivariate resultants are
   `O(deg(P_α)³ · deg_x(f)²)` via subresultant PRS; for typical inputs
   `deg(P_α) ≤ 4`, so the resultant computation is bounded.  High-degree
   `α` (cyclotomic over Q, degree 100+) would need a modular resultant.
2. **Coefficient blow-up in `f(x - sα)`.**  Each shift increases
   coefficient sizes by `O(s · deg P_α)`.  In practice `s ≤ 3` for
   nearly all inputs; the paper conjectures (and we observe) that
   `s ∈ {0, 1}` suffices for most real-world cases.
3. **Squarefree input pre-processing.**  `alg_factor` requires
   squarefree input.  The reduction is `gcd(f, f')` over `Q(α)[x]`,
   which Phase G2 supports directly.
4. **Real vs. complex extensions.**  Trager's algorithm is field-
   agnostic — `Q(i) = Q[y]/(y²+1)` is handled identically to `Q(√2)`.
   Picocas's `Sqrt[-c]` lifts to this directly.
5. **Multivariate input `f(x_1, ..., x_n, α)` over `Q(α)`.**  Falls
   out of Phase G4 *if* the existing bivariate Hensel is taught to
   accept `Q(α)`-valued `cy[]` coefficients.  Significant additional
   scope, deferred to a Phase G8 if user demand materialises.
6. **Interaction with `simp_algebraic`.**  The Simplify-level algebraic-
   reduction logic should not double-process inputs: when `Factor` is
   called with an explicit `Extension`, the algebraic generator is
   resolved at parse time and the input is normalised over `Q(α)`
   *before* entering `alg_factor`.

### 14.7 Estimated total budget

| Phase | LOC source | LOC tests |
|---|---|---|
| G1 — `Q(α)` element type | 250 | 200 |
| G2 — `Q(α)[x]` univariate | 300 | 250 |
| G3 — Norm via resultant | 150 | 100 |
| G4 — sqfr_norm + alg_factor | 250 | 200 |
| G5 — picocas API | 200 | 100 |
| G6 — tower of extensions | 250 | 150 |
| G7 — splitting fields (deferred) | 300 | 200 |
| **MVP total (G1-G5)** | **1150** | **850** |
| **Full Phase G (G1-G6)** | **1400** | **1000** |

Approximately **2-3 person-weeks** for the MVP `Factor[..., Extension
-> Sqrt[c]]`; **3-4 person-weeks** for the full Phase G including
towers; **+1 week** for G7.

### 14.8 Where Phase G slots into the overall implementation order

Phase G has **no hard dependencies** on outstanding F-work: it can
land in parallel with F6 Tier 3/4 if there is appetite.  F6 Tier 3
(full ZUPoly Berlekamp-Zassenhaus) does benefit Phase G indirectly
because the norm `R(x)` may have large coefficients when `deg P_α` is
large.

Updated implementation order (extending §12's table):

| Order | Phase | Reason | Status |
|---|---|---|---|
| 1-6 | (existing F1..F6 Tier 1+2) | (see §12 implementation order) | Done |
| 7 | F6 Tier 3 | Bigint univariate BZ | Not started |
| 8 | F6 Tier 4 | Sparse multivariate Hensel | Not started |
| 9 | F6 cleanup | Widen `extract_monomial` to mpz_t; sign-absorber heuristic | Not started |
| 10 | 5c | Retire `factor_roots` | Not started |
| 11 | **G1** (Q(α) element type) | Self-contained; parallelisable with F6 cleanup. | **Done 2026-05-07** — `src/qa.{c,h}` + 18 tests in `tests/test_qa.c`, all passing. |
| 11b | **G2** (Q(α)[x] univariate) | Polynomial-in-x with QA coefficients. | **Done 2026-05-07** — `src/qaupoly.{c,h}` + 17 tests in `tests/test_qaupoly.c`, all passing.  Headline: `gcd(x²-2, x-√2) = x-√2` over Q(√2) (the Trager lift step). |
| 11c | **G3** (norm via resultant) | Field-norm of `f ∈ Q(α)[x]` via `Resultant_y(P_α, f)`. | **Done 2026-05-07** — `src/qafactor.{c,h}::qaupoly_norm` + 11 tests in `tests/test_qafactor.c`, all passing.  Headlines: `Norm(x − √2) = x² − 2`, `Norm(x − ∛2) = x³ − 2`, `Norm(x² + √2 x + 1) = x⁴ + 1` (cyclotomic Φ_8). |
| 12 | **G4** (sqfr_norm + alg_factor MVP) | First user-visible factoring over `Q(√c)`. | **Done 2026-05-07** — `src/qafactor.{c,h}::qa_sqfr_norm` + `qa_alg_factor` + 9 new tests (20 total in `test_qafactor.c`).  Headlines: `x²−2` over `Q(√2)` → `(x−√2)(x+√2)`; `x⁴+1` over `Q(√2)` → `(x²−√2 x+1)(x²+√2 x+1)`; `x⁴−5x²+6` over `Q(√2)` → `(x−√2)(x+√2)(x²−3)`; `x³−2` over `Q(∛2)` → `(x−∛2)(x²+∛2 x+∛4)`; `x²+1` over `Q(i)`, `x²+x+1` over `Q(√−3)`, plus the `x²−3` over `Q(√2)` irreducible-input edge case. |
| 13 | **G5** (picocas API) | Wire `Extension -> ...` into `Factor`. | **Done 2026-05-07** — hook in `src/facpoly.c::builtin_factor` reads the `Rule[Extension, α]` option; `qa_factor_with_extension` in `src/qafactor.c` dispatches to Trager and renders results.  11 new tests (31 total in `test_qafactor.c`).  Headlines: `Factor[x^2-2, Extension->Sqrt[2]]` → `(x-Sqrt[2])(x+Sqrt[2])`; `Factor[x^4+1, Extension->Sqrt[2]]` → `(x^2-Sqrt[2] x+1)(x^2+Sqrt[2] x+1)`; `Factor[x^3-2, Extension->2^(1/3)]` → `(x-2^(1/3))(x^2+2^(1/3) x+2^(2/3))`; `Factor[x^2+1, Extension->I]` → `(x-I)(x+I)`; `Factor[x^2+x+1, Extension->Sqrt[-3]]` → `(x+1/2-I Sqrt[3]/2)(x+1/2+I Sqrt[3]/2)`; α-bearing inputs (`x^2 - 2 Sqrt[2] x + 2` → `(x-Sqrt[2])^2`); repeated factors with multiplicity ((x²-2)² → (x−√2)²(x+√2)²). |
| 14 | **G6** (tower of extensions) | Multi-radical inputs. | Not started |
| 15 | **G7** (splitting fields) | Foundation for symbolic integration over algebraic extensions (Trager §5). | Deferred |

Recommended near-term sequence:
1. **F6 cleanup (item 9)** — small, immediate value; widens
   `extract_monomial` to handle bigint-coefficient monomials.
2. **G1 + G2** (Trager substrate) — isolated, well-tested; ~5 days.
3. **G3 + G4** (norm + factor) — the first user-visible Phase G
   capability; depends only on G1+G2 + the existing `poly_resultant`
   and `bz_factor_to_expr` substrate.
4. **G5** (`Extension -> ...` syntax) — small wire-up step.

That pre-G6 sequence delivers `Factor[x^4 + 1, Extension -> Sqrt[2]]`
and friends in approximately two weeks of focused work, ahead of the
larger F6 Tier 3/4 efforts.
