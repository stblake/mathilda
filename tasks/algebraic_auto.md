---
title: Automatic algebraic-extension detection across Together, Cancel, and polynomial functions
date_started: 2026-05-12
status: phases A + B + C + E + F + multi-gen shipped 2026-05-13; phase D pending
---

## Status (2026-05-13)

Phases A, B (full polynomial-function coverage), C (Simplify-side
cube-root generalisation), E (nested-radical walker), F
(canonical-form dedup), and multi-generator tower routing shipped.

**Done:**
- Phase A — `extension_autodetect` + `extension_autodetect_args` in
  `src/qafactor.{c,h}`, plus `extract_extension_option_full` in
  `src/options.{c,h}`.  18 unit tests pass
  (`tests/test_extension_autodetect.c`).
- Phase B — auto-detect wired into `Cancel`, `Together`,
  `PolynomialGCD`, `PolynomialLCM`, `PolynomialQuotient`,
  `PolynomialRemainder`, `PolynomialQuotientRemainder`, `Apart`,
  `Factor`.  Tier-1 routes single-generator detections through the
  existing `Extension -> α` path; multi-generator towers (n ≥ 2)
  detect but fall back to no-extension.  17 integration tests pass
  (`tests/test_extension_auto_builtins.c`).
- Phase C — `simp_algebraic` (`src/simp.c`) cube-root and higher
  shortcut.  When the input has exactly one rational-base radical
  generator (any q ≥ 2), simp_algebraic routes through
  `Together[expr, Extension -> α]` (qaupoly substrate), instead of
  bailing on q != 2 as it did before.  12 unit tests pass
  (`tests/test_simp_algebraic_cuberoot.c`, includes 3 nested-radical
  cases added under Phase E).
- Phase E — nested-radical walker.  `extension_autodetect` now
  surfaces `Sqrt[non-integer-base]` and `Power[non-integer-base, p/q]`
  as a separate `GEN_NESTED` generator class, dispatched at tower-
  build time to `qa_resolve_nested_radical` (Phase G8).  Previously
  these caused the entire walker to bail.  20 unit tests pass
  (`tests/test_extension_autodetect.c`, includes 2 new nested cases).
- Multi-generator tower routing — `qa_cancel_with_tower`
  (`src/qafactor.{c,h}`) lifts n ≥ 2 detections into the qaupoly
  substrate via γ-substitution.  Wired into `builtin_cancel` and
  `builtin_together`; previously these fell through on n ≥ 2.  2 new
  integration tests (`tests/test_extension_auto_builtins.c`).
- Phase F — canonical-form dedup of nested radicals.
  `extension_autodetect` runs a post-walker pass that canonicalises
  each nested radicand via `Together[..., Extension -> α]` and then
  rewrites it back into atomic-radical form
  (`autodetect_rewrite_to_atomic`).  Structurally-distinct-but-
  mathematically-equal `Sqrt[u]` terms now collapse to a single
  generator.  An auxiliary `autodetect_walk_intbase_only` walker
  surfaces integer-base generators inside nested radicands while
  skipping inner atomic radicals (which G8 already handles), avoiding
  the double-count regression on `Sqrt[1+Sqrt[2]]`.  1 new test
  (`tests/test_simp_algebraic_cuberoot.c`).
- Regression status: 94 pass / 16 fail (was 92 / 16 — identical
  failure set; new passing test binaries for the auto-builtins and
  cube-root tests).

**Not started:**
- Phase D — `Resultant`, `Discriminant`, `FactorSquareFree`,
  `PolynomialQ` Extension threading.  (`PolynomialMod` is out of
  scope: it's a coefficient-modulus operation, not a polynomial-
  modulus operation, so `Extension` has no semantic for it.)
- Phase F — documentation in `docs/spec/builtins/*.md` and per-builtin
  docstring updates.
- Extending G8 + `qa_resolve_extension` + `qa_expr_to_qaupoly_with_alpha`
  to accept `Power[c, p/q]` for general integer `p` (not just `p == 1`).
  The blocker is Mathilda's Times canonicaliser at `src/times.c:481`,
  which deliberately absorbs coefficient factors into Power exponents
  (`(1/6)·2^(1/3) → (1/3)·2^(-2/3)`).  Phase F's canonicalise-dedup
  produces the equalised radicand structurally, but the canonical
  form contains `Power[c, p/q]` with `p != 1`, which G8 rejects.
  Fixing this requires either changing the canonicaliser (wide impact)
  or extending the qa substrate to accept arbitrary `p` and lift via
  `qa_inverse` for negative exponents (multi-hundred-LoC qa change
  across substitution and render-back paths).  Out of scope for the
  current session.
- Render-back cosmetic minimisation — `qa_cancel_with_tower` outputs
  expressions in terms of γ (the primitive element), not the user's
  original α_i basis.  Mathematical equivalence is preserved but the
  surface form can be unwieldy.

**Motivating example status.**  `Together[D[Integrate[...], x],
Extension -> Automatic]` now produces an intermediate form that
matches Mathematica's `Together[..., Extension -> Automatic]` output
exactly.  Multi-generator routing now ships but doesn't close the
gap — the tower built from `2^(1/3)` plus two distinct (but
mathematically equal) `Sqrt[...]` generators either fails to
construct, or returns a result expressed in a γ-basis that doesn't
visibly collapse to `a x^2 / (x^3 + 2)`.

Closing the gap likely requires *canonical-form dedup before tower
construction*: pre-normalise each nested `Sqrt[u]` radicand against
the already-collected integer-base generators (e.g.
`Sqrt[-1/(9·2^(2/3)) + (2/9)·2^(1/3)]` and `Sqrt[1/(3·2^(2/3))]`
should both reduce to `Sqrt[2^(1/3)/6]`).  Once both Sqrts collapse
to the same generator, the tower has 2 distinct generators (2^(1/3)
and Sqrt[2^(1/3)/6]) instead of 3, and the resulting compositum is
small enough for the surface form to simplify.

Mathilda's `Simplify` on this input separately hits a pre-existing
`1/0` Power-error in a downstream transformation (not related to
this work).

# Automatic algebraic-extension detection (Tier 1)

## Motivating example

```
In[1]:= Integrate[a x / (x^3 + 2), x]
        (* user antiderivative — contains 2^(1/3), 2^(2/3),
           Sqrt[-1/(9·2^(2/3)) + 2/9·2^(1/3)], Sqrt[1/(3·2^(2/3))] *)

In[2]:= D[%, x]                                (* ≈ 3-term sum, mixed radicals *)

In[3]:= Simplify[%]    (* Mathematica: ~1 ms.  Mathilda: hangs / no-op. *)
```

The reduced antiderivative lives in `Q(2^(1/6))(x)` (degree 6).  Once both
`Sqrt`s are recognised as elements of `Q(2^(1/6))` and `2^(1/3) · 2^(2/3)`
reduces to 2 via the minimal polynomial `θ^6 − 2 = 0`, the entire
expression simplifies to `a x² / (x³ + 2)` — the derivative of the input.

Mathilda already has the algebraic substrate (`qa.c`, `qaupoly.c`,
`qafactor.c`, `QATower`).  The single missing piece is **auto-detection**:
when `Together`/`Cancel`/`Factor`/`Apart`/`PolynomialGCD`/`PolynomialLCM`/
`PolynomialQuotient`/`PolynomialRemainder` are called without an explicit
`Extension → α` option (or with `Extension → Automatic`), scan the input
for algebraic generators and feed them to the existing
`qa_resolve_extension_tower` entry point.

## What already exists (no work needed)

| Layer | Component | Status |
|---|---|---|
| Q(α) arithmetic mod minpoly | `src/qa.c` (`QAExt`, `QANum`, `qa_inverse`) | works, any minpoly |
| Q(α)[x] (univariate over the extension) | `src/qaupoly.c` (`qaupoly_gcd`, `qaupoly_divrem`) | works |
| α-Resolver: `Sqrt[c]`, `c^(1/n)`, `I` → `QAExt` | `qafactor.c::qa_resolve_extension` | works |
| Multi-generator tower (compositum, Trager §3) | `qafactor.c::qa_resolve_extension_tower` (`QATower`) | works |
| `Extension → α` plumbing on PolynomialGCD/LCM/Quotient/Remainder/Factor | `poly.c`, `facpoly.c` | works |
| `Extension → α` plumbing on Cancel/Together | `rat.c::cancel_with_extension`, `together_recursive_ext` | works |
| `Extension → α` plumbing on Apart | `parfrac.c` | works |
| Sqrt-only multi-generator rationalisation in Simplify | `simp.c::simp_algebraic` (q == 2 only) | works for Sqrt; bails on cube roots |

## What's missing

1. **Auto-detect generators** when no `Extension` (or `Extension → Automatic`)
   is passed.  Today `extract_extension_option` returns NULL → the
   non-extension fallback path runs, which on the motivating example is
   either a no-op or a multivariate Q[α, x] GCD blowup.

2. **q > 2 in `simp_algebraic`** (Simplify pipeline).  Line 3610 in
   `src/simp.c` bails on any `Power[u, p/q]` with q != 2.  The qa.c
   substrate handles q > 2 fine; only the Simplify-level wrapper hasn't
   been generalised.

3. **Connect the dots**: a single `extension_autodetect(expr) → QATower*`
   helper that lives next to `qa_resolve_extension_tower` and is invoked
   by every public extension-aware entry point when `alpha == NULL`.

## Concrete plan

### Phase A — `extension_autodetect`: collect generators from an expression

**New helper** in `qafactor.c` (declared in `qafactor.h`):

```c
/* Walk `e` collecting every algebraic generator: a sub-expression of
 * the form Power[r, p/q] with rational r and q > 1, or Sqrt[u] where
 * u is rational.  Returns the unique set (by structural equality on
 * the canonical surface form, e.g. 2^(2/3) is the same generator-class
 * as 2^(1/3) — both contribute the single generator 2^(1/3)).
 *
 * Returns NULL on:
 *   - no algebraic generators in the input (caller stays in Q),
 *   - generator outside the qa_resolve_extension recognised set
 *     (e.g. Sqrt[polynomial-in-x]),
 *   - more than QA_AUTODETECT_MAX_GENS (default: 4) distinct generators.
 *
 * Caller owns: qa_tower_free. */
QATower* extension_autodetect(const struct Expr* e);
```

**Algorithm:**

1. Walk `e` collecting every `Power[base, p/q]` with `q != 1` and every
   `Sqrt[u]` (which is `Power[u, 1/2]`).  Skip occurrences inside
   `Hold`/`HoldForm` and inside the function head.
2. For each occurrence, derive the canonical generator surface form:
   - `Power[r, p/q]` with `r ∈ Q`, `gcd(p, q) = 1`, `q > 1`:
     canonical generator is `Power[r, 1/q]` (or `r^(1/q')` where q' = q/gcd(p,q)).
     When the same base `r` appears with multiple denominators
     `q_1, ..., q_k`, use `q = lcm(q_i)` and the generator is
     `Power[r, 1/q]`.  All occurrences are then integer powers of this
     single generator.
   - `Sqrt[u]` with `u ∈ Q`: generator is `Sqrt[u]`.
   - `Sqrt[u]` with `u` non-rational but matching a generator already in
     the set after reduction (e.g. `Sqrt[2^(1/3)/6]` is `2^(1/6)/Sqrt[6]`):
     treat as inducing a *new* generator `Sqrt[u]` only after attempting
     a one-pass reduction of `u` against the generators collected so far.
     Tier-1 simplification: if `u` reduces to `(rational) · θ^k` for
     existing generator `θ`, fold into the existing generator's degree
     (LCM-bump) instead of adding a new generator.
3. Deduplicate: same `r` → bump the radical-degree LCM; literally
   identical surface forms collapse to one.
4. Build the generator list `[α_1, ..., α_n]` and delegate to
   `qa_resolve_extension_tower`.

**Bail conditions** (return NULL, caller uses Q-only path):

- Any `Power[u, p/q]` with `u` non-rational and not a previously-seen
  generator that reduces to a polynomial in θ.
- Any `Power[base, exp]` with `base` containing a free variable
  (`Sqrt[x+1]` is for `simp_algebraic`'s Sqrt-multi-generator path, not
  for the rational-coefficient auto-detect).
- `Complex[..., ...]` or the symbol `I` (folded into the generator list
  as `Sqrt[-1]` already supported by `qa_resolve_extension`).
- More than `QA_AUTODETECT_MAX_GENS` distinct rational-base generators
  (combinatorial blowup in tower construction).

**Implementation skeleton** (~150 LoC):

```c
typedef struct {
    Expr* base;       /* rational base, copy-owned */
    int64_t q_lcm;    /* current LCM of exponent denominators */
    Expr* render;     /* surface form Power[base, 1/q_lcm], copy-owned */
} GenAccum;

static bool collect_one(const Expr* e, GenAccum** gens, size_t* n,
                        size_t* cap);

QATower* extension_autodetect(const Expr* e) {
    GenAccum* gens = NULL;
    size_t n = 0, cap = 0;
    if (!collect_one(e, &gens, &n, &cap) || n == 0) {
        free_gens(gens, n);
        return NULL;
    }
    Expr** renders = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) renders[i] = gens[i].render;
    QATower* tower = qa_resolve_extension_tower(renders, (int)n);
    free(renders);
    free_gens(gens, n);
    return tower;  /* NULL on tower-build failure */
}
```

### Phase B — Wire auto-detect into every Extension-aware builtin

For each builtin that currently accepts `Extension → α`:

```c
/* Pseudo-diff, applied to PolynomialGCD/LCM/Quotient/Remainder,
 * Cancel, Together, Apart, Factor. */
- const Expr* alpha = extract_extension_option(res, &argc);
+ Expr* alpha_owned = NULL;
+ const Expr* alpha = extract_extension_option(res, &argc);
+ bool alpha_is_automatic = (alpha != NULL
+     && alpha->type == EXPR_SYMBOL
+     && alpha->data.symbol == SYM_Automatic);
+ if (alpha == NULL || alpha_is_automatic) {
+     QATower* t = extension_autodetect(<the actual polynomial arg>);
+     if (t) {
+         alpha_owned = qa_tower_export_alpha(t);  /* gen list as
+                                                    one Expr* — single
+                                                    generator for n=1,
+                                                    List[...] for n>1 */
+         alpha = alpha_owned;
+         qa_tower_free(t);
+     }
+ }
  /* ... existing extension-aware path unchanged ... */
+ if (alpha_owned) expr_free(alpha_owned);
```

The single-generator case (`n == 1`) feeds the existing `Extension → α`
paths in `poly.c`/`facpoly.c`/`rat.c` unchanged.

For multi-generator (`n ≥ 2`) we have two options:

- **B.1 (tier-1):** Lower to a single primitive element γ via the existing
  `QATower`, then use the existing single-α path with α = γ.  Surface-form
  output replaces γ back with the user's generators using
  `tower->alpha_renders`.

- **B.2 (deferred):** Native Q(α₁, ..., αₙ)[x] arithmetic that doesn't
  collapse to γ.  Bigger refactor of `qaupoly.c`; not needed for the
  motivating example, which has effectively one generator (2^(1/6)).

Tier-1 ships B.1.  B.2 is a follow-up.

**Auto-detect must NOT fire when:**

- The input is already a single fraction with rational coefficients
  (avoid pointless field construction).  `extension_autodetect` returns
  NULL by construction in that case (no generators found).
- The user explicitly passed `Extension → None`.

### Phase C — Generalise `simp_algebraic` from q == 2 to general q

**File:** `src/simp.c`, around line 3610.

Replace the sign-flip rationalisation
(`alg_sigma_negate(den_r, gens[i])`) with general extended-Euclidean
rationalisation against the relation `g_i^(q_i) − u_i`:

```c
/* For each generator g_i with relation g_i^(q_i) = u_i, treat den_r as
 * a univariate polynomial in g_i with coefficients in Q[vars,
 * g_1, ..., g_{i-1}, g_{i+1}, ..., g_n].  Apply qaupoly-style extended
 * Euclidean to find s(g_i), t(g_i) with
 *     s · den_r + t · (g_i^(q_i) − u_i) = c
 * for some c free of g_i.  Then 1/den_r ≡ s(g_i) / c mod (g_i^(q_i) − u_i),
 * so num_r/den_r = (num_r · s) / c after reduction. */
static Expr* alg_rationalise_one(Expr* num_r, Expr* den_r,
                                 const char* gen, int64_t q,
                                 const Expr* u);
```

For `q == 2` this collapses to the existing
`alg_sigma_negate` ( `s = a − b·g_i` when `den_r = a + b·g_i` ).  For
`q == 3`, `s` is the product of two cyclotomic conjugates etc., computed
by extended Euclidean in the polynomial ring (no need for cyclotomic
ω, the Bezout coefficients are real-rational).

Then update the generator collector to accept any `q ≥ 2`:

```c
- if (q != 2) return false;            // simp.c:3610
+ if (q < 2 || q > ALG_MAX_DEG) return false;
```

with `ALG_MAX_DEG = 6` as a sane cap (degree-6 covers `2^(1/6)`, the
motivating example).

This step removes the cube-root bail and lets `simp_algebraic` handle
the entire motivating expression at the Simplify level (independent of
whether Together/Cancel's auto-detect path also fires).

### Phase D — Polynomial-function coverage

In addition to GCD/LCM/Quot/Rem/Cancel/Together/Apart/Factor, audit
these polynomial builtins for Extension threading:

| Builtin | File | Action |
|---|---|---|
| `PolynomialMod` | `poly.c` | Add `Extension` option + auto-detect path (mod-extension just reduces coefficients in Q(α)). |
| `Resultant`, `Discriminant` | `poly.c` | Same. Resultant over Q(α)[x] is just resultant of QAUPolys. |
| `FactorSquareFree` | `facpoly.c` | Same.  Square-free factorisation over Q(α) via `qa_factor_with_extension`'s square-free routine (already implemented internally). |
| `Variables`, `Coefficient`, `CoefficientList` | `poly.c` | No Extension semantics needed — these are syntactic.  Audit only to confirm they treat algebraic-number sub-expressions correctly when those are *atoms* of the working ring (i.e. don't list `2^(1/3)` as a "variable"). |
| `HornerForm` | `poly.c` | Coefficient-side only; no Extension semantics needed. |
| `IntegerPart`, `FractionalPart` of algebraic numbers | `core.c` | Not in scope (numeric only). |

`PolynomialQ[poly, Extension → α]` should also accept the option (return
True iff `poly ∈ Q(α)[vars]`); add as part of Phase D.

### Phase E — Tests

Three test surfaces:

**E.1 — Auto-detect unit tests** (`tests/test_extension_autodetect.c`,
new file).  No evaluation involved — direct exercise of
`extension_autodetect`:

```c
ASSERT_GEN_COUNT("Sqrt[2]",                      1);
ASSERT_GEN_COUNT("2^(1/3)",                       1);
ASSERT_GEN_COUNT("Sqrt[2] + Sqrt[3]",             2);
ASSERT_GEN_COUNT("2^(1/3) + 2^(2/3)",             1);  /* LCM merge */
ASSERT_GEN_COUNT("2^(1/3) + Sqrt[3]",             2);
ASSERT_GEN_COUNT("3 x^2 + 2 x + 1",               0);  /* no extension */
ASSERT_GEN_NULL ("Sqrt[x]");                            /* polynomial radicand */
ASSERT_GEN_NULL ("Sqrt[x^2 + 1]");                      /* same */
ASSERT_GEN_DEG  ("Power[2, 1/2] + Power[2, 1/3]", 6);  /* tower deg = lcm */
```

**E.2 — Builtin integration tests** (`tests/test_extension_auto_*.c`,
one file per builtin or one shared, TBD):

```c
/* Cancel */
EVAL_EQ("Cancel[(2^(1/3))^3]",                            "2");
EVAL_EQ("Cancel[2^(2/3) 2^(1/3)]",                        "2");
EVAL_EQ("Cancel[(x^3 - 2) / (x - 2^(1/3))]",
        "2^(2/3) + 2^(1/3) x + x^2");
EVAL_EQ("Cancel[(Sqrt[2] + Sqrt[3]) / (Sqrt[2] - Sqrt[3])]",
        "-5 - 2 Sqrt[6]");

/* Together — the motivating example.  Build the radical-laden RHS
 * directly to avoid depending on Integrate.  Then differentiate and
 * Together: must collapse to a x^2 / (x^3 + 2). */
EVAL_EQ("Together[D[<antiderivative-form>, x]]",
        "(a x^2) / (2 + x^3)");

/* PolynomialGCD */
EVAL_EQ("PolynomialGCD[x^3 - 2, x - 2^(1/3)]", "x - 2^(1/3)");

/* Factor */
EVAL_EQ("Factor[x^3 - 2]",
        "(x - 2^(1/3)) (2^(2/3) + 2^(1/3) x + x^2)");
/* (only when Extension auto-detects an external 2^(1/3) in context;
 * standalone x^3 - 2 still factors over Q.  Use a context-aware test:
 * Factor[(x^3 - 2)(y - 2^(1/3))] should split the cubic.) */

/* Apart */
EVAL_EQ("Apart[1 / (x^3 - 2), x]",
        "1/(3 2^(2/3) (x - 2^(1/3))) - (2^(1/3) + 2 x)/(6 (2^(2/3) + 2^(1/3) x + x^2))");
```

**E.3 — Simplify integration** (`tests/test_simp_algebraic_cuberoot.c`,
extension of existing `simp_algebraic` tests):

```c
EVAL_EQ("Simplify[D[Integrate[a x / (x^3 + 2), x], x]]",
        "(a x^2) / (x^3 + 2)");
EVAL_EQ("Simplify[(2^(1/3) + 1) (2^(2/3) - 2^(1/3) + 1)]", "3");
EVAL_EQ("Simplify[1 / (1 + 2^(1/3))]",
        "(1 - 2^(1/3) + 2^(2/3)) / 3");
```

**E.4 — Regression**: full existing test suite must pass.  The
auto-detect path must be NULL-returning (and therefore inert) on every
input that doesn't carry a `Power[rational, p/q]` or `Sqrt[rational]`
sub-expression.

### Phase F — Documentation and Information strings

- Update `docs/spec/builtins/algebra.md` (and `polynomial.md`,
  `rational.md`) with the new `Extension → Automatic` behaviour.
- Update each affected builtin's `symtab_set_docstring` in `info.c` to
  describe auto-detection.
- Add a `docs/spec/changelog/2026-05.md` entry summarising the work.

## Ordering and gating

```
A (extension_autodetect)                      ← foundational, no behaviour change
   │
   ├── B (wire into Cancel, Together, Apart, GCD, LCM, Quot/Rem)  ← user-visible
   │      │
   │      ├── B.1 single-generator tier-1
   │      └── B.2 multi-generator tier-1 (collapse to primitive elt)
   │
   ├── C (simp_algebraic: q > 2)                ← Simplify-side fallback
   │
   ├── D (PolynomialMod, Resultant, FactorSquareFree, PolynomialQ)
   │
   └── E (tests for A–D)
```

Phases A → B → E.1/E.2 unblock the motivating example via Together /
Cancel even without Simplify.  Phase C makes the Simplify-side path
generalise the same problem.  Phase D rounds out the polynomial-function
surface for parity with Mathematica's Extension → Automatic family.

## Risks and open questions

1. **Surface-form ambiguity.**  Mathilda normalises `Power[2, 2/3]` and
   `Power[2, 1/3] · Power[2, 1/3]` differently depending on which Times
   the canonicaliser sees first.  Auto-detect must collect generators
   from the *post-normalised* form to avoid duplicate enrollment.  Plan:
   run a one-shot `Together` (or `Expand` then `Together`) on the
   collected-generator working copy before calling
   `qa_resolve_extension_tower`.

2. **Detection vs. evaluation cost.**  The walk in `extension_autodetect`
   is O(|e|).  Each top-level `Together`/`Cancel` call now pays a single
   walk; on inputs with no algebraic generators (the common case) the
   walk is a fast no-op (early-return on first non-rational radical).
   Worst case: a deep `Hold`ed expression triggers the walk on every
   sub-call inside a tight loop.  Memoise on `simp_memo_wrap` style
   (`$ExtensionAutoDetect`) inside Simplify to amortise.

3. **`Sqrt[u]` with `u` itself algebraic** (e.g. the motivating
   `Sqrt[α/6]` after reduction).  Tier-1 strategy: detect that `u`
   reduces against the current generator set to `(rational) · θ^k` and
   bump the generator-degree LCM accordingly (e.g. `Sqrt[2^(1/3)]`
   bumps θ = 2^(1/3) of degree 3 to θ = 2^(1/6) of degree 6).  When `u`
   doesn't reduce that cleanly, leave the `Sqrt` as an opaque generator
   and rely on multi-generator tower construction.

4. **`Extension → Automatic` vs. `None` semantics.**  Mathematica
   distinguishes:
   - `Extension → None` (default in some functions): treat radicals as
     opaque.
   - `Extension → Automatic`: detect.
   Mathilda's current "default" is no-extension.  After this work, the
   semantics shift to "default == Automatic" because the auto-detect
   walk is cheap.  Document this clearly (it's a behaviour change for
   users who relied on the no-extension fallback).  Provide
   `Extension → None` as an opt-out.

5. **`simp_algebraic` complexity gate**.  Today it accepts results that
   are ≤ input complexity.  For q > 2, the Bezout rationalisation can
   produce an intermediate expression that is briefly larger than the
   input before Cancel kicks in.  Gate decision: gate on the *final*
   form after the trailing `Cancel`, not on intermediates.  The current
   code already does this — keep as-is.

## Out of scope (Tier 2+)

- Transcendental extensions (`Log`, `Exp`) for the Risch integrator.
- `RootOf`-style abstract algebraic numbers (Maple-style).
- Algebraic extensions of complete fields (`Q_p`, etc.).
- Multi-generator native Q(α₁, ..., αₙ)[x] arithmetic without primitive-
  element collapse (relevant only when collapse blows up the degree;
  rare in CAS workloads).
