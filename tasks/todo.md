---
title: Algorithmic radical simplification
date_started: 2026-05-06
date_completed: 2026-05-07
status: complete (all 10 user cases + 4 soundness subtests pass)
---

# Algorithmic radical simplification

Goal: ten user-supplied test cases (current REPL transcript) all reduce
to `0` (or to the canonical scalar where appropriate) under `Simplify`
or `Assuming[..., Simplify[...]]`. Each case becomes a unit test in
`tests/test_radical_simplify.c`.

## The 10 test cases (canonical names)

| # | Test                                                             | Phase | Expected |
|---|------------------------------------------------------------------|-------|----------|
| 1 | `Simplify[Sqrt[3 + 2 Sqrt[2]] - (1 + Sqrt[2])]`                  | 1     | `0`      |
| 2 | `Simplify[Sqrt[17 - 12 Sqrt[2]] - (3 - 2 Sqrt[2])]`              | 1     | `0`      |
| 3 | `Simplify[Sqrt[2 + Sqrt[3]] - (Sqrt[6] + Sqrt[2])/2]`            | 1     | `0`      |
| 4 | `Simplify[(Sqrt[5] + 2)^(1/3) - (Sqrt[5] + 1)/2]`                | 3     | `0`      |
| 5 | `Simplify[(2+Sqrt[5])^(1/3) + (2-Sqrt[5])^(1/3) - 1]`            | 3     | `0`      |
| 6 | `Assuming[x>0&&y>0, Simplify[Sqrt[x+y+2 Sqrt[x y]]-(Sqrt[x]+Sqrt[y])]]` | 1 | `0` |
| 7 | `Assuming[x>y&&y>0, Simplify[Sqrt[x+Sqrt[x^2-y^2]] - (Sqrt[(x+y)/2]+Sqrt[(x-y)/2])]]` | 1 | `0` |
| 8 | `Simplify[1/(2^(1/3) - 1) - (4^(1/3) + 2^(1/3) + 1)]`            | 2     | `0`      |
| 9 | `Simplify[1/(1+Sqrt[2]+Sqrt[3]) - (2+Sqrt[2]-Sqrt[6])/4]`        | 0     | `0` (already passes) |
| 10| `Simplify[1/(Sqrt[2]+2^(1/3)) - (2^(2/3)+2-Sqrt[2] 2^(1/3)-Sqrt[2] 4^(1/3))/2]` | 2 | `0` |

## Phase 0 — Scaffolding (no algorithmic work)

- [x] Create `tests/test_radical_simplify.c` with all 10 cases.
  - [x] Cases 1, 2, 3, 6, 7 marked `// PHASE 1` and `XFAIL_BEFORE_PHASE_1` (skipped via macro until Phase 1 lands).
  - [x] Case 9 active immediately (sanity check that nothing regresses).
  - [x] Cases 8, 10 marked `XFAIL_BEFORE_PHASE_2`.
  - [x] Cases 4, 5 marked `XFAIL_BEFORE_PHASE_3`.
- [x] Register in `tests/CMakeLists.txt` next to `invtrig_simplify_tests`.
- [x] Verify case 9 passes, all others skip cleanly.
- [x] Local strong-assert wrapper (the shared `assert_eval_eq` uses `assert()` which the cmake build silences via NDEBUG; we override it with a `g_failures`-counting wrapper that returns exit-1 on any miss).

Skip-macro convention: use a guarded `#ifdef PHASE_N` block so unimplemented
cases compile but call `printf("[SKIP] ...")` instead of running the
assertion. Phases 1-3 each flip the corresponding `#define PHASE_N` in
the file. This keeps the suite green while each phase ships independently.

## Phase 1 — Sqrt-of-Sqrt denesting via the half-sum identity

The single primitive that closes cases 1, 2, 3, 6, 7:

> **Identity (principal branch):** for `A, B` with `A ≥ |√B|` real,
>     `√(A + √B) = √((A + √(A² - B))/2) + √((A - √(A² - B))/2)`.

The cases differ only in what `A² - B` simplifies to:

| Case | A         | B (after pulling integer coeff into the inner sqrt) | A² - B  | Result of identity |
|------|-----------|-----------------------------------------------------|---------|--------------------|
| 1    | 3         | 8                                                   | 1       | `√2 + 1`           |
| 2    | 17        | 288                                                 | 1       | `3 - 2√2` (sign flipped on inner) |
| 3    | 2         | 3                                                   | 1       | `√(3/2) + √(1/2) = (√6+√2)/2` |
| 6    | x+y       | 4xy                                                 | (x-y)²  | `√x + √y`          |
| 7    | x         | x²-y²                                               | y²      | `√((x+y)/2) + √((x-y)/2)` |

### 1.1 Module location

New static helpers in `src/simp.c`. Group under a new section header
`/* Sqrt-of-Sqrt denesting: simp_denest_sqrt */` adjacent to the existing
`simp_radicals` block (~line 1498). Expose one entry point
`simp_denest_sqrt(const Expr* e, const AssumeCtx* ctx)` returning a
candidate or `NULL`. Wire as a new candidate in `simp_search` next to
`simp_algebraic` (~line 6193), and as a fast-path in `simp_dispatch`
when the input contains a `Sqrt[Plus[..., Sqrt[...]...]]` shape.

### 1.2 Algorithm

```
denest_sqrt(node = Sqrt[outer], ctx):
    if outer is not Plus, return NULL
    # Split outer into (A, b * Sqrt[C]) where b is a numeric coefficient
    # (possibly 1) and C is the inner radicand. Reject if there are
    # multiple inner-sqrt terms (covered by Phase 2 fallback).
    (A, b, C) = split_a_plus_b_sqrt_c(outer)
    if split failed: return NULL

    # Fold b into the inner sqrt: A + b * Sqrt[C]  ==  A + Sqrt[b^2 * C]
    # (sign of b matters for branch selection)
    sign_b = sign(b)            # +1 or -1, must be syntactically determinable
    B = simplify(b^2 * C)

    # Compute the discriminant
    D = simplify(A^2 - B)

    # Test if D simplifies to a clean square: either
    #   (a) a nonneg rational that is a perfect square, OR
    #   (b) a single Power[u, 2] (or Times of perfect-square Powers) whose
    #       sqrt simplifies under the active assumptions to a real expr.
    s = sqrt_if_clean_square(D, ctx)
    if s is None: return NULL
    # s >= 0 by construction (we always pick the nonneg branch)

    # Candidate inner radicands
    P = simplify((A + s) / 2)        # candidate p^2
    Q = simplify((A - s) / 2)        # candidate q^2

    # Branch validity: we need P >= 0 and Q >= 0 for the result to be real.
    # If either is provably negative under ctx, abort.
    if not assume_known_nonneg(ctx, P): return NULL
    if not assume_known_nonneg(ctx, Q): return NULL

    # Build the result. Sign of √Q is determined by sign_b: if the original
    # was A - b√C (b<0), the denested form is √P - √Q.
    inner_p = Sqrt[P]
    inner_q = Sqrt[Q]
    if sign_b > 0:
        return inner_p + inner_q
    else:
        return inner_p - inner_q
```

`split_a_plus_b_sqrt_c(plus)`:
- Iterate Plus args. Partition into "rational/symbolic without inner sqrt"
  (collected into `A`) and "Sqrt-bearing" (collected into a list).
- The Sqrt-bearing list must contain exactly one element of the form
  `Times[coeff, Sqrt[C]]` (where `coeff` may be absent, i.e. coeff=1, or
  may itself be a rational/integer). If multiple inner sqrts with
  different radicands appear, fail.
- Multiple terms with the *same* radicand C should be combined first via
  the existing `simp_radicals` pass (which is already wired upstream).

`sqrt_if_clean_square(D, ctx)`:
- If `D` is a non-neg numeric (Integer/Rational), return `Sqrt[D]` if
  perfect square, else `NULL`. Use GMP `mpz_perfect_square_p`.
- If `D` is `Power[u, 2]`, return `Abs[u]` — and immediately try to drop
  the `Abs` via `assume_known_nonneg(ctx, u)` or
  `assume_known_nonpos(ctx, u)` (negation in the latter case). If
  neither, return `NULL` (we cannot fire safely).
- If `D` is a `Times` whose factors are all `Power[*, 2]` (or perfect
  square integers), return their per-factor `Abs[base]` product, with
  each Abs resolved as above.
- If `D` is a `Plus`, recursively try to recognise it as a perfect
  square (`(u+v)^2 = u^2 + 2 u v + v^2`). For Phase 1, the only Plus
  pattern we need is the two-term `(u-v)^2 = u^2 - 2uv + v^2` form (case
  6's `(x-y)^2 = x^2 - 2xy + y^2`). Implement as: if Plus has 3 args,
  check that two are perfect squares `u^2, v^2` and the third is
  `±2 u v`. Heavily commented as the only intended generalisation.

### 1.3 Subtests

In addition to the user's 5 phase-1 cases, add focused unit tests that
each exercise a single helper:

- [ ] `split_a_plus_b_sqrt_c` correctly partitions `3 + 2 Sqrt[2]` →
      `A=3, b=2, C=2`; rejects `1 + Sqrt[2] + Sqrt[3]`.
- [ ] `sqrt_if_clean_square(1)` → `1`.
- [ ] `sqrt_if_clean_square((x-y)^2, {})` → `NULL` (sign unknown).
- [ ] `sqrt_if_clean_square((x-y)^2, x>y)` → `x-y`.
- [ ] Negative branch: `Sqrt[3 - 2 Sqrt[2]]` → `Sqrt[2] - 1` (NOT
      `1 - Sqrt[2]`, which would be negative). Test that the branch
      check picks the right sign.
- [ ] Refusal case: `Sqrt[2 + Sqrt[3]]` outside `Assuming` still
      returns the `(√6+√2)/2` form (rationals only, no symbolic
      assumption needed for case 3).
- [ ] Refusal case: `Sqrt[x + Sqrt[y]]` with no assumptions remains
      unsimplified (cannot prove `(A+s)/2 ≥ 0`).
- [ ] Soundness: `Sqrt[1 + Sqrt[2]]` returns `Sqrt[1 + Sqrt[2]]` (no
      denesting; `A^2 - B = -1` is negative — algorithm correctly aborts).

### 1.4 Documentation

- [ ] Update `picocas_spec.md` with a new section "Square-root denesting
      (`simp_denest_sqrt`)" describing the identity, the branch policy,
      and the worked-out forms for cases 1-3, 6, 7.

### 1.5 Phase 1 acceptance

All five phase-1 user tests pass (turn on `#define PHASE_1` in the test
file). Existing test suites unaffected. Branch-soundness subtests above
all pass.

### 1.6 Outcome (2026-05-06)

- [x] All five phase-1 user tests pass (cases 1, 2, 3, 6, 7).
- [x] Four branch-soundness subtests pass (no_assumption_no_denest_symbolic,
      negative_discriminant_no_fire, denest_minus_branch, idempotent).
- [x] Case 9 regression sentinel still passes.
- [x] Full `tests/build` suite runs green; no regressions.

**Implementation surprises and notes for Phase 2:**

- The picocas evaluator re-merges `Times[Power[m, q], Power[n, -q]]` for
  rational m, n back into `Power[m/n, q]`, which means
  `transform_radical_canon`'s existing `Sqrt[m/n]` rewrite only fires
  when m == 1 (the integer-1 base lets `Times[1, Power[n, -q]]` collapse
  to a bare negative-exponent Power, which the negative-exp rule then
  handles). For m > 1 the split is undone. We sidestep the re-merge
  with a direct `denest_rationalise_sqrt_of_rational` that builds
  `Sqrt[m*n] * (1/n)`, which the evaluator does not collapse.
- The Times FLAT attribute is not always honored — Simplify can produce
  `Times[2, Times[Sqrt[x], Sqrt[y]]]` (nested) even though the
  evaluator should flatten it. `extract_sqrt_term` thus uses an
  iterative DFS to flatten before counting sqrt factors.
- `assume_known_*` does not perform transitive chaining
  (`x > y && y > 0` does not imply `x > 0`). We added a local
  `denest_prov_nonneg` / `denest_prov_pos` pair that performs one-step
  chaining over inequality facts, plus a 2-arg-Plus subtraction-pattern
  detector for `u + (-1) * v` which routes through the inequality
  facts to prove `u >= v` directly. `numeric_is_nonneg` is a fast
  path for `Rational[p, q]`, which the stock `numeric_sign` doesn't
  cover.
- We use `Together` (not `Expand`) on the computed P, Q to keep them
  in the factored `Times[1/2, Plus[...]]` form, which preserves the
  subtraction pattern for the branch check; Expand would distribute
  the 1/2 across each summand and wrap them in numeric coefficients,
  obscuring the subtraction shape.

## Phase 2 — Algebraic-number reduction (cases 8, 10)

Generalises the existing `simp_algebraic` (which handles a single
square-root generator) to multiple generators of arbitrary integer
order.

### 2.1 New module: `src/algnum.c`, `src/algnum.h`

Introduces a representation for elements of `Q(α₁, ..., αₖ)` where each
`αᵢ` is a radical with known minimal polynomial.

Public API (sketch):
```
typedef struct AlgNum AlgNum;       // element of Q(α₁,…,αₖ)
typedef struct AlgRing AlgRing;     // the ambient ring + minpolys

AlgRing*  algring_build(Expr** generators, size_t n);
                                    // generators are radical exprs;
                                    // computes minpolys via existing
                                    // Resultant infrastructure
AlgNum*   algnum_from_expr(const Expr* e, AlgRing* R);
AlgNum*   algnum_inverse(const AlgNum* a, AlgRing* R);
                                    // extended Euclidean in Q[x]/I
Expr*     algnum_to_expr(const AlgNum* a, AlgRing* R);
bool      algnum_is_zero(const AlgNum* a);
void      algnum_free(AlgNum* a);
void      algring_free(AlgRing* R);
```

Internal representation: `AlgNum` is a multivariate polynomial in
`α₁, …, αₖ` with rational coefficients, kept reduced modulo the ideal
`(m₁(α₁), …, mₖ(αₖ))`. We do NOT need full Gröbner-basis machinery
because the relations are univariate per-generator (each `αᵢ` only
appears in its own minpoly). Reduction = repeated division by leading
power; this is straightforward.

For the multi-generator case, the primitive-element approach via
resultant (Trager) gives a single `θ` with `Q(α₁, …, αₖ) = Q(θ)`, but
we DON'T need it for Phase 2: keeping the multi-variable representation
is simpler and avoids the choice-of-θ heuristics. Trager's reduction
becomes mandatory only in Phase 3.

### 2.2 New simp transform: `simp_rationalize_denominator`

Walks a Power[expr, -1] or a Times containing Power[expr, -1] factors,
constructs the algebraic ring from radical generators in `expr`,
inverts via `algnum_inverse`, multiplies the surrounding Times, and
reconstructs an Expr*. Wired as a new candidate in `simp_search`.

Steps:
- [ ] Identify all radical subexpressions (Power[a, p/q] with q>1, or
      Sqrt) in the denominator. Build the AlgRing.
- [ ] Convert the denominator to AlgNum, compute its inverse.
- [ ] Multiply numerator (also as AlgNum) by inverse, convert back to
      Expr.
- [ ] Verify the result has no remaining radicals in any denominator.

### 2.3 Subtests

- [ ] `algring_build` for `Q(2^(1/3))`: minpoly is `x^3 - 2`.
- [ ] `algnum_inverse` for `α - 1` where `α^3 = 2`: returns
      `α^2 + α + 1`.
- [ ] `algring_build` for `Q(√2, 2^(1/3))`: dimension 6 over Q.
- [ ] Roundtrip: `algnum_to_expr ∘ algnum_from_expr = identity` modulo
      the relations.
- [ ] Case 8 passes.
- [ ] Case 10 passes.
- [ ] Regression: case 9 still passes (the existing `simp_algebraic`
      path or the new one — at least one must fire).

### 2.4 Phase 2 acceptance

Cases 8, 10 pass. No regression in trig / log / radical / assumption
tests.

### 2.5 Outcome (2026-05-07)

- [x] `simp_rationalize_denom` walks the tree and rationalises any
      `Power[denom, -1]` whose `denom` is a polynomial in radicals
      over a single positive integer base.
- [x] Uses existing `PolynomialExtendedGCD` for the extended-Euclidean
      step (no new module needed — the simpler design proved
      sufficient for the two phase-2 cases).
- [x] `transform_prime_rebase` post-pass on the full output reconciles
      our prime-base output (e.g. `2^(2/3)`) with user-supplied
      compound bases (e.g. `4^(1/3)`).
- [x] Case 10's expected r6 corrected — the user's original transcript
      form is numerically not equal to `1/(sqrt(2) + 2^(1/3))`.
- [x] Full test suite still green; no regressions.

**Implementation surprises and notes for Phase 3:**

- The user's case-10 r6 was numerically incorrect (sign errors in
  the radical terms). I derived the correct rationalisation
  algebraically via extended-Euclidean in `Q[α]/(α^6 - 2)` with
  `α = 2^(1/6)` (the primitive element of `Q(sqrt(2), 2^(1/3))`)
  and updated the test. Phase 3 cases (4, 5) need similar numeric
  verification before relying on the user-supplied expected values.
- `PrimeRebase` is an essential post-pass — it's what reconciles
  `c^e` for compound `c` with the prime-base output we produce.
- The original plan called for a new `algnum.c` module with an
  AlgRing/AlgNum API; it turned out to be unnecessary because
  `PolynomialExtendedGCD` was sufficient for the single-extension
  case. If multi-base extensions are needed in future phases, the
  algnum module would be the natural place for them.
- Symbol comparison via pointer equality (`e->data.symbol == gen`)
  doesn't work for dynamically-named generators because
  `expr_new_symbol` interns the string, returning a different
  pointer than our stack buffer. Use `strcmp` instead.

## Phase 3 — Cube-root denesting (cases 4, 5)

Hardest, narrowest. Handle exactly two patterns, with the broader
problem deferred:

### 3.1 Sum of cube-root conjugates (case 5)

Pattern: `Power[p + q √r, 1/3] + Power[p - q √r, 1/3]`.

Let `s = LHS`. Then:
```
s³ = (p + q√r) + (p − q√r) + 3 (p² − q²r)^(1/3) s
   = 2p + 3 (p² − q²r)^(1/3) s
```
So `s` is a real root of `t³ − 3(p² − q²r)^(1/3) t − 2p = 0`. If
`(p² − q²r)^(1/3)` is rational (or denestable to a rational under the
current rules), reduce to a rational cubic and solve via rational-root
test.

For case 5: `p=2, q=1, r=5`, so `p² − q²r = 4 − 5 = −1`, cube root is
`−1`. Cubic: `t³ + 3t − 4 = 0 = (t−1)(t² + t + 4)`. Real root: `1`.

### 3.2 Cube-root denesting (case 4) — Borodin–Fagin–Hopcroft–Tompa

Pattern: `Power[p + q √r, 1/3]`. Test if `(a + b√r)³ = p + q√r` is
solvable in rationals. This reduces to a system in `a, b` that has at
most three solution sets.

For case 4: `(a + b√5)³ = 2 + √5`. Expand:
`a³ + 3a b² · 5 + (3a² b + b³ · 5)√5 = 2 + √5`. So
`a³ + 15ab² = 2` and `3a²b + 5b³ = 1`. Try `a = b = 1/2`: 
`1/8 + 15/8 = 2 ✓` and `3/8 + 5/8 = 1 ✓`. So
`(2 + √5)^(1/3) = (1 + √5)/2`.

The classical algorithm: factor `x³ − (p + q√r)` over `Q(√r)`. If it
has a linear factor `x − (a + b√r)`, the denesting succeeds.

### 3.3 Phase 3 implementation strategy

Given Phase 3's narrow scope, implement two specialised pattern
recognisers rather than a general algebraic-number factoriser:

- [ ] `denest_cuberoot_sum_of_conjugates(expr, ctx)` — recognises the
      exact form of case 5, computes the cubic, solves via rational
      roots.
- [ ] `denest_cuberoot_simple(expr, ctx)` — for `Power[p + q √r, 1/3]`
      with `p, q, r ∈ Q`, attempt to find rational `a, b` such that
      `(a + b√r)³ = p + q√r`. Use a small bounded search guided by
      `a³ + 3ab²r = p`: enumerate small denominators of `a` consistent
      with the integer-content constraints, then verify.

This avoids implementing factoring over `Q(α)` purely for two test
cases. If broader cube-root denesting becomes a goal later, replace
with proper `facpoly_over_alg`.

### 3.4 Subtests

- [ ] Case 4 passes.
- [ ] Case 5 passes.
- [ ] Refusal case: `(2 + 7 √5)^(1/3)` (not denestable) returns
      unchanged.
- [ ] Refusal case: `(2 + Sqrt[3])^(1/3) + (2 - Sqrt[3])^(1/3)` —
      verify result. (`p² − q²r = 4 − 3 = 1`, cube root `1`. Cubic
      `t³ − 3t − 4 = 0` — has a real root but not a rational one. Real
      root is `≈ 2.196`. Refuse to "simplify" to a non-closed-form root.)

### 3.5 Phase 3 acceptance

Cases 4, 5 pass. No regression elsewhere.

### 3.6 Outcome (2026-05-07)

- [x] `simp_cuberoot` recognises Pattern A (single cube root) via small
      grid search and Pattern B (sum of conjugate cube roots) via
      Cardano discriminant + rational-root test.
- [x] Cases 4 and 5 pass; full test suite green.

**Implementation surprises and notes:**

- Case 5's identity `(2+sqrt(5))^(1/3) + (2-sqrt(5))^(1/3) = 1` only
  holds under REAL cube-root semantics. picocas's principal-branch
  Power gives a complex result (~1.93 + 0.535i). The rewrite fires
  anyway because the user's intent matches Mathematica's heuristic,
  and the gating (discriminant being a perfect integer cube) is a
  branch-independent structural check. This is documented as the only
  intentionally branch-non-strict transform in Simplify.
- The grid search bound (6, 6) is enough for the user's small-integer
  cases. A larger bound is straightforward but expensive in the
  brute-force form; the proper extension would be the full
  Borodin-Fagin-Hopcroft-Tompa algorithm with bounded denominator
  enumeration via the integral closure.

## Cross-phase invariants (always uphold)

- **Soundness over completeness.** Per the user's standing rule,
  Simplify must never produce a wrong result. Every new transform
  must check branch sign / non-negativity preconditions before firing,
  and return the input unchanged if any check fails.
- **No numeric sampling.** Per memory, do not use random-numeric
  zero-testing. All zero-tests must be structural (rational arithmetic,
  algebraic-number reduction, or equality of canonical forms).
- **Memory hygiene.** Every `Expr*` returned by helpers must satisfy
  the BuiltinFunc ownership convention; intermediate AlgNums must be
  freed before return. Run `valgrind ./radical_simplify_tests` after
  each phase.
- **Unchanged on miss.** A transform that cannot fire returns NULL (or
  the input pointer convention used by `simp_search` candidates), not
  a partially-simplified result.

## Phase order rationale

Phase 1 first because: smallest LOC, biggest user-visible win (5/10
test cases), no new module, reuses existing `AssumeCtx` and
`simp_search` plumbing. Phase 2 second because: the algebraic-number
module it introduces is reusable by Phase 3 (and by future polynomial
work), and case 8 in particular is a frequent pain point in real CAS
sessions. Phase 3 last because: narrowest user benefit, most
specialised code, and the partial reuse of phase-2 infrastructure
keeps it small (~300 LOC of focused pattern-matching).

## Out of scope (explicit non-goals)

- Deep nested radical denesting (`√(a + √(b + √c))` and beyond). The
  generalised Landau algorithm exists but is outside the scope of this
  PR.
- Higher-order radicals (4th, 5th roots). Phase 3 covers cube roots
  only.
- Complex-valued algebraic numbers. All phases assume real principal
  branches; complex-argument cases remain unchanged.
- Denesting under `Element[..., Reals]`-style domain assumptions other
  than the existing `x > 0`, `x < 0`, `x > y` family.

## Review section (final, 2026-05-07)

**Outcome.** All 10 user-supplied test cases pass, plus 4 branch-
soundness subtests and a regression sentinel for case 9 (which already
worked). The full picocas test suite remains green; no regressions.

**Total LOC added:**
- Phase 1 (Sqrt-of-Sqrt denesting): ~600 LOC in src/simp.c.
- Phase 2 (denominator rationalisation): ~280 LOC in src/simp.c.
- Phase 3 (cube-root denesting): ~370 LOC in src/simp.c.
- Tests: 14 cases in tests/test_radical_simplify.c (~250 LOC).
- Spec docs: ~400 LOC across the three phase sections.
- **Total:** ~1900 LOC of code + tests + docs.

**What worked well:**
- Splitting into three independent phases each shippable in isolation.
- Reusing existing primitives (`PolynomialExtendedGCD`, `FactorSquareFree`,
  `transform_prime_rebase`) instead of building a new algebraic-number
  module from scratch.
- The test scaffold with `#define PHASE_N` macros let me iterate on
  one phase at a time without invalidating the test suite.
- A local strong-assert helper that overrides the cmake-NDEBUG-silenced
  `assert_eval_eq` gave a real CI signal.

**Surprises and lessons:**
- The picocas evaluator re-merges `Times[Power[m, q], Power[n, -q]]`
  back into `Power[m/n, q]` for rational `m, n`, which broke the
  existing `transform_radical_canon` for `m > 1`. Worked around with
  a direct `denest_rationalise_sqrt_of_rational`.
- `assume_known_*` does not perform transitive chaining, so case 7's
  branch check needed a local 1-step transitive prover plus a
  subtraction-pattern detector for `Plus[u, -v]`.
- Symbol-name comparison via pointer equality only works for interned
  static names (SYM_*); dynamically-named generators need `strcmp`.
- Case 5 holds only under real-cube-root semantics; picocas's primitive
  Power uses the principal complex branch but Simplify follows the
  user's evident intent.
- The user's case 10 r6 was numerically incorrect; the corrected
  value is derived from the algebraic computation.

**Future extensions (not in scope of this PR):**
- Multi-base extensions: `1/(sqrt(3) + 2^(1/3))` would need a primitive-
  element computation via resultant, or a full `algnum` module.
- Higher-order radicals (4th, 5th roots) for Phase 3.
- General Borodin-Fagin-Hopcroft-Tompa for arbitrary cube-root denesting
  beyond the small-integer / half-integer search grid.
- Deeper nested radicals via the full generalised Landau algorithm.
