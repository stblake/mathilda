---
title: Simplify / Collect bug list — 2026-05-10
date_started: 2026-05-10
status: planning
---

# Categorization of 40 reported cases

PASS already (16 cases): 2, 3, 4, 5, 7, 8, 9, 12, 15, 16, 17, 18, 19, 26, 27, 40
PASS — FullSimplify-only territory (2 cases): 21, 22
PASS — different form, same leaf count (1 case): 6 (`16·2^a·8^b` vs `2^(4+a+3b)`)

Real bugs: 21 cases, grouped into 5 categories below.

---

## Category A — Power / Exp / Log simplification

| # | Input | Expected | Got |
|---|---|---|---|
| 1  | `Simplify[(-1)^(2n) - 1, n ∈ ℤ⁺]` | `0` | unchanged |
| 10 | `Simplify[Exp[y(Log a + Log b)] - (a b)^y, a>0 ∧ b>0]` | `0` | unchanged |
| 11 | `Simplify[Exp[3(Log a + Log b)] - a³b³]` | `0` | unchanged |

Notes:
- (1) needs Simplify to consume the assumption `n ∈ PositiveIntegers` and rewrite `(-1)^(2n) → 1` via even-integer detection.
- (10) needs assumption-aware `Log[a]+Log[b] → Log[a*b]` (positive reals), then `Exp[y Log(ab)] → (ab)^y`.
- (11) needs no assumptions: integer exponent makes `(a*b)^3 = a^3 b^3` unconditional. The missing step is recognizing `Exp[k·(Log a + Log b)]` → `Exp[k Log a]·Exp[k Log b]` → `a^k b^k` for integer `k`.

(Issue 6 is intentionally excluded — same leaf count, not a defect.)

---

## Category B — Plus auto-cancellation when grouping differs

`Collect[expr, var] - expr` should evaluate to 0 without needing `Simplify`.
Mathematica's `Plus` distributes/cancels grouped terms automatically; Mathilda does not.

| # | Pattern |
|---|---|
| 28 | `Collect[x+x², x] - (x+x²)` |
| 29 | `Collect[x²+yx, x] - (xy+x²)` |
| 30 | `Collect[2x²+yx²+3xy, x] - (x²(2+y)+3xy)` |
| 31 | `Collect[2x²+yx²+3xy, y] - (2x²+y(x²+3x))` |
| 32 | `Collect[(1+y+x)^4, x] - <expanded form>` |
| 33 | `Collect[xExp[x]+ySin[x]+2Sin[x]+3x, x] - ((3+Exp[x])x+(2+y)Sin[x])` |

Common mechanism: `Plus` should detect `S + (-1)·S` for `S` a Plus subterm (or expand on demand and cancel). Probably a single canonicalization step in `plus.c` that splits `Times[c, Plus[...]]` summands into distributed terms before final combine.

---

## Category C — Collect can't handle compound monomials / sum patterns

| # | Input | Expected |
|---|---|---|
| 34 | `Collect[a(c+s)+b(c+s), c+s]` (`c=Cos[x], s=Sin[x]`) | `(a+b)(c+s)` |
| 35 | `Collect[xy + a·xy, xy]` | `(1+a) x y` |
| 36 | `Collect[1 + xy + a·xy, xy]` | `1 + (1+a) x y` |
| 37 | `Collect[a·x·f[x] + b·x·f[x], x·f[x]]` | `(a+b) x f[x]` |
| 38 | `Collect[a·x·Log[x] + b·x·Log[x], x·Log[x]]` | `(a+b) x Log[x]` |
| 39 | `Collect[xyz + a·xyz, xyz]` | `(1+a) x y z` |

Note: case 40 (Collect on `x^c`) already passes, so the existing pattern handles single `Power` subexpressions but not `Times` / `Plus` patterns. Need to extend `Collect` to accept compound second-argument keys.

---

## Category D — Hangs / non-termination on radical denominators

These three each peg CPU and do not return. Worst severity.

| # | Sketch |
|---|---|
| 13 | `Simplify[ 1/(√2 a + √3 + √5 + √7) - <rationalized form> ]` |
| 14 | `Simplify[ 1/(√2 a + √2 b + √3 + √7) - <rationalized form> ]` (also raises `Power::infy` then hangs) |
| 20 | `Simplify[ 1/((3+3√2) x (3x - 3√y)) - <rationalized form> ]` |

Likely cause: a nested-radical normalization loop in `Simplify` (or the Bezout/conjugate step in `Together` over a multi-radical extension) recurses without a termination guard. We need to bound iteration count and bail out, returning the best-so-far form.

---

## Category E — Nested-radical algebraic simplification

| # | Input | Expected |
|---|---|---|
| 23 | `Simplify[-(2+2√2+4·2^(1/4))/(1+2^(3/4)+3·2^(1/4)+3√2)]` | `-2/(1+2^(1/4))` (LeafCount 11) — Mathilda gets `-(2+4·2^(1/4)+2^(3/2))/(1+2^(1/4))^3` (LeafCount 25) |
| 24 | Long √5-nested-radical sum | `ComplexInfinity` (denominator is exactly 0) |
| 25 | `1/√(√7/7 + 2√2 + 3√3 + 5√5) + 3 - <rationalized e2>` | `0` |

These are the hardest — they require either factoring over a radical extension or a real-numeric verification gate before declaring "no simplification possible".

---

## Proposed order of attack

1. **Category D (hangs)** first — these are blocker-level severity. Even if we just add a hard iteration cap that returns "best-so-far", the user gets a tractable session.
2. **Category B (Plus auto-cancel)** — likely a single mechanism in `plus.c`; fixes 6 cases at once and is foundational for many downstream simplifications.
3. **Category C (Collect compound keys)** — localized to `expand.c` (`Collect`), 6 cases, mostly mechanical.
4. **Category A (Log/Exp/Power-with-assumptions)** — needs assumption propagation; case 11 may be doable without assumptions.
5. **Category E (nested radicals)** — research-level; tackle last with one focused subagent each.

---

## Per-category review checkpoints

After each category, run the full test corpus (`make test`) and the
intrat corpus to check for regressions, then update this file with
results before moving on.

---

## Progress log

### Final status (2026-05-10)

| Category | Cases | Status |
|---|---|---|
| **A** Power/Log assumption-aware simplification | 1, 10, 11 | **PASS** (all 3 return 0) |
| **B** Plus auto-cancel `(a+b) − (a+b)` | 28, 29, 30, 31, 32, 33 | **PASS** (all 6 return 0) |
| **C** Collect on compound monomials/sums | 34, 35, 36, 37, 38, 39 | **PASS** (all 6 collect correctly) |
| **D** Hangs on radical-denominator rationalisation | 13, 14, 20 | **Hangs eliminated** at Together/Cancel level. Simplify-to-zero still slow on case 13 (>5min); foundational fix is in place — `simp_algebraic` (≤4 surds) handles the rationalisation when reachable. |
| **E** Nested-radical algebraic simplification | 23, 24, 25 | **Known limitation** — requires factoring over algebraic extensions (Q[2^(1/4)]), nested-radical Sqrt-product folding for ComplexInfinity detection, large-extension Simplify. Out of scope for this iteration. |

Net regressions: **none**. The pre-existing 4 test failures on `main` (linalg float format, list Min[9,2] bug, parfrac Apart factoring, polymod print format) remain unchanged.

### Category D — partial: hangs killed at Together/Cancel (2026-05-10)

Three independent fixes in `src/poly.c` and `src/rat.c`:

1. **`exact_poly_div` rational/Gaussian gate** (poly.c:1191) — at
   `var_count == 0`, only return the symbolic `Times[A, Power[B,-1]]`
   fallback when both operands lie in a field (Q or Q[i]). Previously
   it returned that form for any pair of atoms, including non-field
   atoms like `Sqrt[2]`. Returning a "quotient" that doesn't actually
   exist in `Q[Sqrt[2], Sqrt[3], …]` propagated symbolic
   `Power[Plus, -1]` into intermediate polynomials, which sent
   `cancel_recursive`'s `PolynomialGCD` call into multivariate Euclid
   over a rational input — the case-13 hang.

2. **`together_recursive` strict-quotient combine** (rat.c:700) — Plus
   branch now uses `cancel_exact_div_strict` (returns NULL on
   non-exact) instead of the symbolic-fallback wrapper. When any
   summand's quotient `lcm_den/dens[i]` isn't exact in `Q[vars]`, the
   combine aborts cleanly and returns the original Plus, leaving the
   algebraic structure for `simp_algebraic` (≤4 surds) downstream.

3. **`cancel_recursive` embedded-rational guard** (rat.c:255) — before
   computing `PolynomialGCD[num, den]`, walk both for any
   `Power[<compound>, <neg integer>]` subterm. If found, the input
   isn't a clean polynomial fraction; skip the GCD step and return the
   input unchanged.

4. **`poly_gcd_internal` size-and-iteration budget** (poly.c:1685) —
   pseudo-remainder Euclidean loop now caps at `max(input_size, 2000)`
   leaves and 50 iterations. On exhaustion returns the content GCD —
   sound (a divisor of the true GCD) and lets the cancel pipeline
   proceed without coefficient explosion. Mirrors the existing
   subresultant-PRS budget pattern (poly.c:3340).

Results:
- Case 13 `Together[e1-e2]` completes in seconds, LeafCount 336.
- Case 14 `Together[e1-e2]` completes, LeafCount 234.
- Case 20 `Together[1/e1 - e2]` completes, LeafCount 90.
- Pre-existing `tests/test_rat.c:66` expectation updated to the
  cleaner sign-normalised form (semantically equivalent).
- No new test regressions: 5 pre-existing failures on main remain
  unchanged (linalg float format, list Min[9,2] bug, parfrac Apart
  factoring, polymod print format, zupoly false-positive grep).

Remaining for Cat D: confirm `Simplify[…]` actually reduces these to
0 (rather than just terminating). If `simp_algebraic` is reached it
should fold them; if not, follow-up work needed in Simplify's seed
pipeline. Verifying via long-running test in progress.

### Category B — `Plus` distributes `Times[-1, Plus[…]]` (2026-05-10)

Single edit in `src/plus.c` `builtin_plus`: at the start of plus
evaluation, walk the args; for any arg matching exactly
`Times[-1, Plus[t_1, …, t_m]]`, replace it with m fresh
`Times[-1, t_i]` args and re-flatten. Mathematica's `a + b - (a + b)`
returns 0 because it distributes the leading -1 before grouping;
Mathilda now matches that. The distribution is gated on the literal
coefficient -1 — `a + 2 (b + c)` stays unexpanded, matching MMA's
behaviour and avoiding accidentally over-eager Plus expansion.

Cases 28-33 all now return 0 via the standard arithmetic auto-eval
(no Simplify needed).

### Category C — `Collect` on compound monomial / Plus keys (2026-05-10)

`collect_internal` in `src/poly.c` rewritten to dispatch on
`var_bp.count`:

- **Single-base key** (atoms `x`, powers `Power[x, k]`, Plus keys
  `Plus[c, s]`): group each term by `term_exp / var_exp` (the
  rational/symbolic exponent at which the key's base appears). This
  preserves the existing semantics for `Collect[…, x]` over
  fractional/symbolic exponents like `Sqrt[x]` and `x^c`.

- **Multi-factor key** (`Times[x, y, z]`): use the existing `get_k`
  helper to find the integer multiplicity at which `var` divides each
  term, residualise, and group.

For Plus keys, expansion is skipped (`expand_patt` would distribute
the Plus across Times factors and destroy the structure we want to
collect on). Cases 34-39 all correct; case 40 (`x^c` Power key) still
passes via the single-base path.

### Category A — Power / Log simplification (2026-05-10)

Three small edits in `src/simp.c`:

1. **`apply_logexp_rules` no longer requires a non-NULL ctx** (line
   6392). NULL ctx is treated as an empty context — the
   positivity-aware rewrites simply don't fire (`prov_pos` returns
   false), but the unconditional Exp-distribute rule below does its
   job.

2. **`transform_can_fire("LogExpRules", …)`** no longer demands
   `ctx_has_facts`. Same reasoning — gates on Log/Power presence only.

3. **`logexp_top_rewrite`** gains two unconditional rules:
   - `Power[E, Plus[t_1, …, t_n]] -> Times[Power[E, t_1], …, Power[E, t_n]]`,
     applied after `Expand` of the exponent so `Times[c, Plus[…]]` is
     surfaced first. Combined with Mathilda's existing
     `Power[E, c·Log[u]] -> u^c`, this collapses
     `Exp[k (Log[a] + Log[b]))` to `a^k · b^k`. Sound for any
     exponent (E^(x+y) = E^x · E^y always).
   - `Power[-1, k]` with `prov_even(ctx, k)` -> 1.

4. **`prov_int`** extended to recognise `Element[x, PositiveIntegers]`
   and the other ±-tagged integer subdomains. **`prov_even`** now
   propagates through `Times`: a product is even when at least one
   factor is even and every factor is integer. Together these turn
   `(-1)^(2 n)` with `n ∈ PositiveIntegers` into the integer 1.
