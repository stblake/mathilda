# Factor

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Factor[poly] factors a polynomial over the integers.
Factor[poly, Extension -> alpha] factors over Q(alpha), where alpha is
Sqrt[c], c^(1/n) (rational c), or I.  Implements Trager's algebraic-
factoring algorithm via norm + sqfr_norm + alg_factor.
Factor[poly, Extension -> {alpha_1, ..., alpha_n}] factors over the
compositum Q(alpha_1, ..., alpha_n).  The tower is reduced to a single
primitive element gamma = alpha_1 + s_2 alpha_2 + ... via Trager's
primitive-element algorithm (Phase G6).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Factor[1 + 2x + x^2]
Out[1]= (1 + x)^2

In[2]:= Factor[x^10 - 1]
Out[2]= (-1 + x) (1 + x) (1 + x + x^2 + x^3 + x^4) (1 - x + x^2 - x^3 + x^4)

In[3]:= Factor[x^10 - y^10]
Out[3]= (x + y) (x - y) (x^4 - x^3 y + x^2 y^2 - x y^3 + y^4) (x^4 + x^3 y + x^2 y^2 + x y^3 + y^4)

In[4]:= Factor[2x^3 y - 2a^2 x y - 3a^2 x^2 + 3a^4]
Out[4]= (a + x) (a - x) (3 a^2 - 2 x y)

In[5]:= Factor[(x^3 + 2x^2)/(x^2 - 4y^2) - (x + 2)/(x^2 - 4y^2)]
Out[5]= ((-1 + x) (1 + x) (2 + x))/((x - 2 y) (x + 2 y))
```

## Implementation notes

**Algorithm.** `builtin_factor` lives in `src/poly/facpoly.c`, which is assembled from a set of `.inc` fragments (squarefree, Berlekamp–Zassenhaus univariate, bivariate/trivariate Hensel, the heuristic dispatcher, the memo cache, and the builtin entry point). The entry point (`facpoly_factor_builtin.inc`) is a cascade:

1. **Options / algebraic extensions.** A trailing `Extension -> α` (or `-> Automatic`, autodetected via `extension_autodetect`, or a `List` tower) routes to the algebraic-number factorer in `qafactor.c` (Trager's algorithm over `Q(α)`). The Simplify-scoped result memo (`factor_memo_*`) is bypassed in this branch.
2. **Memo + inexact handling.** A per-`Simplify` result cache is consulted; inexact coefficients are sent through `internal_rationalize_then_numericalize`.
3. **Threading.** `Factor` is `LISTABLE`; comparison/logic heads (`Less`, `Equal`, `Inequality`, `And`, …) are hand-threaded elementwise.
4. **Radical generators.** If the input contains a fractional-power sub-expression `u^{p/q}`, it substitutes `u → g^m` (m = lcm of the q's), factors as a polynomial in `g`, and back-substitutes.
5. **Rational normalisation.** `Together` → split into `Numerator`/`Denominator`; each is factored independently (in a direct call) or with a shared variable scope (inside Simplify).

The numerator/denominator are factored by variable count: **univariate** → `bz_factor_to_expr` (the Berlekamp–Zassenhaus pipeline, `facpoly_bz_uni.inc`); **bivariate** → a direct bivariate Hensel lift fast path (`factor_bivariate_via_hensel`), falling back to `FactorSquareFree` + `heuristic_factor`; **≥ 3 variables** → `FactorSquareFree` followed by `heuristic_factor`'s recursive per-piece dispatch.

The univariate **Berlekamp–Zassenhaus** core (`factor_zassenhaus`): take the primitive part, pick a prime `p` (starting at 13) under which the polynomial stays square-free, factor mod `p` by **distinct-degree** (`cz_ddf`) then **equal-degree** (`cz_edf`) Cantor–Zassenhaus splitting, **Hensel-lift** the modular factors to `p^k` (`multifactor_hensel_lift`/`hensel_lift`) with `p^k` chosen large, then **recombine** by brute-force subset enumeration over the lifted factors — each candidate true factor is reconstructed with leading-coefficient and content correction and accepted when `upoly_div_exact_z` divides the remaining polynomial exactly.

**Data structures.** The univariate path uses a packed `UPoly { int deg; int64_t* c; }` with `__int128_t`-guarded modular arithmetic (`upoly_mul_mod`, `upoly_div_rem_mod`, `upoly_gcd_mod`, `mod_inverse_int`) and `UPolyList` collections. Multivariate paths use `Expr*` polynomials plus the `MPoly`/`BPoly` representations from `mpoly.c`/`bpoly.c` for Hensel lifting. The algebraic path uses `QATower` from `qafactor.c`.

**Complexity / limits.** Modular coefficients are carried in `int64`/`int128`, so the chosen `p^k` is bounded below `10^15`; the Zassenhaus recombination is exponential in the number of modular factors (the classical worst case). The integer-arithmetic UPoly path caps degree at 10000. Multivariate factoring is heuristic (`heuristic_factor`) rather than a complete algorithm.

- `Listable`, `Protected`.
- **FLINT acceleration** (when built with FLINT): a univariate polynomial over Z is factored via `fmpz_poly_factor`, and a genuine multivariate polynomial over Q via `fmpq_mpoly_factor`, both from a guarded fast path at the top of `builtin_factor` (plain single-argument form only). This turns the exponential-in-variable-count classical multivariate factoriser — which hangs on inputs like `Factor[x^99 - y^99]` (> 20 s) — into a few-millisecond call. Each irreducible factor is normalised to a positive leading coefficient (highest total degree, deglex tie-break — Mathematica's convention), with the discarded sign folded into a separated rational content, so e.g. `Factor[y^2 - x^2]` → `-(x - y)(x + y)`. Inputs that are not a polynomial over Q (a denominator, fractional/symbolic exponent, or a non-polynomial head such as `Sqrt`/`Sin`) fall through to the classical path unchanged; `Extension`/`GaussianIntegers`/`Modulus` option forms travel their own branch and are untouched. See also the `` FLINT`Factor `` builtin.
- When given a rational expression, first resolves dependencies over `Together` before factoring.
- Uses exact root isolation (Rational Root Theorem limits) and binomial descents structured identically to Zassenhaus recombination, evaluating combinations exact and memory safe.
- Threads natively across lists, logic structures, and numeric groupings perfectly.
- Bivariate inputs whose leading coefficient (in some variable) is the constant `-1` are handled via Wang's leading-coefficient correction, Stage 1: the input is pre-negated to make it monic, the existing monic Hensel pipeline runs on `-P`, and the overall sign is absorbed into the highest-degree factor via `Expand`.  This unlocks inputs of shape `Factor[(1 - x^k)(x - y^m)]` (and similar non-monic cases with constant `±1` LC) that previously fell back to the legacy linear-trial-division loop.
- Bivariate inputs whose leading coefficient (in some variable) is a non-unit integer constant `a` (with `|a| > 1`) are handled via Wang's leading-coefficient correction, Stage 2: the monic substitution `Q(x, y) = a^(d-1) · P(x/a, y)` makes the lift's input monic in x with integer coefficients.  After lifting `Q = G_1 · ... · G_r` via the existing pipeline, the true factors of `P` are recovered as `F_i = G_i(a·x, y) / cont_Z(G_i(a·x, y))`, where the integer content collects exactly the share of `a^(d-1)` redistributed into G_i by the substitution.  Stages 1 and 2 compose, so inputs with negative non-unit LC (e.g. `lc_x(P) = -6`) also enter the structured pipeline.  This unlocks inputs like `Factor[Expand[(2x+3y)(3x+5y)]]`, `Factor[2 a^2 - 5 a b + 3 b^2]`, and three-factor non-monic forms whose LCs are constant in y (or constant in x).
- Bivariate inputs whose leading coefficient (in some variable) is a non-constant polynomial in the other variable are handled via Wang's leading-coefficient correction, Stage 3 (predicted-LC two-factor Hensel).  When `lc_x(P)(y) = A(y)`, Mathilda factors `A` over `Z[y]`, finds `α` with `A(α) = +1` so the squarefree univariate image `P(x, α)` factors into monic Z[x] pieces `u`, `v`, then enumerates distributions of `A`'s irreducible factors between two predicted leading coefficients `q_u, q_v` (with `q_u · q_v = A` and `q_u(α) = q_v(α) = +1`).  The Hensel iteration is modified so each `Δu` correction has its leading-x coefficient PINNED to the y^k coefficient of `q_u`, keeping `lc_x(U)(y) = q_u(y)` invariant across the lift.  This unlocks inputs like `Factor[Expand[(xy+1)(xy+2)]]`, `Factor[Expand[((y²+1)x+1)(x+3)]]`, `Factor[Expand[((y+1)x+1)((y+1)x+2)]]`.  MVP scope: r = 2 (two univariate factors), both monic, `|cont(A)| = 1`, and inputs with non-trivial monomial content fall through so `heuristic_factor`'s Phase 0 path produces the canonical fully-factored form.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- B. M. Trager, "Algebraic factoring and rational function integration", SYMSAC 1976 — the norm / sqfr_norm / alg_factor approach used for the Extension path.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 8 (polynomial factorization).
- H. Zassenhaus, "On Hensel factorization, I", J. Number Theory 1969.
- D. G. Cantor, H. Zassenhaus, "A new algorithm for factoring polynomials over finite fields", Math. Comp. 1981.
- D. Y. Y. Yun, "On square-free decomposition algorithms", SYMSAC 1976.
- K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992).
- B. M. Trager, "Algebraic factoring and rational function integration", SYMSAC 1976.
- Source: [`src/poly/facpoly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/facpoly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Factor[x^4 - 1]
Out[1]= (-1 + x) (1 + x) (1 + x^2)
```

```mathematica
In[1]:= Factor[6 x^2 + 7 x + 2]
Out[1]= (1 + 2 x) (2 + 3 x)
```

```mathematica
In[1]:= Factor[x^2 + 1, Extension -> I]
Out[1]= (-I + x) (I + x)
```

```mathematica
In[1]:= Factor[x^2 - 2, Extension -> Sqrt[2]]
Out[1]= (Sqrt[2] + x) (-Sqrt[2] + x)
```

```mathematica
In[1]:= Factor[x^10 - 1]
Out[1]= (-1 + x) (1 + x) (1 + x + x^2 + x^3 + x^4) (1 - x + x^2 - x^3 + x^4)
```

```mathematica
In[1]:= Factor[x^4 + 1, Extension -> Sqrt[2]]
Out[1]= (1 - Sqrt[2] x + x^2) (1 + Sqrt[2] x + x^2)
```

```mathematica
In[1]:= Factor[x^4 - 5 x^2 + 6, Extension -> {Sqrt[2], Sqrt[3]}]
Out[1]= (Sqrt[2] + x) (Sqrt[3] + x) (-Sqrt[2] + x) (-Sqrt[3] + x)
```

### Notes

Over the integers, `Factor` returns irreducible factors and keeps an
integer content split out front; `x^2 + 1` stays irreducible because it has
no rational roots. Supplying `Extension -> I`, `Extension -> Sqrt[c]`, or
`Extension -> c^(1/n)` factors over the corresponding algebraic number field
via Trager's algorithm, so `x^2 + 1` splits over `Q(I)` and `x^2 - 2` over
`Q(Sqrt[2])`. Factors are printed in canonical term order, which places the
constant term before the leading term (`-1 + x` rather than `x - 1`).

The cyclotomic factorisation of `x^10 - 1` recovers the degree-1, degree-4
cyclotomic polynomials Φ₁, Φ₂, Φ₅, Φ₁₀ over the integers. `x^4 + 1` (the 8th
cyclotomic, irreducible over Q) splits into two real quadratics over
`Q(Sqrt[2])`. Passing a list to `Extension` factors over the compositum: over
`Q(Sqrt[2], Sqrt[3])` the biquadratic `x^4 - 5 x^2 + 6` splits completely into
four linear factors, the tower being reduced to a single primitive element via
Trager's primitive-element algorithm.
