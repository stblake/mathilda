---
source: src/poly/facpoly.c
references:
  - "H. Zassenhaus, \"On Hensel factorization, I\", J. Number Theory 1969."
  - "D. G. Cantor, H. Zassenhaus, \"A new algorithm for factoring polynomials over finite fields\", Math. Comp. 1981."
  - "D. Y. Y. Yun, \"On square-free decomposition algorithms\", SYMSAC 1976."
  - "K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992)."
  - "B. M. Trager, \"Algebraic factoring and rational function integration\", SYMSAC 1976."
---
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
