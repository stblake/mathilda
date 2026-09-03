# Pre-existing FLINT-bridge / rat_canon per-call leak (56 bytes / 1 node)

**Status:** open. Discovered 2026-09-02 while landing `DSolve\`LieSymmetry` `bivariate`
(M10 L3). **Orthogonal to bivariate** — filed for a dedicated fix.

## Symptom
A `DSolve` operation whose `Simplify`/`Together` reduces a certain rational function
leaks exactly **one Expr node (56 bytes) per call**. Reproduces reliably (heavyweight)
via, e.g.:
```
DSolve`LieSymmetry[y'[x] == -1/x + y[x]/x + y[x]^2/x^3, y, x]   (* leaks 56 B/call *)
```
Measured N=1 vs N=8 under `valgrind --leak-check=full`: `definitely lost` grows by exactly
56 bytes per extra call (indirectly-lost flat).

## Precise localization (via -O0 builds of ratcanon.c + flint_bridge.c)
The orphaned node is a `Times` term of the **denominator** produced by the FLINT
round-trip:
```
expr_new_function            expr.c:92
fmpz_mpoly_to_expr           flint_bridge.c:2434   (Times term of the denominator)
flint_tower_reduce           flint_bridge.c:3596   (den = fmpz_mpoly_to_expr(denref(q)))
rat_canon_reduce             ratcanon.c:800        (reduced = flint_tower_reduce(frac,…))
rat_canon_normalize          ratcanon.c:878
builtin_together             rat.c:2120
… (evaluate) …
call_unary_copy              simp_util.c:41        (Simplify → simp_classify)
builtin_simplify             simp_builtins.c:794
ds_simplify                  (from dsolve_lie.c, during lie_first_integral)
```
`fmpz_mpoly_to_expr` and `flint_tower_reduce`'s `out` handling are each internally
balanced; the orphan is in the downstream refcount flow
`rat_canon_reduce:803-804` → `rat_canon_subst_back` → `rco_sign_normalize` →
`extract_num_den` (`src/rat.c:76`), where `reduced` is consumed and its num/den
re-extracted — one denominator term is dropped for this input shape.

## Why it's pre-existing and NOT introduced by bivariate
On `main` with the bivariate change stashed, the SAME ODE (which then **declines** —
no symmetry found) leaks **zero** bytes (13,440 flat at N=1 and N=8). The leak appears
only when a symmetry is **found and integrated** (`lie_first_integral` → `ds_simplify`),
i.e. through the shared `Simplify`→`Together`→FLINT path — not through bivariate's
determining-system code (verified balanced via three independent refactors: fresh-copy,
flat `lie_S_expr`, consume-fresh ownership — none changed the leak). Any `DSolve` /
`Simplify` / `Together` caller that reduces a rational of this shape can hit it.

## The hard part
**No standalone reproducer.** `Together` / `CoefficientList` / `Simplify` on the *parsed*
FullForm of the triggering expression do NOT leak — only the internally-C-built tree
(shared subnodes via `expr_copy`'s refcount bump, ground-stamped) triggers it. So the
trigger is an in-memory node property, not the printed structure.

## Suggested fix approach
1. Build a **C unit-test harness** that constructs the triggering expression via the
   `Expr` API the way the internal callers do (shared `expr_copy` subnodes), then calls
   `builtin_together` / `builtin_simplify` in a loop under valgrind → a fast standalone
   reproducer.
2. Bisect the `rat_canon_reduce` → `rat_canon_subst_back` → `rco_sign_normalize` →
   `extract_num_den` refcount flow to find the dropped denominator term.
3. **Regression-gate** `Together` / `Cancel` / `Simplify` / `Numerator` / `Denominator`
   (hot path — used everywhere) before/after; this is not a surface to edit blind.
