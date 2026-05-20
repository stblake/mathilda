# Risch–Norman (Parallel Risch / pmint) Implementation Plan

## Context

Mathilda currently has a strong rational-function integrator in `src/intrat.c`
(`Integrate\`BronsteinRational`), but anything containing a transcendental atom
(Sin, Cos, Exp, Log, Tan, ...) returns unevaluated. We will add the **parallel
Risch / Risch–Norman heuristic**, faithfully porting Bronstein's *Poor Man's
Integrator* (pmint) — the modern reformulation of Norman's 1976 heuristic — to
C as a new module `src/intrischnorman.c`. The new entry `Integrate\`RischNorman[f, x]`
will be wired in as a fall-through from `Integrate`, called after
`Integrate\`BronsteinRational` declines a non-rational integrand.

References (under `parallel_risch/`):

- `bronstein 2004 - parallel integration.pdf` — the algorithm we are porting.
- `pmint.maple` — the reference 99-line Maple source. Our C will be a near
  line-for-line port.
- `Geddes - Risch Norman Algorithm in MAPLE.pdf` — Norman's original 1976
  heuristic; the predecessor of pmint. Source of much of the test corpus.
- `Davenport - On the Parallel Risch Algorithm (I,II,III).pdf` — context.

The algorithm: rewrite trigs to half-angle Tan; introduce fresh variables for
each transcendental atom; lift `d/dx` to a polynomial vector field on those
variables (scaled by an `lcm` denominator `q`); construct a rational ansatz
with undetermined coefficients `_A_i` plus a sum `Σ _B_j Log[g_j]`; equate
`f − d(candidate)/q` to zero; extract the resulting linear system in
`{_A_i, _B_j}` over Q, solve via the existing `RowReduce`; on failure retry
over Q(I); back-substitute and return.

## Algorithm at a glance (pmint.maple line ↔ C helper)

| Step | Maple line | C helper |
|---|---|---|
| 1. Trig → half-angle Tan | 5 | `convert_to_tan` |
| 2. Collect indets, close under one diff | 6–8 | `collect_indets_closed` |
| 3. Build forward/reverse subst maps | 9–10 | `build_substitution_maps` |
| 4. Compute derivatives, lcm `q`, scale `l` | 11–15 | `build_vector_field` |
| 5. Darboux specials (Tan, Tanh, LambertW) | 16, 22–28 | `get_special` |
| 6. `splitFactor`, `deflation`, `enumerate_monoms` | 80–98, 69–78 | three siblings |
| 7. Candidate ansatz `(Σ _A_i m_i)/cden + Σ _B_j Log[g_j]` | 40–46 | `build_candidate` |
| 8. `numer(f − d(cand)/q)`, extract linear system | 54–62 | `try_integral` + `solve_linear_undet` |
| 9. K=0 then K=I retry | 49–50 | `try_integral` outer loop |
| 10. Back-substitute, zero free `_`-vars | 17–19 | `finalize_result` |

## File layout

- `src/intrischnorman.h` — public surface (~30 LOC): `void intrischnorman_init(void);`
  and the single builtin `Expr* builtin_rischnorman(Expr* res);`.
- `src/intrischnorman.c` — all helpers private (`static`). ~1600 LOC total.
- `src/integrate.c` — extend dispatcher (~25 LOC change). Add the
  `intrischnorman_init()` call after `intrat_init()`.
- `tests/test_intrischnorman.c` — new test binary, registered in
  `tests/CMakeLists.txt`. ~60 integrand corpus.
- `docs/spec/builtins/calculus.md` — document `Integrate\`RischNorman`.
- `docs/spec/changelog/2026-05.md` — add per-phase entries.
- `tasks/risch_norman.md` — task tracker per the project's `CLAUDE.md`
  "Task Management" convention.

No `sym_names.c` change needed: `Integrate\`RischNorman` is a string-keyed
package name (same convention as `Integrate\`BronsteinRational` — see
`intrat.c:3544`).

## Primitives reused (no new infrastructure)

All from Mathilda; no library code needs to be added.

| Need | C entry | File |
|---|---|---|
| Polynomial derivative | `D` builtin (evaluate Function) | (eval) |
| Together, Numerator, Denominator | `internal_together`, `internal_numerator`, `internal_denominator` | `internal.h:92,93,95` |
| Cancel, Expand, Collect | `internal_cancel`, `internal_expand`, `internal_collect` | `internal.h` |
| Coefficient (multivariate-by-list) | `internal_coefficientlist`, `internal_coefficient` | `internal.h:31,32` |
| Polynomial GCD / LCM | `internal_polynomialgcd` (LCM = ab/gcd) | `internal.h:33` |
| Factor over Z | `internal_factor` | `internal.h:53` |
| Factor over Q(I) | `Factor[p, Extension -> I]` (Trager) | `facpoly.c:2796-2879` |
| Replace / ReplaceAll | `internal_replace_all` | `internal.h:160` |
| Variables[expr] | `internal_variables` | `internal.h:30` |
| **Linear system solve** | `builtin_rowreduce` (Bareiss-fraction-free) | `linalg.c:516` |

The single tool we don't have is `Solve` for an arbitrary linear system in
named unknowns, but `RowReduce` on an augmented matrix is strictly more
powerful and is the right primitive — see Phase 4 below.

## Six implementation phases

Each phase is independently compilable, independently testable, and ships a
visible feature increment. Test integrands are accumulated via
`assert_integral_correct` (mirroring the existing helper in
`tests/test_intrat.c:56-72`).

### Phase 1 — Plumbing and skeleton (~150 LOC, ~3h)

**Deliverable.** `src/intrischnorman.{c,h}` with a stub
`builtin_rischnorman` that always returns `NULL`. `intrischnorman_init()`
wired into `integrate_init()` (`integrate.c:171`) right after `intrat_init()`.
Test scaffolding in place, dispatcher unchanged.

**Files.** `src/intrischnorman.c`, `src/intrischnorman.h` (new); `src/integrate.c`,
`src/integrate.h`; `tests/test_intrischnorman.c` (new); `tests/CMakeLists.txt`;
`tasks/risch_norman.md` (new); `docs/spec/changelog/2026-05.md`.

**Acceptance.** `Integrate\`RischNorman[Sin[x], x]` evaluates to itself
(unevaluated), the symbol carries `ATTR_PROTECTED`, all existing
`tests/` binaries still pass.

### Phase 2 — Trig → Tan, indet collection, substitution maps (~350 LOC, ~5h)

**Deliverable.** Three private helpers:

```c
static Expr* convert_to_tan(Expr* f, Expr* x);
static int   collect_indets_closed(Expr* ff, Expr* x,
                                   Expr*** out_si, size_t* out_n);
static int   build_substitution_maps(Expr** si, size_t n,
                                     PMSubMap* out_lin, PMSubMap* out_lout,
                                     Expr*** out_vars);
```

`convert_to_tan` is a recursive rewrite over the trig family and hyperbolic
siblings — Sin[u] → 2T/(1+T²) with T = Tan[u/2]; similarly Cos, Tan, Cot,
Sec, Csc; hyperbolics via Tanh[u/2]. FreeQ short-circuit on subtrees not
involving `x`. Depth-cap at 4; returns NULL on giveup.

`collect_indets_closed` walks the tree gathering Tan/Tanh/Log/Exp/LambertW/
ArcTan/ArcTanh atoms with `D[atom, x] ≠ 0`, then one closure pass: differentiate
each atom and harvest *its* function subexpressions. Dedupe by `expr_hash` then
`expr_eq`. The lone `x` symbol is the first indet.

`build_substitution_maps` allocates `pmint$v_1, pmint$v_2, ...` fresh
symbols via the existing `intern_symbol` helper. `PMSubMap` is the
C-struct introduced below.

**Data structure.**

```c
typedef struct { Expr* lhs; Expr* rhs; } PMSubEntry;
typedef struct { PMSubEntry* items; size_t n; } PMSubMap;
```

`lin` maps original-term → fresh-var; `lout` is the reverse. Materialised on
demand as `List[Rule[lhs, rhs], ...]` for `ReplaceAll` calls.

**Acceptance.** Unit tests for each helper via temporary
`Integrate\`Helpers\`PM*` testable surfaces (gated `ATTR_READPROTECTED`):
`PMConvertToTan[Sin[x]+Cos[x], x]` returns the expected `2T/(1+T²) + (1-T²)/(1+T²)`
form; cardinality of indets is 4 for `Exp[Sin[x]]`; subst maps round-trip
(`subs[lout, subs[lin, ff]] == ff` for several test inputs).

### Phase 3 — Vector field, splitFactor, deflation, enumerate_monoms (~450 LOC, ~5h)

**Deliverable.** Port pmint.maple lines 11–15, 30–35, 69–98:

```c
static int   build_vector_field(Expr** li, size_t n, PMSubMap* lin,
                                Expr* x,
                                Expr*** out_l, Expr** out_q);
static Expr* apply_d(Expr* f,
                     Expr** vars, Expr** l, size_t n);
static int   split_factor(Expr* p, Expr** vars, Expr** l, size_t n,
                          Expr** out_s, Expr** out_h);
static Expr* deflation(Expr* p, Expr** vars, Expr** l, size_t n);
static int   enumerate_monomials(Expr** vars, size_t nv, int total_degree,
                                 Expr*** out_monoms, size_t* out_n);
```

`apply_d(f) = Σ l[k] · D[f, vars[k]]`, expanded and put over a common
denominator. `split_factor` faithfully ports the multivariate recursion from
Maple's lines 80–90 — pick a var with non-zero derivation, extract
content+primitive part wrt it via `CoefficientList`, compute
`s = PolynomialGCD(q, apply_d(q)) / PolynomialGCD(q, D[q, x_k])`, recurse on
`q/s`. `deflation` is the simpler sibling at lines 92–98. Monomial enumeration
is the obvious recursive product (with a `PMINT_MAX_MONOMIALS = 5000` hard
cap; on overflow return error code).

**Acceptance.** Toy differential field K(x, V) with V = Exp[x] (so
`apply_d(V) = V`): assert `split_factor[V·(V+1)^2, dx] == {V, (V+1)^2}`;
`deflation[(V+1)^2(V−1), dx] == (V+1)(V−1)`; `enumerate_monomials[{x,V}, 3]`
returns 10 distinct monomials.

### Phase 4 — Candidate ansatz, linear-system extraction, solve (~700 LOC, ~6h, the bulk)

**Deliverable.** The closing pipeline — minus log-candidates — that lets us
integrate the *Q-rational-only* test cases:

```c
static Expr* build_candidate(Expr** vars, size_t nv, int degree_bound,
                             Expr* cden,
                             char*** out_A_names, size_t* out_nA);
static int   try_integral(Expr* f, Expr** vars, Expr** l, size_t nv,
                          Expr* q,
                          Expr* cand_main,
                          char** A_names, size_t nA,
                          Expr** candlogs, size_t nlog,
                          char*** B_names, size_t* out_nB,
                          PMSubMap* lout,
                          int K_use_extension_I,
                          Expr** out_result);
static int   solve_linear_undet(Expr* equation_numer,
                                Expr** vars, size_t nvars,
                                Expr** unknowns, size_t nunk,
                                PMSubMap* out_solution,
                                int* out_status);
```

**Linear-system extraction algorithm** (the crux):

1. `equation_numer = Expand[Numerator[Together[f − apply_d(cand)/q]]]`.
   This polynomial is *polynomial in `vars`* and *linear in `unknowns`*.
2. Call `internal_coefficientlist(equation_numer, List[vars...])` once —
   this is the bulk extractor `coeff_list_rec` (`poly.c:2942`) and is the
   single fastest path. Yields an `nv`-deep nested list whose leaves are
   the per-monomial coefficients (each a linear form in the unknowns).
3. Recursively flatten the nested list, tracking the exponent vector.
   Skip leaves where `is_zero_poly(leaf)` is true. Each surviving leaf
   becomes one row of the linear system.
4. For each surviving leaf and each unknown `_A_i`, call
   `internal_coefficient(leaf, _A_i)` — must be a pure rational
   (`is_number` assertion catches non-linear bugs). The constant term
   (RHS) is the leaf with all unknowns set to 0; sign-flip to move it
   to the RHS column.
5. Assemble `List[List[q_11, ..., q_1m, b_1], ..., List[..., b_R]]` as
   an Expr matrix. Leaves stay in `EXPR_RATIONAL` / `EXPR_INTEGER` /
   `EXPR_BIGINT` form so `is_zero_poly` works.
6. Call `builtin_rowreduce` on the augmented matrix. The result is
   reduced row-echelon form, fraction-free Bareiss-style
   (`linalg.c:537-589`), so no symbolic blow-up.
7. Decode RREF:
   - Row `[0 ... 0 | c]` with `c ≠ 0` → **infeasible** (status -1; caller
     should try K=I or give up).
   - Pivot in column `j < nunk` → `unknowns[j] -> matrix[r][nunk]`
     (after zeroing all free unknowns; per pmint line 18-19, free
     `_`-vars are pinned to 0).
   - No pivot in column `j` → `unknowns[j] -> 0`.

**Phase 4 corpus** (all Q-rational, no log-candidates yet):
`Exp[x]`, `Exp[a x]`, `x Exp[x]`, `x^2 Exp[x]`, `Exp[2x] Sin[x]`,
`Sin[x] Exp[x]`, `Cos[x] Exp[x]`, `1/(1+x^2)` (cross-check with
BronsteinRational), `Log[x]`, `x Log[x]`.

**Acceptance.** All ten integrands pass `assert_integral_correct`.

### Phase 5 — Log-candidate sums, getSpecial, K=I retry (~450 LOC, ~6h)

**Deliverable.** Logarithm part of the candidate and the K=I retry:

```c
static int   my_factors(Expr* p, int over_Qi,
                        Expr*** out_factors, size_t* out_n);
static int   get_special(Expr* atom, PMSubMap* lin,
                         Expr** out_darboux, int* out_is_special);
```

`my_factors(p, over_Qi=0)` calls `internal_factor({p}, 1)`; for
`over_Qi=1` it builds `Factor[p, Extension -> I]` and dispatches through
`builtin_factor` (`facpoly.c:2796`). Univariate Q(I) factoring via Trager
(`qa_factor_with_extension`) handles the typical case where each candlog
is univariate in one fresh variable. Multivariate Q(I) factoring is
flagged unsupported (returns NULL → caller skips this candidate).

`get_special` is the Maple `getSpecial` (lines 22–28): Tan[u] → Darboux
`1 + Tan[u]^2` (non-integral); Tanh[u] → `1 + Tanh[u]` and `1 - Tanh[u]`
(non-integral); LambertW[u] → `LambertW[u]` (integral, flag true).

`try_integral` is extended to assemble `candlog = my_factors(l1) ∪
my_factors(l2) ∪ my_factors(l3) ∪ {getSpecial integral-flagged polys}`,
build `candidate = cand_main + Σ _B_j Log[candlog_j]`, then run the same
extract-and-solve from Phase 4. If K=0 returns infeasible, retry with
`K_use_extension_I = 1`.

**Phase 5 corpus** (log-bearing and Q(I) cases):
`Tan[x]`, `Cot[x]`, `Sin[x]/Cos[x]`, `Tan[x]^2`, `Sec[x]^2`,
`1/(1+Exp[x])`, `Exp[x]/(1+Exp[x])`, `Log[x]^2`, `Log[x]^3`,
`1/(x Log[x])`, `(1+Log[x])/(x Log[x])`, `1/(1+Sin[x]^2)` (needs Q(I)),
`x/(1+Exp[x])` (Geddes Table II), `Sin[x]/(1+Cos[x])`,
`Log[Log[x]]`, plus the LambertW cases from Bronstein's paper if
LambertW is available.

**Acceptance.** All 15+ log-bearing integrands pass; K=I retry
demonstrably activates on `1/(1+Sin[x]^2)`.

### Phase 6 — Dispatcher hook, post-hoc verification, edge cases (~300 LOC, ~4h)

**Deliverable.** Final wiring.

Extend `src/integrate.c`'s `builtin_integrate` (currently bails at lines
129–133 when the integrand is neither polynomial nor rational): replace
the bail with a fall-through call to `Integrate\`RischNorman[f, x]`.
Detect the unevaluated package head exactly as the existing
`BronsteinRational` detection at lines 147–154 and bubble back as
`Integrate[f, x]`.

Post-hoc verifier (recommended, defensive):

```c
static bool verify_result(Expr* result, Expr* f, Expr* x);
```

After producing a candidate antiderivative, compute
`Cancel[Together[Expand[D[result, x] − f]]]`; if it isn't structurally 0,
return NULL. This catches any silent bugs in coefficient extraction.

Wall-clock budget: `PMINT_WALL_CLOCK_SEC = 30` per integrand, enforced
via `setitimer`. Document `PMINT_MAX_MONOMIALS = 5000` and
`PMINT_WALL_CLOCK_SEC` as `#define`s at the top of `intrischnorman.c`.

Polish: docstrings via `symtab_set_docstring`, `ATTR_PROTECTED`,
`docs/spec/builtins/calculus.md` entry, `docs/spec/changelog/2026-05.md`
phase-by-phase summaries, `Mathilda_spec.md` overview untouched (no new
top-level category).

**Phase 6 corpus** (combined: full Geddes Table II + Bronstein 2004 +
known-fail cases):
- Geddes Table II (~20 cases): includes the rows for
  `Exp[1/(Log[x]+1)] (2 Log[x] + 1)/(Log[x]+1)^2`, `1/(x Log[x])`,
  `x/(1+Exp[x])`, `Tan[a x]`, `Tan[a x]^2`, `Cos[2x] Exp[x]`,
  `Exp[x] Log[x^2-1]/Log[x]^2`, plus the known-fail cases
  `Exp[x^2]`, `1/Log[x]`, `Exp[x^-10]`.
- Bronstein 2004 Bessel-J example (deferred behind a `BesselJ`
  feature gate if not yet implemented).
- Norman 1976 / Davenport corpus (additional ~15 integrands).
- Cross-validation: every test_intrat.c integrand also passes when
  dispatched through the RischNorman path (regression check).

**Acceptance.** Full ~60-integrand corpus passes; known-fail cases
return NULL cleanly within budget; existing `tests/` suite has no
regressions; `Integrate[f, x]` end-to-end works for both rational and
transcendental integrands.

## Memory ownership (across all phases)

Per the `feedback_builtin_ownership` memory: `builtin_rischnorman` must
NOT call `expr_free(res)` — the evaluator owns `res` and frees it
itself. All internal helpers follow the convention:

- **Borrowed input Exprs** (the caller retains ownership).
- **Owned output Exprs** (caller must eventually free).
- `PMSubMap` items are owned by the constructor; `pm_sub_map_free()` is
  the destructor.
- `Expr**` arrays returned by `collect_indets_closed`,
  `build_vector_field`, `enumerate_monomials`: caller frees each Expr,
  then `free()` the array itself.

Match the existing `intrat.c` style exactly.

## Performance characteristics

- **Candidate space**: with the typical degree bound `dg ≈ 5` and
  `nv ≈ 3-4` fresh vars, ~50-150 monomials; the linear system stays
  under 200×200.
- **RowReduce cost**: Bareiss fraction-free at `linalg.c:537-589`; per-
  pivot cost is dominated by `mpq` arithmetic on rational leaves.
  Practical bound: a 200×200 matrix in well under 1 second.
- **Coefficient extraction**: one `Expand`, one `CoefficientList` per
  `try_integral` invocation. Linear in the size of the expanded
  numerator.
- **Hard caps**: `PMINT_MAX_MONOMIALS = 5000`,
  `PMINT_WALL_CLOCK_SEC = 30`, both checked inside
  `enumerate_monomials` and the solver. Exceeding either returns NULL
  (informational message `Integrate\`RischNorman::overlarge`).

## Risk register

1. **Linear-system extraction silently drops cross-terms.** Mitigation:
   aggressive `Expand` + `Together` at every interface; unit-test
   extraction in isolation (Phase 4) via a testable
   `Integrate\`Helpers\`PMExtractLinearSystem[poly, vars, unknowns]`
   surface before wiring into the solver.

2. **`split_factor` multivariate recursion misbehaves.** Mitigation:
   line-by-line port of `pmint.maple:80-90`; defensive `Expand` before
   each `content` step; toy-field acceptance tests in Phase 3 before
   relying on it in Phase 4.

3. **Trig → Tan rewrite doesn't reach nested radicals like
   `Sqrt[1+Sin[x]^2]`.** Mitigation: scoped to non-algebraic trigs in
   Phase 1; documented limitation; depth-cap at 4 returns NULL on giveup.

4. **`Factor[_, Extension -> I]` is univariate-only.** Mitigation:
   multivariate candlog factoring over Q(I) is documented as out of
   scope. The typical pmint case where each candlog is univariate in a
   fresh var works; multivariate Q(I) becomes a future follow-up.

5. **`feedback_builtin_ownership` already cost us before** — confirm
   no `expr_free(res)` in `builtin_rischnorman` via a code review pass
   before each phase merges.

## Verification

Every Phase 4+ test uses `assert_integral_correct(integrand)`, a verbatim
copy of `tests/test_intrat.c:56-72`, asserting
`Cancel[Together[Expand[D[Integrate[f, x], x] − f]]] === 0`. This is the
universal correctness predicate.

Per-phase verification commands:

```sh
# Build (auto-picks up new src/*.c via wildcard in src/makefile:35-37)
make -j$(nproc)

# Run new test binary
cd tests && mkdir -p build && cd build && cmake .. && make -j$(nproc)
./test_intrischnorman_tests

# Full regression
for t in *_tests; do ./$t; done

# Memory check on the new path (Linux)
valgrind --leak-check=full --error-exitcode=1 ./test_intrischnorman_tests
```

End-to-end manual check:

```mathematica
Integrate[Tan[x], x]                 (* → -Log[Cos[x]] or equivalent *)
Integrate[1/(1 + Exp[x]), x]         (* → x - Log[1 + Exp[x]] *)
Integrate[x/(1 + Exp[x]), x]         (* Geddes Table II *)
Integrate[1/(1 + Sin[x]^2), x]       (* Exercises K=I retry *)
Integrate[Exp[x^2], x]               (* Returns unevaluated cleanly *)
```

## Critical files

- `src/intrischnorman.c` (new — ~1600 LOC, bulk of the algorithm)
- `src/intrischnorman.h` (new — public surface)
- `src/integrate.c` — dispatcher hook at lines 129–133 and 171
- `src/intrat.c:3540-3656` — template for `install()` pattern,
  package-context name registration, fraction-free conventions
- `src/linalg.c:516-589` — `RowReduce` (Bareiss-fraction-free, the
  linear-system solver)
- `src/facpoly.c:2796-2879` — `Factor[poly, Extension -> I]` (Trager,
  for K=I retry)
- `src/poly.c:2942-2966` — `coeff_list_rec` / `builtin_coefficientlist`
  (multivariate coefficient extraction)
- `src/internal.h` — all `internal_*` helpers we wire through
- `parallel_risch/pmint.maple` — line-by-line port reference
- `tests/test_intrischnorman.c` (new — ~60-integrand corpus)
- `tests/test_intrat.c:56-72` — template for `assert_integral_correct`
- `tasks/risch_norman.md` (new — per-phase task tracker)
