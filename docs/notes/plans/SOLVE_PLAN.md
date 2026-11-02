# Plan: `Solve` — router + initial polynomial-equality implementation

## Context

Mathilda needs `Solve` to be feature-complete with respect to Mathematica's
core equation-solving capability. Several existing modules already plan for
it as a dependency:

- `src/intrat.c:1393`, `:1532`, `:1613`, `:2895` — `Solve / ToRadicals` is
  referenced as a future replacement for hardcoded rational-integration
  heuristics.
- `src/root.c:14` — the held-form `Root[]` / `RootSum[]` infrastructure
  was split out specifically so "Solve, Reduce, factorisation over algebraic
  extensions" could emit those nodes.

So `Solve` has a defined home in the system; this plan establishes the
router and ships the first specialist (polynomial equality in one variable).

User-confirmed design points:
- Default domain is `Complexes` (return all complex roots).
- Domain `Reals` returns only real roots via **discriminant-aware** per-degree
  branching (not generate-then-filter).
- Binomial `a*x^n + b == 0` (n ≥ 3) returns all n complex roots.
- **n-quadratics** `a*x^(2n) + b*x^n + c == 0` (n ≥ 1) — biquadratic
  (n=2), bicubic-in-x^3 (n=3), etc. — return all 2n radical solutions
  by default, independent of `Cubics` / `Quartics` options. They are not
  routed through Cardano/Ferrari because the substitution u = x^n
  reduces them to a quadratic in u, after which two binomial sub-solves
  recover all the original roots in closed form.
- Irreducible factor of degree ≥ 5: emit one held `Root[Function[t, p], k]`
  per root. Solvable factors emit closed-form rules; this is decided
  **per factor**, not globally.

---

## Architecture

```
Solve[...]                      src/solve.c           (router)
  └── dispatches to one of:
      ├── Solve`SolvePolynomialEquality[...]   src/solvepoly.c   (THIS PR)
      ├── Solve`SolveLinearSystem[...]         src/solvelinsys.c (future)
      └── …                                     (future)
```

All specialists are also reachable directly by their context-qualified
names (precedent: `Integrate\`RischNorman`, `Integrate\`BronsteinRational`).
A single `solve_init()` registers everything and is called once from
`core_init()` in `src/core.c`.

---

## New files

| File | Purpose |
|------|---------|
| `src/solve.h` / `src/solve.c` | `builtin_solve` router. Argument classification, option parsing, dispatch. Owns `solve_init()`. |
| `src/solvepoly.h` / `src/solvepoly.c` | `builtin_solve_polynomial_equality` and all per-degree sub-routines (linear, quadratic, binomial, cubic, quartic, high-degree Root emission). |
| `tests/test_solve.c` | Unit tests for both router and specialist. |

(Defer `src/solvelinsys.c` — listed only so the router's dispatch table is
forward-compatible; no code in this PR.)

---

## Symbols to add

In `src/sym_names.h` (declarations, alphabetical) and `src/sym_names.c`
(matching `const char* SYM_… = NULL;` + `intern_symbol(...)` lines in
`sym_names_init()`):

- `SYM_Solve` → `"Solve"`
- `SYM_SolvePolynomialEquality` → `"Solve\`SolvePolynomialEquality"`
  (the backtick is part of the string; `intern_symbol` is opaque to
  punctuation, and `src/context.c:54-197` resolves the qualified form
  transparently)
- `SYM_Cubics` → `"Cubics"`
- `SYM_Quartics` → `"Quartics"`
- `SYM_GeneratedParameters` → `"GeneratedParameters"`
- `SYM_VerifySolutions` → `"VerifySolutions"`
- `SYM_Assumptions` → `"Assumptions"`  (verify not already present —
  `SYM_AssumptionRules` exists at sym_names.h:48 but `SYM_Assumptions` does not)
- `SYM_InverseFunctions` → `"InverseFunctions"`

Already present and reusable: `SYM_Equal`, `SYM_Unequal`, `SYM_Greater`,
`SYM_GreaterEqual`, `SYM_Less`, `SYM_LessEqual`, `SYM_And`, `SYM_Or`,
`SYM_List`, `SYM_Rule`, `SYM_RuleDelayed`, `SYM_Plus`, `SYM_Times`,
`SYM_Power`, `SYM_Rational`, `SYM_Complex`, `SYM_I`, `SYM_Pi`, `SYM_E`,
`SYM_Function`, `SYM_Slot`, `SYM_Root`, `SYM_Reals`, `SYM_Complexes`,
`SYM_Integers`, `SYM_Rationals`, `SYM_Method`, `SYM_Modulus`,
`SYM_Automatic`.

---

## Router: `builtin_solve` (`src/solve.c`)

Signature: `Expr* builtin_solve(Expr* res)`. Return `NULL` to leave the
call unevaluated.

Attributes: `ATTR_PROTECTED | ATTR_HOLDALL`.  (HoldAll matches Mathematica
and lets the user write `Solve[x^2 == 4, x]` even when `x` already has a
value bound globally.)

```
1. argc = res->data.function.arg_count;  argc >= 2 required.
2. Walk trailing args left-to-right; an arg is an option iff it is
   Rule[sym, _] or RuleDelayed[sym, _] with sym one of the recognized
   option names (Cubics, Quartics, GeneratedParameters, VerifySolutions,
   Assumptions, InverseFunctions, Method, Modulus). Stop at the first
   non-option arg.
3. positional[0] = expr;  positional[1] = vars;  positional[2] = dom (optional).
4. Evaluate `expr` and `vars` once via evaluate() to expand assignments
   that the HoldAll deferred. (Mirrors how integrate.c handles its args.)
5. Normalize `vars`:
   - EXPR_SYMBOL → list_of_one.
   - List[v1, …]  → list_of_n. Each must be EXPR_SYMBOL or return NULL.
6. Normalize `expr` (single conjunction for now):
   - Equal[lhs, rhs] → one polynomial equation.
   - List[e1, …]   → require every e_i to be Equal[]; treat as system.
   - And[e1, …]    → same as List form.
   - Anything else → NULL.
7. Dispatch:
   - 1 variable + 1 Equal[] → forward to
     solvepoly_dispatch(equation, var, dom, &opts).
   - Otherwise → NULL (future: linear systems, inequalities, transcendental).
8. On success: return List[ List[Rule[v, sol1], …], List[Rule[v, sol2], …], … ].
9. On NULL return, do NOT free res (evaluator owns it).
   On non-NULL, expr_free(res) before returning.
```

The qualified specialist `Solve\`SolvePolynomialEquality` accepts the same
positional shape minus the routing logic — it requires its caller to have
already classified the input. It's registered via `symtab_add_builtin` on
the full context-qualified name (`context.c` handles parsing & display).

### Option parsing (mirrors `src/integrate.c:205-235`)

```c
typedef struct {
    bool cubics_radical;     // default false  (-> emit Root[] for cubics)
    bool quartics_radical;   // default false
    bool verify_solutions;   // default false  (Automatic ≈ false here)
    Expr* dom;               // default Complexes; borrowed pointer
    Expr* generated_param;   // default = symbol C (for future parametric)
} SolveOpts;

static bool parse_option(Expr* opt, SolveOpts* o);   // returns false on unknown opt
```

Unknown options trip a one-shot stderr warning (the `last_warned_hash`
pattern at `integrate.c:254-262`) and the whole call returns NULL.

---

## Specialist: `builtin_solve_polynomial_equality` (`src/solvepoly.c`)

Input: `Solve\`SolvePolynomialEquality[Equal[lhs, rhs], var, dom, opts…]`.

Algorithm:

```
1. Move to one side:
     poly_raw = Plus[lhs, Times[-1, rhs]]
     poly     = evaluate(poly_raw)
   (Plus already normalizes a + (-b) ; no special "subtract" helper exists,
    per src/parse.c:958-965.)

2. Check is_polynomial(poly, &var, 1) from src/poly.h:36.
   If not a polynomial in `var`, return NULL (Solve stays unevaluated).

3. degree d = get_degree_poly(poly, var)  (src/poly.h:38).

4. Fast paths (no factoring needed):
     d == 0 with poly != 0  → {}                  (no solutions)
     d == 0 with poly == 0  → {{}}                (full-dim / tautology)
     d == 1                  → solve_linear()
     d == 2                  → solve_quadratic()
     Two-term test (only coefficients at degree d and 0 nonzero)
                             → solve_binomial(d)
     n-quadratic test (only coefficients at degrees 2n, n, 0 nonzero,
       with n ≥ 1, and the degrees form a consistent stride)
                             → solve_nquadratic(a, b, c, n)
       This catches biquadratics (n=2) and higher n-quadratics. Output
       has 2n radical solutions; the Cubics/Quartics options have no
       effect on this path.

5. Otherwise, factor:
     collected = internal_collect({poly, var}, 2);          // poly.h via internal.h:38
     sqfree    = internal_factorsquarefree({collected}, 1); // internal.h:52
   Walk `sqfree` as Times[..., Power[f_i, m_i], ...] (m_i = 1 if bare).
   For each (f_i, m_i):
     factored = internal_factor({f_i}, 1);                  // internal.h:53
     Walk `factored` similarly. For each irreducible (g_j, n_j):
       dg = get_degree_poly(g_j, var)
       mult = m_i * n_j
       per-factor solutions (FIRST MATCH WINS):
         dg == 1: solve_linear(...) replicated `mult` times
         dg == 2: solve_quadratic(...) — each root replicated `mult` times
         g_j is binomial in var: solve_binomial(...) replicated `mult` times
         g_j is n-quadratic in var (degrees 2n, n, 0; n ≥ 2):
             solve_nquadratic(...) replicated `mult` times
             (checked BEFORE the cubic/quartic branches so it captures
             biquadratics even though dg == 4.)
         dg == 3: opts.cubics_radical
                     ? solve_cubic_radical(...)
                     : root_objects(g_j, var, 1..3)
         dg == 4: opts.quartics_radical
                     ? solve_quartic_radical(...)
                     : root_objects(g_j, var, 1..4)
         dg >= 5: root_objects(g_j, var, 1..dg)

6. Build the outer List of singleton-rule Lists:
     return List[ List[Rule[var, sol_k]] for each sol_k ]

7. Memory: every intermediate (poly, collected, sqfree, factored, opts.dom
   borrows) is freed exactly once. The accumulator that holds per-factor
   Expr** arrays is freed inclusive of unconsumed entries on early return.
```

### Sub-routines (all `static` in `src/solvepoly.c`)

```c
// All return freshly allocated Expr* (owned by caller). Inputs are
// deep-copied where needed; caller owns originals.

static Expr* solve_linear(Expr* a, Expr* b);
//   a x + b == 0  →  -b / a

static Expr** solve_quadratic(Expr* a, Expr* b, Expr* c,
                              bool reals_only, size_t* out_n);
//   Returns 0/1/2 solutions. Discriminant D = b^2 - 4 a c built symbolically.
//   reals_only:
//     - If Sign[D] proves D < 0 → out_n = 0.
//     - If Sign[D] proves D == 0 → out_n = 1 (single root, NOT replicated;
//       the caller multiplies by multiplicity).
//     - If Sign[D] proves D > 0 → out_n = 2 (real roots).
//     - If Sign[D] is unresolved → out_n = 2, symbolic Sqrt[D] form;
//       caller decides. (For dom = Reals with unresolved sign, conservative
//       fallback: emit both, leave the user with symbolic Re-test.)

static Expr** solve_binomial(Expr* a, Expr* b, int64_t n,
                             bool reals_only, size_t* out_n);
//   a x^n + b == 0  →  r = (-b/a)^(1/n)
//   Complexes: { r * Exp[2 Pi I k / n] : k = 0..n-1 }
//                (E^(2 Pi I k / n) collapses for k that is a multiple
//                 of n via existing trig/Exp rules — produced as
//                 literals via expr_new_function("Exp"...).)
//   Reals: discriminant-aware:
//     n odd  → 1 real root: { r }  (sign of -b/a doesn't matter; Mathilda's
//                                    Sign/Power handles -ve real bases for
//                                    odd integer roots).
//     n even, Sign[-b/a] > 0 → 2 roots: { ±r }
//     n even, Sign[-b/a] < 0 → 0 roots
//     unresolved sign       → fall through to Complexes (conservative).

static Expr** solve_nquadratic(Expr* a, Expr* b, Expr* c, int64_t n,
                               bool reals_only, size_t* out_n);
//   a x^(2n) + b x^n + c == 0
//   Internally:  let u = x^n.  Solve  a u^2 + b u + c == 0  via
//   solve_quadratic(...) to get {u_1, u_2}.  Then for each u_i, build
//   the binomial  x^n − u_i == 0  and call solve_binomial(1, -u_i, n, ...)
//   to recover n roots.  Concatenated: 2n roots in Complexes,
//   fewer in Reals (per the binomial Reals branch).
//   When reals_only is true and the inner quadratic in u has 0 real
//   roots, returns 0 solutions immediately.  When one or both u_i are
//   real, the binomial step's even/odd-n logic decides per u_i.

static Expr** solve_cubic_radical(Expr* a, Expr* b, Expr* c, Expr* d,
                                  bool reals_only, size_t* out_n);
//   Cardano. Reals branch uses sign of cubic discriminant
//     Δ = 18 a b c d − 4 b³ d + b² c² − 4 a c³ − 27 a² d²
//   Δ > 0: 3 distinct real roots (casus irreducibilis — uses trig form).
//   Δ == 0: multiple root, ≤ 2 distinct real.
//   Δ < 0: 1 real, 2 complex conjugate (return 1 if reals_only).

static Expr** solve_quartic_radical(Expr* a, Expr* b, Expr* c, Expr* d, Expr* e,
                                    bool reals_only, size_t* out_n);
//   Ferrari (depressed quartic + resolvent cubic). Reals filter via
//   resolvent cubic + per-pair discriminants.

static Expr* root_object(Expr* poly_in_var, Expr* var, int k);
//   Build Root[Function[poly_in_slot1], k] following the canonical
//   Slot[1] form at src/root.c:231-262. Reuses substitute_bvar_with_slot
//   if exposed; otherwise replicate the 18-line helper locally.
```

### Helpers (likely also `static`)

- `static bool is_binomial(Expr* poly, Expr* var, int d, Expr** a_out, Expr** b_out, int64_t* n_out);`
  — quick check: only the leading coefficient and the constant term are nonzero.
- `static bool is_nquadratic(Expr* poly, Expr* var, int d, Expr** a_out, Expr** b_out, Expr** c_out, int64_t* n_out);`
  — quick check: exactly three nonzero coefficients at degrees 2n, n, 0
    with n = d/2 ≥ 1. Falls through to the binomial path if b is zero
    (i.e., a pure binomial is handled by `is_binomial` instead, so n ≥ 2
    here is the meaningful case — though n == 1 just means "ordinary
    quadratic" and is harmless if it slips through).
- `static int try_sign(Expr* e);` — `1` / `0` / `-1` / `INT_MIN` for
  unresolved, via `evaluate(Sign[e])`.
- `static Expr* sqrt_expr(Expr* e);` — builds `Power[e, Rational[1,2]]`.

---

## Memory safety

Mathilda's standard contract (`SPEC.md` §4):

- `builtin_solve` and `builtin_solve_polynomial_equality` only call
  `expr_free(res)` along the success path. On `NULL` return, `res` is
  untouched.
- Every sub-routine returns a freshly allocated `Expr*` or `Expr**`
  arrays; the dispatcher accumulates and concatenates them, never freeing
  a borrowed pointer.
- Use `expr_copy()` whenever the same coefficient appears in two
  independently-owned output expressions.
- Test plan includes a valgrind pass over the full `test_solve.c` suite
  (per CLAUDE.md "Efficient and careful memory management").

---

## Registration

### `src/core.c`

Add (alphabetical with neighbors):
```c
#include "solve.h"   // near other includes around line ~30
…
solve_init();        // inside core_init(), near rat_init() / facpoly_init()
```

### `solve_init()` (in `src/solve.c`)

```c
void solve_init(void) {
    symtab_add_builtin("Solve", builtin_solve);
    symtab_get_def("Solve")->attributes |= ATTR_PROTECTED | ATTR_HOLDALL;
    symtab_set_docstring("Solve", "<see docstring below>");

    /* Option-name symbols are interned via sym_names_init; nothing to
       register here beyond the documentation. */
    symtab_set_docstring("Cubics", "<…>");
    symtab_set_docstring("Quartics", "<…>");
    symtab_set_docstring("GeneratedParameters", "<…>");
    symtab_set_docstring("VerifySolutions", "<…>");

    solvepoly_init();   /* declared in src/solvepoly.h */
}
```

### `solvepoly_init()` (in `src/solvepoly.c`)

```c
void solvepoly_init(void) {
    symtab_add_builtin("Solve`SolvePolynomialEquality",
                       builtin_solve_polynomial_equality);
    SymbolDef* def = symtab_get_def("Solve`SolvePolynomialEquality");
    def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Solve`SolvePolynomialEquality", "<…>");
}
```

### Docstring text

Follow the inline-in-init pattern (most modules do this; `info.c` is
reserved for Plus/Times/Power/Subtract that have no module init).

---

## Spec & changelog

- `docs/spec/builtins/arithmetic-and-algebra.md`: add `## Solve`,
  `## Cubics`, `## Quartics`, `## VerifySolutions`,
  `## Solve\`SolvePolynomialEquality` sections in the existing format
  (signature, options table, examples).
- `docs/spec/changelog/2026-05.md`: append `## Solve` with a 1-paragraph
  summary citing scope (router + polynomial equality), supported degree
  ranges, and the dom semantics.
- `Mathilda_spec.md`: only touch if adding a new top-level category
  (per CLAUDE.md guidance). "Polynomial Manipulation" already exists, so
  the per-category file change above is sufficient.

---

## Tests (`tests/test_solve.c`)

Use the `run_test(input, fullform_expected)` helper from
`tests/test_facpoly.c:13-30`.  Register in `tests/CMakeLists.txt`:

```cmake
add_executable(solve_tests test_solve.c ${COMMON_SRC})
target_include_directories(solve_tests PRIVATE ../src)
```

Coverage matrix:

| Category | Case | Expected |
|----------|------|----------|
| Linear | `Solve[2 x + 3 == 0, x]` | `{{x -> -3/2}}` |
| Linear LHS != 0 | `Solve[2 x + 3 == 1, x]` | `{{x -> -1}}` |
| Quadratic real | `Solve[x^2 - 5 x + 6 == 0, x]` | 2 rules: x→2, x→3 |
| Quadratic complex | `Solve[x^2 + 1 == 0, x]` | `{{x -> -I}, {x -> I}}` |
| Quad Reals filter | `Solve[x^2 + 1 == 0, x, Reals]` | `{}` |
| Multiplicity | `Solve[(x-1)^2 == 0, x]` | `{{x -> 1}, {x -> 1}}` |
| Pre-factored | `Solve[(x-1)(x-2)(x-3) == 0, x]` | 3 rules |
| Binomial Cx | `Solve[x^4 - 1 == 0, x]` | 4 complex roots |
| Binomial Reals | `Solve[x^4 - 1 == 0, x, Reals]` | 2 real roots: 1, -1 |
| Cubic default | `Solve[x^3 + x + 1 == 0, x]` | 3 Root[] objects |
| Cubic radical | `Solve[x^3 + x + 1 == 0, x, Cubics -> True]` | 3 radical exprs |
| Biquadratic | `Solve[x^4 - 5 x^2 + 4 == 0, x]` | 4 radical rules: ±1, ±2 (no Quartics opt needed) |
| Biquadratic Cx | `Solve[x^4 + 3 x^2 + 2 == 0, x]` | 4 complex radical rules: ±I, ±I Sqrt[2] |
| Biquadratic Reals | `Solve[x^4 + 3 x^2 + 2 == 0, x, Reals]` | `{}` |
| n-quadratic n=3 | `Solve[x^6 - 9 x^3 + 8 == 0, x]` | 6 radical rules (cube roots of 1 and 8) |
| n-quad does not need Quartics | `Solve[x^4 - 5 x^2 + 4 == 0, x, Quartics -> False]` | same 4 radical rules as above |
| Quartic non-biquad | `Solve[x^4 + x + 1 == 0, x]` | 4 Root[] objects (no n-quadratic shape) |
| Quartic radical | `Solve[x^4 + x + 1 == 0, x, Quartics -> True]` | 4 Ferrari radical rules |
| Quintic | `Solve[x^5 - x - 1 == 0, x]` | 5 Root[] objects |
| Mixed | `Solve[(x-1)(x^5 - x - 1) == 0, x]` | 1 rule + 5 Root[] |
| Trivial | `Solve[1 == 0, x]` | `{}` |
| Tautology | `Solve[0 == 0, x]` | `{{}}` |
| Non-poly | `Solve[Sin[x] == 0, x]` | unevaluated (Solve[…] returns) |
| Free of var | `Solve[y == 0, x]` | `{}` (constant non-zero in x) |
| Held var | `x = 5; Solve[x^2 == 4, x]` | `{{x -> -2}, {x -> 2}}` (HoldAll works) |
| Negative test | `Solve[x + 1 == 0, {x, y}]` | NULL (system-of-1 in 2 vars: router declines for now) |

All tests must pass under `make && ./solve_tests`, and the binary must
report zero leaks under `valgrind --leak-check=full --error-exitcode=1
./solve_tests`.

---

## Reused functions (no reinvention)

| Need | Function | File |
|------|----------|------|
| Factor input | `internal_factor(args, n)` | `src/internal.h:53` |
| Square-free split | `internal_factorsquarefree` | `src/internal.h:52` |
| Group by var | `internal_collect` | `src/internal.h:38` |
| Polynomial test | `is_polynomial` | `src/poly.h:36` |
| Degree | `get_degree_poly` | `src/poly.h:38` |
| Coefficient | `get_coeff` | `src/poly.h:39` |
| Free-of-var | `contains_any_symbol_from` | `src/poly.h:45` |
| Force-eval intermediate | `evaluate` | `src/eval.h` |
| Held `Root[Function[…], k]` | (replicate `substitute_bvar_with_slot`) | `src/root.c:212-228` |
| Option dispatch idiom | parse_method_option pattern | `src/integrate.c:215-235` |
| Build expressions | `expr_new_function`, `expr_new_symbol`, `expr_new_integer`, `expr_copy`, `expr_free` | `src/expr.h` |

---

## Files modified / added

| File | Change |
|------|--------|
| `src/solve.h` | **NEW** — declarations for router + `solve_init` |
| `src/solve.c` | **NEW** — router + option parser |
| `src/solvepoly.h` | **NEW** — declarations for specialist |
| `src/solvepoly.c` | **NEW** — polynomial-equality solver + per-degree subroutines |
| `src/sym_names.h` | add SYM_* declarations listed above |
| `src/sym_names.c` | add SYM_* NULL defs + intern_symbol() lines |
| `src/core.c` | `#include "solve.h"` and call `solve_init()` |
| `tests/test_solve.c` | **NEW** |
| `tests/CMakeLists.txt` | add `solve_tests` executable |
| `docs/spec/builtins/arithmetic-and-algebra.md` | add 5 sections |
| `docs/spec/changelog/2026-05.md` | append `## Solve` entry |
| `SOLVE_PLAN.md` | **NEW** — checked into repo root as the canonical design doc (per user request) |

---

## Verification (end-to-end)

1. `make -j$(nproc)` from repo root — must compile cleanly under
   `gcc -std=c99 -Wall -Wextra`; no POSIX-only types (`ssize_t`,
   `strdup`, `M_PI`, etc.) per CLAUDE.md guardrails.
2. `cd tests && mkdir -p build && cd build && cmake .. && make -j` — builds
   the new `solve_tests` target alongside existing ones.
3. `./solve_tests` — every row in the test matrix passes.
4. `valgrind --leak-check=full --error-exitcode=1 ./solve_tests` — zero
   leaks, zero errors.
5. `for t in *_tests; do "$t" || echo "REGRESSION in $t"; done` — no
   pre-existing test regresses (Solve doesn't touch their code paths, but
   verify nothing in sym_names.c, core.c, or info.c changes broke them).
6. REPL smoke check: `./Mathilda` then manually evaluate `Solve[x^2 - 5 x + 6 == 0, x]`,
   `Solve[x^4 - 1 == 0, x, Reals]`, `Solve[x^5 - x - 1 == 0, x]`, `?Solve`.
   Confirm output and docstring render.
