# `Product` — Development Plan (build sheet)

> Executable build sheet for the symbolic `Product[]` subsystem, structured so
> each stage can be picked up independently. Companion to the algorithm survey in
> [`PRODUCT_PLAN.md`](PRODUCT_PLAN.md) (read that first for the math and
> references) and modelled on [`SUM_DEVEL_PLAN.md`](SUM_DEVEL_PLAN.md).
>
> **Distinct from `NProduct`** (numerical, shipped — `NPRODUCT_PLAN.md`).
> `Product` returns an *exact closed form* (`n!`, `Pochhammer[a,n]`, `Gamma`
> ratios, `QPochhammer`, …). The naïve `Exp[NSum[Log f]]` reduction that powers
> `NProduct` does **not** work symbolically (`Log f` is not hypergeometric even
> when `f` is) — products are handled in their own multiplicative terms.
>
> **Two headline requirements (from the brief):**
> 1. **Every algorithm is reachable via `Method`** — `Method -> "Telescoping" |
>    "Rational" | "Geometric" | "QProduct" | …`.
> 2. **`Method -> Automatic` (the default) is a polyalgorithm** that selects the
>    best algorithm for the input via a cascade, exactly as `Sum`/`Integrate` do.
> 3. **Extensive unit *and* stress tests** (§9).

---

## 1. Goal & scope

Add a builtin `Product[f, {i, imin, imax (, di)}, opts]` — the multiplicative
analogue of the existing `Sum` (`src/sum/`). It must:

- evaluate finite explicit products (`Product[k,{k,1,5}] → 120`);
- give closed forms for symbolic-bound products of hypergeometric terms
  (`Product[k,{k,1,n}] → n!`, `Product[k+a,{k,1,n}] → Pochhammer[a+1,n]`);
- give indefinite products / anti-quotients (`Product[f, i]`);
- evaluate convergent infinite products (`Product[1-1/k^2,{k,2,∞}] → 1/2`,
  Wallis → `π/2`) and flag divergent ones;
- support multidimensional `Product[f, {i,…}, {j,…}, …]` with inner bounds that
  may depend on an outer index;
- expose each algorithm through `Method` with an `Automatic` polyalgorithm
  default;
- carry attributes `HoldAll | Protected` (the index is `Block`/own-value
  localised exactly as `Sum`/`Table`/`Do`);
- be valgrind-clean and add **zero** new leaks.

`Product` mirrors `Sum`'s architecture so closely that **`src/sum/` is the
reference implementation** — read each `sum_*.c` alongside its `product_*.c`
target.

---

## 2. Architecture — a multiplicative mirror of `src/sum/`

Create `src/product/`, one dispatcher + one file per algorithm, each algorithm
also exposed as a context-qualified builtin (``Product`Rational`` etc.), exactly
as `Sum` exposes ``Sum`Gosper`` and `Integrate` exposes
``Integrate`BronsteinRational``.

```
src/product/
  product.h            mirror of sum.h
  product.c            dispatcher (Stage 0): surface, finite expansion,
                       multi-iterator nesting, Method cascade
  product_internal.h   shared-helper API (mirror of sum_internal.h)
  product_util.c       prod_eval / prod_subst / prod_factor / prod_free_of /
                       prod_int / prod_div / product_stage_args
  product_telescoping.c   Product`Telescoping   (Stage 1)
  product_rational.c      Product`Rational      (Stage 2 — the workhorse)
  product_geometric.c     Product`Geometric     (Stage 3)
  product_infinite.c      Product`Infinite      (Stage 4; or fold into try_*)
  product_qproduct.c      Product`QProduct      (Stage 5 — needs QPochhammer)
```

**Adding a stage stays additive** (the `Sum` invariant): a new
`src/product/product_<name>.c`, one `try_*` line in `dispatch_def`/
`dispatch_indef` (`product.c`), one `product_<name>_init()` call in
`product_init()`. The cascade tolerates an absent stage (an unresolved
context-builtin call counts as "fall through").

**Build wiring** (the top-level `makefile` auto-discovers `src/**/*.c`, so
`./Mathilda` picks the files up automatically; only the *tests* need explicit
listing — see §8):
- register from `core_init()` in `src/core.c` by calling `product_init()` right
  after the existing `sum_init();` line.

---

## 3. Hard-won constraints (read before coding any stage)

Inherited from the `Sum` build (`SUM_DEVEL_PLAN.md` §2) and project memory; the
product-specific ones are flagged **[P]**.

1. **`Together`/`Factor` infinite-loop on a symbolic exponential** — any
   `Power[sym, var]` (e.g. `a^i`, `q^k`). `Cancel` and plain `evaluate` are safe.
   **Rule:** never run `Together`/`Factor` on an expression carrying `base^var`;
   normalise only the rational *coefficient* and multiply the exponential back
   in. **[P] This is central here**, not a corner case: Stage 3 (geometric) and
   Stage 5 (`q`-products) live on `base^k` factors. See `sum_geometric.c` and
   memory `project_together_factor_hang_exponential`.
2. **`get_degree_poly` needs `Expand` first** (it does not expand `(i+3)^5`).
3. **`Solve` cannot solve symbolic linear *systems*** — use triangular
   back-substitution or `SolveAlways`; for `Product` the root extraction in
   Stage 2 uses `Factor` + factor-walking, not `Solve`, sidestepping this.
4. **Finite numeric ranges never reach the cascade** — Stage 0 expands them
   directly, so the closed-form stages only ever see symbolic bounds, `Infinity`,
   or the indefinite form. They may assume non-numeric limits.
5. **Output need not be maximally combined.** The test oracle is mathematical
   (§9), not string equality.
6. **[P] Empty / reversed range is `1`, not `0`** (Karr convention): the
   multiplicative identity. `fold_times` of zero terms returns `1`; the
   telescoping shortcut must guard `imin <= imax+1` and fall back to expansion
   (→ `1`) on an empty range, exactly as `Sum` guards `min<=max` to fold to `0`.
7. **[P] A single zero factor forces the whole product to `0`** regardless of the
   closed form — relevant for definite products whose range straddles a root, and
   for `Gamma`-pole bookkeeping. Finite expansion handles this for free; the
   closed-form stages assume the range avoids exact integer roots of the
   denominator (true for symbolic bounds).
8. **[P] `evaluate(X)` does not free `X`** — bind every
   `evaluate(expr_new_function(...))` node to a temp and `expr_free` it. The
   evaluator owns `res`; **never `expr_free(res)`** in the builtin (memory:
   `feedback_builtin_res_ownership`, `project_simplify_zero_leak`).
9. **[P] `Gamma`/`Pochhammer`/`FactorialPower` now evaluate** (added since the
   Sum plan; verified `src/special_functions/gamma.c`,
   `src/special_functions/pochhammer.c`, `src/numbertheory.c`). So closed forms
   built from them reduce at integer bounds — the finite-expansion oracle (§9)
   is viable. `QPochhammer`, `BarnesG`, `Hyperfactorial`, Glaisher are **missing**
   and gate Stages 5–6 (§7).

---

## 4. The `Method` polyalgorithm (core requirement)

`product.c` parses a `Method` option and drives either a single named stage or
the `Automatic` cascade — the exact shape of `sum.c`'s `SumMethod` /
`parse_method_option` / `dispatch_def` / `dispatch_indef`.

```c
typedef enum {
    PROD_METHOD_AUTOMATIC = 0,
    PROD_METHOD_TELESCOPING,
    PROD_METHOD_RATIONAL,      /* aka "Hypergeometric"            */
    PROD_METHOD_GEOMETRIC,     /* aka "PolynomialExponential"     */
    PROD_METHOD_QPRODUCT,
    PROD_METHOD_INVALID
} ProdMethod;
```

`Method -> "Name"` string map (mirror `parse_method_option`; an unrecognised
string → `PROD_METHOD_INVALID` → `Product` stays unevaluated, matching `Sum`):

| `Method` value | Stage | Algorithm |
|---|---|---|
| `Automatic` (default) | — | polyalgorithm cascade (below) |
| `"Telescoping"` | 1 | rational anti-quotient `g(k+1)/g(k)` |
| `"Rational"` / `"Hypergeometric"` | 2 | rational → Pochhammer/Gamma |
| `"Geometric"` / `"PolynomialExponential"` | 3 | `base^k` factors via `Sum` exponent |
| `"QProduct"` | 5 | `q`-rational → QPochhammer |

**The `Automatic` polyalgorithm cascade** (cheapest/most-specific first, so the
*nicest* closed form wins — pure rational before Gamma before `q`):

```c
static Expr* dispatch_def(ProdMethod m, Expr* f, Expr* var, Expr* lo, Expr* hi) {
    switch (m) {
    case PROD_METHOD_AUTOMATIC: {
        Expr* r = try_def("Product`Telescoping", f, var, lo, hi);   /* clean rational */
        if (!r) r = try_def("Product`Rational",  f, var, lo, hi);   /* Pochhammer/Gamma */
        if (!r) r = try_def("Product`Geometric",  f, var, lo, hi);  /* base^k */
        if (!r) r = try_def("Product`QProduct",   f, var, lo, hi);  /* q-rational */
        return r;
    }
    case PROD_METHOD_TELESCOPING: return try_def("Product`Telescoping", f, var, lo, hi);
    case PROD_METHOD_RATIONAL:    return try_def("Product`Rational",    f, var, lo, hi);
    case PROD_METHOD_GEOMETRIC:   return try_def("Product`Geometric",   f, var, lo, hi);
    case PROD_METHOD_QPRODUCT:    return try_def("Product`QProduct",    f, var, lo, hi);
    default: return NULL;
    }
}
```

`dispatch_indef` is the same over `head[f,i]`. **Cascade ordering rationale:**
Telescoping first because it yields a Gamma-free rational answer; `Rational`
second because it subsumes telescoping but introduces `Pochhammer`/`Gamma`;
`Geometric` and `QProduct` handle the `base^k` shapes the rational stages
reject. (`Product`Infinite` of Stage 4 is *not* a cascade entry — it is a
wrapper that, when `imax == Infinity`, runs the cascade for the finite form then
takes a limit / consults the identity table; see §6.4.)

---

## 5. Shared helpers (`product_internal.h` / `product_util.c`)

Direct multiplicative analogue of `sum_internal.h` / `sum_util.c`. Identical
except `prod_div` replaces `sum_sub`, and `prod_factor` carries the symbolic-power
guard.

```c
Expr* prod_eval(const char* head, Expr** args, size_t n);      /* = sum_eval   */
Expr* prod_subst(Expr* e, Expr* var, Expr* val);               /* = sum_subst  */
Expr* prod_factor(Expr* e);    /* Factor with fallback; caller must ensure e   */
                               /* carries no base^var (constraint 3.1)         */
bool  prod_free_of(Expr* e, Expr* var);                        /* = sum_free_of */
Expr* prod_int(int64_t v);
Expr* prod_div(Expr* a, Expr* b);   /* a / b  =  Times[a, Power[b,-1]], eval'd  */
bool  product_stage_args(Expr* res, Expr** f, Expr** var,
                         Expr** imin, Expr** imax, bool* definite); /* = sum_stage_args */
```

`product_stage_args` is identical to `sum_stage_args`: `head[f,i]` (argc 2,
indefinite) or `head[f,i,imin,imax]` (argc 4, definite); out-pointers alias
`res`'s args.

Two product-specific helpers used by Stages 1–2 (put them in `product_util.c`):

```c
/* True if e contains Power[base, expr-involving var] — gates prod_factor/Together. */
bool prod_has_symbolic_power(Expr* e, Expr* var);

/* Factor a polynomial in var into  c, [{root_i, mult_i}]  over Q (rational roots
 * exact; irreducible-quadratic+ roots reported so the caller can choose to bail
 * or emit RootSum-style output). Built on the existing Factor / factor-walk +
 * get_degree_poly/get_coeff. Returns false if e is not polynomial in var. */
bool prod_linear_factors(Expr* e, Expr* var,
                         Expr** lead_out, /* c */
                         Expr*** roots_out, int** mults_out, size_t* n_out,
                         bool* all_linear_out);
```

---

## 6. Stages

### Stage 0 — dispatcher `product.c`

A line-for-line port of `sum.c` with four changes:

1. **`fold_times` replaces `fold_plus`** — empty product → `expr_new_integer(1)`;
   otherwise `Times[terms…]` evaluated. (`Sum`'s `fold_plus` returns `0` /
   `Plus[…]`.)
2. **Finite expansion** (`expand_list`, `expand_range`) — identical iteration,
   folding with `fold_times`. Reuse `iter.h` (`iter_spec_parse`,
   `iter_spec_resolve_numeric`, `iter_spec_shadow`/`_restore`) verbatim.
3. **Telescoping shortcut guard** — for a unit-step integer-bounded **non-empty**
   range (`!is_real && di==1 && min<=max`), call `dispatch_def` first (the
   closed form is cheaper than per-term expansion). Empty range must fall through
   to expansion → `1` (constraint 3.6). Same guard structure as `sum_one_spec`.
4. **Multi-iterator nesting** — `Product[f, s1,…,sk, opts]` →
   `Product[Product[f, sk, opts], s1,…,s_{k-1}, opts]`, re-evaluated. Identical
   to `sum.c`'s rewrite; the inner product is evaluated under each outer binding,
   so inner bounds may depend on outer indices.

`builtin_product` is `sum.c`'s `builtin_sum` with `SYM_Sum`→`SYM_Product`,
`fold_plus`→`fold_times`, and the `SumMethod` type swapped for `ProdMethod`.

```c
void product_init(void) {
    symtab_add_builtin("Product", builtin_product);
    symtab_get_def("Product")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Product", /* terse — examples live in docs/spec */
        "Product[f, {i, imax}] gives the product of f for i from 1 to imax. "
        "Product[f, {i, imin, imax}], Product[f, {i, imin, imax, di}] and "
        "Product[f, {i, {i1, i2, ...}}] use the standard iterator forms; multiple "
        "iterators give nested products. Product[f, i] gives the indefinite "
        "product. Symbolic and infinite products are evaluated in closed form via "
        "Method -> \"Telescoping\" | \"Rational\" | \"Geometric\" | \"QProduct\".");

    void product_telescoping_init(void); product_telescoping_init();
    void product_rational_init(void);    product_rational_init();
    void product_geometric_init(void);   product_geometric_init();
    /* product_qproduct_init(); product_infinite_init();  — later stages */
}
```

**Done when:** `Product[k,{k,1,5}]→120`, `Product[k,{k,1,0}]→1` (empty),
`Product[2,{k,1,n}]` stays held (no stage yet) but `Product[2,{k,1,3}]→8`,
multi-iterator `Product[i j,{i,1,2},{j,1,2}]→4`, attributes
`HoldAll|Protected`, and `Method->"Bogus"` stays unevaluated.

### Stage 1 — `Product`Telescoping` (`product_telescoping.c`)

Catches the products whose anti-quotient is itself **rational** (no
`Gamma`/`Pochhammer` needed), giving the cleanest answers:
`Product[1+1/k,{k,1,n}]→n+1`, `Product[k/(k+1),…]→1/(n+1)`,
`Product[1-1/k^2,{k,2,n}]→(n+1)/(2n)`.

**Algorithm.** For rational `f`, find rational `g` with `f(k) = g(k+1)/g(k)`;
then `∏_{k=imin}^{imax} f = g(imax+1)/g(imin)` (indefinite: return `g(i)` up to
the free constant, mirroring `Sum`'s indefinite convention).

Construction via factor-chains (the dispersion idea in miniature):
1. Bail if `prod_has_symbolic_power(f, var)` (constraint 3.1).
2. `Together[f]`, split numerator/denominator, `prod_linear_factors` each →
   a root→exponent map `e(r)` (numerator roots positive, denominator negative).
   Bail unless all factors are linear over ℚ (`all_linear`).
3. Group roots into **integer-spaced chains** (same fractional part). `f`
   telescopes iff, within every chain, `e` is the discrete difference of some
   integer exponent function `d`: i.e. `e(r) = d(r-1) - d(r)` solved by the
   prefix sum `d(r) = -∑_{s ≤ r in chain} e(s)`, and the chain's exponents sum to
   `0` (so the ladder closes). Each generator contributes `(k - r)^{d(r)}` to
   `g`.
4. If every chain closes **and** the leading constant is `1`, emit
   `g = ∏ (k - r)^{d(r)}`; return `prod_div(prod_subst(g, var, imax+1),
   prod_subst(g, var, imin))`, run through `Cancel` (not `Factor` — cheaper and
   no hang risk).
5. Otherwise return `NULL` (fall through to Stage 2, which always succeeds for
   rational `f` but introduces Gammas).

**Implementation note.** Stage 1 and Stage 2 share `prod_linear_factors`; Stage 1
is precisely "the all-cancel case of Stage 2." An acceptable simpler first cut:
implement Stage 2, then implement Stage 1 as "run Stage 2, and if the result is
`Gamma`/`Pochhammer`-free after `Simplify`, keep it; else fall through" — but the
standalone chain construction above is preferred (no dependence on `Gamma`
simplification, and it is what makes `Method->"Telescoping"` meaningfully
distinct).

**Done when:** the three targets above plus `Product[(k+1)/(k+2),{k,1,n}]` give
Gamma-free rational closed forms, oracle-verified at `n→7`.

### Stage 2 — `Product`Rational` (`product_rational.c`) — the workhorse

Closed form for **any** rational `f` of `k` (the Maple/SymPy `_eval_product`
core). Output in `Pochhammer` (primary) / `Gamma` (fallback) / `Factorial` /
`FactorialPower`.

**Atomic identity** (definite, unit step):
```
∏_{k=imin}^{imax} (k + a) = Pochhammer[imin + a, imax - imin + 1]
                          = Gamma[imax + a + 1] / Gamma[imin + a].
```

**Pipeline** for `∏_{k=imin}^{imax} f`:
1. Bail if `prod_has_symbolic_power(f, var)` (→ Stage 3 territory).
2. `Together[f]` → `num/den`; `prod_linear_factors` num and den.
   `all_linear == false` (irreducible quadratic+) → **gate**: handle real-linear
   first; for an irreducible factor emit the `Gamma`-over-roots form only behind
   a follow-up flag (like `Sum`Rational`'s `RootSum` gate). First cut: bail.
3. Leading constant `c` (= `lead_num/lead_den`): contributes `c^(imax-imin+1)`.
4. Each numerator root `-a` of multiplicity `m`:
   `Pochhammer[imin+a, imax-imin+1]^m`; each denominator root divides.
5. **Minimal-output pass (optional, deferrable):** match a numerator root and a
   denominator root differing by an integer → telescope that pair to a rational
   factor instead of two Gammas (the "minimal number of factorials" refinement).
   Deferring this just yields a less-reduced but correct `Gamma` ratio that
   `Simplify` can often finish.
6. Indefinite form `Product[f, i]`: `∏_{k=1}^{i}(k+a) = Pochhammer[a+1, i]`
   (mirror `Sum`'s indefinite lower-bound convention; document it).
7. Assemble `Times[…]`, evaluate, return. **Do not `Factor`** the assembled
   Gamma/Pochhammer product (it isn't polynomial); a light `Simplify` or
   `FunctionExpand`-free `evaluate` suffices.

`f` is rational in `k` ⇒ no `k` in any exponent ⇒ `Factor`/root extraction are
hang-safe (constraint 3.1 does not bite here).

**Targets:** `Product[k,{k,1,n}]→n!`, `Product[k^2,…]→(n!)^2`,
`Product[k+a,…]→Pochhammer[1+a,n]`, `Product[(2k-1)/(2k),{k,1,n}]`
(→ `Binomial[2n,n]/4^n` after Simplify), `Product[x+i,{i,0,n-1}]→Pochhammer[x,n]`,
`Product[x-i,{i,0,n-1}]→FactorialPower[x,n]`.

**Done when:** the targets oracle-verify at `n→7` (Pochhammer/Gamma evaluate to
exact rationals/integers), and prior stages still pass.

### Stage 3 — `Product`Geometric` (`product_geometric.c`)

Handles factors `base^k` with `base` free of `k`, and `base^(linear in k)`:
`∏_{k} c^k = c^(∑ k)` — route the **exponent** through the shipped `Sum`.

**Pipeline:**
1. Split `f` into the part free of `var` and the rest (`prod_free_of`); a
   constant factor `c` → `c^(imax-imin+1)`.
2. For each factor `base^(p(k))` with `base` free of `var` and `p` polynomial in
   `var`: emit `base^(Sum[p(k), {k, imin, imax}])` — `Sum` returns the polynomial
   series closed form (Stages 1–2 of `Sum`). **Keep `base^…` factored; never
   `Together`/`Factor` it** (constraint 3.1).
3. A residual rational cofactor → hand to `Product`Rational` (Stage 2) and
   multiply.

**Targets:** `Product[2^k,{k,1,n}]→2^(n(n+1)/2)`, `Product[a^k,{k,0,n}]`,
`Product[k 2^k,{k,1,n}]` (geometric × rational). **Done when** these
oracle-verify and `Together`/`Factor` are never reached with a `base^var` present
(add a regression that would hang if violated, run under a timeout).

### Stage 4 — `Product`Infinite` + convergence (`product_infinite.c`)

`imax == Infinity`. Two routes, tried in order:

1. **Identity table (`.m` rules).** A new internal file
   `src/internal/product_infinite.m` (loaded from `init.m`, mirroring `deriv.m` /
   the integral tables) with the recognised convergent infinite products as
   `Product` DownValues / a recognizer:
   - Wallis `∏ 4k²/(4k²-1) = π/2`;
   - sine/cosine: `∏(1 - z²/k²) = Sin[π z]/(π z)`, cosine analogue;
   - the family `∏_{k≥1}(k²+a²)/(k²+b²) = b Sinh[π a]/(a Sinh[π b])`, special
     cases `∏(1+1/k²)=Sinh[π]/π`, telescoping `∏_{k≥2}(1-1/k²)=1/2`.
   This is "C for performance, rules for mathematics" (SPEC §13) — table-driven
   recognition belongs in `.m`.
2. **Limit of the finite closed form.** Run the finite cascade (Stages 1–3) for
   symbolic `imax`, then `Limit[closedform, imax → Infinity]` (`src/calculus/`).
   Recognise residual `Gamma`/`Pochhammer` limits.

**Convergence gate (must run before asserting any infinite value).** From
`PRODUCT_PLAN.md` §6.3:
- no zero factor in range;
- **term → 1** (`Limit[f, k→∞] == 1`), else divergent → emit `Product::div`,
  stay unevaluated (matches `Sum`'s `Sum::div`);
- for `f = P/Q`: require `deg P == deg Q`, **equal leading coefficients**, **and
  equal next-to-leading coefficients** (root-sum P = root-sum Q) — the decisive
  test that separates convergent `a_k ~ c/k²` from divergent `a_k ~ c/k`. This is
  pure polynomial inspection (`get_degree_poly` + `get_coeff`), GMP-exact.

**Options** (wire into `product.c`'s option parser, extending the `Method`
parse): `VerifyConvergence -> True|False` (default `True`),
`GenerateConditions -> True|False`.

**Targets:** `Product[1-1/k^2,{k,2,∞}]→1/2`, `Product[1+1/k^2,{k,1,∞}]→Sinh[π]/π`,
Wallis → `π/2`, and `Product[1+1/k,{k,1,∞}]` → divergent (held, `Product::div`).

### Stage 5 — `Product`QProduct` (`product_qproduct.c`) — gated on `QPochhammer`

`q`-rational factors (rational in `q^k`): `∏_{k=0}^{n-1}(1 - a q^k) =
QPochhammer[a, q, n]`. Mirror of Stage 2 with `Pochhammer→QPochhammer`,
`Gamma→QGamma`. **Prerequisite:** a `QPochhammer` builtin (and `QGamma`,
`[n]_q!`) — none exist (§7). Keep `q^k` factored throughout (constraint 3.1).
**Targets:** `Product[1-q^k,{k,1,n}]→QPochhammer[q,q,n]`,
`Product[1-a q^k,{k,0,n-1}]→QPochhammer[a,q,n]`.

### Stage 6 — Barnes-G / Hyperfactorial / Glaisher results — gated

`Product[k^k,{k,1,n}]→Hyperfactorial[n]`,
`Product[Gamma[k],{k,1,n-1}]→BarnesG[n]`. **Prerequisites:** `BarnesG`,
`Hyperfactorial` builtins + the Glaisher–Kinkelin constant (§7). Best as `.m`
recognizer rules once the special functions exist.

### Stage 7 — GP-form unification / minimality (refinement)

Rebuild Stages 1–2 on an explicit Gosper–Petkovšek normal form of the term-ratio
(`PRODUCT_PLAN.md` §5.1), and/or the Abramov–Petkovšek minimal decomposition
(§5.2) for provably-minimal `Gamma` output. Reuses the additive Gosper engine
(`sum_gosper.c`) + `src/poly/` factorisation + dispersion. Lowest priority; the
§5–6 root-matching already gives near-minimal output for common inputs.

---

## 7. Prerequisite builtins (general, not in `src/product/`)

Each is a standalone builtin (`src/special_functions/…` or `src/numbertheory.c`),
registered from `core_init()`, with attributes, a terse `symtab_set_docstring`,
exact/symbolic rules, a `D` rule, and machine **and** MPFR `N` evaluation (cover
`EXPR_REAL` *and* `EXPR_MPFR`, add an `N[x,35]` test — memory
`feedback_numeric_builtins_cover_mpfr`).

| Builtin | Status (verified 2026-06-18) | Needed by |
|---|---|---|
| `Pochhammer`, `Gamma`, `FactorialPower`, `Factorial`, `Binomial`, `PolyGamma` | **exist & evaluate** | Stages 1–4 |
| `Sin`, `Sinh`, `Pi`, `Limit` | exist | Stage 4 |
| `QPochhammer[a,q,n]`, `QGamma[z,q]` | **missing** | Stage 5 |
| `BarnesG[z]`, `Hyperfactorial[n]`, `Glaisher` | **missing** | Stage 6 |

**`QPochhammer[a, q, n]` and `QPochhammer[a, q]`** (`= (a;q)_∞`, `|q|<1`).
`Listable | NumericFunction | Protected`. Exact finite product for non-negative
integer `n`; `QPochhammer[q,q,n] = (q;q)_n`; `N` via the defining product /
`mpfr`. Foundation for the whole `q`-stage (`QGamma`, `[n]_q!`, Gaussian
binomials are thin wrappers).

**`BarnesG[z]`** — `G(z+1) = Gamma[z] G(z)`, `G(1)=G(2)=1`, integer
`G(n+1) = ∏_{k=1}^{n-1} k!` (superfactorial, exact with GMP), non-integer via a
LogGamma/`ζ'(-1)` asymptotic. **`Hyperfactorial[n]` = `∏_{k=1}^{n} k^k`** (exact
GMP for integer `n`; relates to `BarnesG` and Glaisher `A = exp(1/12 - ζ'(-1)) ≈
1.2824271291`). **`Glaisher`** — protected constant (`Constant` attribute, MPFR
value), like `EulerGamma`.

---

## 8. Exact edit sites

**Create** `src/product/{product.h, product.c, product_internal.h,
product_util.c, product_telescoping.c, product_rational.c, product_geometric.c}`
(Stages 0–3); add `product_infinite.c`, `product_qproduct.c` with their stages.

**`src/sym_names.h` / `src/sym_names.c`** — intern `Product` (next to
`SYM_Sum`, `sym_names.c:1027`). The `Method` values are *strings*, not symbols,
so no per-method symbols. Later stages also intern `QProduct` (if exposed),
`BarnesG`, `Hyperfactorial`, `Glaisher`, `GenerateConditions`. (`VerifyConvergence`
already exists — used by `NProduct`/`NSum`.)

**`src/core.c`** — add `void product_init(void); product_init();` immediately
after the existing `sum_init();` call.

**`src/info.c`** — docstrings for `Product` and each ``Product`*`` stage (terse;
examples go in `docs/spec` — memory `feedback_no_examples_in_docstrings`).

**`tests/CMakeLists.txt`** — two sites (mirror the `Sum` block):
- in `COMMON_SRC`, after the `src/sum/…` lines (~261–266), add the six/seven
  `../src/product/product*.c` files. **A new `src/*.c` MUST be in `COMMON_SRC`
  or every `*_tests` binary fails to link `_product_init`** (memory
  `project_tests_common_src_list`).
- after the `sum_tests` block (~340–342), add:
  ```cmake
  add_executable(product_tests test_product.c ${COMMON_SRC})
  target_link_libraries(product_tests m)
  target_include_directories(product_tests PRIVATE ../src ../src/product)
  add_executable(product_stress_tests test_product_stress.c ${COMMON_SRC})
  target_link_libraries(product_stress_tests m)
  target_include_directories(product_stress_tests PRIVATE ../src ../src/product)
  ```
  (and `add_dependencies(product_tests build_ecm)` etc. if ECM is linked, as the
  other suites do ~1600).

**`docs/spec/builtins/calculus.md`** — a `## Product` section beside `## Sum` /
`## NProduct`.

**`docs/spec/changelog/2026-06-15.md`** — the Monday-of-ISO-week file for
2026-06-18 (Thu) is **2026-06-15**. Add the feature note (create the file with a
`# Changelog: week of …` heading if absent, and add a row to `Mathilda_spec.md`'s
changelog table).

---

## 9. Test plan — unit **and** stress (the brief's emphasis)

### 9.1 Unit tests — `tests/test_product.c`

Port `tests/test_sum.c`'s harness verbatim: `eval_str`, the aborting `check()`
(exact `InputForm` string), and the oracle. **The product oracle is
multiplicative** — add alongside `same()`:

```c
/* Mathematical equality of a closed form against direct finite expansion:
 * substitute the symbolic bound n -> several concrete integers and compare to
 * the Stage-0 expanded product. Robust to Pochhammer/Gamma vs rational forms,
 * since both sides evaluate to exact numbers. */
static void oracle(const char* closed_form_in_n, const char* product_in_n) {
    for (int n = 0; n <= 8; n++) {            /* include n=0: empty product = 1 */
        char buf[512];
        snprintf(buf, sizeof buf, "Simplify[(%s) - (%s) /. n -> %d]",
                 closed_form_in_n, product_in_n, n);
        char* s = eval_str(buf);
        if (strcmp(s, "0") != 0) { /* FAIL: report n, both sides */ }
        free(s);
    }
}
```

Use `oracle("n!", "Product[k,{k,1,n}]")` etc. Prefer the oracle over brittle
string `check()` for any Gamma/Pochhammer answer (output form is not canonical —
constraint 3.5); reserve `check()` for the canonical cases (`n!`, `n+1`,
`1/(1+n)`, `120`, `1`).

**Unit case checklist** (one block per stage):
- *Stage 0:* `Product[k,{k,1,5}]→120`; empty `Product[k,{k,1,0}]→1` and reversed
  `Product[k,{k,5,1}]→1`; `Product[3,{k,1,4}]→81`; list spec
  `Product[k,{k,{2,3,5}}]→30`; multi-iterator `Product[i+j,{i,1,2},{j,1,2}]`;
  dependent inner bound `Product[i,{i,1,n},{}]`… ; `HoldAll`/`Protected` via
  `MemberQ[Attributes[Product], …]`; `Method->"Bogus"` and 1-arg `Product[f]`
  stay unevaluated.
- *Stage 1:* `Product[1+1/k,{k,1,n}]→1+n`; `Product[k/(k+1),{k,1,n}]`;
  `Product[1-1/k^2,{k,2,n}]`; oracle each. `Method->"Telescoping"` forces this
  stage; verify it returns Gamma-free output.
- *Stage 2:* `Product[k,{k,1,n}]→n!`; `Product[k^2,…]`; `Product[k+a,…]`;
  `Product[(2k-1)/(2k),{k,1,n}]`; `Product[x+i,{i,0,n-1}]→Pochhammer[x,n]`;
  indefinite `Product[k,k]`; `Method->"Rational"` forcing.
- *Stage 3:* `Product[2^k,{k,1,n}]`; `Product[a^k,{k,0,n}]`;
  `Product[k 2^k,{k,1,n}]`; a **hang regression** (`Product[a^k,{k,1,n}]` must
  return promptly — guards constraint 3.1) run with a wall-clock guard.
- *Stage 4:* `Product[1-1/k^2,{k,2,∞}]→1/2`; `Product[1+1/k^2,{k,1,∞}]→Sinh[π]/π`;
  Wallis→`π/2`; divergent `Product[1+1/k,{k,1,∞}]` stays held;
  `VerifyConvergence->False` behaviour.
- *Methods:* for each algorithm, assert `Method->"Name"` yields the same value as
  `Automatic` on an in-class input, and falls through (unevaluated) on an
  out-of-class input — proving the polyalgorithm and the per-method routing.
- *Per new special function* (Stages 5–6): an `N[x,35]` MPFR test.

### 9.2 Stress tests — `tests/test_product_stress.c` (+ optional `.m` corpus)

Model on the recent **NSolve stress suite** (`git log`: `test(nsolve)`) and
`IntegrateRationalTests.m`. Two layers:

1. **Generated families, oracle-verified.** Programmatically build large batches
   and check each by finite expansion at several bounds (no hand-written
   expected strings — the oracle is ground truth):
   - rational products `∏ (k+a)/(k+b)` over a grid of small integer/rational
     `a,b` (hundreds of cases) — Stage 1/2;
   - polynomial products `∏ P(k)` for random low-degree integer `P` — Stage 2;
   - geometric `∏ p(k) r^k` over small `r` and low-degree `p` — Stage 3;
   - telescoping ladders `∏ h(k+1)/h(k)` for random factored `h` — Stage 1
     (answer must be Gamma-free);
   - `q`-products once Stage 5 lands.
   For each: assert `closed_form /. n->{0..8}` equals the direct product, and
   assert the closed form is **stable** (`Method->Name` == `Automatic`).
2. **Known-hard / regression corpus** — a curated list with expected canonical
   forms: the `PRODUCT_PLAN.md` §7 benchmark table, plus divergence cases, empty
   ranges, half-integer shifts, and the symbolic-exponent hang guards.

**Stress-run protocol** (memory `project_crc_corpus_run_crashes_harness`,
`project_pre_existing_test_crashes`, `feedback_no_background_mathilda_pollers`):
build and run the stress binary **foreground only**, redirect to `/tmp`, with a
long timeout (macOS has no `timeout` — use a generous foreground run); **never
background** it or write poller loops (prior sessions OOM'd / flooded tmpfs that
way). A correct-but-unsimplifiable answer should surface as a reported oracle
mismatch with both sides printed, not a hang.

---

## 10. Cross-cutting `Product` options (wire in `product.c`)

| Option | Default | Owner stage |
|---|---|---|
| `Method -> "..."` | `Automatic` | Stage 0 (extend names per stage) |
| `VerifyConvergence -> True\|False` | `True` | Stage 4 |
| `GenerateConditions -> True\|False` | `False` | Stage 4 |
| `Assumptions -> ...` | `$Assumptions` | Stage 2/4 (passed to Simplify/Limit) |
| `GeneratedParameters -> C` | `None` | Stage 1/2 (indefinite-product constant) |

`N[unevaluated Product]` should route to `NProduct` (already shipped) — note here
so the `N` interaction is designed for, exactly as `Sum`/`NSum`.

---

## 11. Verification protocol (per stage)

A stage is **done** only when:
1. its new `test_product.c` cases pass (`check()` for canonical forms; `oracle()`
   for Gamma/Pochhammer forms at `n ∈ {0..8}`);
2. the **stress** families for that stage pass (oracle-verified, foreground run);
3. prior stages' suites still pass — run only the affected `*_tests` binaries
   (`product_tests`, `product_stress_tests`, and `sum_tests` if `Sum` was touched
   by Stage 3's exponent routing) — memory `feedback_scope_tests_to_change`;
4. `valgrind --leak-check=full` is clean on the new paths — diff against the
   known macOS baseline noise (~12.8 kB / 400 blocks from dyld/Accelerate;
   memory `project_macos_valgrind_baseline_noise`); watch the
   `evaluate`-doesn't-free pitfall (constraint 3.8);
5. docs updated: `docs/spec/builtins/calculus.md` + the weekly changelog
   `docs/spec/changelog/2026-06-15.md` (and the `Mathilda_spec.md` changelog-table
   row).

**Build & run (foreground, one process at a time):**
```bash
make -j$(sysctl -n hw.ncpu)                      # ./Mathilda (auto-discovers src/product/*.c)
cd tests && mkdir -p build && cd build
cmake .. >/dev/null && make product_tests product_stress_tests sum_tests -j$(sysctl -n hw.ncpu)
./product_tests
./product_stress_tests > /tmp/prod_stress.out 2>&1   # foreground; never background
./sum_tests                                          # regression (Stage 3 touches Sum)
```

---

## 12. Risks & mitigations

- **`Together`/`Factor` hang on `base^var`** (constraint 3.1) — the dominant
  risk, since Stages 3/5 live on such factors. Mitigation: the
  `prod_has_symbolic_power` guard at the top of every stage that might `Factor`,
  plus a hang-regression test under a wall-clock guard.
- **Empty / reversed / zero-straddling ranges** — `fold_times`→`1`, telescoping
  shortcut guarded `imin<=imax+1`; a zero factor in a finite range is caught by
  expansion. Test `n=0` in every oracle.
- **Non-canonical Gamma/Pochhammer output** — use the multiplicative oracle, not
  string equality (constraint 3.5).
- **Infinite-product false convergence** — the degree/leading/next-coefficient
  gate (Stage 4) rejects `a_k ~ c/k`; mirror `Sum::div`. Don't claim a value
  without the gate.
- **Leaks** — follow the `builtin_sum` ownership pattern exactly; every
  constructed node paired with `expr_free`/`eval`-and-free; never `expr_free(res)`.
- **Stress-run harness crashes** — foreground only, redirect to `/tmp`, no
  backgrounding/pollers (§9.2).

---

## 13. Suggested landing order

Stages **0 → 1 → 2 → 3** are self-contained on existing machinery (`Pochhammer`/
`Gamma`/`FactorialPower` all evaluate today) and deliver the bulk of practical
`Product`. **Stage 4** adds infinite products (needs the `.m` table + `Limit` +
the convergence gate). **Stage 5** (`q`) and **Stage 6** (Barnes-G/Glaisher) each
land only after their prerequisite builtins (§7). **Stage 7** is the
provably-minimal refinement. Each stage ships with its unit + stress cases and a
changelog note before the next begins.
