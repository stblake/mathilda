# Sum — Development Plan (remaining stages)

Status as of 2026-05-30: **Stages 0–3 are implemented and merged**
(branch `feature/sum-stages-0-3`). This document specifies the remaining
work — the deferred summation stages and the special-function builtins they
require — so each can be picked up independently.

---

## 1. Where we are

`Sum` lives in `src/sum/`: a dispatcher (`sum.c`) plus one file per algorithm,
each algorithm also exposed as a context-qualified builtin (`Sum`Polynomial`,
`Sum`Geometric`, `Sum`Gosper`) — the same shape `Integrate` uses for
`Integrate`BronsteinRational`. `Sum` is `HoldAll | Protected`.

**Shipped (Stages 0–3, no new special functions required):**

| Stage | File | Handles |
|------|------|---------|
| 0 | `sum.c` | iterator surface, finite/list expansion, multi-iterator nesting, Method cascade, indefinite routing |
| 1 | `sum_polynomial.c` | polynomial sums (Newton falling-factorial) |
| 2 | `sum_geometric.c` | `p(i) r^i` (geometric / polynomial-exponential) |
| 3 | `sum_gosper.c` | Gosper indefinite hypergeometric summation + `DifferenceDelta` |

**Shared helpers already available** (`src/sum/sum_internal.h`, `sum_util.c`):
`sum_eval(head,args,n)`, `sum_subst(e,var,val)`, `sum_factor(e)`,
`sum_free_of(e,var)`, `sum_sub(a,b)`, `sum_int(v)`,
`sum_stage_args(res,&f,&var,&imin,&imax,&definite)`.

**Adding a stage is additive**: new `src/sum/sum_<name>.c`, one `try_*` line in
`dispatch_def`/`dispatch_indef` (`sum.c`), one `sum_<name>_init()` call in
`sum_init()`. The cascade tolerates an absent stage (an unresolved
context-builtin call counts as "fall through").

---

## 2. Hard-won constraints (read before coding any stage)

These were discovered building Stages 0–3 and apply to every later stage:

1. **`Together` / `Factor` infinite-loop on a symbolic exponential** —
   any expression containing `Power[sym, var]` (e.g. `a^i`). `Cancel` and plain
   `evaluate` are safe. **Rule:** never run Together/Factor on an expression
   carrying `r^var`; normalise only the rational *coefficient* and multiply the
   exponential back in. See `sum_geometric.c`. (Memory:
   `project_together_factor_hang_exponential`.)
2. **`Solve` cannot solve symbolic linear *systems*** (e.g. 3 equations in
   `{c0,c1,c2}` with a symbolic parameter). Single-equation `Solve` and
   `SolveAlways` work. Stages 1–3 instead use either triangular back-substitution
   (`sum_geometric.c`) or `SolveAlways` (`sum_gosper.c`).
3. **`SolveAlways[eqn, var]` solves for *all* non-`var` symbols** — including any
   stray symbolic parameters in the coefficients, not just the unknowns you
   intend. Fine for parameter-free inputs; for parametric inputs prefer building
   the coefficient system explicitly (`CoefficientList`) and the triangular
   elimination pattern from `sum_geometric.c:solve_q`.
4. **`get_degree_poly` needs `Expand` first** (it does not expand `(i+3)^5`).
5. **`Simplify` reduces factorial ratios** (`(k+1)!/k! → k+1`) but **not
   binomial ratios** (`Binomial[n,k+1]/Binomial[n,k]` is left alone). This blocks
   Gosper/Zeilberger on binomial terms until a hypergeometric term-ratio
   normaliser or `FunctionExpand` exists (see §4 Stage 9 and §3 `FunctionExpand`).
6. **Finite numeric ranges never reach the cascade** — Stage 0 expands them
   directly. Closed-form stages only see symbolic bounds, `Infinity`, or the
   indefinite form, so they may assume non-numeric limits.
7. **Output need not be maximally combined.** Definite results may print as a
   difference of two terms (e.g. `1 - 1/(1+n)` rather than `n/(1+n)`); the test
   oracle checks `Simplify[lhs - rhs] == 0`, not string equality.

---

## 3. Prerequisite special-function builtins

These are **general builtins**, not part of `src/sum/` — give each its own
`src/<name>.c` (or a shared `src/specialfn/`) registered from `core_init()`.
Each needs: forms, attributes, `symtab_set_docstring` (terse — examples go in
`docs/spec/`), exact/symbolic simplifications, derivative (`D`) rules (add to
`src/internal/deriv.m` or native `deriv.c`), and machine/MPFR numeric evaluation
(`N`) following the pattern in `feedback_numeric_builtins_cover_mpfr` (cover
`EXPR_REAL` *and* `EXPR_MPFR`, add an `N[x, 35]` test).

Inventory (verified 2026-05-30):

| Builtin | Status | Needed by |
|---------|--------|-----------|
| `TrigToExp`, `ExpToTrig` | **exist, functional** | Stage 6 |
| `EulerGamma`, `Gamma` | symbol only (no evaluation) | Stages 4, 5, 8 |
| `PolyGamma`, `HarmonicNumber` | **missing** | Stage 5 |
| `Zeta` (Hurwitz 3-arg), `LerchPhi`, `BernoulliB`, `Pochhammer`, `PolyLog`, `FunctionExpand` | **missing** | Stages 8, 9, 10 |

### Group A — for Stage 5 (rational sums)

**`EulerGamma`** — Euler–Mascheroni constant γ. Already an interned symbol;
make it a *protected constant* (like `Pi`/`E`): `Constant` attribute, numeric
value under `N` (MPFR: `mpfr_const_euler`). No arguments.

**`PolyGamma[z]` and `PolyGamma[n, z]`** — `ψ(z) = Γ'(z)/Γ(z)` and its `n`-th
derivative `ψ⁽ⁿ⁾(z)`. Attributes `Listable | NumericFunction | Protected`.
Key relations to implement:
- Recurrence: `PolyGamma[n, z+1] == PolyGamma[n, z] + (-1)^n n! / z^(n+1)`
  (used to normalise integer shifts; do **not** auto-apply — only on demand).
- Special values: `PolyGamma[0, 1] == -EulerGamma`;
  `PolyGamma[n, 1] == (-1)^(n+1) n! Zeta[n+1]` (needs `Zeta`, Group B);
  `PolyGamma[0, 1/2] == -EulerGamma - 2 Log[2]`.
- Derivative: `D[PolyGamma[n, z], z] == PolyGamma[n+1, z]`.
- `N`: digamma/polygamma via MPFR (`mpfr_digamma` for `n==0`; series/reflection
  for `n>=1`).

**`HarmonicNumber[n]` and `HarmonicNumber[n, r]`** — `H_n = Σ_{k=1}^n 1/k` and
the generalised `H_n^(r) = Σ_{k=1}^n 1/k^r`. Attributes
`Listable | NumericFunction | Protected`. Implement as thin wrappers once
`PolyGamma`/`Zeta` exist:
- `HarmonicNumber[n] == EulerGamma + PolyGamma[0, n+1]`.
- `HarmonicNumber[n, r] == Zeta[r] - Zeta[r, n+1]` (needs Hurwitz `Zeta`).
- Exact rationals for integer `n` (`HarmonicNumber[4] == 25/12`).

### Group B — for Stages 8, 9, 10

**`Zeta[s]` and `Zeta[s, a]`** — Riemann ζ and Hurwitz `ζ(s,a) = Σ_{k≥0}
(k+a)^{-s}`. Attributes `Listable | NumericFunction | Protected`. Needed:
- `Zeta[2] == Pi^2/6`, `Zeta[4] == Pi^4/90`, even `Zeta[2n]` via `BernoulliB`:
  `Zeta[2n] == (-1)^(n+1) BernoulliB[2n] (2 Pi)^(2n) / (2 (2n)!)`.
- `Zeta[-n] == -BernoulliB[n+1]/(n+1)` (regularization, Stage 10).
- `Zeta[s, 1] == Zeta[s]`; recurrence `Zeta[s, a+1] == Zeta[s, a] - a^(-s)`.
- `D[Zeta[s, a], a] == -s Zeta[s+1, a]`.
- `N`: MPFR `mpfr_zeta` (Riemann); Hurwitz via Euler–Maclaurin.

**`BernoulliB[n]` and `BernoulliB[n, x]`** — Bernoulli numbers (exact rationals)
and Bernoulli polynomials. Attributes `Listable | NumericFunction | Protected`.
`BernoulliB[0]==1`, `BernoulliB[1]==-1/2`, `BernoulliB[2]==1/6`, odd `>1` are 0.
Compute via the recurrence `Σ_{k=0}^{n} Binomial[n+1,k] B_k == 0`. Enables a
tidy Faulhaber form for Stage 1 and the `Zeta` even-value table.

**`Pochhammer[a, n]`** — rising factorial `a(a+1)…(a+n-1) == Gamma[a+n]/Gamma[a]`.
Attributes `Listable | NumericFunction | Protected`. Exact for integer `n`;
needed for hypergeometric output (Stage 9) and as the WL-canonical form for
several closed forms.

**`Gamma[z]`** (make the existing symbol evaluate) — `Gamma[n] == (n-1)!`,
recurrence `Gamma[z+1] == z Gamma[z]`, `Gamma[1/2] == Sqrt[Pi]`,
`D[Gamma[z], z] == Gamma[z] PolyGamma[0, z]`, `N` via `mpfr_gamma`. Needed by the
logarithmic sums (`Sum[Log[i+1], i] == Log[Gamma[1+i]]`, Stage 5/8 boundary) and
Pochhammer.

**`LerchPhi[z, s, a]`** — `Φ(z,s,a) = Σ_{k≥0} z^k/(k+a)^s`. Attributes
`Listable | NumericFunction | Protected`. Generalises `Zeta` (`z==1`) and
`PolyLog` (`a==1`). Needed for rational-exponential sums and many infinite sums.
`N` via the defining series with convergence acceleration.

**`PolyLog[n, z]`** — `Σ_{k≥1} z^k/k^n == z LerchPhi[z, n, 1]`. Needed for some
infinite rational-exponential and Euler-sum results.

**`FunctionExpand`** — at minimum, reduce ratios of `Binomial`/`Pochhammer`/
`Gamma` to rational functions (e.g. `Binomial[n,k+1]/Binomial[n,k] →
(n-k)/(k+1)`). This is the missing piece that unblocks Gosper/Zeilberger on
binomial terms (§2.5). Could alternatively be a focused internal
"hypergeometric term ratio" normaliser used only inside `sum_gosper.c` /
`sum_zeilberger.c`.

---

## 4. Remaining summation stages

Dependency graph:

```
Stage 4 (PolyGamma, HarmonicNumber, EulerGamma) ─┐
                                                  ├─► Stage 5 (Sum`Rational)
Stage 3 (Gosper, shipped) ───────────────────────┘
Stage 2/3 + TrigToExp/ExpToTrig ───────────────────► Stage 6 (Sum`Trig)
Stage 7 (Zeta, LerchPhi, BernoulliB, Pochhammer, Gamma) ┐
Stage 5 + Limit ────────────────────────────────────────├─► Stage 8 (Sum`Infinite)
Stage 3 + Stage 7 ──────────────────────────────────────┴─► Stage 9 (Sum`Zeilberger)
Stage 7 + Stage 8 ─────────────────────────────────────────► Stage 10 (Regularization)
```

### Stage 5 — `Sum`Rational` (`sum_rational.c`)

**Status (2026-06-19): infinite path shipped.** `src/sum/sum_rational.c`
(`Method -> "Rational"`, in the cascade before `Sum`Hypergeometric`) closes
`Sum[p(i)/q(i), {i, imin, Infinity}]` for concrete integer `imin`,
`deg q >= deg p + 2`: linear/rational poles → `Zeta`/`PolyGamma`,
complex-conjugate quadratics → `Coth` + conjugate digamma, real radical roots →
`PolyGamma` over an auto-detected field extension. See
[`tasks/sum_rational.md`](tasks/sum_rational.md) and the 2026-06-15 changelog.
Still open: symbolic-coefficient denominators (direct-residue route, Phase D) and
the **finite/indefinite** rational path below (HarmonicNumber targets, Phase E).

Finite/indefinite rational `f = p(i)/q(i)`. Pipeline:
1. Let Gosper (Stage 3) try first — it already returns the rational
   antidifference when one exists; only the non-Gosper-summable remainder needs
   `PolyGamma`.
2. Partial-fraction decompose with the existing `Apart`. Each term
   `c/(i+a)^m` integrates to a polygamma:
   - `Σ_i 1/(i+a) = PolyGamma[0, i+a]` (+ the indefinite constant);
   - `Σ_i 1/(i+a)^m = (-1)^m/(m-1)! · PolyGamma[m-1, i+a]` for `m ≥ 2`;
   - integer `a` collapses to `HarmonicNumber` (`HarmonicNumber[i-1]`,
     `HarmonicNumber[i-1, m]`).
3. Definite = `F(imax+1) - F(imin)`; reuse the `sum_stage_args` /
   `sum_sub` / `Cancel` output pattern.
4. **Irreducible denominators of degree ≥ 2** (e.g. `1/(i^2+1)`): WL emits
   `RootSum[... PolyGamma ...]`. Handle the linear and (real-root) quadratic
   cases first; gate the general `RootSum` form behind a follow-up — Mathilda
   already has `RootSum` (`src/root.c`), so the eventual form is reachable.

Targets: `Sum[1/(i(i+6)), {i,1,n}]`, `Sum[1/i, {i,1,n}] → HarmonicNumber[n]`,
`Sum[1/i^2, {i,1,n}] → HarmonicNumber[n,2]`. Oracle as in Stages 1–3.

### Stage 6 — `Sum`Trig` (`sum_trig.c`)

Trigonometric polynomials in `Sin[k i]`, `Cos[k i]`. `TrigToExp`/`ExpToTrig`
already exist. Pipeline: `TrigToExp[f]` → a sum of `p(i) e^{c i}` terms →
Stage 2 (`Sum`Geometric`) / Stage 3 over base `e^{c}` → `ExpToTrig` and
`ComplexExpand` to recollect real trig. Watch §2.1 (keep the exponential factored;
no Together/Factor). Targets: `Sum[Sin[a i + b], i]`,
`Sum[Sin[i+1] Cos[i], i]`, `Sum[i Sin[i], {i,1,n}]`.

### Stage 8 — `Sum`Infinite` (`sum_infinite.c` or integrated into each `try_*`)

`imax == Infinity`. Take the symbolic `Limit` (`src/calculus/limit.c`) of the
finite closed form found by Stages 1–6 as `imax → ∞`, then recognise the
residual:
- geometric tails (`|r| < 1`) → `r^a/(1-r)`;
- `Σ 1/i^s → Zeta[s]`, `Σ 1/(i+a)^s → Zeta[s, a]`;
- `Σ z^i/i^s → LerchPhi[z, s, 1]` / `PolyLog`.
Add the options `VerifyConvergence` (default `True`: detect divergence, emit
`Sum::div`, stay unevaluated) and `GenerateConditions` (return
`expr if condition`). Targets: `Sum[1/i^2, {i,1,Infinity}] → Pi^2/6`,
`Sum[x^i, {i,0,Infinity}] → 1/(1-x)`, `Sum[x^i/i, {i,1,Infinity}] → -Log[1-x]`.

### Stage 9 — `Sum`Zeilberger` (`sum_zeilberger.c`)

Creative telescoping for parametric definite hypergeometric sums
`Σ_i F(n,i)` (binomial identities). Reuses the Stage-3 Gosper core,
parametrised in `n`: find a recurrence `Σ_j a_j(n) S(n+j) = 0` via a Gosper run
on `F(n,i)`, then solve/telescope to closed form (often `Pochhammer`/`Gamma`).
**Blocked on** the binomial term-ratio reduction (§2.5 / `FunctionExpand`) and
Group B specials. Targets: `Sum[Binomial[n,k], {k,0,n}] → 2^n`,
`Sum[Binomial[n,k]^2, {k,0,n}] → Binomial[2n,n]`.

### Stage 10 — Regularization (`sum_regularize.c`)

`Regularization -> "Abel" | "Borel" | "Cesaro" | "Dirichlet" | "Euler"` for
divergent sums. Mostly a wrapper that swaps which special-value/limit procedure
is consulted (`Zeta`/`BernoulliB` analytic continuation). Lowest priority.
Targets: `Sum[(-1)^i, {i,0,Infinity}, Regularization->"Abel"] → 1/2`.

---

## 5. Cross-cutting `Sum` options (wire in `sum.c`)

The Stage-0 dispatcher already strips trailing `Rule`/`RuleDelayed` options and
parses `Method`. Extend the option parser as the stages that need them land:

| Option | Default | Owner stage |
|--------|---------|-------------|
| `Method -> "..."` | `Automatic` | done (extend names per stage) |
| `Assumptions -> ...` | `$Assumptions` | Stage 5/8 (passed to Simplify/Limit) |
| `GenerateConditions -> True\|False` | `False` | Stage 8 |
| `GeneratedParameters -> C` | `None` | Stage 1/5 (indefinite constant naming) |
| `Regularization -> "..."` | `None` | Stage 10 |
| `VerifyConvergence -> True\|False` | `True` | Stage 8 |

`N[unevaluated Sum]` should route to numeric summation (`NSum`) — a separate
future builtin, noted here so the option/`N` interaction is designed for.

---

## 6. Verification (per stage)

Mirror `tests/test_sum.c` (add cases there; it already links the whole
`src/sum/` set via `tests/CMakeLists.txt` `COMMON_SRC`):
- exact-string `check()` for canonical forms;
- `same(a, b)` oracle = `Simplify[(a) - (b)] == 0`, robust to output form;
- **finite-expansion oracle**: every closed form must equal the Stage-0 direct
  expansion at several integer bounds (substitute `n -> 7`, compare);
- for each new special function, an `N[x, 35]` MPFR test
  (`feedback_numeric_builtins_cover_mpfr`);
- run only the affected `*_tests` binaries (`feedback_scope_tests_to_change`);
- `valgrind --leak-check=full` the new paths — watch the
  `evaluate(expr_new_function(...))` pitfall: `evaluate` does **not** free its
  input, so bind the node to a temp and `expr_free` it.

Each stage is "done" when its new cases pass, the prior stages' suites still
pass, and the doc + weekly changelog are updated (`docs/spec/builtins/calculus.md`,
`docs/spec/changelog/<Monday>.md`).
