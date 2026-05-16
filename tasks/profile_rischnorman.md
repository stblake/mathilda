# Profiling `Integrate`RischNorman` — 2026-05-15

## What was done

1. **Benchmark harness** at `tests/bench_intrischnorman.c` — 52 integrands
   spanning Phases 4–7 of the test corpus.  Build target
   `bench_intrischnorman` in `tests/CMakeLists.txt`.
2. **Phase timers** added to `src/intrischnorman.c`, gated on
   `PMINT_PROFILE=1` env var.  Top-level (16 phases of
   `rischnorman_integrate`) plus sub-phases of `try_integral_full`
   and `solve_linear_undet`.  Always-on time accumulation;
   `clock_gettime(CLOCK_MONOTONIC)` overhead is ~30 ns / phase boundary
   so the cost when un-dumped is negligible.
3. **macOS `sample`** of `bench_intrischnorman --only
   p7-exp_x_plus_inv_log --reps 2` — the slowest case (~12 s) — gave
   function-level confirmation of where the cycles go.

## Headline numbers (52-integrand corpus, total ≈ 34 s)

```
top-level
  try_integral          32.82 s   97.1 %   <-- dominant
  output_cleanup         0.44 s    1.3 %
  Together(ff)           0.16 s    0.5 %
  ff_decompose           0.15 s    0.5 %
  build_candidate        0.07 s    0.2 %
  my_factors             0.06 s    0.2 %
  ... all others combined ~0.16 s, <0.5 %

try_integral breakdown (32.82 s)
  ti.solve_rowred       23.36 s   71.2 %   <-- RowReduce
  ti.term_build          4.53 s   13.8 %   <-- big Expand
  ti.solve_clist         1.76 s    5.4 %   <-- eval_expand+CoefficientList
  ti.apply_d             1.55 s    4.7 %
  ti.solve_walk          1.50 s    4.6 %   <-- walk_coefficient_table
  ti.substitute          0.09 s    0.3 %
  ti.solve_decode        0.01 s    0.0 %
```

Slowest 10 integrands (best-of-1 wall time):
- `p7-exp_x_plus_inv_log`     11.7 s
- `p7-sin_sin_2x`              6.8 s
- `p7-sin_over_x_plus_log_cos` 5.9 s
- `p7-exp_inv_log`             4.5 s
- `p7-exp_over_log_shifted`    1.2 s
- `p7-tan_div`                 0.66 s
- everything else < 0.34 s.

The slow tail (the 4 worst cases) accounts for ~85 % of total time;
each one hits a hard RowReduce.

## Function-level confirmation (sample output, `p7-exp_x_plus_inv_log`)

Top-of-stack hot leaves (out of ~14 500 samples):

| samples | function                                    |
|--------:|---------------------------------------------|
|   3 069 | `nanov2_find_block_and_allocate` (malloc)   |
|   2 757 | `_nanov2_free`                              |
|   2 498 | `nanov2_malloc_type`                        |
|   1 222 | `expr_contains_symbol`                      |
|     386 | `intern_symbol`                             |
|     384 | `_platform_strcmp`                          |
|     344 | `nanov2_calloc_type`                        |
|     251 | `evaluate_step`                             |
|     242 | `collect_symbols_in`                        |
|     233 | `expr_free`                                 |
|     220 | `expr_copy`                                 |

→ **~55 % of leaf time is allocator traffic.**

Call chain dominating that time:

```
builtin_rowreduce
  evaluate_step
    builtin_times          (per-cell arithmetic)
      multiply_numbers
        expr_bigint_normalize
          __gmpz_fits_slong_p / __gmpz_clear / nanov2_malloc / _nanov2_free
        expr_new_integer
          nanov2_malloc_type
```

Every matrix-cell operation goes through:
`evaluate_step` → attribute lookup → arg evaluation loop → Orderless
flatten/sort → `builtin_times` → numeric-contagion → `multiply_numbers`
→ bigint normalize → 2 allocations + 1 free, *per scalar multiply*.

## Concrete optimization targets, ranked by expected payoff

### 1. Bypass the Expr/evaluator path inside RowReduce (largest win)

The augmented matrix passed to `RowReduce` (called from
`solve_linear_undet` at intrischnorman.c:1737) has rational entries.
The current RowReduce in `linalg.c` performs every scalar op through
the symbolic evaluator — that is the source of the 23 s.

Two paths:

- **Quick win:** detect "all entries are integers or `Rational[]`",
  copy the matrix into `mpq_t[][]` (or `mpz_t[][]` if integer-only),
  run a fraction-free Bareiss elimination, then convert the RREF back
  to `Rational[]` Expr nodes.  pmint matrices are integer-coefficient
  by construction (the equation is built from polynomials over Q with
  symbolic A/B unknowns), so this short-circuit should fire on every
  rischnorman call.  Expected 10–50× on the slow tail.

- **Bigger lift:** a sparse RREF.  Inspecting an actual pmint matrix
  for the slow cases would tell us the density, but the equations are
  per-monomial × per-unknown — for `p7-exp_x_plus_inv_log` with
  several A's and several B's spread across many monomials, density
  is likely 10–30 %.  Sparse representation skips zero-pivot ops.

### 2. Skip the redundant Expand inside `solve_linear_undet` (small win)

intrischnorman.c:1619 calls `eval_expand(equation_numer)` again even
though `try_integral_full` already produced the equation as a fully
expanded polynomial at line 2068 (the final `eval_expand(mk_plus2(…))`
of term1/term2/term3).  Drop this — it shows up as 5.4 % of
`try_integral` (1.8 s).

### 3. Extract coefficients without round-tripping through Expr

`solve_linear_undet` builds a `CoefficientList` Expr (line 1620) and
then `walk_coefficient_table` recursively unboxes it back into rows
(line 1634).  Combined that is 10 % of `try_integral`.  A direct
multivariate-coefficient extractor that emits straight into the
`mpq_t**` matrix would eliminate both phases.  Pairs naturally with
optimization #1.

### 4. Pool / arena allocate short-lived Expr nodes inside try_integral

Allocator overhead is ~55 % of leaf time and most of those nodes die
within a single iterate.  Either:
- a per-call arena allocator (free-all-at-once when leaving
  `try_integral_full`), or
- a freelist for small fixed-size Expr nodes (integers + small
  Function nodes).

This is more invasive (touches `expr.c`) but would also help every
other heavy CAS path (Together, Cancel, Expand on big polynomials).

### 5. `expr_contains_symbol` micro-optimization (small)

8 % of leaf time, called from inside `builtin_times` / `collect_symbols_in`.
It walks the tree on every call; the symbols of interest are fixed
across an evaluation, so a `set<Symbol>` membership test instead of a
linear walk would help — but this is a small win compared to #1.

## Files touched

- `src/intrischnorman.c` — phase-timer infrastructure +
  instrumentation of `rischnorman_integrate`, `try_integral_full`,
  `solve_linear_undet`.  Always on; output gated by
  `PMINT_PROFILE=1`.
- `tests/bench_intrischnorman.c` — new bench harness.
- `tests/CMakeLists.txt` — wired in `bench_intrischnorman` target.

## How to re-run

```
cd tests/build
make -j8 bench_intrischnorman
./bench_intrischnorman                          # corpus, baseline
PMINT_PROFILE=1 ./bench_intrischnorman          # corpus + phase breakdown
./bench_intrischnorman --only p7-exp_x_plus_inv_log --reps 3
PMINT_PROFILE=1 ./bench_intrischnorman --only p7-exp_x_plus_inv_log
# for sample:
./bench_intrischnorman --only <case> --reps 3 &
sample $! 18 -mayDie -file /tmp/sample.txt
```

---

## Implementation #1 — direct mpq Gauss-Jordan (landed 2026-05-15)

`solve_linear_undet` now classifies the augmented matrix once
populated; if every entry is `Integer | BigInt | Rational[int, int]`
the elimination is performed directly on `mpq_t**` via Gauss-Jordan
with row-index indirection, bypassing the symbolic evaluator entirely.
BigInt-bearing numerators / denominators flow through `mpq_t`
unchanged.  Non-Q matrices (symbol-bearing entries like `a` from
`Exp[a x]`) fall through to the existing symbolic `RowReduce` path.

**All 94 `test_intrischnorman` cases pass.  `intrat_tests` clean.**

### Speedup (52-case bench, best-of-1, total wall time)

|                                | before  | after  | ratio   |
|--------------------------------|--------:|-------:|--------:|
| **Total**                      | 34.1 s  | 10.2 s | **3.35×** |
| p7-exp_x_plus_inv_log          | 11.74 s | 2.27 s | **5.2×**  |
| p7-sin_sin_2x                  |  6.79 s | 1.93 s | **3.5×**  |
| p7-sin_over_x_plus_log_cos     |  5.95 s | 1.74 s | **3.4×**  |
| p7-exp_inv_log                 |  4.44 s | 1.42 s | **3.1×**  |
| p7-exp_over_log_shifted        |  1.17 s | 0.37 s | **3.1×**  |
| p7-tan_div                     |  0.66 s | 0.23 s | **2.9×**  |

### Post-optimization `try_integral` breakdown (9.2 s, 54 calls)

```
ti.term_build      4.40 s   47.9 %      <-- new top (Expand-heavy)
ti.solve_clist     1.63 s   17.7 %      <-- CoefficientList + Expand
ti.apply_d         1.56 s   17.0 %
ti.solve_walk      1.47 s   15.9 %      <-- walk_coefficient_table
ti.solve_qfast     0.03 s    0.3 %      <-- the new fast solver
ti.solve_rowred    0.01 s    0.1 %      <-- 3/54 fall-throughs only
ti.substitute      0.09 s    1.0 %
```

RowReduce went from **71.2 % → 0.1 %** of `try_integral`.  The fast
path fires on 51/54 calls (94 %); the 3 fall-throughs are the
parameter-bearing integrand (`Exp[a x]`) whose matrix carries `a` as
a symbol.

### Next bottleneck

`ti.term_build` (eval_expand on `term1 / term2 / term3 /
equation_numer`) and `ti.solve_clist + ti.solve_walk` are now the
hot phases — these are the remaining items #2 and #3 from the
original ranked list above.

---

## Implementation #2+#3 — direct-from-expanded extractor (landed 2026-05-15)

Combined item #2 (drop redundant `eval_expand` inside
`solve_linear_undet`) and item #3 (direct multivariate coefficient
extractor bypassing `CoefficientList` + `walk_coefficient_table`).

The new path `try_solve_direct_q`:

1. Walks `Plus` children of the (already-expanded) `equation_numer`.
2. Decomposes each term into `(mpq scalar, var exponent vector,
   at-most-one unknown index)` via a small recursive `Times` walker.
3. Groups terms by monomial key into a sparse `QMatBuild` of `mpq_t`
   rows, accumulating column entries on the fly.
4. Prunes zero rows, flattens to a dense flat `mpq_t*` matrix and
   runs Gauss-Jordan (shared with #1 via the factored-out
   `gauss_jordan_mpq` + `emit_q_solution` helpers).

Bails to the previous `CoefficientList` path (which still has its
own Q-fast detector on `Expr**` rows) when any term carries a
foreign symbol — preserving the symbolic fallback for
parameter-bearing integrands.

**All 94 `test_intrischnorman` cases pass.  `intrat_tests` clean.**

### Speedup (52-case bench, best-of-1, total wall time)

|                                | baseline | after #1 | after #1+#2+#3 | total ratio |
|--------------------------------|---------:|---------:|---------------:|------------:|
| **Total**                      | 34.1 s   | 10.2 s   | **7.0 s**      | **4.87×**   |
| p7-exp_x_plus_inv_log          | 11.74 s  |  2.27 s  | **1.77 s**     | **6.6×**    |
| p7-sin_over_x_plus_log_cos     |  5.95 s  |  1.74 s  | **1.34 s**     | **4.4×**    |
| p7-exp_inv_log                 |  4.44 s  |  1.42 s  | **1.31 s**     | **3.4×**    |
| p7-sin_sin_2x                  |  6.79 s  |  1.93 s  | **0.84 s**     | **8.1×**    |
| p7-exp_over_log_shifted        |  1.17 s  |  0.37 s  | **0.17 s**     | **6.9×**    |
| p7-tan_div                     |  0.66 s  |  0.23 s  | **0.13 s**     | **5.1×**    |

### Post-#1+#2+#3 `try_integral` breakdown (6.02 s, 54 calls)

```
ti.term_build      4.30 s   71.5 %     <-- new top: Expand on term1/2/3
ti.apply_d         1.56 s   26.0 %     <-- derivatives
ti.substitute      0.09 s    1.5 %
ti.solve_directq   0.04 s    0.7 %     <-- new direct extractor (51/54)
ti.solve_clist     0.004 s   0.1 %     <-- 3/54 fall-throughs only
ti.solve_walk      0.004 s   0.1 %     <-- 3/54 fall-throughs only
ti.solve_rowred    0.008 s   0.1 %     <-- 3/54 fall-throughs only
ti.solve_qfast     0.0   s   0.0 %
```

Combined `solve_*` cost is now **~1%** of `try_integral`, down from
the original **92%** baseline (71.2 RowReduce + 13.8 + 5.4 + 4.6 +
…).  The fast path fires on 51/54 calls; the same 3 parameter-bearing
fall-throughs as #1.

### Next bottleneck

`ti.term_build` dominates at 71.5 % — the four `eval_expand` calls
on `term1 / term2 / term3 / equation_numer` inside
`try_integral_full`.  These are full Expand passes on products of
already-expanded polynomials.  Two avenues:

- **Targeted polynomial multiply**: implement
  `expand_poly_mult(a, b, vars)` that walks two expanded sums and
  emits the cross-product as a single Plus, bypassing the generic
  Expand pipeline (which routes every multiply through
  `evaluate_step → builtin_times → multiply_numbers`).
- **Combine term2 / term3**: share `prod_g` and `cden_sq` factors so
  fewer multiplications are needed.

`ti.apply_d` at 26 % is next; mostly Expand inside `apply_d` itself.

After those, the leaf-time allocator share (~55 % of sample
profile) suggests an arena allocator for short-lived Expr nodes
inside `try_integral_full` would still be a meaningful win.

---

## Implementation #4 — Structural distributor replaces term_build's Expand chain

### Drill-down

A second profile pass added sub-sub-phase timers to `term_build`:

```
ti.term_build           4.32 s   71.5 %   (try_integral)
  tb.dpart              2.69 s   44.5 %   ← d_num·cden − cand_num·d_cden
    tb.dpart_a          2.36 s   39.1 %   ← eval_expand(d_num · cden)  *the* hotspot
    tb.dpart_b          0.08 s    1.3 %
    tb.dpart_s          0.25 s    4.1 %
  tb.equation           0.97 s   16.0 %   ← eval_expand(term1 − term2 − term3)
  tb.term2              0.61 s   10.0 %   ← eval_expand(ff_denom · prod_g · dpart)
  tb.term3              0.045 s   0.7 %
  tb.term1              0.008 s   0.1 %
  tb.cden_sq            0.002 s   0.0 %
  tb.prodg              0.001 s   0.0 %
```

A single `eval_expand(d_num · cden)` was 39 % of `try_integral`.  All
the per-product Expand calls share the same root cause: the generic
evaluator pipeline canonicalises (sort, like-term merge, attribute
dispatch) that the downstream Q-linear solver does not need.

### Change

Replaced the entire term-build block in `try_integral_full` with a
**single structural distribution pass**:

1. Build the equation as one unevaluated Plus/Times tree
   `Plus[+term1, −term2a, +term2b, −term3]`, where each `term*` is a
   `Times[…]` of factors and `prod_g`, `cden_sq` are themselves Times
   chains (not pre-expanded).
2. A new `distribute_to_plus(e)` recursively distributes the tree:
   - `Plus`: concatenate distributed children.
   - `Times`: cartesian product of distributed factors via
     `mk_times_concat` (flat Times-of-atoms, no nested Plus).
   - `Power[Plus, k]` (small `k ≥ 0`): `k`-fold cartesian self-product.
   - Anything else: opaque atom — emitted as-is.
3. **No evaluator calls** — pure structural malloc/expr_copy.
4. The result is a flat `Plus[Times[…]…]` fed straight to
   `try_solve_direct_q`.  Same-base exponents (including negatives
   from logarithmic-derivative `l_k = 1/x` factors) combine naturally
   in `decompose_factor`'s `exp_vec[]` accumulator; the qmb keying
   tolerates negative entries.
5. Relaxed `decompose_factor` to accept negative integer `Power`
   exponents (previously bailed) — required because the structural
   path keeps `Power[x, -1] · Power[x, 2]` as two separate factors
   that combine *during* decomposition.

### Result (52-case corpus, best-of-3)

| | baseline | after #1 | after #1+2+3 | after #1+2+3+4 | total speedup |
|---|---:|---:|---:|---:|---:|
| Total | 34.1 s | 10.2 s | 7.0 s | **2.73 s** | **12.5×** |
| exp_x_plus_inv_log | 11.74 s | 2.27 s | 1.77 s | **0.67 s** | 17.5× |
| sin_over_x_plus_log_cos | 5.9 s | … | 1.37 s | **0.52 s** | 11.4× |
| exp_inv_log | 4.5 s | … | 1.25 s | **0.47 s** | 9.6× |
| sin_sin_2x | 6.79 s | 1.93 s | 0.85 s | **0.073 s** | 93× |

`term_build` collapsed from 4.32 s (71.5 % of try_integral) to **0.085 s
(1.6 %)** — a **51× per-call speedup** on the hottest phase.

### Updated phase breakdown (52 cases, reps=3, 8.33 s cumulative)

```
top-level                                                 share
  try_integral          5.39 s    64.7 %
  output_cleanup        1.31 s    15.7 %
  Together(ff)          0.47 s     5.6 %
  ff_decompose          0.47 s     5.7 %
  splitFactor(ff)       0.13 s     1.5 %
  my_factors            0.18 s     2.2 %
  ... rest < 1 %

try_integral breakdown (5.39 s)
  ti.apply_d            4.65 s    86.3 %   ← NEW dominant phase
  ti.solve_directq      0.27 s     5.1 %
  ti.substitute         0.27 s     5.1 %
  ti.term_build         0.085 s    1.6 %
  ti.solve_clist        0.046 s    0.9 %    (9/162 fall-throughs)
  ti.solve_walk         0.013 s    0.2 %    (9/162)
  ti.solve_rowred       0.023 s    0.4 %    (9/162)
  ti.solve_qfast        0      s   0.0 %
```

The fast path fires on **153/162** call sites; same 9 parameter-bearing
fall-throughs as #2+#3.  No correctness regressions: full 116-suite
test corpus (`intrischnorman_tests`, `intrat_tests`, plus the broader
Mathilda test suite) is green.

### Next bottleneck

`ti.apply_d` now dominates at 86 % of `try_integral`.  Its work is:

```c
sum = Σ_k l[k] · D[f, vars[k]];
return eval_and_free(sum);
```

Per call: `nv` invocations of `call_d` (which goes through the full
`evaluate` pipeline for `D[…]`), plus a final `eval_and_free`.  The
heavy callers are `apply_d(cand_num, …)` where `cand_num` has tens of
unknowns × tens of monomials.

Avenues for #5:

- **Targeted polynomial differentiator**: walk `cand_num` as a Plus
  of monomials directly, applying the power rule per factor instead
  of routing through `evaluate → builtin_d → DownValue lookup` for
  every term.
- **Skip the final `eval_and_free`**: `apply_d`'s consumers
  (`term_build`'s structural distributor, the `dpart` builder) no
  longer require canonical form — they can consume an unflattened
  `Plus[Times[l_k, …]…]` directly.

---

## Implementation #5 — Power-rule fast path + skip canonicalisation

Both #5 avenues implemented together:

1. **`diff_monomial_pow(m, var_k)`** / **`diff_poly_pow(f, var_k)`** —
   direct polynomial differentiator.  For each monomial (Times of
   integer/rational scalars, atoms, `Power[var, int]`):
   - Find the unique factor containing `var_k` (must be `var_k`
     itself or `Power[var_k, int]`).
   - All other factors must be free of `var_k`; otherwise bail and
     fall back to the generic `call_d` path.
   - Emit `c · e · v_k^(e-1) · (other factors)` directly.
   - Returns 0 when `var_k` is absent from `m`.

2. **`apply_d_raw`** — same fast-path as `apply_d` but **skips the
   final `eval_and_free`** canonicalisation pass.  Safe because
   `try_integral_full`'s structural distributor walks Plus / Times
   itself and does not require canonical input.

The generic `apply_d` (still called from `split_factor`'s
content / primitive-part recursion, where `eval_poly_gcd` needs
canonical input) retains the canonicalisation pass.

### Result (52-case corpus, best-of-3)

| | after #4 | after #5 | total speedup |
|---|---:|---:|---:|
| Total | 2.73 s | **1.23 s** | **27.7× vs baseline 34.1 s** |
| ti.apply_d | 4.65 s | **0.003 s** | >1000× |
| try_integral | 5.39 s | **0.76 s** | 7.1× |
| exp_x_plus_inv_log | 0.67 s | **0.081 s** | 8.3× |
| exp_inv_log | 0.47 s | **0.040 s** | 11.6× |

### Updated phase breakdown (52 cases, reps=3, 3.88 s cumulative)

```
top-level                                                 share
  output_cleanup        1.44 s    37.1 %   ← NEW dominant
  try_integral          0.76 s    19.6 %
  Together(ff)          0.49 s    12.6 %
  ff_decompose          0.48 s    12.4 %
  build_candidate       0.20 s     5.2 %
  my_factors            0.19 s     4.8 %

try_integral breakdown (0.76 s)
  ti.solve_directq      0.28 s    37.2 %
  ti.substitute         0.28 s    37.3 %
  ti.term_build         0.086 s   11.3 %
  ti.solve_clist        0.045 s    6.0 %   (9 fall-throughs)
  ti.solve_rowred       0.024 s    3.2 %
  ti.solve_walk         0.012 s    1.6 %
  ti.apply_d            0.003 s    0.3 %   ← essentially gone
```

After #5 the `try_integral` budget is well-balanced across solve /
substitute / term_build, and the bottleneck shifts out to the
top-level `output_cleanup` (per-summand `eval_cancel` + Sin/Cos
fresh-substitution round-trip).  No correctness regressions; all 116
test suites green.
