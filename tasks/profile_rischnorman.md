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
