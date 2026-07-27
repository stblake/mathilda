/* Compile[] engine performance gate.
 *
 * Two things are measured, and only one of them is a pass/fail gate:
 *
 *  1. **Did the body compile at all?**  This is the gate that matters most.  A
 *     `compile_expr` returning NULL is indistinguishable from success at the
 *     call site — the caller just quietly interprets — so a coverage gap reads
 *     as a 10-40x slowdown that looks exactly like working code.  Every body
 *     here asserts it compiled; a regression in the compilable subset fails the
 *     build rather than silently costing an order of magnitude.
 *
 *  2. **Optimiser effect**, as the ratio (unoptimised time / optimised time) for
 *     the same body.  A ratio is machine-independent in a way absolute timings
 *     are not, so it can be gated: the optimiser must never make a body slower
 *     than the code it started from (allowing 10% for measurement noise).
 *
 * Absolute ns/call figures are printed for information only — they are wildly
 * machine-dependent, and on this project's hardware end-to-end wall clock is
 * +-40% run to run, so every timing here is a MINIMUM over repeated trials
 * (a minimum is the statistic that noise cannot inflate) and warm.
 *
 * Run: ./bench_compile
 */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "parse.h"
#include "sym_intern.h"
#include "ndarray.h"
#include "compile/compile.h"
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int failures = 0;

/* WALL clock, not `clock()`.
 *
 * `clock()` returns CPU time summed over every thread, so a parallel region that
 * scales PERFECTLY reports as N times SLOWER.  That is exactly what happened
 * here: the threaded fused map measured 0.56-0.83x against the serial one while
 * the region itself was running 6.4x faster.  Any benchmark that a threaded code
 * path can reach has to be wall-clock, and the rest of tests/bench_*.c already
 * is — this file was the exception. */
static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1.0e-9;
}

/* Minimum-of-`reps` wall time for `iters` calls of a compiled all-real program. */
static double time_real(const CompiledProgram* p, size_t nargs, int iters, int reps) {
    double best = 1e300, a[4] = { 0.37, 0.91, 1.23, 0.55 };
    for (int r = 0; r < reps; r++) {
        double t0 = now_s(), acc = 0;
        for (int i = 0; i < iters; i++) {
            a[0] = 0.37 + 1e-9 * i;
            double o;
            if (compiled_eval_real(p, a, &o)) acc += o;
        }
        double t = now_s() - t0;
        if (t < best) best = t;
        (void)acc; (void)nargs;
    }
    return best / iters;
}

/* Compile one body twice — with and without the optimiser — and report the
 * ratio.  `raw` parses without evaluating, which loop bodies need: with the
 * arguments still free symbols the evaluator would close-form the Sum before the
 * compiler ever saw it. */
static void bench(const char* name, const char* src, const char* const* names,
                  size_t nargs, bool raw, int iters) {
    const char* inm[4];
    for (size_t k = 0; k < nargs; k++) inm[k] = intern_symbol(names[k]);
    const CompileType RR[4] = { CT_REAL, CT_REAL, CT_REAL, CT_REAL };

    Expr* b = raw ? parse_expression(src) : eval_and_free(parse_expression(src));
    CompiledProgram* po = compile_expr_ex(b, inm, RR, nargs, 0u);
    CompiledProgram* pr = compile_expr_ex(b, inm, RR, nargs, COMPILE_NO_OPT);

    if (!po || !pr) {
        /* The gate. Not "slow" — not compiled, which is far worse and silent. */
        printf("FAIL: %-26s DID NOT COMPILE (body outside the compilable subset)\n", name);
        failures++;
        compiled_free(po); compiled_free(pr); expr_free(b);
        return;
    }

    double to = time_real(po, nargs, iters, 5);
    double tr = time_real(pr, nargs, iters, 5);
    double ratio = tr / to;
    size_t no = compiled_num_instructions(po), nr = compiled_num_instructions(pr);

    /* Only gate bodies the optimiser actually changed: where it removed nothing
     * the ratio is measuring run-to-run noise, and gating on noise produces a
     * test that fails at random, which is worse than no test. */
    if (no < nr && ratio < 0.85) {
        printf("FAIL: %-26s optimiser made it SLOWER (%.2fx)\n", name, ratio);
        failures++;
    } else {
        printf("ok:   %-26s %7.1f ns/call  opt %.2fx  instrs %zu -> %zu  cse %zu\n",
               name, to * 1e9, ratio, nr, no, compiled_num_cse(po));
    }
    compiled_free(po); compiled_free(pr); expr_free(b);
}

/* Array bodies: the interesting axis is length, because the delegated array path
 * makes one full pass and one temporary buffer per operation, so its cost per
 * element does not fall with length the way a fused loop's would. */
static double time_arr(CompiledProgram* p, Expr* v, CompileType at, int iters) {
    double best = 1e300;
    for (int r = 0; r < 5; r++) {
        double t0 = now_s();
        for (int i = 0; i < iters; i++) {
            CompileValue av, out;
            av.type = at; av.v.a = v;
            if (compiled_eval(p, &av, &out)) {
                if (CT_IS_ARRAY(out.type)) expr_free(out.v.a);
            }
        }
        double t = now_s() - t0;
        if (t < best) best = t;
    }
    return best / iters;
}

/* Threaded fused map vs the SAME fused map on one thread.  Separate from
 * bench_arr because that one measures fusion against delegation, and the
 * delegated ND path threads too — so at large lengths it would compare two
 * threaded implementations and report ~1x while hiding whether either scaled. */
static void bench_par(const char* name, const char* src, size_t len, int iters) {
    const char* inm[1] = { intern_symbol("v") };
    const CompileType AT[1] = { CT_ARRAY(CT_REAL, 1) };
    Expr* b = parse_expression(src);
    CompiledProgram* pp = compile_expr_ex(b, inm, AT, 1, 0u);
    CompiledProgram* ps = compile_expr_ex(b, inm, AT, 1, COMPILE_NO_PAR);
    if (!pp || !ps) {
        printf("FAIL: %-26s DID NOT COMPILE\n", name);
        failures++; compiled_free(pp); compiled_free(ps); expr_free(b); return;
    }
    /* Identical programs mean the fan-out marker never got emitted, i.e. this
     * benchmark is timing the same thing twice and its ratio means nothing. */
    if (compiled_num_instructions(pp) <= compiled_num_instructions(ps)) {
        printf("FAIL: %-26s len=%zu: parallel fan-out did not engage\n", name, len);
        failures++;
    }
    int64_t dims[1]; dims[0] = (int64_t)len;
    double* buf = malloc(len * sizeof(double));
    for (size_t i = 0; i < len; i++) buf[i] = 0.3 + 1.7 * (double)i / (double)len;
    Expr* v = expr_new_ndarray(1, dims, buf, NDT_FLOAT64);

    double ts = time_arr(ps, v, AT[0], iters);
    double tp = time_arr(pp, v, AT[0], iters);
    /* No absolute speed gate: the ratio is bounded by the core count, which is a
     * property of the machine, not of the compiler.  But threading must never
     * make a body SLOWER — that is machine-independent, and it is what a
     * reintroduced CPU-time clock, a lock in the worker path or a contention bug
     * would each look like.  (On a single-core box no fan-out happens at all and
     * the ratio sits at 1.0, so the gate stays quiet rather than flaky.) */
    double ratio = ts / tp;
    if (ratio < 0.9) {
        printf("FAIL: %-26s len=%zu: threaded is SLOWER (%.2fx)\n", name, len, ratio);
        failures++;
    }
    printf("ok:   %-26s len=%-8zu threaded %8.1f us (%5.2f ns/el)  serial %8.1f us  -> %.2fx\n",
           name, len, tp * 1e6, tp / (double)len * 1e9, ts * 1e6, ratio);
    expr_free(v);
    compiled_free(pp); compiled_free(ps); expr_free(b);
}

/* Strip-mined fusion vs the delegated NDArray path, same body, same data.  This
 * ratio is the whole question for array work: delegation makes one full-length
 * pass and one temporary buffer per operation, fusion makes one pass total. */
static void bench_arr(const char* name, const char* src, size_t len, int iters) {
    const char* inm[1] = { intern_symbol("v") };
    const CompileType AT[1] = { CT_ARRAY(CT_REAL, 1) };
    Expr* b = parse_expression(src);                 /* never pre-evaluate */
    CompiledProgram* pf = compile_expr_ex(b, inm, AT, 1, 0u);
    CompiledProgram* pd = compile_expr_ex(b, inm, AT, 1, COMPILE_NO_FUSE);
    if (!pf || !pd) {
        printf("FAIL: %-26s DID NOT COMPILE (%s)\n", name, !pf ? "fused" : "delegated");
        failures++;
        compiled_free(pf); compiled_free(pd); expr_free(b);
        return;
    }
    int64_t dims[1]; dims[0] = (int64_t)len;
    double* buf = malloc(len * sizeof(double));
    for (size_t i = 0; i < len; i++) buf[i] = 0.3 + 1.7 * (double)i / (double)len;
    Expr* v = expr_new_ndarray(1, dims, buf, NDT_FLOAT64);

    double tf = time_arr(pf, v, AT[0], iters);
    double td = time_arr(pd, v, AT[0], iters);
    size_t nf = compiled_num_instructions(pf), nd = compiled_num_instructions(pd);
    /* Identical programs mean fusion silently did NOT engage — the same class of
     * invisible failure as a body that never compiled, so it is a gate. */
    if (nf == nd) {
        printf("FAIL: %-26s len=%zu: fusion did not engage (%zu instrs both ways)\n",
               name, len, nf);
        failures++;
    }
    printf("ok:   %-26s len=%-6zu fused %7.2f us (%5.2f ns/el)  delegated %7.2f us  -> %.2fx  [%zu vs %zu instrs]\n",
           name, len, tf * 1e6, tf / (double)len * 1e9, td * 1e6, td / tf, nf, nd);
    expr_free(v);
    compiled_free(pf); compiled_free(pd); expr_free(b);
}

int main(void) {
    core_init();
    const char* xyz[] = { "x", "y", "z" };

    printf("=== Compile[] scalar bodies (optimiser A/B) ===\n");
    /* dispatch-bound: pure add/mul, so VM instruction dispatch dominates and
     * this is the body that sees encoding and superinstruction work. */
    {
        char buf[4096]; size_t pos = 0;
        const int DEG = 40;
        for (int i = 0; i < DEG; i++) pos += (size_t)snprintf(buf + pos, sizeof buf - pos, "(");
        pos += (size_t)snprintf(buf + pos, sizeof buf - pos, "x");
        for (int i = 0; i < DEG; i++)
            pos += (size_t)snprintf(buf + pos, sizeof buf - pos, " + %d) x", (i % 9) + 1);
        bench("horner deg40 (dispatch)", buf, xyz, 1, true, 500000);
    }
    bench("mixed libm", "Sin[x y] + Exp[-x^2/2] Cos[y] - Sqrt[Abs[x]] + x^3 - 2 y^2",
          xyz, 2, false, 400000);
    bench("shared subexpressions", "Sin[x y] + Cos[x y] Sin[x y] + Sin[x y]^3",
          xyz, 2, false, 400000);
    /* loop bodies: the only place loop-invariant code motion has anything to do.
     * The invariant work here is deliberately expensive (Exp, Sqrt) so hoisting
     * it out of 30 iterations is visible above the noise. */
    bench("loop, invariant libm", "Sum[Exp[-x^2] Sqrt[Abs[y]] + i, {i, 1, 30}]",
          xyz, 2, true, 60000);
    bench("loop, invariant arith", "Sum[(x + y) (x - y) i, {i, 1, 30}]",
          xyz, 2, true, 100000);
    bench("nested loops", "Sum[Sum[Sin[x] y i j, {j, 1, 10}], {i, 1, 10}]",
          xyz, 2, true, 20000);
    bench("Newton (While)",
          "Module[{t = x, k = 0}, While[k < 20, t = (t + x/t)/2; k = k + 1]; t]",
          xyz, 1, true, 150000);
    bench("Nest", "Nest[Function[u, (u + x/u)/2], x, 20]", xyz, 1, true, 150000);

    printf("\n=== Compile[] rank-1 array bodies ===\n");
    bench_arr("Total[Sin v Exp -v]", "Total[Sin[v] Exp[-v] + Sqrt[v]]", 16,    200000);
    bench_arr("Total[Sin v Exp -v]", "Total[Sin[v] Exp[-v] + Sqrt[v]]", 1024,  20000);
    bench_arr("Total[Sin v Exp -v]", "Total[Sin[v] Exp[-v] + Sqrt[v]]", 65536, 300);
    bench_arr("v^2 + 2 v",           "v^2 + 2 v + 1",                   65536, 500);
    /* Reducing bodies allocate NO result buffer, which separates the cost of
     * the elementwise pass itself from the cost of allocating and faulting in a
     * full-length output array. */
    bench_arr("Total[v^2 + 2v + 1]",  "Total[v^2 + 2 v + 1]",           65536, 500);
    bench_arr("Total[v w-ish]",       "Total[v + v v]",                 65536, 500);

    printf("\n=== Compile[] threaded fused map (OP_APAR) ===\n");
    bench_par("Sqrt v + v^2",         "Sqrt[v] + v^2",                   1000000, 60);
    bench_par("Sin v Exp -v",         "Sin[v] Exp[-v] + Sqrt[v]",        1000000, 40);
    bench_par("Gamma v",              "Gamma[v] + Erf[v]",               1000000, 20);

    if (failures == 0) printf("\nAll Compile benchmarks within gate.\n");
    else printf("\n%d Compile benchmark gate(s) FAILED.\n", failures);
    return failures ? 1 : 0;
}
