/* Numeric RHS compiler for NDSolve: the compiled stack-machine must reproduce
 * the symbolic evaluator to rounding, and its colored finite-difference Jacobian
 * must match the analytic (symbolic-D) Jacobian.  We build small systems whose
 * components are a battery of expressions in the reduced-state symbols
 * NDSolve`w<k> and the time variable t, compile them, and compare compiled
 * evaluation / Jacobian against references over many random points.
 *
 * Soft asserts: prints FAIL and keeps going. Run: ./ndsolve_compile_tests */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "parse.h"
#include "sym_intern.h"
#include "numerical_calculus/ndsolve_common.h"
#include "numerical_calculus/ndsolve_compile.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

static int failures = 0;

/* deterministic xorshift RNG for reproducibility */
static uint64_t rng_state = 0x9e3779b97f4a7c15ULL;
static double urand(double lo, double hi) {
    uint64_t x = rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    rng_state = x;
    double u = (double)(x >> 11) * (1.0 / 9007199254740992.0);
    return lo + (hi - lo) * u;
}

/* build the state symbol NDSolve`w<k> */
static Expr* wsym(size_t k) {
    char buf[32];
    snprintf(buf, sizeof buf, "NDSolve`w%zu", k);
    return expr_new_symbol(intern_symbol(buf));
}

/* Reference: substitute w_k -> Y[k], t -> tval into f and evaluate to a double.
 * Returns false if the value is not a finite real (domain error / complex). */
static bool ref_eval(Expr* f, const double* Y, size_t d, double t, double* out) {
    Expr** rules = malloc((d + 1) * sizeof(Expr*));
    for (size_t k = 0; k < d; k++) {
        Expr* r[2] = { wsym(k), expr_new_real(Y[k]) };
        rules[k] = expr_new_function(expr_new_symbol("Rule"), r, 2);
    }
    Expr* rt[2] = { expr_new_symbol("t"), expr_new_real(t) };
    rules[d] = expr_new_function(expr_new_symbol("Rule"), rt, 2);
    Expr* rl = expr_new_function(expr_new_symbol("List"), rules, d + 1);
    free(rules);
    Expr* ra[2] = { expr_copy(f), rl };
    Expr* call = expr_new_function(expr_new_symbol("ReplaceAll"), ra, 2);
    Expr* res = eval_and_free(call);
    bool ok = false;
    if (res) {
        if (res->type == EXPR_REAL) { *out = res->data.real; ok = isfinite(*out); }
        else if (res->type == EXPR_INTEGER) { *out = (double)res->data.integer; ok = true; }
    }
    expr_free(res);
    return ok;
}

/* symbolic Jacobian entry D[f_i, w_j] evaluated at (Y, t) */
static bool ref_jac(Expr* fi, size_t j, const double* Y, size_t d, double t, double* out) {
    Expr* da[2] = { expr_copy(fi), wsym(j) };
    Expr* dcall = expr_new_function(expr_new_symbol("D"), da, 2);
    Expr* dexpr = eval_and_free(dcall);
    bool ok = ref_eval(dexpr, Y, d, t, out);
    expr_free(dexpr);
    return ok;
}

/* Compile a square system of `d` components (fstr[i] parsed) and compare
 * compiled eval + Jacobian to references over `trials` random points in
 * [lo, hi]. */
static void check_system(const char* name, const char** fstr, size_t d,
                         double lo, double hi, int trials, double jtol) {
    Expr** f = malloc(d * sizeof(Expr*));
    for (size_t i = 0; i < d; i++) f[i] = eval_and_free(parse_expression(fstr[i]));

    NdProblem P; memset(&P, 0, sizeof(P));
    P.d = d; P.tvar = intern_symbol("t"); P.f = f;
    P.ysym = malloc(d * sizeof(Expr*));
    for (size_t k = 0; k < d; k++) P.ysym[k] = wsym(k);

    NdCompiled* C = nd_compile_rhs(&P);
    if (!C) {
        printf("FAIL: %s -> nd_compile_rhs returned NULL\n", name); failures++;
        goto cleanup;
    }

    double *Y = malloc(d * sizeof(double)), *fc = malloc(d * sizeof(double));
    double *J = malloc(d * d * sizeof(double));
    int eval_ok = 0, eval_cmp = 0, jac_cmp = 0;
    double max_eval_err = 0.0, max_jac_err = 0.0;
    for (int tr = 0; tr < trials; tr++) {
        for (size_t k = 0; k < d; k++) Y[k] = urand(lo, hi);
        double t = urand(lo, hi);
        if (!nd_compiled_eval(C, t, Y, fc)) continue;     /* domain: skip */
        eval_ok++;
        /* compare each component to the symbolic reference */
        for (size_t i = 0; i < d; i++) {
            double ref;
            if (!ref_eval(f[i], Y, d, t, &ref)) continue;
            eval_cmp++;
            double err = fabs(fc[i] - ref) / (1.0 + fabs(ref));
            if (err > max_eval_err) max_eval_err = err;
        }
        /* Jacobian: colored FD vs symbolic D at a subset of points */
        if (tr % 8 == 0 && nd_compiled_jacobian(C, t, Y, J)) {
            for (size_t i = 0; i < d; i++)
                for (size_t j = 0; j < d; j++) {
                    double ref;
                    if (!ref_jac(f[i], j, Y, d, t, &ref)) continue;
                    jac_cmp++;
                    double err = fabs(J[i*d + j] - ref) / (1.0 + fabs(ref));
                    if (err > max_jac_err) max_jac_err = err;
                }
        }
    }
    if (eval_cmp == 0) { printf("FAIL: %s -> no comparable eval points\n", name); failures++; }
    else if (max_eval_err > 1e-9) {
        printf("FAIL: %s -> eval mismatch max_rel=%.2e (%d cmps)\n", name, max_eval_err, eval_cmp);
        failures++;
    } else if (max_jac_err > jtol) {
        printf("FAIL: %s -> jacobian mismatch max_rel=%.2e (%d cmps)\n", name, max_jac_err, jac_cmp);
        failures++;
    } else {
        printf("ok:   %-34s eval<=%.1e (%d) jac<=%.1e (%d) colors=%d/%zu\n",
               name, max_eval_err, eval_cmp, max_jac_err, jac_cmp, nd_compiled_ncolor(C), d);
    }
    free(Y); free(fc); free(J);
cleanup:
    nd_compiled_free(C);
    for (size_t i = 0; i < d; i++) { expr_free(f[i]); expr_free(P.ysym[i]); }
    free(f); free(P.ysym);
}

/* A single-component system that must fail to compile (unsupported construct)
 * -> NULL, so the solver keeps the symbolic sampler. `expr` references only the
 * state symbol NDSolve`w0 plus the offending construct. */
static void check_bail(const char* name, const char* expr) {
    Expr* f0 = eval_and_free(parse_expression(expr));
    NdProblem P; memset(&P, 0, sizeof(P));
    P.d = 1; P.tvar = intern_symbol("t");
    P.f = &f0;
    Expr* y0 = wsym(0);
    P.ysym = &y0;
    NdCompiled* C = nd_compile_rhs(&P);
    if (C) { printf("FAIL: %s -> compiled but should bail\n", name); failures++; nd_compiled_free(C); }
    else printf("ok:   %-34s bailed to symbolic fallback\n", name);
    expr_free(f0); expr_free(y0);
}

/* Regression guard: the compiler must recognize state symbols by matching the
 * problem's actual reduced-state names, NOT a hardcoded convention.  The ODE
 * front-end names its state NDSolve`y<k> (the PDE front-end uses NDSolve`w<k>);
 * both must compile. */
static void check_naming(void) {
    Expr** f = malloc(2 * sizeof(Expr*));
    f[0] = eval_and_free(parse_expression("NDSolve`y0 NDSolve`y1 + Sin[t]"));
    f[1] = eval_and_free(parse_expression("-NDSolve`y0 + NDSolve`y1^2"));
    NdProblem P; memset(&P, 0, sizeof(P));
    P.d = 2; P.tvar = intern_symbol("t"); P.f = f;
    P.ysym = malloc(2 * sizeof(Expr*));
    P.ysym[0] = expr_new_symbol(intern_symbol("NDSolve`y0"));
    P.ysym[1] = expr_new_symbol(intern_symbol("NDSolve`y1"));
    NdCompiled* C = nd_compile_rhs(&P);
    double Y[2] = { 1.5, 2.0 }, out[2] = { 0, 0 };
    bool ok = C && nd_compiled_eval(C, 0.5, Y, out);
    double e0 = 1.5 * 2.0 + sin(0.5), e1 = -1.5 + 4.0;
    if (ok && fabs(out[0] - e0) < 1e-12 && fabs(out[1] - e1) < 1e-12)
        printf("ok:   %-34s ODE-named state compiles+evals\n", "state naming (NDSolve`y<k>)");
    else { printf("FAIL: ODE-named state (NDSolve`y<k>) not recognized by compiler\n"); failures++; }
    nd_compiled_free(C);
    for (size_t i = 0; i < 2; i++) { expr_free(f[i]); expr_free(P.ysym[i]); }
    free(f); free(P.ysym);
}

int main(void) {
    core_init();
    check_naming();

    /* ---- broad arithmetic/elementary battery (positive domain) ---- */
    const char* battery[] = {
        "3 NDSolve`w0 + 2 NDSolve`w1 - NDSolve`w2/2 + 7",
        "NDSolve`w0 NDSolve`w1 NDSolve`w2",
        "NDSolve`w0^3 - 2 NDSolve`w1^2 + NDSolve`w2",
        "NDSolve`w0^2/NDSolve`w1 + NDSolve`w2^4",
        "Sqrt[NDSolve`w0] + Sqrt[NDSolve`w1 NDSolve`w2]",
        "Exp[-NDSolve`w0] + Log[NDSolve`w1] + Log[2, NDSolve`w2]",
        "Sin[NDSolve`w0] Cos[NDSolve`w1] + Tan[NDSolve`w2/3]",
        "Sinh[NDSolve`w0] - Cosh[NDSolve`w1] + Tanh[NDSolve`w2]",
        "ArcTan[NDSolve`w0] + ArcTan[NDSolve`w1, NDSolve`w2]",
        "Abs[NDSolve`w0 - NDSolve`w1] + Sign[NDSolve`w2 - 1]",
        "Max[NDSolve`w0, NDSolve`w1, NDSolve`w2] - Min[NDSolve`w0, NDSolve`w1]",
        "NDSolve`w0 Exp[NDSolve`w1] Sin[t] + t^2 NDSolve`w2",
        "(1 + NDSolve`w0^2)^(3/2) + NDSolve`w1^(1/2)",
        "1/(1 + NDSolve`w0^2 + NDSolve`w1^2) + NDSolve`w2^(-1)",
        "Erf[NDSolve`w0] + Erfc[NDSolve`w1] + Sqrt[Pi] NDSolve`w2",
        "Pi NDSolve`w0 + E NDSolve`w1 - NDSolve`w2",
    };
    check_system("arithmetic/elementary battery", battery,
                 sizeof battery / sizeof *battery, 0.25, 1.75, 400, 1e-5);

    /* ---- nonlinear-PDE-like couplings (shallow-water / Burgers / reaction) --- */
    const char* pde[] = {
        /* Burgers-like: -u u_x  with a diffusion term */
        "-NDSolve`w0 (NDSolve`w1 - NDSolve`w2) + (NDSolve`w1 - 2 NDSolve`w0 + NDSolve`w2)",
        /* shallow-water momentum-like: -u u_x - g h_x */
        "-NDSolve`w1 (NDSolve`w1 - NDSolve`w0) - (981/100)(NDSolve`w2 - NDSolve`w0)",
        /* mass-like: -(h u)_x */
        "-(NDSolve`w2 NDSolve`w1 - NDSolve`w0 NDSolve`w2)",
        /* cubic reaction */
        "NDSolve`w0 (1 - NDSolve`w0^2) + NDSolve`w1 NDSolve`w2",
    };
    check_system("nonlinear PDE couplings", pde, 4, 0.4, 1.6, 400, 1e-4);

    /* ---- tridiagonal (banded) system: coloring must be O(1) colors ---- */
    {
        const size_t D = 24;
        char** fstr = malloc(D * sizeof(char*));
        for (size_t i = 0; i < D; i++) {
            fstr[i] = malloc(160);
            size_t lo = i ? i - 1 : 0, hi = i + 1 < D ? i + 1 : D - 1;
            snprintf(fstr[i], 160,
                     "NDSolve`w%zu - 2 NDSolve`w%zu + NDSolve`w%zu + Sin[NDSolve`w%zu]",
                     lo, i, hi, i);
        }
        check_system("tridiagonal banded system", (const char**)fstr, D, 0.3, 1.5, 200, 1e-4);
        for (size_t i = 0; i < D; i++) free(fstr[i]);
        free(fstr);
    }

    /* ---- graceful bail on unsupported constructs ---- */
    /* Functions WITH a machine kernel (Gamma, BesselJ, ...) now compile through
     * the shared engine's generic-kernel path; genuine bails are functions with
     * no ndkernels entry (Zeta, PolyLog) or a free symbol. */
    check_bail("no kernel (Zeta)",       "NDSolve`w0 + Zeta[NDSolve`w0]");
    check_bail("no kernel (PolyLog)",    "NDSolve`w0 + PolyLog[2, NDSolve`w0]");
    check_bail("free parameter symbol",  "NDSolve`w0 + freeParameter");

    if (failures == 0) printf("\nAll NDSolve compile tests passed.\n");
    else printf("\n%d NDSolve compile test(s) FAILED.\n", failures);
    return failures ? 1 : 0;
}
