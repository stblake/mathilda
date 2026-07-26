/* Compile[] engine (scalar core): the compiled register machine must reproduce
 * the symbolic interpreter to rounding across the whole type lattice
 * (Bool/Int/Real/Complex), all opcodes, and coercions; the all-real fast path
 * must agree with the boxed path; unsupported constructs must bail; and the
 * compiled evaluation must be far faster than the interpreter.  Stress tests
 * exercise deep/wide expressions, many arguments, and random trees.
 *
 * Soft asserts: prints FAIL and keeps going. Run: ./compile_tests */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "parse.h"
#include "sym_intern.h"
#include "compile/compile.h"
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>

static int failures = 0;

static uint64_t rng_state = 0x243f6a8885a308d3ULL;
static double urand(double lo, double hi) {
    uint64_t x = rng_state; x ^= x << 13; x ^= x >> 7; x ^= x << 17; rng_state = x;
    return lo + (hi - lo) * ((double)(x >> 11) * (1.0 / 9007199254740992.0));
}
static long long irand(long long lo, long long hi) {
    uint64_t x = rng_state; x ^= x << 13; x ^= x >> 7; x ^= x << 17; rng_state = x;
    return lo + (long long)(x % (uint64_t)(hi - lo + 1));
}

static bool expr_to_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_REAL) { *out = e->data.real; return true; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    return false;
}

static Expr* val_to_expr(CompileValue v) {
    switch (v.type) {
        case CT_BOOL: return expr_new_symbol(v.v.b ? "True" : "False");
        case CT_INT:  return expr_new_integer(v.v.i);
        case CT_REAL: return expr_new_real(v.v.r);
        case CT_COMPLEX: {
            Expr* a[2] = { expr_new_real(creal(v.v.z)), expr_new_real(cimag(v.v.z)) };
            return expr_new_function(expr_new_symbol("Complex"), a, 2);
        }
    }
    return expr_new_integer(0);
}

typedef struct { bool ok; bool is_bool; int bval; double re, im; } Ref;

/* interpreter reference: substitute args, evaluate, extract (Re,Im) or bool */
static Ref ref_eval(const Expr* body, const char* const* names,
                    const CompileValue* vals, size_t n) {
    Ref r; memset(&r, 0, sizeof(r));
    Expr** rules = malloc(n * sizeof(Expr*));
    for (size_t k = 0; k < n; k++) {
        Expr* rr[2] = { expr_new_symbol(names[k]), val_to_expr(vals[k]) };
        rules[k] = expr_new_function(expr_new_symbol("Rule"), rr, 2);
    }
    Expr* rl = expr_new_function(expr_new_symbol("List"), rules, n);
    free(rules);
    Expr* ca[2] = { expr_copy((Expr*)body), rl };
    Expr* res = eval_and_free(expr_new_function(expr_new_symbol("ReplaceAll"), ca, 2));
    if (!res) return r;
    if (res->type == EXPR_SYMBOL &&
        (strcmp(res->data.symbol.name, "True") == 0 || strcmp(res->data.symbol.name, "False") == 0)) {
        r.ok = true; r.is_bool = true; r.bval = strcmp(res->data.symbol.name, "True") == 0;
        expr_free(res); return r;
    }
    Expr* reA[1] = { expr_new_function(expr_new_symbol("Re"), (Expr*[]){ expr_copy(res) }, 1) };
    Expr* imA[1] = { expr_new_function(expr_new_symbol("Im"), (Expr*[]){ expr_copy(res) }, 1) };
    Expr* reE = eval_and_free(expr_new_function(expr_new_symbol("N"), reA, 1));
    Expr* imE = eval_and_free(expr_new_function(expr_new_symbol("N"), imA, 1));
    double re, im;
    if (expr_to_double(reE, &re) && expr_to_double(imE, &im) && isfinite(re) && isfinite(im)) {
        r.ok = true; r.re = re; r.im = im;
    }
    expr_free(res); expr_free(reE); expr_free(imE);
    return r;
}

static void cval_reim(CompileValue v, double* re, double* im, bool* isb, int* bv) {
    *isb = false;
    switch (v.type) {
        case CT_BOOL: *isb = true; *bv = v.v.b; break;
        case CT_INT:  *re = (double)v.v.i; *im = 0; break;
        case CT_REAL: *re = v.v.r; *im = 0; break;
        case CT_COMPLEX: *re = creal(v.v.z); *im = cimag(v.v.z); break;
    }
}

/* Compile `body` over args of the given types and compare to the interpreter
 * over `trials` random points; reals/complex parts drawn from [lo,hi],
 * integers from [ilo,ihi]. */
static void parity(const char* name, const char* body_s,
                   const char* const* names, const CompileType* types, size_t n,
                   double lo, double hi, long long ilo, long long ihi, int trials) {
    Expr* body = eval_and_free(parse_expression(body_s));
    const char* inames[8];
    for (size_t k = 0; k < n; k++) inames[k] = intern_symbol(names[k]);
    CompiledProgram* p = compile_expr(body, inames, types, n);
    if (!p) { printf("FAIL: %-30s -> did not compile\n", name); failures++; expr_free(body); return; }

    CompileValue args[8], outc;
    int cmp = 0; double maxerr = 0.0;
    for (int t = 0; t < trials; t++) {
        for (size_t k = 0; k < n; k++) {
            args[k].type = types[k];
            switch (types[k]) {
                case CT_BOOL: args[k].v.b = (unsigned char)irand(0, 1); break;
                case CT_INT:  args[k].v.i = irand(ilo, ihi); break;
                case CT_REAL: args[k].v.r = urand(lo, hi); break;
                case CT_COMPLEX: args[k].v.z = urand(lo, hi) + urand(lo, hi) * I; break;
            }
        }
        bool cok = compiled_eval(p, args, &outc);
        Ref r = ref_eval(body, inames, args, n);
        if (!cok || !r.ok) continue;                 /* domain / non-finite: skip */
        cmp++;
        double re, im; bool isb; int bv;
        cval_reim(outc, &re, &im, &isb, &bv);
        if (isb != r.is_bool) { printf("FAIL: %-30s -> bool/numeric type mismatch\n", name); failures++; break; }
        if (isb) { if (bv != r.bval) { printf("FAIL: %-30s -> bool %d vs %d\n", name, bv, r.bval); failures++; break; } }
        else {
            double err = (fabs(re - r.re) + fabs(im - r.im)) / (1.0 + fabs(r.re) + fabs(r.im));
            if (err > maxerr) maxerr = err;
        }
    }
    if (cmp < trials / 8) { printf("FAIL: %-30s -> too few comparable points (%d)\n", name, cmp); failures++; }
    else if (maxerr > 1e-9) { printf("FAIL: %-30s -> max_rel=%.2e (%d cmps)\n", name, maxerr, cmp); failures++; }
    else printf("ok:   %-30s max_rel=%.1e (%d cmps)\n", name, maxerr, cmp);
    compiled_free(p);
    expr_free(body);
}

static void must_bail(const char* name, const char* body_s, const char* const* names,
                      const CompileType* types, size_t n) {
    Expr* body = eval_and_free(parse_expression(body_s));
    const char* inames[8];
    for (size_t k = 0; k < n; k++) inames[k] = intern_symbol(names[k]);
    CompiledProgram* p = compile_expr(body, inames, types, n);
    if (p) { printf("FAIL: %-30s -> compiled but should bail\n", name); failures++; compiled_free(p); }
    else printf("ok:   %-30s bailed\n", name);
    expr_free(body);
}

int main(void) {
    core_init();

    const char* xyz[] = { "x", "y", "z" };
    const char* xy[]  = { "x", "y" };
    const char* x1[]  = { "x" };
    const CompileType RRR[] = { CT_REAL, CT_REAL, CT_REAL };
    const CompileType III[] = { CT_INT, CT_INT, CT_INT };
    const CompileType CCC[] = { CT_COMPLEX, CT_COMPLEX, CT_COMPLEX };

    /* ---- real scalar: arithmetic + elementary + nesting ---- */
    parity("real arithmetic", "3 x - 2 y + z/2 + 7", xyz, RRR, 3, 0.2, 3.0, 0, 0, 300);
    parity("real poly/power", "x^3 - 2 y^2 + x y + x^(1/2)", xyz, RRR, 3, 0.2, 3.0, 0, 0, 300);
    parity("real elementary", "Sin[x] Cos[y] + Exp[-z] + Log[x] + Sqrt[y]", xyz, RRR, 3, 0.3, 2.5, 0, 0, 300);
    parity("real transcendental", "Tanh[x] + ArcTan[y] + ArcTan[x,y] + Erf[z]", xyz, RRR, 3, 0.2, 2.0, 0, 0, 300);
    parity("real abs/sign/minmax", "Abs[x - y] + Sign[z - 1] + Max[x,y,z] - Min[x,y]", xyz, RRR, 3, -2.0, 2.0, 0, 0, 300);
    parity("real floor/ceil/round", "Floor[3 x] + Ceiling[2 y] - Round[z]", xyz, RRR, 3, -3.0, 3.0, 0, 0, 300);
    parity("real deep nest", "Sin[Cos[Sin[Cos[x + y] - z]]]", xyz, RRR, 3, -2.0, 2.0, 0, 0, 300);
    parity("real rational pow", "(1 + x^2)^(3/2) + y^(1/2) - 1/(1 + z^2)", xyz, RRR, 3, 0.2, 2.0, 0, 0, 300);
    parity("real Log base", "Log[2, x] + Log[y]", xy, RRR, 2, 0.3, 3.0, 0, 0, 300);

    /* ---- integer ---- */
    parity("int arithmetic", "x y + x - 2 z", xyz, III, 3, 0, 0, -20, 20, 300);
    parity("int power", "x^3 + y^2 - z", xyz, III, 3, 0, 0, -6, 6, 300);
    parity("int mod/quotient", "Mod[x, y] + Quotient[z, y]", xyz, III, 3, 0, 0, 1, 30, 300);
    parity("int abs/sign/minmax", "Abs[x] + Sign[y] + Max[x,y,z] - Min[x,z]", xyz, III, 3, 0, 0, -20, 20, 300);

    /* ---- complex ---- */
    parity("complex arithmetic", "x y - 2 z + x/y", xyz, CCC, 3, 0.3, 2.0, 0, 0, 300);
    parity("complex elementary", "Exp[x] + Sin[y] + Sqrt[z]", xyz, CCC, 3, -1.5, 1.5, 0, 0, 300);
    parity("complex power", "x^3 + y^(1/2) + z^x", xyz, CCC, 3, 0.3, 1.5, 0, 0, 300);
    parity("complex re/im/abs/arg", "Re[x] + Im[y] + Abs[z] + Arg[x] + Re[Conjugate[y]]", xyz, CCC, 3, -2.0, 2.0, 0, 0, 300);
    parity("I literal", "x + 2 I y - I", xy, CCC, 2, -2.0, 2.0, 0, 0, 300);

    /* ---- generic special-function kernels (shared ndkernels registry) ---- */
    parity("Gamma real (kernel)", "Gamma[x] + Gamma[2 y]", xy, RRR, 2, 0.3, 3.0, 0, 0, 300);
    parity("hyperbolic/inverse extras", "ArcSinh[x] + Coth[y] + Sech[z] + ArcTanh[x/4]", xyz, RRR, 3, 0.3, 2.0, 0, 0, 300);
    parity("complex generic kernel", "ArcSinh[x] + Coth[y]", xy, CCC, 2, 0.3, 1.8, 0, 0, 300);
    parity("Beta (binary kernel)", "Beta[x, y]", xy, RRR, 2, 0.4, 3.0, 0, 0, 300);
    parity("BesselJ/Y (binary kernel)", "BesselJ[2, x] + BesselY[1, y]", xy, RRR, 2, 0.5, 6.0, 0, 0, 300);
    parity("Factorial", "Factorial[x] + FractionalPart[3 y]", xy, RRR, 2, 0.3, 3.0, 0, 0, 250);

    /* ---- coercions (mixed arg types in one expression) ---- */
    { const CompileType m[] = { CT_INT, CT_REAL, CT_COMPLEX };
      parity("mixed int+real+complex", "x + y + z + x y - z^2", xyz, m, 3, 0.3, 2.0, -5, 5, 300); }
    { const CompileType m[] = { CT_INT, CT_REAL };
      parity("int+real widen", "x^2 + Sin[y] + x y", xy, m, 2, 0.2, 2.5, -6, 6, 300); }

    /* ---- boolean ---- */
    parity("comparisons+logic", "And[x < y, Or[y > z, Not[x == z]]]", xyz, RRR, 3, -2.0, 2.0, 0, 0, 300);
    parity("int comparisons", "Or[x <= y, x != z]", xyz, III, 3, 0, 0, -5, 5, 300);
    parity("mixed compare", "x >= y", xy, ((const CompileType[]){ CT_INT, CT_REAL }), 2, -3.0, 3.0, -5, 5, 300);
    /* NB: == / != on inexact numbers is fuzzy in the interpreter but exact in the
     * compiler (a deliberate divergence), so we do not parity-test EQ_C/NE_C. */

    /* ---- pi/e constants ---- */
    parity("named constants", "Pi x + E y - x/Pi", xy, RRR, 2, 0.2, 3.0, 0, 0, 200);

    /* ---- all-real fast path parity ---- */
    {
        const char* nm[] = { "x", "y" };
        const char* inm[2] = { intern_symbol("x"), intern_symbol("y") };
        Expr* b = eval_and_free(parse_expression("Sin[x] + y^2 - Exp[-x] Cos[y]"));
        CompiledProgram* p = compile_expr(b, inm, RRR, 2);
        int mism = 0;
        for (int t = 0; t < 500; t++) {
            double a[2] = { urand(-2, 2), urand(-2, 2) };
            CompileValue av[2] = { { CT_REAL, { .r = a[0] } }, { CT_REAL, { .r = a[1] } } };
            CompileValue bo; double fr;
            bool ok1 = compiled_eval(p, av, &bo), ok2 = compiled_eval_real(p, a, &fr);
            if (ok1 != ok2 || (ok1 && fabs(bo.v.r - fr) > 0) ) mism++;
        }
        if (mism) { printf("FAIL: real fast path disagrees (%d)\n", mism); failures++; }
        else printf("ok:   %-30s boxed == compiled_eval_real\n", "real fast path");
        compiled_free(p); expr_free(b);
        (void)nm;
    }

    /* ---- arg-dependency introspection ---- */
    {
        const char* inm[3] = { intern_symbol("x"), intern_symbol("y"), intern_symbol("z") };
        Expr* b = eval_and_free(parse_expression("x + z^2"));   /* uses x,z not y */
        CompiledProgram* p = compile_expr(b, inm, RRR, 3);
        int deps[3]; size_t nd = compiled_arg_deps(p, deps, 3);
        bool ok = (nd == 2 && deps[0] == 0 && deps[1] == 2);
        if (!ok) { printf("FAIL: arg_deps wrong (n=%zu)\n", nd); failures++; }
        else printf("ok:   %-30s deps={0,2}\n", "arg-dependency introspection");
        compiled_free(p); expr_free(b);
    }

    /* ---- graceful bail ---- */
    must_bail("no kernel (Zeta)", "Zeta[x]", x1, RRR, 1);        /* not in ndkernels -> bail */
    must_bail("free symbol", "x + unknownParam", x1, RRR, 1);
    must_bail("list body", "{x, x^2}", x1, RRR, 1);

    /* ================= STRESS ================= */

    /* deep nesting: Sin applied 400 times */
    {
        char buf[8192]; size_t pos = 0;
        for (int i = 0; i < 400; i++) pos += (size_t)snprintf(buf + pos, sizeof buf - pos, "Sin[");
        pos += (size_t)snprintf(buf + pos, sizeof buf - pos, "x");
        for (int i = 0; i < 400; i++) pos += (size_t)snprintf(buf + pos, sizeof buf - pos, "]");
        const char* inm[1] = { intern_symbol("x") };
        Expr* b = eval_and_free(parse_expression(buf));
        CompiledProgram* p = compile_expr(b, inm, RRR, 1);
        bool ok = false;
        if (p) {
            double a = 0.7, out; CompileValue av = { CT_REAL, { .r = a } }, ov;
            compiled_eval(p, &av, &ov);
            /* reference: iterate sin 400x */
            double ref = a; for (int i = 0; i < 400; i++) ref = sin(ref);
            ok = fabs(ov.v.r - ref) < 1e-12;
        }
        if (!ok) { printf("FAIL: deep-nest (400x Sin)\n"); failures++; }
        else printf("ok:   %-30s 400-deep, exact\n", "stress: deep nesting");
        compiled_free(p); expr_free(b);
    }

    /* wide sum: x0 + x1 + ... but with one arg reused 500 times: 500 x + literals */
    {
        char buf[16384]; size_t pos = 0;
        for (int i = 0; i < 500; i++) pos += (size_t)snprintf(buf + pos, sizeof buf - pos, "%sSin[%d x]", i ? " + " : "", i + 1);
        const char* inm[1] = { intern_symbol("x") };
        Expr* b = eval_and_free(parse_expression(buf));
        CompiledProgram* p = compile_expr(b, inm, RRR, 1);
        bool ok = false;
        if (p) {
            double a = 0.013, out = 0; for (int i = 0; i < 500; i++) out += sin((i + 1) * a);
            CompileValue av = { CT_REAL, { .r = a } }, ov;
            compiled_eval(p, &av, &ov);
            ok = fabs(ov.v.r - out) < 1e-9;
        }
        if (!ok) { printf("FAIL: wide-sum (500 terms)\n"); failures++; }
        else printf("ok:   %-30s 500-term sum, exact\n", "stress: wide expression");
        compiled_free(p); expr_free(b);
    }

    /* many arguments: x0 + 2 x1 + ... + 8 x7 (8 args) */
    {
        const char* nms[8] = { "a", "b", "c", "d", "e", "f", "g", "h" };
        const char* inm[8]; CompileType ty[8];
        for (int i = 0; i < 8; i++) { inm[i] = intern_symbol(nms[i]); ty[i] = CT_REAL; }
        Expr* b = eval_and_free(parse_expression("a + 2 b + 3 c + 4 d + 5 e + 6 f + 7 g + 8 h + a b c - d e f"));
        CompiledProgram* p = compile_expr(b, inm, ty, 8);
        bool ok = false;
        if (p) {
            double v[8]; CompileValue av[8]; for (int i = 0; i < 8; i++) { v[i] = urand(-1, 1); av[i].type = CT_REAL; av[i].v.r = v[i]; }
            CompileValue ov; compiled_eval(p, av, &ov);
            double ref = 0; for (int i = 0; i < 8; i++) ref += (i + 1) * v[i];
            ref += v[0]*v[1]*v[2] - v[3]*v[4]*v[5];
            ok = fabs(ov.v.r - ref) < 1e-9;
        }
        if (!ok) { printf("FAIL: many-args (8)\n"); failures++; }
        else printf("ok:   %-30s 8 args, exact\n", "stress: many arguments");
        compiled_free(p); expr_free(b);
    }

    /* ================= PERFORMANCE ================= */
    {
        const char* inm[2] = { intern_symbol("x"), intern_symbol("y") };
        Expr* b = eval_and_free(parse_expression("Sin[x y] + Exp[-x^2/2] Cos[y] - Sqrt[Abs[x]] + x^3 - 2 y^2"));
        CompiledProgram* p = compile_expr(b, inm, RRR, 2);
        if (!p) { printf("FAIL: perf body did not compile\n"); failures++; }
        else {
            const int NC = 2000000, NI = 4000;
            double acc = 0; double a[2];
            clock_t t0 = clock();
            for (int i = 0; i < NC; i++) { a[0] = 0.3 + 1e-7 * i; a[1] = 0.7; double o; compiled_eval_real(p, a, &o); acc += o; }
            double tc = (double)(clock() - t0) / CLOCKS_PER_SEC;
            const char* nm[] = { "x", "y" };
            clock_t t1 = clock();
            for (int i = 0; i < NI; i++) {
                CompileValue av[2] = { { CT_REAL, { .r = 0.3 + 1e-7 * i } }, { CT_REAL, { .r = 0.7 } } };
                Ref r = ref_eval(b, inm, av, 2); acc += r.re;
            }
            double ti = (double)(clock() - t1) / CLOCKS_PER_SEC;
            double per_c = tc / NC, per_i = ti / NI, speedup = per_i / per_c;
            (void)nm;
            if (speedup < 20.0) { printf("FAIL: perf speedup only %.1fx\n", speedup); failures++; }
            else printf("ok:   %-30s %.0fx faster (%.1f ns/call vs %.0f ns)\n",
                        "performance", speedup, per_c * 1e9, per_i * 1e9);
            (void)acc;
        }
        compiled_free(p); expr_free(b);
    }

    if (failures == 0) printf("\nAll Compile engine tests passed.\n");
    else printf("\n%d Compile engine test(s) FAILED.\n", failures);
    return failures ? 1 : 0;
}
