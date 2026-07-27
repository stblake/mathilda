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
#include "match.h"
#include "ndarray.h"
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

/* ------------------------------------------------------------------ *
 *  Arrays (M3a): rank-1 machine vectors                               *
 * ------------------------------------------------------------------ */

/* A rank-1 float64 NDArray of `n` entries drawn from [lo,hi]. */
static Expr* make_vec(size_t n, double lo, double hi) {
    double* buf = malloc(n * sizeof(double));
    for (size_t k = 0; k < n; k++) buf[k] = urand(lo, hi);
    int64_t dims[1]; dims[0] = (int64_t)n;
    return expr_new_ndarray(1, dims, buf, NDT_FLOAT64);   /* adopts buf */
}

/* A rank-1 complex64 NDArray (interleaved re,im) with entries from [lo,hi]^2. */
static Expr* make_cvec(size_t n, double lo, double hi) {
    double* buf = malloc(2 * n * sizeof(double));
    for (size_t k = 0; k < n; k++) { buf[2 * k] = urand(lo, hi); buf[2 * k + 1] = urand(lo, hi); }
    int64_t dims[1]; dims[0] = (int64_t)n;
    return expr_new_ndarray(1, dims, buf, NDT_COMPLEX64);   /* adopts buf */
}

/* Interpreter reference: substitute the NDArray arguments and evaluate, so the
 * comparison is against the very same NDArray fast paths the VM delegates to.
 * Substitution goes through the binder rather than ReplaceAll, because
 * ReplaceAll evaluates its first argument first — which would collapse the
 * unevaluated body (Length[v] -> 0) before v was ever bound. */
static Expr* ref_eval_arr(const Expr* body, const char* const* names,
                          Expr* const* vals, size_t n) {
    MatchEnv* env = env_new();
    for (size_t k = 0; k < n; k++) env_set(env, names[k], vals[k]);
    Expr* sub = replace_bindings((Expr*)body, env);
    env_free(env);
    return eval_and_free(sub);
}

/* (Re, Im) of a numeric scalar Expr — Real, Integer, or Complex[re, im]. */
static bool scalar_reim(const Expr* e, double* re, double* im) {
    if (expr_to_double(e, re)) { *im = 0.0; return true; }
    Expr* reA[1] = { expr_new_function(expr_new_symbol("Re"), (Expr*[]){ expr_copy((Expr*)e) }, 1) };
    Expr* imA[1] = { expr_new_function(expr_new_symbol("Im"), (Expr*[]){ expr_copy((Expr*)e) }, 1) };
    Expr* reE = eval_and_free(expr_new_function(expr_new_symbol("N"), reA, 1));
    Expr* imE = eval_and_free(expr_new_function(expr_new_symbol("N"), imA, 1));
    bool ok = expr_to_double(reE, re) && expr_to_double(imE, im);
    expr_free(reE); expr_free(imE);
    return ok;
}

/* Accumulate the max relative difference between two results (array or scalar).
 * Returns false when the two do not even agree on shape/kind. */
static bool arr_cmp(const Expr* got, const Expr* want, double* maxerr) {
    if (!got || !want) return false;
    if (got->type == EXPR_NDARRAY || want->type == EXPR_NDARRAY) {
        if (got->type != EXPR_NDARRAY || want->type != EXPR_NDARRAY) return false;
        size_t n = ndarray_size(got);
        if (n != ndarray_size(want) || got->data.ndarray.rank != want->data.ndarray.rank) return false;
        for (size_t k = 0; k < n; k++) {
            double ar, ai, br, bi;
            ndt_get(got->data.ndarray.data, k, got->data.ndarray.dtype, &ar, &ai);
            ndt_get(want->data.ndarray.data, k, want->data.ndarray.dtype, &br, &bi);
            double e = (fabs(ar - br) + fabs(ai - bi)) / (1.0 + fabs(br) + fabs(bi));
            if (e > *maxerr) *maxerr = e;
        }
        return true;
    }
    double gr, gi, wr, wi;
    if (!scalar_reim(got, &gr, &gi) || !scalar_reim(want, &wr, &wi)) return false;
    double e = (fabs(gr - wr) + fabs(gi - wi)) / (1.0 + fabs(wr) + fabs(wi));
    if (e > *maxerr) *maxerr = e;
    return true;
}

/* Box a compiled array/scalar result as an Expr the comparison can read.
 * Takes ownership of an array result. */
static Expr* aval_to_expr(CompileValue v) {
    if (CT_IS_ARRAY(v.type)) return v.v.a;
    return val_to_expr(v);
}

/* Compile `body` over array-typed args and compare to the interpreter over
 * `trials` freshly built vectors.  The body is parsed but NOT evaluated: with
 * the array parameters still free symbols the evaluator would rewrite the very
 * constructs under test (Total[v] collapses to v, With[{u=v},...] inlines), so
 * only the raw parse tree exercises the array lowering. */
static void parity_arr(const char* name, const char* body_s,
                       const char* const* names, const CompileType* types, size_t n,
                       size_t len, double lo, double hi, int trials) {
    Expr* body = parse_expression(body_s);
    const char* inames[4];
    for (size_t k = 0; k < n; k++) inames[k] = intern_symbol(names[k]);
    CompiledProgram* p = compile_expr(body, inames, types, n);
    if (!p) { printf("FAIL: %-30s -> did not compile\n", name); failures++; expr_free(body); return; }

    int cmp = 0; double maxerr = 0.0; bool shape_ok = true;
    for (int t = 0; t < trials; t++) {
        Expr* vecs[4]; CompileValue args[4], outc;
        for (size_t k = 0; k < n; k++) {
            vecs[k] = CT_ELEM(types[k]) == CT_COMPLEX ? make_cvec(len, lo, hi)
                                                      : make_vec(len, lo, hi);
            args[k].type = types[k]; args[k].v.a = vecs[k];
        }
        bool cok = compiled_eval(p, args, &outc);
        if (cok) {
            Expr* want = ref_eval_arr(body, inames, vecs, n);
            Expr* got  = aval_to_expr(outc);
            if (!arr_cmp(got, want, &maxerr)) shape_ok = false;
            else cmp++;
            expr_free(got); expr_free(want);
        }
        for (size_t k = 0; k < n; k++) expr_free(vecs[k]);   /* args are borrowed */
        if (!shape_ok) break;
    }
    if (!shape_ok) { printf("FAIL: %-30s -> result shape/kind mismatch\n", name); failures++; }
    else if (cmp < trials) { printf("FAIL: %-30s -> only %d/%d evaluated\n", name, cmp, trials); failures++; }
    else if (maxerr > 1e-12) { printf("FAIL: %-30s -> max_rel=%.2e\n", name, maxerr); failures++; }
    else printf("ok:   %-30s max_rel=%.1e (%d cmps)\n", name, maxerr, cmp);
    compiled_free(p);
    expr_free(body);
}

static void bail_body(const char* name, Expr* body, const char* const* names,
                      const CompileType* types, size_t n) {
    const char* inames[8];
    for (size_t k = 0; k < n; k++) inames[k] = intern_symbol(names[k]);
    CompiledProgram* p = compile_expr(body, inames, types, n);
    if (p) { printf("FAIL: %-30s -> compiled but should bail\n", name); failures++; compiled_free(p); }
    else printf("ok:   %-30s bailed\n", name);
    expr_free(body);
}
static void must_bail(const char* name, const char* body_s, const char* const* names,
                      const CompileType* types, size_t n) {
    bail_body(name, eval_and_free(parse_expression(body_s)), names, types, n);
}
/* Same, on the raw parse tree — see parity_arr on why array bodies must not be
 * pre-evaluated. */
static void must_bail_raw(const char* name, const char* body_s, const char* const* names,
                          const CompileType* types, size_t n) {
    bail_body(name, parse_expression(body_s), names, types, n);
}

/* ---- optimiser A/B ------------------------------------------------------
 * Compile the same body twice, with and without the optimiser, and require the
 * two programs to agree BITWISE over a randomised argument sweep.  Bitwise (not
 * "to rounding") is the right gate: a pass that reassociated a sum or contracted
 * a multiply-add would still look accurate but would break the engine's stated
 * parity contract with the interpreter. */
static long long* ab_opt_tot = NULL;
static long long* ab_raw_tot = NULL;
static int*       ab_count   = NULL;
static void ab_opt_init(long long* o, long long* r, int* n) {
    ab_opt_tot = o; ab_raw_tot = r; ab_count = n;
}

static void ab_opt(const char* name, const char* body_s, const char* const* names,
                   const CompileType* types, size_t n, bool raw_parse) {
    const char* inames[4];
    for (size_t k = 0; k < n; k++) inames[k] = intern_symbol(names[k]);
    Expr* body = raw_parse ? parse_expression(body_s)
                           : eval_and_free(parse_expression(body_s));
    CompiledProgram* po = compile_expr_ex(body, inames, types, n, 0u);
    CompiledProgram* pr = compile_expr_ex(body, inames, types, n, COMPILE_NO_OPT);

    if (!po != !pr) {
        printf("FAIL: %-30s optimiser changed whether the body compiles\n", name);
        failures++;
        compiled_free(po); compiled_free(pr); expr_free(body);
        return;
    }
    if (!po) {   /* both bailed: nothing to compare, and that is consistent */
        expr_free(body);
        return;
    }

    int bad = 0;
    for (int trial = 0; trial < 60 && !bad; trial++) {
        CompileValue av[4];
        for (size_t k = 0; k < n; k++) {
            av[k].type = types[k];
            switch (types[k]) {
                case CT_INT:     av[k].v.i = irand(1, 12); break;
                case CT_COMPLEX: av[k].v.z = urand(0.3, 2.5) + urand(-2.0, 2.0) * I; break;
                default:         av[k].v.r = urand(0.35, 3.0); break;
            }
        }
        CompileValue oo, rr;
        bool so = compiled_eval(po, av, &oo);
        bool sr = compiled_eval(pr, av, &rr);
        if (so != sr) { bad = 1; break; }
        if (!so) continue;                       /* both declined, consistently */
        if (oo.type != rr.type) { bad = 1; break; }
        /* memcmp, not ==: NaN must equal NaN and -0.0 must differ from +0.0,
         * because either would mean a pass changed the arithmetic. */
        switch (oo.type) {
            case CT_BOOL:    if (oo.v.b != rr.v.b) bad = 1; break;
            case CT_INT:     if (oo.v.i != rr.v.i) bad = 1; break;
            case CT_REAL:    if (memcmp(&oo.v.r, &rr.v.r, sizeof(double)) != 0) bad = 1; break;
            case CT_COMPLEX: if (memcmp(&oo.v.z, &rr.v.z, sizeof(double _Complex)) != 0) bad = 1; break;
            default: break;
        }
    }

    size_t no = compiled_num_instructions(po), nr = compiled_num_instructions(pr);
    if (bad) {
        printf("FAIL: %-30s optimised result differs from unoptimised\n", name);
        failures++;
    } else if (no > nr) {
        printf("FAIL: %-30s optimiser GREW the program (%zu -> %zu)\n", name, nr, no);
        failures++;
    } else {
        if (ab_opt_tot) { *ab_opt_tot += (long long)no; *ab_raw_tot += (long long)nr; (*ab_count)++; }
    }
    compiled_free(po); compiled_free(pr); expr_free(body);
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

    /* ---- control flow: If (branches) ---- */
    parity("If real branches", "If[x < y, x^2 + 1, Sin[y] - z]", xyz, RRR, 3, -2.0, 2.0, 0, 0, 300);
    parity("If nested", "If[x < 0, If[y < 0, x + y, x - y], x y]", xyz, RRR, 3, -2.0, 2.0, 0, 0, 300);
    parity("If int branches", "If[x > y, x + z, x - z]", xyz, III, 3, 0, 0, -10, 10, 300);
    parity("If widen branches", "If[x < y, x, 5 y/2]", xy, ((const CompileType[]){ CT_INT, CT_REAL }), 2, -3.0, 3.0, -5, 5, 300);
    parity("If complex branches", "If[Re[x] < 1, x^2, Exp[y]]", xy, CCC, 2, -1.5, 1.5, 0, 0, 300);
    parity("If bool branches", "If[x < y, y < z, x < z]", xyz, RRR, 3, -2.0, 2.0, 0, 0, 300);
    parity("If in arithmetic", "3 If[x < y, x, y] + Abs[z]", xyz, RRR, 3, -2.0, 2.0, 0, 0, 300);

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

    /* ---- arrays (M3a): rank-1 machine vectors ----
     * The compiled array ops delegate to the same NDArray fast paths the
     * interpreter uses, so parity here is exact, not just to rounding. */
    {
        const char* vw[] = { "v", "w" };
        const CompileType AA[] = { CT_ARRAY(CT_REAL, 1), CT_ARRAY(CT_REAL, 1) };
        const int N = 40;

        /* array -> array */
        parity_arr("vec + scalar",      "v + 1",        vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("vec * scalar",      "2 v",          vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("vec + vec",         "v + w",        vw, AA, 2, 64, 0.3, 4.0, N);
        parity_arr("vec * vec",         "v w",          vw, AA, 2, 64, 0.3, 4.0, N);
        parity_arr("vec - vec",         "v - w",        vw, AA, 2, 64, 0.3, 4.0, N);
        parity_arr("vec - scalar",      "v - 2",        vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("scalar - vec",      "3 - v",        vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("vec negate",        "-v",           vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("vec / scalar",      "v / 3",        vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("vec / vec",         "v / w",        vw, AA, 2, 64, 0.3, 4.0, N);
        parity_arr("vec ^ integer",     "v^3",          vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("vec ^ vec",         "v^w",          vw, AA, 2, 64, 0.3, 2.0, N);
        parity_arr("scalar ^ vec",      "2^v",          vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("Sqrt[vec]",         "Sqrt[v]",      vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("Sin[vec]",          "Sin[v]",       vw, AA, 1, 64, -3.0, 3.0, N);
        parity_arr("Exp/Log[vec]",      "Exp[-v] + Log[v]", vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("Abs[vec]",          "Abs[v]",       vw, AA, 1, 64, -4.0, 4.0, N);
        parity_arr("Gamma[vec]",        "Gamma[v]",     vw, AA, 1, 64, 0.5, 4.0, N);
        parity_arr("Log[b, vec]",       "Log[2, v]",    vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("ArcTan[vec, y]",    "ArcTan[v, 2]", vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("nested vec expr",   "Sin[v w] + Exp[-v] (w + 2)", vw, AA, 2, 64, 0.3, 3.0, N);

        /* array -> scalar */
        parity_arr("Total[vec]",        "Total[v]",     vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("Total[Sin[vec]]",   "Total[Sin[v]]", vw, AA, 1, 64, -3.0, 3.0, N);
        parity_arr("Length[vec]",       "Length[v]",    vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("mean",              "Total[v] / Length[v]", vw, AA, 1, 64, 0.3, 4.0, N);
        parity_arr("dot product",       "Total[v w]",   vw, AA, 2, 64, 0.3, 4.0, N);
        parity_arr("scalar in array expr", "Total[(v - Total[v]/Length[v])^2] / Length[v]",
                   vw, AA, 1, 64, 0.3, 4.0, N);

        /* complex-element vectors: the complex64 buffer and the complex scalar
         * broadcast, both of which read the other half of the register slot. */
        {
            const CompileType CA[] = { CT_ARRAY(CT_COMPLEX, 1), CT_ARRAY(CT_COMPLEX, 1) };
            parity_arr("cvec + cvec",   "v + w",      vw, CA, 2, 32, 0.3, 3.0, N);
            parity_arr("cvec * cvec",   "v w",        vw, CA, 2, 32, 0.3, 3.0, N);
            parity_arr("cvec + I",      "v + 2 I",    vw, CA, 1, 32, 0.3, 3.0, N);
            parity_arr("Exp[cvec]",     "Exp[v]",     vw, CA, 1, 32, 0.3, 2.0, N);
            parity_arr("Abs[cvec]",     "Abs[v]",     vw, CA, 1, 32, 0.3, 3.0, N);
            parity_arr("Total[cvec]",   "Total[v]",   vw, CA, 1, 32, 0.3, 3.0, N);
            parity_arr("cvec ^ 2",      "v^2",        vw, CA, 1, 32, 0.3, 3.0, N);
        }

        /* Constructs that would need to COPY an array handle rather than move
         * it — deferred to M3b, and each must bail rather than alias. */
        must_bail_raw("array identity",     "v",                 vw, AA, 1);
        must_bail_raw("array If branch",    "If[1 < 2, v, v+1]", vw, AA, 1);
        must_bail_raw("array With local",   "With[{u = v}, u + 1]", vw, AA, 1);
        must_bail_raw("array Sum body",     "Sum[v, {i, 1, 3}]", vw, AA, 1);
        must_bail_raw("array comparison",   "v < w",             vw, AA, 2);
        must_bail_raw("array Max",          "Max[v, w]",         vw, AA, 2);
        must_bail_raw("array Mod",          "Mod[v, w]",         vw, AA, 2);

        /* Runtime contract: a real-typed program that would have to leave the
         * real axis fails the call (caller falls back), it does not lie. */
        {
            const char* inm[1] = { intern_symbol("v") };
            Expr* b = parse_expression("Sqrt[v]");
            CompiledProgram* p = compile_expr(b, inm, AA, 1);
            double neg[3] = { 1.0, -4.0, 9.0 }, pos[3] = { 1.0, 4.0, 9.0 };
            int bad = 0;
            for (int t = 0; t < 2; t++) {
                double* src = t ? pos : neg;
                double* buf = malloc(3 * sizeof(double));
                memcpy(buf, src, 3 * sizeof(double));
                int64_t dims[1] = { 3 };
                Expr* v = expr_new_ndarray(1, dims, buf, NDT_FLOAT64);
                CompileValue a, o; a.type = AA[0]; a.v.a = v;
                bool ok = p && compiled_eval(p, &a, &o);
                if (ok != (t == 1)) bad++;
                if (ok) expr_free(o.v.a);
                expr_free(v);
            }
            if (!p || bad) { printf("FAIL: %-30s -> complex-promotion not rejected\n", "Sqrt[vec] real contract"); failures++; }
            else printf("ok:   %-30s negative entry -> fallback\n", "Sqrt[vec] real contract");
            compiled_free(p); expr_free(b);
        }

        /* Array temporaries must be released every call (the frame is reused),
         * and released again when a later call aborts partway through. */
        {
            const char* inm[1] = { intern_symbol("v") };
            Expr* b = parse_expression("Total[Sin[v] Exp[-v] + Sqrt[v]]");
            CompiledProgram* p = compile_expr(b, inm, AA, 1);
            double acc = 0; int ok = 1;
            for (int t = 0; t < 5000 && p; t++) {
                Expr* v = make_vec(32, 0.3, 3.0);
                CompileValue a, o; a.type = AA[0]; a.v.a = v;
                if (!compiled_eval(p, &a, &o)) ok = 0; else acc += o.v.r;
                expr_free(v);
            }
            if (!p || !ok) { printf("FAIL: %-30s -> repeated array calls failed\n", "array temp lifetime"); failures++; }
            else printf("ok:   %-30s 5000 calls, acc=%.3f\n", "array temp lifetime", acc);
            compiled_free(p); expr_free(b);
        }

        /* Where the array path pays.  The buffer work is the SAME code the
         * interpreter runs (both call the ND kernels), so the win is purely the
         * per-operation evaluator round-trip — Expr build, attribute lookup,
         * dispatch — that the compiled program does not pay.  It is therefore
         * largest for short vectors, where the buffer pass is cheap relative to
         * that fixed cost, and shrinks as the vector grows. */
        for (int which = 0; which < 2; which++) {
            size_t len = which ? 4096 : 16;
            const char* inm[1] = { intern_symbol("v") };
            Expr* b = parse_expression("Total[Sin[v] Exp[-v] + Sqrt[v]]");
            CompiledProgram* p = compile_expr(b, inm, AA, 1);
            Expr* v = make_vec(len, 0.3, 3.0);
            CompileValue a, o; a.type = AA[0]; a.v.a = v;
            const int NC = which ? 300 : 20000, NI = which ? 300 : 20000;
            double acc = 0;
            clock_t t0 = clock();
            for (int t = 0; t < NC; t++) { if (compiled_eval(p, &a, &o)) acc += o.v.r; }
            double tc = (double)(clock() - t0) / CLOCKS_PER_SEC / NC;
            t0 = clock();
            for (int t = 0; t < NI; t++) { Expr* r = ref_eval_arr(b, inm, &v, 1); expr_free(r); }
            double ti = (double)(clock() - t0) / CLOCKS_PER_SEC / NI;
            printf("ok:   %-30s len=%-5zu %.1fx faster (%.2f us vs %.2f us)\n",
                   "array performance", len, ti / tc, tc * 1e6, ti * 1e6);
            (void)acc;
            expr_free(v); compiled_free(p); expr_free(b);
        }

        /* An array temporary produced inside a Do loop must be freed per
         * iteration, not accumulated until teardown. */
        {
            const char* inm[1] = { intern_symbol("v") };
            Expr* b = parse_expression("Module[{s = 0.0}, Do[s = s + Total[Sin[v]], {i, 1, 2000}]; s]");
            CompiledProgram* p = compile_expr(b, inm, AA, 1);
            Expr* v = make_vec(16, 0.3, 3.0);
            CompileValue a, o; a.type = AA[0]; a.v.a = v;
            bool ok = p && compiled_eval(p, &a, &o);
            if (!ok) { printf("FAIL: %-30s -> loop over array temps failed\n", "array temp freed in loop"); failures++; }
            else printf("ok:   %-30s s=%.3f\n", "array temp freed in loop", o.v.r);
            expr_free(v); compiled_free(p); expr_free(b);
        }
    }

    /* ---- control flow: Sum / Product (integer-counted loops) ----
     * Parsed WITHOUT evaluation so the loop stays symbolic (the future Compile[]
     * builtin is HoldAll); references are computed directly in C. */
    {
        const char* inm[2] = { intern_symbol("x"), intern_symbol("n") };
        CompileType ty[2] = { CT_REAL, CT_INT };
        struct { const char* nm; const char* body; } cases[] = {
            { "Sum Sin", "Sum[Sin[i x], {i, 1, n}]" },
            { "Sum i^2 x", "Sum[i^2 x, {i, 1, n}]" },
            { "Sum with If", "Sum[If[Mod[i,2] == 0, i x, -i x], {i, 1, n}]" },
            { "Product 1+x/i", "Product[1 + x/i, {i, 1, n}]" },
            { "Sum Gamma[i]", "Sum[Gamma[i] x, {i, 1, n}]" },
            { "nested Sum", "Sum[Sum[i j x, {j, 1, i}], {i, 1, n}]" },
        };
        for (size_t k = 0; k < sizeof cases / sizeof cases[0]; k++) {
            Expr* b = parse_expression(cases[k].body);   /* NOT evaluated */
            CompiledProgram* p = compile_expr(b, inm, ty, 2);
            if (!p) { printf("FAIL: %-30s -> did not compile\n", cases[k].nm); failures++; expr_free(b); continue; }
            double maxerr = 0; int cmp = 0;
            for (int t = 0; t < 200; t++) {
                double x = urand(0.2, 1.5); long long n = irand(1, 9);
                CompileValue av[2] = { { CT_REAL, { .r = x } }, { CT_INT, { .i = n } } }, o;
                if (!compiled_eval(p, av, &o)) continue;
                double ref;
                if (k == 0) { ref = 0; for (long long i = 1; i <= n; i++) ref += sin(i * x); }
                else if (k == 1) { ref = 0; for (long long i = 1; i <= n; i++) ref += (double)(i * i) * x; }
                else if (k == 2) { ref = 0; for (long long i = 1; i <= n; i++) ref += (i % 2 == 0 ? 1 : -1) * i * x; }
                else if (k == 3) { ref = 1; for (long long i = 1; i <= n; i++) ref *= 1 + x / (double)i; }
                else if (k == 4) { ref = 0; for (long long i = 1; i <= n; i++) ref += tgamma((double)i) * x; }
                else { ref = 0; for (long long i = 1; i <= n; i++) for (long long j = 1; j <= i; j++) ref += (double)(i * j) * x; }
                double err = fabs(o.v.r - ref) / (1.0 + fabs(ref));
                if (err > maxerr) maxerr = err;
                cmp++;
            }
            if (cmp < 50 || maxerr > 1e-9) { printf("FAIL: %-30s -> max_rel=%.2e (%d)\n", cases[k].nm, maxerr, cmp); failures++; }
            else printf("ok:   %-30s max_rel=%.1e (%d cmps)\n", cases[k].nm, maxerr, cmp);
            compiled_free(p); expr_free(b);
        }
    }

    /* ---- procedural: Module/With locals, Set/AddTo/Increment, Do/While/For,
     * CompoundExpression.  Parsed unevaluated; C references. ---- */
    {
        const char* in1[1] = { intern_symbol("x") };
        const char* in2[2] = { intern_symbol("x"), intern_symbol("n") };
        CompileType ty1[1] = { CT_REAL };
        CompileType ty2[2] = { CT_REAL, CT_INT };
        struct { const char* nm; const char* body; int na; } cases[] = {
            { "Module+Set+CompoundExpr", "Module[{s = 0., t = 1.}, s = s + x; t = t x; s + t]", 1 },
            { "Do accumulation",         "Module[{s = 0.}, Do[s = s + 1/i^2, {i, 1, n}]; s]", 2 },
            { "While Newton sqrt",       "Module[{r = x, k = 0}, While[k < 25, r = (r + x/r)/2; k = k + 1]; r]", 1 },
            { "For sum 1..n",            "Module[{s = 0, i = 0}, For[i = 1, i <= n, i = i + 1, s = s + i]; s]", 2 },
            { "Increment + i^2 sum",     "Module[{s = 0, i = 0}, While[i < n, i++; s = s + i i]; s]", 2 },
            { "With + AddTo/TimesBy",    "With[{a = x}, Module[{s = 0.}, s += a; TimesBy[s, 3]; s]]", 1 },
        };
        for (size_t k = 0; k < sizeof cases / sizeof cases[0]; k++) {
            Expr* b = parse_expression(cases[k].body);
            const char* const* nm = cases[k].na == 1 ? in1 : in2;
            const CompileType* ty = cases[k].na == 1 ? ty1 : ty2;
            CompiledProgram* p = compile_expr(b, nm, ty, (size_t)cases[k].na);
            if (!p) { printf("FAIL: %-30s -> did not compile\n", cases[k].nm); failures++; expr_free(b); continue; }
            double maxerr = 0; int cmp = 0;
            for (int t = 0; t < 150; t++) {
                double x = urand(0.5, 3.0); long long n = irand(1, 12);
                CompileValue av[2] = { { CT_REAL, { .r = x } }, { CT_INT, { .i = n } } }, o;
                if (!compiled_eval(p, av, &o)) continue;
                double got = (o.type == CT_INT) ? (double)o.v.i : o.v.r, ref = 0;
                if (k == 0) ref = 2 * x;
                else if (k == 1) { for (long long i = 1; i <= n; i++) ref += 1.0 / (double)(i * i); }
                else if (k == 2) ref = sqrt(x);
                else if (k == 3) { for (long long i = 1; i <= n; i++) ref += (double)i; }
                else if (k == 4) { for (long long i = 1; i <= n; i++) ref += (double)(i * i); }
                else if (k == 5) ref = 3.0 * x;   /* With+AddTo/TimesBy: (0+x)*3 */
                else if (k == 8) ref = 3.0 * x;    /* diag TimesBy int */
                else ref = x;                       /* diag: return local == x */
                double err = fabs(got - ref) / (1.0 + fabs(ref));
                if (err > maxerr) maxerr = err;
                cmp++;
            }
            if (cmp < 40 || maxerr > 1e-9) { printf("FAIL: %-30s -> max_rel=%.2e (%d)\n", cases[k].nm, maxerr, cmp); failures++; }
            else printf("ok:   %-30s max_rel=%.1e (%d cmps)\n", cases[k].nm, maxerr, cmp);
            compiled_free(p); expr_free(b);
        }
    }

    /* ---- Counted iterator spellings.  Do accepts every form the interpreter's
     * iter_spec_parse does; Sum/Product require a named iterator because the
     * interpreter rejects a bare count for them.  A missing form here is not a
     * slow path — it makes the WHOLE surrounding Compile[] bail, which is how
     * Do[body, {n}] cost the Newton-fractal benchmark ~37x. ---- */
    {
        const char* in1[1] = { intern_symbol("n") };
        CompileType ty1[1] = { CT_INT };
        struct { const char* nm; const char* body; int compiles; } fc[] = {
            { "Do {n} bare count",    "Module[{s = 0}, Do[s = s + 1, {n}]; s]",           1 },
            { "Do n unbraced count",  "Module[{s = 0}, Do[s = s + 1, n]; s]",             1 },
            { "Do {i, hi}",           "Module[{s = 0}, Do[s = s + i, {i, n}]; s]",        1 },
            { "Do {i, lo, hi}",       "Module[{s = 0}, Do[s = s + i, {i, 1, n}]; s]",     1 },
            { "Do {i,lo,hi,+di}",     "Module[{s = 0}, Do[s = s + 1, {i, 1, n, 2}]; s]",  1 },
            { "Do {i,hi,lo,-di}",     "Module[{s = 0}, Do[s = s + 1, {i, n, 1, -1}]; s]", 1 },
            { "Sum {i, hi}",          "Sum[i, {i, n}]",                                   1 },
            { "Sum {n} rejected",     "Sum[2, {n}]",                                       0 },
            { "Product {n} rejected", "Product[2, {n}]",                                   0 },
        };
        for (size_t k = 0; k < sizeof fc / sizeof fc[0]; k++) {
            Expr* b = parse_expression(fc[k].body);
            CompiledProgram* p = compile_expr(b, in1, ty1, 1);
            if (!fc[k].compiles) {
                /* Must bail: the interpreter leaves these unevaluated, so a
                 * compiled answer would be a divergence, not a bonus. */
                if (p) { printf("FAIL: %-30s -> compiled but must bail\n", fc[k].nm); failures++; compiled_free(p); }
                else printf("ok:   %-30s bails (matches interpreter)\n", fc[k].nm);
                expr_free(b); continue;
            }
            if (!p) { printf("FAIL: %-30s -> did not compile\n", fc[k].nm); failures++; expr_free(b); continue; }
            int bad = 0, cmp = 0;
            for (long long n = 0; n <= 12; n++) {
                CompileValue av[1] = { { CT_INT, { .i = n } } }, o;
                if (!compiled_eval(p, av, &o) || o.type != CT_INT) { bad++; continue; }
                long long ref = 0;
                switch (k) {
                    case 0: case 1: ref = n > 0 ? n : 0; break;                       /* count */
                    case 2: case 3: for (long long i = 1; i <= n; i++) ref += i; break;
                    case 4: for (long long i = 1; i <= n; i += 2) ref++; break;
                    case 5: for (long long i = n; i >= 1; i--) ref++; break;
                    case 6: for (long long i = 1; i <= n; i++) ref += i; break;
                    default: break;
                }
                if (o.v.i != ref) { bad++; printf("      n=%lld got=%lld want=%lld\n", n, o.v.i, ref); }
                cmp++;
            }
            if (bad || cmp < 13) { printf("FAIL: %-30s -> %d bad of %d\n", fc[k].nm, bad, cmp); failures++; }
            else printf("ok:   %-30s %d counts exact\n", fc[k].nm, cmp);
            compiled_free(p); expr_free(b);
        }
    }

    /* ---- Multi-statement loop body, written three ways.  A one-liner only
     * exercises the loop scaffolding; this body carries four locals of three
     * different types (Complex u/du, Real acc, Integer cnt), five statements per
     * iteration, and a two-armed If that mutates state in BOTH arms — so the
     * accumulator/register discipline is under real pressure across the back
     * edge.  Do / While / For must agree BITWISE (same arithmetic, only the loop
     * control differs) and must match an independent C reference. ---- */
    {
        const char* inzn[2] = { intern_symbol("z"), intern_symbol("n") };
        CompileType tyzn[2] = { CT_COMPLEX, CT_INT };
        /* the shared five statements, spliced into each loop form */
        #define ML_STEP "du = (2 u + 1/u^2)/3 - u; u = u + du; acc = acc + Abs[du]; " \
                        "If[Re[u] > 0, cnt = cnt + 1, cnt = cnt - 1]"
        #define ML_DECL "u = z, du = 0. + 0. I, acc = 0., cnt = 0"
        #define ML_TAIL "acc + cnt + Re[u]"
        struct { const char* nm; const char* body; } ml[] = {
            { "multi-line Do",    "Module[{" ML_DECL "}, Do[" ML_STEP ", {n}]; " ML_TAIL "]" },
            { "multi-line While", "Module[{" ML_DECL ", k = 0}, While[k < n, " ML_STEP "; k = k + 1]; " ML_TAIL "]" },
            { "multi-line For",   "Module[{" ML_DECL ", k = 0}, For[k = 0, k < n, k = k + 1, " ML_STEP "]; " ML_TAIL "]" },
        };
        #undef ML_STEP
        #undef ML_DECL
        #undef ML_TAIL
        CompiledProgram* mp[3] = { NULL, NULL, NULL };
        for (size_t k = 0; k < 3; k++) {
            Expr* b = parse_expression(ml[k].body);   /* raw: evaluating would rewrite the loop */
            mp[k] = compile_expr(b, inzn, tyzn, 2);
            if (!mp[k]) { printf("FAIL: %-30s -> did not compile\n", ml[k].nm); failures++; }
            expr_free(b);
        }
        if (mp[0] && mp[1] && mp[2]) {
            double maxrel = 0; int cmp = 0, disagree = 0;
            for (int t = 0; t < 200; t++) {
                /* stay off the Julia boundary: the iteration is chaotic there and
                 * would amplify the reference's own rounding without telling us
                 * anything about the lowering */
                double re = urand(0.35, 1.4), im = urand(0.35, 1.4);
                if (t & 1) re = -re;
                if (t & 2) im = -im;
                long long n = irand(0, 14);
                double _Complex z = re + im * I;
                CompileValue av[2] = { { CT_COMPLEX, { .z = 0 } }, { CT_INT, { .i = n } } };
                av[0].v.z = z;
                double got[3];
                int ok3 = 1;
                for (size_t k = 0; k < 3; k++) {
                    CompileValue o;
                    if (!compiled_eval(mp[k], av, &o) || o.type != CT_REAL) { ok3 = 0; break; }
                    got[k] = o.v.r;
                }
                if (!ok3) continue;
                /* the three spellings are the same computation -> bitwise equal */
                if (got[0] != got[1] || got[0] != got[2]) {
                    if (disagree++ < 3)
                        printf("      n=%lld z=%g%+gi  Do=%.17g While=%.17g For=%.17g\n",
                               n, re, im, got[0], got[1], got[2]);
                }
                /* independent C reference */
                double _Complex u = z, du = 0;
                double acc = 0; long long cnt = 0;
                for (long long q = 0; q < n; q++) {
                    du = (2.0 * u + 1.0 / (u * u)) / 3.0 - u;
                    u = u + du;
                    acc += cabs(du);
                    if (creal(u) > 0) cnt++; else cnt--;
                }
                double ref = acc + (double)cnt + creal(u);
                double rel = fabs(got[0] - ref) / (1.0 + fabs(ref));
                if (rel > maxrel) maxrel = rel;
                cmp++;
            }
            if (disagree) { printf("FAIL: multi-line loops disagree (%d of %d)\n", disagree, cmp); failures++; }
            else if (cmp < 60 || maxrel > 1e-9) {
                printf("FAIL: multi-line loop body -> max_rel=%.2e (%d cmps)\n", maxrel, cmp); failures++;
            } else {
                printf("ok:   %-30s Do==While==For bitwise, max_rel=%.1e (%d cmps)\n",
                       "multi-line loop body", maxrel, cmp);
            }
        }
        for (size_t k = 0; k < 3; k++) if (mp[k]) compiled_free(mp[k]);
    }

    /* ---- Nest[Function[u, body], x, n]: parsed unevaluated; C references. ---- */
    {
        const char* in1[1] = { intern_symbol("x") };
        const char* in2[2] = { intern_symbol("x"), intern_symbol("n") };
        CompileType ty1[1] = { CT_REAL };
        CompileType ty2[2] = { CT_REAL, CT_INT };
        struct { const char* nm; const char* body; int na; } nc[] = {
            { "Nest linear map",       "Nest[Function[u, u/2 + 1], x, n]", 2 },
            { "Nest Newton sqrt",      "Nest[Function[{u}, (u + x/u)/2], x, 8]", 1 },
            { "Nest 2^n (int)",        "Nest[Function[u, 2 u], 1, n]", 2 },
            { "Nest widen int->real",  "Nest[Function[u, u + 0.5], 0, n]", 2 },
        };
        for (size_t k = 0; k < sizeof nc / sizeof nc[0]; k++) {
            Expr* b = parse_expression(nc[k].body);
            const char* const* nm = nc[k].na == 1 ? in1 : in2;
            const CompileType* ty = nc[k].na == 1 ? ty1 : ty2;
            CompiledProgram* p = compile_expr(b, nm, ty, (size_t)nc[k].na);
            if (!p) { printf("FAIL: %-30s -> did not compile\n", nc[k].nm); failures++; expr_free(b); continue; }
            double maxerr = 0; int cmp = 0;
            for (int t = 0; t < 150; t++) {
                double x = urand(1.0, 3.0); long long n = irand(0, 10);
                CompileValue av[2] = { { CT_REAL, { .r = x } }, { CT_INT, { .i = n } } }, o;
                if (!compiled_eval(p, av, &o)) continue;
                double got = (o.type == CT_INT) ? (double)o.v.i : o.v.r, ref = 0;
                if (k == 0) { ref = x; for (long long i = 0; i < n; i++) ref = ref / 2 + 1; }
                else if (k == 1) { ref = x; for (int i = 0; i < 8; i++) ref = (ref + x / ref) / 2; }
                else if (k == 2) { long long r = 1; for (long long i = 0; i < n; i++) r *= 2; ref = (double)r; }
                else { ref = 0; for (long long i = 0; i < n; i++) ref += 0.5; }
                double err = fabs(got - ref) / (1.0 + fabs(ref));
                if (err > maxerr) maxerr = err;
                cmp++;
            }
            if (cmp < 40 || maxerr > 1e-9) { printf("FAIL: %-30s -> max_rel=%.2e (%d)\n", nc[k].nm, maxerr, cmp); failures++; }
            else printf("ok:   %-30s max_rel=%.1e (%d cmps)\n", nc[k].nm, maxerr, cmp);
            compiled_free(p); expr_free(b);
        }
    }

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

    /* dispatch-bound micro-benchmark: a deep Horner polynomial is pure add/mul,
     * so VM instruction dispatch (not libm) dominates — the regime where a
     * threaded interpreter helps and which matches a stencil RHS. */
    {
        char buf[4096]; size_t pos = 0;
        const int DEG = 40;
        for (int i = 0; i < DEG; i++) pos += (size_t)snprintf(buf + pos, sizeof buf - pos, "(");
        pos += (size_t)snprintf(buf + pos, sizeof buf - pos, "x");
        for (int i = 0; i < DEG; i++) pos += (size_t)snprintf(buf + pos, sizeof buf - pos, " + %d) x", (i % 9) + 1);
        const char* inm[1] = { intern_symbol("x") };
        Expr* b = parse_expression(buf);           /* NOT evaluated: keep Horner nesting */
        CompiledProgram* p = compile_expr(b, inm, RRR, 1);
        if (!p) { printf("FAIL: dispatch-perf body did not compile\n"); failures++; }
        else {
            const int NC = 5000000; double acc = 0, a;
            clock_t t0 = clock();
            for (int i = 0; i < NC; i++) { a = 0.3 + 1e-9 * i; double o; compiled_eval_real(p, &a, &o); acc += o; }
            double per = (double)(clock() - t0) / CLOCKS_PER_SEC / NC;
            printf("ok:   %-30s %.2f ns/call (%d arith ops)\n", "dispatch perf", per * 1e9, 2 * DEG);
            (void)acc;
        }
        compiled_free(p); expr_free(b);
    }

    /* ================= OPTIMISER A/B =================
     * The bytecode optimiser (constant folding, CSE, copy propagation, DCE,
     * LICM) is required to be RESULT-PRESERVING, not merely accurate: it must
     * never reassociate floating point or contract a multiply-add.  So the gate
     * is bitwise identity between a body compiled with and without the passes,
     * over a randomised argument sweep — plus bail parity, since a pass must not
     * change which bodies compile at all. */
    {
        long long tot_opt = 0, tot_raw = 0;
        int ab_bodies = 0;
        ab_opt_init(&tot_opt, &tot_raw, &ab_bodies);
        #define AB(nm, src, na, ty, raw) ab_opt(nm, src, xyz, ty, na, raw)

        AB("straight-line arith",   "x + 2 y - x y + 3",            2, RRR, false);
        AB("shared subexpression",  "Sin[x y] + Cos[x y] + Sin[x y]^2", 2, RRR, false);
        AB("repeated power",        "x^2 + x^2 y + (x^2)^3",        2, RRR, false);
        AB("constant subtree",      "x + 2 Pi + Sqrt[2] Exp[1]",    1, RRR, false);
        AB("nested constants",      "x (1 + 2) (3 - 1) / 4",        1, RRR, false);
        AB("libm chain",            "Sin[Exp[-x^2/2]] + Sqrt[Abs[x]] - Log[1 + x^2]", 1, RRR, false);
        AB("complex arith",         "(x + I y)^3 + Exp[x + I y]",   2, CCC, false);
        AB("int arith",             "x y + Mod[x, 7] - Quotient[y, 3]", 2, III, false);
        AB("comparisons",           "If[x < y, x^2 + Sin[y], y^2 - Cos[x]]", 2, RRR, false);
        AB("kernel heads",          "Gamma[x] + BesselJ[2, y] + Erf[x]", 2, RRR, false);
        AB("Max/Min",               "Max[x, y] - Min[x, y] + Max[x, 2]", 2, RRR, false);
        /* Loop bodies: the LICM pass only has anything to do here.  Parsed
         * UNEVALUATED, or the interpreter closed-forms the Sum before we see it. */
        AB("Do loop",               "Module[{s = 0.}, Do[s = s + Sin[x] i, {i, 1, 20}]; s]", 1, RRR, true);
        AB("Sum invariant body",    "Sum[Sin[x] Cos[y] i, {i, 1, 30}]", 2, RRR, true);
        AB("nested loops",          "Sum[Sum[x y i j, {j, 1, 8}], {i, 1, 8}]", 2, RRR, true);
        AB("While loop",            "Module[{t = x, k = 0}, While[k < 12, t = (t + x/t)/2; k = k + 1]; t]", 1, RRR, true);
        AB("For loop",              "Module[{s = 0.}, For[i = 1, i <= 15, i = i + 1, s = s + x^2 i]; s]", 1, RRR, true);
        AB("Nest",                  "Nest[Function[u, (u + x/u)/2], x, 14]", 1, RRR, true);
        AB("With locals",           "With[{a = Sin[x], b = Cos[x]}, a b + a/b + a^2]", 1, RRR, true);
        AB("loop-invariant heavy",  "Sum[Exp[-x^2] Sqrt[Abs[y]] + i, {i, 1, 25}]", 2, RRR, true);

        #undef AB
        if (ab_bodies)
            printf("ok:   %-30s %d bodies bitwise-identical, %lld -> %lld instrs (%.0f%% removed)\n",
                   "optimiser A/B", ab_bodies, tot_raw, tot_opt,
                   100.0 * (double)(tot_raw - tot_opt) / (double)(tot_raw ? tot_raw : 1));
    }

    /* Arrays must survive the optimiser untouched — their ownership lives in
     * instruction flags, not in the dataflow the scalar passes reason about. */
    {
        const char* vw[] = { "v", "w" };
        const CompileType AA[] = { CT_ARRAY(CT_REAL, 1), CT_ARRAY(CT_REAL, 1) };
        parity_arr("opt: vec chain",  "Total[Sin[v] Exp[-v] + Sqrt[v]]", vw, AA, 1, 64, 0.3, 4.0, 20);
        parity_arr("opt: vec + vec",  "v w + v - w",  vw, AA, 2, 64, 0.3, 4.0, 20);
    }

    if (failures == 0) printf("\nAll Compile engine tests passed.\n");
    else printf("\n%d Compile engine test(s) FAILED.\n", failures);
    return failures ? 1 : 0;
}
