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

/* When set, a compiled array call that DECLINES (returns false) is accepted
 * rather than failing the test.  Declining is the documented contract for a
 * real-typed program whose buffer would have gone complex, or that hit a pole —
 * the caller falls back to the interpreter — so randomised bodies, which hit
 * those cases by construction, must not treat it as an error. */
static bool allow_decline = false;


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
/* Rank-N real array with `rank` dims taken from `dims`. */
static Expr* make_nd(int rank, const int64_t* dims, double lo, double hi) {
    size_t n = 1;
    for (int i = 0; i < rank; i++) n *= (size_t)dims[i];
    double* buf = malloc(n * sizeof(double));
    for (size_t k = 0; k < n; k++) buf[k] = urand(lo, hi);
    return expr_new_ndarray(rank, dims, buf, NDT_FLOAT64);   /* adopts buf */
}

/* Parity for rank >= 2.  The delegated NDArray path is already rank-general, so
 * lifting the compiler's rank-1 front gate is all that higher rank needs; this
 * checks that claim against the interpreter for the same body. */
static void parity_nd(const char* name, const char* body_s,
                      const char* const* names, size_t n, int rank,
                      const int64_t* dims, unsigned flags, int trials) {
    Expr* body = parse_expression(body_s);
    const char* inames[4];
    CompileType types[4];
    for (size_t k = 0; k < n; k++) {
        inames[k] = intern_symbol(names[k]);
        types[k] = CT_ARRAY(CT_REAL, rank);
    }
    CompiledProgram* p = compile_expr_ex(body, inames, types, n, flags);
    if (!p) { printf("FAIL: %-30s -> did not compile (rank %d)\n", name, rank); failures++; expr_free(body); return; }

    int cmp = 0; double maxerr = 0.0; bool shape_ok = true;
    for (int t = 0; t < trials; t++) {
        Expr* vecs[4]; CompileValue args[4], outc;
        for (size_t k = 0; k < n; k++) {
            vecs[k] = make_nd(rank, dims, 0.4, 3.0);
            args[k].type = types[k]; args[k].v.a = vecs[k];
        }
        if (compiled_eval(p, args, &outc)) {
            Expr* want = ref_eval_arr(body, inames, vecs, n);
            Expr* got  = aval_to_expr(outc);
            if (!arr_cmp(got, want, &maxerr)) shape_ok = false; else cmp++;
            expr_free(got); expr_free(want);
        }
        for (size_t k = 0; k < n; k++) expr_free(vecs[k]);
        if (!shape_ok) break;
    }
    if (!shape_ok) { printf("FAIL: %-30s -> rank-%d shape/kind mismatch\n", name, rank); failures++; }
    else if (cmp < trials && !allow_decline) { printf("FAIL: %-30s -> only %d/%d evaluated\n", name, cmp, trials); failures++; }
    else if (cmp == 0 && allow_decline) { /* declined every trial: the contract, not a bug */ }
    else if (maxerr > 1e-12) { printf("FAIL: %-30s -> max_rel=%.2e\n", name, maxerr); failures++; }
    else printf("ok:   %-30s rank %d, max_rel=%.1e (%d cmps)\n", name, rank, maxerr, cmp);
    compiled_free(p); expr_free(body);
}

static void parity_arr_ex(const char* name, const char* body_s,
                       const char* const* names, const CompileType* types, size_t n,
                       size_t len, double lo, double hi, int trials, unsigned flags) {
    Expr* body = parse_expression(body_s);
    const char* inames[4];
    for (size_t k = 0; k < n; k++) inames[k] = intern_symbol(names[k]);
    CompiledProgram* p = compile_expr_ex(body, inames, types, n, flags);
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
    else if (cmp < trials && !allow_decline) { printf("FAIL: %-30s -> only %d/%d evaluated\n", name, cmp, trials); failures++; }
    else if (cmp == 0 && allow_decline) { /* declined every trial: the contract, not a bug */ }
    else if (maxerr > 1e-12) { printf("FAIL: %-30s -> max_rel=%.2e\n", name, maxerr); failures++; }
    else printf("ok:   %-30s max_rel=%.1e (%d cmps)\n", name, maxerr, cmp);
    compiled_free(p);
    expr_free(body);
}

static void parity_arr(const char* name, const char* body_s,
                       const char* const* names, const CompileType* types, size_t n,
                       size_t len, double lo, double hi, int trials) {
    parity_arr_ex(name, body_s, names, types, n, len, lo, hi, trials, 0u);
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

/* Evaluate `src` (which mentions the symbol `xq`) at x = xv.
 *
 * Substitution, NOT string formatting: printing a double with %.17g yields a
 * 17-digit literal, which Mathilda reads as an ARBITRARY-PRECISION number, so
 * the whole reference computes in MPFR and comes back as EXPR_MPFR rather than a
 * machine Real.  Binding a genuine machine Real sidesteps that entirely. */
static Expr* ref_at(const char* src, double xv) {
    Expr* tmpl = parse_expression(src);
    MatchEnv* env = env_new();
    Expr* val = expr_new_real(xv);
    env_set(env, intern_symbol("xq"), val);
    Expr* sub = replace_bindings(tmpl, env);
    env_free(env);
    expr_free(val);
    expr_free(tmpl);
    return eval_and_free(sub);
}

/* Does the compiled program's DECLARED RESULT TYPE agree with the HEAD the
 * interpreter produces for the same call?
 *
 * The numeric parity tests cannot see this.  They compare values, and a value
 * compares equal whether it came back as the Integer -1 or the Real -1. — but
 * the two are not interchangeable downstream (`IntegerQ`, exact arithmetic,
 * printing all differ), so a compiled path that answers with the wrong head is
 * answering DIFFERENTLY from the interpreter, which the engine forbids.
 *
 * This is not hypothetical: it is why `UnitStep` is lowered by hand rather than
 * registered as a kernel, and `Sign[-2.5]` shipped as `-1.` against the
 * interpreter's `-1` until this check existed. */
static void parity_head(const char* name, const char* src, double xv) {
    Expr* body = parse_expression(src);
    const char* inm[1] = { intern_symbol("xq") };
    const CompileType RR[1] = { CT_REAL };
    CompiledProgram* p = compile_expr(body, inm, RR, 1);
    if (!p) { printf("FAIL: head %-24s did not compile\n", name); failures++; expr_free(body); return; }

    Expr* want = ref_at(src, xv);
    const char* wh = want->type == EXPR_INTEGER ? "Integer"
                   : want->type == EXPR_REAL    ? "Real"
                   : (want->type == EXPR_FUNCTION
                      && want->data.function.head->type == EXPR_SYMBOL
                      && strcmp(want->data.function.head->data.symbol.name, "Complex") == 0)
                     ? "Complex" : "other";
    CompileType rt = compiled_result_type(p);
    const char* gh = rt == CT_INT ? "Integer" : rt == CT_REAL ? "Real"
                   : rt == CT_COMPLEX ? "Complex" : "other";
    /* A Complex whose imaginary part is exactly zero is unboxed back to a Real,
     * so CT_COMPLEX legitimately produces a Real head at some arguments. */
    bool ok = strcmp(wh, gh) == 0
              || (rt == CT_COMPLEX && strcmp(wh, "Real") == 0);
    if (!ok) {
        printf("FAIL: head %-24s at x=%g: interpreter gives %s, compiled declares %s\n",
               name, xv, wh, gh);
        failures++;
    } else {
        printf("ok:   head %-24s %s\n", name, gh);
    }
    expr_free(want); compiled_free(p); expr_free(body);
}

/* ---- strip-mining stress helpers ---------------------------------------- */

/* Compile the same array body fused and delegated, and require agreement.  Not
 * bitwise: the delegated reduction sums pairwise while the fused one accumulates
 * a tile at a time, so the two legitimately differ in the last bits. */
static bool fused_matches_delegated(const char* body_s, const char* const* names,
                                    const CompileType* types, size_t n, size_t len) {
    const char* inames[4];
    for (size_t k = 0; k < n; k++) inames[k] = intern_symbol(names[k]);
    Expr* body = parse_expression(body_s);
    CompiledProgram* pf = compile_expr_ex(body, inames, types, n, 0u);
    CompiledProgram* pd = compile_expr_ex(body, inames, types, n, COMPILE_NO_FUSE);
    bool ok = true;
    if (!pf) ok = (pd == NULL);          /* fusion must not be the weaker path */
    else if (!pd) ok = true;             /* fusion compiles strictly more: fine */
    else {
        Expr* vecs[4]; CompileValue args[4], of, od;
        for (size_t k = 0; k < n; k++) {
            vecs[k] = make_vec(len, 0.4, 3.0);
            args[k].type = types[k]; args[k].v.a = vecs[k];
        }
        bool sf = compiled_eval(pf, args, &of), sd = compiled_eval(pd, args, &od);
        if (sf != sd) ok = false;
        else if (sf) {
            if (CT_IS_ARRAY(of.type) != CT_IS_ARRAY(od.type)) ok = false;
            else if (CT_IS_ARRAY(of.type)) {
                double err = 0;
                Expr* a = aval_to_expr(of); Expr* b2 = aval_to_expr(od);
                if (!arr_cmp(a, b2, &err) || err > 1e-12) ok = false;
                expr_free(a); expr_free(b2);
            } else if (fabs(of.v.r - od.v.r) > 1e-11 * (fabs(od.v.r) + 1.0)) ok = false;
            if (CT_IS_ARRAY(of.type)) { } /* aval_to_expr already consumed copies */
        }
        for (size_t k = 0; k < n; k++) expr_free(vecs[k]);
    }
    compiled_free(pf); compiled_free(pd); expr_free(body);
    return ok;
}

/* A random elementwise body over `v` and `w`.
 *
 * ONLY Listable heads appear — that is precisely the set fusion may thread over,
 * so the generator cannot produce something the compiler is right to refuse.
 * Two further constraints keep every generated body a legitimate test rather
 * than a known-uncompilable one: at least one ARRAY leaf must reach the root
 * (`want_arr`), because a body of pure scalars is not an array program at all
 * and a bare argument array cannot be the result (the caller would end up
 * freeing a value it does not own); and `Total` only ever wraps the finished
 * tree, because Total of a scalar is not compilable and never should be. */
static void rand_arr_tree(char* buf, size_t cap, int depth, bool want_arr) {
    if (depth <= 0) {
        if (want_arr) snprintf(buf, cap, "%s", irand(0, 1) ? "v" : "w");
        else {
            /* Scalar leaves must be INEXACT.  With an exact integer the
             * interpreter keeps Cos[3] symbolic while the compiler folds it to a
             * machine number, so the reference and the compiled result would
             * differ in kind for reasons that have nothing to do with fusion. */
            static const char* leaf[] = { "v", "w", "2.", "0.5", "3.25" };
            snprintf(buf, cap, "%s", leaf[irand(0, 4)]);
        }
        return;
    }
    char a[256], b[256];
    int pick = (int)irand(0, 8);
    if (pick <= 4) {                       /* binary: one side carries the array */
        static const char* op[] = { "+", "-", "*", "/", "+" };
        bool left = irand(0, 1) != 0;
        rand_arr_tree(a, sizeof a, depth - 1, want_arr && left);
        rand_arr_tree(b, sizeof b, depth - 1, want_arr && !left);
        snprintf(buf, cap, "(%s) %s (%s)", a, op[pick], b);
    } else if (pick <= 6) {                /* unary Listable head */
        static const char* fn[] = { "Sin", "Cos", "Exp", "Sqrt", "Abs", "Tanh" };
        rand_arr_tree(a, sizeof a, depth - 1, want_arr);
        snprintf(buf, cap, "%s[%s]", fn[irand(0, 5)], a);
    } else {                               /* integer power */
        rand_arr_tree(a, sizeof a, depth - 1, want_arr);
        snprintf(buf, cap, "(%s)^%d", a, (int)irand(2, 4));
    }
}

static void rand_arr_body(char* buf, size_t cap, int depth) {
    char inner[400];
    rand_arr_tree(inner, sizeof inner, depth, true);
    if (irand(0, 3) == 0) snprintf(buf, cap, "Total[%s]", inner);
    else                  snprintf(buf, cap, "(%s) + (v)", inner);
}

/* ---- optimiser A/B ------------------------------------------------------
 * Compile the same body twice, with and without the optimiser, and require the
 * two programs to agree BITWISE over a randomised argument sweep.  Bitwise (not
 * "to rounding") is the right gate: a pass that reassociated a sum or contracted
 * a multiply-add would still look accurate but would break the engine's stated
 * parity contract with the interpreter. */
/* Which pass the A/B is currently disabling, so one harness gates them all. */
static unsigned    ab_off_flag = 0;
static const char* ab_what = "the optimiser";

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
    CompiledProgram* pr = compile_expr_ex(body, inames, types, n, ab_off_flag);

    if (!po != !pr) {
        printf("FAIL: %-30s %s changed whether the body compiles\n", name, ab_what);
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
        printf("FAIL: %-30s result differs with %s\n", name, ab_what);
        failures++;
    } else if (no > nr) {
        printf("FAIL: %-30s %s GREW the program (%zu -> %zu)\n", name, ab_what, nr, no);
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
    const char* x1i[] = { intern_symbol("x") };
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
    /* Heads the coverage audit caught bailing even though the interpreter
     * evaluates them: real Mod/Quotient and Arg of a real.  Each one bailing
     * cost the WHOLE surrounding body, silently. */
    parity("real Mod", "Mod[x, y]", xyz, RRR, 2, 0.3, 6.0, 0, 0, 300);
    parity("real Mod negative", "Mod[-x, y]", xyz, RRR, 2, 0.3, 6.0, 0, 0, 300);
    parity("real Quotient", "Quotient[x, y]", xyz, RRR, 2, 0.3, 6.0, 0, 0, 300);
    parity("real Arg (positive)", "Arg[x]", xyz, RRR, 1, 0.3, 6.0, 0, 0, 200);
    parity("real Arg (negative)", "Arg[-x]", xyz, RRR, 1, 0.3, 6.0, 0, 0, 200);
    parity("Mod/Arg in a chain", "Sin[Mod[x, y]] + Arg[-x] Quotient[x, y]", xyz, RRR, 2, 0.4, 5.0, 0, 0, 200);

    /* Exponential-integral family: new machine kernels (sf_machine.c) for
     * modules that previously computed only in MPFR, so these heads used to bail
     * and take the whole surrounding body with them.  Ranges stay inside each
     * function's real domain; the out-of-domain behaviour is checked separately
     * below, where the kernel must DECLINE rather than invent a real answer. */
    parity("ExpIntegralEi",   "ExpIntegralEi[x]",  xyz, RRR, 1, 0.3, 12.0, 0, 0, 300);
    parity("ExpIntegralEi -",  "ExpIntegralEi[-x]", xyz, RRR, 1, 0.3, 12.0, 0, 0, 300);
    parity("ExpIntegralEi big","ExpIntegralEi[x]",  xyz, RRR, 1, 45.0, 90.0, 0, 0, 200);
    parity("LogIntegral",     "LogIntegral[x]",    xyz, RRR, 1, 1.5, 500.0, 0, 0, 300);
    parity("LogIntegral <1",  "LogIntegral[x]",    xyz, RRR, 1, 0.05, 0.9, 0, 0, 200);
    parity("SinIntegral",     "SinIntegral[x]",    xyz, RRR, 1, -40.0, 40.0, 0, 0, 400);
    parity("CosIntegral",     "CosIntegral[x]",    xyz, RRR, 1, 0.05, 40.0, 0, 0, 400);
    parity("SinhIntegral",    "SinhIntegral[x]",   xyz, RRR, 1, -25.0, 25.0, 0, 0, 300);
    parity("CoshIntegral",    "CoshIntegral[x]",   xyz, RRR, 1, 0.05, 25.0, 0, 0, 300);
    parity("Sinc",            "Sinc[x]",           xyz, RRR, 1, -20.0, 20.0, 0, 0, 300);
    parity("InverseErf",      "InverseErf[x]",     xyz, RRR, 1, -0.98, 0.98, 0, 0, 300);
    parity("InverseErfc",     "InverseErfc[x]",    xyz, RRR, 1, 0.02, 1.98, 0, 0, 300);
    parity("expint in a chain",
           "Sin[SinIntegral[x]] + CosIntegral[x] Sinc[y] + ExpIntegralEi[-x]",
           xyz, RRR, 2, 0.4, 8.0, 0, 0, 300);

    /* Second kernel batch: Erfi, ProductLog, Fresnel, digamma, Zeta and the
     * real continuations of Fibonacci/LucasL.  Ranges stay inside each real
     * domain; the poles and branch points are covered by the decline list. */
    parity("Erfi",            "Erfi[x]",           xyz, RRR, 1, -6.0, 6.0, 0, 0, 300);
    parity("Erfi large",      "Erfi[x]",           xyz, RRR, 1, 6.0, 20.0, 0, 0, 200);
    parity("ProductLog",      "ProductLog[x]",     xyz, RRR, 1, -0.36, 40.0, 0, 0, 300);
    parity("ProductLog small","ProductLog[x]",     xyz, RRR, 1, 0.001, 1.0, 0, 0, 200);
    parity("FresnelC",        "FresnelC[x]",       xyz, RRR, 1, -12.0, 12.0, 0, 0, 400);
    parity("FresnelS",        "FresnelS[x]",       xyz, RRR, 1, -12.0, 12.0, 0, 0, 400);
    parity("PolyGamma",       "PolyGamma[x]",      xyz, RRR, 1, 0.2, 30.0, 0, 0, 300);
    parity("PolyGamma neg",   "PolyGamma[-x - 0.3]", xyz, RRR, 1, 0.1, 8.0, 0, 0, 200);
    parity("HarmonicNumber",  "HarmonicNumber[x]", xyz, RRR, 1, 0.1, 60.0, 0, 0, 300);
    parity("Zeta",            "Zeta[x]",           xyz, RRR, 1, 1.2, 30.0, 0, 0, 300);
    parity("Zeta strip",      "Zeta[x]",           xyz, RRR, 1, 0.55, 0.95, 0, 0, 200);
    parity("Zeta negative",   "Zeta[-x]",          xyz, RRR, 1, 0.2, 12.0, 0, 0, 300);
    parity("Fibonacci real",  "Fibonacci[x]",      xyz, RRR, 1, -6.0, 25.0, 0, 0, 300);
    parity("LucasL real",     "LucasL[x]",         xyz, RRR, 1, -6.0, 25.0, 0, 0, 300);
    parity("sf batch chain",  "Erfi[x]/100 + FresnelC[y] + PolyGamma[x] Zeta[y + 2]",
           xyz, RRR, 2, 0.4, 5.0, 0, 0, 300);
    /* PolyGamma matters as a BINARY kernel even written unary: the evaluator
     * canonicalises PolyGamma[x] to PolyGamma[0, x] before the compiler sees it. */
    parity("PolyGamma[n,x]",  "PolyGamma[2, x]",   xyz, RRR, 1, 0.3, 20.0, 0, 0, 300);
    parity("PolyGamma[5,x]",  "PolyGamma[5, x]",   xyz, RRR, 1, 0.5, 12.0, 0, 0, 200);
    parity("HurwitzZeta",     "HurwitzZeta[x + 1.2, y]", xyz, RRR, 2, 0.3, 9.0, 0, 0, 300);
    parity("Pochhammer",      "Pochhammer[x, 3]",  xyz, RRR, 1, 0.3, 9.0, 0, 0, 300);
    parity("Pochhammer real", "Pochhammer[x, 2.5]", xyz, RRR, 1, 0.4, 9.0, 0, 0, 300);
    parity("Binomial",        "Binomial[x, 2]",    xyz, RRR, 1, 0.3, 12.0, 0, 0, 300);
    parity("Binomial real",   "Binomial[x, 2.5]",  xyz, RRR, 1, 3.0, 12.0, 0, 0, 300);
    parity("LegendreP",       "LegendreP[3, x]",   xyz, RRR, 1, -0.99, 0.99, 0, 0, 300);
    parity("LegendreP high",  "LegendreP[12, x]",  xyz, RRR, 1, -0.9, 0.9, 0, 0, 300);
    /* Airy: three methods, not two.  The ascending series and the asymptotic
     * expansion do not meet in double precision, and the band between them is
     * covered by Taylor marching of y'' = x y from whichever end is exact.  The
     * band ranges below are the point of these tests — they used to be a decline. */
    parity("AiryAi small",    "AiryAi[x]",         xyz, RRR, 1, -2.5, 2.5, 0, 0, 300);
    parity("AiryBi small",    "AiryBi[x]",         xyz, RRR, 1, -2.5, 2.5, 0, 0, 300);
    parity("AiryAiPrime",     "AiryAiPrime[x]",    xyz, RRR, 1, -2.5, 2.5, 0, 0, 300);
    parity("AiryBiPrime",     "AiryBiPrime[x]",    xyz, RRR, 1, -2.5, 2.5, 0, 0, 300);
    parity("AiryAi large",    "AiryAi[x]",         xyz, RRR, 1, 8.0, 30.0, 0, 0, 200);
    parity("AiryAi large -",  "AiryAi[-x]",        xyz, RRR, 1, 8.0, 30.0, 0, 0, 200);
    parity("AiryBi large",    "AiryBi[x]",         xyz, RRR, 1, 8.0, 40.0, 0, 0, 200);
    /* The formerly-uncovered band, both signs.  Ai and Bi are marched in
     * OPPOSITE directions here (each in the one where it dominates), so a test
     * that only walked x > 0 would miss half the mechanism. */
    parity("AiryAi band",     "AiryAi[x]",         xyz, RRR, 1, 2.5, 8.0, 0, 0, 300);
    parity("AiryAi band -",   "AiryAi[-x]",        xyz, RRR, 1, 2.5, 8.0, 0, 0, 300);
    parity("AiryBi band",     "AiryBi[x]",         xyz, RRR, 1, 2.5, 8.0, 0, 0, 300);
    parity("AiryBi band -",   "AiryBi[-x]",        xyz, RRR, 1, 2.5, 8.0, 0, 0, 300);
    parity("AiryAiPrime band", "AiryAiPrime[x]",   xyz, RRR, 1, 2.5, 8.0, 0, 0, 300);
    parity("AiryAiPrime band -", "AiryAiPrime[-x]", xyz, RRR, 1, 2.5, 8.0, 0, 0, 300);
    parity("AiryBiPrime band", "AiryBiPrime[x]",   xyz, RRR, 1, 2.5, 8.0, 0, 0, 300);
    parity("AiryBiPrime band -", "AiryBiPrime[-x]", xyz, RRR, 1, 2.5, 8.0, 0, 0, 300);
    /* Across both seams, where a discontinuity between methods would show. */
    parity("AiryAi all",      "AiryAi[x]",         xyz, RRR, 1, -12.0, 12.0, 0, 0, 400);
    parity("AiryBi all",      "AiryBi[x]",         xyz, RRR, 1, -12.0, 9.0, 0, 0, 400);

    /* Bessel I/K, PolyLog, QPochhammer and the hypergeometrics.  BesselK needed
     * a continued fraction for x > 2: both small-x forms compute a decaying K
     * from quantities growing like e^x and lose ~2x/ln(10) digits. */
    parity("BesselI int order",  "BesselI[2, x]",   xyz, RRR, 1, 0.05, 18.0, 0, 0, 300);
    parity("BesselI real order", "BesselI[2.5, x]", xyz, RRR, 1, 0.05, 18.0, 0, 0, 300);
    parity("BesselI order 0",    "BesselI[0, x]",   xyz, RRR, 1, 0.05, 25.0, 0, 0, 300);
    parity("BesselK int order",  "BesselK[2, x]",   xyz, RRR, 1, 0.05, 40.0, 0, 0, 400);
    parity("BesselK real order", "BesselK[2.5, x]", xyz, RRR, 1, 0.05, 40.0, 0, 0, 400);
    parity("BesselK order 0",    "BesselK[0, x]",   xyz, RRR, 1, 0.05, 40.0, 0, 0, 400);
    parity("BesselK high order", "BesselK[7, x]",   xyz, RRR, 1, 0.5, 40.0, 0, 0, 300);
    parity("PolyLog 2",       "PolyLog[2, x]",     xyz, RRR, 1, -0.99, 0.99, 0, 0, 400);
    parity("PolyLog 3",       "PolyLog[3, x]",     xyz, RRR, 1, -0.99, 0.99, 0, 0, 300);
    parity("PolyLog real s",  "PolyLog[2.5, x]",   xyz, RRR, 1, -0.95, 0.95, 0, 0, 300);
    parity("QPochhammer",     "QPochhammer[x, 0.3]", xyz, RRR, 1, -2.0, 2.0, 0, 0, 300);
    parity("QPochhammer q",   "QPochhammer[0.5, x]", xyz, RRR, 1, -0.95, 0.95, 0, 0, 300);
    parity("Hypergeometric0F1", "Hypergeometric0F1[2., x]", xyz, RRR, 1, -20.0, 20.0, 0, 0, 300);
    parity("LerchPhi",        "LerchPhi[x, 2., 1.]", xyz, RRR, 1, -0.95, 0.95, 0, 0, 300);
    parity("LerchPhi s",      "LerchPhi[0.4, x, 2.]", xyz, RRR, 1, 0.5, 6.0, 0, 0, 300);
    parity("Hypergeometric1F1", "Hypergeometric1F1[1., 2., x]", xyz, RRR, 1, -40.0, 20.0, 0, 0, 300);
    parity("Hypergeometric1F1 b", "Hypergeometric1F1[0.5, x, 1.5]", xyz, RRR, 1, 0.3, 9.0, 0, 0, 300);
    parity("Hypergeometric2F1", "Hypergeometric2F1[1., 2., 3., x]", xyz, RRR, 1, -0.95, 0.95, 0, 0, 300);
    parity("Zeta at 0",       "Zeta[x - 1.]",      xyz, RRR, 1, 0.5, 1.5, 0, 0, 200);
    parity("sf batch 3 chain",
           "BesselK[1, x] + PolyLog[2, y/10.] Hypergeometric0F1[2., x] + BesselI[0, y]",
           xyz, RRR, 2, 0.4, 6.0, 0, 0, 300);

    /* UnitStep / Clip / Rescale are lowered by hand, not registered as kernels,
     * because their result TYPE is the difficulty: UnitStep must stay an
     * Integer, and Clip/Rescale take a list of bounds. */
    parity("UnitStep",        "UnitStep[x]",       xyz, RRR, 1, -3.0, 3.0, 0, 0, 300);
    parity("UnitStep 2-arg",  "UnitStep[x, y]",    xyz, RRR, 2, -3.0, 3.0, 0, 0, 300);
    parity("UnitStep in sum", "UnitStep[x - 1.] x + UnitStep[1. - x] x^2",
           xyz, RRR, 1, -3.0, 3.0, 0, 0, 300);
    parity("Clip",            "Clip[x, {1., 3.}]", xyz, RRR, 1, -5.0, 8.0, 0, 0, 300);
    parity("Rescale",         "Rescale[x, {0., 4.}]", xyz, RRR, 1, -5.0, 9.0, 0, 0, 300);
    parity("Clip+Rescale",    "Rescale[x, {1., 5.}] + Clip[x, {2., 4.}]",
           xyz, RRR, 1, -5.0, 9.0, 0, 0, 300);

    /* Out of domain, the machine kernel must DECLINE so the caller falls back —
     * these are exactly the arguments where the interpreter leaves the real
     * axis (Ci, Chi and li of a negative) or hits a pole (Ei at 0, li at 1).
     * Returning a real number here would be the compiled path answering
     * something the interpreter does not. */
    {
        static const struct { const char* body; double x; } DECLINE[] = {
            { "CosIntegral[x]",   -1.0 }, { "CosIntegral[x]",    0.0 },
            { "CoshIntegral[x]",  -2.0 }, { "CoshIntegral[x]",   0.0 },
            { "LogIntegral[x]",   -1.0 }, { "LogIntegral[x]",    1.0 },
            { "LogIntegral[x]",    0.0 }, { "ExpIntegralEi[x]",  0.0 },
            { "InverseErf[x]",     1.0 }, { "InverseErf[x]",    -1.0 },
            { "InverseErfc[x]",    0.0 }, { "InverseErfc[x]",    2.0 },
            { "ProductLog[x]",    -0.5 },          /* below the branch point */
            { "Zeta[x]",           1.0 },          /* pole */
            { "PolyGamma[x]",      0.0 }, { "PolyGamma[x]", -3.0 },   /* poles */
            { "BesselK[1, x]",    -1.0 },          /* complex for x < 0 */
            { "PolyLog[2, x]",     1.5 },          /* past the branch cut */
            { "QPochhammer[0.5, x]", 1.5 },        /* |q| >= 1: no convergence */
            { "Hypergeometric2F1[1., 2., 3., x]", 1.5 },  /* outside the disc */
        };
        int bad = 0, n = (int)(sizeof DECLINE / sizeof DECLINE[0]);
        for (int i = 0; i < n; i++) {
            Expr* b = parse_expression(DECLINE[i].body);
            CompiledProgram* p = compile_expr(b, x1i, RRR, 1);
            if (!p) { bad++; expr_free(b); continue; }
            double xv = DECLINE[i].x, got;
            if (compiled_eval_real(p, &xv, &got)) {
                printf("      %s at x=%g returned %g instead of declining\n",
                       DECLINE[i].body, xv, got);
                bad++;
            }
            compiled_free(p); expr_free(b);
        }
        if (bad) { printf("FAIL: %-30s %d/%d did not decline\n", "expint out-of-domain", bad, n); failures++; }
        else printf("ok:   %-30s %d arguments decline to the interpreter\n",
                    "expint out-of-domain", n);
    }
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
    /* Complex Sign, FractionalPart and Rescale.  These bailed for _Complex
     * until the audit was given a complex column — Sign's kernel was already in
     * the registry and the compiler's own inline lowering was shadowing it. */
    parity("complex Sign",    "Sign[x]",               x1i, CCC, 1, 0.3, 2.0, 0, 0, 300);
    parity("complex Sign sum", "Sign[x] + Sign[y]",    xy,  CCC, 2, -2.0, 2.0, 0, 0, 300);
    parity("complex FracPart", "FractionalPart[x]",    x1i, CCC, 1, -4.0, 4.0, 0, 0, 300);
    parity("complex Rescale", "Rescale[x, {0., 4.}]",  x1i, CCC, 1, -3.0, 5.0, 0, 0, 300);
    parity("complex Rescale 2", "Rescale[x, {1., 5.}] + Sign[y]", xy, CCC, 2, -3.0, 5.0, 0, 0, 300);
    /* Complex Gamma and LogGamma, sharing the interpreter's own Lanczos series.
     * Both half-planes and both sides of the Re = 1/2 reflection boundary: the
     * reflection is where the branch structure lives, and LogGamma's continued
     * branch there was wrong (principal, so short by a multiple of 2 pi i) until
     * this work.  `x - 3.` shifts a positive-quadrant sample into the upper half
     * plane with mixed real part; Conjugate puts the mirror image in the lower. */
    parity("complex Gamma",       "Gamma[x]",                 x1i, CCC, 1, 0.3, 3.0, 0, 0, 300);
    parity("complex Gamma refl",  "Gamma[x - 3.]",            x1i, CCC, 1, 0.2, 4.0, 0, 0, 300);
    parity("complex Gamma refl-", "Gamma[Conjugate[x] - 3.]", x1i, CCC, 1, 0.2, 4.0, 0, 0, 300);
    parity("complex LogGamma",      "LogGamma[x]",                 x1i, CCC, 1, 0.3, 3.0, 0, 0, 300);
    parity("complex LogGamma refl", "LogGamma[x - 3.]",            x1i, CCC, 1, 0.2, 4.0, 0, 0, 300);
    parity("complex LogGamma refl-","LogGamma[Conjugate[x] - 3.]", x1i, CCC, 1, 0.2, 4.0, 0, 0, 300);
    parity("complex LogGamma far",  "LogGamma[x + 12. I]",         x1i, CCC, 1, 0.3, 3.0, 0, 0, 300);
    parity("complex Gamma chain", "Gamma[x] Gamma[y] / Gamma[x + y]", xy, CCC, 2, 0.3, 2.5, 0, 0, 300);
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
    must_bail("no kernel (BarnesG)", "BarnesG[x]", x1, RRR, 1); /* no machine kernel -> bail */
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
        for (int pass = 0; pass < 2; pass++) {
        ab_off_flag = pass ? COMPILE_NO_CSE : COMPILE_NO_OPT;
        ab_what     = pass ? "CSE" : "the optimiser";

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
        /* Bodies with genuinely repeated subtrees — what Expr-level CSE exists
         * for, and the shapes where hoisting could go wrong. */
        AB("repeat: same call thrice", "Sin[x y] + Cos[x y] Sin[x y] + Sin[x y]^3", 2, RRR, false);
        AB("repeat: nested repeats",   "Exp[Sin[x] + Sin[x]] + Sin[x] (Sin[x] + 1)", 1, RRR, false);
        AB("repeat: deep shared",      "Sqrt[Abs[x y]] + Log[1 + Sqrt[Abs[x y]]] + Sqrt[Abs[x y]]^2", 2, RRR, false);
        AB("repeat: under a loop",     "Sum[Sin[x] Cos[y] + Sin[x] i, {i, 1, 12}]", 2, RRR, true);
        AB("repeat: inside If",        "If[x < y, Sin[x y] + Sin[x y], Cos[x y] - Cos[x y]^2]", 2, RRR, false);
        AB("repeat: with locals",      "With[{a = Sin[x y]}, a + Sin[x y] + a Sin[x y]]", 2, RRR, true);
        AB("repeat: complex",          "(x + I y)^2 + Exp[(x + I y)^2]", 2, CCC, false);
        }

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

    /* ================= ANY RANK =================
     * The delegated NDArray path was already rank-general; the compiler's own
     * rank-1 front gate was the only thing keeping matrices and higher tensors
     * out.  With it lifted, an elementwise body works at any rank the packed
     * type encoding can name.  (`Total` deliberately stays rank-1: at rank >= 2
     * it reduces only the LEADING axis, which is a different operation.) */
    {
        const char* vw[] = { "v", "w" };
        const int64_t d2[2] = { 7, 5 };
        const int64_t d3[3] = { 4, 3, 5 };
        const int64_t d4[4] = { 3, 2, 4, 3 };
        parity_nd("rank-2 elementwise",  "v + 2 w",              vw, 2, 2, d2, COMPILE_NO_FUSE, 20);
        parity_nd("rank-2 libm chain",   "Sin[v] Exp[-v] + Sqrt[v]", vw, 1, 2, d2, COMPILE_NO_FUSE, 20);
        parity_nd("rank-2 power",        "v^3 + w^2",            vw, 2, 2, d2, COMPILE_NO_FUSE, 20);
        parity_nd("rank-3 elementwise",  "Sin[v w] + Exp[-v]",   vw, 2, 3, d3, COMPILE_NO_FUSE, 20);
        parity_nd("rank-3 Gamma",        "Gamma[v] + Log[2, w]", vw, 2, 3, d3, COMPILE_NO_FUSE, 20);
        parity_nd("rank-4 elementwise",  "v w - v / w",          vw, 2, 4, d4, COMPILE_NO_FUSE, 20);
    }

    /* ================= ELEMENTWISE FUSION (opt-in) =================
     * Fusion strip-mines an elementwise chain into ONE pass over the buffers
     * rather than one NDArray pass (and one temporary buffer) per operation.  It
     * is ON by default (COMPILE_NO_FUSE disables it), so these run the fused
     * path; the block above runs the same shapes delegated. */
    {
        const char* vw[] = { "v", "w" };
        const CompileType AA[] = { CT_ARRAY(CT_REAL, 1), CT_ARRAY(CT_REAL, 1) };
        const CompileType CA[] = { CT_ARRAY(CT_COMPLEX, 1), CT_ARRAY(CT_COMPLEX, 1) };
        const int64_t d2[2] = { 7, 5 };
        const int64_t d3[3] = { 4, 3, 5 };
        parity_arr_ex("fuse: vec chain",  "Sin[v] Exp[-v] + Sqrt[v]", vw, AA, 1, 64, 0.3, 4.0, 20, 0u);
        parity_arr_ex("fuse: Total chain","Total[Sin[v] Exp[-v] + Sqrt[v]]", vw, AA, 1, 64, 0.3, 4.0, 20, 0u);
        parity_arr_ex("fuse: vec + vec",  "v w + v - w",   vw, AA, 2, 64, 0.3, 4.0, 20, 0u);
        parity_arr_ex("fuse: power",      "v^3 + 2 v + 1", vw, AA, 1, 64, 0.3, 4.0, 20, 0u);
        parity_arr_ex("fuse: Gamma",      "Gamma[v] + Erf[w]", vw, AA, 2, 64, 0.5, 3.0, 20, 0u);
        parity_arr_ex("fuse: complex vec","v w + Exp[v]",  vw, CA, 2, 32, 0.3, 2.0, 20, 0u);
        parity_nd    ("fuse: rank-2",     "Sin[v] + v w",  vw, 2, 2, d2, 0u, 20);
        parity_nd    ("fuse: rank-3",     "v^2 - Exp[-w]", vw, 2, 3, d3, 0u, 20);

        /* Non-Listable heads must NOT fuse: the interpreter does not thread them,
         * so an elementwise loop would quietly answer something different.
         * Max[v,w] is the largest single element, not an elementwise maximum. */
        must_bail_raw("fuse: Max stays scalar", "Max[v, w]", vw, AA, 2);
        must_bail_raw("fuse: If stays scalar",  "If[v > 0, v, -v]", vw, AA, 1);
        /* The interpreter leaves ArcTan[NDArray, NDArray] unevaluated because the
         * ND binary kernels are one-array-plus-scalar only.  Fusion could compute
         * it, and must not: answering where the interpreter declines is the same
         * class of bug as answering differently. */
        must_bail_raw("fuse: ArcTan[v,w] declines", "ArcTan[v, w]", vw, AA, 2);
        must_bail_raw("fuse: BesselJ[v,w] declines", "BesselJ[v, w]", vw, AA, 2);
    }

    /* ================= STRIP-MINING STRESS =================
     * Fusion processes VBLOCK elements per opcode, so the tail of the last tile
     * is the one thing a length-64 test can never exercise.  Every body is run
     * at lengths straddling the tile boundary — and at lengths shorter than one
     * tile, where the very first tile is already partial. */
    {
        static const size_t LENS[] = { 1, 2, 3, 7, 31, 63, 64, 65, 66, 127, 128, 129,
                                       255, 256, 257, 1000, 4097 };
        static const size_t NLEN = sizeof LENS / sizeof LENS[0];
        static const char* BODIES[] = {
            "v + 1", "2 v", "v w", "v + w", "v - w", "v / w", "v^2", "v^3",
            "v^2 + 2 v + 1", "Sqrt[v]", "Sin[v]", "Exp[-v] + Log[v]",
            "Sin[v w] + Exp[-v] (w + 2)", "Gamma[v] + Erf[w]",
            "Total[v]", "Total[v w]", "Total[Sin[v] Exp[-v] + Sqrt[v]]",
            /* Log[b,x] lowers to the arithmetic identity Log[x]/Log[b], so it
             * threads array-by-array in both paths.  ArcTan does NOT (see
             * fuse_listable): the interpreter leaves ArcTan[nd, nd] unevaluated,
             * so it is checked as a bail below rather than for parity. */
            "Log[v, w]", "ArcTan[v, 2]",
            "Total[(v - Total[v]/Length[v])^2] / Length[v]",
        };
        static const size_t NBODY = sizeof BODIES / sizeof BODIES[0];
        const char* vw[] = { "v", "w" };
        const CompileType AA[] = { CT_ARRAY(CT_REAL, 1), CT_ARRAY(CT_REAL, 1) };
        int checked = 0, bad = 0;
        for (size_t bi = 0; bi < NBODY; bi++) {
            for (size_t li = 0; li < NLEN; li++) {
                int before = failures;
                parity_arr("(stress)", BODIES[bi], vw, AA, 2, LENS[li], 0.4, 3.0, 2);
                if (failures > before) {
                    printf("      ^ body \"%s\" at length %zu\n", BODIES[bi], LENS[li]);
                    bad++;
                }
                checked++;
            }
        }
        /* parity_arr prints a line per call; the summary is what matters. */
        printf("ok:   %-30s %d body x length combinations, %d bad\n",
               "stress: tile boundaries", checked, bad);
    }

    /* Fused and delegated must agree.  Not bitwise — the delegated reduction
     * sums pairwise while the fused one accumulates a tile at a time — but far
     * inside any tolerance that could hide a real disagreement. */
    {
        const char* vw[] = { "v", "w" };
        const CompileType AA[] = { CT_ARRAY(CT_REAL, 1), CT_ARRAY(CT_REAL, 1) };
        static const char* BODIES[] = {
            "v w + v - w", "v^3 + 2 v + 1", "Sin[v] Exp[-v] + Sqrt[v]",
            "Gamma[v] + Log[2, w]", "Abs[v - w] + ArcTan[v, w]",
            "Total[v w]", "Total[Sqrt[v] + w^2]",
        };
        int bad = 0, n = 0;
        for (size_t i = 0; i < sizeof BODIES / sizeof BODIES[0]; i++) {
            for (size_t len = 63; len <= 129; len += 33) {
                if (!fused_matches_delegated(BODIES[i], vw, AA, 2, len)) {
                    printf("FAIL: fused != delegated for \"%s\" at length %zu\n", BODIES[i], len);
                    failures++; bad++;
                }
                n++;
            }
        }
        if (!bad) printf("ok:   %-30s %d bodies agree with the delegated path\n",
                         "stress: fused vs delegated", n);
    }

    /* Randomised elementwise trees.  Only Listable heads are generated, because
     * that is exactly the set fusion is allowed to thread over; the point is to
     * cover shapes the hand-written bodies do not. */
    {
        const char* vw[] = { "v", "w" };
        const CompileType AA[] = { CT_ARRAY(CT_REAL, 1), CT_ARRAY(CT_REAL, 1) };
        int bad = 0, ncompiled = 0;
        allow_decline = true;
        for (int t = 0; t < 120; t++) {
            char buf[512];
            rand_arr_body(buf, sizeof buf, 3);
            size_t len = (size_t)irand(1, 200);
            int before = failures;
            parity_arr("(random)", buf, vw, AA, 2, len, 0.4, 3.0, 1);
            if (failures > before) { printf("      ^ random body \"%s\" len %zu\n", buf, len); bad++; }
            else ncompiled++;
        }
        allow_decline = false;
        printf("ok:   %-30s %d random trees, %d bad\n",
               "stress: random array bodies", ncompiled, bad);
    }

    /* The frame and the tile bank are reused across calls, so a stale handle or
     * a tile left holding the previous call's data would only show up on the
     * SECOND call.  Same program, many calls, alternating lengths. */
    {
        const char* inm[1] = { intern_symbol("v") };
        const CompileType AT[1] = { CT_ARRAY(CT_REAL, 1) };
        Expr* b = parse_expression("Total[Sqrt[v] + v^2]");
        CompiledProgram* p = compile_expr(b, inm, AT, 1);
        int bad = 0;
        if (!p) { printf("FAIL: reuse body did not compile\n"); failures++; }
        else {
            for (int it = 0; it < 200; it++) {
                size_t len = (size_t)(1 + (it * 37) % 300);
                Expr* v = make_vec(len, 0.4, 3.0);
                CompileValue av, out;
                av.type = AT[0]; av.v.a = v;
                double want = 0;
                const double* raw = (const double*)v->data.ndarray.data;
                for (size_t k = 0; k < len; k++) want += sqrt(raw[k]) + raw[k] * raw[k];
                if (!compiled_eval(p, &av, &out) || fabs(out.v.r - want) > 1e-9 * (fabs(want) + 1))
                    bad++;
                expr_free(v);
            }
            if (bad) { printf("FAIL: %d/200 repeated calls wrong (frame/tile reuse)\n", bad); failures++; }
            else printf("ok:   %-30s 200 calls, varying lengths\n", "stress: frame + tile reuse");
        }
        compiled_free(p); expr_free(b);
    }

    /* ================= RESULT HEAD, NOT JUST RESULT VALUE =================
     * Every head whose result TYPE differs from its argument type, checked
     * against the interpreter.  See parity_head on why the numeric tests miss
     * this entirely. */
    parity_head("Sign real +",      "Sign[xq]",              2.5);
    parity_head("Sign real -",      "Sign[xq]",             -2.5);
    parity_head("Sign real 0",      "Sign[xq]",              0.0);
    parity_head("Floor",            "Floor[xq]",             2.5);
    parity_head("Ceiling",          "Ceiling[xq]",           2.5);
    parity_head("Round",            "Round[xq]",             2.5);
    parity_head("UnitStep",         "UnitStep[xq]",          0.5);
    parity_head("IntegerPart",      "IntegerPart[xq]",       2.5);
    parity_head("IntegerPart -",    "IntegerPart[xq]",      -2.5);
    parity_head("Floor -",          "Floor[xq]",            -2.5);
    parity_head("FractionalPart",   "FractionalPart[xq]",    2.5);
    parity_head("FractionalPart -", "FractionalPart[xq]",   -2.5);
    parity_head("Quotient",         "Quotient[xq, 3.]",      5.5);
    parity_head("Abs",              "Abs[xq]",              -2.5);
    parity_head("Clip",             "Clip[xq, {1., 3.}]",    2.0);
    parity_head("Rescale",          "Rescale[xq, {0., 4.}]", 2.0);
    parity_head("Max",              "Max[xq, 1.]",           2.5);
    parity_head("Sqrt",             "Sqrt[xq]",              2.5);

    /* ================= THREADED FUSED MAP (OP_APAR) =================
     * A fused MAP is split across threads above NDARRAY_THREAD_THRESHOLD.  Each
     * output element is computed by the same operations on the same inputs
     * whichever core runs it, so the threaded and single-threaded programs must
     * agree BITWISE — not "to 1e-12".  memcmp is the assertion precisely because
     * a race, a shared tile buffer or an off-by-one chunk boundary would show up
     * as a handful of wrong elements that any tolerance-based check would pass.
     *
     * The lengths straddle the fan-out threshold and are deliberately NOT
     * multiples of VBLOCK, so the last tile of every chunk is short. */
    {
        const char* inm[2] = { intern_symbol("v"), intern_symbol("w") };
        const CompileType AT[2] = { CT_ARRAY(CT_REAL, 1), CT_ARRAY(CT_REAL, 1) };
        static const char* BODIES[] = {
            "Sqrt[v] + v^2",
            "Sin[v] Exp[-v] + Sqrt[v]",
            "v w + Log[v] Cos[w]",
            "Gamma[v] + Erf[w]",              /* via the machine-kernel path */
        };
        static const size_t LENS[] = { 99999, 100000, 100001, 262144, 1000003 };
        int bad = 0, ran = 0, threaded_seen = 0;
        for (size_t bi = 0; bi < sizeof BODIES / sizeof BODIES[0]; bi++) {
            Expr* b = parse_expression(BODIES[bi]);
            CompiledProgram* pp = compile_expr_ex(b, inm, AT, 2, 0u);
            CompiledProgram* ps = compile_expr_ex(b, inm, AT, 2, COMPILE_NO_PAR);
            if (!pp || !ps) { printf("FAIL: par A/B body %zu did not compile\n", bi); failures++; }
            else {
                /* The two programs must differ by exactly the APAR marker; if
                 * they are identical the fan-out never got emitted and this
                 * whole section is silently testing nothing. */
                if (compiled_num_instructions(pp) <= compiled_num_instructions(ps)) bad++;
                else threaded_seen++;
                for (size_t li = 0; li < sizeof LENS / sizeof LENS[0]; li++) {
                    size_t len = LENS[li];
                    Expr* v = make_vec(len, 0.4, 3.0);
                    Expr* w = make_vec(len, 0.4, 3.0);
                    CompileValue args[2], op_, os_;
                    args[0].type = AT[0]; args[0].v.a = v;
                    args[1].type = AT[1]; args[1].v.a = w;
                    bool sp = compiled_eval(pp, args, &op_);
                    bool ss = compiled_eval(ps, args, &os_);
                    ran++;
                    if (sp != ss) bad++;
                    else if (sp) {
                        Expr* ap = aval_to_expr(op_);
                        Expr* as = aval_to_expr(os_);
                        if (ndarray_size(ap) != ndarray_size(as)
                            || memcmp(ap->data.ndarray.data, as->data.ndarray.data,
                                      len * sizeof(double)) != 0)
                            bad++;
                        expr_free(ap); expr_free(as);
                    }
                    expr_free(v); expr_free(w);
                }
            }
            compiled_free(pp); compiled_free(ps); expr_free(b);
        }
        if (bad) { printf("FAIL: threaded fused map: %d/%d mismatched\n", bad, ran); failures++; }
        else if (!threaded_seen) { printf("FAIL: threaded fused map: fan-out never emitted\n"); failures++; }
        else printf("ok:   %-30s %d runs bitwise identical\n", "threaded fused map (APAR)", ran);
    }

    /* ================= COMPILED -> COMPILED CALLS =================
     * A CompiledFunction callee is INLINED up to a depth cap, beyond which the
     * whole body used to bail — a chain of eleven compiled functions dropped
     * entirely to the interpreter.  OP_CALL is the fallback: the callee runs on
     * its own frame with machine values passed in registers, no Expr and no
     * evaluator round-trip, so the chain compiles instead. */
    {
        /* ELEVEN nested applications of one compiled callee.  The inliner
         * pastes the first eight in; past its depth cap the rest used to bail
         * the entire body, and are now CALLed.
         *
         * (A chain of DISTINCT Compile[] objects would not work here, by design:
         * user Compile[] does not fold globals, because the object outlives its
         * defining scope — so gc2's body cannot resolve gc1 at the time gc2 is
         * compiled.  Nesting one already-compiled callee is the shape that
         * actually arises.) */
        const int DEPTH = 11;
        /* The callee is deliberately LARGER than INLINE_MAX_INSTRS, so the
         * compiler CALLs it rather than pasting it in — which is what exercises
         * OP_CALL.  Same body twice: once as a Compile[] object, once as a
         * DownValue, which never compiles and is therefore an honest reference. */
        /* `Set` evaluates to the assigned value, so the result must be freed —
         * eval_and_free only consumes its ARGUMENT. */
        expr_free(eval_and_free(parse_expression(
            "gcf = Compile[{x}, Sin[x]/2 + Cos[x]/3 + Sqrt[Abs[x]]/5 + Exp[-x x]/7"
            " + Log[1 + x x]/11 + Tanh[x]/13 + ArcTan[x]/17 + x/19]")));
        expr_free(eval_and_free(parse_expression(
            "icf[a_] := Sin[a]/2 + Cos[a]/3 + Sqrt[Abs[a]]/5 + Exp[-a a]/7"
            " + Log[1 + a a]/11 + Tanh[a]/13 + ArcTan[a]/17 + a/19")));

        char call[512], ref0[512];
        { size_t q = 0, r = 0;
          for (int k = 0; k < DEPTH; k++) q += (size_t)snprintf(call + q, sizeof call - q, "gcf[");
          q += (size_t)snprintf(call + q, sizeof call - q, "x");
          for (int k = 0; k < DEPTH; k++) q += (size_t)snprintf(call + q, sizeof call - q, "]");
          for (int k = 0; k < DEPTH; k++) r += (size_t)snprintf(ref0 + r, sizeof ref0 - r, "icf[");
          r += (size_t)snprintf(ref0 + r, sizeof ref0 - r, "xq");
          for (int k = 0; k < DEPTH; k++) r += (size_t)snprintf(ref0 + r, sizeof ref0 - r, "]");
        }

        const char* inm[1] = { intern_symbol("x") };
        const CompileType RR1[1] = { CT_REAL };
        Expr* b = parse_expression(call);
        CompiledProgram* p = compile_expr_ex(b, inm, RR1, 1, COMPILE_FOLD_GLOBALS);
        if (!p) {
            printf("FAIL: %-30s depth-%d nesting did not compile\n", "compiled->compiled call", DEPTH);
            failures++;
        } else {
            int bad = 0, nofin = 0, noref = 0; double maxerr = 0;
            for (int t = 0; t < 100; t++) {
                double xv = urand(0.2, 3.0), got;
                if (!compiled_eval_real(p, &xv, &got)) { bad++; nofin++; continue; }
                Expr* r = ref_at(ref0, xv);
                double want;
                if (!expr_to_double(r, &want)) { bad++; noref++; }
                else {
                    double e = fabs(got - want) / (fabs(want) + 1e-30);
                    if (e > maxerr) { maxerr = e;
                        if (e > 1e-12) printf("      x=%.17g got=%.17g want=%.17g\n", xv, got, want); }
                }
                expr_free(r);
            }
            if (bad || maxerr > 1e-12) {
                printf("FAIL: %-30s %d bad (%d eval-fail, %d ref-fail), max_rel=%.2e\n",
                       "compiled->compiled call", bad, nofin, noref, maxerr);
                failures++;
            } else printf("ok:   %-30s depth %d, max_rel=%.1e (100 pts)\n",
                          "compiled->compiled call", DEPTH, maxerr);
            compiled_free(p);
        }
        expr_free(b);

        /* The callee gets its OWN frame.  Inside a loop the caller's frame is
         * live across every re-entry, so a shared frame would corrupt the
         * accumulator — the test the old per-program frame could not have passed. */
        {
            char loop[700];
            snprintf(loop, sizeof loop,
                     "Module[{s = 0.}, Do[s = s + %s, {i, 1, 8}]; s]", call);
            Expr* b2 = parse_expression(loop);
            CompiledProgram* p2 = compile_expr_ex(b2, inm, RR1, 1, COMPILE_FOLD_GLOBALS);
            if (!p2) { printf("FAIL: call inside a loop did not compile\n"); failures++; }
            else {
                int bad = 0; double maxerr = 0;
                for (int t = 0; t < 40; t++) {
                    double xv = urand(0.2, 2.0), got;
                    if (!compiled_eval_real(p2, &xv, &got)) { bad++; continue; }
                    Expr* r = ref_at(ref0, xv);
                    double v = 0;
                    if (!expr_to_double(r, &v)) bad++;
                    expr_free(r);
                    double want = 8.0 * v;           /* body does not depend on i */
                    double e = fabs(got - want) / (fabs(want) + 1e-30);
                    if (e > maxerr) maxerr = e;
                }
                if (bad || maxerr > 1e-12) {
                    printf("FAIL: %-30s %d bad, max_rel=%.2e\n", "call inside a loop", bad, maxerr);
                    failures++;
                } else printf("ok:   %-30s max_rel=%.1e (40 pts x 8 iters)\n",
                              "call inside a loop", maxerr);
                compiled_free(p2);
            }
            expr_free(b2);
        }
    }

    if (failures == 0) printf("\nAll Compile engine tests passed.\n");
    else printf("\n%d Compile engine test(s) FAILED.\n", failures);
    return failures ? 1 : 0;
}
