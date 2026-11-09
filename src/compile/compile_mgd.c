/* Mathilda — Compile[]: the arbitrary-precision managed-scalar emitter.
 *
 * A self-contained lowering for arbitrary-precision scalar bodies (MPFR reals /
 * complex under WorkingPrecision -> n, GMP integers under "BigIntegers" ->
 * True), kept ENTIRELY SEPARATE from the machine emit_node so the machine path
 * is provably unchanged.  A managed register's Slot holds a POINTER to its
 * per-call container (mpfr_t / mpz_t / ncpx), bound once at frame entry; the
 * VM-side arena and frame entry/exit live in compile_vm.c.  Managed programs
 * force COMPILE_NO_OPT | COMPILE_NO_CSE so the optimiser never sees a managed
 * register.  See docs/design/compile-arbitrary-precision.md. */
#include "compile_emit.h"        /* Ctx / Val / ins / arg_find */
#include "../arithmetic.h"       /* is_rational, expr_new_* */
#include "../sym_names.h"         /* SYM_Plus, SYM_Times, ...: used in every build */
#ifdef USE_MPFR
#include "../numeric_complex.h"  /* ncpx, numeric_mpfr_* — the managed containers */
#include "../numeric.h"          /* numeric_mpfr_pow_fast */
#include "../special_functions/productlog.h"   /* productlog_mpfr_real */
#include "../symtab.h"           /* attributes: NumericFunction (constant folding) */
#include "../attr.h"
#include "../eval.h"
#endif
#include "compile_mgd_special.h"  /* MG_SPECIAL: the special-function registry */
#include <math.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>

/* ================================================================== *
 *  Arbitrary precision — the managed-scalar emitter (emit_mgd)        *
 * ================================================================== *
 * A self-contained lowering for arbitrary-precision scalar bodies, kept
 * ENTIRELY SEPARATE from the machine emit_node so the machine path is provably
 * unchanged (the machine emitter is not touched at all).  A program that reaches
 * here has prec_bits > 0 (MPFR reals/complex at that precision) or allow_bigint
 * (GMP integers), and the whole body is lowered in that one domain.  v1 scope:
 * straight-line arithmetic — literals, declared args, Plus/Times/Power/Subtract
 * and the unary elementary functions.  Anything else bails to the interpreter,
 * exactly like an unsupported head on the machine path.  Managed programs force
 * COMPILE_NO_OPT | COMPILE_NO_CSE, so the optimiser never sees a managed
 * register.  See docs/design/compile-arbitrary-precision.md. */

/* MGC_* complex-op selectors are declared in compile_internal.h — the VM reads
 * the same imm.i selector the emitter writes here. */

/* Allocate a fresh managed register of container type `t`.  Bump-only: a managed
 * register is never popped, so each physical managed slot has exactly ONE
 * container type for the whole program, which is what lets frame entry allocate
 * the right container by register number.  Managed bodies are scalar loops with
 * modest temp counts, so the wasted registers are cheap. */
static int alloc_mgd(Ctx* c, CompileType t) {
    if (c->mgd_top >= c->mgd_type_cap) {
        int nc = c->mgd_type_cap ? c->mgd_type_cap * 2 : 32;
        CompileType* nt = realloc(c->mgd_type, (size_t)nc * sizeof *nt);
        if (!nt) { c->ok = false; return MGD_VREG; }
        c->mgd_type = nt; c->mgd_type_cap = nc;
    }
    int idx = c->mgd_top;
    c->mgd_type[idx] = t;
    c->mgd_top++;
    if (c->mgd_top > c->mgd_max) c->mgd_max = c->mgd_top;
    return MGD_VREG + idx;
}

/* Free one program-owned literal container (used on the bail path and by
 * compiled_free). */
void mgd_const_free(MgdConst* mc) {
    if (!mc || !mc->ptr) return;
    switch (mc->kind) {
        case 1: mpz_clear((mpz_ptr)mc->ptr); break;
#ifdef USE_MPFR
        case 0: mpfr_clear((mpfr_ptr)mc->ptr); break;
        case 2: ncpx_clear((ncpx*)mc->ptr); break;
#endif
        case MGD_CONST_SPECIAL:   /* owns its literals; frees itself */
            mgd_special_call_free((MgdSpecialCall*)mc->ptr);
            mc->ptr = NULL;
            return;
        default: break;
    }
    free(mc->ptr);
    mc->ptr = NULL;
}

/* Record a program-owned literal container; its stable pointer rides in an
 * instruction immediate and ownership transfers to the program at finalize. */
static bool ctx_add_mgd_const(Ctx* c, int kind, void* ptr) {
    if (c->nmgd_consts == c->mgd_consts_cap) {
        int nc = c->mgd_consts_cap ? c->mgd_consts_cap * 2 : 8;
        MgdConst* nb = realloc(c->mgd_consts, (size_t)nc * sizeof *nb);
        if (!nb) { c->ok = false; return false; }
        c->mgd_consts = nb; c->mgd_consts_cap = nc;
    }
    c->mgd_consts[c->nmgd_consts].kind = kind;
    c->mgd_consts[c->nmgd_consts].ptr  = ptr;
    c->nmgd_consts++;
    return true;
}

/* An integer-valued exponent that fits a machine long (Power fast path). */
static bool mgd_small_int(const Expr* e, long* n) {
    if (e->type == EXPR_INTEGER) { *n = (long)e->data.integer; return true; }
    if (e->type == EXPR_BIGINT && mpz_fits_slong_p(e->data.bigint)) {
        *n = mpz_get_si(e->data.bigint); return true;
    }
    return false;
}

static bool emit_mgd_int(Ctx* c, const Expr* e, Val* out);

/* ---- GMP (bignum) constant materialisation ---- */
/* z = v for any int64 v.  mpz_set_si takes a `long`, which is 32 bits on LLP64
 * Windows, so a value outside it goes through mpz_import of its magnitude. */
static void mgd_mpz_set_i64(mpz_ptr z, int64_t v) {
    if (v >= LONG_MIN && v <= LONG_MAX) { mpz_set_si(z, (long)v); return; }
    uint64_t mag = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
    mpz_import(z, 1, 1, sizeof mag, 0, 0, &mag);
    if (v < 0) mpz_neg(z, z);
}
bool mgd_mpz_from_expr(mpz_ptr z, const Expr* e) {
    if (e->type == EXPR_INTEGER) { mgd_mpz_set_i64(z, e->data.integer); return true; }
    if (e->type == EXPR_BIGINT)  { mpz_set(z, e->data.bigint); return true; }
    return false;   /* Real / Rational are not exact integers */
}
static int mgd_new_mpz_const(Ctx* c, const Expr* lit) {
    mpz_ptr z = malloc(sizeof *z);
    if (!z) { c->ok = false; return -1; }
    mpz_init(z);
    if (!mgd_mpz_from_expr(z, lit) || !ctx_add_mgd_const(c, 1, z)) {
        mpz_clear(z); free(z); return -1;
    }
    Slot k; memset(&k, 0, sizeof k); k.p = z;
    int dst = alloc_mgd(c, CT_BIGINT);
    ins(c, OP_MG_ZCONST, (uint32_t)dst, 0, 0, k);
    return dst;
}

#ifdef USE_MPFR
/* Rational[p, q] with integer (machine or bignum) parts. */
static bool mgd_is_exact_rational(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2
        && expr_is_integer_like(e->data.function.args[0])
        && expr_is_integer_like(e->data.function.args[1]);
}

/* ---- MPFR (real) constant materialisation ---- */
bool mgd_mpfr_from_expr(mpfr_ptr m, const Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER:
            if (e->data.integer >= LONG_MIN && e->data.integer <= LONG_MAX) {
                mpfr_set_si(m, (long)e->data.integer, MPFR_RNDN);
            } else {               /* `long` is 32 bits on LLP64 Windows */
                mpz_t z;
                mpz_init(z);
                mgd_mpz_set_i64(z, e->data.integer);
                mpfr_set_z(m, z, MPFR_RNDN);
                mpz_clear(z);
            }
            return true;
        case EXPR_REAL:    mpfr_set_d (m, e->data.real,          MPFR_RNDN); return true;
        case EXPR_BIGINT:  mpfr_set_z (m, e->data.bigint,        MPFR_RNDN); return true;
        case EXPR_MPFR:    mpfr_set   (m, e->data.mpfr,          MPFR_RNDN); return true;
        default: break;
    }
    if (mgd_is_exact_rational(e)) {
        /* exact quotient rounded once (a bignum numerator or denominator --
         * (99/100)^11 -- is as much a constant as 1/3) */
        mpq_t q;
        mpq_init(q);
        mgd_mpz_from_expr(mpq_numref(q), e->data.function.args[0]);
        mgd_mpz_from_expr(mpq_denref(q), e->data.function.args[1]);
        bool ok = mpz_sgn(mpq_denref(q)) != 0;
        if (ok) { mpq_canonicalize(q); mpfr_set_q(m, q, MPFR_RNDN); }
        mpq_clear(q);
        return ok;
    }
    return false;
}
/* A named real constant (Pi, E, ...) at working precision. */
static bool mgd_mpfr_from_const(mpfr_ptr m, const char* nm) {
    if (nm == SYM_Pi)         { mpfr_const_pi(m, MPFR_RNDN); return true; }
    if (nm == SYM_E)          { mpfr_set_ui(m, 1, MPFR_RNDN); mpfr_exp(m, m, MPFR_RNDN); return true; }
    if (nm == SYM_EulerGamma) { mpfr_const_euler(m, MPFR_RNDN); return true; }
    if (nm == SYM_Catalan)    { mpfr_const_catalan(m, MPFR_RNDN); return true; }
    if (nm == SYM_Degree)     { mpfr_const_pi(m, MPFR_RNDN); mpfr_div_ui(m, m, 180, MPFR_RNDN); return true; }
    return false;
}
static int mgd_new_mpfr_const_expr(Ctx* c, const Expr* lit) {
    mpfr_ptr m = malloc(sizeof *m);
    if (!m) { c->ok = false; return -1; }
    mpfr_init2(m, c->prec_bits);
    if (!mgd_mpfr_from_expr(m, lit) || !ctx_add_mgd_const(c, 0, m)) {
        mpfr_clear(m); free(m); return -1;
    }
    Slot k; memset(&k, 0, sizeof k); k.p = m;
    int dst = alloc_mgd(c, CT_BIGREAL);
    ins(c, OP_MG_RCONST, (uint32_t)dst, 0, 0, k);
    return dst;
}
static int mgd_new_mpfr_const_named(Ctx* c, const char* nm) {
    mpfr_ptr m = malloc(sizeof *m);
    if (!m) { c->ok = false; return -1; }
    mpfr_init2(m, c->prec_bits);
    if (!mgd_mpfr_from_const(m, nm) || !ctx_add_mgd_const(c, 0, m)) {
        mpfr_clear(m); free(m); return -1;
    }
    Slot k; memset(&k, 0, sizeof k); k.p = m;
    int dst = alloc_mgd(c, CT_BIGREAL);
    ins(c, OP_MG_RCONST, (uint32_t)dst, 0, 0, k);
    return dst;
}
/* ---- ncpx (complex) constant materialisation ---- */
static int mgd_new_ncpx_const(Ctx* c, const Expr* re, const Expr* im) {
    ncpx* z = malloc(sizeof *z);
    if (!z) { c->ok = false; return -1; }
    ncpx_init(z, c->prec_bits);
    if (!mgd_mpfr_from_expr(z->re, re) || !mgd_mpfr_from_expr(z->im, im)
        || !ctx_add_mgd_const(c, 2, z)) {
        ncpx_clear(z); free(z); return -1;
    }
    Slot k; memset(&k, 0, sizeof k); k.p = z;
    int dst = alloc_mgd(c, CT_BIGCOMPLEX);
    ins(c, OP_MG_CCONST, (uint32_t)dst, 0, 0, k);
    return dst;
}

/* Lift a managed real value to a managed complex one (zero imaginary part). */
static Val mgd_to_complex(Ctx* c, Val v) {
    if (v.type == CT_BIGCOMPLEX) return v;
    Slot k; memset(&k, 0, sizeof k); k.i = 0;
    int dst = alloc_mgd(c, CT_BIGCOMPLEX);
    ins(c, OP_MG_CFROM_R, (uint32_t)dst, (uint32_t)v.reg, 0, k);
    Val r = { dst, true, CT_BIGCOMPLEX, false };
    return r;
}

/* Real elementary-function head -> mpfr unary function pointer, or NULL. */
static const void* mgd_real_unary_fn(const char* h) {
    if (h == SYM_Sqrt) return (const void*)mpfr_sqrt;
    if (h == SYM_Exp) return (const void*)mpfr_exp;
    if (h == SYM_Log) return (const void*)mpfr_log;
    if (h == SYM_Sin) return (const void*)mpfr_sin;
    if (h == SYM_Cos) return (const void*)mpfr_cos;
    if (h == SYM_Tan) return (const void*)mpfr_tan;
    if (h == SYM_Csc) return (const void*)mpfr_csc;
    if (h == SYM_Sec) return (const void*)mpfr_sec;
    if (h == SYM_Cot) return (const void*)mpfr_cot;
    if (h == SYM_Sinh) return (const void*)mpfr_sinh;
    if (h == SYM_Cosh) return (const void*)mpfr_cosh;
    if (h == SYM_Tanh) return (const void*)mpfr_tanh;
    /* the reciprocal hyperbolic functions: a separated (Fourier) series such as
     * Sum (-1)^k/(2k+1) Sech[(2k+1) Pi/2] was interpreted per term without them */
    if (h == SYM_Sech) return (const void*)mpfr_sech;
    if (h == SYM_Csch) return (const void*)mpfr_csch;
    if (h == SYM_Coth) return (const void*)mpfr_coth;
    /* NaN outside the real domain (ArcCosh x < 1, ArcTanh |x| > 1): declined */
    if (h == SYM_ArcSinh) return (const void*)mpfr_asinh;
    if (h == SYM_ArcCosh) return (const void*)mpfr_acosh;
    if (h == SYM_ArcTanh) return (const void*)mpfr_atanh;
    if (h == SYM_ArcSin) return (const void*)mpfr_asin;
    if (h == SYM_ArcCos) return (const void*)mpfr_acos;
    if (h == SYM_ArcTan) return (const void*)mpfr_atan;
    if (h == SYM_Abs) return (const void*)mpfr_abs;
    /* NaN outside its real domain: the VM declines and the point is interpreted */
    if (h == SYM_ProductLog) return (const void*)productlog_mpfr_real;
    return NULL;
}
/* Complex elementary-function head -> MG_CUN selector, or -1.  The result of Abs
 * is real, reported via `*to_real`. */
static int mgd_complex_unary_sel(const char* h, bool* to_real) {
    *to_real = false;
    if (h == SYM_Exp) return MGC_EXP;
    if (h == SYM_Log) return MGC_LOG;
    if (h == SYM_Sin) return MGC_SIN;
    if (h == SYM_Cos) return MGC_COS;
    if (h == SYM_Sqrt) return MGC_SQRT;
    /* the reciprocal trig functions: 1/Sin[x]^201 evaluates to Csc[x]^201, so
     * without them a residue contour on |x| = 1 was interpreted per sample */
    if (h == SYM_Tan) return MGC_TAN;
    if (h == SYM_Csc) return MGC_CSC;
    if (h == SYM_Sec) return MGC_SEC;
    if (h == SYM_Cot) return MGC_COT;
    if (h == SYM_Abs) { *to_real = true; return MGC_ABS; }
    if (h == SYM_Arg) { *to_real = true; return MGC_ARG; }
    return -1;
}

static bool emit_mgd_real(Ctx* c, const Expr* e, Val* out);

/* Combine two managed-real/complex operands with a real mpfr binary fn or its
 * complex selector, promoting a real operand to complex when the other is. */
static Val mgd_bin(Ctx* c, Val a, Val b, const void* real_fn, int cplx_sel) {
    if (a.type == CT_BIGCOMPLEX || b.type == CT_BIGCOMPLEX) {
        a = mgd_to_complex(c, a); b = mgd_to_complex(c, b);
        Slot k; memset(&k, 0, sizeof k); k.i = cplx_sel;
        int dst = alloc_mgd(c, CT_BIGCOMPLEX);
        ins(c, OP_MG_CBIN, (uint32_t)dst, (uint32_t)a.reg, (uint32_t)b.reg, k);
        Val r = { dst, true, CT_BIGCOMPLEX, false };
        return r;
    }
    Slot k; memset(&k, 0, sizeof k); k.p = real_fn;
    int dst = alloc_mgd(c, CT_BIGREAL);
    ins(c, OP_MG_RBIN, (uint32_t)dst, (uint32_t)a.reg, (uint32_t)b.reg, k);
    Val r = { dst, true, CT_BIGREAL, false };
    return r;
}

/* r = x^(1/q) for x >= 0 (q, a small positive integer, carried as an MPFR
 * constant so the RBIN signature serves); NaN for x < 0. */
static int mgd_real_rootn(mpfr_ptr r, mpfr_srcptr x, mpfr_srcptr q, mpfr_rnd_t rnd) {
    if (mpfr_sgn(x) < 0) { mpfr_set_nan(r); return 0; }
    return mpfr_rootn_ui(r, x, mpfr_get_ui(q, MPFR_RNDN), rnd);
}

/* base^n for a machine-int exponent n. */
static Val mgd_powi(Ctx* c, Val base, long n) {
    Slot k; memset(&k, 0, sizeof k); k.i = (long long)n;
    if (base.type == CT_BIGCOMPLEX) {
        int dst = alloc_mgd(c, CT_BIGCOMPLEX);
        ins(c, OP_MG_CPOWI, (uint32_t)dst, (uint32_t)base.reg, 0, k);
        Val r = { dst, true, CT_BIGCOMPLEX, false };
        return r;
    }
    int dst = alloc_mgd(c, CT_BIGREAL);
    ins(c, OP_MG_RPOWI, (uint32_t)dst, (uint32_t)base.reg, 0, k);
    Val r = { dst, true, CT_BIGREAL, false };
    return r;
}

static bool emit_mgd_real_node(Ctx* c, const Expr* e, Val* out);

/* x = k * unit * rest with a numeric leading factor k*unit, k an integer,
 * 2 <= |k| <= 64 or k = -1, unit = 1 or a primitive Gaussian integer: returns
 * Power[Exp[unit rest], k] (caller frees), else NULL. */
static Expr* mgd_exp_harmonic(const Expr* x) {
    if (x->type != EXPR_FUNCTION || x->data.function.head->type != EXPR_SYMBOL
        || x->data.function.head->data.symbol.name != SYM_Times
        || x->data.function.arg_count < 2) return NULL;
    const Expr* c0 = x->data.function.args[0];
    int64_t re = 0, im = 0;
    Expr *cr, *ci;
    if (c0->type == EXPR_INTEGER) re = c0->data.integer;
    else if (is_complex((Expr*)c0, &cr, &ci) && cr->type == EXPR_INTEGER
             && ci->type == EXPR_INTEGER) { re = cr->data.integer; im = ci->data.integer; }
    else return NULL;
    int64_t g = re < 0 ? -re : re, b = im < 0 ? -im : im;
    while (b) { int64_t r = g % b; g = b; b = r; }      /* gcd(|re|, |im|) */
    if (g == 0) return NULL;
    /* the sign goes into k so that unit's first nonzero part is positive */
    int64_t k = ((re != 0 ? re : im) < 0) ? -g : g;
    if (!(k == -1 || (k >= 2 && k <= 64) || (k <= -2 && k >= -64))) return NULL;
    int64_t ur = re / k, ui = im / k;
    size_t n = x->data.function.arg_count;
    Expr** args = malloc(n * sizeof *args);
    if (!args) return NULL;
    size_t m = 0;
    if (!(ur == 1 && ui == 0))
        args[m++] = (ui == 0) ? expr_new_integer(ur)
                              : make_complex(expr_new_integer(ur), expr_new_integer(ui));
    for (size_t i = 1; i < n; i++) args[m++] = expr_copy(x->data.function.args[i]);
    Expr* u = (m == 1) ? args[0] : expr_new_function(expr_new_symbol(SYM_Times), args, m);
    free(args);
    Expr* ev[1] = { u };
    Expr* ex = expr_new_function(expr_new_symbol(SYM_Exp), ev, 1);
    Expr* pv[2] = { ex, expr_new_integer(k) };
    return expr_new_function(expr_new_symbol(SYM_Power), pv, 2);
}

/* The MPFR real/complex domain, with value numbering: a compound subtree that
 * was already lowered is not lowered again.  The managed programs skip the
 * optimiser's CSE (it never sees a managed register), so before this a contour
 * body (1 + E^(I t))^2000 (1 + E^(I t) + E^(2 I t))^1000 ... evaluated E^(I t)
 * -- an MPFR sin_cos at the working precision -- once per occurrence. */
static bool emit_mgd_real(Ctx* c, const Expr* e, Val* out) {
    if (!c->ok) return false;
    if (e->type != EXPR_FUNCTION) return emit_mgd_real_node(c, e, out);
    uint64_t h = (uint64_t)expr_hash(e);
    for (int i = 0; i < c->nmgd_memo; i++)
        if (c->mgd_memo_hash[i] == h && expr_eq(c->mgd_memo_key[i], e)) {
            *out = c->mgd_memo_val[i];
            return true;
        }
    if (!emit_mgd_real_node(c, e, out) || !c->ok) return false;
    if (c->nmgd_memo == c->mgd_memo_cap) {
        int nc = c->mgd_memo_cap ? 2 * c->mgd_memo_cap : 32;
        Expr** nk = realloc(c->mgd_memo_key, (size_t)nc * sizeof *nk);
        if (nk) c->mgd_memo_key = nk;
        uint64_t* nh = realloc(c->mgd_memo_hash, (size_t)nc * sizeof *nh);
        if (nh) c->mgd_memo_hash = nh;
        Val* nv = realloc(c->mgd_memo_val, (size_t)nc * sizeof *nv);
        if (nv) c->mgd_memo_val = nv;
        if (!nk || !nh || !nv) return true;          /* no memo: still correct */
        c->mgd_memo_cap = nc;
    }
    c->mgd_memo_key[c->nmgd_memo] = expr_copy((Expr*)e);
    c->mgd_memo_hash[c->nmgd_memo] = h;
    c->mgd_memo_val[c->nmgd_memo] = *out;
    c->nmgd_memo++;
    return true;
}

/* Copy of e with each symbol names[i] replaced by vals[i] (copied).  With is
 * substitution, so this IS its semantics; an inner construct that rebinds one
 * of the names gets a non-symbol in its binding position, which no managed
 * lowering accepts, so it bails rather than being mis-scoped. */
static Expr* mgd_subst(const Expr* e, const char* const* names, Expr* const* vals,
                       size_t n) {
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < n; i++)
            if (e->data.symbol.name == names[i]) return expr_copy(vals[i]);
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t na = e->data.function.arg_count;
    Expr** av = malloc((na ? na : 1) * sizeof *av);
    if (!av) return NULL;
    Expr* hd = mgd_subst(e->data.function.head, names, vals, n);
    bool ok = hd != NULL;
    for (size_t i = 0; i < na; i++) {
        av[i] = ok ? mgd_subst(e->data.function.args[i], names, vals, n) : NULL;
        if (!av[i]) ok = false;
    }
    if (!ok) {
        expr_free(hd);
        for (size_t i = 0; i < na; i++) expr_free(av[i]);
        free(av);
        return NULL;
    }
    Expr* r = expr_new_function(hd, av, na);
    free(av);
    return r;
}

/* With[{v1 = e1, ...}, body]: lowered as the substituted body.  Value numbering
 * (emit_mgd_real's memo) then evaluates each e_i once however often v_i occurs,
 * which is exactly what a register-bound local would give, without a managed
 * scope table that the memo would have to be keyed on.  Returns -1 when e is
 * not a With, else the lowering result (0/1). */
static int emit_mgd_with(Ctx* c, const Expr* e, Val* out) {
    if (e->data.function.arg_count != 2) return -1;
    const Expr* L = e->data.function.args[0];
    if (L->type != EXPR_FUNCTION || L->data.function.head->type != EXPR_SYMBOL
        || L->data.function.head->data.symbol.name != SYM_List) return 0;
    size_t nl = L->data.function.arg_count;
    if (nl == 0 || nl > 32) return 0;
    const char* names[32];
    Expr* vals[32];
    for (size_t i = 0; i < nl; i++) {
        const Expr* sp = L->data.function.args[i];
        if (sp->type != EXPR_FUNCTION || sp->data.function.head->type != EXPR_SYMBOL
            || sp->data.function.arg_count != 2
            || (sp->data.function.head->data.symbol.name != SYM_Set
                && sp->data.function.head->data.symbol.name != SYM_SetDelayed)
            || sp->data.function.args[0]->type != EXPR_SYMBOL) return 0;
        names[i] = sp->data.function.args[0]->data.symbol.name;
        vals[i] = sp->data.function.args[1];
        /* a local shadowing a declared argument is still substitution: the
         * body's occurrences are the local, the argument is unreachable */
    }
    Expr* body = mgd_subst(e->data.function.args[1], names, vals, nl);
    if (!body) { c->ok = false; return 0; }
    bool ok = emit_mgd_real(c, body, out);
    expr_free(body);
    return (ok && c->ok) ? 1 : 0;
}

/* ---- special functions (MG_SPECIAL) and closed constants ---- */

/* A named constant numericalize knows (the literal-parameter test below). */
static bool mgd_named_constant(const char* nm) {
    double d;
    return nm == SYM_Pi || nm == SYM_E || nm == SYM_EulerGamma || nm == SYM_Catalan
        || nm == SYM_Degree || numeric_constant_machine_value(nm, &d);
}

/* Is `e` CLOSED — free of every declared argument and built only from numbers,
 * named constants, Lists and NumericFunction heads — so that its value is the
 * same at every call?  Such a subtree is a literal parameter of a special
 * function (the 2 of BesselJ[2, x], the {1/3} of HypergeometricPFQ[{1/3}, ...])
 * or, as a whole call, a constant.  The NumericFunction requirement is what
 * keeps RandomReal[] or a user function with side effects out. */
static bool mgd_closed(Ctx* c, const Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_REAL: case EXPR_BIGINT:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_SYMBOL:
            return arg_find(c, e->data.symbol.name) < 0
                && mgd_named_constant(e->data.symbol.name);
        case EXPR_FUNCTION: {
            const Expr* hd = e->data.function.head;
            if (hd->type != EXPR_SYMBOL) return false;
            const char* nm = hd->data.symbol.name;
            if (nm != SYM_List && nm != SYM_Rational && nm != SYM_Complex) {
                SymbolDef* d = symtab_lookup(nm);
                if (!d || !(d->attributes & ATTR_NUMERICFUNCTION)) return false;
            }
            for (size_t i = 0; i < e->data.function.arg_count; i++)
                if (!mgd_closed(c, e->data.function.args[i])) return false;
            return true;
        }
        default:
            return false;
    }
}

/* A closed call as a managed constant: N[e] at the working precision, once, at
 * compile time — what the interpreter computes at every sample. */
static bool mgd_fold_const(Ctx* c, const Expr* e, Val* out) {
    NumericSpec sp; sp.mode = NUMERIC_MODE_MPFR; sp.bits = c->prec_bits;
    Expr* v = numericalize(e, sp);
    if (!v) return false;
    Expr *re, *im;
    int cst = -1;
    CompileType t = CT_BIGREAL;
    if (v->type == EXPR_MPFR) cst = mgd_new_mpfr_const_expr(c, v);
    else if (is_complex(v, &re, &im) && (re->type == EXPR_MPFR || im->type == EXPR_MPFR)) {
        cst = mgd_new_ncpx_const(c, re, im);
        t = CT_BIGCOMPLEX;
    }
    expr_free(v);
    if (cst < 0) return false;
    Val r = { cst, true, t, false }; *out = r;
    return c->ok;
}

/* Lower f[args] for a head in the special-function registry: every closed
 * argument becomes a literal parameter (evaluated once, passed exactly as the
 * interpreter sees it — the integer order of BesselJ[2, x] stays an integer),
 * every other argument a managed register, moved into one consecutive operand
 * range.  The result is complex when any operand is; a real program meeting a
 * complex value at run time declines that point. */
static bool emit_mgd_special(Ctx* c, const Expr* e, const char* h,
                             Expr* const* a, size_t n, Val* out) {
    (void)e;
    const MgdSpecialEntry* ent = mgd_special_lookup(h, n);
    if (!ent || n > 16) return false;
    Val ops[16];
    size_t nreg = 0;
    bool any_cplx = false;
    MgdSpecialCall* sc = calloc(1, sizeof *sc);
    if (!sc) { c->ok = false; return false; }
    sc->head = h; sc->nargs = n;
    sc->lit = calloc(n ? n : 1, sizeof *sc->lit);
    sc->cplx = calloc(n ? n : 1, 1);
    if (!sc->lit || !sc->cplx) { mgd_special_call_free(sc); c->ok = false; return false; }
    for (size_t i = 0; i < n; i++) {
        if (mgd_closed(c, a[i])) {
            Expr* cp = expr_copy(a[i]);
            sc->lit[i] = evaluate(cp);
            expr_free(cp);
            if (!sc->lit[i]) { mgd_special_call_free(sc); return false; }
            continue;
        }
        if (!emit_mgd_real(c, a[i], &ops[nreg])) { mgd_special_call_free(sc); return false; }
        sc->cplx[i] = ops[nreg].type == CT_BIGCOMPLEX;
        any_cplx = any_cplx || sc->cplx[i];
        nreg++;
    }
    if (nreg == 0 || (any_cplx && !(ent->flags & MGS_CPLX))
        || !ctx_add_mgd_const(c, MGD_CONST_SPECIAL, sc)) {
        mgd_special_call_free(sc);
        return false;
    }
    /* ownership of sc is the program's from here (freed with the consts) */
    int first = -1;
    int regs[16];
    for (size_t j = 0; j < nreg; j++) {
        regs[j] = alloc_mgd(c, ops[j].type);
        if (j == 0) first = regs[j];
    }
    Slot z; memset(&z, 0, sizeof z);
    for (size_t j = 0; j < nreg; j++)
        ins(c, ops[j].type == CT_BIGCOMPLEX ? OP_MG_CMOV : OP_MG_RMOV,
            (uint32_t)regs[j], (uint32_t)ops[j].reg, 0, z);
    /* complex when an operand is, or when the head is complex-valued on the
     * real line (HankelH1 = J + I Y) */
    bool cout = any_cplx || (ent->flags & MGS_COUT);
    CompileType rt = cout ? CT_BIGCOMPLEX : CT_BIGREAL;
    sc->out_cplx = cout;
    int dst = alloc_mgd(c, rt);
    Slot k; memset(&k, 0, sizeof k); k.p = sc;
    ins_f(c, OP_MG_SPECIAL, (uint16_t)nreg, (uint32_t)dst, (uint32_t)first, 0, k);
    Val r = { dst, true, rt, false }; *out = r;
    return c->ok;
}

static bool emit_mgd_real_node(Ctx* c, const Expr* e, Val* out) {
    if (!c->ok) return false;
    /* leaf: declared argument */
    if (e->type == EXPR_SYMBOL) {
        int reg = arg_find(c, e->data.symbol.name);
        if (reg >= 0) {
            CompileType t = c->arg_types[reg];
            if (t == CT_BIGREAL || t == CT_BIGCOMPLEX) {
                Val v = { reg, false, t, false }; *out = v; return true;
            }
            return false;
        }
        int cst = mgd_new_mpfr_const_named(c, e->data.symbol.name);
        if (cst < 0) return false;
        Val v = { cst, true, CT_BIGREAL, false }; *out = v; return true;
    }
    /* leaf: numeric literal */
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT
#ifdef USE_MPFR
        || e->type == EXPR_MPFR
#endif
       ) {
        int cst = mgd_new_mpfr_const_expr(c, e);
        if (cst < 0) return false;
        Val v = { cst, true, CT_BIGREAL, false }; *out = v; return true;
    }
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    Expr* const* a = e->data.function.args;
    size_t n = e->data.function.arg_count;

    /* Complex[re, im] literal (both parts numeric) */
    if (h == SYM_Complex && n == 2) {
        int cst = mgd_new_ncpx_const(c, a[0], a[1]);
        if (cst < 0) return false;
        Val v = { cst, true, CT_BIGCOMPLEX, false }; *out = v; return true;
    }
    /* rational literal Rational[p, q] */
    if (mgd_is_exact_rational(e)) {
        int cst = mgd_new_mpfr_const_expr(c, e);
        if (cst < 0) return false;
        Val v = { cst, true, CT_BIGREAL, false }; *out = v; return true;
    }

    if (h == SYM_With) {
        int w = emit_mgd_with(c, e, out);
        if (w >= 0) return w == 1;
    }
    if (h == SYM_Plus && n >= 1) {
        Val acc; if (!emit_mgd_real(c, a[0], &acc)) return false;
        for (size_t i = 1; i < n; i++) {
            Val v; if (!emit_mgd_real(c, a[i], &v)) return false;
            acc = mgd_bin(c, acc, v, (const void*)mpfr_add, MGC_ADD);
        }
        *out = acc; return c->ok;
    }
    if (h == SYM_Times && n >= 1) {
        Val acc; if (!emit_mgd_real(c, a[0], &acc)) return false;
        for (size_t i = 1; i < n; i++) {
            Val v; if (!emit_mgd_real(c, a[i], &v)) return false;
            acc = mgd_bin(c, acc, v, (const void*)mpfr_mul, MGC_MUL);
        }
        *out = acc; return c->ok;
    }
    if (h == SYM_Subtract && n == 2) {
        Val x, y; if (!emit_mgd_real(c, a[0], &x) || !emit_mgd_real(c, a[1], &y)) return false;
        *out = mgd_bin(c, x, y, (const void*)mpfr_sub, MGC_SUB); return c->ok;
    }
    if (h == SYM_Divide && n == 2) {
        Val x, y; if (!emit_mgd_real(c, a[0], &x) || !emit_mgd_real(c, a[1], &y)) return false;
        *out = mgd_bin(c, x, y, (const void*)mpfr_div, MGC_DIV); return c->ok;
    }
    if (h == SYM_Power && n == 2) {
        /* E^x is exp(x): lowering it as a general power would take log E (an
         * AGM at the working precision) and, for a complex x, a complex log and
         * exp, on every evaluation. */
        if (a[0]->type == EXPR_SYMBOL && a[0]->data.symbol.name == SYM_E) {
            /* E^(k u) for an integer multiple k of a unit coefficient (1 or a
             * Gaussian-integer direction such as I) is lowered as (E^u)^k: with
             * value numbering, the harmonics E^(2 I t), E^(3 I t), ... of a
             * contour body share one exponential, each costing a few products
             * instead of an MPFR sin_cos at the working precision. */
            Expr* t = mgd_exp_harmonic(a[1]);
            if (!t) {
                Expr* av[1] = { expr_copy(a[1]) };
                t = expr_new_function(expr_new_symbol(SYM_Exp), av, 1);
            }
            bool ok = emit_mgd_real(c, t, out);
            expr_free(t);
            return ok && c->ok;
        }
        /* b^(p/2) = Sqrt[b]^p exactly on the principal branch (both are
         * exp((p/2) Log b)): a square root and a few multiplications instead
         * of a log and an exp. */
        {
            int64_t pn, pd;
            if (is_rational(a[1], &pn, &pd) && pd == 2 && pn >= -64 && pn <= 64) {
                Expr* av[1] = { expr_copy(a[0]) };
                Expr* t = expr_new_function(expr_new_symbol(SYM_Sqrt), av, 1);
                Val sq;
                bool ok = emit_mgd_real(c, t, &sq);
                expr_free(t);
                if (!ok || !c->ok) return false;
                *out = (pn == 1) ? sq : mgd_powi(c, sq, (long)pn);
                return c->ok;
            }
        }
        Val base; if (!emit_mgd_real(c, a[0], &base)) return false;
        long ni;
        if (mgd_small_int(a[1], &ni)) { *out = mgd_powi(c, base, ni); return c->ok; }
        /* b^(p/q), real b, 3 <= q <= 64, |p| <= 64: the q-th root (a Newton
         * iteration, correctly rounded) and p products instead of exp(p/q log b)
         * -- an AGM log and an exp at the working precision.  A negative base
         * yields NaN, so the point is interpreted (its principal power is
         * complex), exactly as the log/exp route declines it. */
        {
            int64_t pn, pd;
            if (base.type == CT_BIGREAL && is_rational(a[1], &pn, &pd)
                && pd >= 3 && pd <= 64 && pn >= -64 && pn <= 64) {
                Expr* qe = expr_new_integer(pd);
                int qc = mgd_new_mpfr_const_expr(c, qe);
                expr_free(qe);
                if (qc < 0) return false;
                Val q = { qc, true, CT_BIGREAL, false };
                Val r = mgd_bin(c, base, q, (const void*)mgd_real_rootn, MGC_POW);
                *out = (pn == 1) ? r : mgd_powi(c, r, (long)pn);
                return c->ok;
            }
        }
        Val ex; if (!emit_mgd_real(c, a[1], &ex)) return false;
        *out = mgd_bin(c, base, ex, (const void*)numeric_mpfr_pow_fast, MGC_POW); return c->ok;
    }
    if (h == SYM_Minus && n == 1) {
        Val x; if (!emit_mgd_real(c, a[0], &x)) return false;
        if (x.type == CT_BIGCOMPLEX) {
            Slot k; memset(&k, 0, sizeof k); k.i = MGC_NEG;
            int dst = alloc_mgd(c, CT_BIGCOMPLEX);
            ins(c, OP_MG_CUN, (uint32_t)dst, (uint32_t)x.reg, 0, k);
            Val r = { dst, true, CT_BIGCOMPLEX, false }; *out = r; return c->ok;
        }
        Slot k; memset(&k, 0, sizeof k); k.p = (const void*)mpfr_neg;
        int dst = alloc_mgd(c, CT_BIGREAL);
        ins(c, OP_MG_RUN, (uint32_t)dst, (uint32_t)x.reg, 0, k);
        Val r = { dst, true, CT_BIGREAL, false }; *out = r; return c->ok;
    }
    /* unary elementary functions: a direct mpfr_* / ncpx_* opcode where one
     * exists for the operand's domain; otherwise (Sinh of a complex, Erf, ...)
     * fall through to the special-function registry below.  The operand is
     * value-numbered, so lowering it here and again there emits it once. */
    bool tr_;
    if (n == 1 && (mgd_real_unary_fn(h) || mgd_complex_unary_sel(h, &tr_) >= 0)) {
        Val x; if (!emit_mgd_real(c, a[0], &x)) return false;
        if (x.type == CT_BIGCOMPLEX) {
            bool to_real; int sel = mgd_complex_unary_sel(h, &to_real);
            if (sel >= 0) {
                Slot k; memset(&k, 0, sizeof k); k.i = sel;
                CompileType rt = to_real ? CT_BIGREAL : CT_BIGCOMPLEX;
                int dst = alloc_mgd(c, rt);
                ins(c, OP_MG_CUN, (uint32_t)dst, (uint32_t)x.reg, 0, k);
                Val r = { dst, true, rt, false }; *out = r; return c->ok;
            }
        } else {
            const void* fn = mgd_real_unary_fn(h);
            if (fn) {
                Slot k; memset(&k, 0, sizeof k); k.p = fn;
                int dst = alloc_mgd(c, CT_BIGREAL);
                ins(c, OP_MG_RUN, (uint32_t)dst, (uint32_t)x.reg, 0, k);
                Val r = { dst, true, CT_BIGREAL, false }; *out = r; return c->ok;
            }
        }
    }
    /* a closed numeric call (BesselJZero[0, 3], Gamma[1/3]) is a constant at
     * the working precision; anything else of a registered head is lowered to
     * its own arbitrary-precision evaluator */
    if (mgd_closed(c, e)) return mgd_fold_const(c, e, out);
    return emit_mgd_special(c, e, h, a, n, out);
}
#else
static bool emit_mgd_real(Ctx* c, const Expr* e, Val* out) {
    (void)c; (void)e; (void)out; return false;
}
#endif /* USE_MPFR */

/* The GMP (bignum) integer domain: exact Plus/Times/Subtract/Minus and
 * non-negative integer Power.  Division is a Rational, which no managed-int
 * container holds, so it bails. */
static bool emit_mgd_int(Ctx* c, const Expr* e, Val* out) {
    if (!c->ok) return false;
    if (e->type == EXPR_SYMBOL) {
        int reg = arg_find(c, e->data.symbol.name);
        if (reg >= 0 && c->arg_types[reg] == CT_BIGINT) {
            Val v = { reg, false, CT_BIGINT, false }; *out = v; return true;
        }
        return false;
    }
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) {
        int cst = mgd_new_mpz_const(c, e);
        if (cst < 0) return false;
        Val v = { cst, true, CT_BIGINT, false }; *out = v; return true;
    }
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    Expr* const* a = e->data.function.args;
    size_t n = e->data.function.arg_count;

    if ((h == SYM_Plus || h == SYM_Times) && n >= 1) {
        uint16_t op = (h[0] == 'P') ? OP_MG_ZADD : OP_MG_ZMUL;
        Val acc; if (!emit_mgd_int(c, a[0], &acc)) return false;
        for (size_t i = 1; i < n; i++) {
            Val v; if (!emit_mgd_int(c, a[i], &v)) return false;
            Slot z; memset(&z, 0, sizeof z);
            int dst = alloc_mgd(c, CT_BIGINT);
            ins(c, op, (uint32_t)dst, (uint32_t)acc.reg, (uint32_t)v.reg, z);
            Val r = { dst, true, CT_BIGINT, false }; acc = r;
        }
        *out = acc; return c->ok;
    }
    if (h == SYM_Subtract && n == 2) {
        Val x, y; if (!emit_mgd_int(c, a[0], &x) || !emit_mgd_int(c, a[1], &y)) return false;
        Slot z; memset(&z, 0, sizeof z);
        int dst = alloc_mgd(c, CT_BIGINT);
        ins(c, OP_MG_ZSUB, (uint32_t)dst, (uint32_t)x.reg, (uint32_t)y.reg, z);
        Val r = { dst, true, CT_BIGINT, false }; *out = r; return c->ok;
    }
    if (h == SYM_Minus && n == 1) {
        Val x; if (!emit_mgd_int(c, a[0], &x)) return false;
        Slot z; memset(&z, 0, sizeof z);
        int dst = alloc_mgd(c, CT_BIGINT);
        ins(c, OP_MG_ZNEG, (uint32_t)dst, (uint32_t)x.reg, 0, z);
        Val r = { dst, true, CT_BIGINT, false }; *out = r; return c->ok;
    }
    if (h == SYM_Power && n == 2) {
        long ni;
        if (!mgd_small_int(a[1], &ni) || ni < 0) return false;   /* negative -> Rational */
        Val base; if (!emit_mgd_int(c, a[0], &base)) return false;
        Slot k; memset(&k, 0, sizeof k); k.i = (long long)ni;
        int dst = alloc_mgd(c, CT_BIGINT);
        ins(c, OP_MG_ZPOWI, (uint32_t)dst, (uint32_t)base.reg, 0, k);
        Val r = { dst, true, CT_BIGINT, false }; *out = r; return c->ok;
    }
    return false;
}

/* Managed emit entry: route to the domain the program was compiled for, and
 * record the innermost failing node for the bail diagnostic (first writer wins,
 * mirroring the machine emit wrapper). */
bool emit_mgd(Ctx* c, const Expr* e, Val* out) {
    bool ok = (c->prec_bits > 0) ? emit_mgd_real(c, e, out)
            : c->allow_bigint    ? emit_mgd_int(c, e, out)
            : false;
    for (int i = 0; i < c->nmgd_memo; i++) expr_free(c->mgd_memo_key[i]);
    free(c->mgd_memo_key); free(c->mgd_memo_hash); free(c->mgd_memo_val);
    c->mgd_memo_key = NULL; c->mgd_memo_hash = NULL; c->mgd_memo_val = NULL;
    c->nmgd_memo = c->mgd_memo_cap = 0;
    if ((!ok || !c->ok) && !c->bail_node) c->bail_node = e;
    return ok && c->ok;
}
