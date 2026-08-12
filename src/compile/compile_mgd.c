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
#ifdef USE_MPFR
#include "../numeric_complex.h"  /* ncpx, numeric_mpfr_* — the managed containers */
#endif
#include <math.h>
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
bool mgd_mpz_from_expr(mpz_ptr z, const Expr* e) {
    if (e->type == EXPR_INTEGER) { mpz_set_si(z, (long)e->data.integer); return true; }
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
/* ---- MPFR (real) constant materialisation ---- */
bool mgd_mpfr_from_expr(mpfr_ptr m, const Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(m, (long)e->data.integer, MPFR_RNDN); return true;
        case EXPR_REAL:    mpfr_set_d (m, e->data.real,          MPFR_RNDN); return true;
        case EXPR_BIGINT:  mpfr_set_z (m, e->data.bigint,        MPFR_RNDN); return true;
        case EXPR_MPFR:    mpfr_set   (m, e->data.mpfr,          MPFR_RNDN); return true;
        default: break;
    }
    int64_t nn, dd;
    if (is_rational(e, &nn, &dd) && dd != 0) {
        mpfr_t den; mpfr_init2(den, mpfr_get_prec(m));
        mpfr_set_si(m, (long)nn, MPFR_RNDN);
        mpfr_set_si(den, (long)dd, MPFR_RNDN);
        mpfr_div(m, m, den, MPFR_RNDN);
        mpfr_clear(den);
        return true;
    }
    return false;
}
/* A named real constant (Pi, E, ...) at working precision. */
static bool mgd_mpfr_from_const(mpfr_ptr m, const char* nm) {
    if (strcmp(nm, "Pi") == 0)         { mpfr_const_pi(m, MPFR_RNDN); return true; }
    if (strcmp(nm, "E") == 0)          { mpfr_set_ui(m, 1, MPFR_RNDN); mpfr_exp(m, m, MPFR_RNDN); return true; }
    if (strcmp(nm, "EulerGamma") == 0) { mpfr_const_euler(m, MPFR_RNDN); return true; }
    if (strcmp(nm, "Catalan") == 0)    { mpfr_const_catalan(m, MPFR_RNDN); return true; }
    if (strcmp(nm, "Degree") == 0)     { mpfr_const_pi(m, MPFR_RNDN); mpfr_div_ui(m, m, 180, MPFR_RNDN); return true; }
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
    if (strcmp(h, "Sqrt") == 0) return (const void*)mpfr_sqrt;
    if (strcmp(h, "Exp")  == 0) return (const void*)mpfr_exp;
    if (strcmp(h, "Log")  == 0) return (const void*)mpfr_log;
    if (strcmp(h, "Sin")  == 0) return (const void*)mpfr_sin;
    if (strcmp(h, "Cos")  == 0) return (const void*)mpfr_cos;
    if (strcmp(h, "Tan")  == 0) return (const void*)mpfr_tan;
    if (strcmp(h, "Sinh") == 0) return (const void*)mpfr_sinh;
    if (strcmp(h, "Cosh") == 0) return (const void*)mpfr_cosh;
    if (strcmp(h, "Tanh") == 0) return (const void*)mpfr_tanh;
    if (strcmp(h, "ArcSin") == 0) return (const void*)mpfr_asin;
    if (strcmp(h, "ArcCos") == 0) return (const void*)mpfr_acos;
    if (strcmp(h, "ArcTan") == 0) return (const void*)mpfr_atan;
    if (strcmp(h, "Abs")  == 0) return (const void*)mpfr_abs;
    return NULL;
}
/* Complex elementary-function head -> MG_CUN selector, or -1.  The result of Abs
 * is real, reported via `*to_real`. */
static int mgd_complex_unary_sel(const char* h, bool* to_real) {
    *to_real = false;
    if (strcmp(h, "Exp")  == 0) return MGC_EXP;
    if (strcmp(h, "Log")  == 0) return MGC_LOG;
    if (strcmp(h, "Sin")  == 0) return MGC_SIN;
    if (strcmp(h, "Cos")  == 0) return MGC_COS;
    if (strcmp(h, "Sqrt") == 0) return MGC_SQRT;
    if (strcmp(h, "Abs")  == 0) { *to_real = true; return MGC_ABS; }
    if (strcmp(h, "Arg")  == 0) { *to_real = true; return MGC_ARG; }
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

/* The MPFR real/complex domain. */
static bool emit_mgd_real(Ctx* c, const Expr* e, Val* out) {
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
    if (strcmp(h, "Complex") == 0 && n == 2) {
        int cst = mgd_new_ncpx_const(c, a[0], a[1]);
        if (cst < 0) return false;
        Val v = { cst, true, CT_BIGCOMPLEX, false }; *out = v; return true;
    }
    /* rational literal Rational[p, q] */
    { int64_t nn, dd; if (is_rational(e, &nn, &dd)) {
        int cst = mgd_new_mpfr_const_expr(c, e);
        if (cst < 0) return false;
        Val v = { cst, true, CT_BIGREAL, false }; *out = v; return true;
    } }

    if (strcmp(h, "Plus") == 0 && n >= 1) {
        Val acc; if (!emit_mgd_real(c, a[0], &acc)) return false;
        for (size_t i = 1; i < n; i++) {
            Val v; if (!emit_mgd_real(c, a[i], &v)) return false;
            acc = mgd_bin(c, acc, v, (const void*)mpfr_add, MGC_ADD);
        }
        *out = acc; return c->ok;
    }
    if (strcmp(h, "Times") == 0 && n >= 1) {
        Val acc; if (!emit_mgd_real(c, a[0], &acc)) return false;
        for (size_t i = 1; i < n; i++) {
            Val v; if (!emit_mgd_real(c, a[i], &v)) return false;
            acc = mgd_bin(c, acc, v, (const void*)mpfr_mul, MGC_MUL);
        }
        *out = acc; return c->ok;
    }
    if (strcmp(h, "Subtract") == 0 && n == 2) {
        Val x, y; if (!emit_mgd_real(c, a[0], &x) || !emit_mgd_real(c, a[1], &y)) return false;
        *out = mgd_bin(c, x, y, (const void*)mpfr_sub, MGC_SUB); return c->ok;
    }
    if (strcmp(h, "Divide") == 0 && n == 2) {
        Val x, y; if (!emit_mgd_real(c, a[0], &x) || !emit_mgd_real(c, a[1], &y)) return false;
        *out = mgd_bin(c, x, y, (const void*)mpfr_div, MGC_DIV); return c->ok;
    }
    if (strcmp(h, "Power") == 0 && n == 2) {
        Val base; if (!emit_mgd_real(c, a[0], &base)) return false;
        long ni;
        if (mgd_small_int(a[1], &ni)) { *out = mgd_powi(c, base, ni); return c->ok; }
        Val ex; if (!emit_mgd_real(c, a[1], &ex)) return false;
        *out = mgd_bin(c, base, ex, (const void*)mpfr_pow, MGC_POW); return c->ok;
    }
    if (strcmp(h, "Minus") == 0 && n == 1) {
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
    /* unary elementary functions */
    if (n == 1) {
        Val x; if (!emit_mgd_real(c, a[0], &x)) return false;
        if (x.type == CT_BIGCOMPLEX) {
            bool to_real; int sel = mgd_complex_unary_sel(h, &to_real);
            if (sel < 0) return false;
            Slot k; memset(&k, 0, sizeof k); k.i = sel;
            CompileType rt = to_real ? CT_BIGREAL : CT_BIGCOMPLEX;
            int dst = alloc_mgd(c, rt);
            ins(c, OP_MG_CUN, (uint32_t)dst, (uint32_t)x.reg, 0, k);
            Val r = { dst, true, rt, false }; *out = r; return c->ok;
        }
        const void* fn = mgd_real_unary_fn(h);
        if (!fn) return false;
        Slot k; memset(&k, 0, sizeof k); k.p = fn;
        int dst = alloc_mgd(c, CT_BIGREAL);
        ins(c, OP_MG_RUN, (uint32_t)dst, (uint32_t)x.reg, 0, k);
        Val r = { dst, true, CT_BIGREAL, false }; *out = r; return c->ok;
    }
    return false;
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

    if ((strcmp(h, "Plus") == 0 || strcmp(h, "Times") == 0) && n >= 1) {
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
    if (strcmp(h, "Subtract") == 0 && n == 2) {
        Val x, y; if (!emit_mgd_int(c, a[0], &x) || !emit_mgd_int(c, a[1], &y)) return false;
        Slot z; memset(&z, 0, sizeof z);
        int dst = alloc_mgd(c, CT_BIGINT);
        ins(c, OP_MG_ZSUB, (uint32_t)dst, (uint32_t)x.reg, (uint32_t)y.reg, z);
        Val r = { dst, true, CT_BIGINT, false }; *out = r; return c->ok;
    }
    if (strcmp(h, "Minus") == 0 && n == 1) {
        Val x; if (!emit_mgd_int(c, a[0], &x)) return false;
        Slot z; memset(&z, 0, sizeof z);
        int dst = alloc_mgd(c, CT_BIGINT);
        ins(c, OP_MG_ZNEG, (uint32_t)dst, (uint32_t)x.reg, 0, z);
        Val r = { dst, true, CT_BIGINT, false }; *out = r; return c->ok;
    }
    if (strcmp(h, "Power") == 0 && n == 2) {
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
    if ((!ok || !c->ok) && !c->bail_node) c->bail_node = e;
    return ok && c->ok;
}
