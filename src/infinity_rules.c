/* infinity_rules.c -- numeric special values for NumericFunction heads.
 *
 * One mechanism instead of a hand patch per head (audit Wave 3, shared fix
 * S11). The evaluator classifies the arguments of every Listable call anyway
 * (eval.c, listable_arg_kinds); that pass also notes an Indeterminate argument
 * and an infinite one (Infinity, ComplexInfinity, c Infinity), and only then
 * calls numeric_special_value_rule. So Sin[1.2] pays nothing new, and the rule
 * runs before UpValues, DownValues and the builtin, for C and .m heads alike.
 *
 * 1. Indeterminate. A NumericFunction is a function of the numeric value of
 *    its arguments; an argument whose value is undetermined leaves the value
 *    undetermined, so f[..., Indeterminate, ...] is Indeterminate. That holds
 *    for every non-Flat NumericFunction head (the Flat ones -- Plus, Times,
 *    Max, Min -- carry their own rules and are not routed here). The two
 *    list-valued heads, ReIm and QuotientRemainder, give a pair of them.
 *
 * 2. Infinities. Where the limit along a direction is systematic, a row of
 *    INF_ROWS gives it for each direction class (Infinity, -Infinity,
 *    I Infinity, -I Infinity, ComplexInfinity, any other direction). Every
 *    entry is the limit of the function along the ray (an essential
 *    singularity -- ComplexInfinity into Sin, Exp, Erf -- or an oscillation
 *    without a limit is Indeterminate; a bounded oscillation of a real
 *    function along the real axis is Interval[{-1, 1}]). The branch of each
 *    inverse function follows Mathilda's own numerics (ArcSin[10.^8] =
 *    1.5708 - 19.11 I, so ArcSin[Infinity] = -I Infinity); the table was
 *    checked against them at |z| = 10^8.
 *    Heads that READ the direction (Abs, Re, Im, ReIm, Sign, Arg, Conjugate),
 *    Power with an infinite base or exponent, and the classical orthogonal /
 *    Bernoulli polynomials of explicit degree (whose value at infinity is the
 *    leading term) are rules rather than rows. */

#include "infinity_rules.h"
#include "eval.h"
#include "arithmetic.h"
#include "numeric.h"
#include "sym_names.h"
#include "sym_intern.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ---- small predicates and constructors --------------------------------- */

static bool sym_is(const Expr* e, const char* s) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == s;
}

static bool head_sym(const Expr* e, const char* s) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == s;
}

/* Sign (-1, 0, 1) of a real number leaf (Integer, BigInt, Rational, Real,
 * MPFR); 2 for anything else. */
static int real_number_sign(const Expr* e) {
    if (!e) return 2;
    switch (e->type) {
        case EXPR_INTEGER: return e->data.integer > 0 ? 1 : e->data.integer < 0 ? -1 : 0;
        case EXPR_BIGINT:  return mpz_sgn(e->data.bigint);
        case EXPR_REAL:    return e->data.real > 0.0 ? 1 : e->data.real < 0.0 ? -1 :
                                  (e->data.real == 0.0 ? 0 : 2);   /* NaN: 2 */
#ifdef USE_MPFR
        case EXPR_MPFR:    return mpfr_nan_p(e->data.mpfr) ? 2 : mpfr_sgn(e->data.mpfr);
#endif
        default: break;
    }
    if (head_sym(e, SYM_Rational) && e->data.function.arg_count == 2) {
        int a = real_number_sign(e->data.function.args[0]);
        int b = real_number_sign(e->data.function.args[1]);
        if (a == 2 || b == 2) return 2;
        return a * b;
    }
    return 2;
}

static Expr* mk2(const char* h, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(h), args, 2);
}
static Expr* mk1(const char* h, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(h), args, 1);
}

static Expr* mk_inf(void)  { return expr_new_symbol(SYM_Infinity); }
static Expr* mk_ninf(void) { return mk2(SYM_Times, expr_new_integer(-1), mk_inf()); }
static Expr* mk_cinf(void) { return expr_new_symbol(SYM_ComplexInfinity); }
static Expr* mk_indet(void) { return expr_new_symbol(SYM_Indeterminate); }
static Expr* mk_unit_i(int s) { return make_complex(expr_new_integer(0), expr_new_integer(s)); }
/* c Infinity, evaluated (Times normalises the direction). */
static Expr* mk_dir_inf(Expr* c) { return eval_and_free(mk2(SYM_Times, c, mk_inf())); }

/* ---- direction classification ------------------------------------------ */

InfDir infinity_direction(const Expr* e, Expr** dir) {
    if (dir) *dir = NULL;
    if (!e) return INF_DIR_NONE;
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol.name == SYM_Infinity) return INF_DIR_POS;
        if (e->data.symbol.name == SYM_ComplexInfinity) return INF_DIR_CINF;
        return INF_DIR_NONE;
    }
    if (!head_sym(e, SYM_Times)) return INF_DIR_NONE;
    size_t ac = e->data.function.arg_count;
    if (ac < 2 || !sym_is(e->data.function.args[ac - 1], SYM_Infinity)) return INF_DIR_NONE;
    if (ac == 2) {
        Expr* c = e->data.function.args[0];
        int s = real_number_sign(c);
        if (s == 1) return INF_DIR_POS;
        if (s == -1) return INF_DIR_NEG;
        if (s == 0) return INF_DIR_NONE;
        Expr *re, *im;
        if (is_complex(c, &re, &im) && real_number_sign(re) == 0) {
            int si = real_number_sign(im);
            if (si == 1) return INF_DIR_I;
            if (si == -1) return INF_DIR_NEG_I;
        }
    }
    for (size_t i = 0; i + 1 < ac; i++) {
        const Expr* f = e->data.function.args[i];
        if (sym_is(f, SYM_Infinity) || sym_is(f, SYM_ComplexInfinity) ||
            sym_is(f, SYM_Indeterminate) || !expr_is_numeric_quantity(f))
            return INF_DIR_NONE;
    }
    if (dir) {
        if (ac == 2) *dir = expr_copy(e->data.function.args[0]);
        else {
            Expr** a = malloc(sizeof(Expr*) * (ac - 1));
            for (size_t i = 0; i + 1 < ac; i++) a[i] = expr_copy(e->data.function.args[i]);
            *dir = expr_new_function(expr_new_symbol(SYM_Times), a, ac - 1);
            free(a);
        }
    }
    return INF_DIR_OTHER;
}

/* The direction of a classified infinity as a (new) expression. */
static Expr* direction_expr(InfDir d, const Expr* e) {
    switch (d) {
        case INF_DIR_POS:   return expr_new_integer(1);
        case INF_DIR_NEG:   return expr_new_integer(-1);
        case INF_DIR_I:     return mk_unit_i(1);
        case INF_DIR_NEG_I: return mk_unit_i(-1);
        case INF_DIR_OTHER: { Expr* dir = NULL; infinity_direction(e, &dir); return dir; }
        default: return NULL;
    }
}

int numeric_special_arg_kind(const Expr* arg) {
    if (!arg) return 0;
    if (arg->type == EXPR_SYMBOL) {
        const char* s = arg->data.symbol.name;
        if (s == SYM_Indeterminate) return INFR_ARG_INDET;
        if (s == SYM_Infinity || s == SYM_ComplexInfinity) return INFR_ARG_INF;
        return 0;
    }
    if (head_sym(arg, SYM_Times)) {
        size_t ac = arg->data.function.arg_count;
        if (ac && sym_is(arg->data.function.args[ac - 1], SYM_Infinity)) return INFR_ARG_INF;
    }
    return 0;
}

/* ---- the table ---------------------------------------------------------- */

enum {
    R_NONE = 0,   /* no rule: the builtin decides */
    R_INDET, R_ZERO, R_INF, R_NINF, R_IINF, R_NIINF, R_CINF,
    R_IVL,        /* Interval[{-1, 1}] */
    R_IVL01,      /* Interval[{0, 1}] */
    R_IVLM10,     /* Interval[{-1, 0}] */
    R_PI_2, R_NPI_2,
    R_PARITY_INF  /* (-1)^n Infinity, n the integer first argument */
};

/* Conditions on the arguments other than the infinite one. */
enum {
    C_NONE = 0,   /* unary */
    C_NUM,        /* every other argument a finite numeric quantity */
    C_REAL,       /* every other argument a real number */
    C_INT         /* every other argument an integer (R_PARITY_INF) */
};

/* Direction slots of InfRow.v. */
enum { V_POS, V_NEG, V_I, V_NEG_I, V_CINF, V_OTHER, V_N };

typedef struct {
    const char* name;
    unsigned char arity, pos, cond;
    unsigned char v[V_N];
} InfRow;

/*                                          +Inf     -Inf     I Inf    -I Inf   CInf     other */
static const InfRow INF_ROWS[] = {
    { "Sin",        1, 0, C_NONE, { R_IVL,   R_IVL,   R_NONE,  R_NONE,  R_INDET, R_NONE } },
    { "Cos",        1, 0, C_NONE, { R_IVL,   R_IVL,   R_NONE,  R_NONE,  R_INDET, R_NONE } },
    { "Tan",        1, 0, C_NONE, { R_INDET, R_INDET, R_NONE,  R_NONE,  R_INDET, R_NONE } },
    { "Cot",        1, 0, C_NONE, { R_INDET, R_INDET, R_NONE,  R_NONE,  R_INDET, R_NONE } },
    { "Sec",        1, 0, C_NONE, { R_INDET, R_INDET, R_ZERO,  R_ZERO,  R_INDET, R_NONE } },
    { "Csc",        1, 0, C_NONE, { R_INDET, R_INDET, R_ZERO,  R_ZERO,  R_INDET, R_NONE } },
    { "Sinh",       1, 0, C_NONE, { R_NONE,  R_NONE,  R_NONE,  R_NONE,  R_INDET, R_NONE } },
    { "Cosh",       1, 0, C_NONE, { R_NONE,  R_NONE,  R_IVL,   R_IVL,   R_INDET, R_NONE } },
    { "Tanh",       1, 0, C_NONE, { R_NONE,  R_NONE,  R_INDET, R_INDET, R_INDET, R_NONE } },
    { "Coth",       1, 0, C_NONE, { R_NONE,  R_NONE,  R_INDET, R_INDET, R_INDET, R_NONE } },
    { "Sech",       1, 0, C_NONE, { R_NONE,  R_NONE,  R_INDET, R_INDET, R_INDET, R_NONE } },
    { "Csch",       1, 0, C_NONE, { R_NONE,  R_NONE,  R_INDET, R_INDET, R_INDET, R_NONE } },
    { "ArcSin",     1, 0, C_NONE, { R_NIINF, R_IINF,  R_NONE,  R_NIINF, R_NONE,  R_NONE } },
    { "ArcCos",     1, 0, C_NONE, { R_IINF,  R_NIINF, R_NIINF, R_IINF,  R_NONE,  R_NONE } },
    { "ArcSinh",    1, 0, C_NONE, { R_NONE,  R_NONE,  R_INF,   R_NINF,  R_NONE,  R_NONE } },
    { "ArcCosh",    1, 0, C_NONE, { R_NONE,  R_INF,   R_INF,   R_INF,   R_NONE,  R_NONE } },
    { "Log",        1, 0, C_NONE, { R_NONE,  R_INF,   R_INF,   R_INF,   R_INF,   R_INF  } },
    { "Erf",        1, 0, C_NONE, { R_NONE,  R_NONE,  R_NONE,  R_NONE,  R_INDET, R_NONE } },
    { "Erfc",       1, 0, C_NONE, { R_NONE,  R_NONE,  R_NONE,  R_NONE,  R_INDET, R_NONE } },
    { "Erfi",       1, 0, C_NONE, { R_NONE,  R_NONE,  R_NONE,  R_NONE,  R_INDET, R_NONE } },
    { "Factorial",  1, 0, C_NONE, { R_INF,   R_NONE,  R_NONE,  R_NONE,  R_NONE,  R_NONE } },
    { "EllipticK",  1, 0, C_NONE, { R_NONE,  R_ZERO,  R_ZERO,  R_ZERO,  R_NONE,  R_ZERO } },
    { "EllipticE",  1, 0, C_NONE, { R_NONE,  R_INF,   R_NONE,  R_NONE,  R_NONE,  R_NONE } },
    { "BesselJ",    2, 1, C_NUM,  { R_ZERO,  R_ZERO,  R_NONE,  R_NONE,  R_NONE,  R_NONE } },
    { "BesselY",    2, 1, C_NUM,  { R_ZERO,  R_ZERO,  R_NONE,  R_NONE,  R_NONE,  R_NONE } },
    { "BesselK",    2, 1, C_NUM,  { R_ZERO,  R_NONE,  R_NONE,  R_NONE,  R_NONE,  R_NONE } },
    { "BesselI",    2, 1, C_NUM,  { R_INF,   R_NONE,  R_NONE,  R_NONE,  R_NONE,  R_NONE } },
    { "BesselI",    2, 1, C_INT,  { R_INF,   R_PARITY_INF, R_NONE, R_NONE, R_NONE, R_NONE } },
    { "ExpIntegralE", 2, 1, C_NUM, { R_ZERO, R_NONE,  R_NONE,  R_NONE,  R_NONE,  R_NONE } },
    { "ArcTan",     2, 1, C_REAL, { R_PI_2,  R_NPI_2, R_NONE,  R_NONE,  R_NONE,  R_NONE } },
    { "ArcTan",     2, 0, C_REAL, { R_ZERO,  R_NONE,  R_NONE,  R_NONE,  R_NONE,  R_NONE } },
    /* x - Floor[x] ranges over [0, 1) however large x grows */
    { "FractionalPart", 1, 0, C_NONE, { R_IVL01, R_IVLM10, R_NONE, R_NONE, R_INDET, R_NONE } },
};
#define N_INF_ROWS (sizeof(INF_ROWS) / sizeof(INF_ROWS[0]))
static const char* g_row_sym[N_INF_ROWS];
static bool g_rows_ready = false;

static void rows_init(void) {
    for (size_t i = 0; i < N_INF_ROWS; i++) g_row_sym[i] = intern_symbol(INF_ROWS[i].name);
    g_rows_ready = true;
}

static bool finite_numeric(const Expr* e) {
    return numeric_special_arg_kind(e) == 0 && !head_sym(e, SYM_DirectedInfinity) &&
           expr_is_numeric_quantity(e);
}

static bool row_cond_holds(const InfRow* r, const Expr* call) {
    for (size_t i = 0; i < r->arity; i++) {
        if (i == r->pos) continue;
        const Expr* a = call->data.function.args[i];
        switch (r->cond) {
            case C_NUM:  if (!finite_numeric(a)) return false; break;
            case C_REAL: { int s = real_number_sign(a); if (s == 2) return false; break; }
            case C_INT:  if (a->type != EXPR_INTEGER && a->type != EXPR_BIGINT) return false; break;
            default: break;
        }
    }
    return true;
}

static Expr* build_code(int code, const Expr* call) {
    switch (code) {
        case R_INDET: return mk_indet();
        case R_ZERO:  return expr_new_integer(0);
        case R_INF:   return mk_inf();
        case R_NINF:  return mk_ninf();
        case R_IINF:  return mk2(SYM_Times, mk_unit_i(1), mk_inf());
        case R_NIINF: return mk2(SYM_Times, mk_unit_i(-1), mk_inf());
        case R_CINF:  return mk_cinf();
        case R_IVL:   return eval_and_free(mk1(SYM_Interval,
                          mk2(SYM_List, expr_new_integer(-1), expr_new_integer(1))));
        case R_IVL01: return eval_and_free(mk1(SYM_Interval,
                          mk2(SYM_List, expr_new_integer(0), expr_new_integer(1))));
        case R_IVLM10: return eval_and_free(mk1(SYM_Interval,
                          mk2(SYM_List, expr_new_integer(-1), expr_new_integer(0))));
        case R_PI_2:  return mk2(SYM_Times, make_rational(1, 2), expr_new_symbol(SYM_Pi));
        case R_NPI_2: return mk2(SYM_Times, make_rational(-1, 2), expr_new_symbol(SYM_Pi));
        case R_PARITY_INF: {
            const Expr* n = call->data.function.args[0];
            bool odd = n->type == EXPR_INTEGER ? (n->data.integer & 1) != 0
                                               : mpz_odd_p(n->data.bigint) != 0;
            return odd ? mk_ninf() : mk_inf();
        }
        default: return NULL;
    }
}

static int dir_slot(InfDir d) {
    switch (d) {
        case INF_DIR_POS: return V_POS;
        case INF_DIR_NEG: return V_NEG;
        case INF_DIR_I: return V_I;
        case INF_DIR_NEG_I: return V_NEG_I;
        case INF_DIR_CINF: return V_CINF;
        case INF_DIR_OTHER: return V_OTHER;
        default: return -1;
    }
}

/* ---- heads that read the direction ------------------------------------- */

/* Sign of the real or imaginary part of the direction `u` (owned): -1/0/1,
 * or 2 when undecided. */
static int part_sign(const char* part, Expr* u) {
    Expr* s = eval_and_free(mk1(SYM_Sign, mk1(part, u)));
    int r = (s && s->type == EXPR_INTEGER) ? (int)s->data.integer : 2;
    expr_free(s);
    return r;
}

static Expr* part_of_infinity(const char* part, InfDir d, const Expr* arg) {
    bool re = (part == SYM_Re);
    switch (d) {
        case INF_DIR_POS:   return re ? mk_inf() : expr_new_integer(0);
        case INF_DIR_NEG:   return re ? mk_ninf() : expr_new_integer(0);
        case INF_DIR_I:     return re ? expr_new_integer(0) : mk_inf();
        case INF_DIR_NEG_I: return re ? expr_new_integer(0) : mk_ninf();
        case INF_DIR_CINF:  return mk_indet();
        case INF_DIR_OTHER: {
            int s = part_sign(part, direction_expr(d, arg));
            if (s == 0) return expr_new_integer(0);
            if (s == 1) return mk_inf();
            if (s == -1) return mk_ninf();
            return NULL;
        }
        default: return NULL;
    }
}

static Expr* directional_rule(const char* h, const Expr* arg) {
    InfDir d = infinity_direction(arg, NULL);
    if (d == INF_DIR_NONE) return NULL;
    if (h == SYM_Abs) return mk_inf();
    if (h == SYM_Re || h == SYM_Im) return part_of_infinity(h, d, arg);
    if (h == SYM_ReIm) {
        Expr* r = part_of_infinity(SYM_Re, d, arg);
        Expr* i = part_of_infinity(SYM_Im, d, arg);
        if (!r || !i) { expr_free(r); expr_free(i); return NULL; }
        return mk2(SYM_List, r, i);
    }
    if (h == SYM_Conjugate) {
        if (d == INF_DIR_CINF) return mk_cinf();
        if (d == INF_DIR_POS) return mk_inf();
        if (d == INF_DIR_NEG) return mk_ninf();
        return mk_dir_inf(eval_and_free(mk1(SYM_Conjugate, direction_expr(d, arg))));
    }
    if (h == SYM_Sign) {
        if (d == INF_DIR_CINF) return mk_indet();
        return eval_and_free(mk1(SYM_Sign, direction_expr(d, arg)));
    }
    if (h == SYM_Arg) {
        if (d == INF_DIR_CINF) return NULL;
        return eval_and_free(mk1(SYM_Arg, direction_expr(d, arg)));
    }
    return NULL;
}

/* ---- Power ---------------------------------------------------------------- */

/* A positive real base other than 1 (a number, E, or Pi). */
static bool positive_real_base(const Expr* b) {
    if (sym_is(b, SYM_E) || sym_is(b, SYM_Pi)) return true;
    if (real_number_sign(b) != 1) return false;
    if (b->type == EXPR_INTEGER && b->data.integer == 1) return false;
    if (b->type == EXPR_REAL && b->data.real == 1.0) return false;
#ifdef USE_MPFR
    if (b->type == EXPR_MPFR && mpfr_cmp_ui(b->data.mpfr, 1) == 0) return false;
#endif
    return true;
}

static Expr* power_rule(const Expr* call) {
    const Expr* b = call->data.function.args[0];
    const Expr* x = call->data.function.args[1];
    InfDir db = infinity_direction(b, NULL);
    InfDir dx = infinity_direction(x, NULL);
    if (db != INF_DIR_NONE && dx == INF_DIR_NONE) {
        /* (u Infinity)^p, p a nonzero exact real: u^p Infinity for p > 0
         * (Sqrt[-Infinity] = I Infinity), 0 for p < 0. Infinity^p itself and
         * an inexact or symbolic exponent stay with the builtin. */
        if (x->type != EXPR_INTEGER && x->type != EXPR_BIGINT && !head_sym(x, SYM_Rational))
            return NULL;
        int s = real_number_sign(x);
        if (s == 0 || s == 2) return NULL;
        if (s < 0) return expr_new_integer(0);
        if (db == INF_DIR_CINF) return mk_cinf();
        if (db == INF_DIR_POS) return NULL;
        Expr* u = direction_expr(db, b);
        return mk_dir_inf(eval_and_free(mk2(SYM_Power, u, expr_copy((Expr*)x))));
    }
    if (db == INF_DIR_NONE && (dx == INF_DIR_CINF || dx == INF_DIR_I || dx == INF_DIR_NEG_I)) {
        /* b^z for a positive base b != 1: |b^z| = b^Re[z] and the phase is
         * Im[z] Log[b], so along the imaginary axis b^z circles forever and
         * at ComplexInfinity it has an essential singularity -- no value.
         * (Exp[I Infinity], E^ComplexInfinity, 2^ComplexInfinity.) */
        if (positive_real_base(b)) return mk_indet();
    }
    return NULL;
}

/* ---- polynomials of explicit degree --------------------------------------- */

/* p_n(x) for the classical families: p_0 = 1, and for n >= 1 the value at
 * u Infinity is the leading term lc (u Infinity)^n = sign(lc) u^n Infinity.
 * Every leading coefficient is positive except LaguerreL's (-1)^n / n!. */
static Expr* poly_rule(const char* h, const Expr* call) {
    if (call->data.function.arg_count != 2) return NULL;
    const Expr* n = call->data.function.args[0];
    const Expr* x = call->data.function.args[1];
    if (n->type != EXPR_INTEGER || n->data.integer < 0 || n->data.integer > 100000) return NULL;
    InfDir d = infinity_direction(x, NULL);
    if (d == INF_DIR_NONE) return NULL;
    int64_t deg = n->data.integer;
    if (deg == 0) return expr_new_integer(1);
    if (d == INF_DIR_CINF) return mk_cinf();
    int lc_sign = (h == SYM_LaguerreL && (deg & 1)) ? -1 : 1;
    Expr* u = direction_expr(d, x);
    Expr* un = eval_and_free(mk2(SYM_Power, u, expr_new_integer(deg)));
    if (lc_sign < 0) un = eval_and_free(mk2(SYM_Times, expr_new_integer(-1), un));
    return mk_dir_inf(un);
}

static bool is_poly_head(const char* h) {
    return h == SYM_ChebyshevT || h == SYM_ChebyshevU || h == SYM_HermiteH ||
           h == SYM_LegendreP || h == SYM_LaguerreL || h == SYM_BernoulliB;
}

/* HurwitzZeta[Infinity, a] = lim Sum[(k + a)^-s]: every term with k + a > 1
 * vanishes, a term equal to 1 stays 1, and one in (0, 1) blows up -- so for
 * a real number a: 0 when a > 1, 1 when a = 1, Infinity when 0 < a < 1
 * (audit HurwitzZeta#2). */
static Expr* hurwitz_rule(const Expr* call) {
    if (infinity_direction(call->data.function.args[0], NULL) != INF_DIR_POS) return NULL;
    const Expr* a = call->data.function.args[1];
    if (real_number_sign(a) != 1) return NULL;
    Expr* d = eval_and_free(mk2(SYM_Plus, expr_copy((Expr*)a), expr_new_integer(-1)));
    int s = real_number_sign(d);
    expr_free(d);
    if (s == 1) return expr_new_integer(0);
    if (s == 0) return expr_new_integer(1);
    if (s == -1) return mk_inf();
    return NULL;
}

/* BernoulliB and the other polynomial heads that are not NumericFunction
 * reach the leading-term rule through this (infinity_rules.h). */
Expr* infinity_poly_value(const Expr* call) {
    if (!call || call->type != EXPR_FUNCTION || call->data.function.head->type != EXPR_SYMBOL)
        return NULL;
    const char* h = call->data.function.head->data.symbol.name;
    return is_poly_head(h) ? poly_rule(h, call) : NULL;
}

/* ---- entry -------------------------------------------------------------- */

Expr* numeric_special_value_rule(const Expr* call, bool has_indet) {
    if (!call || call->type != EXPR_FUNCTION || call->data.function.head->type != EXPR_SYMBOL)
        return NULL;
    const char* h = call->data.function.head->data.symbol.name;
    size_t argc = call->data.function.arg_count;
    if (argc == 0) return NULL;

    if (has_indet) {
        if (h == SYM_ReIm || h == SYM_QuotientRemainder)
            return mk2(SYM_List, mk_indet(), mk_indet());
        return mk_indet();
    }

    if (argc == 1 && (h == SYM_Abs || h == SYM_Re || h == SYM_Im || h == SYM_ReIm ||
                      h == SYM_Sign || h == SYM_Arg || h == SYM_Conjugate))
        return directional_rule(h, call->data.function.args[0]);
    if (h == SYM_Power && argc == 2) return power_rule(call);
    if (is_poly_head(h)) return poly_rule(h, call);
    if (argc == 2 && h == intern_symbol("HurwitzZeta")) return hurwitz_rule(call);
    if (argc == 3 && h == intern_symbol("QPochhammer") &&
        infinity_direction(call->data.function.args[2], NULL) == INF_DIR_POS)
        /* (a; q)_Infinity is the two-argument QPochhammer (QPochhammer#2) */
        return mk2(intern_symbol("QPochhammer"), expr_copy(call->data.function.args[0]),
                   expr_copy(call->data.function.args[1]));

    if (!g_rows_ready) rows_init();
    for (size_t i = 0; i < N_INF_ROWS; i++) {
        if (g_row_sym[i] != h) continue;
        const InfRow* r = &INF_ROWS[i];
        if (r->arity != argc) continue;
        int slot = dir_slot(infinity_direction(call->data.function.args[r->pos], NULL));
        if (slot < 0 || r->v[slot] == R_NONE) continue;
        if (!row_cond_holds(r, call)) continue;
        Expr* out = build_code(r->v[slot], call);
        if (out) return out;
    }
    return NULL;
}

/* ---- DirectedInfinity ------------------------------------------------------
 *
 * Mathilda's canonical infinities are Infinity, ComplexInfinity and
 * c Infinity (c a unit number), so DirectedInfinity[z] with a numeric z
 * evaluates to that form: DirectedInfinity[1] is Infinity (and SameQ to it),
 * DirectedInfinity[-1] is -Infinity, DirectedInfinity[2 I] is I Infinity,
 * DirectedInfinity[1 + I] is (1 + I)/Sqrt[2] Infinity, and DirectedInfinity[]
 * and DirectedInfinity[0] are ComplexInfinity. A symbolic direction stays. */
Expr* builtin_directedinfinity(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return mk_cinf();
    if (n != 1) return NULL;
    Expr* z = res->data.function.args[0];
    if (sym_is(z, SYM_Indeterminate)) return mk_indet();
    if (!finite_numeric(z)) return NULL;
    Expr* s = eval_and_free(mk1(SYM_Sign, expr_copy(z)));
    if (!s) return NULL;
    if (s->type == EXPR_INTEGER) {
        int64_t v = s->data.integer;
        expr_free(s);
        if (v == 0) return mk_cinf();
        return v > 0 ? mk_inf() : mk_ninf();
    }
    if (head_sym(s, SYM_Sign)) { expr_free(s); return NULL; }   /* undecided */
    return mk2(SYM_Times, s, mk_inf());
}

/* ---- poles hit by an inexact argument --------------------------------------
 *
 * A machine or MPFR kernel that lands exactly on a pole -- Cot[0.],
 * ArcTanh[1.], Factorial[-1.], Csc[0.`20] -- produces inf (or NaN), which the
 * evaluator's net would report as Overflow[] (General::ovfl) though nothing
 * overflowed. A pole sits at a simple exact point, and the exact call there
 * says what the infinity is (ComplexInfinity, Infinity, -I Infinity), so the
 * inexact call answers the same: the result agrees with f[exact] at the
 * point the inexact argument names. The exact re-evaluation is confined to
 * arguments that are small dyadic rationals (|x| <= 2^20, denominator <= 2^8),
 * so it is cheap and can never start a huge exact computation (an overflow
 * past MPFR's range needs a far larger argument, and stays Overflow[]). */

#define POLE_MAX_MAG 1048576.0   /* 2^20 */
#define POLE_DEN 256.0           /* 2^8 */

/* The exact dyadic value of a finite inexact leaf, or NULL. Sets *inexact. */
static Expr* pole_exact_leaf(const Expr* e, bool* inexact) {
    double v;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return expr_copy((Expr*)e);
    if (head_sym(e, SYM_Rational)) return expr_copy((Expr*)e);
    if (e->type == EXPR_REAL) v = e->data.real;
#ifdef USE_MPFR
    else if (e->type == EXPR_MPFR) {
        if (!mpfr_number_p(e->data.mpfr)) return NULL;
        v = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
        if (mpfr_cmp_d(e->data.mpfr, v) != 0) return NULL;
    }
#endif
    else return NULL;
    if (!(fabs(v) <= POLE_MAX_MAG)) return NULL;           /* also rejects NaN/inf */
    double s = v * POLE_DEN;
    if (s != floor(s)) return NULL;
    *inexact = true;
    return make_rational((int64_t)s, (int64_t)POLE_DEN);
}

static Expr* pole_exact_arg(const Expr* e, bool* inexact) {
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) {
        Expr* r = pole_exact_leaf(re, inexact);
        Expr* i = r ? pole_exact_leaf(im, inexact) : NULL;
        if (!r || !i) { expr_free(r); expr_free(i); return NULL; }
        return eval_and_free(mk2(SYM_Plus, r, mk2(SYM_Times, i, mk_unit_i(1))));
    }
    return pole_exact_leaf(e, inexact);
}

/* Cheap pre-scan, no allocation: every argument a number (or a Complex of
 * them) and at least one inexact. A symbolic call is rejected at its first
 * argument's type test. */
static bool pole_leaf_ok(const Expr* e, bool* inexact) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (e->type == EXPR_REAL
#ifdef USE_MPFR
        || e->type == EXPR_MPFR
#endif
        ) { *inexact = true; return true; }
    return head_sym(e, SYM_Rational);
}
static bool pole_args_ok(const Expr* call) {
    bool inexact = false;
    for (size_t i = 0; i < call->data.function.arg_count; i++) {
        const Expr* a = call->data.function.args[i];
        Expr *re, *im;
        if (!a) return false;   /* an argument the builtin took over */
        if (pole_leaf_ok(a, &inexact)) continue;
        if (is_complex((Expr*)a, &re, &im) && pole_leaf_ok(re, &inexact) &&
            pole_leaf_ok(im, &inexact))
            continue;
        return false;
    }
    return inexact;
}

static int g_pole_depth = 0;

Expr* numeric_pole_value(const Expr* call) {
    if (!call || call->type != EXPR_FUNCTION || g_pole_depth > 0) return NULL;
    /* Power keeps its own division-by-zero semantics (1/0. is ComplexInfinity
     * with Power::infy, and an internal probe's muted 0.^-1 must not turn into
     * a certified pole: the zero tester reads a vanished MPFR denominator as
     * lost precision and climbs its ladder). */
    if (call->data.function.head->type == EXPR_SYMBOL &&
        call->data.function.head->data.symbol.name == SYM_Power)
        return NULL;
    size_t n = call->data.function.arg_count;
    if (n == 0 || n > 8 || !pole_args_ok(call)) return NULL;
    Expr* args[8];
    bool inexact = false;
    for (size_t i = 0; i < n; i++) {
        args[i] = pole_exact_arg(call->data.function.args[i], &inexact);
        if (!args[i]) {
            for (size_t j = 0; j < i; j++) expr_free(args[j]);
            return NULL;
        }
    }
    if (!inexact) {
        for (size_t j = 0; j < n; j++) expr_free(args[j]);
        return NULL;
    }
    Expr* exact = expr_new_function(expr_copy(call->data.function.head), args, n);
    g_pole_depth++;
    g_arith_warnings_muted++;
    Expr* r = eval_and_free(exact);
    g_arith_warnings_muted--;
    g_pole_depth--;
    if (r && infinity_direction(r, NULL) != INF_DIR_NONE) return r;
    expr_free(r);
    return NULL;
}
