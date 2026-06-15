/*
 * nresidue.c — NResidue[expr, {z, z0}, opts]
 *
 * Numerically finds the residue of `expr` at z = z0 by integrating around a
 * small circle in the complex plane (the periodic-trapezoidal Cauchy
 * integral; see quadrature.h). Unlike the symbolic Residue, this works for
 * essential singularities (Exp[1/x], Sin[1/x]) where the power series does
 * not exist.
 *
 *   NResidue[expr, {z, z0}]               residue near z = z0
 *   NResidue[{e1, e2, ...}, {z, z0}]      threads element-wise over arg 1
 *
 * Options (trailing Rule[...] in any order):
 *   Radius           -> r | Automatic   contour radius (default 1/100)
 *   WorkingPrecision -> MachinePrecision | digits
 *   PrecisionGoal    -> digits | Automatic | Infinity
 *   MaxRecursion     -> n               max N-doublings (default 10)
 *   Method           -> "Trapezoidal"   (only method in this version)
 *
 * The sampler binds the variable to each complex sample point (Block-style,
 * via temporary OwnValues) and numerically evaluates the integrand, so the
 * user's global symbol table is left untouched. NResidue cannot tell a tiny
 * spurious residual from a true zero — Chop the result when needed.
 *
 * Memory: receives `res` owned by the evaluator. Returns a fresh Expr* on
 * success or NULL (unevaluated). Never frees `res`. All temporary OwnValues
 * are removed before returning, on every path.
 */

#include "nresidue.h"
#include "quadrature.h"

#include <complex.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#  include <mpfr.h>
#endif

#include "arithmetic.h"   /* is_complex, make_complex, is_rational */
#include "attr.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ *
 *  Options                                                            *
 * ------------------------------------------------------------------ */

typedef struct {
    bool   radius_auto;        /* Radius -> Automatic                  */
    double radius;             /* fixed radius (default 1/100)         */
    bool   prec_mpfr;          /* WorkingPrecision selects MPFR        */
    long   bits;               /* MPFR working precision in bits       */
    int    max_recursion;      /* default 10                           */
    double prec_goal_digits;   /* < 0 => derive from WorkingPrecision  */
} NrOpts;

static void nr_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    fprintf(stderr, "NResidue::%s: ", tag);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* ------------------------------------------------------------------ *
 *  Numeric leaf <-> C value helpers                                   *
 * ------------------------------------------------------------------ */

static bool nr_to_double_real(Expr* e, double* out) {
    if (!e) return false;
    int64_t rn, rd;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;   return true;
        case EXPR_REAL:    *out = e->data.real;              return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        case EXPR_FUNCTION:
            if (is_rational(e, &rn, &rd)) { *out = (double)rn / (double)rd; return true; }
            return false;
        default: return false;
    }
}

static bool nr_to_complex(Expr* e, double _Complex* out) {
    if (!e) return false;
    double rv;
    if (nr_to_double_real(e, &rv)) { *out = rv + 0.0 * I; return true; }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        double r, i;
        if (nr_to_double_real(re, &r) && nr_to_double_real(im, &i)) {
            *out = r + i * I;
            return true;
        }
    }
    return false;
}

static Expr* nr_from_complex_d(double _Complex c) {
    double r = creal(c), i = cimag(c);
    if (i == 0.0) return expr_new_real(r);
    return make_complex(expr_new_real(r), expr_new_real(i));
}

#ifdef USE_MPFR
static Expr* nr_from_complex_mpfr(const mpfr_t re, const mpfr_t im) {
    if (mpfr_zero_p(im)) return expr_new_mpfr_copy(re);
    return make_complex(expr_new_mpfr_copy(re), expr_new_mpfr_copy(im));
}
#endif

/* ------------------------------------------------------------------ *
 *  Block-style variable binding (mirrors FindRoot)                    *
 * ------------------------------------------------------------------ */

typedef struct {
    const char* name;      /* interned                           */
    Rule*       saved_own; /* prior OwnValue chain (borrowed)    */
    uint32_t    saved_attrs;
    bool        valid;
} NrBind;

static void nr_bind_snapshot(NrBind* b, const char* name) {
    SymbolDef* def = symtab_get_def(name);
    b->name = name;
    b->saved_own = def->own_values;
    b->saved_attrs = def->attributes;
    def->own_values = NULL;
    b->valid = true;
}

static void nr_bind_clear_temp(SymbolDef* def) {
    Rule* curr = def->own_values;
    while (curr) {
        Rule* next = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        curr = next;
    }
    def->own_values = NULL;
}

static void nr_bind_set(NrBind* b, Expr* value) {
    SymbolDef* def = symtab_get_def(b->name);
    nr_bind_clear_temp(def);
    Expr* sym = expr_new_symbol(b->name);
    symtab_add_own_value(b->name, sym, value);
    expr_free(sym);
}

static void nr_bind_restore(NrBind* b) {
    if (!b->valid) return;
    SymbolDef* def = symtab_get_def(b->name);
    nr_bind_clear_temp(def);
    def->own_values = b->saved_own;
    def->attributes = b->saved_attrs;
    b->valid = false;
    eval_clock_bump();
}

/* ------------------------------------------------------------------ *
 *  Sampler context + callbacks                                        *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr*       f;        /* borrowed integrand (variable free)  */
    NrBind*     bind;     /* the active variable binding         */
    NumericSpec spec;     /* machine or MPFR                     */
    long        bits;     /* MPFR bits (0 => machine)            */
} NrCtx;

/* Evaluate f after binding the variable to `value` (consumed). Returns the
 * numericalized result (caller frees) or NULL. */
static Expr* nr_eval_at(NrCtx* c, Expr* value) {
    nr_bind_set(c->bind, value);   /* takes a copy internally */
    expr_free(value);
    eval_clock_bump();
    Expr* raw = eval_and_free(expr_copy(c->f));
    if (!raw) return NULL;
    Expr* num = numericalize(raw, c->spec);
    expr_free(raw);
    return num;
}

static bool nr_sample_machine(void* vctx, double _Complex z, double _Complex* out) {
    NrCtx* c = (NrCtx*)vctx;
    Expr* zv = nr_from_complex_d(z);
    Expr* num = nr_eval_at(c, zv);
    if (!num) return false;
    bool ok = nr_to_complex(num, out);
    expr_free(num);
    return ok;
}

#ifdef USE_MPFR
static bool nr_sample_mpfr(void* vctx, const mpfr_t z_re, const mpfr_t z_im,
                           mpfr_t out_re, mpfr_t out_im) {
    NrCtx* c = (NrCtx*)vctx;
    Expr* zv = nr_from_complex_mpfr(z_re, z_im);
    Expr* num = nr_eval_at(c, zv);
    if (!num) return false;
    bool inexact;
    bool ok = get_approx_mpfr(num, out_re, out_im, &inexact);
    expr_free(num);
    return ok;
}
#endif

/* ------------------------------------------------------------------ *
 *  Option parsing                                                     *
 * ------------------------------------------------------------------ */

static bool nr_is_known_option(const char* s) {
    return s == SYM_Radius
        || s == SYM_WorkingPrecision
        || s == SYM_PrecisionGoal
        || s == SYM_MaxRecursion
        || s == SYM_Method;
}

static bool nr_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && nr_is_known_option(lhs->data.symbol);
}

static bool nr_parse_working_precision(Expr* val, bool* mpfr, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol == SYM_MachinePrecision) {
        *mpfr = false; *bits = 0; return true;
    }
    double digits;
    if (!nr_to_double_real(val, &digits) || digits <= 0.0) return false;
#ifdef USE_MPFR
    if (digits <= NUMERIC_MACHINE_PRECISION_DIGITS) { *mpfr = false; *bits = 0; }
    else { *mpfr = true; *bits = numeric_digits_to_bits(digits); }
#else
    *mpfr = false; *bits = 0;
#endif
    return true;
}

static bool nr_parse_goal(Expr* val, double* digits_out) {
    if (val->type == EXPR_SYMBOL) {
        if (val->data.symbol == SYM_Automatic) { *digits_out = -1.0;      return true; }
        if (val->data.symbol == SYM_Infinity)  { *digits_out = INFINITY;  return true; }
        return false;
    }
    return nr_to_double_real(val, digits_out);
}

static bool nr_apply_option(Expr* rule, NrOpts* o) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol;

    if (name == SYM_Radius) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) {
            o->radius_auto = true; return true;
        }
        double r;
        if (!nr_to_double_real(rhs, &r) || r <= 0.0) {
            nr_warn("badopt", "Radius must be a positive number or Automatic");
            return false;
        }
        o->radius_auto = false; o->radius = r; return true;
    }
    if (name == SYM_WorkingPrecision) {
        if (!nr_parse_working_precision(rhs, &o->prec_mpfr, &o->bits)) {
            nr_warn("badopt", "invalid WorkingPrecision value");
            return false;
        }
        return true;
    }
    if (name == SYM_PrecisionGoal) return nr_parse_goal(rhs, &o->prec_goal_digits);
    if (name == SYM_MaxRecursion) {
        if (rhs->type == EXPR_INTEGER && rhs->data.integer >= 0) {
            o->max_recursion = (int)rhs->data.integer; return true;
        }
        nr_warn("badopt", "MaxRecursion must be a non-negative integer");
        return false;
    }
    if (name == SYM_Method) {
        /* Only the trapezoidal rule is implemented; accept Automatic and the
         * "Trapezoidal" name, warn (but proceed) on anything else. */
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) return true;
        if (rhs->type == EXPR_STRING && strcmp(rhs->data.string, "Trapezoidal") == 0) return true;
        if (rhs->type == EXPR_SYMBOL && strcmp(rhs->data.symbol, "Trapezoidal") == 0) return true;
        nr_warn("badmeth", "only Method -> \"Trapezoidal\" is supported; using it");
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Manual list-threading over arg 1                                   *
 * ------------------------------------------------------------------ */

static Expr* nr_thread_over_list(Expr* res) {
    Expr* lst = res->data.function.args[0];
    size_t n = lst->data.function.arg_count;
    size_t argc = res->data.function.arg_count;
    Expr** items = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr** ca = malloc(sizeof(Expr*) * argc);
        ca[0] = expr_copy(lst->data.function.args[i]);
        for (size_t j = 1; j < argc; j++) ca[j] = expr_copy(res->data.function.args[j]);
        items[i] = expr_new_function(expr_new_symbol(SYM_NResidue), ca, argc);
        free(ca);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return eval_and_free(out);
}

/* ------------------------------------------------------------------ *
 *  Result-status reporting shared by both precision paths             *
 * ------------------------------------------------------------------ */

static Expr* nr_report(QdStatus status, int n_used, Expr* value /*owned*/) {
    switch (status) {
        case QD_OK:
            return value;
        case QD_NOCONV:
            nr_warn("ncvi", "failed to converge to prescribed accuracy after %d "
                            "sample points; result may be inaccurate", n_used);
            return value;
        case QD_BRANCHCUT:
            nr_warn("bcut", "the integrand does not appear analytic on the contour "
                            "(possible branch cut); result is unreliable");
            return value;
        case QD_NONNUMERIC:
        default:
            nr_warn("nnum", "the integrand could not be evaluated to a number on "
                            "the integration contour");
            expr_free(value);
            return NULL;
    }
}

/* ------------------------------------------------------------------ *
 *  Entry point                                                        *
 * ------------------------------------------------------------------ */

Expr* builtin_nresidue(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    /* Manual threading: if arg 1 is a List, thread NResidue over it and keep
     * the {z, z0} spec + options fixed. (NResidue is intentionally NOT
     * ATTR_LISTABLE: generic threading would split the {z, z0} spec.) */
    Expr* arg0 = res->data.function.args[0];
    if (arg0->type == EXPR_FUNCTION
        && arg0->data.function.head->type == EXPR_SYMBOL
        && arg0->data.function.head->data.symbol == SYM_List) {
        return nr_thread_over_list(res);
    }

    /* Peel trailing options; exactly two positional args must remain. */
    size_t pos_end = argc;
    while (pos_end > 0 && nr_is_option_arg(res->data.function.args[pos_end - 1])) pos_end--;
    for (size_t i = pos_end; i < argc; i++) {
        if (!nr_is_option_arg(res->data.function.args[i])) {
            nr_warn("badopt", "unrecognised option in trailing position");
            return NULL;
        }
    }
    if (pos_end != 2) return NULL;

    NrOpts opts;
    opts.radius_auto = false;
    opts.radius = 0.01;            /* 1/100 default */
    opts.prec_mpfr = false;
    opts.bits = 0;
    opts.max_recursion = 10;
    opts.prec_goal_digits = -1.0;  /* derive */
    for (size_t i = pos_end; i < argc; i++) {
        if (!nr_apply_option(res->data.function.args[i], &opts)) return NULL;
    }

    /* Parse the {z, z0} spec. */
    Expr* spec = res->data.function.args[1];
    if (spec->type != EXPR_FUNCTION
        || spec->data.function.head->type != EXPR_SYMBOL
        || spec->data.function.head->data.symbol != SYM_List
        || spec->data.function.arg_count != 2) {
        nr_warn("ivar", "second argument must be {z, z0}");
        return NULL;
    }
    Expr* var = spec->data.function.args[0];
    if (var->type != EXPR_SYMBOL) {
        nr_warn("ivar", "the variable in {z, z0} must be a symbol");
        return NULL;
    }

    /* Working precision / goal. */
    double wp_digits = opts.prec_mpfr ? numeric_bits_to_digits(opts.bits)
                                      : NUMERIC_MACHINE_PRECISION_DIGITS;
    double prec_goal = opts.prec_goal_digits;
    if (prec_goal < 0.0)        prec_goal = wp_digits - 2.0;   /* two guard digits */
    if (prec_goal > wp_digits)  prec_goal = wp_digits;
    if (prec_goal < 1.0)        prec_goal = 1.0;

    NumericSpec spec_num;
#ifdef USE_MPFR
    if (opts.prec_mpfr) { spec_num.mode = NUMERIC_MODE_MPFR; spec_num.bits = opts.bits;
                          spec_num.preserve_inexact = false; }
    else                  spec_num = numeric_machine_spec();
#else
    spec_num = numeric_machine_spec();
#endif

    /* Numericalise z0. */
    Expr* z0_raw = eval_and_free(expr_copy(spec->data.function.args[1]));
    Expr* z0 = z0_raw ? numericalize(z0_raw, spec_num) : NULL;
    expr_free(z0_raw);
    if (!z0) { nr_warn("nnum", "z0 is not numeric"); return NULL; }

    NrBind bind;
    nr_bind_snapshot(&bind, var->data.symbol);

    NrCtx ctx;
    ctx.f = res->data.function.args[0];
    ctx.bind = &bind;
    ctx.spec = spec_num;
    ctx.bits = opts.bits;

    double radius_arg = opts.radius_auto ? -1.0 : opts.radius;
    Expr* result = NULL;

#ifdef USE_MPFR
    if (opts.prec_mpfr) {
        mpfr_t z0r, z0i, outr, outi;
        mpfr_init2(z0r, opts.bits); mpfr_init2(z0i, opts.bits);
        mpfr_init2(outr, opts.bits); mpfr_init2(outi, opts.bits);
        bool inexact;
        if (!get_approx_mpfr(z0, z0r, z0i, &inexact)) {
            nr_warn("nnum", "z0 is not numeric");
        } else {
            QdResultMpfr R = qd_contour_residue_mpfr(nr_sample_mpfr, &ctx,
                                                     z0r, z0i, radius_arg, opts.bits,
                                                     prec_goal, opts.max_recursion,
                                                     outr, outi);
            Expr* value = (R.status == QD_NONNUMERIC) ? NULL
                                                      : nr_from_complex_mpfr(outr, outi);
            result = nr_report(R.status, R.n_used, value);
        }
        mpfr_clears(z0r, z0i, outr, outi, (mpfr_ptr)0);
    } else
#endif
    {
        double _Complex z0c;
        if (!nr_to_complex(z0, &z0c)) {
            nr_warn("nnum", "z0 is not numeric");
        } else {
            QdResultMachine R = qd_contour_residue_machine(nr_sample_machine, &ctx,
                                                           z0c, radius_arg,
                                                           prec_goal, opts.max_recursion);
            Expr* value = (R.status == QD_NONNUMERIC) ? NULL : nr_from_complex_d(R.value);
            result = nr_report(R.status, R.n_used, value);
        }
    }

    nr_bind_restore(&bind);
    expr_free(z0);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void nresidue_init(void) {
    symtab_add_builtin("NResidue", builtin_nresidue);
    /* Protected only. List-threading is manual (see builtin_nresidue), so
     * ATTR_LISTABLE is deliberately NOT set — it would split the {z, z0}
     * spec across two bogus calls. */
    symtab_get_def("NResidue")->attributes |= ATTR_PROTECTED;
}
