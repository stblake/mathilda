/*
 * sum_logrational.c -- Sum`LogRational: log-of-rational telescoping sums.
 *
 * Sums whose summand is a rational function plus a log of a rational function,
 *
 *   term(k) = P(k) + Log R(k),
 *
 * where P and Log R may individually diverge (like +c/k and -c/k) but combine
 * to a convergent series.  The canonical example:
 *
 *   Sum[1/k + Log[(k-1)/k], {k, 2, Inf}] = EulerGamma - 1,
 *
 * which unblocks Product[E^(1/k)(1-1/k), {k,2,Inf}] = E^(EulerGamma-1) via the
 * Product`LogSum bridge.
 *
 * Method (rigorous).  Write P via Apart as sum_i c_i/(k-rho_i)^{m_i}, and
 * R(k) = L * prod_num (k-a)^{e} / prod_den (k-b)^{e} (roots over Q).  Using
 *   sum_{k=imin}^N Log(k-rho)   = LogGamma(N+1-rho) - LogGamma(imin-rho),
 *   sum_{k=imin}^N c/(k-rho)    = c(PolyGamma(0,N+1-rho) - PolyGamma(0,imin-rho)),
 * the N-dependent parts cancel (Stirling / digamma asymptotics) precisely when
 *   L == 1,  #num roots == #den roots,  and
 *   sum c_i  -  sum_num e*a  +  sum_den e*b  == 0    (the 1/k coefficient),
 * leaving the exact finite value
 *
 *   V =  sum_{m>=2} c_i Zeta[m_i, imin-rho_i]           (convergent pole parts)
 *      - sum_{m==1} c_i PolyGamma[0, imin-rho_i]         (simple poles)
 *      - sum_num e * LogGamma[imin-a]                    (numerator log-roots)
 *      + sum_den e * LogGamma[imin-b].                   (denominator log-roots)
 *
 * Divergent inputs (nonzero polynomial part, L != 1, unequal counts, failed 1/k
 * cancellation) and non-rational log arguments / irrational roots fall through
 * (NULL).
 *
 * Memory contract: takes ownership of res but must not free it; returns an owned
 * closed form or NULL to fall through.
 */

#include "sum_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "sym_names.h"
#include "poly.h"
#include <string.h>
#include <stdlib.h>

/* Reused from the Product module (product_util.c); always linked in. */
bool prod_linear_factors(Expr* e, Expr* var,
                         Expr** lead_out,
                         Expr*** roots_out, int** mults_out, size_t* n_out,
                         bool* all_linear_out);

static bool is_infinity_sym(const Expr* e) {
    return e->type == EXPR_SYMBOL && e->data.symbol == SYM_Infinity;
}

static bool head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, name) == 0;
}

/* True if e contains a Log[...] anywhere. */
static bool has_log(const Expr* e) {
    if (!e) return false;
    if (head_is(e, "Log")) return true;
    if (e->type == EXPR_FUNCTION) {
        if (has_log(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (has_log(e->data.function.args[i])) return true;
    }
    return false;
}

/* e1 + e2 (adopts both). */
static Expr* add2(Expr* a, Expr* b) {
    Expr* s = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ a, b }, 2);
    Expr* r = evaluate(s);
    expr_free(s);
    return r;
}

/* imin - rho (both copied). */
static Expr* anchor_minus(Expr* imin, Expr* rho) {
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(rho) }, 2);
    Expr* s = expr_new_function(expr_new_symbol(SYM_Plus),
                  (Expr*[]){ expr_copy(imin), neg }, 2);
    Expr* r = evaluate(s);
    expr_free(s);
    return r;
}

/* PossibleZeroQ[e] == True (consumes e). */
static bool is_zero(Expr* e) {
    Expr* r = sum_eval("PossibleZeroQ", (Expr*[]){ e }, 1);
    bool yes = (r && r->type == EXPR_SYMBOL && r->data.symbol == SYM_True);
    if (r) expr_free(r);
    return yes;
}

Expr* builtin_sum_logrational(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite || !is_infinity_sym(imax)) return NULL;
    if (imin->type != EXPR_INTEGER) return NULL;
    if (!has_log(f)) return NULL;   /* pure rational belongs to Sum`Rational */

    /* --- Split term into rational part P and log-argument product Racc. --- */
    bool is_plus = head_is(f, "Plus");
    size_t nadd = is_plus ? f->data.function.arg_count : 1;

    Expr* P = expr_new_integer(0);
    Expr* Racc = expr_new_integer(1);
    bool bail = false;

    for (size_t i = 0; i < nadd && !bail; i++) {
        Expr* addend = is_plus ? f->data.function.args[i] : f;
        if (!has_log(addend)) {
            P = add2(P, expr_copy(addend));
            continue;
        }
        /* addend must be coeff * Log[arg], coeff free of var. */
        Expr* logarg = NULL;
        Expr* coeff = NULL;
        if (head_is(addend, "Log") && addend->data.function.arg_count == 1) {
            logarg = expr_copy(addend->data.function.args[0]);
            coeff = expr_new_integer(1);
        } else if (head_is(addend, "Times")) {
            Expr* c = expr_new_integer(1);
            for (size_t j = 0; j < addend->data.function.arg_count; j++) {
                Expr* fac = addend->data.function.args[j];
                if (head_is(fac, "Log") && fac->data.function.arg_count == 1
                        && !logarg) {
                    logarg = expr_copy(fac->data.function.args[0]);
                } else if (!has_log(fac) && sum_free_of(fac, var)) {
                    Expr* nc = expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ c, expr_copy(fac) }, 2); /* adopts c */
                    c = evaluate(nc);
                    expr_free(nc);
                } else {
                    /* a var-dependent non-log factor, or a second Log: unsupported */
                    expr_free(c); c = NULL; break;
                }
            }
            coeff = c;
        }
        if (!logarg || !coeff) { if (logarg) expr_free(logarg); if (coeff) expr_free(coeff); bail = true; break; }

        /* Racc *= arg^coeff. */
        Expr* p = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ logarg, coeff }, 2);   /* adopts logarg, coeff */
        Expr* nr = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ Racc, p }, 2);         /* adopts Racc, p */
        Racc = evaluate(nr);
        expr_free(nr);
    }
    if (bail) { expr_free(P); expr_free(Racc); return NULL; }

    /* --- Factor R = Racc into numerator/denominator roots over Q. --- */
    Expr* Rt = sum_eval("Together", (Expr*[]){ Racc }, 1);   /* adopts Racc */
    Expr* Rnum = sum_eval("Numerator", (Expr*[]){ expr_copy(Rt) }, 1);
    Expr* Rden = sum_eval("Denominator", (Expr*[]){ Rt }, 1);/* adopts Rt */

    Expr *leadN = NULL, *leadD = NULL;
    Expr **rootsN = NULL, **rootsD = NULL;
    int *multsN = NULL, *multsD = NULL;
    size_t nN = 0, nD = 0;
    bool linN = false, linD = false;
    bool okN = prod_linear_factors(Rnum, var, &leadN, &rootsN, &multsN, &nN, &linN);
    bool okD = prod_linear_factors(Rden, var, &leadD, &rootsD, &multsD, &nD, &linD);
    expr_free(Rnum); expr_free(Rden);

    /* Cleanup lambda (manual). */
    #define LR_FREE_ROOTS() do { \
        if (rootsN) { for (size_t _i=0;_i<nN;_i++) expr_free(rootsN[_i]); free(rootsN); } \
        if (rootsD) { for (size_t _i=0;_i<nD;_i++) expr_free(rootsD[_i]); free(rootsD); } \
        free(multsN); free(multsD); \
        if (leadN) expr_free(leadN); \
        if (leadD) expr_free(leadD); \
    } while (0)

    if (!okN || !okD || !linN || !linD) { LR_FREE_ROOTS(); expr_free(P); return NULL; }

    /* Convergence: leading ratio L == 1, and equal root counts (with mult). */
    Expr* Lratio = sum_eval("Simplify", (Expr*[]){
        expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_copy(leadN),
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ expr_copy(leadD), expr_new_integer(-1) }, 2) }, 2) }, 1);
    Expr* Lm1 = add2(Lratio, expr_new_integer(-1));
    bool L_ok = is_zero(Lm1);   /* consumes Lm1 */
    if (!L_ok) { LR_FREE_ROOTS(); expr_free(P); return NULL; }

    int cntN = 0, cntD = 0;
    for (size_t i = 0; i < nN; i++) cntN += multsN[i];
    for (size_t i = 0; i < nD; i++) cntD += multsD[i];
    if (cntN != cntD) { LR_FREE_ROOTS(); expr_free(P); return NULL; }

    /* --- Decompose the rational part via Apart and assemble V. --- */
    Expr* Pap = sum_eval("Apart", (Expr*[]){ P, expr_copy(var) }, 2);  /* adopts P */

    Expr* V = expr_new_integer(0);
    Expr* sumC = expr_new_integer(0);   /* sum of simple-pole residues */

    bool is_pplus = head_is(Pap, "Plus");
    size_t npa = is_pplus ? Pap->data.function.arg_count : 1;
    for (size_t i = 0; i < npa && !bail; i++) {
        Expr* t = is_pplus ? Pap->data.function.args[i] : Pap;
        /* Skip an exact zero. */
        if (t->type == EXPR_INTEGER && t->data.integer == 0) continue;

        Expr* tden = sum_eval("Denominator", (Expr*[]){ expr_copy(t) }, 1);
        int degd = get_degree_poly(tden, var);
        if (degd <= 0) {
            /* Polynomial / constant part -> the series diverges. */
            expr_free(tden); bail = true; break;
        }
        Expr* tnum = sum_eval("Numerator", (Expr*[]){ expr_copy(t) }, 1);
        if (!sum_free_of(tnum, var)) { expr_free(tnum); expr_free(tden); bail = true; break; }

        Expr *dlead = NULL, **drts = NULL; int *dmul = NULL; size_t dn = 0; bool dlin = false;
        bool dok = prod_linear_factors(tden, var, &dlead, &drts, &dmul, &dn, &dlin);
        expr_free(tden);
        if (!dok || !dlin || dn != 1) {
            if (drts) { for (size_t _i=0;_i<dn;_i++) expr_free(drts[_i]); free(drts); }
            free(dmul); if (dlead) expr_free(dlead);
            expr_free(tnum); bail = true; break;
        }
        Expr* rho = drts[0];
        int m = dmul[0];
        /* c = tnum / dlead. */
        Expr* c = sum_eval("Simplify", (Expr*[]){
            expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ tnum,   /* adopts tnum */
                           expr_new_function(expr_new_symbol(SYM_Power),
                               (Expr*[]){ expr_copy(dlead), expr_new_integer(-1) }, 2) }, 2) }, 1);
        Expr* anchor = anchor_minus(imin, rho);   /* imin - rho */

        if (m == 1) {
            /* - c * PolyGamma[0, imin-rho];  sumC += c. */
            Expr* pg = expr_new_function(expr_new_symbol("PolyGamma"),
                           (Expr*[]){ expr_new_integer(0), expr_copy(anchor) }, 2);
            Expr* contrib = expr_new_function(expr_new_symbol(SYM_Times),
                                (Expr*[]){ expr_new_integer(-1), expr_copy(c), pg }, 3); /* adopts pg */
            V = add2(V, contrib);
            sumC = add2(sumC, expr_copy(c));
        } else {
            /* + c * Zeta[m, imin-rho]. */
            Expr* z = expr_new_function(expr_new_symbol("Zeta"),
                          (Expr*[]){ expr_new_integer(m), expr_copy(anchor) }, 2);
            Expr* contrib = expr_new_function(expr_new_symbol(SYM_Times),
                                (Expr*[]){ expr_copy(c), z }, 2);   /* adopts z */
            V = add2(V, contrib);
        }
        expr_free(anchor); expr_free(c);
        expr_free(rho); free(drts); free(dmul); expr_free(dlead);
    }
    expr_free(Pap);
    if (bail) { expr_free(V); expr_free(sumC); LR_FREE_ROOTS(); return NULL; }

    /* 1/k cancellation: sumC - sum_num e*a + sum_den e*b == 0. */
    Expr* onek = expr_copy(sumC);
    for (size_t i = 0; i < nN; i++) {
        Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_integer(-multsN[i]), expr_copy(rootsN[i]) }, 2);
        onek = add2(onek, t);
    }
    for (size_t i = 0; i < nD; i++) {
        Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_integer(multsD[i]), expr_copy(rootsD[i]) }, 2);
        onek = add2(onek, t);
    }
    bool conv = is_zero(onek);   /* consumes onek */
    expr_free(sumC);
    if (!conv) { expr_free(V); LR_FREE_ROOTS(); return NULL; }

    /* Log-root contributions: -sum_num e*LogGamma[imin-a] + sum_den e*LogGamma[imin-b]. */
    for (size_t i = 0; i < nN; i++) {
        Expr* anchor = anchor_minus(imin, rootsN[i]);
        Expr* lg = expr_new_function(expr_new_symbol("LogGamma"),
                       (Expr*[]){ anchor }, 1);   /* adopts anchor */
        Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_integer(-multsN[i]), lg }, 2);  /* adopts lg */
        V = add2(V, t);
    }
    for (size_t i = 0; i < nD; i++) {
        Expr* anchor = anchor_minus(imin, rootsD[i]);
        Expr* lg = expr_new_function(expr_new_symbol("LogGamma"),
                       (Expr*[]){ anchor }, 1);
        Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_integer(multsD[i]), lg }, 2);
        V = add2(V, t);
    }

    LR_FREE_ROOTS();
    #undef LR_FREE_ROOTS

    Expr* out = sum_eval("Simplify", (Expr*[]){ V }, 1);   /* adopts V */
    return out;
}

void sum_logrational_init(void) {
    symtab_add_builtin("Sum`LogRational", builtin_sum_logrational);
    symtab_get_def("Sum`LogRational")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Sum`LogRational",
        "Sum`LogRational[f, i, imin, Infinity] evaluates convergent sums of a "
        "rational function plus a log of a rational function (e.g. "
        "Sum[1/k+Log[(k-1)/k]] = EulerGamma-1) via matched digamma/LogGamma "
        "asymptotics. Returns unevaluated for divergent or unrecognised forms.");
}
