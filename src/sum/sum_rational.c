/*
 * sum_rational.c -- Sum`Rational: summation of rational functions.
 *
 * For a summand p(i)/q(i) that is a rational function of the index i, with the
 * convergence condition deg q >= deg p + 2, the infinite sum from a concrete
 * integer lower bound is a finite combination of Hurwitz-zeta / polygamma
 * values.  Decompose into linear partial fractions over the denominator's
 * splitting field,
 *
 *     p(i)/q(i)  =  sum_j sum_{k=1}^{m_j}  c_{j,k} / (i - rho_j)^k ,
 *
 * and sum each term from i = imin to infinity with the master identity
 * (n = i - rho ranges over imin-rho, imin-rho+1, ...):
 *
 *   k >= 2:  sum_{i>=imin} 1/(i-rho)^k  =  Zeta[k, imin - rho]
 *            ( = ((-1)^k/(k-1)!) PolyGamma[k-1, imin - rho] ).
 *   k == 1:  individually divergent; convergence forces sum_j c_{j,1} = 0, and
 *            sum_{i>=imin} c/(i-rho)  ->  -c (PolyGamma[0, imin - rho] + EulerGamma).
 *            The EulerGamma is kept per term (it cancels because sum c = 0) so
 *            the printed form matches Mathematica.
 *
 *   Sum`Rational[f, i, imin, Infinity]  ->  closed form
 *
 * Phase A handles rational poles only: it calls Apart[f, i] (no field
 * extension), which splits over Q.  A term whose denominator is an irreducible
 * higher-degree factor (complex / radical pole) is left for the extension path;
 * encountering one makes this stage fall through (NULL) so the Sum stays held.
 *
 * The infinite case needs a concrete integer lower bound.  A FINITE range with
 * exact integer bounds is summed by sr_finite (polygamma differences, no
 * convergence condition; see there).  Every other shape (indefinite, symbolic
 * finite bounds, symbolic/non-integer imin, a non-rational summand, or a
 * divergent deg q < deg p + 2 series) falls through unevaluated.
 *
 * Memory contract: builtin_sum_rational takes ownership of res but must not
 * free it; every Expr* allocated here is freed on all paths.
 */

#include "sum_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "poly.h"
#include "arithmetic.h"   /* is_infinity_sym, make_rational_mpz */
#include "sym_names.h"
#include "simp.h"         /* simp_default_complexity */
#include "rationalize.h"  /* internal_contains_inexact */
#include <gmp.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

/* e is a polynomial in var? */
static bool sr_poly_in(Expr* e, Expr* var) {
    Expr* vars[1] = { var };
    return is_polynomial(e, vars, 1);
}

/* degree of e in var (expanding first); -1 if not a polynomial / zero. */
static int sr_pdeg(Expr* e, Expr* var) {
    Expr* ea[1] = { expr_copy(e) };
    Expr* ex = sum_eval("Expand", ea, 1);
    int d = get_degree_poly(ex, var);
    expr_free(ex);
    return d;
}

/* True if e is Power[base, exp] with exp a negative integer; sets *base (alias,
 * not copied) and *k = -exp. */
static bool sr_neg_power(Expr* e, Expr** base, int* k) {
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol.name != SYM_Power
        || e->data.function.arg_count != 2) return false;
    Expr* ex = e->data.function.args[1];
    if (ex->type != EXPR_INTEGER || ex->data.integer >= 0) return false;
    *base = e->data.function.args[0];
    *k = (int)(-ex->data.integer);
    return true;
}

/*
 * Count the var-dependent negative-power (pole) factors in `term`, descending
 * into nested Times (Apart's extension output is not guaranteed to flatten).
 * The last one found is recorded in *fb (alias) / *fk.  A proper partial-
 * fraction term has exactly one; a numerator factor that depends on var (e.g.
 * the `i` in i/(i^2+1)) is ignored here and recovered later via Cancel.
 */
static int sr_find_pole(Expr* e, Expr* var, Expr** fb, int* fk) {
    Expr* b; int kk;
    if (sr_neg_power(e, &b, &kk) && !sum_free_of(b, var)) { *fb = b; *fk = kk; return 1; }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Times) {
        int c = 0;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            c += sr_find_pole(e->data.function.args[i], var, fb, fk);
        return c;
    }
    return 0;
}

/* Residue coefficient of a pole term: Cancel[term * base^k], i.e. the term with
 * its denominator power divided out. */
static Expr* sr_residue(Expr* term, Expr* base, int k, Expr* var) {
    (void)var;
    Expr* basek = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_copy(base), sum_int(k) }, 2);
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ expr_copy(term), basek }, 2);
    return sum_eval("Cancel", (Expr*[]){ prod }, 1);
}

/* Build the contribution of a single pole term c0/base^k to the infinite sum.
 * base is linear in var: base = b1 var + b0, rho = -b0/b1, a = imin - rho, and
 * the effective residue is c = c0 * b1^(-k).  Returns an owned (unevaluated)
 * Expr*, or NULL if base is not linear in var (e.g. an irreducible quadratic
 * pole -- defer to the extension path). c0 is consumed; base is borrowed (used
 * only via copies here, so the caller retains ownership). */
static Expr* sr_term_contribution(Expr* c0, Expr* base, int k,
                                  Expr* var, Expr* imin) {
    if (sr_pdeg(base, var) != 1) { expr_free(c0); return NULL; }

    Expr* zero = sum_int(0);
    Expr* one  = sum_int(1);
    Expr* b0 = sum_subst(base, var, zero);          /* base(0) */
    Expr* b_at1 = sum_subst(base, var, one);        /* base(1) */
    expr_free(zero); expr_free(one);
    Expr* b1 = sum_sub(b_at1, b0);                  /* b1 = base(1)-base(0) */
    expr_free(b_at1);

    /* rho = -b0/b1 */
    Expr* rho = sum_eval("Times",
                    (Expr*[]){ sum_int(-1), b0,
                               expr_new_function(expr_new_symbol(SYM_Power),
                                   (Expr*[]){ expr_copy(b1), sum_int(-1) }, 2) }, 3);
    /* a = imin - rho */
    Expr* a = sum_sub(imin, rho);
    expr_free(rho);

    /* c = c0 * b1^(-k) */
    Expr* c = sum_eval("Times",
                  (Expr*[]){ c0,
                             expr_new_function(expr_new_symbol(SYM_Power),
                                 (Expr*[]){ b1, sum_int(-k) }, 2) }, 2);

    Expr* contrib;
    if (k == 1) {
        /* -c (PolyGamma[0, a] + EulerGamma) */
        Expr* pg = expr_new_function(expr_new_symbol(SYM_PolyGamma),
                       (Expr*[]){ sum_int(0), a }, 2);
        Expr* eg = expr_new_symbol(SYM_EulerGamma);
        Expr* bracket = expr_new_function(expr_new_symbol(SYM_Plus),
                            (Expr*[]){ pg, eg }, 2);
        contrib = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ sum_int(-1), c, bracket }, 3);
    } else {
        /* c Zeta[k, a] */
        Expr* zeta = expr_new_function(expr_new_symbol(SYM_Zeta),
                         (Expr*[]){ sum_int(k), a }, 2);
        contrib = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ c, zeta }, 2);
    }
    return contrib;
}

/*
 * Closed form of Sum_{i>=imin} N(i)/base(i) for a single irreducible quadratic
 * denominator with complex-conjugate roots (discriminant < 0) and multiplicity
 * k == 1.  Writing base = b2 ((i - alpha)^2 + beta^2) and N = A i + B, split the
 * numerator about i = alpha:
 *
 *   N/b2 = Atil (i - alpha) + S,   Atil = A/b2,  S = Atil alpha + B/b2,
 *
 * and sum the two pieces (n = i - alpha, t = imin - alpha):
 *
 *   odd part   Atil (i-alpha)/((i-alpha)^2+beta^2):
 *       -> -(Atil/2)(PolyGamma[0, t-beta I] + PolyGamma[0, t+beta I] + 2 gamma)
 *          (a conjugate digamma sum -- no elementary form, kept as PolyGamma).
 *   even part  S /((i-alpha)^2+beta^2):
 *       -> S ((pi beta Coth[pi beta] - 1)/(2 beta^2) - sum_{j=1}^{t-1} 1/(j^2+beta^2))
 *          (the symmetric series collapses to Coth).
 *
 * Returns NULL -- so the caller falls back to the radical/complex extension
 * route -- when the discriminant is not (numerically) negative (real radical
 * roots), or t = imin - alpha is not a positive integer.
 */
static Expr* sr_quadratic_contribution(Expr* term, Expr* base, Expr* var, Expr* imin) {
    Expr* b2 = sum_eval("Coefficient", (Expr*[]){ expr_copy(base), expr_copy(var), sum_int(2) }, 3);
    Expr* b1 = sum_eval("Coefficient", (Expr*[]){ expr_copy(base), expr_copy(var), sum_int(1) }, 3);
    Expr* b0 = sum_eval("Coefficient", (Expr*[]){ expr_copy(base), expr_copy(var), sum_int(0) }, 3);

    /* disc = b1^2 - 4 b2 b0; require numerically negative (complex roots). */
    Expr* b1sq   = sum_eval("Times", (Expr*[]){ expr_copy(b1), expr_copy(b1) }, 2);
    Expr* fourbb = sum_eval("Times", (Expr*[]){ sum_int(-4), expr_copy(b2), expr_copy(b0) }, 3);
    Expr* disc   = sum_eval("Together", (Expr*[]){
                       sum_eval("Plus", (Expr*[]){ b1sq, fourbb }, 2) }, 1);
    Expr* discN  = sum_eval("N", (Expr*[]){ expr_copy(disc) }, 1);
    bool neg = (discN->type == EXPR_REAL && discN->data.real < 0.0)
            || (discN->type == EXPR_INTEGER && discN->data.integer < 0);
    expr_free(discN);
    if (!neg) { expr_free(b2); expr_free(b1); expr_free(b0); expr_free(disc); return NULL; }

    /* alpha = -b1/(2 b2). */
    Expr* twob2 = sum_eval("Times", (Expr*[]){ sum_int(2), expr_copy(b2) }, 2);
    Expr* alpha = sum_eval("Together", (Expr*[]){
                      sum_eval("Times", (Expr*[]){ sum_int(-1), expr_copy(b1),
                          expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ twob2, sum_int(-1) }, 2) }, 3) }, 1);

    /* t = imin - alpha must be a positive integer for the Coth collapse. */
    Expr* t = sum_sub(imin, alpha);
    if (t->type != EXPR_INTEGER || t->data.integer < 1) {
        expr_free(b2); expr_free(b1); expr_free(b0); expr_free(disc);
        expr_free(alpha); expr_free(t);
        return NULL;
    }
    int64_t tval = t->data.integer;
    expr_free(t);

    /* beta2 = -disc/(4 b2^2); beta = Sqrt[beta2]. */
    Expr* b2sq     = sum_eval("Times", (Expr*[]){ expr_copy(b2), expr_copy(b2) }, 2);
    Expr* fourb2sq = sum_eval("Times", (Expr*[]){ sum_int(4), b2sq }, 2);
    Expr* beta2 = sum_eval("Together", (Expr*[]){
                      sum_eval("Times", (Expr*[]){ sum_int(-1), disc,
                          expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ fourb2sq, sum_int(-1) }, 2) }, 3) }, 1);
    Expr* beta = sum_eval("Sqrt", (Expr*[]){ expr_copy(beta2) }, 1);

    /* numerator N = Cancel[term*base]; A = [i^1], B = [i^0]; Atil=A/b2, S=Atil*alpha+B/b2. */
    Expr* N = sum_eval("Cancel", (Expr*[]){
                  sum_eval("Times", (Expr*[]){ expr_copy(term), expr_copy(base) }, 2) }, 1);
    Expr* A = sum_eval("Coefficient", (Expr*[]){ expr_copy(N), expr_copy(var), sum_int(1) }, 3);
    Expr* B = sum_eval("Coefficient", (Expr*[]){ N, expr_copy(var), sum_int(0) }, 3);
    Expr* invb2  = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(b2), sum_int(-1) }, 2);
    Expr* invb2b = expr_copy(invb2);
    Expr* Atil = sum_eval("Times", (Expr*[]){ A, invb2 }, 2);
    Expr* Btil = sum_eval("Times", (Expr*[]){ B, invb2b }, 2);
    expr_free(b2); expr_free(b1); expr_free(b0);
    Expr* S = sum_eval("Plus", (Expr*[]){
                  sum_eval("Times", (Expr*[]){ expr_copy(Atil), expr_copy(alpha) }, 2), Btil }, 2);
    expr_free(alpha);

    /* a = t - beta I, abar = t + beta I. */
    Expr* bI = sum_eval("Times", (Expr*[]){ expr_copy(beta), expr_new_symbol(SYM_I) }, 2);
    Expr* a    = sum_eval("Plus", (Expr*[]){ sum_int(tval),
                     sum_eval("Times", (Expr*[]){ sum_int(-1), expr_copy(bI) }, 2) }, 2);
    Expr* abar = sum_eval("Plus", (Expr*[]){ sum_int(tval), bI }, 2);

    /* odd part: -(Atil/2)(PolyGamma[0,a] + PolyGamma[0,abar] + 2 EulerGamma). */
    Expr* pg_a  = expr_new_function(expr_new_symbol(SYM_PolyGamma), (Expr*[]){ sum_int(0), a }, 2);
    Expr* pg_ab = expr_new_function(expr_new_symbol(SYM_PolyGamma), (Expr*[]){ sum_int(0), abar }, 2);
    Expr* twog  = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ sum_int(2), expr_new_symbol(SYM_EulerGamma) }, 2);
    Expr* brk = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ pg_a, pg_ab, twog }, 3);
    Expr* halfA = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                      sum_int(-1), Atil,
                      expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ sum_int(2), sum_int(-1) }, 2) }, 3);
    Expr* odd = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ halfA, brk }, 2);

    /* even part: S ((Pi beta Coth[Pi beta]-1)/(2 beta2) - sum_{j=1}^{t-1} 1/(j^2+beta2)). */
    Expr* pibeta = sum_eval("Times", (Expr*[]){ expr_new_symbol("Pi"), beta }, 2);
    Expr* coth   = expr_new_function(expr_new_symbol("Coth"), (Expr*[]){ expr_copy(pibeta) }, 1);
    Expr* num_main = sum_eval("Plus", (Expr*[]){
                         sum_eval("Times", (Expr*[]){ pibeta, coth }, 2), sum_int(-1) }, 2);
    Expr* twobeta2 = sum_eval("Times", (Expr*[]){ sum_int(2), expr_copy(beta2) }, 2);
    Expr* coth_main = sum_eval("Times", (Expr*[]){ num_main,
                          expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ twobeta2, sum_int(-1) }, 2) }, 2);
    Expr* corr = sum_int(0);
    for (int64_t j = 1; j < tval; j++) {
        Expr* denj  = sum_eval("Plus", (Expr*[]){ sum_int(j * j), expr_copy(beta2) }, 2);
        Expr* termj = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ denj, sum_int(-1) }, 2);
        corr = sum_eval("Plus", (Expr*[]){ corr, termj }, 2);
    }
    expr_free(beta2);
    Expr* cothpart = sum_eval("Plus", (Expr*[]){ coth_main,
                         sum_eval("Times", (Expr*[]){ sum_int(-1), corr }, 2) }, 2);
    Expr* even = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ S, cothpart }, 2);

    Expr* contrib = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ odd, even }, 2);
    Expr* r = evaluate(contrib);
    expr_free(contrib);
    return r;
}

/*
 * Decompose f via Apart -- optionally over an Extension field -- and sum each
 * linear-pole term with the master identity.  When `ext` is non-NULL it is an
 * owned List of field generators ({I, Sqrt[5], ...}) and is consumed into the
 * Apart[..., Extension -> ext] call.  Returns the owned closed form, or NULL if
 * any additive term fails to classify as a linear pole (e.g. an unfactored
 * irreducible quadratic when no/insufficient extension was supplied).
 */
static Expr* sr_decompose(Expr* f, Expr* var, Expr* imin, Expr* ext,
                          bool* needs_ext) {
    if (needs_ext) *needs_ext = false;
    Expr* apart;
    if (ext) {
        Expr* opt = expr_new_function(expr_new_symbol(SYM_Rule),
                        (Expr*[]){ expr_new_symbol(SYM_Extension), ext }, 2);
        apart = sum_eval("Apart",
                    (Expr*[]){ expr_copy(f), expr_copy(var), opt }, 3);
    } else {
        apart = sum_eval("Apart", (Expr*[]){ expr_copy(f), expr_copy(var) }, 2);
    }

    /* Normalize to a list of additive terms. */
    Expr** terms;
    size_t nterms;
    bool is_plus = apart->type == EXPR_FUNCTION
                && apart->data.function.head->type == EXPR_SYMBOL
                && apart->data.function.head->data.symbol.name == SYM_Plus;
    if (is_plus) {
        nterms = apart->data.function.arg_count;
        terms = apart->data.function.args;
    } else {
        nterms = 1;
        terms = &apart;
    }

    Expr** contribs = malloc(sizeof(Expr*) * (nterms ? nterms : 1));
    size_t nc = 0;
    bool ok = true;
    for (size_t i = 0; i < nterms; i++) {
        Expr *base = NULL; int k = 0;
        if (sr_find_pole(terms[i], var, &base, &k) != 1) { ok = false; break; }
        int d = sr_pdeg(base, var);
        Expr* contrib = NULL;
        if (d == 1) {
            /* Linear pole: master identity (Hurwitz Zeta / digamma). */
            Expr* c0 = sr_residue(terms[i], base, k, var);
            contrib = sr_term_contribution(c0, base, k, var, imin);
        } else if (d == 2 && !ext) {
            /* Irreducible quadratic over Q: complex-conjugate roots collapse to
             * a Coth / conjugate-digamma form (only k == 1).  Real radical roots
             * (disc >= 0) or higher multiplicity return NULL -> extension path. */
            if (k == 1)
                contrib = sr_quadratic_contribution(terms[i], base, var, imin);
            if (!contrib) { ok = false; if (needs_ext) *needs_ext = true; break; }
        }
        if (!contrib) {
            ok = false;
            /* Over Q a non-linear leftover means we still need the extension. */
            if (needs_ext && !ext) *needs_ext = true;
            break;
        }
        contribs[nc++] = contrib;
    }

    if (!ok) {
        for (size_t i = 0; i < nc; i++) expr_free(contribs[i]);
        free(contribs);
        expr_free(apart);
        return NULL;
    }

    Expr* total = expr_new_function(expr_new_symbol(SYM_Plus), contribs, nc);
    free(contribs);
    Expr* evaled = evaluate(total);
    expr_free(total);
    expr_free(apart);

    /* Simplify for a clean WL-matching form: it rationalises the complex /
     * radical residue coefficients (cancelling the spurious overall I that a
     * bare Together leaves) and groups the shared EulerGamma in the conjugate
     * case.  Taken only when strictly simpler: on a tie the evaluated sum
     * keeps its form, as the reference's Sum does -- Sum[1/i^2, {i, 2,
     * Infinity}] is -1 + Pi^2/6, where Simplify takes the content out. */
    Expr* plain = expr_copy(evaled);
    Expr* result = sum_eval("Simplify", (Expr*[]){ evaled }, 1);
    if (result && simp_default_complexity(result) >= simp_default_complexity(plain)) {
        expr_free(result);
        return plain;
    }
    expr_free(plain);
    return result;
}

/*
 * Recursively gather the algebraic generators implicit in a Solve root `e`:
 * sets *need_i when a nonzero imaginary part appears, and appends each distinct
 * surd Power[c, p/q] (non-integer rational exponent) to *surds.  Sets
 * *unsupported on a Root[] object, whose generator cannot be extracted as a
 * radical (defer to the direct-residue path / leave held).
 */
static void sr_walk_gens(Expr* e, bool* need_i, Expr*** surds, size_t* ns,
                         size_t* cap, bool* unsupported) {
    if (!e || *unsupported) return;
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol.name == SYM_I) *need_i = true;
        return;
    }
    if (e->type != EXPR_FUNCTION) return;

    Expr* head = e->data.function.head;
    size_t argc = e->data.function.arg_count;
    if (head->type == EXPR_SYMBOL) {
        if (head->data.symbol.name == SYM_Complex && argc == 2) {
            Expr* im = e->data.function.args[1];
            bool zero = (im->type == EXPR_INTEGER && im->data.integer == 0)
                     || (im->type == EXPR_REAL && im->data.real == 0.0);
            if (!zero) *need_i = true;
            return;
        }
        if (head->data.symbol.name == SYM_Root) { *unsupported = true; return; }
        if (head->data.symbol.name == SYM_Power && argc == 2) {
            Expr* ex = e->data.function.args[1];
            if (ex->type == EXPR_FUNCTION
                && ex->data.function.head->type == EXPR_SYMBOL
                && ex->data.function.head->data.symbol.name == SYM_Rational) {
                for (size_t i = 0; i < *ns; i++)
                    if (expr_eq((*surds)[i], e)) return;   /* dedup */
                if (*ns == *cap) {
                    *cap *= 2;
                    *surds = realloc(*surds, sizeof(Expr*) * (*cap));
                }
                (*surds)[(*ns)++] = expr_copy(e);
                return;
            }
        }
    }
    sr_walk_gens(head, need_i, surds, ns, cap, unsupported);
    for (size_t i = 0; i < argc; i++)
        sr_walk_gens(e->data.function.args[i], need_i, surds, ns, cap, unsupported);
}

/*
 * Build the Extension generator list for f's denominator: solve den == 0,
 * walk the roots, and assemble {I?, surds...}.  Returns an owned List, or NULL
 * if there are no algebraic generators (pure rational, already handled) or the
 * roots involve an inextricable Root[] object.
 */
static Expr* sr_generators(Expr* den, Expr* var) {
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                   (Expr*[]){ expr_copy(den), sum_int(0) }, 2);
    Expr* sol = sum_eval("Solve", (Expr*[]){ eq, expr_copy(var) }, 2);
    /* roots = var /. sol */
    Expr* roots = sum_eval("ReplaceAll",
                      (Expr*[]){ expr_copy(var), sol }, 2);
    if (!roots || roots->type != EXPR_FUNCTION
        || roots->data.function.head->type != EXPR_SYMBOL
        || roots->data.function.head->data.symbol.name != SYM_List) {
        if (roots) expr_free(roots);
        return NULL;
    }

    bool need_i = false, unsupported = false;
    size_t cap = 4, ns = 0;
    Expr** surds = malloc(sizeof(Expr*) * cap);
    sr_walk_gens(roots, &need_i, &surds, &ns, &cap, &unsupported);
    expr_free(roots);

    if (unsupported || (!need_i && ns == 0)) {
        for (size_t i = 0; i < ns; i++) expr_free(surds[i]);
        free(surds);
        return NULL;
    }

    size_t total = ns + (need_i ? 1 : 0);
    Expr** gens = malloc(sizeof(Expr*) * total);
    size_t g = 0;
    if (need_i) gens[g++] = expr_new_symbol(SYM_I);
    for (size_t i = 0; i < ns; i++) gens[g++] = surds[i];
    free(surds);
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), gens, total);
    free(gens);   /* expr_new_function copies the array; free our copy */
    return list;
}

/* ------------------------------------------------------------------ */
/*  Finite ranges: Sum[p(i)/q(i), {i, imin, imax}] with integer bounds  */
/* ------------------------------------------------------------------ */

/* PolyGamma[m, a] (owned args). */
static Expr* sr_polygamma(int m, Expr* a) {
    return expr_new_function(expr_new_symbol(SYM_PolyGamma),
                             (Expr*[]){ sum_int(m), a }, 2);
}

/* Head test without evaluation. */
static bool sr_head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

/* The finite sum  sum_{i=imin}^{imax} 1/(i - rho)^k  as a difference of two
 * polygamma values at A = imin - rho and B = imax + 1 - rho:
 *
 *   k == 1:  PolyGamma[0, B] - PolyGamma[0, A]
 *   k >= 2:  (-1)^k/(k-1)! (PolyGamma[k-1, A] - PolyGamma[k-1, B])
 *
 * (Zeta[k, A] - Zeta[k, B] with Zeta[k, x] = (-1)^k PolyGamma[k-1, x]/(k-1)!;
 * both are the telescoped digamma recurrence psi(x+1) - psi(x) = 1/x.)  These
 * are the reference's forms; at a small positive-integer A or B the PolyGamma
 * evaluates to its exact harmonic value, elsewhere it stays symbolic.  Returns
 * an owned, unevaluated expression.  A, B are consumed. */
static Expr* sr_finite_pole_sum(int k, Expr* A, Expr* B) {
    if (k == 1) {
        Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ sum_int(-1), sr_polygamma(0, A) }, 2);
        return expr_new_function(expr_new_symbol(SYM_Plus),
                   (Expr*[]){ sr_polygamma(0, B), neg }, 2);
    }
    mpz_t fact;
    mpz_init(fact);
    mpz_fac_ui(fact, (unsigned long)(k - 1));
    if (k % 2) mpz_neg(fact, fact);              /* (-1)^k (k-1)! */
    mpz_t one;
    mpz_init_set_ui(one, 1);
    Expr* coef = make_rational_mpz(one, fact);   /* (-1)^k/(k-1)! */
    mpz_clears(fact, one, (mpz_ptr)0);
    Expr* negB = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ sum_int(-1), sr_polygamma(k - 1, B) }, 2);
    Expr* diff = expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ sr_polygamma(k - 1, A), negB }, 2);
    return expr_new_function(expr_new_symbol(SYM_Times),
               (Expr*[]){ coef, diff }, 2);
}

/* Contribution of one linear pole term c0/base^k (base = b1 i + b0) to the
 * finite sum: c0 b1^(-k) times sr_finite_pole_sum at rho = -b0/b1.  Sets
 * *in_range when rho is a real rational in [imin, imax] -- the summand itself
 * is singular at an index of the range, so there is no finite value to give.
 * Returns NULL when base is not linear.  c0 is consumed, base borrowed. */
static Expr* sr_finite_linear_term(Expr* c0, Expr* base, int k, Expr* var,
                                   Expr* imin, Expr* imax, bool* in_range) {
    if (sr_pdeg(base, var) != 1) { expr_free(c0); return NULL; }
    Expr* zero = sum_int(0);
    Expr* one  = sum_int(1);
    Expr* b0 = sum_subst(base, var, zero);
    Expr* b_at1 = sum_subst(base, var, one);
    expr_free(zero); expr_free(one);
    Expr* b1 = sum_sub(b_at1, b0);
    expr_free(b_at1);

    Expr* rho = sum_eval("Times",
                    (Expr*[]){ sum_int(-1), b0,
                               expr_new_function(expr_new_symbol(SYM_Power),
                                   (Expr*[]){ expr_copy(b1), sum_int(-1) }, 2) }, 3);
    /* A pole at a non-integer rational is never hit by an integer index. */
    if (expr_is_integer_like(rho)) {
        Expr* c = sum_eval("LessEqual", (Expr*[]){ expr_copy(imin), expr_copy(rho),
                                                    expr_copy(imax) }, 3);
        if (c->type == EXPR_SYMBOL && c->data.symbol.name == SYM_True) *in_range = true;
        expr_free(c);
    }
    Expr* A = sum_sub(imin, rho);
    Expr* imax1 = sum_eval("Plus", (Expr*[]){ expr_copy(imax), sum_int(1) }, 2);
    Expr* B = sum_sub(imax1, rho);
    expr_free(imax1);
    expr_free(rho);

    Expr* c = sum_eval("Times",
                  (Expr*[]){ c0,
                             expr_new_function(expr_new_symbol(SYM_Power),
                                 (Expr*[]){ b1, sum_int(-k) }, 2) }, 2);
    return expr_new_function(expr_new_symbol(SYM_Times),
               (Expr*[]){ c, sr_finite_pole_sum(k, A, B) }, 2);
}

/* Contribution of a simple pole term A(i)/q(i), q irreducible over Q of degree
 * >= 2 with roots the field-generator search cannot express by radicals
 * (Root objects): sum over the roots r of q of the residue R(r) = A(r)/q'(r)
 * times the digamma difference, as a RootSum (the reference's form):
 *
 *   RootSum[q(#)&, R(#) PolyGamma[0, imax + 1 - #]&]
 *     - RootSum[q(#)&, R(#) PolyGamma[0, imin - #]&]
 *
 * R is reduced to a polynomial of degree < deg q:  R = A s mod q, where
 * s q' + t q = g (PolynomialExtendedGCD, g a nonzero constant since q is
 * squarefree).  Returns NULL when q is not squarefree-coprime to q' or the
 * pieces fail to come back polynomial.  term and base are borrowed. */
static Expr* sr_finite_rootsum_term(Expr* term, Expr* base, Expr* var,
                                    Expr* imin, Expr* imax) {
    Expr* numer = sr_residue(term, base, 1, var);           /* A(i) */
    if (!sr_poly_in(numer, var)) { expr_free(numer); return NULL; }
    Expr* dq = sum_eval("D", (Expr*[]){ expr_copy(base), expr_copy(var) }, 2);
    Expr* eg = sum_eval("PolynomialExtendedGCD",
                   (Expr*[]){ dq, expr_copy(base), expr_copy(var) }, 3);
    /* eg = {g, {s, t}} */
    if (!sr_head_is(eg, SYM_List) || eg->data.function.arg_count != 2
        || !sr_head_is(eg->data.function.args[1], SYM_List)
        || eg->data.function.args[1]->data.function.arg_count != 2
        || !sum_free_of(eg->data.function.args[0], var)) {
        expr_free(eg); expr_free(numer); return NULL;
    }
    Expr* g = expr_copy(eg->data.function.args[0]);
    Expr* s = expr_copy(eg->data.function.args[1]->data.function.args[0]);
    expr_free(eg);
    Expr* as = sum_eval("Expand", (Expr*[]){
                   sum_eval("Times", (Expr*[]){ numer, s,
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ g, sum_int(-1) }, 2) }, 3) }, 1);
    Expr* R = sum_eval("PolynomialRemainder",
                  (Expr*[]){ as, expr_copy(base), expr_copy(var) }, 3);
    if (!sr_poly_in(R, var)) { expr_free(R); return NULL; }

    Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot),
                                   (Expr*[]){ sum_int(1) }, 1);
    Expr* q_s = sum_subst(base, var, slot);
    Expr* R_s = sum_subst(R, var, slot);
    expr_free(R);

    Expr* imax1 = sum_eval("Plus", (Expr*[]){ expr_copy(imax), sum_int(1) }, 2);
    Expr* ends[2] = { imax1, expr_copy(imin) };
    Expr* rs[2];
    for (int e = 0; e < 2; e++) {
        Expr* arg = expr_new_function(expr_new_symbol(SYM_Plus),
                        (Expr*[]){ ends[e],
                                   expr_new_function(expr_new_symbol(SYM_Times),
                                       (Expr*[]){ sum_int(-1), expr_copy(slot) }, 2) }, 2);
        Expr* body = expr_new_function(expr_new_symbol(SYM_Times),
                         (Expr*[]){ expr_copy(R_s), sr_polygamma(0, arg) }, 2);
        Expr* fq = expr_new_function(expr_new_symbol(SYM_Function),
                       (Expr*[]){ expr_copy(q_s) }, 1);
        Expr* fb = expr_new_function(expr_new_symbol(SYM_Function),
                       (Expr*[]){ body }, 1);
        rs[e] = expr_new_function(expr_new_symbol(SYM_RootSum),
                    (Expr*[]){ fq, fb }, 2);
    }
    expr_free(slot); expr_free(q_s); expr_free(R_s);
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ sum_int(-1), rs[1] }, 2);
    return expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ rs[0], neg }, 2);
}

/*
 * Decompose f with Apart (over Q, or over `ext` when non-NULL -- consumed) and
 * sum every additive piece over the finite range: the polynomial part through
 * Sum`Polynomial, each linear pole through sr_finite_linear_term, and -- only
 * when allow_rootsum -- each simple pole at an irreducible non-linear factor
 * through sr_finite_rootsum_term.  On a non-linear pole that is not handled,
 * sets *needs_ext (when non-NULL) and returns NULL.  Sets *in_range when a pole
 * of the summand lies on an index of the range.
 */
static Expr* sr_finite_decompose(Expr* f, Expr* var, Expr* imin, Expr* imax,
                                 Expr* ext, bool allow_rootsum,
                                 bool* needs_ext, bool* in_range) {
    Expr* apart;
    if (ext) {
        Expr* opt = expr_new_function(expr_new_symbol(SYM_Rule),
                        (Expr*[]){ expr_new_symbol(SYM_Extension), ext }, 2);
        apart = sum_eval("Apart", (Expr*[]){ expr_copy(f), expr_copy(var), opt }, 3);
    } else {
        apart = sum_eval("Apart", (Expr*[]){ expr_copy(f), expr_copy(var) }, 2);
    }
    Expr** terms;
    size_t nterms;
    if (sr_head_is(apart, SYM_Plus)) {
        nterms = apart->data.function.arg_count;
        terms = apart->data.function.args;
    } else {
        nterms = 1;
        terms = &apart;
    }

    Expr** contribs = malloc(sizeof(Expr*) * (nterms + 1));
    Expr** polys = malloc(sizeof(Expr*) * (nterms + 1));
    size_t nc = 0, np = 0;
    bool ok = true;
    for (size_t i = 0; i < nterms && ok; i++) {
        Expr *base = NULL; int k = 0;
        int npoles = sr_find_pole(terms[i], var, &base, &k);
        if (npoles == 0) {
            if (!sr_poly_in(terms[i], var)) { ok = false; break; }
            polys[np++] = expr_copy(terms[i]);
            continue;
        }
        if (npoles != 1) { ok = false; break; }
        Expr* contrib = NULL;
        int d = sr_pdeg(base, var);
        if (d == 1) {
            Expr* c0 = sr_residue(terms[i], base, k, var);
            contrib = sr_finite_linear_term(c0, base, k, var, imin, imax, in_range);
        } else if (d >= 2 && k == 1 && allow_rootsum) {
            contrib = sr_finite_rootsum_term(terms[i], base, var, imin, imax);
        } else if (needs_ext && !ext) {
            *needs_ext = true;
        }
        if (!contrib) { ok = false; break; }
        contribs[nc++] = contrib;
    }
    expr_free(apart);

    if (ok && np > 0) {
        /* The polynomial part sums in closed form (Faulhaber). */
        Expr* p = expr_new_function(expr_new_symbol(SYM_Plus), polys, np);
        np = 0;
        Expr* ps = sum_eval("Sum`Polynomial",
                       (Expr*[]){ p, expr_copy(var), expr_copy(imin), expr_copy(imax) }, 4);
        if (!ps || sr_head_is(ps, "Sum`Polynomial")) {
            if (ps) expr_free(ps);
            ok = false;
        } else {
            contribs[nc++] = ps;
        }
    }
    for (size_t i = 0; i < np; i++) expr_free(polys[i]);
    free(polys);
    if (!ok || *in_range) {
        for (size_t i = 0; i < nc; i++) expr_free(contribs[i]);
        free(contribs);
        return NULL;
    }
    Expr* total = expr_new_function(expr_new_symbol(SYM_Plus), contribs, nc);
    free(contribs);
    Expr* r = evaluate(total);
    expr_free(total);
    return r;
}

/*
 * Finite Sum`Rational: imin and imax exact integers, f a rational function of
 * var.  Unlike the infinite case there is no convergence condition: every
 * rational summand without a pole on the range has a closed form in
 * polygamma values, and its cost does not depend on the width of the range --
 * the term-by-term sum of Sum[1/(k + 3/7), {k, 12, 3*10^6}] is a rational of
 * millions of digits, built by millions of big-number additions.
 *
 * Over Q first; then, for a denominator whose roots are radicals or complex
 * (k^2 + 1, k^2 - 2), Apart over the generated extension so every pole is
 * linear (PolyGamma at algebraic arguments, the reference's form); finally,
 * for roots only expressible as Root objects (k^3 + k + 1), a RootSum per
 * irreducible simple factor.  A pole on an index of the range, a repeated
 * non-linear irreducible factor, or a non-rational summand returns NULL.
 */
static Expr* sr_finite(Expr* f, Expr* var, Expr* imin, Expr* imax) {
    if (!expr_is_integer_like(imin) || !expr_is_integer_like(imax)) return NULL;
    /* An inexact summand (1.0/k, N[1/k]) is summed numerically by the caller's
     * compiled loop: Apart would rationalise it and return an exact form. */
    if (internal_contains_inexact(f)) return NULL;
    Expr* tog = sum_eval("Together", (Expr*[]){ expr_copy(f) }, 1);
    Expr* num = sum_eval("Numerator", (Expr*[]){ expr_copy(tog) }, 1);
    Expr* den = sum_eval("Denominator", (Expr*[]){ tog }, 1);
    bool rational = sr_poly_in(num, var) && sr_poly_in(den, var)
                 && sr_pdeg(den, var) >= 1;
    expr_free(num);
    if (!rational) { expr_free(den); return NULL; }

    bool needs_ext = false, in_range = false;
    Expr* result = sr_finite_decompose(f, var, imin, imax, NULL, false,
                                       &needs_ext, &in_range);
    if (!result && needs_ext && !in_range) {
        Expr* ext = sr_generators(den, var);
        if (ext)
            result = sr_finite_decompose(f, var, imin, imax, ext, false,
                                         NULL, &in_range);
        else
            result = sr_finite_decompose(f, var, imin, imax, NULL, true,
                                         NULL, &in_range);
    }
    expr_free(den);
    return result;
}

Expr* builtin_sum_rational(Expr* res);

Expr* builtin_sum_rational(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite) return NULL;

    /* A finite range with exact integer bounds: the polygamma closed form. */
    if (!is_infinity_sym(imax)) return sr_finite(f, var, imin, imax);

    /* Infinite sums with a concrete integer lower bound. */
    if (imin->type != EXPR_INTEGER) return NULL;

    /* --- Convergence gate (before the expensive Apart). ---
     * Together, then split into Numerator/Denominator; both must be polynomial
     * in var, and deg q >= deg p + 2 (else the series diverges -> stay held). */
    Expr* tog = sum_eval("Together", (Expr*[]){ expr_copy(f) }, 1);
    Expr* num = sum_eval("Numerator", (Expr*[]){ expr_copy(tog) }, 1);
    Expr* den = sum_eval("Denominator", (Expr*[]){ tog }, 1);
    if (!sr_poly_in(num, var) || !sr_poly_in(den, var)) {
        expr_free(num); expr_free(den);
        return NULL;
    }
    int degp = sr_pdeg(num, var);
    int degq = sr_pdeg(den, var);
    expr_free(num);
    if (degq < 0 || degp < 0 || degq < degp + 2) { expr_free(den); return NULL; }

    /* Phases A & C over Q: linear poles -> Hurwitz Zeta / digamma; irreducible
     * quadratics with complex-conjugate roots -> Coth / conjugate-digamma form. */
    bool needs_ext = false;
    Expr* result = sr_decompose(f, var, imin, NULL, &needs_ext);

    /* Phase B: a factor with real radical roots (or higher degree) did not split
     * over Q -- discover the field generators and retry Apart over the
     * extension, which fully splits everything into linear poles. */
    if (!result && needs_ext) {
        Expr* ext = sr_generators(den, var);
        if (ext) result = sr_decompose(f, var, imin, ext, NULL);
    }
    expr_free(den);
    return result;
}

void sum_rational_init(void) {
    symtab_add_builtin("Sum`Rational", builtin_sum_rational);
    symtab_get_def("Sum`Rational")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Sum`Rational",
        "Sum`Rational[f, i, imin, Infinity] gives the closed form of an infinite "
        "sum of a rational function f of i (deg of denominator >= deg of "
        "numerator + 2); Sum`Rational[f, i, imin, imax] with integer bounds gives "
        "the finite sum. Decomposes into partial fractions and sums each term via "
        "Hurwitz Zeta / PolyGamma (RootSum for poles at Root objects). Returns "
        "unevaluated for symbolic finite bounds, indefinite, non-rational, or "
        "divergent inputs, or a pole on the range.");
}
