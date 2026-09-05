/*
 * dsolve_undetcoeff.c — DSolve`UndeterminedCoefficients.
 *
 * A tidy particular solution of a constant-coefficient linear ODE
 *     a_n y^(n) + ... + a_0 y == g(x)
 * whose forcing g is a "UC function" — a sum of terms x^m Exp[a x] {1|Cos[b x]|Sin[b x]}.
 * The homogeneous part is the usual characteristic-root fundamental set
 * (dsolve_homog_basis); the particular is found by superposition over the additive
 * terms of g, each with the ansatz
 *     y_p = x^s Exp[a x] (Cos[b x] Sum A_k x^k + Sin[b x] Sum B_k x^k)
 * (drop the trig factor when b == 0), where s (the resonance shift) is found by
 * trying s = 0, 1, 2, ... until the linear system for the A_k, B_k is solvable —
 * so no separate root-multiplicity bookkeeping is needed.  Mirrors SymPy's
 * `nth_linear_constant_coeff_undetermined_coefficients`.
 *
 * Runs immediately BEFORE LinearConstantCoefficients in the cascade: UC forcing
 * gets this tidy particular, everything else falls through to that method's
 * variation-of-parameters (which handles any forcing).  Constant coefficients
 * only in this first cut; g == 0 (homogeneous) is left to constcoeff.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "../parse.h"
#include <stdlib.h>

/* First match of `ruleStr` (with placeholder XV for the independent variable, uu_
 * for the pattern var) anywhere in T, else the integer `deflt`. */
static Expr* uc_first_case(const Expr* T, const char* ruleStr, const char* xvar, int deflt) {
    Expr* rule = parse_expression(ruleStr);
    if (!rule) return expr_new_integer(deflt);
    rule = ds_subst(rule, expr_new_symbol("XV"), expr_new_symbol(xvar));
    Expr* lvl = expr_new_function(expr_new_symbol(SYM_List),
                    (Expr*[]){ expr_new_integer(0), expr_new_symbol("Infinity") }, 2);
    Expr* cs = eval_and_free(expr_new_function(expr_new_symbol("Cases"),
                   (Expr*[]){ expr_copy((Expr*)T), rule, lvl }, 3));
    Expr* out = expr_new_integer(deflt);
    if (head_is(cs, SYM_List) && cs->data.function.arg_count >= 1)
        out = expr_copy(cs->data.function.args[0]);
    expr_free(cs);
    return out;
}

/* Sum_{k=0}^{d} sym[k] x^k, and append the sym[k] funcapps to *vars (grown). */
static Expr* uc_poly(const char* sym, int d, const char* xvar, Expr*** vars, size_t* nv) {
    Expr** terms = malloc((size_t)(d + 1) * sizeof(Expr*));
    for (int k = 0; k <= d; k++) {
        Expr* v = expr_new_function(expr_new_symbol(sym), (Expr*[]){ expr_new_integer(k) }, 1);
        Expr* xk = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_new_symbol(xvar), expr_new_integer(k) }, 2);
        terms[k] = ds_call2(SYM_Times, expr_copy(v), xk);
        *vars = realloc(*vars, (*nv + 1) * sizeof(Expr*));
        (*vars)[(*nv)++] = v;
    }
    Expr* s = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)(d + 1));
    free(terms);
    return s;
}

/* Coefficient[e, v] with v a full expression (e.g. Cos[b x]); e,v borrowed. */
static Expr* coeff_of(const Expr* e, const Expr* v) {
    return eval_and_free(ds_call2("Coefficient", expr_copy((Expr*)e), expr_copy((Expr*)v)));
}

/* Particular solution for a single UC term T (in xvar), given the constant
 * coefficients a[0..n].  Returns owned y_p or NULL if T is not UC / not solvable. */
static Expr* uc_particular(const Expr* T, Expr** a, int n, const char* xvar) {
    /* extract the exponential rate alpha and trig frequency beta */
    Expr* alpha = uc_first_case(T, "Power[E, uu_] :> Coefficient[uu, XV]", xvar, 0);
    Expr* beta  = uc_first_case(T, "(Cos[uu_] | Sin[uu_]) :> Coefficient[uu, XV]", xvar, 0);
    bool has_trig = !ds_is_zero(beta);

    /* polynomial degree d: strip the Exp / Cos / Sin factors */
    Expr* stripRules = parse_expression("{Power[E, _] -> 1, Cos[_] -> 1, Sin[_] -> 1}");
    Expr* polyPart = eval_and_free(ds_call2("ReplaceAll", expr_copy((Expr*)T), stripRules));
    Expr* degE = eval_and_free(ds_call2("Exponent", polyPart, expr_new_symbol(xvar)));
    long d = (degE->type == EXPR_INTEGER) ? degE->data.integer : -1;
    expr_free(degE);
    if (d < 0) { expr_free(alpha); expr_free(beta); return NULL; }

    const char* uca = intern_symbol("DSolve`uca");
    const char* ucb = intern_symbol("DSolve`ucb");
    Expr* expfac = eval_and_free(ds_call1("Exp",
                       ds_call2(SYM_Times, expr_copy(alpha), expr_new_symbol(xvar))));
    Expr* cosfac = has_trig ? ds_call1("Cos", ds_call2(SYM_Times, expr_copy(beta), expr_new_symbol(xvar))) : NULL;
    Expr* sinfac = has_trig ? ds_call1("Sin", ds_call2(SYM_Times, expr_copy(beta), expr_new_symbol(xvar))) : NULL;

    Expr* yp = NULL;
    for (int s = 0; s <= n && !yp; s++) {
        Expr** vars = NULL; size_t nv = 0;
        Expr* pa = uc_poly(uca, (int)d, xvar, &vars, &nv);
        Expr* trig_combo;
        if (has_trig) {
            Expr* pb = uc_poly(ucb, (int)d, xvar, &vars, &nv);
            trig_combo = ds_call2(SYM_Plus,
                             ds_call2(SYM_Times, expr_copy(cosfac), pa),
                             ds_call2(SYM_Times, expr_copy(sinfac), pb));
        } else {
            trig_combo = pa;
        }
        Expr* xs = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_new_symbol(xvar), expr_new_integer(s) }, 2);
        Expr* trial = eval_and_free(ds_call2(SYM_Times, xs,
                          ds_call2(SYM_Times, expr_copy(expfac), trig_combo)));

        /* L[trial] = Sum_j a_j D[trial, {x, j}] */
        Expr* Lt = expr_new_integer(0);
        for (int j = 0; j <= n; j++) {
            Expr* dj = expr_copy(trial);
            for (int i = 0; i < j; i++) dj = ds_d(dj, expr_new_symbol(xvar));
            Lt = eval_and_free(ds_call2(SYM_Plus, Lt, ds_call2(SYM_Times, expr_copy(a[j]), dj)));
        }
        /* E = Expand[(L[trial] - T) / Exp[alpha x]].  Expand (not Simplify): Simplify
         * would rewrite Cos[b x]/Sin[b x] (e.g. Cos[2x] -> 1-2 Sin[x]^2), which
         * breaks the Coefficient[.,Cos[b x]] extraction below; the Exp[a x]/Exp[-a x]
         * cancellation is automatic in Times, and Expand keeps the trig atoms intact. */
        Expr* Ediff = ds_call2(SYM_Subtract, Lt, expr_copy((Expr*)T));
        Expr* Eexpr = eval_and_free(ds_call1("Expand",
                          eval_and_free(ds_call2(SYM_Times, Ediff,
                              expr_new_function(expr_new_symbol(SYM_Power),
                                  (Expr*[]){ expr_copy(expfac), expr_new_integer(-1) }, 2)))));

        /* coefficient equations: CoefficientList of the Cos- and Sin-parts in x */
        Expr* eqsList = NULL;
        if (has_trig) {
            Expr* Pc = coeff_of(Eexpr, cosfac);
            Expr* Ps = coeff_of(Eexpr, sinfac);
            /* leftover (a constant/other part) must vanish for the extraction to be valid */
            Expr* leftover = eval_and_free(ds_call1("Expand",
                                 eval_and_free(ds_call2(SYM_Subtract, expr_copy(Eexpr),
                                     ds_call2(SYM_Plus,
                                         ds_call2(SYM_Times, expr_copy(cosfac), expr_copy(Pc)),
                                         ds_call2(SYM_Times, expr_copy(sinfac), expr_copy(Ps)))))));
            bool clean = ds_is_zero(leftover);
            expr_free(leftover);
            if (clean) {
                Expr* cl1 = eval_and_free(ds_call2("CoefficientList", Pc, expr_new_symbol(xvar)));
                Expr* cl2 = eval_and_free(ds_call2("CoefficientList", Ps, expr_new_symbol(xvar)));
                eqsList = eval_and_free(ds_call2("Join", cl1, cl2));
            } else { expr_free(Pc); expr_free(Ps); }
        } else {
            eqsList = eval_and_free(ds_call2("CoefficientList", expr_copy(Eexpr), expr_new_symbol(xvar)));
        }
        expr_free(Eexpr);

        if (eqsList && head_is(eqsList, SYM_List)) {
            /* Map[# == 0 &, eqsList] */
            size_t m = eqsList->data.function.arg_count;
            Expr** eqs = malloc((m ? m : 1) * sizeof(Expr*));
            for (size_t i = 0; i < m; i++)
                eqs[i] = expr_new_function(expr_new_symbol(SYM_Equal),
                             (Expr*[]){ expr_copy(eqsList->data.function.args[i]),
                                        expr_new_integer(0) }, 2);
            Expr* eqsAnd = expr_new_function(expr_new_symbol(SYM_List), eqs, m);
            free(eqs);
            Expr** vcopy = malloc((nv ? nv : 1) * sizeof(Expr*));
            for (size_t i = 0; i < nv; i++) vcopy[i] = expr_copy(vars[i]);
            Expr* varList = expr_new_function(expr_new_symbol(SYM_List), vcopy, nv);
            free(vcopy);
            Expr* sol = eval_and_free(ds_call2("Solve", eqsAnd, varList));
            /* substitute the first solution branch into the trial */
            if (head_is(sol, SYM_List) && sol->data.function.arg_count >= 1
                && head_is(sol->data.function.args[0], SYM_List)) {
                Expr* cand = eval_and_free(ds_call2("ReplaceAll", expr_copy(trial),
                                 expr_copy(sol->data.function.args[0])));
                /* accept only if every undetermined coefficient was determined AND
                 * L[cand] == T is a DECIDABLE zero.  The undetermined-coefficient
                 * linear solve only proves the ansatz template matched; a forcing
                 * term that is NOT a UC function (Q[x], Log[x], Exp[x^2], Cos[x^2],
                 * ...) collapses to degree 0 and is "solved" as y_p = T itself,
                 * whose residual (e.g. T'(x)) the permissive substrate verify
                 * cannot reject.  Requiring a decidable zero here makes those
                 * decline, so linear1's integral form / constcoeff variation of
                 * parameters handles the equation instead of a wrong y_p. */
                if (!ds_contains(cand, uca) && !ds_contains(cand, ucb)) {
                    Expr* cand_s = ds_simplify(cand);
                    Expr* Lc = expr_new_integer(0);
                    for (int j = 0; j <= n; j++) {
                        Expr* dj = expr_copy(cand_s);
                        for (int i = 0; i < j; i++) dj = ds_d(dj, expr_new_symbol(xvar));
                        Lc = eval_and_free(ds_call2(SYM_Plus, Lc,
                                 ds_call2(SYM_Times, expr_copy(a[j]), dj)));
                    }
                    Expr* resid = ds_call2(SYM_Subtract, Lc, expr_copy((Expr*)T));
                    bool solves = ds_is_zero(resid);
                    expr_free(resid);
                    if (solves) yp = cand_s; else expr_free(cand_s);
                } else expr_free(cand);
            }
            expr_free(sol);
        }
        if (eqsList) expr_free(eqsList);
        for (size_t i = 0; i < nv; i++) expr_free(vars[i]);
        free(vars);
        expr_free(trial);
    }

    expr_free(alpha); expr_free(beta); expr_free(expfac);
    if (cosfac) expr_free(cosfac);
    if (sinfac) expr_free(sinfac);
    return yp;
}

Expr** dsolve_undetcoeff_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    const char* xvar = P->ind_names[0];

    Expr** a; Expr* g; int n;
    if (!dsolve_linear_coeffs(P, &a, &g, &n)) return NULL;
    dsolve_linear_normalize(a, &g, n, P->ind_names[0]);

    /* constant coefficients, nonzero leading, nonzero forcing */
    bool ok = !ds_is_zero(g);
    for (int k = 0; k <= n && ok; k++) if (!ds_free_of(a[k], xvar)) ok = false;
    if (ok && ds_is_zero(a[n])) ok = false;
    if (!ok) { for (int k = 0; k <= n; k++) expr_free(a[k]); free(a); expr_free(g); return NULL; }

    /* homogeneous fundamental set from the characteristic polynomial */
    const char* lam = intern_symbol("DSolve`uclam");
    Expr** ch = malloc((size_t)(n + 1) * sizeof(Expr*));
    for (int k = 0; k <= n; k++)
        ch[k] = ds_call2(SYM_Times, expr_copy(a[k]),
                    expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_new_symbol(lam), expr_new_integer(k) }, 2));
    Expr* charpoly = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), ch, (size_t)(n + 1)));
    free(ch);
    size_t nb = 0;
    Expr** basis = dsolve_homog_basis(charpoly, lam, xvar, n, &nb);
    expr_free(charpoly);

    Expr* general = NULL;
    if (basis && nb == (size_t)n) {
        /* particular by superposition over the additive terms of Expand[g] */
        Expr* gE = eval_and_free(ds_call1("Expand", expr_copy(g)));
        size_t nterms; Expr** terms;
        if (head_is(gE, SYM_Plus)) { nterms = gE->data.function.arg_count; terms = gE->data.function.args; }
        else { nterms = 1; terms = &gE; }

        Expr* yp = expr_new_integer(0);
        bool all_ok = true;
        for (size_t t = 0; t < nterms && all_ok; t++) {
            Expr* ypt = uc_particular(terms[t], a, n, xvar);
            if (!ypt) all_ok = false;
            else yp = eval_and_free(ds_call2(SYM_Plus, yp, ypt));
        }
        expr_free(gE);

        if (all_ok) {
            Expr** hterms = malloc((size_t)n * sizeof(Expr*));
            for (int k = 0; k < n; k++)
                hterms[k] = ds_call2(SYM_Times, ds_const(k + 1), expr_copy(basis[k]));
            Expr* Hgen = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), hterms, (size_t)n));
            free(hterms);
            general = eval_and_free(ds_call2(SYM_Plus, Hgen, yp));
        } else {
            expr_free(yp);
        }
    }
    if (basis) { for (size_t i = 0; i < nb; i++) expr_free(basis[i]); free(basis); }
    for (int k = 0; k <= n; k++) expr_free(a[k]);
    free(a); expr_free(g);

    if (!general) return NULL;
    Expr** out = malloc(sizeof(Expr*));
    out[0] = general;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_undetcoeff(Expr* res) {
    return dsolve_method_builtin(res, dsolve_undetcoeff_try);
}

void dsolve_undetcoeff_init(void) {
    symtab_add_builtin("DSolve`UndeterminedCoefficients", builtin_dsolve_undetcoeff);
    symtab_get_def("DSolve`UndeterminedCoefficients")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`UndeterminedCoefficients",
        "DSolve`UndeterminedCoefficients[eqn, y, x] solves a constant-coefficient linear "
        "ODE with polynomial/exponential/sinusoidal forcing by the method of undetermined "
        "coefficients (characteristic-root homogeneous set + a matched trial particular).");
}
