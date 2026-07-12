/* risch_structure.c — the Risch structure theorems (Bronstein §9.3).
 *
 * See risch_structure.h.  The Q-span decision (risch_rational_span) sets up the
 * identity theta == Sum c_i g_i with fresh constant unknowns c_i and hands the
 * coefficient-matching to SolveAlways (which is RowReduce/FLINT-backed); a
 * rational solution certifies membership, its absence certifies independence.
 * The structure-theorem front-ends build the monomial generators (Dt_i for a
 * logarithm, Dt_i/t_i for an exponential) and the target (Da/a or Db) and defer
 * to the span decision.
 */

#include "risch_structure.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"
#include "arithmetic.h"
#include "risch_field.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Evaluation helpers.                                                 */
/* ------------------------------------------------------------------ */

static Expr* rs_cp(const Expr* e) { return expr_copy((Expr*)e); }
static Expr* rs_eval_adopt(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}
static Expr* rs_fn(const char* head, Expr** args, size_t n) {
    return rs_eval_adopt(expr_new_function(expr_new_symbol(head), args, n));
}
static Expr* rs_call1(const char* head, Expr* a) { return rs_fn(head, (Expr*[]){ a }, 1); }
static Expr* rs_call2(const char* head, Expr* a, Expr* b) { return rs_fn(head, (Expr*[]){ a, b }, 2); }
static Expr* rs_times(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){ a, b }, 2);
}
static Expr* rs_pow(Expr* a, long n) {
    return expr_new_function(expr_new_symbol("Power"), (Expr*[]){ a, expr_new_integer(n) }, 2);
}

/* True iff e is a rational-number constant (Integer, BigInt, or Rational[p,q]). */
static bool rs_is_rational(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    return e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == intern_symbol("Rational");
}
/* True iff e is the structural integer 0. */
static bool rs_is_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}
static bool rs_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == intern_symbol("List");
}

/* ------------------------------------------------------------------ */
/* The rational Q-span decision.                                       */
/* ------------------------------------------------------------------ */

Expr* risch_rational_span(const Expr* theta, Expr* const* gens, size_t m,
                          const Expr* vars) {
    /* theta == 0 is in every span with all-zero coefficients. */
    Expr* th0 = rs_call1("Together", rs_cp(theta));
    bool theta_zero = rs_is_zero(th0);
    expr_free(th0);
    if (theta_zero) {
        Expr** zs = (m ? malloc(m * sizeof(Expr*)) : NULL);
        for (size_t i = 0; i < m; i++) zs[i] = expr_new_integer(0);
        Expr* r = expr_new_function(expr_new_symbol("List"), zs, m);
        free(zs);
        return r;
    }
    if (m == 0) return NULL;   /* nonzero theta cannot be an empty combination */

    /* Fresh constant unknowns c_i and the combination Sum c_i g_i. */
    Expr** cs = malloc(m * sizeof(Expr*));
    Expr** terms = malloc(m * sizeof(Expr*));
    for (size_t i = 0; i < m; i++) {
        char name[32];
        snprintf(name, sizeof name, "$RSpan$%zu", i);
        cs[i] = expr_new_symbol(name);
        terms[i] = rs_times(rs_cp(cs[i]), rs_cp(gens[i]));
    }
    Expr* sum = expr_new_function(expr_new_symbol("Plus"), terms, m);
    free(terms);
    /* Clear denominators so SolveAlways sees a polynomial identity:
     * Numerator[Together[theta - Sum c_i g_i]] == 0. */
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ rs_cp(theta),
                   expr_new_function(expr_new_symbol("Times"),
                       (Expr*[]){ expr_new_integer(-1), sum }, 2) }, 2);
    Expr* num = rs_call1("Numerator", rs_call1("Together", diff));
    Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
        (Expr*[]){ num, expr_new_integer(0) }, 2);   /* num == 0 */
    /* Dependent generators make the system underdetermined by design; mute the
     * Solve::svars informational message during this internal probe. */
    arith_warnings_mute_push();
    Expr* sol = rs_call2("SolveAlways", eqn, rs_cp(vars));
    arith_warnings_mute_pop();

    /* sol is a List of solution-rule-lists; {} means no solution (theta != 0 so
     * it is a genuine contradiction, not the always-true degenerate case). */
    Expr* result = NULL;
    if (rs_is_list(sol) && sol->data.function.arg_count >= 1) {
        Expr* first = sol->data.function.args[0];   /* List of Rule[c_i, val] */
        Expr** res = malloc(m * sizeof(Expr*));
        for (size_t i = 0; i < m; i++) res[i] = NULL;
        if (rs_is_list(first)) {
            for (size_t k = 0; k < first->data.function.arg_count; k++) {
                Expr* rule = first->data.function.args[k];
                if (rule->type == EXPR_FUNCTION &&
                    rule->data.function.head->type == EXPR_SYMBOL &&
                    rule->data.function.head->data.symbol == intern_symbol("Rule") &&
                    rule->data.function.arg_count == 2) {
                    Expr* lhs = rule->data.function.args[0];
                    for (size_t i = 0; i < m; i++)
                        if (lhs->type == EXPR_SYMBOL && cs[i]->type == EXPR_SYMBOL &&
                            lhs->data.symbol == cs[i]->data.symbol) {
                            res[i] = rs_cp(rule->data.function.args[1]);
                        }
                }
            }
        }
        /* Free parameters (unconstrained c_i) get pinned to 0; also substitute
         * any residual c_j appearing in another coefficient's value. */
        Expr** subrules = malloc(m * sizeof(Expr*));
        for (size_t i = 0; i < m; i++)
            subrules[i] = expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ rs_cp(cs[i]), expr_new_integer(0) }, 2);
        Expr* srlist = expr_new_function(expr_new_symbol("List"), subrules, m);
        free(subrules);

        bool ok = true;
        for (size_t i = 0; i < m; i++) {
            if (!res[i]) res[i] = expr_new_integer(0);
            Expr* v = rs_eval_adopt(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ res[i], rs_cp(srlist) }, 2));   /* adopts res[i] */
            res[i] = v;
            if (!rs_is_rational(res[i])) ok = false;
        }
        expr_free(srlist);
        if (ok) {
            result = expr_new_function(expr_new_symbol("List"), res, m);
        } else {
            for (size_t i = 0; i < m; i++) expr_free(res[i]);
        }
        free(res);
    }
    expr_free(sol);
    for (size_t i = 0; i < m; i++) expr_free(cs[i]);
    free(cs);
    return result;
}

/* ------------------------------------------------------------------ */
/* Structure-theorem front-ends.                                       */
/* ------------------------------------------------------------------ */

/* Decode {{t_i, "Exp"|"Log", Dt_i}, ...} into parallel arrays and build both the
 * derivation rules {x->1, t_i->Dt_i} and the generator list.  Returns the number
 * of monomials, or (size_t)-1 on malformed input.  On success the caller owns
 * *gens (array + elements), *vars (a List Expr), and *deriv_rules (a List Expr). */
static size_t rs_decode_tower(const Expr* x, const Expr* mons,
                              Expr*** gens_out, Expr** vars_out, Expr** deriv_out) {
    if (!rs_is_list(mons)) return (size_t)-1;
    size_t n = mons->data.function.arg_count;

    Expr** gens = (n ? malloc(n * sizeof(Expr*)) : NULL);
    /* derivation rules: {x->1} plus one per monomial. */
    Expr** drules = malloc((n + 1) * sizeof(Expr*));
    Expr** vars = malloc((n + 1) * sizeof(Expr*));
    drules[0] = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ rs_cp(x), expr_new_integer(1) }, 2);
    vars[0] = rs_cp(x);

    size_t done = 0;   /* number of monomials fully processed */
    bool ok = true;
    for (size_t i = 0; i < n; i++) {
        Expr* mi = mons->data.function.args[i];
        if (!rs_is_list(mi) || mi->data.function.arg_count != 3) { ok = false; break; }
        Expr* ti = mi->data.function.args[0];
        Expr* kind = mi->data.function.args[1];
        Expr* Dti = mi->data.function.args[2];
        if (kind->type != EXPR_STRING) { ok = false; break; }
        const char* ks = kind->data.string;

        /* generator: Log -> Dt_i ; Exp -> Dt_i / t_i ; Tan -> Dt_i / (t_i^2+1) */
        Expr* g;
        if (strcmp(ks, "Log") == 0) {
            g = rs_cp(Dti);
        } else if (strcmp(ks, "Exp") == 0) {
            g = rs_call1("Cancel", rs_times(rs_cp(Dti), rs_pow(rs_cp(ti), -1)));
        } else if (strcmp(ks, "Tan") == 0) {
            Expr* tsq1 = rs_eval_adopt(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ rs_pow(rs_cp(ti), 2), expr_new_integer(1) }, 2));
            g = rs_call1("Cancel", rs_times(rs_cp(Dti), rs_pow(tsq1, -1)));
        } else { ok = false; break; }

        gens[i] = g;
        drules[i + 1] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ rs_cp(ti), rs_cp(Dti) }, 2);
        vars[i + 1] = rs_cp(ti);
        done++;
    }

    if (!ok) {
        for (size_t i = 0; i < done; i++) expr_free(gens[i]);
        for (size_t i = 0; i <= done; i++) { expr_free(drules[i]); expr_free(vars[i]); }
        free(gens); free(drules); free(vars);
        return (size_t)-1;
    }

    *gens_out = gens;
    *vars_out = expr_new_function(expr_new_symbol("List"), vars, n + 1);
    free(vars);
    *deriv_out = expr_new_function(expr_new_symbol("List"), drules, n + 1);
    free(drules);
    return n;
}

/* Shared body for LogReducible / ExpReducible.  is_log selects the target:
 * log test -> Da/a ; exp test -> Db. */
static Expr* rs_reducible(const Expr* arg, const Expr* x, const Expr* mons, bool is_log) {
    Expr** gens = NULL; Expr* vars = NULL; Expr* drules = NULL;
    size_t m = rs_decode_tower(x, mons, &gens, &vars, &drules);
    if (m == (size_t)-1) return NULL;

    RischDeriv d;
    if (!risch_deriv_from_rules(drules, &d)) {
        for (size_t i = 0; i < m; i++) expr_free(gens[i]);
        free(gens); expr_free(vars); expr_free(drules);
        return NULL;
    }

    Expr* Darg = risch_field_deriv(arg, &d);   /* D_tower[arg] */
    Expr* theta = is_log
        ? rs_call1("Cancel", rs_times(Darg, rs_pow(rs_cp(arg), -1)))  /* Da/a */
        : Darg;                                                       /* Db */

    Expr* coeffs = risch_rational_span(theta, gens, m, vars);

    risch_deriv_free(&d);
    expr_free(theta);
    for (size_t i = 0; i < m; i++) expr_free(gens[i]);
    free(gens); expr_free(vars); expr_free(drules);

    return coeffs ? coeffs : expr_new_symbol("False");
}

/* ------------------------------------------------------------------ */
/* Builtins.                                                           */
/* ------------------------------------------------------------------ */

/* Risch`RationalSpan[theta, {g...}, {vars}] -> {r...} or False. */
static Expr* builtin_risch_rationalspan(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* theta = res->data.function.args[0];
    Expr* glist = res->data.function.args[1];
    Expr* vars = res->data.function.args[2];
    if (!rs_is_list(glist) || !rs_is_list(vars)) return NULL;
    size_t m = glist->data.function.arg_count;
    Expr* r = risch_rational_span(theta, glist->data.function.args, m, vars);
    return r ? r : expr_new_symbol("False");
}

/* Risch`LogReducible[a, x, mons] / Risch`ExpReducible[b, x, mons]. */
static Expr* builtin_risch_logreducible(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    return rs_reducible(res->data.function.args[0], res->data.function.args[1],
                        res->data.function.args[2], true);
}
static Expr* builtin_risch_expreducible(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    return rs_reducible(res->data.function.args[0], res->data.function.args[1],
                        res->data.function.args[2], false);
}

/* ------------------------------------------------------------------ */
/* Registration.                                                       */
/* ------------------------------------------------------------------ */

static void rs_install(const char* name, Expr* (*fn)(Expr*), const char* doc) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    if (doc) symtab_set_docstring(name, doc);
}

void integrate_risch_structure_init(void) {
    rs_install("Risch`RationalSpan", builtin_risch_rationalspan,
        "Risch`RationalSpan[theta, {g1,...}, {vars}] returns {r1,...} of rational\n"
        "numbers with theta = Sum r_i g_i (coefficients matched over vars), or\n"
        "False if theta is not in the rational span of the generators.");
    rs_install("Risch`LogReducible", builtin_risch_logreducible,
        "Risch`LogReducible[a, x, {{t, \"Exp\"|\"Log\", Dt}, ...}] applies the Risch\n"
        "structure theorem: returns the rational coefficients expressing Da/a in the\n"
        "monomial generators when Log[a] is reducible over the tower, else False\n"
        "(Log[a] is a new, algebraically independent monomial).");
    rs_install("Risch`ExpReducible", builtin_risch_expreducible,
        "Risch`ExpReducible[b, x, {{t, \"Exp\"|\"Log\", Dt}, ...}] applies the Risch\n"
        "structure theorem: returns the rational coefficients expressing Db in the\n"
        "monomial generators when Exp[b] is reducible over the tower, else False\n"
        "(Exp[b] is a new, algebraically independent monomial).");
}
