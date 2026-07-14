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
           e->data.function.head->data.symbol.name == intern_symbol("Rational");
}
/* True iff e is the structural integer 0. */
static bool rs_is_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}
static bool rs_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == intern_symbol("List");
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
                    rule->data.function.head->data.symbol.name == intern_symbol("Rule") &&
                    rule->data.function.arg_count == 2) {
                    Expr* lhs = rule->data.function.args[0];
                    for (size_t i = 0; i < m; i++)
                        if (lhs->type == EXPR_SYMBOL && cs[i]->type == EXPR_SYMBOL &&
                            lhs->data.symbol.name == cs[i]->data.symbol.name) {
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

/* Monomial kinds and their structure-theorem index-set group.  The REAL Risch
 * structure theorem (Bronstein Cor. 9.3.2) partitions the tower monomials into
 * four DISJOINT index sets — E (exponentials), L (logarithms), T (tangents),
 * A (arc-tangents) — and the four reducibility decisions each use only the
 * relevant PAIR of generator sets:
 *   log/exp reducibility (eq. 9.12/9.13):     generators from E u L (group 0);
 *   tan/arctan reducibility (eq. 9.14/9.15):  generators from T u A (group 1).
 * The generator of a monomial t_i is  Dt_i (Log, ArcTan),  Dt_i/t_i (Exp),
 * Dt_i/(t_i^2+1) (Tan).  Lumping all generators (as the complex-only code did)
 * is correct ONLY while no tangent monomials are present; the partitioning below
 * makes the decision faithful to the disjoint-index real theorem. */
enum { RS_LOG = 0, RS_EXP = 1, RS_TAN = 2, RS_ARCTAN = 3 };
static int rs_kind_group(int k) { return (k == RS_TAN || k == RS_ARCTAN) ? 1 : 0; }

/* Target-decision codes (which of the four Cor. 9.3.1/9.3.2 questions to answer). */
enum { RS_LOG_R = 0, RS_EXP_R = 1, RS_TAN_R = 2, RS_ARCTAN_R = 3 };
static int rs_target_group(int t) { return (t == RS_TAN_R || t == RS_ARCTAN_R) ? 1 : 0; }

/* Decode {{t_i, "Exp"|"Log"|"Tan"|"ArcTan", Dt_i [, arg_i]}, ...} into parallel
 * arrays and build the derivation rules {x->1, t_i->Dt_i} and the generator list.
 * Returns the number of monomials, or (size_t)-1 on malformed input.  On success
 * the caller owns *gens (array + elements), *kinds (int array), *vars (List Expr),
 * *deriv_out (List Expr), and — when bases_out != NULL — *bases (array + elements,
 * the witness base for each monomial: the log ARGUMENT arg_i when a 4-tuple gives
 * it, else the monomial t_i itself).  The optional 4th element carries a_i for a
 * logarithm t_i = log(a_i), used only for radical-witness reconstruction. */
static size_t rs_decode_tower(const Expr* x, const Expr* mons,
                              Expr*** gens_out, int** kinds_out, Expr*** bases_out,
                              Expr** vars_out, Expr** deriv_out) {
    if (!rs_is_list(mons)) return (size_t)-1;
    size_t n = mons->data.function.arg_count;

    Expr** gens = (n ? malloc(n * sizeof(Expr*)) : NULL);
    int*   kinds = (n ? malloc(n * sizeof(int)) : NULL);
    Expr** bases = (bases_out && n ? malloc(n * sizeof(Expr*)) : NULL);
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
        if (!rs_is_list(mi) ||
            (mi->data.function.arg_count != 3 && mi->data.function.arg_count != 4)) {
            ok = false; break;
        }
        Expr* ti = mi->data.function.args[0];
        Expr* kind = mi->data.function.args[1];
        Expr* Dti = mi->data.function.args[2];
        Expr* argi = (mi->data.function.arg_count == 4) ? mi->data.function.args[3] : NULL;
        if (kind->type != EXPR_STRING) { ok = false; break; }
        const char* ks = kind->data.string;

        /* generator by kind; group per the disjoint index sets. */
        Expr* g; int kc;
        if (strcmp(ks, "Log") == 0) {
            kc = RS_LOG; g = rs_cp(Dti);
        } else if (strcmp(ks, "Exp") == 0) {
            kc = RS_EXP; g = rs_call1("Cancel", rs_times(rs_cp(Dti), rs_pow(rs_cp(ti), -1)));
        } else if (strcmp(ks, "Tan") == 0) {
            kc = RS_TAN;
            Expr* tsq1 = rs_eval_adopt(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ rs_pow(rs_cp(ti), 2), expr_new_integer(1) }, 2));
            g = rs_call1("Cancel", rs_times(rs_cp(Dti), rs_pow(tsq1, -1)));
        } else if (strcmp(ks, "ArcTan") == 0) {
            kc = RS_ARCTAN; g = rs_cp(Dti);
        } else { ok = false; break; }

        gens[i] = g;
        kinds[i] = kc;
        if (bases) bases[i] = (kc == RS_LOG && argi) ? rs_cp(argi) : rs_cp(ti);
        drules[i + 1] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ rs_cp(ti), rs_cp(Dti) }, 2);
        vars[i + 1] = rs_cp(ti);
        done++;
    }

    if (!ok) {
        for (size_t i = 0; i < done; i++) { expr_free(gens[i]); if (bases) expr_free(bases[i]); }
        for (size_t i = 0; i <= done; i++) { expr_free(drules[i]); expr_free(vars[i]); }
        free(gens); free(kinds); free(bases); free(drules); free(vars);
        return (size_t)-1;
    }

    *gens_out = gens;
    *kinds_out = kinds;
    if (bases_out) *bases_out = bases;
    *vars_out = expr_new_function(expr_new_symbol("List"), vars, n + 1);
    free(vars);
    *deriv_out = expr_new_function(expr_new_symbol("List"), drules, n + 1);
    free(drules);
    return n;
}

/* Shared body for the four structure-theorem reducibility decisions.  `target`
 * (RS_*_R) selects both the target theta and the generator group:
 *   RS_LOG_R:    theta = Da/a       gens E u L   (Cor. i,   eq. 9.8/9.12)
 *   RS_EXP_R:    theta = Db         gens E u L   (Cor. ii,  eq. 9.9/9.13)
 *   RS_TAN_R:    theta = Db         gens T u A   (Cor. iv,  eq. 9.15)
 *   RS_ARCTAN_R: theta = Db/(b^2+1) gens T u A   (Cor. iii, eq. 9.14)
 * Returns an owned List of rational coefficients in FULL monomial-list order
 * (0 at monomials outside the target's group), or NULL to signal "False". */
static Expr* rs_structure_decide(const Expr* arg, const Expr* x, const Expr* mons,
                                 int target) {
    Expr** gens = NULL; int* kinds = NULL; Expr* vars = NULL; Expr* drules = NULL;
    size_t m = rs_decode_tower(x, mons, &gens, &kinds, NULL, &vars, &drules);
    if (m == (size_t)-1) return NULL;

    RischDeriv d;
    if (!risch_deriv_from_rules(drules, &d)) {
        for (size_t i = 0; i < m; i++) expr_free(gens[i]);
        free(gens); free(kinds); expr_free(vars); expr_free(drules);
        return NULL;
    }

    Expr* Darg = risch_field_deriv(arg, &d);   /* D_tower[arg] */
    Expr* theta;
    switch (target) {
        case RS_LOG_R:
            theta = rs_call1("Cancel", rs_times(Darg, rs_pow(rs_cp(arg), -1)));   /* Da/a */
            break;
        case RS_ARCTAN_R: {                                                       /* Db/(b^2+1) */
            Expr* bsq1 = rs_eval_adopt(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ rs_pow(rs_cp(arg), 2), expr_new_integer(1) }, 2));
            theta = rs_call1("Cancel", rs_times(Darg, rs_pow(bsq1, -1)));
            break;
        }
        default: /* RS_EXP_R, RS_TAN_R */
            theta = Darg;                                                         /* Db */
            break;
    }

    /* Restrict to the generators of the target's index-set group. */
    int grp = rs_target_group(target);
    Expr** fg = (m ? malloc(m * sizeof(Expr*)) : NULL);
    size_t* idx = (m ? malloc(m * sizeof(size_t)) : NULL);
    size_t fm = 0;
    for (size_t i = 0; i < m; i++)
        if (rs_kind_group(kinds[i]) == grp) { fg[fm] = gens[i]; idx[fm] = i; fm++; }

    Expr* coeffs = risch_rational_span(theta, fg, fm, vars);

    Expr* result = NULL;
    if (coeffs) {
        Expr** full = (m ? malloc(m * sizeof(Expr*)) : NULL);
        for (size_t i = 0; i < m; i++) full[i] = expr_new_integer(0);
        for (size_t j = 0; j < fm; j++) {
            expr_free(full[idx[j]]);
            full[idx[j]] = rs_cp(coeffs->data.function.args[j]);
        }
        result = expr_new_function(expr_new_symbol("List"), full, m);
        free(full);
        expr_free(coeffs);
    }

    risch_deriv_free(&d);
    expr_free(theta);
    for (size_t i = 0; i < m; i++) expr_free(gens[i]);
    free(gens); free(kinds); free(fg); free(idx);
    expr_free(vars); expr_free(drules);

    return result;
}

/* ------------------------------------------------------------------ */
/* Logarithmic derivative of a radical (Bronstein §5.12 / eq. 7.44).   */
/* ------------------------------------------------------------------ */

/* Decide whether f is the logarithmic derivative of a K-radical over the tower,
 * via the structure-theorem test (Cor. 9.3.1/9.3.2 (ii), eq. 9.13/7.44): f is
 * such a logarithmic derivative iff there are r_i in Q with
 *     f = Sum_{i in L} r_i Dt_i + Sum_{i in E} r_i Dt_i/t_i
 * (only the E u L generators participate).  On success, put the r_i over a common
 * denominator n>0, set e_i = n*r_i in Z, and reconstruct the witness radical
 *     u = Prod base_i^{e_i}     (base_i = arg_i for a logarithm t_i = log(a_i),
 *                                = t_i for an exponential),
 * so that D[u]/u = n*f exactly.  Returns an owned List[n, u], or NULL for "False".
 *
 * This is the exact, complete decision when f = Db is the derivative of a field
 * element (the form arising from integration); the full §5.12 recursive method
 * (which also admits radicals built from the log MONOMIALS themselves, e.g.
 * Bronstein's u = x^5 log(x)) is a strictly broader completeness item. */
static Expr* rs_logderiv_radical(const Expr* f, const Expr* x, const Expr* mons) {
    Expr** gens = NULL; int* kinds = NULL; Expr** bases = NULL;
    Expr* vars = NULL; Expr* drules = NULL;
    size_t m = rs_decode_tower(x, mons, &gens, &kinds, &bases, &vars, &drules);
    if (m == (size_t)-1) return NULL;

    /* E u L generators only (group 0). */
    Expr** fg = (m ? malloc(m * sizeof(Expr*)) : NULL);
    size_t* idx = (m ? malloc(m * sizeof(size_t)) : NULL);
    size_t fm = 0;
    for (size_t i = 0; i < m; i++)
        if (rs_kind_group(kinds[i]) == 0) { fg[fm] = gens[i]; idx[fm] = i; fm++; }

    Expr* coeffs = risch_rational_span(f, fg, fm, vars);

    Expr* result = NULL;
    if (coeffs) {
        /* n = LCM of the coefficient denominators (n = 1 when all integral / empty). */
        Expr* n;
        if (fm == 0) {
            n = expr_new_integer(1);
        } else {
            Expr** ds = malloc(fm * sizeof(Expr*));
            for (size_t j = 0; j < fm; j++)
                ds[j] = rs_call1("Denominator", rs_cp(coeffs->data.function.args[j]));
            Expr* dl = expr_new_function(expr_new_symbol("List"), ds, fm);
            free(ds);
            n = rs_fn("Apply", (Expr*[]){ expr_new_symbol("LCM"), dl }, 2);
        }

        /* u = Prod base_i^{n r_i}, skipping zero exponents. */
        Expr** facs = (fm ? malloc(fm * sizeof(Expr*)) : NULL);
        size_t nf = 0;
        for (size_t j = 0; j < fm; j++) {
            Expr* ej = rs_eval_adopt(rs_times(rs_cp(n), rs_cp(coeffs->data.function.args[j])));
            if (rs_is_zero(ej)) { expr_free(ej); continue; }
            facs[nf++] = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ rs_cp(bases[idx[j]]), ej }, 2);
        }
        Expr* u;
        if (nf == 0) u = expr_new_integer(1);
        else if (nf == 1) u = rs_eval_adopt(facs[0]);
        else u = rs_eval_adopt(expr_new_function(expr_new_symbol("Times"), facs, nf));
        free(facs);

        result = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){ n, u }, 2);
    }

    if (coeffs) expr_free(coeffs);
    for (size_t i = 0; i < m; i++) { expr_free(gens[i]); if (bases) expr_free(bases[i]); }
    free(gens); free(kinds); free(bases); free(fg); free(idx);
    expr_free(vars); expr_free(drules);

    return result;
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

/* The four structure-theorem reducibility builtins share one shape. */
static Expr* rs_reducible_builtin(Expr* res, int target) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* r = rs_structure_decide(res->data.function.args[0], res->data.function.args[1],
                                  res->data.function.args[2], target);
    return r ? r : expr_new_symbol("False");
}
/* Risch`LogReducible[a, x, mons] / Risch`ExpReducible[b, x, mons]
 * / Risch`TanReducible[b, x, mons] / Risch`ArcTanReducible[b, x, mons]. */
static Expr* builtin_risch_logreducible(Expr* res)    { return rs_reducible_builtin(res, RS_LOG_R); }
static Expr* builtin_risch_expreducible(Expr* res)    { return rs_reducible_builtin(res, RS_EXP_R); }
static Expr* builtin_risch_tanreducible(Expr* res)    { return rs_reducible_builtin(res, RS_TAN_R); }
static Expr* builtin_risch_arctanreducible(Expr* res) { return rs_reducible_builtin(res, RS_ARCTAN_R); }

/* Risch`LogarithmicDerivativeOfRadical[f, x, mons] -> {n, u} or False. */
static Expr* builtin_risch_logderiv_radical(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* r = rs_logderiv_radical(res->data.function.args[0], res->data.function.args[1],
                                  res->data.function.args[2]);
    return r ? r : expr_new_symbol("False");
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
    rs_install("Risch`TanReducible", builtin_risch_tanreducible,
        "Risch`TanReducible[b, x, {{t, \"Tan\"|\"ArcTan\"|\"Exp\"|\"Log\", Dt}, ...}]\n"
        "applies the REAL Risch structure theorem (Cor. 9.3.2 iv, eq. 9.15): returns\n"
        "the rational coefficients over the tangent/arc-tangent generators when Tan[b]\n"
        "is reducible over the tower, else False. Only T u A monomials participate.");
    rs_install("Risch`ArcTanReducible", builtin_risch_arctanreducible,
        "Risch`ArcTanReducible[b, x, {{t, \"Tan\"|\"ArcTan\"|\"Exp\"|\"Log\", Dt}, ...}]\n"
        "applies the REAL Risch structure theorem (Cor. 9.3.2 iii, eq. 9.14): returns\n"
        "the rational coefficients expressing Db/(b^2+1) over the tangent/arc-tangent\n"
        "generators when ArcTan[b] is reducible over the tower, else False.");
    rs_install("Risch`LogarithmicDerivativeOfRadical", builtin_risch_logderiv_radical,
        "Risch`LogarithmicDerivativeOfRadical[f, x, {{t, kind, Dt[, arg]}, ...}]\n"
        "decides whether f is the logarithmic derivative of a K-radical over the tower\n"
        "(Bronstein §5.12, structure-theorem test eq. 7.44/9.13): returns {n, u} with\n"
        "n a positive integer and u the witness radical such that D[u]/u == n f, else\n"
        "False. Only the Exp/Log (E u L) monomials participate; supply the log argument\n"
        "as the optional 4th monomial element for the witness u.");
}
