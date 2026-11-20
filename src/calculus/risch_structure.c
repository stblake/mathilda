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
#include "risch_util.h"
#include "context.h"

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
/* An explicit algebraic number: integers, rationals, Complex of those, and
 * Plus / Times / Power with a rational exponent / Root[poly &, k] over them. */
static bool rs_is_algebraic_number(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    size_t n = e->data.function.arg_count;
    if (h == intern_symbol("Rational") || h == intern_symbol("Complex")) {
        for (size_t i = 0; i < n; i++)
            if (!rs_is_algebraic_number(e->data.function.args[i])) return false;
        return true;
    }
    if (h == intern_symbol("Root")) return true;   /* Root[f &, k]: algebraic by construction */
    if (h == intern_symbol("Power") && n == 2) {
        const Expr* ex = e->data.function.args[1];
        bool ratexp = ex->type == EXPR_INTEGER || ex->type == EXPR_BIGINT
            || (ex->type == EXPR_FUNCTION && ex->data.function.head->type == EXPR_SYMBOL
                && ex->data.function.head->data.symbol.name == intern_symbol("Rational"));
        return ratexp && rs_is_algebraic_number(e->data.function.args[0]);
    }
    if (h == intern_symbol("Plus") || h == intern_symbol("Times")) {
        for (size_t i = 0; i < n; i++)
            if (!rs_is_algebraic_number(e->data.function.args[i])) return false;
        return true;
    }
    return false;
}

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

/* Symbolic parameters: rt_collect_params (risch_util.c), shared with the
 * ansatz sites' rt_solve_ansatz. */

/* ------------------------------------------------------------------ */
/* The rational Q-span decision.                                       */
/* ------------------------------------------------------------------ */

/* The Gaussian constants of e.  Returns 1 when e carries a Complex[a, b] with
 * rational a and b, -1 when it carries a Complex of anything else (a machine
 * or MPFR part), 0 when it carries none. */
static int rs_gauss_kind(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    int k = 0;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Complex")) {
        if (e->data.function.arg_count != 2) return -1;
        return (rs_is_rational(e->data.function.args[0])
                && rs_is_rational(e->data.function.args[1])) ? 1 : -1;
    }
    int kh = rs_gauss_kind(e->data.function.head);
    if (kh < 0) return -1;
    if (kh) k = 1;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int ki = rs_gauss_kind(e->data.function.args[i]);
        if (ki < 0) return -1;
        if (ki) k = 1;
    }
    return k;
}
/* True when every constant of e is real, so that iota^2 + 1 stays
 * irreducible over the field they generate: numbers, symbols (variables,
 * parameters, Pi, E), Complex (spelled with iota), and Plus / Times / Power
 * / Rational / Equal over those, a Power with a non-integer exponent only of
 * a positive rational base.  Anything else -- Sqrt[-3] = I Sqrt[3] in another
 * spelling, a Root, Log[-2], an unknown head -- is not known to be real, and
 * splitting at iota there could miss a relation. */
static bool rs_real_constants(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    size_t n = e->data.function.arg_count;
    if (h == intern_symbol("Complex") || h == intern_symbol("Rational")) return true;
    if (h == intern_symbol("Power") && n == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* ex = e->data.function.args[1];
        bool intexp = ex->type == EXPR_INTEGER || ex->type == EXPR_BIGINT;
        /* E^(I x): iota inside an exponent is not polynomial arithmetic */
        if (!intexp && rs_gauss_kind(ex) != 0) return false;
        if (!intexp && (b->type == EXPR_INTEGER || b->type == EXPR_BIGINT || rs_is_rational(b))) {
            Expr* pos = rs_call1("Positive", rs_cp(b));
            bool ok = pos && pos->type == EXPR_SYMBOL && pos->data.symbol.name == intern_symbol("True");
            if (pos) expr_free(pos);
            if (!ok) return false;
            return rs_real_constants(ex);
        }
        if (!intexp && b->type != EXPR_SYMBOL) {
            /* a non-integer power of a compound base is real only for a
             * positive base, which is not checked here: a variable base
             * (a tower symbol, a parameter) is the identity's business */
            return false;
        }
        return rs_real_constants(b) && rs_real_constants(ex);
    }
    if (h == intern_symbol("Plus") || h == intern_symbol("Times") || h == intern_symbol("Equal")
        || h == intern_symbol("List")) {
        for (size_t i = 0; i < n; i++)
            if (!rs_real_constants(e->data.function.args[i])) return false;
        return true;
    }
    return false;
}
/* e with every Gaussian constant Complex[a, b] spelled a + b iota (owned,
 * unevaluated). */
static Expr* rs_gauss_subst(const Expr* e, const Expr* iota) {
    if (!e || e->type != EXPR_FUNCTION) return rs_cp(e);
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Complex")
        && e->data.function.arg_count == 2)
        return expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
            rs_cp(e->data.function.args[0]),
            rs_times(rs_cp(e->data.function.args[1]), rs_cp(iota)) }, 2);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = rs_gauss_subst(e->data.function.args[i], iota);
    Expr* r = expr_new_function(rs_gauss_subst(e->data.function.head, iota), a, n);
    free(a);
    return r;
}

Expr* risch_rational_span(const Expr* theta, Expr* const* gens, size_t m,
                          const Expr* vars, bool* failed) {
    *failed = false;
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
    /* Gaussian constants.  Matching coefficients over C(vars) finds the
     * C-span, not the Q-span: -I over {I, 1} is c_0 I + c_1 = -I, a line
     * of complex solutions whose pinned point c_1 = -I is not rational, so
     * the decision gave up (gauss-elem-27: E^(-I (1 + x)) beside E^x and
     * E^(-2 I + I x), every exponential tower over Q(i) stayed "dependent,
     * undecided").  With the constants in Q(i), spell I as a fresh symbol
     * iota: an identity in Q(i)(vars) with rational c_i holds iff the
     * numerator, reduced mod iota^2 + 1, vanishes as a polynomial in vars
     * and iota -- its iota^0 and iota^1 parts separately.  (The denominator
     * at iota = I is the original nonzero one.)  The system is then over Q,
     * its solution space is Q-defined, and an empty one is a proof. */
    Expr* iota = NULL;
    if (rs_is_list(vars)) {
        int gk = rs_gauss_kind(theta);
        for (size_t i = 0; i < m && gk >= 0; i++) {
            int g = rs_gauss_kind(gens[i]);
            if (g < 0) gk = -1; else if (g) gk = 1;
        }
        bool real = gk == 1 && rs_real_constants(theta);
        for (size_t i = 0; i < m && real; i++) real = rs_real_constants(gens[i]);
        if (real) iota = expr_new_symbol("$RSpanIota");
    }
    Expr* num;
    bool split_bad = false;   /* the reduction mod iota^2 + 1 did not come out linear */
    if (iota) {
        Expr* dg = rs_eval_adopt(rs_gauss_subst(diff, iota));
        expr_free(diff);
        Expr* nm = rs_call1("Numerator", rs_call1("Together", dg));
        Expr* mp = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
            rs_pow(rs_cp(iota), 2), expr_new_integer(1) }, 2);
        num = rs_call1("Expand", rs_fn("PolynomialRemainder",
            (Expr*[]){ rs_call1("Expand", nm), mp, rs_cp(iota) }, 3));
        /* only a remainder that is a polynomial of degree <= 1 in iota is
         * the split; an empty solution set of anything else is no proof */
        Expr* pq = num ? rs_call2("PolynomialQ", rs_cp(num), rs_cp(iota)) : NULL;
        Expr* ex = num ? rs_call2("Exponent", rs_cp(num), rs_cp(iota)) : NULL;
        split_bad = !(pq && pq->type == EXPR_SYMBOL && pq->data.symbol.name == intern_symbol("True")
                      && ex && ex->type == EXPR_INTEGER && ex->data.integer <= 1);
        if (pq) expr_free(pq);
        if (ex) expr_free(ex);
        if (!num) num = expr_new_symbol("$Failed");
    } else {
        num = rs_call1("Numerator", rs_call1("Together", diff));
    }
    Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
        (Expr*[]){ num, expr_new_integer(0) }, 2);   /* num == 0 */
    /* Dependent generators make the system underdetermined by design; mute the
     * Solve::svars informational message during this internal probe. */
    /* The identity must hold for ALL values of the symbolic parameters too, so
     * they join the SolveAlways variables.  The structure-theorem multipliers
     * r_i are rational NUMBERS (Bronstein Thm 9.3.x), never functions of the
     * parameters: left out of `vars`, a parameter c is solved for as one more
     * unknown beside the c_i, so E^(-c x) over E^(c x) gave the underdetermined
     * -c - r c == 0 (Solve::nsdim, "not decided") instead of r = -1, and x over
     * E^(c x) could equally have produced the bogus r = 1/c. */
    Expr* allvars = rs_cp(vars);
    bool nonrat = false;   /* a parameter outside a rational position */
    if (rs_is_list(vars)) {
        const char** ps = NULL; size_t np = 0, pcap = 0;
        rt_collect_params(theta, vars, true, &nonrat, &ps, &np, &pcap);
        for (size_t i = 0; i < m; i++)
            rt_collect_params(gens[i], vars, true, &nonrat, &ps, &np, &pcap);
        size_t ni = iota ? 1 : 0;
        if (np || ni) {
            size_t nv = vars->data.function.arg_count;
            Expr** vv = malloc((nv + np + ni) * sizeof(Expr*));
            for (size_t i = 0; i < nv; i++) vv[i] = rs_cp(vars->data.function.args[i]);
            for (size_t i = 0; i < np; i++) vv[nv + i] = expr_new_symbol(ps[i]);
            if (ni) vv[nv + np] = rs_cp(iota);
            expr_free(allvars);
            allvars = expr_new_function(expr_new_symbol("List"), vv, nv + np + ni);
            free(vv);
        }
        free(ps);
    }
    if (iota) expr_free(iota);
    arith_warnings_mute_push();
    Expr* sol = rs_call2("SolveAlways", eqn, allvars);
    arith_warnings_mute_pop();

    /* sol is a List of solution-rule-lists; {} means no solution (theta != 0 so
     * it is a genuine contradiction, not the always-true degenerate case) -- the
     * only NULL that is a proof.  Anything that is not a List (SolveAlways left
     * unevaluated or gave up) is an undecided answer. */
    Expr* result = NULL;
    if (!rs_is_list(sol)) *failed = true;
    /* ...unless a parameter sits outside a rational position, where SolveAlways
     * treats Sqrt[c] / Log[c] as extra unknowns: then {} is not a proof. */
    if (rs_is_list(sol) && sol->data.function.arg_count == 0 && (nonrat || split_bad)) *failed = true;
    if (rs_is_list(sol) && sol->data.function.arg_count >= 1) {
        Expr* first = sol->data.function.args[0];   /* List of Rule[c_i, val] */
        Expr** res = malloc(m * sizeof(Expr*));
        for (size_t i = 0; i < m; i++) res[i] = NULL;
        /* `unique`: one solution, every rule binds some c_i, every c_i is
         * bound to a numeric constant free of the c's -- the solution point is
         * then the only one (the system in the c's has full rank), so a
         * non-rational value at it PROVES there is no rational relation
         * (review F9: 2 over {Sqrt[2]}, E^(Sqrt[2] u) beside E^(2 u)). */
        bool unique = sol->data.function.arg_count == 1 && rs_is_list(first) && !nonrat;
        if (rs_is_list(first)) {
            for (size_t k = 0; k < first->data.function.arg_count; k++) {
                Expr* rule = first->data.function.args[k];
                bool bound = false;
                if (rule->type == EXPR_FUNCTION &&
                    rule->data.function.head->type == EXPR_SYMBOL &&
                    rule->data.function.head->data.symbol.name == intern_symbol("Rule") &&
                    rule->data.function.arg_count == 2) {
                    Expr* lhs = rule->data.function.args[0];
                    for (size_t i = 0; i < m; i++)
                        if (lhs->type == EXPR_SYMBOL && cs[i]->type == EXPR_SYMBOL &&
                            lhs->data.symbol.name == cs[i]->data.symbol.name) {
                            res[i] = rs_cp(rule->data.function.args[1]);
                            bound = true;
                        }
                }
                if (!bound) unique = false;
            }
        }
        /* Each value must be an explicit ALGEBRAIC number (only numbers, +, x,
         * rational powers, Root): its rationality is then decidable, by
         * RootReduce to the canonical form.  A transcendental value (Pi,
         * Log[2]/Log[4] unsimplified) leaves the question open.  The
         * reduction is also applied to every algebraic value before the
         * rationality test below: (3 + 3 Sqrt[2])/(1 + Sqrt[2]) is 3, and the
         * syntactic test alone read it as irrational -- a false certificate
         * found by the hunt during development. */
        for (size_t i = 0; i < m; i++) {
            if (!res[i]) { unique = false; continue; }
            if (!rs_is_algebraic_number(res[i])) { unique = false; continue; }
            if (!rs_is_rational(res[i])) {
                Expr* rr = rs_call1("RootReduce", rs_cp(res[i]));
                if (rr && rs_is_algebraic_number(rr)) { expr_free(res[i]); res[i] = rr; }
                else { if (rr) expr_free(rr); unique = false; }
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
            /* Certify the relation outside the solver: theta - Sum r_i g_i must
             * cancel to the exact 0.  SolveAlways solves for every non-variable
             * atom too (Log[c] -> 0, Sin[x] -> 0) and the rules for those are
             * dropped above, so its r_i alone are not evidence. */
            Expr** ts = malloc((m + 1) * sizeof(Expr*));
            ts[0] = rs_cp(theta);
            for (size_t i = 0; i < m; i++)
                ts[i + 1] = rs_times(expr_new_integer(-1), rs_times(rs_cp(res[i]), rs_cp(gens[i])));
            Expr* chk = rs_call1("Together",
                expr_new_function(expr_new_symbol("Plus"), ts, m + 1));
            free(ts);
            ok = rs_is_zero(chk);
            expr_free(chk);
        }
        bool nonrat_proof = false;
        if (!ok && unique) {
            /* Certify the exact (non-rational) point before reading it as a
             * proof: theta - Sum v_i g_i must cancel to 0. */
            Expr** ts = malloc((m + 1) * sizeof(Expr*));
            ts[0] = rs_cp(theta);
            for (size_t i = 0; i < m; i++)
                ts[i + 1] = rs_times(expr_new_integer(-1), rs_times(rs_cp(res[i]), rs_cp(gens[i])));
            Expr* chk = rs_call1("Together", rs_call1("Expand",
                expr_new_function(expr_new_symbol("Plus"), ts, m + 1)));
            free(ts);
            nonrat_proof = rs_is_zero(chk);
            expr_free(chk);
        }
        if (ok) {
            result = expr_new_function(expr_new_symbol("List"), res, m);
        } else if (nonrat_proof) {
            /* the unique relation has a non-rational coefficient: False */
            for (size_t i = 0; i < m; i++) expr_free(res[i]);
        } else {
            /* A non-rational particular solution: over a rational coefficient
             * field the affine solution space is defined over Q and its pinned
             * point is rational, so this happens only with algebraic or
             * symbolic constants, where a rational point may still exist.
             * Likewise a solution that fails the certificate above.
             * Undecided, not a proof. */
            for (size_t i = 0; i < m; i++) expr_free(res[i]);
            *failed = true;
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
                                 int target, bool* failed) {
    *failed = true;   /* until the span decision runs, a NULL is not a proof */
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

    Expr* coeffs = risch_rational_span(theta, fg, fm, vars, failed);

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
static Expr* rs_logderiv_radical(const Expr* f, const Expr* x, const Expr* mons,
                                 bool* failed) {
    *failed = true;   /* until the span decision runs, a NULL is not a proof */
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

    Expr* coeffs = risch_rational_span(f, fg, fm, vars, failed);

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
    bool failed;
    Expr* r = risch_rational_span(theta, glist->data.function.args, m, vars, &failed);
    return r ? r : expr_new_symbol(failed ? "$Failed" : "False");
}

/* The four structure-theorem reducibility builtins share one shape. */
static Expr* rs_reducible_builtin(Expr* res, int target) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    bool failed;
    Expr* r = rs_structure_decide(res->data.function.args[0], res->data.function.args[1],
                                  res->data.function.args[2], target, &failed);
    return r ? r : expr_new_symbol(failed ? "$Failed" : "False");
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
    bool failed;
    Expr* r = rs_logderiv_radical(res->data.function.args[0], res->data.function.args[1],
                                  res->data.function.args[2], &failed);
    return r ? r : expr_new_symbol(failed ? "$Failed" : "False");
}

/* ------------------------------------------------------------------ */
/* Registration.                                                       */
/* ------------------------------------------------------------------ */

static void rs_install(const char* name, Expr* (*fn)(Expr*), const char* doc) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    if (doc) symtab_set_docstring(name, doc);
}

void integrate_risch_structure_init(void) {
    rs_install("Risch`RationalSpan", builtin_risch_rationalspan,
        "Risch`RationalSpan[theta, {g1,...}, {vars}] returns {r1,...} of rational\n"
        "numbers with theta = Sum r_i g_i (coefficients matched over vars), or\n"
        "False if theta is not in the rational span of the generators, or\n"
        " $Failed when the decision could not be made (a malformed tower, or the linear\n"
        "solve gave up); False is returned only as a proof.");
    rs_install("Risch`LogReducible", builtin_risch_logreducible,
        "Risch`LogReducible[a, x, {{t, \"Exp\"|\"Log\", Dt}, ...}] applies the Risch\n"
        "structure theorem: returns the rational coefficients expressing Da/a in the\n"
        "monomial generators when Log[a] is reducible over the tower, else False\n"
        "(Log[a] is a new, algebraically independent monomial), or\n"
        " $Failed when the decision could not be made (a malformed tower, or the linear\n"
        "solve gave up); False is returned only as a proof.");
    rs_install("Risch`ExpReducible", builtin_risch_expreducible,
        "Risch`ExpReducible[b, x, {{t, \"Exp\"|\"Log\", Dt}, ...}] applies the Risch\n"
        "structure theorem: returns the rational coefficients expressing Db in the\n"
        "monomial generators when Exp[b] is reducible over the tower, else False\n"
        "(Exp[b] is a new, algebraically independent monomial), or\n"
        " $Failed when the decision could not be made (a malformed tower, or the linear\n"
        "solve gave up); False is returned only as a proof.");
    rs_install("Risch`TanReducible", builtin_risch_tanreducible,
        "Risch`TanReducible[b, x, {{t, \"Tan\"|\"ArcTan\"|\"Exp\"|\"Log\", Dt}, ...}]\n"
        "applies the REAL Risch structure theorem (Cor. 9.3.2 iv, eq. 9.15): returns\n"
        "the rational coefficients over the tangent/arc-tangent generators when Tan[b]\n"
        "is reducible over the tower, else False. Only T u A monomials participate.\n"
        " $Failed when the decision could not be made (a malformed tower, or the linear\n"
        "solve gave up); False is returned only as a proof.");
    rs_install("Risch`ArcTanReducible", builtin_risch_arctanreducible,
        "Risch`ArcTanReducible[b, x, {{t, \"Tan\"|\"ArcTan\"|\"Exp\"|\"Log\", Dt}, ...}]\n"
        "applies the REAL Risch structure theorem (Cor. 9.3.2 iii, eq. 9.14): returns\n"
        "the rational coefficients expressing Db/(b^2+1) over the tangent/arc-tangent\n"
        "generators when ArcTan[b] is reducible over the tower, else False, or\n"
        " $Failed when the decision could not be made (a malformed tower, or the linear\n"
        "solve gave up); False is returned only as a proof.");
    rs_install("Risch`LogarithmicDerivativeOfRadical", builtin_risch_logderiv_radical,
        "Risch`LogarithmicDerivativeOfRadical[f, x, {{t, kind, Dt[, arg]}, ...}]\n"
        "decides whether f is the logarithmic derivative of a K-radical over the tower\n"
        "(Bronstein §5.12, structure-theorem test eq. 7.44/9.13): returns {n, u} with\n"
        "n a positive integer and u the witness radical such that D[u]/u == n f, else\n"
        "False. Only the Exp/Log (E u L) monomials participate; supply the log argument\n"
        "as the optional 4th monomial element for the witness u.\n"
        " $Failed when the decision could not be made (a malformed tower, or the linear\n"
        "solve gave up); False is returned only as a proof.");
}
