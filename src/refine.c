/* Mathilda -- Refine implementation.
 *
 * Refine[expr, assum] gives the form of `expr` that would be obtained if the
 * symbols in it were replaced by explicit values satisfying `assum`. It is a
 * thin orchestrator over machinery that already exists and was deliberately
 * shared for it (simp.h reserves AssumeCtx "so future modules (Refine, ...)
 * can share it"):
 *
 *   1. Argument / option parsing mirrors builtin_possible_zero_q
 *      (src/zero_test.c): one positional assumption, an Assumptions -> X option
 *      that overrides $Assumptions, and a TimeConstraint option. The effective
 *      assumption is And-combined exactly as PossibleZeroQ/Simplify do.
 *
 *   2. Predicate expressions are decided to True/False:
 *        Element[x, dom]              -> element_decide (assumption-aware,
 *                                        including compound-expression domain
 *                                        inference via prov_re/prov_int).
 *        Equal / Unequal              -> zero_test_decide_assuming on lhs - rhs,
 *                                        with a Reduce fallback.
 *        Less / ... / logic           -> Reduce/CAD entailment over the reals:
 *                                        P is True  iff Reduce[A && !P] is False,
 *                                        P is False iff Reduce[A &&  P] is False.
 *
 *   3. Every other expression is rewritten by apply_assumption_rules
 *      (Sqrt[x^2]->x/-x/Abs[x], Log[x^p]->p Log[x], integer-k trig, Abs/Sign/
 *      Conjugate under sign facts, Floor/Ceiling/Mod under integer/interval
 *      facts, ...) followed by a deep-positivity post-pass that resolves
 *      Sign[p]/Abs[p]/Sqrt[p^2] for symbolic polynomials p that the fast
 *      sign prover cannot settle, using the same Reduce entailment.
 *
 * TimeConstraint (default 30s) is enforced with Simplify's cooperative
 * wall-clock budget (simp_set_time_budget / simp_mono_seconds): Reduce/CAD has
 * no interior abort hook, so the deadline is checked BEFORE each entailment
 * call and a running Reduce is never preempted (best-effort, and never the
 * malloc-lock-prone async TimeConstrained[]).
 *
 * Refine is a purely symbolic/structural head: it has no element-wise numeric
 * semantics, so it deliberately carries no NDArray/packed kernel and no
 * Compile[] lowering.
 *
 * Ownership follows the builtin contract: return a new tree or steal from
 * `res`; never expr_free(res).
 */

#include "refine.h"

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "sym_names.h"

#include "simp.h"           /* AssumeCtx, assume_ctx_*, assume_known_*, read_dollar_assumptions, apply_assumption_rules */
#include "simp_internal.h"  /* is_rule_with_lhs, element_decide, prov_*, is_real_constant_symbol, simp time budget */
#include "numeric.h"          /* expr_is_numeric_quantity */
#include "zero_test.h"      /* zero_test_decide_assuming, ZeroTestResult */
#include "reduce.h"         /* builtin_reduce (via evaluate of Reduce[...]) */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "message.h"   /* mth_message: Quiet/Check funnel */

/* Default TimeConstraint (seconds) for a single condition check, matching
 * Mathematica's Refine default. */
#define REFINE_DEFAULT_TIMECONSTRAINT 30.0

/* Cap on the number of Reduce entailment calls per Refine, so the deep
 * positivity walk over a large expression cannot fan out into an unbounded
 * number of CAD runs. */
#define REFINE_MAX_ENTAILMENT_CALLS 24

/* ------------------------------------------------------------------------- */
/* Small helpers                                                             */
/* ------------------------------------------------------------------------- */

static int head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head && e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

static const char* head_name(const Expr* e) {
    return (e && e->type == EXPR_FUNCTION && e->data.function.head &&
            e->data.function.head->type == EXPR_SYMBOL)
           ? e->data.function.head->data.symbol.name : NULL;
}

static Expr* bool_sym(int v) { return expr_new_symbol(v ? SYM_True : SYM_False); }

static int is_true_sym(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True;
}
static int is_false_sym(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_False;
}

/* A relational / logical head that Refine can try to decide via Reduce. */
static int is_relational_head(const char* h) {
    if (!h) return 0;
    static const char* heads[] = {
        "Equal", "Unequal", "Less", "LessEqual", "Greater", "GreaterEqual",
        "Inequality", "And", "Or", "Not", "Implies", "Xor", NULL };
    for (int i = 0; heads[i]; i++) if (strcmp(h, heads[i]) == 0) return 1;
    return 0;
}

/* True if `name` is a registered option name of Refine (Assumptions,
 * TimeConstraint). Consulted so that a positional assumption written as a
 * rule (e.g. Refine[expr, a -> b]) is not mistaken for a trailing option. */
static int is_refine_option_name(const char* name) {
    Expr* opts = symtab_get_options("Refine");   /* borrowed List of rules */
    if (!opts || opts->type != EXPR_FUNCTION) return 0;
    for (size_t i = 0; i < opts->data.function.arg_count; i++) {
        const Expr* r = opts->data.function.args[i];
        if ((head_is(r, "Rule") || head_is(r, "RuleDelayed")) &&
            r->data.function.arg_count == 2) {
            const Expr* lhs = r->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL &&
                strcmp(lhs->data.symbol.name, name) == 0) return 1;
        }
    }
    return 0;
}

/* A Rule/RuleDelayed whose LHS symbol names one of Refine's options. */
static int is_option_rule(const Expr* e) {
    if ((head_is(e, "Rule") || head_is(e, "RuleDelayed")) &&
        e->data.function.arg_count == 2) {
        const Expr* lhs = e->data.function.args[0];
        return lhs->type == EXPR_SYMBOL &&
               is_refine_option_name(lhs->data.symbol.name);
    }
    return 0;
}

/* And-combine two owned assumption expressions and evaluate (canonicalises
 * And[True, x] -> x, flattens nested And). Both inputs are consumed. */
static Expr* combine_and(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    Expr* call = expr_new_function(expr_new_symbol(SYM_And), args, 2);
    Expr* out = evaluate(call);
    expr_free(call);
    return out;
}

/* Parse a TimeConstraint value into a wall-clock budget in seconds. Infinity /
 * DirectedInfinity[1] / non-positive / unparsable -> HUGE_VAL ("no limit"). A
 * NULL option means "use the Refine default" (30s). */
static double parse_time_budget(const Expr* e) {
    if (!e) return REFINE_DEFAULT_TIMECONSTRAINT;
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity)
        return HUGE_VAL;
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_DirectedInfinity)
        return HUGE_VAL;
    /* A machine-real positive value is the per-condition budget. */
    if (e->type == EXPR_INTEGER && e->data.integer > 0)
        return (double)e->data.integer;
    if (e->type == EXPR_REAL && e->data.real > 0.0)
        return e->data.real;
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational &&
        e->data.function.arg_count == 2) {
        const Expr* p = e->data.function.args[0];
        const Expr* q = e->data.function.args[1];
        if (p->type == EXPR_INTEGER && q->type == EXPR_INTEGER && q->data.integer != 0) {
            double v = (double)p->data.integer / (double)q->data.integer;
            if (v > 0.0) return v;
        }
    }
    return REFINE_DEFAULT_TIMECONSTRAINT;
}

/* ------------------------------------------------------------------------- */
/* Reduce/CAD entailment                                                     */
/* ------------------------------------------------------------------------- */

/* State threaded through the entailment machinery: the wall-clock deadline
 * (HUGE_VAL for none) and a decrementing call budget. */
typedef struct {
    double deadline;   /* simp_mono_seconds() cutoff; HUGE_VAL means unlimited */
    int    calls_left; /* remaining Reduce invocations */
} RefineBudget;

static int budget_ok(RefineBudget* b) {
    if (b->calls_left <= 0) return 0;
    if (b->deadline < HUGE_VAL && simp_mono_seconds() > b->deadline) return 0;
    return 1;
}

/* Recursively collect distinct bare-symbol leaves of `e` that are not
 * protected real constants (Pi, E, ...). Appends owned copies to *vars. Do NOT
 * reuse the polynomial collect_variables, which atomizes x^p / Floor[x] into
 * whole "variables" that Reduce rejects. */
static void collect_bare_vars(const Expr* e, Expr*** vars, size_t* n, size_t* cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        if (is_real_constant_symbol(nm)) return;
        for (size_t i = 0; i < *n; i++) {
            if ((*vars)[i]->type == EXPR_SYMBOL &&
                strcmp((*vars)[i]->data.symbol.name, nm) == 0) return; /* dup */
        }
        if (*n == *cap) {
            *cap = *cap ? *cap * 2 : 8;
            *vars = (Expr**)realloc(*vars, *cap * sizeof(Expr*));
        }
        (*vars)[(*n)++] = expr_copy((Expr*)e);
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        /* The head of a call is a structural symbol (Plus, Sin, ...), never a
         * solve variable; descend into arguments only. */
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            collect_bare_vars(e->data.function.args[i], vars, n, cap);
    }
}

/* Decide whether `stmt` is unsatisfiable over `domain` (a domain symbol name
 * such as "Reals" or "Complexes"). Returns 1 iff Reduce proves it False
 * (unsatisfiable); 0 otherwise (satisfiable, a residual formula, out of reach,
 * or no usable variables). `stmt` is consumed. */
static int reduce_is_unsat(Expr* stmt, const char* domain, RefineBudget* b);

/* A conjunction whose conjuncts fall into groups sharing no variable is
 * unsatisfiable iff one of the groups is.  Deciding the groups apart keeps an
 * irrelevant fact out of each Reduce call: Refine[Pi/a > 0, a > 0 && s > 0]
 * asked Reduce[a > 0 && s > 0 && Pi/a <= 0, {a, s}, Reals], which the
 * multivariate CAD leaves unevaluated (a transcendental coefficient), although
 * the group {a > 0, Pi/a <= 0} alone is decided False.  Returns 1 (some group
 * unsatisfiable), 0 (none proved), or -1 when the statement is not split
 * (a single group); consumes stmt only when it returns 0 or 1. */
static int reduce_unsat_by_components(Expr* stmt, const char* domain, RefineBudget* b) {
    if (!stmt || stmt->type != EXPR_FUNCTION || head_name(stmt) != SYM_And) return -1;
    size_t n = stmt->data.function.arg_count;
    if (n < 2 || n > 64) return -1;
    Expr*** cv = calloc(n, sizeof(Expr**));
    size_t* cn = calloc(n, sizeof(size_t));
    int* comp = malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) {
        size_t cap = 0;
        collect_bare_vars(stmt->data.function.args[i], &cv[i], &cn[i], &cap);
        comp[i] = (int)i;
    }
    /* union conjuncts that share a variable (union-find, path halving) */
    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++) {
            bool share = false;
            for (size_t p = 0; p < cn[i] && !share; p++)
                for (size_t q = 0; q < cn[j] && !share; q++)
                    if (expr_eq(cv[i][p], cv[j][q])) share = true;
            if (!share) continue;
            int ri = (int)i, rj = (int)j;
            while (comp[ri] != ri) ri = comp[ri] = comp[comp[ri]];
            while (comp[rj] != rj) rj = comp[rj] = comp[comp[rj]];
            if (ri != rj) comp[rj] = ri;
        }
    int nroots = 0;
    for (size_t i = 0; i < n; i++) {
        int r = (int)i;
        while (comp[r] != r) r = comp[r];
        comp[i] = r;
        if (r == (int)i) nroots++;
    }
    int result = -1;
    if (nroots >= 2) {
        result = 0;
        for (size_t r = 0; r < n && result == 0; r++) {
            if (comp[r] != (int)r) continue;
            size_t m = 0;
            for (size_t i = 0; i < n; i++) if (comp[i] == (int)r) m++;
            Expr** part = malloc(m * sizeof(Expr*));
            size_t k = 0;
            for (size_t i = 0; i < n; i++)
                if (comp[i] == (int)r) part[k++] = expr_copy(stmt->data.function.args[i]);
            Expr* g = (m == 1) ? part[0]
                               : expr_new_function(expr_new_symbol(SYM_And), part, m);
            free(part);
            if (reduce_is_unsat(g, domain, b) == 1) result = 1;
        }
        expr_free(stmt);
    }
    for (size_t i = 0; i < n; i++) {
        for (size_t p = 0; p < cn[i]; p++) expr_free(cv[i][p]);
        free(cv[i]);
    }
    free(cv); free(cn); free(comp);
    return result;
}

static int reduce_is_unsat(Expr* stmt, const char* domain, RefineBudget* b) {
    if (!budget_ok(b)) { expr_free(stmt); return 0; }

    /* A statement that already evaluated to a literal decides directly. */
    if (is_false_sym(stmt)) { expr_free(stmt); return 1; }
    if (is_true_sym(stmt))  { expr_free(stmt); return 0; }

    {
        int c = reduce_unsat_by_components(stmt, domain, b);
        if (c >= 0) return c;
    }

    Expr** vars = NULL; size_t nv = 0, cap = 0;
    collect_bare_vars(stmt, &vars, &nv, &cap);
    if (nv == 0 || nv > 6) {           /* nothing to solve for, or too wide for CAD */
        for (size_t i = 0; i < nv; i++) expr_free(vars[i]);
        free(vars);
        expr_free(stmt);
        return 0;
    }

    Expr* varlist = expr_new_function(expr_new_symbol(SYM_List), vars, nv);
    free(vars);

    Expr* rargs[3] = { stmt, varlist, expr_new_symbol(domain) };
    Expr* call = expr_new_function(expr_new_symbol("Reduce"), rargs, 3);
    b->calls_left--;
    /* An internal probe: whatever Reduce or the Solve under it reports
     * (Solve::svars on a > b && a == b over the complexes) stays silent. */
    mth_msg_suppress_push();
    Expr* out = evaluate(call);
    mth_msg_suppress_pop();
    expr_free(call);

    int unsat = is_false_sym(out);     /* only a literal False proves unsatisfiability */
    if (out) expr_free(out);
    return unsat;
}

/* Build the conjunction of the relational/logical facts in `ctx` (the
 * inequality/equation assumptions), as an owned expression, or the symbol True
 * when there are none. Domain-declaration facts (Element[...]) are omitted: the
 * Reals domain already makes the variables real, which is exactly the
 * "quantities appearing algebraically in inequalities are assumed real" rule. */
static Expr* relational_assumption_conj(const AssumeCtx* ctx) {
    if (!ctx || ctx->count == 0) return expr_new_symbol(SYM_True);
    Expr** kept = (Expr**)calloc(ctx->count, sizeof(Expr*));
    size_t nk = 0;
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (f->type == EXPR_FUNCTION && is_relational_head(head_name(f)))
            kept[nk++] = expr_copy((Expr*)f);
    }
    Expr* out;
    if (nk == 0)      out = expr_new_symbol(SYM_True);
    else if (nk == 1) { out = kept[0]; }
    else              out = expr_new_function(expr_new_symbol(SYM_And), kept, nk);
    free(kept);
    return out;
}

/* Given the relational assumption conjunction `A` (borrowed) and a predicate
 * `P` (borrowed), decide P over `domain`:
 *   returns  1  if A entails P     (Reduce[A && !P] is unsatisfiable)
 *   returns  0  if A entails not-P (Reduce[A &&  P] is unsatisfiable)
 *   returns -1  if undecided. */
static int reduce_entail_predicate(const Expr* A, const Expr* P,
                                   const char* domain, RefineBudget* b) {
    /* A && !P unsatisfiable  =>  P holds everywhere A holds. */
    Expr* notP = eval_and_free(expr_new_function(expr_new_symbol(SYM_Not),
                              (Expr*[]){ expr_copy((Expr*)P) }, 1));
    Expr* stmt_true = combine_and(expr_copy((Expr*)A), notP);
    if (reduce_is_unsat(stmt_true, domain, b)) return 1;

    /* A && P unsatisfiable  =>  P fails everywhere A holds. */
    Expr* stmt_false = combine_and(expr_copy((Expr*)A), expr_copy((Expr*)P));
    if (reduce_is_unsat(stmt_false, domain, b)) return 0;

    return -1;
}

/* ------------------------------------------------------------------------- */
/* Predicate folding (the top-level expression is itself a statement)        */
/* ------------------------------------------------------------------------- */

/* Try to decide the predicate expression `P` to True/False. Returns an owned
 * True/False symbol, or NULL when it cannot be decided (caller falls back to
 * the rewrite path). */
static Expr* refine_decide_predicate(const Expr* P, const AssumeCtx* ctx,
                                     const Expr* Aconj, RefineBudget* b) {
    const char* h = head_name(P);
    if (!h) return NULL;

    /* Element[x, dom]: assumption-aware domain membership (compound x handled
     * inside element_decide via prov_re/prov_int). */
    if (h == SYM_Element && P->data.function.arg_count == 2) {
        const Expr* x   = P->data.function.args[0];
        const Expr* dom = P->data.function.args[1];
        if (dom->type == EXPR_SYMBOL) {
            int d = element_decide(x, dom->data.symbol.name, ctx);
            if (d == 1) return bool_sym(1);
            if (d == 0) return bool_sym(0);
        }
        return NULL;
    }

    /* Equal / Unequal: decide via the assumption-aware zero test on the
     * difference, then fall back to Reduce over the complexes (an equation is
     * not restricted to the reals). */
    if ((h == SYM_Equal || h == SYM_Unequal) &&
        P->data.function.arg_count == 2) {
        int is_eq = (h == SYM_Equal);
        Expr* lhs = P->data.function.args[0];
        Expr* rhs = P->data.function.args[1];
        Expr* diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(lhs),
                       expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ expr_new_integer(-1), expr_copy(rhs) }, 2) }, 2));
        ZeroTestResult z = zero_test_decide_assuming(diff, ctx);
        if (z == ZERO_TEST_TRUE)  { expr_free(diff); return bool_sym(is_eq ? 1 : 0); }
        /* A non-zero verdict says the difference is not identically zero on
         * the assumed region, which decides the equation only when the
         * difference is a constant: m == a is not False, merely not an
         * identity (the reference leaves it alone).  A symbolic difference
         * goes on to the Reduce entailment, which does decide x == -1 under
         * x > 0. */
        if (z == ZERO_TEST_FALSE && expr_is_numeric_quantity(diff)) {
            expr_free(diff);
            return bool_sym(is_eq ? 0 : 1);
        }
        /* Undecided by the zero test: try to prove the difference is exactly
         * zero by rewriting it under the (equality) assumptions -- the same
         * substitution engine that turns Refine[a - b, a == b] into 0. Every
         * rewrite is value-preserving under the assumptions, so if the reduced
         * difference is unconditionally zero, the equation holds. */
        Expr* red = apply_assumption_rules(diff, ctx);
        if (red) {
            Expr* redev = eval_and_free(red);
            ZeroTestResult zr = zero_test_decide(redev);
            expr_free(redev);
            if (zr == ZERO_TEST_TRUE) { expr_free(diff); return bool_sym(is_eq ? 1 : 0); }
        }
        expr_free(diff);
        int r = reduce_entail_predicate(Aconj, P, "Complexes", b);
        if (r == 1) return bool_sym(1);
        if (r == 0) return bool_sym(0);
        /* Both sides real under the assumptions (an inequality assumption
         * makes its variables real: a > b, x > 0 && y < 0): the equation lives
         * on the reals, where Reduce decides a == b and x^2 == -1. */
        if (prov_re(ctx, lhs) && prov_re(ctx, rhs)) {
            r = reduce_entail_predicate(Aconj, P, "Reals", b);
            if (r == 1) return bool_sym(1);
            if (r == 0) return bool_sym(0);
        }
        return NULL;
    }

    /* Inequalities and logical combinations: Reduce/CAD entailment over the
     * reals (algebraic quantities in inequalities are assumed real). */
    if (is_relational_head(h)) {
        int r = reduce_entail_predicate(Aconj, P, "Reals", b);
        if (r == 1) return bool_sym(1);
        if (r == 0) return bool_sym(0);
    }
    return NULL;
}

/* ------------------------------------------------------------------------- */
/* Deep-positivity post-pass                                                 */
/* ------------------------------------------------------------------------- */

/* If p^2 written as Power[p, 2] (the canonical Sqrt-of-square inner form). */
static const Expr* square_base(const Expr* e) {
    if (head_is(e, "Power") && e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* x = e->data.function.args[1];
        if (x->type == EXPR_INTEGER && x->data.integer == 2) return b;
    }
    return NULL;
}

/* Power[p, Rational[1,2]] -> its radicand p, else NULL. */
static const Expr* sqrt_radicand(const Expr* e) {
    if (head_is(e, "Power") && e->data.function.arg_count == 2) {
        const Expr* x = e->data.function.args[1];
        if (x->type == EXPR_FUNCTION && x->data.function.head &&
            x->data.function.head->type == EXPR_SYMBOL &&
            x->data.function.head->data.symbol.name == SYM_Rational &&
            x->data.function.arg_count == 2) {
            const Expr* p = x->data.function.args[0];
            const Expr* q = x->data.function.args[1];
            if (p->type == EXPR_INTEGER && p->data.integer == 1 &&
                q->type == EXPR_INTEGER && q->data.integer == 2)
                return e->data.function.args[0];
        }
    }
    return NULL;
}

/* Prove the sign of a symbolic expression `p` via Reduce over the reals:
 *   +1 if p > 0 everywhere A holds, -1 if p < 0, 0 if undecided. */
static int deep_sign(const Expr* p, const AssumeCtx* ctx, const Expr* Aconj,
                     RefineBudget* b) {
    /* Reduce over the Reals reads every subexpression as real-valued, so for a
     * non-real p BOTH p > 0 and p <= 0 are unsatisfiable and the entailment
     * below "proves" p > 0 vacuously: Arg[(-a - Sqrt[-3 a^2])/2] -> 0 and
     * Abs[-a + I a] -> -a + I a under a > 0.  The sign is only meaningful once
     * p is known to be real. */
    if (!prov_re(ctx, p)) return 0;
    Expr* pos = expr_new_function(expr_new_symbol(SYM_Greater),
                    (Expr*[]){ expr_copy((Expr*)p), expr_new_integer(0) }, 2);
    int r = reduce_entail_predicate(Aconj, pos, "Reals", b);
    expr_free(pos);
    if (r == 1) return 1;
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Less),
                    (Expr*[]){ expr_copy((Expr*)p), expr_new_integer(0) }, 2);
    r = reduce_entail_predicate(Aconj, neg, "Reals", b);
    expr_free(neg);
    if (r == 1) return -1;
    return 0;
}

/* Extract perfect-square factors from a Sqrt radicand under proved signs:
 *   Sqrt[Product base_i^(2 k_i) * R] -> (Product |base_i|^k_i) * Sqrt[R],
 * each |base_i| resolved by deep_sign.  A Plus radicand is Together'd first so a
 * denominator is exposed (1 + b^2/a^2 -> (a^2+b^2)/a^2 -> Sqrt[a^2+b^2]/a).
 * Returns an owned rewritten tree, or NULL when nothing is extractable.  Sound:
 * every pulled |base|^k is non-negative, matching Sqrt's principal branch. */
static Expr* sqrt_extract(const Expr* rad, const AssumeCtx* ctx,
                          const Expr* Aconj, RefineBudget* b) {
    Expr* together = NULL;
    const Expr* work = rad;
    if (head_is(rad, "Plus")) {
        together = eval_and_free(expr_new_function(expr_new_symbol("Together"),
                      (Expr*[]){ expr_copy((Expr*)rad) }, 1));
        if (together) work = together;
    }
    const Expr* one[1];
    Expr* const* fs;
    size_t nf;
    if (head_is(work, "Times")) {
        fs = work->data.function.args; nf = work->data.function.arg_count;
    } else { one[0] = work; fs = (Expr* const*)one; nf = 1; }

    Expr** pulled = (Expr**)calloc(nf ? nf : 1, sizeof(Expr*)); size_t np = 0;
    Expr** rem    = (Expr**)calloc(nf ? nf : 1, sizeof(Expr*)); size_t nr = 0;
    int any = 0;
    for (size_t i = 0; i < nf; i++) {
        const Expr* f = fs[i];
        if (head_is(f, "Power") && f->data.function.arg_count == 2 &&
            f->data.function.args[1]->type == EXPR_INTEGER &&
            (f->data.function.args[1]->data.integer % 2) == 0 &&
            budget_ok(b)) {
            const Expr* base = f->data.function.args[0];
            long k = f->data.function.args[1]->data.integer / 2;
            int s = deep_sign(base, ctx, Aconj, b);
            if (s != 0) {
                Expr* absb = (s > 0)
                    ? expr_copy((Expr*)base)
                    : expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)base) }, 2);
                pulled[np++] = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ absb, expr_new_integer(k) }, 2);
                any = 1;
                continue;
            }
        }
        rem[nr++] = expr_copy((Expr*)f);
    }
    if (together) expr_free(together);
    if (!any) {
        for (size_t i = 0; i < np; i++) expr_free(pulled[i]);
        for (size_t i = 0; i < nr; i++) expr_free(rem[i]);
        free(pulled); free(rem);
        return NULL;
    }
    Expr* pulled_prod = (np == 1) ? pulled[0]
        : expr_new_function(expr_new_symbol(SYM_Times), pulled, np);
    Expr* out;
    if (nr == 0) {
        out = pulled_prod;                       /* Sqrt[1] = 1 */
    } else {
        Expr* rem_prod = (nr == 1) ? rem[0]
            : expr_new_function(expr_new_symbol(SYM_Times), rem, nr);
        Expr* half = expr_new_function(expr_new_symbol(SYM_Rational),
              (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2);
        Expr* sq = expr_new_function(expr_new_symbol(SYM_Power),
              (Expr*[]){ rem_prod, half }, 2);
        out = expr_new_function(expr_new_symbol(SYM_Times),
              (Expr*[]){ pulled_prod, sq }, 2);
    }
    free(pulled); free(rem);
    return eval_and_free(out);                   /* tidy Times[1,..], Power[x,1] */
}

/* Bottom-up walk resolving Sign[p], Abs[p] and Sqrt[p^2] for polynomial p whose
 * sign the fast prover could not settle. Returns an owned rebuilt tree; sets
 * *changed when any node was rewritten. */
static Expr* deep_positivity_walk(const Expr* e, const AssumeCtx* ctx,
                                  const Expr* Aconj, RefineBudget* b,
                                  int* changed) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Rebuild children first. */
    size_t n = e->data.function.arg_count;
    Expr** args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++)
        args[i] = deep_positivity_walk(e->data.function.args[i], ctx, Aconj, b, changed);
    Expr* rebuilt = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);

    if (!budget_ok(b)) return rebuilt;

    const char* h = head_name(rebuilt);
    if (!h) return rebuilt;

    /* Sign[p] -> +/-1 ; Abs[p] -> +/-p. */
    if ((h == SYM_Sign || h == SYM_Abs) &&
        rebuilt->data.function.arg_count == 1) {
        const Expr* p = rebuilt->data.function.args[0];
        if (p->type == EXPR_FUNCTION) {          /* only worth a CAD run for compound p */
            int s = deep_sign(p, ctx, Aconj, b);
            if (s != 0) {
                *changed = 1;
                Expr* out;
                if (h == SYM_Sign) {
                    out = expr_new_integer(s);
                } else if (s > 0) {
                    out = expr_copy((Expr*)p);
                } else {
                    out = expr_new_function(expr_new_symbol(SYM_Times),
                             (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)p) }, 2);
                }
                expr_free(rebuilt);
                return out;
            }
        }
    }

    /* Arg[p] -> 0 (p > 0) / Pi (p < 0) for a provably-signed real compound p.
     * Sound because deep_sign only settles the sign of a REAL expression (the
     * same basis as the Abs/Sqrt collapses above); a complex argument leaves
     * deep_sign == 0 and Arg untouched.  Collapses the keyhole residue forms
     * (e.g. Arg[-a] with a > 0 -> Pi, which closes Log[x]/(x^2-a^2) PV). */
    if (h == SYM_Arg && rebuilt->data.function.arg_count == 1) {
        const Expr* p = rebuilt->data.function.args[0];
        if (p->type == EXPR_FUNCTION) {          /* only worth a CAD run for compound p */
            int s = deep_sign(p, ctx, Aconj, b);
            if (s != 0) {
                *changed = 1;
                Expr* out = (s > 0) ? expr_new_integer(0) : expr_new_symbol(SYM_Pi);
                expr_free(rebuilt);
                return out;
            }
        }
    }

    /* Sqrt[p^2] -> +/-p when the sign of p is known. */
    const Expr* rad = sqrt_radicand(rebuilt);
    if (rad) {
        const Expr* base = square_base(rad);
        if (base && base->type == EXPR_FUNCTION) {
            int s = deep_sign(base, ctx, Aconj, b);
            if (s != 0) {
                *changed = 1;
                Expr* out = (s > 0)
                    ? expr_copy((Expr*)base)
                    : expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)base) }, 2);
                expr_free(rebuilt);
                return out;
            }
        }
    }
    /* General radicand under a `^(1/2)` or `^(-1/2)` power: a product of squared
     * factors, or a quotient whose denominator is a perfect square
     * (Sqrt[(a^2+b^2)/a^2] -> Sqrt[a^2+b^2]/a, and 1/Sqrt[...] = ^(-1/2) too, so
     * 1/(a Sqrt[1+b^2/a^2]) -> 1/Sqrt[a^2+b^2]).  Only for a Times/Plus radicand,
     * so the single-square path above is untouched. */
    {
        int half = 0;                       /* +1 for ^(1/2), -1 for ^(-1/2) */
        const Expr* grad = NULL;
        if (head_is(rebuilt, "Power") && rebuilt->data.function.arg_count == 2) {
            const Expr* ex = rebuilt->data.function.args[1];
            if (ex->type == EXPR_FUNCTION && head_is(ex, "Rational") &&
                ex->data.function.arg_count == 2 &&
                ex->data.function.args[1]->type == EXPR_INTEGER &&
                ex->data.function.args[1]->data.integer == 2 &&
                ex->data.function.args[0]->type == EXPR_INTEGER) {
                long pnum = ex->data.function.args[0]->data.integer;
                if (pnum == 1)  { half = 1;  grad = rebuilt->data.function.args[0]; }
                if (pnum == -1) { half = -1; grad = rebuilt->data.function.args[0]; }
            }
        }
        if (grad && (head_is(grad, "Times") || head_is(grad, "Plus"))) {
            Expr* ext = sqrt_extract(grad, ctx, Aconj, b);
            if (ext) {
                if (half < 0)
                    ext = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ ext, expr_new_integer(-1) }, 2));
                *changed = 1; expr_free(rebuilt); return ext;
            }
        }
    }

    return rebuilt;
}

/* The unit direction a DirectedInfinity[d] collapses to under the facts, as
 * an owned +1 / -1 / I / -I literal, or NULL when it cannot be decided.
 * Refines Sign[d] with the same two passes the rewrite path runs (the fast
 * rule set, then the Reduce-backed deep positivity walk), so the direction
 * is settled exactly as Refine[Sign[d], assum] would settle it. */
static Expr* refine_unit_direction(const Expr* d, const AssumeCtx* ctx,
                                   const Expr* Aconj, RefineBudget* b) {
    Expr* sgn = eval_and_free(expr_new_function(expr_new_symbol(SYM_Sign),
                    (Expr*[]){ expr_copy((Expr*)d) }, 1));
    Expr* r = apply_assumption_rules(sgn, ctx);
    if (r) { expr_free(sgn); sgn = r; }
    int changed = 0;
    Expr* deep = deep_positivity_walk(sgn, ctx, Aconj, b, &changed);
    expr_free(sgn);
    deep = eval_and_free(deep);
    if (deep->type == EXPR_INTEGER &&
        (deep->data.integer == 1 || deep->data.integer == -1))
        return deep;
    if (head_is(deep, "Complex") && deep->data.function.arg_count == 2) {
        const Expr* re = deep->data.function.args[0];
        const Expr* im = deep->data.function.args[1];
        if (re->type == EXPR_INTEGER && re->data.integer == 0 &&
            im->type == EXPR_INTEGER &&
            (im->data.integer == 1 || im->data.integer == -1))
            return deep;
    }
    expr_free(deep);
    return NULL;
}

/* Bottom-up walk canonicalising every DirectedInfinity[d] whose direction is
 * decided by the facts: DirectedInfinity[1/(s - 1)] under s > 1 is Infinity,
 * DirectedInfinity[1 - s] is -Infinity. An undecidable direction is left as
 * it is. A literal direction is already canonical (the evaluator normalises
 * it), so only a symbolic one is examined. Returns an owned tree. */
static Expr* refine_directed_infinity(const Expr* e, const AssumeCtx* ctx,
                                      const Expr* Aconj, RefineBudget* b) {
    if (!e || e->type != EXPR_FUNCTION) return e ? expr_copy((Expr*)e) : NULL;
    size_t n = e->data.function.arg_count;
    if (head_is(e, "DirectedInfinity") && n == 1) {
        const Expr* d = e->data.function.args[0];
        if (d->type == EXPR_SYMBOL || d->type == EXPR_FUNCTION) {
            Expr* u = budget_ok(b) ? refine_unit_direction(d, ctx, Aconj, b) : NULL;
            if (u) {
                /* +-1 take the canonical Infinity / -Infinity spelling. */
                if (u->type == EXPR_INTEGER) {
                    int64_t sg = u->data.integer;
                    expr_free(u);
                    Expr* inf = expr_new_symbol(SYM_Infinity);
                    return sg > 0 ? inf
                        : expr_new_function(expr_new_symbol(SYM_Times),
                              (Expr*[]){ expr_new_integer(-1), inf }, 2);
                }
                return expr_new_function(expr_new_symbol(SYM_DirectedInfinity),
                                         (Expr*[]){ u }, 1);
            }
        }
        return expr_copy((Expr*)e);
    }
    Expr** args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++)
        args[i] = refine_directed_infinity(e->data.function.args[i], ctx, Aconj, b);
    Expr* out = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return out;
}

/* ------------------------------------------------------------------------- */
/* builtin_refine                                                            */
/* ------------------------------------------------------------------------- */

Expr* builtin_refine(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    /* Split positional args from trailing options. */
    Expr* opt_assumptions   = NULL;
    Expr* opt_timeconstraint = NULL;
    Expr* positional[2] = { NULL, NULL };
    size_t npos = 0;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (is_rule_with_lhs(a, "Assumptions")) {
            opt_assumptions = a->data.function.args[1];
        } else if (is_rule_with_lhs(a, "TimeConstraint")) {
            opt_timeconstraint = a->data.function.args[1];
        } else if (is_option_rule(a)) {
            /* a registered option we do not consume here -- ignore */
        } else {
            if (npos < 2) positional[npos] = a;
            npos++;
        }
    }

    if (npos == 0) {
        mth_message("Refine", "argt", "Refine called with 0 arguments; 1 or 2 arguments are expected.");
        return NULL;
    }
    if (npos > 2) {
        mth_message("Refine", "argx", "Refine called with %zu arguments; 1 or 2 arguments are expected.", npos);
        return NULL;
    }

    Expr* expr = positional[0];
    Expr* positional_assum = (npos == 2) ? positional[1] : NULL;

    /* Effective assumption (identical policy to PossibleZeroQ / Simplify):
     *   Assumptions->X with a positional A  -> A && X   ($Assumptions unused)
     *   Assumptions->X only                 -> X        ($Assumptions unused)
     *   positional A only                   -> A && $Assumptions
     *   neither                             -> $Assumptions
     * Then evaluate to canonicalise (And[True, x>0] -> x>0). */
    Expr* effective;
    if (opt_assumptions) {
        if (positional_assum)
            effective = combine_and(expr_copy(positional_assum), expr_copy(opt_assumptions));
        else
            effective = eval_and_free(expr_copy(opt_assumptions));
    } else {
        Expr* dollar = read_dollar_assumptions();
        if (positional_assum)
            effective = combine_and(expr_copy(positional_assum), dollar);
        else
            effective = dollar;
    }

    AssumeCtx* ctx = assume_ctx_from_expr(effective);
    /* Inconsistent assumptions: every statement is vacuously decided. We keep
     * the conservative behaviour of returning the expression rewritten under
     * whatever facts survive rather than fabricating a value. */
    Expr* Aconj = relational_assumption_conj(ctx);

    double budget = parse_time_budget(opt_timeconstraint);
    RefineBudget bud;
    bud.deadline   = (budget < HUGE_VAL) ? simp_mono_seconds() + budget : HUGE_VAL;
    bud.calls_left = REFINE_MAX_ENTAILMENT_CALLS;

    /* Bound any Simplify nested inside Reduce/CAD by the same budget. */
    double saved_budget = simp_current_time_budget();
    simp_set_time_budget(budget);

    Expr* result = NULL;

    /* (1) Predicate expression -> try to decide True/False. */
    const char* eh = head_name(expr);
    if (eh && (eh == SYM_Element || is_relational_head(eh))) {
        result = refine_decide_predicate(expr, ctx, Aconj, &bud);
    }

    /* (2) Rewrite path. */
    if (!result) {
        int have_facts = (ctx && ctx->count > 0);
        if (!have_facts) {
            /* Nothing to assume: identity. Steal the positional expr. */
            simp_set_time_budget(saved_budget);
            expr_free(Aconj);
            assume_ctx_free(ctx);
            expr_free(effective);
            /* Steal the positional expr out of res (the evaluator frees res). */
            for (size_t i = 0; i < argc; i++) {
                if (res->data.function.args[i] == expr) { res->data.function.args[i] = NULL; break; }
            }
            return expr;
        }

        Expr* rewritten = apply_assumption_rules(expr, ctx);
        if (!rewritten) rewritten = expr_copy(expr);

        int changed = 0;
        Expr* deep = deep_positivity_walk(rewritten, ctx, Aconj, &bud, &changed);
        expr_free(rewritten);

        /* Directed infinities: settle a symbolic direction to its unit. */
        Expr* dirs = refine_directed_infinity(deep, ctx, Aconj, &bud);
        expr_free(deep);

        result = eval_and_free(dirs);
    }

    simp_set_time_budget(saved_budget);
    expr_free(Aconj);
    assume_ctx_free(ctx);
    expr_free(effective);
    return result;
}

/* ------------------------------------------------------------------------- */
/* Registration                                                              */
/* ------------------------------------------------------------------------- */

void refine_init(void) {
    symtab_add_builtin("Refine", builtin_refine);
    symtab_get_def("Refine")->attributes |= ATTR_PROTECTED;

    /* Options[Refine] = {Assumptions -> $Assumptions, TimeConstraint -> 30}. */
    Expr* r_assum = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol(SYM_Assumptions),
                   expr_new_symbol(SYM_DollarAssumptions) }, 2);
    Expr* r_time = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol("TimeConstraint"), expr_new_integer(30) }, 2);
    Expr* opts = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ r_assum, r_time }, 2);
    symtab_set_options("Refine", opts);
}
