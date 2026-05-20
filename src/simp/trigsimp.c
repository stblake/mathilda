#include "trigsimp.h"
#include "attr.h"
#include "eval.h"
#include "facpoly.h"
#include "parse.h"
#include "symtab.h"
#include "sym_names.h"
#include "trig_canon.h"
#include <stdlib.h>
#include <string.h>

/* Per-call memo helpers (Phase 9): builtin_trigexpand and
 * builtin_trigfactor are invoked many times during the Simplify
 * candidate-set search.  When Simplify pushes a FactorMemo we cache
 * results here, keyed on `TrigExpand[X]` / `TrigFactor[X]`.  These
 * keys never collide with `Factor[X]` / `TrigRoundtrip[X]` because
 * the head is part of the hash.
 *
 * The wrapper pattern: the public entry point looks up in the memo,
 * delegates to the impl on miss, then stores the result.  The impl
 * does the actual work.  We don't wrap the threading-over-logic-
 * heads short-circuit because that recursively re-invokes the
 * builtin, which itself goes through the memo. */
/* Compute a canonical form of `arg` for memo lookup.  Mathilda's
 * evaluator handles Plus/Times orderless+flat normalisation but does
 * not auto-distribute scalars (`1/8 (-18 a + 6 b)` is structurally
 * different from `-9/4 a + 3/4 b`) or auto-combine fractions
 * (`(s + t)/(c + d)` vs `(s/(c+d)) + (t/(c+d))`).  Running
 * Together(Expand(.)) once unifies these surface forms.
 *
 * Note: an experimental switch to TrigToExp(arg) keying (which folds
 * in trig identities like `Cos[x]^2 = (1 + Cos[2x])/2`) was tried
 * but reverted -- TrigToExp's per-call cost (~5 ms) exceeded the
 * dedup gains for TrigFactor / TrigExpand on typical Simplify
 * intermediates, where most inputs are unique post-canonicalisation.
 * It worked for TrigRoundtrip (Phase 14) only because TrigToExp is
 * the first stage of TrigRoundtrip's pipeline anyway -- the cost
 * was free-rider.  For TrigFactor / TrigExpand the first stage is
 * different, so the TrigToExp computation is purely overhead. */
static Expr* trig_canonicalize(Expr* arg) {
    Expr* expand_args[1] = { expr_copy(arg) };
    Expr* expand_call = expr_new_function(expr_new_symbol("Expand"),
                                          expand_args, 1);
    Expr* expanded = evaluate(expand_call);
    expr_free(expand_call);

    Expr* together_args[1] = { expanded };
    Expr* together_call = expr_new_function(expr_new_symbol("Together"),
                                            together_args, 1);
    Expr* canonical = evaluate(together_call);
    expr_free(together_call);
    return canonical;
}

static Expr* trig_memo_call(Expr* res, const char* name,
                            Expr* (*impl)(Expr*)) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return impl(res);
    }
    FactorMemo* memo = factor_memo_active();
    if (!memo) return impl(res);

    Expr* arg = res->data.function.args[0];
    Expr* canonical = trig_canonicalize(arg);

    Expr* key_args[1] = { canonical };  /* takes ownership */
    Expr* key = expr_new_function(expr_new_symbol(name), key_args, 1);

    const Expr* hit = factor_memo_lookup(memo, key);
    if (hit) {
        Expr* cached = expr_copy((Expr*)hit);
        expr_free(key);
        return cached;
    }

    /* On a cache miss, pass the ORIGINAL input (not the canonical
     * form) to impl: TrigFactor / TrigExpand expect trig atoms, not
     * exponential forms, so calling impl on `TrigToExp(input)` would
     * produce a useless un-factored result.  Caching under the
     * canonical key + computing on the original is correct: any two
     * inputs that share the canonical key are algebraically equal,
     * so impl(input1) and impl(input2) are mathematically equivalent
     * even if they differ in surface form.  We accept the cached
     * result as the canonical output. */
    Expr* result = impl(res);

    if (result) factor_memo_store(memo, key, result);
    expr_free(key);
    return result;
}

/*
 * trigsimp.c
 *
 * Trigonometric simplification builtins:
 *   TrigToExp  -- rewrite circular/hyperbolic (and their inverses) as Exp/Log.
 *   ExpToTrig  -- rewrite Exp/Log forms back into circular/hyperbolic functions.
 *   TrigExpand -- expand sums and integer multiples that appear in arguments of
 *                 trigonometric (circular and hyperbolic) functions.
 */

static Expr* trig_to_exp_rules = NULL;

static Expr* exp_to_trig_rules = NULL;
static Expr* exp_to_trig_simp = NULL;

static Expr* trig_expand_rules = NULL;
static Expr* trig_expand_pythag = NULL;

static Expr* trig_factor_to_sincos = NULL;
static Expr* trig_factor_identities = NULL;
static Expr* trig_factor_from_sincos = NULL;

static Expr* trig_reduce_rules = NULL;
static Expr* trig_reduce_collapse = NULL;

Expr* builtin_exptotrig(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* Disable trig_canon so the intermediate Sin/Cos forms produced by
     * exp_to_trig_rules (e.g. Cosh[x] -> (Exp[x] + Exp[-x])/2 unwound) are
     * not prematurely collapsed into Tan/Sec/Csc shapes that the rule set
     * doesn't recognise. */
    trig_canon_suppress_inc();

    Expr* replace_args[2] = { expr_copy(arg), expr_copy(exp_to_trig_rules) };
    Expr* replace_expr = expr_new_function(expr_new_symbol("ReplaceRepeated"), replace_args, 2);
    Expr* replaced = evaluate(replace_expr);
    expr_free(replace_expr);

    Expr* together_args[1] = { replaced };
    Expr* together_expr = expr_new_function(expr_new_symbol("Together"), together_args, 1);
    Expr* togethered = evaluate(together_expr);
    expr_free(together_expr);

    Expr* cancel_args[1] = { togethered };
    Expr* cancel_expr = expr_new_function(expr_new_symbol("Cancel"), cancel_args, 1);
    Expr* cancelled = evaluate(cancel_expr);
    expr_free(cancel_expr);

    Expr* rep_args[2] = { cancelled, expr_copy(exp_to_trig_simp) };
    Expr* rep_simp = expr_new_function(expr_new_symbol("ReplaceRepeated"), rep_args, 2);
    Expr* result = evaluate(rep_simp);
    expr_free(rep_simp);

    trig_canon_suppress_dec();
    return result;
}

Expr* builtin_trigtoexp(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* Memo lookup keyed on the raw input.  TrigToExp is pure
     * (ReplaceAll + Expand) and shows up both as a Simplify
     * candidate-set transform AND inside transform_trig_roundtrip's
     * miss path.  On Tan-rich inputs each call costs 30-130 ms; the
     * candidate-set search re-enters with structurally identical
     * subexpressions across iterations.  Caching under the raw input
     * lets the second call (and any TrigRoundtrip internal call on
     * the same expression) hit immediately.
     *
     * The result is keyed under a `TrigToExp` head so it never
     * collides with the `TrigRoundtrip` / `TrigFactor` / `Factor`
     * entries that share the same FactorMemo. */
    FactorMemo* memo = factor_memo_active();
    Expr* memo_key = NULL;
    if (memo) {
        Expr* key_args[1] = { expr_copy(arg) };
        memo_key = expr_new_function(expr_new_symbol("TrigToExp"),
                                     key_args, 1);
        const Expr* hit = factor_memo_lookup(memo, memo_key);
        if (hit) {
            Expr* cached = expr_copy((Expr*)hit);
            expr_free(memo_key);
            return cached;
        }
    }

    /* Suppress trig_canon for the duration: Tan[x] under trig_to_exp_rules
     * rewrites to (Sin[x]/Cos[x] expressed via E^(I x)). Without suppression,
     * the Sin/Cos rewrite would prematurely re-collapse to Tan, defeating the
     * rule. */
    trig_canon_suppress_inc();

    Expr* replace_args[2] = { expr_copy(arg), expr_copy(trig_to_exp_rules) };
    Expr* replace_expr = expr_new_function(expr_new_symbol("ReplaceAll"), replace_args, 2);
    Expr* replaced = evaluate(replace_expr);
    expr_free(replace_expr);

    Expr* expand_args[1] = { replaced };
    Expr* expand_expr = expr_new_function(expr_new_symbol("Expand"), expand_args, 1);
    Expr* result = evaluate(expand_expr);
    expr_free(expand_expr);

    trig_canon_suppress_dec();

    if (memo_key) {
        factor_memo_store(memo, memo_key, result);
        expr_free(memo_key);
    }

    return result;
}

/* Heads over which TrigExpand threads (in addition to List, which is handled
 * via the ATTR_LISTABLE attribute). Mathematica's TrigExpand threads over
 * equations, inequalities, and logic functions.
 */
static int trigexpand_threads_over(const char* name) {
    return strcmp(name, "Equal") == 0 ||
           strcmp(name, "Unequal") == 0 ||
           strcmp(name, "Less") == 0 ||
           strcmp(name, "LessEqual") == 0 ||
           strcmp(name, "Greater") == 0 ||
           strcmp(name, "GreaterEqual") == 0 ||
           strcmp(name, "SameQ") == 0 ||
           strcmp(name, "UnsameQ") == 0 ||
           strcmp(name, "And") == 0 ||
           strcmp(name, "Or") == 0 ||
           strcmp(name, "Not") == 0 ||
           strcmp(name, "Xor") == 0 ||
           strcmp(name, "Implies") == 0;
}

/*
 * Detect whether the input contains a Pythagorean-eligible squared structure
 * of the form Sin[a]^k together with Cos[a]^k (or the hyperbolic analog
 * Sinh[a]^k together with Cosh[a]^k) for some k >= 2 with the same argument
 * a. The Factor-based Pythagorean collapse pass below only produces useful
 * work when the input has this shape; without this guard, Factor on a
 * multivariate non-Pythagorean expansion (e.g. TrigExpand[Sin[2 x + 3 y]]
 * after expansion becomes a 4-variable polynomial in Cos[x], Sin[x], Cos[y],
 * Sin[y]) becomes prohibitively slow.
 */
#define TRIG_PAIR_CAP 64
typedef struct {
    int kinds[TRIG_PAIR_CAP];   /* 0=Sin, 1=Cos, 2=Sinh, 3=Cosh */
    const Expr* args[TRIG_PAIR_CAP];
    size_t count;
    int overflow;
} TrigSquareList;

static int trig_pair_kind(const char* name) {
    if (strcmp(name, "Sin") == 0)  return 0;
    if (strcmp(name, "Cos") == 0)  return 1;
    if (strcmp(name, "Sinh") == 0) return 2;
    if (strcmp(name, "Cosh") == 0) return 3;
    return -1;
}

static void collect_trig_squares(const Expr* e, TrigSquareList* L) {
    if (!e || L->overflow) return;
    if (e->type != EXPR_FUNCTION) return;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (base->type == EXPR_FUNCTION &&
            base->data.function.head &&
            base->data.function.head->type == EXPR_SYMBOL &&
            base->data.function.arg_count == 1 &&
            exp->type == EXPR_INTEGER && exp->data.integer >= 2) {
            int k = trig_pair_kind(base->data.function.head->data.symbol);
            if (k >= 0) {
                if (L->count < TRIG_PAIR_CAP) {
                    L->kinds[L->count] = k;
                    L->args[L->count]  = base->data.function.args[0];
                    L->count++;
                } else {
                    L->overflow = 1;
                }
            }
        }
    }
    collect_trig_squares(e->data.function.head, L);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        collect_trig_squares(e->data.function.args[i], L);
    }
}

static int input_has_pythag_pair(const Expr* e) {
    TrigSquareList L = {{0}, {0}, 0, 0};
    collect_trig_squares(e, &L);
    /* Conservatively run Factor when the buffer overflowed -- the input is
     * complex enough that we cannot cheaply rule out a Pythagorean pair. */
    if (L.overflow) return 1;
    for (size_t i = 0; i < L.count; i++) {
        for (size_t j = i + 1; j < L.count; j++) {
            int ki = L.kinds[i], kj = L.kinds[j];
            int circ = (ki == 0 && kj == 1) || (ki == 1 && kj == 0);
            int hyp  = (ki == 2 && kj == 3) || (ki == 3 && kj == 2);
            if ((circ || hyp) && expr_eq((Expr*)L.args[i], (Expr*)L.args[j])) {
                return 1;
            }
        }
    }
    return 0;
}

/*
 * Count the number of distinct trig atoms (Sin/Cos/Sinh/Cosh applied to a
 * specific argument) that appear at integer power >= 2 anywhere in the
 * expression. Two atoms are considered the same if they have the same head
 * (Sin/Cos/Sinh/Cosh) and structurally equal arguments.
 *
 * This is used as a proxy for the cost of running Factor: heuristic_factor's
 * factor_roots routine probes O(v_count * num_high_deg_vars * 10) trial
 * divisions, each of which is itself an expensive multivariate polynomial
 * GCD. With more than two distinct trig atoms appearing at squared-or-higher
 * degree, the polynomial Factor sees is dense enough in trig "variables" that
 * the trial divisions run prohibitively slowly without producing useful
 * factorizations (since the polynomial is irreducible over Z without the
 * Sin^2 + Cos^2 = 1 identity).
 */
static size_t count_distinct_squared_trig_atoms(const Expr* e) {
    TrigSquareList L = {{0}, {0}, 0, 0};
    collect_trig_squares(e, &L);
    /* Conservative on overflow: report a high count so callers skip Factor. */
    if (L.overflow) return TRIG_PAIR_CAP;
    size_t distinct = 0;
    for (size_t i = 0; i < L.count; i++) {
        int already_seen = 0;
        for (size_t j = 0; j < i; j++) {
            if (L.kinds[i] == L.kinds[j] &&
                expr_eq((Expr*)L.args[i], (Expr*)L.args[j])) {
                already_seen = 1;
                break;
            }
        }
        if (!already_seen) distinct++;
    }
    return distinct;
}

/* Count every Sin/Cos/Sinh/Cosh call in the tree (any power or none).
 * Two calls match if same head and structurally equal argument. Used as
 * a gate against feeding Factor a 4+-variable trig polynomial -- when
 * TrigRoundtrip on a sub-expression produces forms like
 *   Cosh[2 x] - Cosh[2 x] Cosh[4 x] - Sinh[2 x] - Sinh[4 x] + ...
 * (4 distinct trig atoms at degree 1, dense cross-products), Factor
 * stalls for tens of seconds even though no atom appears squared. */
static void collect_trig_atoms(const Expr* e, TrigSquareList* L) {
    if (!e || L->overflow) return;
    if (e->type != EXPR_FUNCTION) return;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.arg_count == 1) {
        int k = trig_pair_kind(e->data.function.head->data.symbol);
        if (k >= 0) {
            if (L->count < TRIG_PAIR_CAP) {
                L->kinds[L->count] = k;
                L->args[L->count]  = e->data.function.args[0];
                L->count++;
            } else {
                L->overflow = 1;
                return;
            }
        }
    }
    collect_trig_atoms(e->data.function.head, L);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        collect_trig_atoms(e->data.function.args[i], L);
    }
}

static size_t count_distinct_trig_atoms(const Expr* e) {
    TrigSquareList L = {{0}, {0}, 0, 0};
    collect_trig_atoms(e, &L);
    if (L.overflow) return TRIG_PAIR_CAP;
    size_t distinct = 0;
    for (size_t i = 0; i < L.count; i++) {
        int already_seen = 0;
        for (size_t j = 0; j < i; j++) {
            if (L.kinds[i] == L.kinds[j] &&
                expr_eq((Expr*)L.args[i], (Expr*)L.args[j])) {
                already_seen = 1;
                break;
            }
        }
        if (!already_seen) distinct++;
    }
    return distinct;
}

/* Threshold above which Factor on a polynomial whose variables are trig
 * atoms (any power, not just squared) stalls. The motivating case from
 * bottom-up Simplify is the post-TrigRoundtrip form of
 * Cosh[x]^2 - Sinh[x]^2 which has 4 distinct degree-1 trig atoms
 * (Cosh[2x], Sinh[2x], Cosh[4x], Sinh[4x]) -- Factor took ~150s before
 * this guard was added. With <= 3 distinct atoms Factor remains
 * tractable. */
#define TRIG_FACTOR_TOTAL_ATOM_THRESHOLD 3

/*
 * Threshold above which Factor's multivariate cost becomes prohibitive on
 * trig expansions. With <= 2 distinct squared trig atoms (the typical
 * Sin[x]/Cos[x] or Sinh[x]/Cosh[x] pair plus at most one extra), Factor stays
 * tractable; beyond that the multivariate polynomial densifies and Factor's
 * trial-division loop in factor_roots stalls.
 */
#define TRIG_FACTOR_ATOM_THRESHOLD 2

/*
 * Maximum integer exponent on any squared trig atom that we still feed to
 * Factor. Even when only Sin[x] and Cos[x] appear (atom count == 2), the
 * post-TrigExpand polynomial for Sin[n x] / Cos[n x] grows to total degree n,
 * and Factor's multivariate trial-division loop stalls hard on dense forms
 * of that shape. The known-good Pythagorean cases reach degree 4
 * (e.g. Sin[x]^4 + 2 Sin[x]^2 Cos[x]^2 + Cos[x]^4 -> 1), so a cap of 4
 * preserves them while ruling out the high-degree blow-ups produced by
 * TrigFactor's Path B on inputs like Sin[5 x] + Sin[7 x].
 */
#define TRIG_FACTOR_DEGREE_THRESHOLD 4

/*
 * Maximum integer exponent appearing on any Sin/Cos/Sinh/Cosh atom in `e`.
 * Returns 0 if no such Power[trig_atom, k>=2] subterm exists. Used together
 * with count_distinct_squared_trig_atoms to gate Factor in the TrigFactor
 * pipeline. */
static int64_t max_trig_atom_power(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    int64_t best = 0;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (base->type == EXPR_FUNCTION &&
            base->data.function.head &&
            base->data.function.head->type == EXPR_SYMBOL &&
            base->data.function.arg_count == 1 &&
            exp->type == EXPR_INTEGER && exp->data.integer >= 2 &&
            trig_pair_kind(base->data.function.head->data.symbol) >= 0) {
            best = exp->data.integer;
        }
    }
    int64_t hp = max_trig_atom_power(e->data.function.head);
    if (hp > best) best = hp;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int64_t cp = max_trig_atom_power(e->data.function.args[i]);
        if (cp > best) best = cp;
    }
    return best;
}

/*
 * Test whether an expression tree contains any Power[_, n] subterm with a
 * negative-integer exponent. This is how Mathilda represents reciprocals
 * (e.g. 1/x = Power[x, -1], a/b = Times[a, Power[b, -1]]). Expand does not
 * distribute across negative powers, so the Factor-based Pythagorean
 * collapse must not be applied to expressions that hide structure inside a
 * denominator -- otherwise we would permanently replace a canonical form
 * such as "Cos[x]^2 - Sin[x]^2" with its factored denominator form
 * "(Cos[x] + Sin[x])(Cos[x] - Sin[x])".
 */
static int has_reciprocal_power(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_INTEGER && exp->data.integer < 0) return 1;
    }
    if (has_reciprocal_power(e->data.function.head)) return 1;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_reciprocal_power(e->data.function.args[i])) return 1;
    }
    return 0;
}

static Expr* builtin_trigexpand_impl(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* Thread over equations, inequalities, and logic heads. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        trigexpand_threads_over(arg->data.function.head->data.symbol)) {
        size_t n = arg->data.function.arg_count;
        Expr** new_args = (Expr**)calloc(n, sizeof(Expr*));
        if (!new_args) return NULL;
        for (size_t i = 0; i < n; i++) {
            Expr* inner[1] = { expr_copy(arg->data.function.args[i]) };
            Expr* wrap = expr_new_function(expr_new_symbol("TrigExpand"), inner, 1);
            new_args[i] = evaluate(wrap);
            expr_free(wrap);
        }
        return expr_new_function(expr_copy(arg->data.function.head), new_args, n);
    }

    /* Suppress trig_canon for the duration of the pipeline: the expansion
     * rules rewrite Tan[2x] etc. into ratios of Sin/Cos, which the Times
     * canonicalizer would otherwise immediately re-collapse back to Tan. */
    trig_canon_suppress_inc();

    /* Apply the angle-addition and multiple-angle rules to a fixed point. */
    Expr* replace_args[2] = { expr_copy(arg), expr_copy(trig_expand_rules) };
    Expr* replace_expr = expr_new_function(expr_new_symbol("ReplaceRepeated"), replace_args, 2);
    Expr* replaced = evaluate(replace_expr);
    expr_free(replace_expr);

    /* Distribute products of sums so the result is a flat sum of monomials. */
    Expr* expand_args[1] = { replaced };
    Expr* expand_expr = expr_new_function(expr_new_symbol("Expand"), expand_args, 1);
    Expr* expanded = evaluate(expand_expr);
    expr_free(expand_expr);

    /*
     * Factor-based Pythagorean collapse.
     *
     * Expanding Sin[n*x]^2 + Cos[n*x]^2 produces (Sin[x]^2 + Cos[x]^2)^n in
     * distributed form -- a polynomial Factor recognises. Applying the
     * (Sin[x]^2 + Cos[x]^2)^n -> 1 rule after Factor collapses it, and the
     * final Expand restores the canonical monomial form for any non-collapsing
     * subexpressions (so e.g. "Cos[x]^2 - Sin[x]^2" survives unchanged).
     *
     * Expand does not distribute across negative powers, so applying Factor
     * to an expression with a denominator would leave the denominator in
     * factored form permanently (e.g. Tan[2 x] -> 2 Cos Sin / ((Cos+Sin)(Cos-Sin))
     * instead of 2 Cos Sin / (Cos^2 - Sin^2)). Skip the Factor-based pass for
     * such expressions and only run the direct-sum Pythagorean rules.
     */
    /* Skip the Factor pass entirely when the input has no Pythagorean-eligible
     * squared structure -- Factor on a multivariate non-Pythagorean expansion
     * (e.g. TrigExpand[Sin[2 x + 3 y]]) is prohibitively slow and would not
     * have produced any useful collapse. Also skip when the expanded form has
     * more than two distinct squared trig atoms; even if a Pythagorean pair
     * exists structurally, Factor on the resulting multivariate polynomial
     * stalls without producing any useful collapse (e.g. TrigFactor's Path B
     * may invoke TrigExpand on an already-expanded form like the polynomial
     * for Sin[2 x + 3 y]). */
    Expr* reduced;
    if (has_reciprocal_power(expanded) || !input_has_pythag_pair(arg) ||
        count_distinct_squared_trig_atoms(expanded) > TRIG_FACTOR_ATOM_THRESHOLD) {
        Expr* pyth_args[2] = { expanded, expr_copy(trig_expand_pythag) };
        Expr* pyth_expr = expr_new_function(expr_new_symbol("ReplaceRepeated"), pyth_args, 2);
        reduced = evaluate(pyth_expr);
        expr_free(pyth_expr);
    } else {
        Expr* factor_args[1] = { expanded };
        Expr* factor_expr = expr_new_function(expr_new_symbol("Factor"), factor_args, 1);
        Expr* factored = evaluate(factor_expr);
        expr_free(factor_expr);

        Expr* pyth_args[2] = { factored, expr_copy(trig_expand_pythag) };
        Expr* pyth_expr = expr_new_function(expr_new_symbol("ReplaceRepeated"), pyth_args, 2);
        reduced = evaluate(pyth_expr);
        expr_free(pyth_expr);
    }

    /* Re-expand so that any residual factored products are distributed back
     * into the canonical TrigExpand monomial form. */
    Expr* final_args[1] = { reduced };
    Expr* final_expr = expr_new_function(expr_new_symbol("Expand"), final_args, 1);
    Expr* result = evaluate(final_expr);
    expr_free(final_expr);

    trig_canon_suppress_dec();
    return result;
}

/*
 * Count the total number of atomic leaves (integers, reals, symbols, strings,
 * bigints) in an expression tree. Used as a complexity measure when choosing
 * between alternative TrigFactor results.
 */
static size_t trigfactor_leaf_count(const Expr* e) {
    if (!e) return 0;
    if (e->type != EXPR_FUNCTION) return 1;
    size_t total = trigfactor_leaf_count(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        total += trigfactor_leaf_count(e->data.function.args[i]);
    }
    return total;
}

/*
 * Core TrigFactor pipeline (see builtin_trigfactor for the overall design).
 * Takes ownership of `input` and returns a new Expr*.
 */
static Expr* trigfactor_run_pipeline(Expr* input) {
    /* See builtin_trigreduce_impl for the rationale: the to_sincos rewrite
     * needs trig_canon disabled so the Sin/Cos polynomial form survives
     * long enough for Factor and the identity rules to act on it. */
    trig_canon_suppress_inc();

    /* Rewrite reciprocal heads as Sin/Cos ratios. */
    Expr* ts_args[2] = { input, expr_copy(trig_factor_to_sincos) };
    Expr* ts_expr = expr_new_function(expr_new_symbol("ReplaceRepeated"), ts_args, 2);
    Expr* sincos = evaluate(ts_expr);
    expr_free(ts_expr);

    /* Combine over a common denominator. */
    Expr* tog_args[1] = { sincos };
    Expr* tog_expr = expr_new_function(expr_new_symbol("Together"), tog_args, 1);
    Expr* togethered = evaluate(tog_expr);
    expr_free(tog_expr);

    /* Factor -- Mathilda Factor treats trig atoms as polynomial variables.
     * Skip Factor when the polynomial would be too expensive to factor:
     *
     *   1. More than two distinct squared trig atoms: the multivariate
     *      polynomial has too many high-degree variables and Factor's
     *      trial-division loop in factor_roots stalls (e.g. on the
     *      polynomial that TrigExpand produces for Sin[2 x + 3 y]).
     *   2. Even with only two atoms (Sin[x] and Cos[x], the typical
     *      Pythagorean pair), an exponent above
     *      TRIG_FACTOR_DEGREE_THRESHOLD signals a dense high-degree
     *      polynomial. Path B of TrigFactor on Sin[m x] + Sin[n x]
     *      (m,n >= 3) produces exactly such a form -- TrigExpand turns
     *      Sin[7 x] into a degree-7 polynomial in Sin[x], Cos[x], and
     *      Factor stalls on it.
     *
     * In either case the identity rules applied below still match the
     * Pythagorean structure that survives without Factor (e.g.
     * (Sin[x]^2 + Cos[x]^2)(Cosh[y]^2 - Sinh[y]^2) collapses to 1 directly
     * via the Times-context identity rules). */
    Expr* factored;
    if (count_distinct_squared_trig_atoms(togethered) > TRIG_FACTOR_ATOM_THRESHOLD ||
        max_trig_atom_power(togethered) > TRIG_FACTOR_DEGREE_THRESHOLD ||
        count_distinct_trig_atoms(togethered) > TRIG_FACTOR_TOTAL_ATOM_THRESHOLD) {
        factored = togethered;
    } else {
        Expr* fac_args[1] = { togethered };
        Expr* fac_expr = expr_new_function(expr_new_symbol("Factor"), fac_args, 1);
        factored = evaluate(fac_expr);
        expr_free(fac_expr);
    }

    /* Apply identity collapse rules to a fixed point. */
    Expr* id_args[2] = { factored, expr_copy(trig_factor_identities) };
    Expr* id_expr = expr_new_function(expr_new_symbol("ReplaceRepeated"), id_args, 2);
    Expr* collapsed = evaluate(id_expr);
    expr_free(id_expr);

    /* Restore Tan/Sec/Cot/Csc (and hyperbolic analogs). */
    Expr* fs_args[2] = { collapsed, expr_copy(trig_factor_from_sincos) };
    Expr* fs_expr = expr_new_function(expr_new_symbol("ReplaceRepeated"), fs_args, 2);
    Expr* result = evaluate(fs_expr);
    expr_free(fs_expr);

    trig_canon_suppress_dec();
    return result;
}

/*
 * TrigFactor -- factor trigonometric expressions.
 *
 * TrigFactor is broadly the functional inverse of TrigExpand for the
 * structural identities that both support. It combines sums of trig
 * products into angle-sum forms, collapses Pythagorean identities in both
 * circular and hyperbolic form, recognises the reverse of the double-angle
 * identities, and factors polynomial structure over the trig atoms.
 *
 * Two parallel paths are tried and the shorter result (by leaf count) is
 * returned:
 *
 *   Path A   -- run the core pipeline on the original expression. This
 *               preserves angle-sum structure such as Cos[x+y] and is the
 *               preferred route for expressions like Sin[x+y]^2 + Tan[x+y]
 *               where expanding the angle-sum would make the result much
 *               larger.
 *   Path B   -- first apply TrigExpand to the argument, then run the core
 *               pipeline. This catches cancellations that only become
 *               visible after the angle-sum is expanded, e.g.
 *               Cos[x+y] + Sin[x] Sin[y] -> Cos[x] Cos[y].
 *
 * The core pipeline itself:
 *   1. Rewrite reciprocal heads (Tan, Cot, Sec, Csc and their hyperbolic
 *      analogs) in terms of Sin/Cos/Sinh/Cosh so Factor sees the full
 *      polynomial structure.
 *   2. Together to combine over a common denominator.
 *   3. Factor the rational.
 *   4. Apply identity collapse rules: Pythagorean sums, reverse
 *      angle-addition, reverse double-angle, and factored-form Pythagorean
 *      identities that arise from Factor.
 *   5. Re-collapse Sin[x]/Cos[x] (and friends) back to Tan[x] etc. so
 *      surface-form reciprocals survive the round-trip.
 *
 * Threading over equations, inequalities, and logic heads is applied before
 * the pipeline (List threading is delivered by ATTR_LISTABLE).
 */
/* Names of the unary trig / hyperbolic atoms that TrigFactor and
 * TrigExpand operate on. */
static bool is_trig_head(const char* h) {
    return strcmp(h, "Sin") == 0 || strcmp(h, "Cos") == 0 ||
           strcmp(h, "Tan") == 0 || strcmp(h, "Cot") == 0 ||
           strcmp(h, "Sec") == 0 || strcmp(h, "Csc") == 0 ||
           strcmp(h, "Sinh") == 0 || strcmp(h, "Cosh") == 0 ||
           strcmp(h, "Tanh") == 0 || strcmp(h, "Coth") == 0 ||
           strcmp(h, "Sech") == 0 || strcmp(h, "Csch") == 0;
}

/* True if `e` is a unary trig atom: Sin[arg], Cos[arg], etc. */
static bool is_trig_atom(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && is_trig_head(e->data.function.head->data.symbol)
        && e->data.function.arg_count == 1;
}

/* Does the expression contain a structural pattern that TrigExpand
 * could plausibly turn into a cancellable form?  Specifically: a
 * Power[trig_atom, k] with k >= 2, or a Times of two or more distinct
 * trig atoms (or trig atoms times other things).  Plain
 * sums/products with bare trig atoms (no powers, no trig-trig
 * products) won't gain cancellations from TrigExpand -- the expansion
 * just reinverts after factoring.
 *
 * We walk recursively through Plus, Times, Power, and division
 * structure, returning true at the first cancellable shape we find. */
static bool has_compound_trig_structure(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL) return false;

    const char* h = e->data.function.head->data.symbol;

    /* Power[trig_atom, k] with k != 1: cancellable. */
    if (h == SYM_Power && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        if (is_trig_atom(base) && exp->type == EXPR_INTEGER &&
            (exp->data.integer >= 2 || exp->data.integer <= -2)) {
            return true;
        }
        /* Also recurse: Power[X, n] where X is compound. */
        return has_compound_trig_structure(base);
    }

    if (h == SYM_Times) {
        /* Count distinct trig atoms (ignoring numeric factors).  Two
         * or more trig atoms in a single Times product means
         * TrigExpand may expose product-to-sum cancellations. */
        int trig_count = 0;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            Expr* a = e->data.function.args[i];
            if (is_trig_atom(a)) trig_count++;
            else if (a->type == EXPR_FUNCTION
                     && a->data.function.head
                     && a->data.function.head->type == EXPR_SYMBOL
                     && a->data.function.head->data.symbol == SYM_Power
                     && a->data.function.arg_count == 2
                     && is_trig_atom(a->data.function.args[0])) {
                /* Power of a trig atom -- compound on its own. */
                return true;
            }
            else if (has_compound_trig_structure(a)) return true;
        }
        return trig_count >= 2;
    }

    if (h == SYM_Plus) {
        /* Sums by themselves of bare trig atoms aren't cancellable
         * (Sin[x] + Sin[3x] doesn't gain from TrigExpand→Factor).
         * Recurse only -- we're looking for compound structure
         * inside one of the terms. */
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (has_compound_trig_structure(e->data.function.args[i])) return true;
        }
        return false;
    }

    /* Trig atom on its own is NOT compound. */
    if (is_trig_atom(e)) return false;

    /* For other heads, recurse into args. */
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_compound_trig_structure(e->data.function.args[i])) return true;
    }
    return false;
}

static Expr* builtin_trigfactor_impl(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* 1. Thread over equations, inequalities, and logic heads. Listable is
     * handled automatically by the attribute system. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        trigexpand_threads_over(arg->data.function.head->data.symbol)) {
        size_t n = arg->data.function.arg_count;
        Expr** new_args = (Expr**)calloc(n, sizeof(Expr*));
        if (!new_args) return NULL;
        for (size_t i = 0; i < n; i++) {
            Expr* inner[1] = { expr_copy(arg->data.function.args[i]) };
            Expr* wrap = expr_new_function(expr_new_symbol("TrigFactor"), inner, 1);
            new_args[i] = evaluate(wrap);
            expr_free(wrap);
        }
        return expr_new_function(expr_copy(arg->data.function.head), new_args, n);
    }

    /* Path A: pipeline applied to the argument as-is. */
    Expr* result_a = trigfactor_run_pipeline(expr_copy(arg));

    /* If Path A changed the expression at all, trust that result and skip
     * Path B. Path B (which TrigExpand-expands the input) can be expensive
     * for expressions with angle-sum arguments -- e.g. Sin[x+y]^2 + Tan[x+y]
     * blows up into a multivariate rational that Factor has to reduce over
     * many terms. Path B is only valuable when Path A was unproductive,
     * because TrigExpand may then expose cancellations (for example,
     * Cos[x+y] + Sin[x] Sin[y] -> Cos[x] Cos[y]). */
    if (!expr_eq(result_a, arg)) {
        return result_a;
    }

    /* Cheap-skip heuristic: when Path A is a no-op AND the input has
     * no compound trig structure (no Power[trig, k>=2], no
     * Times[trig, trig, ...], no nested compound forms), Path B's
     * TrigExpand → Factor → identity-rules cycle empirically reduces
     * to the input itself.  These are typical Simplify intermediates
     * like Sin[x] + Sin[3x] over 2(Cos[x]+Cos[3x]), where TrigExpand
     * just reinverts after factoring.  Skipping them saves 50-250 ms
     * per call without losing any genuine cancellation: the user's
     * primary case Sin[x]^3 + Sin[3x] - 3 Sin[x] still benefits
     * because it has Sin[x]^3 (a Power[trig, 3]).
     *
     * The predicate is conservative: only inputs we are confident
     * about get the skip.  Anything with a power, product, or
     * compound substructure stays on Path B. */
    if (!has_compound_trig_structure(arg)) {
        return result_a;
    }

    /* Path B: apply TrigExpand first, then run the pipeline. */
    Expr* te_args[1] = { expr_copy(arg) };
    Expr* te_expr = expr_new_function(expr_new_symbol("TrigExpand"), te_args, 1);
    Expr* expanded = evaluate(te_expr);
    expr_free(te_expr);
    Expr* result_b = trigfactor_run_pipeline(expanded);

    /* Keep Path B only if it actually produced a simpler expression than
     * Path A (which in this branch equals the input). This protects against
     * the TrigExpand round-trip blowing the expression up. */
    if (trigfactor_leaf_count(result_b) < trigfactor_leaf_count(result_a)) {
        expr_free(result_a);
        return result_b;
    }
    expr_free(result_b);
    return result_a;
}

/*
 * TrigReduce -- inverse of TrigExpand for the product/power direction.
 *
 * Rewrites products and integer powers of single-argument circular and
 * hyperbolic trig calls into single trig calls of compound (sum/multiple)
 * arguments. The classical product-to-sum and power-reduction identities:
 *
 *   Sin[a] Cos[b]   = (Sin[a + b] + Sin[a - b]) / 2
 *   Cos[a] Cos[b]   = (Cos[a + b] + Cos[a - b]) / 2
 *   Sin[a] Sin[b]   = (Cos[a - b] - Cos[a + b]) / 2
 *   Sin[x]^2        = (1 - Cos[2 x]) / 2
 *   Cos[x]^2        = (1 + Cos[2 x]) / 2
 * (and hyperbolic analogues).
 *
 * Pipeline:
 *   1. Convert reciprocal heads (Tan/Cot/Sec/Csc and hyperbolic analogues)
 *      to Sin/Cos ratios so the product-to-sum rules can see ratios.
 *   2. Iterate (ReplaceRepeated rules, Expand) until fixed point. The
 *      iteration is necessary because Expand re-exposes Cos[2 x]^2 terms
 *      hidden inside (1 - Cos[2 x])^2 / 4 after a power-reduction pass on
 *      Sin[x]^4, etc.
 *   3. Together to renormalise denominators (e.g. Tan[x] + Tan[y] becomes
 *      a single rational).
 *   4. Apply angle-addition collapse rules so numerators of the form
 *      Sin[a] Cos[b] + Cos[a] Sin[b] become Sin[a + b].
 *   5. Restore reciprocal trig (Tan/Sec/Csc, etc.) where possible.
 *
 * Threading over List is delivered by ATTR_LISTABLE; equations,
 * inequalities, and logical operators are threaded explicitly here using
 * the same helper as TrigExpand / TrigFactor.
 */
static Expr* trigreduce_apply_rules(Expr* input, Expr* rules) {
    Expr* args[2] = { input, expr_copy(rules) };
    Expr* call = expr_new_function(expr_new_symbol("ReplaceRepeated"),
                                   args, 2);
    Expr* result = evaluate(call);
    expr_free(call);
    return result;
}

static Expr* trigreduce_call_unary(const char* name, Expr* input) {
    Expr* args[1] = { input };
    Expr* call = expr_new_function(expr_new_symbol(name), args, 1);
    Expr* result = evaluate(call);
    expr_free(call);
    return result;
}

static Expr* builtin_trigreduce_impl(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;
    Expr* arg = res->data.function.args[0];

    /* Thread over equations, inequalities, and logic heads. List
     * threading is delivered by ATTR_LISTABLE in the attribute
     * system. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        trigexpand_threads_over(arg->data.function.head->data.symbol)) {
        size_t n = arg->data.function.arg_count;
        Expr** new_args = (Expr**)calloc(n, sizeof(Expr*));
        if (!new_args) return NULL;
        for (size_t i = 0; i < n; i++) {
            Expr* inner[1] = { expr_copy(arg->data.function.args[i]) };
            Expr* wrap = expr_new_function(expr_new_symbol("TrigReduce"),
                                           inner, 1);
            new_args[i] = evaluate(wrap);
            expr_free(wrap);
        }
        return expr_new_function(expr_copy(arg->data.function.head),
                                 new_args, n);
    }

    /* Suppress Times/Power-level trig canonicalization for the entire
     * pipeline: step 1's Tan -> Sin/Cos rewrite would otherwise be undone
     * immediately when its result evaluates back through Times. The outer
     * evaluator re-evaluates whatever we return, at which point trig_canon
     * fires on the final Sin/Cos result and (together with step 5's explicit
     * from_sincos rules) yields a clean Tan/Sec/Csc-named answer. */
    trig_canon_suppress_inc();

    /* Step 1: rewrite reciprocal heads so the polynomial structure is
     * visible to the product-to-sum rules. */
    Expr* current = expr_copy(arg);
    current = trigreduce_apply_rules(current, trig_factor_to_sincos);

    /* Step 2: iterate (ReplaceRepeated, Expand) until fixed point. The
     * safety bound caps work on pathological inputs. Sin[x]^(2^k) needs
     * k iterations to fully reduce, so 16 covers exponents up through
     * Sin[x]^65536 -- well beyond any realistic input. */
    Expr* prev = NULL;
    for (int iter = 0; iter < 16; iter++) {
        Expr* before = expr_copy(current);
        current = trigreduce_apply_rules(current, trig_reduce_rules);
        current = trigreduce_call_unary("Expand", current);
        int converged = expr_eq(before, current);
        expr_free(before);
        if (prev) expr_free(prev);
        prev = NULL;
        if (converged) break;
    }
    if (prev) expr_free(prev);

    /* Step 3: combine over a common denominator so collapse rules can
     * see numerators as a single Plus. */
    current = trigreduce_call_unary("Together", current);

    /* Step 4: angle-addition collapse rules. Idempotent on inputs
     * lacking the matching shape. */
    current = trigreduce_apply_rules(current, trig_reduce_collapse);

    /* Step 5: restore Tan/Cot/Sec/Csc (and hyperbolic analogues) where
     * the result has the corresponding ratio or reciprocal shape. */
    current = trigreduce_apply_rules(current, trig_factor_from_sincos);

    /* Step 6: final canonicalisation. Expand distributes any remaining
     * outer scalar (e.g. `1/2 (2 Cos[a+b] + 2 Sin[a+b])` flattens to
     * `Cos[a+b] + Sin[a+b]` which Together leaves untouched, while
     * irreducible fractions like `(3 - 4 Cos[2x] + Cos[4x])/2` are
     * preserved as a single rational by the trailing Together). */
    current = trigreduce_call_unary("Expand", current);
    current = trigreduce_call_unary("Together", current);

    trig_canon_suppress_dec();
    return current;
}

/* Public entry points: route through trig_memo_call so the FactorMemo
 * (when active inside a Simplify call) caches results, then delegate
 * to the implementation on a cache miss. */
Expr* builtin_trigexpand(Expr* res) {
    return trig_memo_call(res, "TrigExpand", builtin_trigexpand_impl);
}

Expr* builtin_trigfactor(Expr* res) {
    return trig_memo_call(res, "TrigFactor", builtin_trigfactor_impl);
}

Expr* builtin_trigreduce(Expr* res) {
    return trig_memo_call(res, "TrigReduce", builtin_trigreduce_impl);
}

void trigsimp_init(void) {
    const char* rules_str = "{ "
        "Sin[x_] :> 1/2 I E^(-I x) - 1/2 I E^(I x), "
        "Cos[x_] :> E^(-I x)/2 + E^(I x)/2, "
        "Tan[x_] :> (I (E^(-I x) - E^(I x)))/(E^(-I x) + E^(I x)), "
        "Cot[x_] :> (-I (E^(-I x) + E^(I x)))/(E^(-I x) - E^(I x)), "
        "Sec[x_] :> 2/(E^(-I x) + E^(I x)), "
        "Csc[x_] :> (2 I)/(E^(I x) - E^(-I x)), "
        "Sinh[x_] :> -(E^-x/2) + E^x/2, "
        "Cosh[x_] :> E^-x/2 + E^x/2, "
        "Tanh[x_] :> -(E^-x/(E^-x+E^x)) + E^x/(E^-x+E^x), "
        "Coth[x_] :> E^-x/(E^-x-E^x) + E^x/(E^-x-E^x), "
        "Sech[x_] :> 2/(E^-x + E^x), "
        "Csch[x_] :> 2/(-E^-x + E^x), "
        "ArcSin[x_] :> -I Log[I x + Sqrt[1 - x^2]], "
        "ArcCos[x_] :> Pi/2 + I Log[I x + Sqrt[1 - x^2]], "
        "ArcTan[x_] :> 1/2 I Log[1 - I x] - 1/2 I Log[1 + I x], "
        "ArcTan[x_, y_] :> -I Log[(x + I y)/Sqrt[x^2 + y^2]], "
        "ArcCot[x_] :> 1/2 I Log[1 - I/x] - 1/2 I Log[1 + I/x], "
        "ArcSec[x_] :> Pi/2 + I Log[I/x + Sqrt[1 - 1/x^2]], "
        "ArcCsc[x_] :> -I Log[I/x + Sqrt[1 - 1/x^2]], "
        "ArcSinh[x_] :> Log[x + Sqrt[1 + x^2]], "
        "ArcCosh[x_] :> Log[x + Sqrt[-1 + x] Sqrt[1 + x]], "
        "ArcTanh[x_] :> -(1/2) Log[1 - x] + 1/2 Log[1 + x], "
        "ArcCoth[x_] :> -(1/2) Log[1 - 1/x] + 1/2 Log[1 + 1/x], "
        "ArcSech[x_] :> Log[1/x + Sqrt[-1 + 1/x] Sqrt[1 + 1/x]], "
        "ArcCsch[x_] :> Log[1/x + Sqrt[1 + 1/x^2]] "
    "}";

    trig_to_exp_rules = parse_expression(rules_str);

    const char* exp_to_trig_rules_str = "{ "
        "E^((I) x_.) :> Cos[x] + I Sin[x], "
        "E^((-I) x_.) :> Cos[x] - I Sin[x], "
        "E^(-x_.) :> Cosh[x] - Sinh[x], "
        "E^x_. :> Cosh[x] + Sinh[x], "
        "Log[a_ + (I) x_] - Log[a_ - (I) x_] :> 2 I ArcTan[x/a], "
        "Log[a_ - (I) x_] - Log[a_ + (I) x_] :> -2 I ArcTan[x/a], "
        "Log[a_ + x_] - Log[a_ - x_] :> 2 ArcTanh[x/a], "
        "Log[a_ - x_] - Log[a_ + x_] :> -2 ArcTanh[x/a], "
        "Log[x_ + Sqrt[1 + x_^2]] :> ArcSinh[x], "
        "Log[x_ + Sqrt[-1 + x_] Sqrt[1 + x_]] :> ArcCosh[x], "
        "Log[(I) x_ + Sqrt[1 - x_^2]] :> (I) ArcSin[x], "
        "Log[1/x_ + Sqrt[1 + 1/x_^2]] :> ArcCsch[x], "
        "Log[1/x_ + Sqrt[-1 + 1/x_] Sqrt[1 + 1/x_]] :> ArcSech[x], "
        "Log[(I)/x_ + Sqrt[1 - 1/x_^2]] :> (I) ArcCsc[x] "
    "}";

    const char* exp_to_trig_simp_str = "{ "
        "Sin[x_] * Cos[x_]^-1 :> Tan[x], "
        "Cos[x_] * Sin[x_]^-1 :> Cot[x], "
        "Sinh[x_] * Cosh[x_]^-1 :> Tanh[x], "
        "Cosh[x_] * Sinh[x_]^-1 :> Coth[x], "
        "Cos[x_]^-1 :> Sec[x], "
        "Sin[x_]^-1 :> Csc[x], "
        "Cosh[x_]^-1 :> Sech[x], "
        "Sinh[x_]^-1 :> Csch[x] "
    "}";

    exp_to_trig_rules = parse_expression(exp_to_trig_rules_str);
    exp_to_trig_simp = parse_expression(exp_to_trig_simp_str);

    /*
     * TrigExpand rules.
     *
     * 1. Angle-addition: split Sin[a + b + ...] into the binary form by
     *    matching the first summand and grouping the rest back up with Plus.
     *    Multi-term arguments are handled by iterated application of the
     *    binary formula via ReplaceRepeated.
     * 2. Multiple-angle: rewrite Sin[n x] in terms of Sin[(n-1) x] and
     *    Cos[(n-1) x] for integer n >= 2, so ReplaceRepeated reduces the
     *    argument down to Sin[x] and Cos[x]. Negative-integer arguments are
     *    folded out by the existing parity simplifications in trig.c / hyperbolic.c
     *    (e.g. Sin[-2 x] auto-evaluates to -Sin[2 x]).
     * 3. Tan/Cot/Sec/Csc (and their hyperbolic analogs) with a sum argument
     *    or an integer-multiple argument are rewritten in terms of Sin/Cos
     *    (resp. Sinh/Cosh) so the base rules can expand them.
     */
    const char* expand_rules_str = "{ "
        "Sin[x_ + y__] :> Sin[x] Cos[Plus[y]] + Cos[x] Sin[Plus[y]], "
        "Cos[x_ + y__] :> Cos[x] Cos[Plus[y]] - Sin[x] Sin[Plus[y]], "
        "Sinh[x_ + y__] :> Sinh[x] Cosh[Plus[y]] + Cosh[x] Sinh[Plus[y]], "
        "Cosh[x_ + y__] :> Cosh[x] Cosh[Plus[y]] + Sinh[x] Sinh[Plus[y]], "
        "Sin[n_Integer x_] /; n > 1 :> Sin[(n-1) x] Cos[x] + Cos[(n-1) x] Sin[x], "
        "Cos[n_Integer x_] /; n > 1 :> Cos[(n-1) x] Cos[x] - Sin[(n-1) x] Sin[x], "
        "Sinh[n_Integer x_] /; n > 1 :> Sinh[(n-1) x] Cosh[x] + Cosh[(n-1) x] Sinh[x], "
        "Cosh[n_Integer x_] /; n > 1 :> Cosh[(n-1) x] Cosh[x] + Sinh[(n-1) x] Sinh[x], "
        "Tan[x_ + y__] :> Sin[x + y] / Cos[x + y], "
        "Cot[x_ + y__] :> Cos[x + y] / Sin[x + y], "
        "Sec[x_ + y__] :> 1 / Cos[x + y], "
        "Csc[x_ + y__] :> 1 / Sin[x + y], "
        "Tanh[x_ + y__] :> Sinh[x + y] / Cosh[x + y], "
        "Coth[x_ + y__] :> Cosh[x + y] / Sinh[x + y], "
        "Sech[x_ + y__] :> 1 / Cosh[x + y], "
        "Csch[x_ + y__] :> 1 / Sinh[x + y], "
        "Tan[n_Integer x_] /; n > 1 :> Sin[n x] / Cos[n x], "
        "Cot[n_Integer x_] /; n > 1 :> Cos[n x] / Sin[n x], "
        "Sec[n_Integer x_] /; n > 1 :> 1 / Cos[n x], "
        "Csc[n_Integer x_] /; n > 1 :> 1 / Sin[n x], "
        "Tanh[n_Integer x_] /; n > 1 :> Sinh[n x] / Cosh[n x], "
        "Coth[n_Integer x_] /; n > 1 :> Cosh[n x] / Sinh[n x], "
        "Sech[n_Integer x_] /; n > 1 :> 1 / Cosh[n x], "
        "Csch[n_Integer x_] /; n > 1 :> 1 / Sinh[n x], "
        /* Inverse-trig compositions. After angle-addition expansion of
         * something like Sin[ArcSin[t] + ArcSin[v]], the residual cross-terms
         * Cos[ArcSin[v]] should reduce to Sqrt[1 - v^2] (and so on) so the
         * result is a polynomial expression in t, v rather than a mixed
         * inverse-trig form. */
        "Cos[ArcSin[x_]] :> Sqrt[1 - x^2], "
        "Sin[ArcCos[x_]] :> Sqrt[1 - x^2], "
        "Tan[ArcSin[x_]] :> x / Sqrt[1 - x^2], "
        "Cot[ArcSin[x_]] :> Sqrt[1 - x^2] / x, "
        "Sec[ArcSin[x_]] :> 1 / Sqrt[1 - x^2], "
        "Csc[ArcSin[x_]] :> 1 / x, "
        "Tan[ArcCos[x_]] :> Sqrt[1 - x^2] / x, "
        "Cot[ArcCos[x_]] :> x / Sqrt[1 - x^2], "
        "Sec[ArcCos[x_]] :> 1 / x, "
        "Csc[ArcCos[x_]] :> 1 / Sqrt[1 - x^2], "
        "Sin[ArcTan[x_]] :> x / Sqrt[1 + x^2], "
        "Cos[ArcTan[x_]] :> 1 / Sqrt[1 + x^2], "
        "Cot[ArcTan[x_]] :> 1 / x, "
        "Sec[ArcTan[x_]] :> Sqrt[1 + x^2], "
        "Csc[ArcTan[x_]] :> Sqrt[1 + x^2] / x, "
        "Sin[ArcCot[x_]] :> 1 / Sqrt[1 + x^2], "
        "Cos[ArcCot[x_]] :> x / Sqrt[1 + x^2], "
        "Tan[ArcCot[x_]] :> 1 / x, "
        "Sec[ArcCot[x_]] :> Sqrt[1 + x^2] / x, "
        "Csc[ArcCot[x_]] :> Sqrt[1 + x^2], "
        "Sin[ArcSec[x_]] :> Sqrt[1 - 1/x^2], "
        "Cos[ArcSec[x_]] :> 1 / x, "
        "Tan[ArcSec[x_]] :> Sqrt[x^2 - 1], "
        "Cot[ArcSec[x_]] :> 1 / Sqrt[x^2 - 1], "
        "Csc[ArcSec[x_]] :> 1 / Sqrt[1 - 1/x^2], "
        "Sin[ArcCsc[x_]] :> 1 / x, "
        "Cos[ArcCsc[x_]] :> Sqrt[1 - 1/x^2], "
        "Tan[ArcCsc[x_]] :> 1 / Sqrt[x^2 - 1], "
        "Cot[ArcCsc[x_]] :> Sqrt[x^2 - 1], "
        "Sec[ArcCsc[x_]] :> 1 / Sqrt[1 - 1/x^2], "
        /* Hyperbolic-of-inverse-hyperbolic compositions (real branch). */
        "Cosh[ArcSinh[x_]] :> Sqrt[1 + x^2], "
        "Sinh[ArcCosh[x_]] :> Sqrt[-1 + x^2], "
        "Tanh[ArcSinh[x_]] :> x / Sqrt[1 + x^2], "
        "Coth[ArcSinh[x_]] :> Sqrt[1 + x^2] / x, "
        "Sech[ArcSinh[x_]] :> 1 / Sqrt[1 + x^2], "
        "Csch[ArcSinh[x_]] :> 1 / x, "
        "Tanh[ArcCosh[x_]] :> Sqrt[-1 + x^2] / x, "
        "Coth[ArcCosh[x_]] :> x / Sqrt[-1 + x^2], "
        "Sech[ArcCosh[x_]] :> 1 / x, "
        "Csch[ArcCosh[x_]] :> 1 / Sqrt[-1 + x^2], "
        "Sinh[ArcTanh[x_]] :> x / Sqrt[1 - x^2], "
        "Cosh[ArcTanh[x_]] :> 1 / Sqrt[1 - x^2], "
        "Coth[ArcTanh[x_]] :> 1 / x, "
        "Sech[ArcTanh[x_]] :> Sqrt[1 - x^2], "
        "Csch[ArcTanh[x_]] :> Sqrt[1 - x^2] / x "
    "}";

    trig_expand_rules = parse_expression(expand_rules_str);

    /*
     * Pythagorean identities applied after Factor.
     *
     * Factor does not produce a canonical sign for these expressions: depending
     * on the input, Factor[Expand[Sin[n*x]^2 + Cos[n*x]^2]] may emerge as
     * (Sin^2 + Cos^2)^n, or (-Sin^2 - Cos^2)^n wrapped in an outer -1, and the
     * hyperbolic analog likewise may factor with the minus attached to either
     * Sinh or Cosh. We enumerate both signs; each rule fires once and the
     * remaining constant factors are picked up by Expand.
     *
     *   (Sin[x]^2 + Cos[x]^2)^n           circular Pythagorean identity,
     *   (-Sin[x]^2 - Cos[x]^2)^n          same with negated base;
     *   (Cosh+Sinh)^n (Cosh-Sinh)^n       hyperbolic Pythagorean via real
     *                                     factorization,
     *   (Cosh+Sinh)^n (-Cosh+Sinh)^n      same with negated second factor;
     *   r___ tails                        allow arbitrary extra factors in the
     *                                     surrounding Times (e.g. a leading
     *                                     constant that Factor pulled out),
     *   <base>+r___ sum rules             catch the identity when it appears as
     *                                     part of a larger sum that Factor left
     *                                     in the original summed form.
     */
    const char* pythag_str = "{ "
        "(Sin[x_]^2 + Cos[x_]^2)^n_. r___ :> 1 r, "
        "(-Sin[x_]^2 - Cos[x_]^2)^n_. r___ :> (-1)^n r, "
        "(Cosh[x_] + Sinh[x_])^n_. (Cosh[x_] - Sinh[x_])^n_. r___ :> 1 r, "
        "(Cosh[x_] + Sinh[x_])^n_. (-Cosh[x_] + Sinh[x_])^n_. r___ :> (-1)^n r, "
        "Sin[x_]^2 + Cos[x_]^2 + r___ :> 1 + r, "
        "Cosh[x_]^2 - Sinh[x_]^2 + r___ :> 1 + r, "
        "Sinh[x_]^2 - Cosh[x_]^2 + r___ :> -1 + r "
    "}";

    trig_expand_pythag = parse_expression(pythag_str);

    /*
     * TrigFactor rules.
     *
     * Split into three phases so we can control ordering precisely:
     *   to_sincos    -- rewrite Tan/Cot/Sec/Csc (and hyperbolic) as Sin/Cos
     *                   ratios before polynomial factoring,
     *   identities   -- Pythagorean collapse, reverse angle-addition, reverse
     *                   double-angle, and factored-form Pythagorean collapses
     *                   that arise from Factor,
     *   from_sincos  -- restore reciprocal trig on the output so Tan, Sec,
     *                   Tanh, Sech etc. are preserved in the factored form.
     */
    const char* to_sincos_str = "{ "
        "Tan[x_]  :> Sin[x] / Cos[x], "
        "Cot[x_]  :> Cos[x] / Sin[x], "
        "Sec[x_]  :> 1 / Cos[x], "
        "Csc[x_]  :> 1 / Sin[x], "
        "Tanh[x_] :> Sinh[x] / Cosh[x], "
        "Coth[x_] :> Cosh[x] / Sinh[x], "
        "Sech[x_] :> 1 / Cosh[x], "
        "Csch[x_] :> 1 / Sinh[x] "
    "}";
    trig_factor_to_sincos = parse_expression(to_sincos_str);

    /*
     * Identity rules. The tail patterns (r___ in Plus, r___ in Times) let each
     * rule fire inside a larger sum or product without requiring the sum /
     * product to match exactly. Orderless attributes on Plus and Times handle
     * the permutation search.
     *
     * Pythagorean collapses (circular Sin^2 + Cos^2 = 1, hyperbolic
     * Cosh^2 - Sinh^2 = 1) appear in both directions because Factor may
     * emerge with either sign depending on which way the leading coefficient
     * went.
     *
     * Factored forms:
     *   (Cos - Sin)(Cos + Sin) = Cos^2 - Sin^2   = Cos[2x]
     *   (Cosh - 1)(Cosh + 1)   = Cosh^2 - 1      = Sinh^2
     *   (Cosh - Sinh)(Cosh + Sinh) = Cosh^2 - Sinh^2 = 1
     *
     * Reverse angle-addition and reverse double-angle fire inside a larger
     * sum so that mixed expressions such as
     *   Sin[a]Cos[b] + Cos[a]Sin[b] + extra_term
     * collapse to Sin[a + b] + extra_term cleanly.
     */
    const char* identities_str = "{ "
        /* Pythagorean identities (Plus context, with optional coefficient). */
        "a_. Sin[x_]^2 + a_. Cos[x_]^2 + r___ :> a + r, "
        "a_. Cosh[x_]^2 - a_. Sinh[x_]^2 + r___ :> a + r, "
        "-a_. Cosh[x_]^2 + a_. Sinh[x_]^2 + r___ :> -a + r, "
        /* Pythagorean identities (Times/factored context). */
        "(Sin[x_]^2 + Cos[x_]^2) r___ :> 1 r, "
        "(Cosh[x_]^2 - Sinh[x_]^2) r___ :> 1 r, "
        "(-Cosh[x_]^2 + Sinh[x_]^2) r___ :> -1 r, "
        /* (Cosh + Sinh)(Cosh - Sinh) = 1 and variant with leading sign. */
        "(Cosh[x_] - Sinh[x_]) (Cosh[x_] + Sinh[x_]) r___ :> 1 r, "
        "(-Cosh[x_] + Sinh[x_]) (Cosh[x_] + Sinh[x_]) r___ :> -1 r, "
        /* (Cosh + 1)(Cosh - 1) = Sinh^2 -- appears from Factor[Cosh^2 - 1]. */
        "(-1 + Cosh[x_]) (1 + Cosh[x_]) r___ :> Sinh[x]^2 r, "
        "(1 - Cosh[x_]) (1 + Cosh[x_]) r___ :> -Sinh[x]^2 r, "
        /* (Cos - Sin)(Cos + Sin) = Cos[2x] -- from Factor[Cos^2 - Sin^2]. */
        "(Cos[x_] - Sin[x_]) (Cos[x_] + Sin[x_]) r___ :> Cos[2 x] r, "
        "(-Cos[x_] + Sin[x_]) (Cos[x_] + Sin[x_]) r___ :> -Cos[2 x] r, "
        /* Reverse angle-addition (circular). */
        "Sin[a_] Cos[b_] + Cos[a_] Sin[b_] + r___ :> Sin[a + b] + r, "
        "-Sin[a_] Cos[b_] - Cos[a_] Sin[b_] + r___ :> -Sin[a + b] + r, "
        "Sin[a_] Cos[b_] - Cos[a_] Sin[b_] + r___ :> Sin[a - b] + r, "
        "-Sin[a_] Cos[b_] + Cos[a_] Sin[b_] + r___ :> -Sin[a - b] + r, "
        "Cos[a_] Cos[b_] - Sin[a_] Sin[b_] + r___ :> Cos[a + b] + r, "
        "-Cos[a_] Cos[b_] + Sin[a_] Sin[b_] + r___ :> -Cos[a + b] + r, "
        "Cos[a_] Cos[b_] + Sin[a_] Sin[b_] + r___ :> Cos[a - b] + r, "
        /* Reverse angle-addition (hyperbolic). */
        "Sinh[a_] Cosh[b_] + Cosh[a_] Sinh[b_] + r___ :> Sinh[a + b] + r, "
        "-Sinh[a_] Cosh[b_] - Cosh[a_] Sinh[b_] + r___ :> -Sinh[a + b] + r, "
        "Sinh[a_] Cosh[b_] - Cosh[a_] Sinh[b_] + r___ :> Sinh[a - b] + r, "
        "-Sinh[a_] Cosh[b_] + Cosh[a_] Sinh[b_] + r___ :> -Sinh[a - b] + r, "
        "Cosh[a_] Cosh[b_] + Sinh[a_] Sinh[b_] + r___ :> Cosh[a + b] + r, "
        "-Cosh[a_] Cosh[b_] - Sinh[a_] Sinh[b_] + r___ :> -Cosh[a + b] + r, "
        "Cosh[a_] Cosh[b_] - Sinh[a_] Sinh[b_] + r___ :> Cosh[a - b] + r, "
        /* Reverse double-angle (Plus contexts). */
        "Cos[x_]^2 - Sin[x_]^2 + r___ :> Cos[2 x] + r, "
        "-Cos[x_]^2 + Sin[x_]^2 + r___ :> -Cos[2 x] + r, "
        "Cosh[x_]^2 + Sinh[x_]^2 + r___ :> Cosh[2 x] + r, "
        "-Cosh[x_]^2 - Sinh[x_]^2 + r___ :> -Cosh[2 x] + r, "
        /* Reverse double-angle (Times contexts -- 2 Sin Cos = Sin[2x]). */
        "2 Sin[x_] Cos[x_] r___ :> Sin[2 x] r, "
        "-2 Sin[x_] Cos[x_] r___ :> -Sin[2 x] r, "
        "2 Sinh[x_] Cosh[x_] r___ :> Sinh[2 x] r, "
        "-2 Sinh[x_] Cosh[x_] r___ :> -Sinh[2 x] r, "
        /* Linear-combination factoring: a Sin[x] + b Cos[x] -> Sqrt[a^2+b^2]
         * Sin[x + ArcTan[a, b]]. Only fires when both coefficients are
         * numeric AND real; for complex coefficients Sqrt[a^2+b^2] can
         * collapse to zero (e.g. a=I, b=1 gives Sqrt[I^2+1] = Sqrt[0] = 0,
         * which would falsely rewrite Cos[x] + I Sin[x] = E^(I x) to 0).
         * For symbolic coefficients the result would be more complex than
         * the input. The 2-arg ArcTan handles quadrant correctly. */
        "a_. Sin[x_] + b_. Cos[x_] + r___ /; (NumberQ[a] && NumberQ[b] && Im[a] === 0 && Im[b] === 0) "
            ":> Sqrt[a^2 + b^2] Sin[x + ArcTan[a, b]] + r "
    "}";
    trig_factor_identities = parse_expression(identities_str);

    /*
     * from_sincos: restore Tan/Sec/Cot/Csc (and hyperbolic analogs) from the
     * Sin/Cos ratio form. The higher-power rules (Sin^n Cos^-n) must appear
     * before the degree-1 rules so ReplaceRepeated chooses the most specific
     * match first.
     */
    const char* from_sincos_str = "{ "
        /* Sin/Cos ratio -> Tan/Cot (must precede the bare reciprocal rules
         * below so a ratio never matches as a bare Sec/Csc factor). */
        "Sin[x_]^2 Cos[x_]^(-2) r___ :> Tan[x]^2 r, "
        "Cos[x_]^2 Sin[x_]^(-2) r___ :> Cot[x]^2 r, "
        "Sin[x_]   Cos[x_]^(-1) r___ :> Tan[x] r, "
        "Cos[x_]   Sin[x_]^(-1) r___ :> Cot[x] r, "
        "Sinh[x_]^2 Cosh[x_]^(-2) r___ :> Tanh[x]^2 r, "
        "Cosh[x_]^2 Sinh[x_]^(-2) r___ :> Coth[x]^2 r, "
        "Sinh[x_]   Cosh[x_]^(-1) r___ :> Tanh[x] r, "
        "Cosh[x_]   Sinh[x_]^(-1) r___ :> Coth[x] r, "
        /* Bare reciprocal -> Sec/Csc. The r___ tail is required so factors
         * can exist alongside without gating the rewrite. */
        "Cos[x_]^(-2) r___ :> Sec[x]^2 r, "
        "Sin[x_]^(-2) r___ :> Csc[x]^2 r, "
        "Cos[x_]^(-1) r___ :> Sec[x] r, "
        "Sin[x_]^(-1) r___ :> Csc[x] r, "
        "Cosh[x_]^(-2) r___ :> Sech[x]^2 r, "
        "Sinh[x_]^(-2) r___ :> Csch[x]^2 r, "
        "Cosh[x_]^(-1) r___ :> Sech[x] r, "
        "Sinh[x_]^(-1) r___ :> Csch[x] r "
    "}";
    trig_factor_from_sincos = parse_expression(from_sincos_str);

    /*
     * TrigReduce rules.
     *
     * Two parallel directions of trig identity:
     *   - Power reduction:   Sin[x]^n / Cos[x]^n / Sinh[x]^n / Cosh[x]^n
     *     for integer n >= 2 are each pulled down by two powers, with
     *     the (1 +/- Cos[2 x])/2 (resp. hyperbolic) factor surfacing.
     *     Repeated application drops Sin[x]^4 -> Sin[x]^2 (1 - Cos[2 x])/2
     *     -> ((1 - Cos[2 x])/2)^2, after which the Expand step in the
     *     pipeline distributes the square and the rule fires again on
     *     any newly-exposed Cos[2 x]^2 term.
     *   - Product-to-sum: pairs of single-argument trig calls multiplied
     *     together collapse into single trig calls of compound arguments.
     *
     * Mathilda's Times is Orderless, so the matcher commutes the factors
     * in `Sin[a_] Cos[b_]` to align with whichever order the canonical
     * sort produced. Only one direction of each pair is needed -- the
     * commuted form is delivered by the matcher.
     */
    /* The Expand wrappers around the constructed arguments are needed
     * because Mathilda does not auto-distribute Times across Plus
     * (e.g. 2 (x + y) does not flatten to 2 x + 2 y, and (a + b) - (a - b)
     * stays as Plus[a, b, Times[-1, Plus[a, -b]]] rather than reducing to
     * 2 b). Without the Expand, the rule body would leave canonical-but-
     * unsimplified arguments inside the Sin/Cos/Sinh/Cosh wrappers, and
     * later Together/collapse passes could not see the actual scalar
     * (e.g. Sin[2 y] hidden inside Sin[(x + y) - (x - y)]). Expand inside
     * the rule body is evaluated as part of the substitution, so the
     * argument is canonicalised before the surrounding trig head sees
     * it. */
    const char* reduce_rules_str = "{ "
        "Sin[x_]^n_Integer /; n >= 2 :> Sin[x]^(n - 2) (1 - Cos[Expand[2 x]]) / 2, "
        "Cos[x_]^n_Integer /; n >= 2 :> Cos[x]^(n - 2) (1 + Cos[Expand[2 x]]) / 2, "
        "Sinh[x_]^n_Integer /; n >= 2 :> Sinh[x]^(n - 2) (-1 + Cosh[Expand[2 x]]) / 2, "
        "Cosh[x_]^n_Integer /; n >= 2 :> Cosh[x]^(n - 2) (1 + Cosh[Expand[2 x]]) / 2, "
        "Sin[a_] Cos[b_] r___ :> r (Sin[Expand[a + b]] + Sin[Expand[a - b]]) / 2, "
        "Sin[a_] Sin[b_] r___ :> r (Cos[Expand[a - b]] - Cos[Expand[a + b]]) / 2, "
        "Cos[a_] Cos[b_] r___ :> r (Cos[Expand[a + b]] + Cos[Expand[a - b]]) / 2, "
        "Sinh[a_] Cosh[b_] r___ :> r (Sinh[Expand[a + b]] + Sinh[Expand[a - b]]) / 2, "
        "Sinh[a_] Sinh[b_] r___ :> r (Cosh[Expand[a + b]] - Cosh[Expand[a - b]]) / 2, "
        "Cosh[a_] Cosh[b_] r___ :> r (Cosh[Expand[a + b]] + Cosh[Expand[a - b]]) / 2 "
    "}";
    trig_reduce_rules = parse_expression(reduce_rules_str);

    /*
     * Angle-addition collapse rules.
     *
     * Run after Together so a numerator of the form
     *   Sin[a] Cos[b] + Cos[a] Sin[b]
     * (often produced by Tan[a] + Tan[b] passing through the to-sin/cos
     * rewrite then Together) collapses to Sin[a + b]. Each rule has
     * sign variants because the sort that Mathilda applies inside the
     * Plus may surface either coefficient orientation.
     *
     * Trailing `r___` lets the rule fire inside a larger sum without
     * requiring an exact two-term match.
     */
    /*
     * The final four cancellation rules pull Sin[-X] = -Sin[X] (and the
     * Sinh analogue) and the even-function combinations Cos[-X] = Cos[X]
     * (and Cosh) out of sum-of-trig terms. Mathilda canonicalises
     * Sin[Times[-1, x]] but does NOT canonicalise Sin[Plus[Times[-1, a],
     * b]] (i.e. Sin[-a + b] is left as-is rather than simplifying to
     * -Sin[a - b]). After applying the product-to-sum rules to a Sin[a]
     * Cos[b] + Cos[a] Sin[b] form whose factors disagree on which
     * variable lands in the Sin slot, the matcher binds the rule's
     * pattern variables in opposite directions for the two terms,
     * producing Sin[a - b] + Sin[b - a] which mathematically equals zero
     * but is not recognised by the auto-evaluator. The condition
     * `SameQ[Expand[a + b], 0]` confirms the two argument expressions
     * are negatives of each other; the rules then combine the
     * coefficients with the appropriate sign for each parity (Sin/Sinh
     * are odd; Cos/Cosh are even).
     */
    const char* collapse_str = "{ "
        /* Circular angle-addition collapses, coefficient-aware.  c_. on
         * both sides defaults to 1, so the four rules subsume the eight
         * literal-coefficient cases (including the +1 / -1 split that
         * was previously enumerated).  Forcing the same `c` on both
         * Plus children is what lets a coefficient-bearing numerator
         * like
         *   k Cos[a] Cos[b] - k Sin[a] Sin[b]   ->  k Cos[a + b]
         * collapse -- the literal-1 form was the only case the prior
         * eight rules handled. */
        "c_. Sin[a_] Cos[b_] + c_. Cos[a_] Sin[b_] + r___ :> c Sin[Expand[a + b]] + r, "
        "c_. Sin[a_] Cos[b_] - c_. Cos[a_] Sin[b_] + r___ :> c Sin[Expand[a - b]] + r, "
        "c_. Cos[a_] Cos[b_] - c_. Sin[a_] Sin[b_] + r___ :> c Cos[Expand[a + b]] + r, "
        "c_. Cos[a_] Cos[b_] + c_. Sin[a_] Sin[b_] + r___ :> c Cos[Expand[a - b]] + r, "
        /* Hyperbolic angle-addition collapses, coefficient-aware. */
        "c_. Sinh[a_] Cosh[b_] + c_. Cosh[a_] Sinh[b_] + r___ :> c Sinh[Expand[a + b]] + r, "
        "c_. Sinh[a_] Cosh[b_] - c_. Cosh[a_] Sinh[b_] + r___ :> c Sinh[Expand[a - b]] + r, "
        "c_. Cosh[a_] Cosh[b_] + c_. Sinh[a_] Sinh[b_] + r___ :> c Cosh[Expand[a + b]] + r, "
        "c_. Cosh[a_] Cosh[b_] - c_. Sinh[a_] Sinh[b_] + r___ :> c Cosh[Expand[a - b]] + r, "
        /* Negative-argument cancellations */
        "c1_. Sin[a_] + c2_. Sin[b_] + r___ /; SameQ[Expand[a + b], 0] :> (c1 - c2) Sin[a] + r, "
        "c1_. Cos[a_] + c2_. Cos[b_] + r___ /; SameQ[Expand[a + b], 0] :> (c1 + c2) Cos[a] + r, "
        "c1_. Sinh[a_] + c2_. Sinh[b_] + r___ /; SameQ[Expand[a + b], 0] :> (c1 - c2) Sinh[a] + r, "
        "c1_. Cosh[a_] + c2_. Cosh[b_] + r___ /; SameQ[Expand[a + b], 0] :> (c1 + c2) Cosh[a] + r "
    "}";
    trig_reduce_collapse = parse_expression(collapse_str);

    symtab_add_builtin("ExpToTrig", builtin_exptotrig);
    symtab_get_def("ExpToTrig")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("TrigToExp", builtin_trigtoexp);
    symtab_get_def("TrigToExp")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("TrigExpand", builtin_trigexpand);
    symtab_get_def("TrigExpand")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("TrigFactor", builtin_trigfactor);
    symtab_get_def("TrigFactor")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("TrigReduce", builtin_trigreduce);
    symtab_get_def("TrigReduce")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);
}
