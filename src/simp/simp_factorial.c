#include "simp.h"
#include "simp_internal.h"
#include "arithmetic.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expand.h"
#include "facpoly.h"
#include "numeric.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "expr.h"
#include "rationalize.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "trigrat.h"
#include "qa.h"
#include "qafactor.h"
#include "simp_log.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif


/* ----------------------------------------------------------------------- */
/* simp_factorial -- general factorial simplification                       */
/* ----------------------------------------------------------------------- */

/*
 * Simplify expressions that contain Factorial[...] atoms via a single
 * principled four-step procedure rather than a per-pattern table.
 *
 *   Step A. Decompose every Factorial[arg] as arg = sym + c, where c is the
 *           integer addend in arg's Plus form (0 if arg has no integer
 *           term) and sym is the rest. Run Expand on the arg first so that
 *           shapes like 2(n-1) become 2n-2 (sym = 2n, c = -2).
 *   Step B. Group factorials whose sym parts are structurally equal. The
 *           group base is b = sym + min(c_i); each member's offset is
 *           k_i = c_i - min(c_i) >= 0.
 *   Step C. Rewrite each Factorial[b + k_i] as
 *               Factorial[b] * (b+1) * (b+2) * ... * (b+k_i).
 *           (k_i = 0 leaves the call unchanged.) Built as a literal Times
 *           product, then evaluated; the unifying Factorial[b] factor lets
 *           algebraic transforms see the cancellation/absorption.
 *   Step D. Run Together -> Cancel -> Expand on the rewritten form so that
 *           common Factorial[b] factors cancel and additive collapses fire
 *           ((n+1)! - n n! -> n!,   1/n! - 1/(n+1)! -> n/(n+1)!).
 *   Step E. Re-fold: walk every Times node holding a Factorial[b] factor
 *           and absorb consecutive (b+1)(b+2)...(b+k) cofactors back into
 *           Factorial[b+k]. A small gap-fill (budget 2) allows shapes like
 *           Factorial[b]*(b+1)(b+3)  ->  Factorial[b+3] / (b+2),
 *           which closes (n^2 - 1)*(n-2)! -> (n+1)!/n.
 *
 * The whole transform is gated on `contains_factorial(e)` and the result
 * is propagated as a Simplify seed: only the strictly-shorter form
 * survives the standard SimplifyCount tiebreak.
 */

bool contains_factorial(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (simp_eq_head_sym(e, "Factorial") && e->data.function.arg_count == 1) {
        return true;
    }
    if (contains_factorial(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_factorial(e->data.function.args[i])) return true;
    }
    return false;
}

/* Decompose `arg` as sym + c with c an int64 offset. Returns true on
 * success. *sym_out is a freshly-allocated Expr* (caller must free).
 * Conservative: any BigInt or Real addend aborts the decomposition (we
 * cannot soundly fold a factorial with a non-int64 offset under our
 * Pochhammer-style expansion). */
static bool simp_fact_decompose(const Expr* arg, Expr** sym_out,
                                int64_t* c_out) {
    *c_out = 0;
    *sym_out = NULL;
    if (arg->type == EXPR_INTEGER) {
        *c_out = arg->data.integer;
        *sym_out = expr_new_integer(0);
        return true;
    }
    if (arg->type == EXPR_BIGINT || arg->type == EXPR_REAL) return false;
    if (simp_eq_head_sym(arg, "Plus")) {
        size_t n = arg->data.function.arg_count;
        int64_t total = 0;
        Expr** rest = (Expr**)calloc(n, sizeof(Expr*));
        size_t rest_count = 0;
        for (size_t i = 0; i < n; i++) {
            const Expr* t = arg->data.function.args[i];
            if (t->type == EXPR_INTEGER) {
                int64_t v = t->data.integer;
                int64_t nt;
                if (__builtin_add_overflow(total, v, &nt)) {
                    for (size_t j = 0; j < rest_count; j++) expr_free(rest[j]);
                    free(rest);
                    return false;
                }
                total = nt;
            } else if (t->type == EXPR_BIGINT || t->type == EXPR_REAL) {
                for (size_t j = 0; j < rest_count; j++) expr_free(rest[j]);
                free(rest);
                return false;
            } else {
                rest[rest_count++] = expr_copy((Expr*)t);
            }
        }
        *c_out = total;
        if (rest_count == 0) {
            *sym_out = expr_new_integer(0);
        } else if (rest_count == 1) {
            *sym_out = rest[0];
        } else {
            *sym_out = expr_new_function(expr_new_symbol(SYM_Plus), rest,
                                         rest_count);
        }
        free(rest);
        return true;
    }
    /* Atom or non-Plus: c = 0, sym = arg. */
    *sym_out = expr_copy((Expr*)arg);
    return true;
}

/* Build the canonical expression sym + c (returns owned Expr*; consumes
 * `sym` ownership). Special-cases sym == 0 (returns just c) and c == 0
 * (returns just sym). */
static Expr* simp_fact_make_offset(Expr* sym, int64_t c) {
    bool sym_is_zero = (sym->type == EXPR_INTEGER && sym->data.integer == 0);
    if (sym_is_zero) {
        expr_free(sym);
        return expr_new_integer(c);
    }
    if (c == 0) return sym;
    Expr* args[2] = { sym, expr_new_integer(c) };
    return expr_new_function(expr_new_symbol(SYM_Plus), args, 2);
}

/* Group registry. Each entry maps a symbolic part to its minimum integer
 * offset across the input. Linear-scan find by structural equality. */
typedef struct {
    Expr*   sym;       /* owned */
    int64_t min_c;
    int64_t max_c;
    size_t  count;
} SimpFactGroup;

typedef struct {
    SimpFactGroup* items;
    size_t count;
    size_t cap;
} SimpFactGroupSet;

static void simp_fgs_init(SimpFactGroupSet* s) {
    s->items = NULL; s->count = 0; s->cap = 0;
}
static void simp_fgs_free(SimpFactGroupSet* s) {
    for (size_t i = 0; i < s->count; i++) expr_free(s->items[i].sym);
    free(s->items);
    s->items = NULL; s->count = 0; s->cap = 0;
}
static SimpFactGroup* simp_fgs_find(SimpFactGroupSet* s, const Expr* sym) {
    for (size_t i = 0; i < s->count; i++) {
        if (expr_eq(s->items[i].sym, sym)) return &s->items[i];
    }
    return NULL;
}
static void simp_fgs_record(SimpFactGroupSet* s, const Expr* sym, int64_t c) {
    SimpFactGroup* g = simp_fgs_find(s, sym);
    if (g) {
        if (c < g->min_c) g->min_c = c;
        if (c > g->max_c) g->max_c = c;
        g->count++;
        return;
    }
    if (s->count >= s->cap) {
        s->cap = s->cap ? s->cap * 2 : 4;
        s->items = (SimpFactGroup*)realloc(s->items, s->cap * sizeof(*s->items));
    }
    s->items[s->count].sym   = expr_copy((Expr*)sym);
    s->items[s->count].min_c = c;
    s->items[s->count].max_c = c;
    s->items[s->count].count = 1;
    s->count++;
}

/* Walk the expression collecting (sym, c) for every Factorial[arg]. Each
 * Factorial argument is first run through Expand so shapes like 2(n-1)
 * decompose under the expanded form 2n-2. */
static void simp_fact_gather(const Expr* e, SimpFactGroupSet* groups) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (simp_eq_head_sym(e, "Factorial") && e->data.function.arg_count == 1) {
        Expr* arg_exp = call_unary_copy("Expand", e->data.function.args[0]);
        if (arg_exp) {
            Expr* sym = NULL; int64_t c = 0;
            if (simp_fact_decompose(arg_exp, &sym, &c)) {
                simp_fgs_record(groups, sym, c);
                expr_free(sym);
            }
            expr_free(arg_exp);
        }
    }
    simp_fact_gather(e->data.function.head, groups);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        simp_fact_gather(e->data.function.args[i], groups);
    }
}

/* Build Factorial[b] * (b+1) * (b+2) * ... * (b + offset) for offset >= 0.
 * `b_sym_template` is borrowed -- the helper deep-copies as needed.
 * `b_const` is the integer addend that combines with `b_sym_template`
 * to form b. */
static Expr* simp_fact_pochhammer_expansion(const Expr* b_sym_template,
                                            int64_t b_const,
                                            int64_t offset) {
    /* Pochhammer guard: keep at most 32 explicit multiplied factors so a
     * pathological Factorial[n + 1000000] does not blow up the rewrite.
     * Beyond this we leave the call alone. */
    if (offset < 0 || offset > 32) return NULL;

    Expr* base = simp_fact_make_offset(expr_copy((Expr*)b_sym_template),
                                       b_const);
    Expr* fact_args[1] = { base };
    Expr* fact_b = expr_new_function(expr_new_symbol("Factorial"), fact_args, 1);
    if (offset == 0) return fact_b;

    /* Times[Factorial[b], (b+1), (b+2), ..., (b+offset)] */
    size_t n_args = (size_t)offset + 1;
    Expr** args = (Expr**)calloc(n_args, sizeof(Expr*));
    args[0] = fact_b;
    for (int64_t k = 1; k <= offset; k++) {
        args[(size_t)k] = simp_fact_make_offset(
            expr_copy((Expr*)b_sym_template), b_const + k);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_Times), args, n_args);
    free(args);
    return out;
}

/* Recursively rewrite each Factorial[arg] in `e` per the group base. The
 * returned tree is a fresh allocation. */
static Expr* simp_fact_rewrite(const Expr* e, const SimpFactGroupSet* groups) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    if (simp_eq_head_sym(e, "Factorial") && e->data.function.arg_count == 1) {
        Expr* arg_exp = call_unary_copy("Expand", e->data.function.args[0]);
        if (arg_exp) {
            Expr* sym = NULL; int64_t c = 0;
            if (simp_fact_decompose(arg_exp, &sym, &c)) {
                /* find the group for this sym */
                for (size_t i = 0; i < groups->count; i++) {
                    if (expr_eq(groups->items[i].sym, sym)) {
                        int64_t mc = groups->items[i].min_c;
                        Expr* exp = simp_fact_pochhammer_expansion(
                            sym, mc, c - mc);
                        if (exp) {
                            expr_free(sym);
                            expr_free(arg_exp);
                            return exp;
                        }
                        break;
                    }
                }
                expr_free(sym);
            }
            expr_free(arg_exp);
        }
        /* fall through: rewrite arg recursively (not a known group) */
    }

    /* Generic function: rewrite head + args. */
    Expr* new_head = simp_fact_rewrite(e->data.function.head, groups);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = simp_fact_rewrite(e->data.function.args[i], groups);
    }
    Expr* out = expr_new_function(new_head, new_args, n);
    free(new_args);
    return out;
}

/* Re-fold helper: given a list of Times args, look at the first Factorial
 * factor and try to absorb consecutive (base+j) cofactors. Returns a
 * new Times-content list (possibly with a Power[..., -1] denominator
 * for gap-fill terms) on success, or NULL when no fold fires.
 *
 * Inputs:
 *   args        -- borrowed Times children (caller still owns).
 *   n           -- count of `args`.
 *   out_args    -- on success, *out_args is a fresh Expr** of *out_n
 *                  freshly-allocated children that the caller owns.
 *   out_n       -- count of children in *out_args.
 */
static bool simp_fact_refold_times(Expr** args, size_t n,
                                   Expr*** out_args, size_t* out_n) {
    if (n < 2) return false;
    /* Find the first Factorial[b] direct child. */
    size_t fact_idx = (size_t)-1;
    Expr* base_sym = NULL;
    int64_t base_c = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* a = args[i];
        if (simp_eq_head_sym(a, "Factorial") &&
            a->data.function.arg_count == 1) {
            Expr* exp = call_unary_copy("Expand", a->data.function.args[0]);
            if (exp) {
                if (simp_fact_decompose(exp, &base_sym, &base_c)) {
                    fact_idx = i;
                    expr_free(exp);
                    break;
                }
                expr_free(exp);
            }
        }
    }
    if (fact_idx == (size_t)-1) return false;

    /* Collect candidate absorbing factors: args of the form sym + c with
     * the SAME sym as `base_sym`. Record (arg_index, j = c - base_c). */
    typedef struct { size_t idx; int64_t j; } AbsCand;
    AbsCand* cands = (AbsCand*)calloc(n, sizeof(AbsCand));
    size_t ncand = 0;
    for (size_t i = 0; i < n; i++) {
        if (i == fact_idx) continue;
        Expr* a = args[i];
        Expr* sym2 = NULL; int64_t c2 = 0;
        if (!simp_fact_decompose(a, &sym2, &c2)) continue;
        if (expr_eq(sym2, base_sym)) {
            int64_t j = c2 - base_c;
            if (j >= 1) {
                cands[ncand].idx = i;
                cands[ncand].j   = j;
                ncand++;
            }
        }
        expr_free(sym2);
    }
    if (ncand == 0) {
        free(cands);
        expr_free(base_sym);
        return false;
    }

    /* Sort cands by j ascending. Small n -- a bubble pass is fine. */
    for (size_t i = 0; i < ncand; i++) {
        for (size_t j = i + 1; j < ncand; j++) {
            if (cands[j].j < cands[i].j) {
                AbsCand t = cands[i]; cands[i] = cands[j]; cands[j] = t;
            }
        }
    }
    /* Drop duplicate j-values (a Times node should not contain duplicate
     * (base+j) factors, but be defensive). */
    {
        size_t w = 0;
        for (size_t i = 0; i < ncand; i++) {
            if (w == 0 || cands[w-1].j != cands[i].j) {
                cands[w++] = cands[i];
            }
        }
        ncand = w;
    }

    /* Find the largest k such that:
     *   - {1, 2, ..., k} subset of cand_offsets union gap_set
     *   - |gap_set| <= GAP_BUDGET
     *   - in_count >= gap_count (else the fold strictly grows the form).
     * The bound k_max = largest cand offset (folding past it adds factors
     * we did not have).
     */
    const int GAP_BUDGET = 2;
    int64_t kmax = cands[ncand - 1].j;
    int64_t best_k = 0;
    int64_t best_gaps = 0;
    /* We iterate k from 1 to kmax. Maintain a pointer into cands for
     * "how many cand offsets are <= k". */
    size_t ci = 0;
    for (int64_t k = 1; k <= kmax; k++) {
        while (ci < ncand && cands[ci].j <= k) ci++;
        int64_t n_in = (int64_t)ci;
        int64_t n_gaps = k - n_in;
        if (n_gaps > GAP_BUDGET) break;
        if (n_in < n_gaps) continue;
        if (k > best_k) { best_k = k; best_gaps = n_gaps; }
    }
    if (best_k == 0) {
        free(cands);
        expr_free(base_sym);
        return false;
    }

    /* Identify which cand indices we are absorbing (j <= best_k). */
    bool* absorbed = (bool*)calloc(n, sizeof(bool));
    bool* present_offset = (bool*)calloc((size_t)best_k + 1, sizeof(bool));
    for (size_t i = 0; i < ncand; i++) {
        if (cands[i].j <= best_k) {
            absorbed[cands[i].idx] = true;
            present_offset[cands[i].j] = true;
        }
    }

    /* Build the new args list:
     *   - replace args[fact_idx] with Factorial[base_sym + base_c + best_k].
     *   - drop args[absorbed].
     *   - keep all other args unchanged.
     *   - append gap_factors as Power[(base_sym + base_c + m), -1] for each
     *     m in 1..best_k missing from present_offset.
     */
    size_t new_n = 0;
    for (size_t i = 0; i < n; i++) {
        if (i == fact_idx) { new_n++; continue; }
        if (!absorbed[i]) new_n++;
    }
    new_n += (size_t)best_gaps;
    Expr** built = (Expr**)calloc(new_n ? new_n : 1, sizeof(Expr*));
    size_t w = 0;
    for (size_t i = 0; i < n; i++) {
        if (i == fact_idx) {
            Expr* nb = simp_fact_make_offset(expr_copy(base_sym),
                                             base_c + best_k);
            Expr* fa[1] = { nb };
            built[w++] = expr_new_function(expr_new_symbol("Factorial"),
                                           fa, 1);
        } else if (!absorbed[i]) {
            built[w++] = expr_copy(args[i]);
        }
    }
    for (int64_t m = 1; m <= best_k; m++) {
        if (!present_offset[m]) {
            Expr* fac = simp_fact_make_offset(expr_copy(base_sym),
                                              base_c + m);
            Expr* pa[2] = { fac, expr_new_integer(-1) };
            built[w++] = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
        }
    }

    free(absorbed);
    free(present_offset);
    free(cands);
    expr_free(base_sym);
    *out_args = built;
    *out_n = new_n;
    return true;
}

/* Combine Times children that all carry exponent -1 into a single
 * Power[Times[..., -1]]. Mathilda's evaluator does NOT auto-coalesce
 * separate Power[a, -1] * Power[b, -1] into Power[Times[a, b], -1]
 * (Mathematica does, Mathilda doesn't), so the auto-cancelled output of
 * a factorial rewrite sits at SimplifyCount 12 even when the same
 * expression printed as `1/(a*b)` would score 9. Lifting all negative-
 * exponent factors into a shared denominator gets us back to the lower
 * canonical score so the FactorialRules seed wins the tiebreak.
 *
 * Conservative: only combines exponent -1 (the trivial "denominator"
 * case). Mixed-sign exponents are left untouched. */
static Expr* simp_fact_combine_inverses(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (!simp_eq_head_sym(e, "Times")) {
        /* Recurse into children. */
        Expr* new_head = simp_fact_combine_inverses(e->data.function.head);
        size_t n = e->data.function.arg_count;
        Expr** new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            new_args[i] = simp_fact_combine_inverses(e->data.function.args[i]);
        }
        Expr* out = expr_new_function(new_head, new_args, n);
        free(new_args);
        return out;
    }
    size_t n = e->data.function.arg_count;
    /* Recurse into each Times child first. */
    Expr** child = (Expr**)calloc(n, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        child[i] = simp_fact_combine_inverses(e->data.function.args[i]);
    }
    /* Partition children: numerator (no -1 power), denominator (Power[X, -1]). */
    Expr** num = (Expr**)calloc(n, sizeof(Expr*));
    Expr** den = (Expr**)calloc(n, sizeof(Expr*));
    size_t nn = 0, dn = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* a = child[i];
        if (simp_eq_head_sym(a, "Power") &&
            a->data.function.arg_count == 2 &&
            a->data.function.args[1]->type == EXPR_INTEGER &&
            a->data.function.args[1]->data.integer == -1) {
            den[dn++] = expr_copy(a->data.function.args[0]);
        } else {
            num[nn++] = expr_copy(a);
        }
    }
    /* If 0 or 1 denominator factors, nothing to combine. Restore original. */
    if (dn <= 1) {
        for (size_t i = 0; i < nn; i++) expr_free(num[i]);
        for (size_t i = 0; i < dn; i++) expr_free(den[i]);
        free(num); free(den);
        Expr* out = expr_new_function(
            expr_copy(e->data.function.head), child, n);
        /* expr_new_function takes ownership of `child` slot but we still
         * need to free our wrapper allocation. */
        free(child);
        return out;
    }
    /* Build combined denominator: Times[den[0], ..., den[dn-1]] */
    Expr* den_times;
    if (dn == 1) {
        den_times = den[0];
    } else {
        Expr** dargs = (Expr**)calloc(dn, sizeof(Expr*));
        for (size_t i = 0; i < dn; i++) dargs[i] = den[i];
        den_times = expr_new_function(expr_new_symbol(SYM_Times), dargs, dn);
        free(dargs);
    }
    Expr* den_pow_args[2] = { den_times, expr_new_integer(-1) };
    Expr* den_pow = expr_new_function(expr_new_symbol(SYM_Power),
                                      den_pow_args, 2);
    /* Build new Times: numerator factors + den_pow. */
    if (nn == 0) {
        for (size_t i = 0; i < n; i++) expr_free(child[i]);
        free(child); free(num); free(den);
        return den_pow;
    }
    Expr** new_args = (Expr**)calloc(nn + 1, sizeof(Expr*));
    for (size_t i = 0; i < nn; i++) new_args[i] = num[i];
    new_args[nn] = den_pow;
    Expr* out = expr_new_function(expr_new_symbol(SYM_Times), new_args, nn + 1);
    free(new_args);
    for (size_t i = 0; i < n; i++) expr_free(child[i]);
    free(child); free(num); free(den);
    return out;
}

/* Re-fold tree walk. Bottom-up rebuild; at each Times node, attempt the
 * absorption pass repeatedly until it stops firing. */
static Expr* simp_fact_refold(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Recurse first. */
    Expr* new_head = simp_fact_refold(e->data.function.head);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = simp_fact_refold(e->data.function.args[i]);
    }

    if (simp_eq_head_sym(e, "Times")) {
        /* Iterate refold until no more folds. Bounded by initial n: each
         * fold reduces the cofactor count by at least 1 (one absorbed
         * arg replaced; possibly +gaps but those are appended once). */
        for (size_t guard = 0; guard < 16; guard++) {
            Expr** out = NULL; size_t out_n = 0;
            if (!simp_fact_refold_times(new_args, n, &out, &out_n)) break;
            for (size_t i = 0; i < n; i++) expr_free(new_args[i]);
            free(new_args);
            new_args = out;
            n = out_n;
        }
    }

    Expr* out = expr_new_function(new_head, new_args, n);
    free(new_args);
    return out;
}

/* Recognize the canonical (2 v)! / (2^v v!) shape and rewrite to
 * Factorial2[2 v - 1]. The check is purely structural: input must be a
 * Times whose direct children include
 *    Factorial[Times[2, v]]                    (numerator factorial)
 *    Power[Factorial[v], -1]                   (denominator factorial)
 *    Power[2, Times[-1, v]]                    (1 / 2^v)
 * with the same v expression appearing in all three. Other Times
 * children are kept as a residual cofactor on the result. */
static Expr* simp_fact_double_factorial(const Expr* e) {
    if (!simp_eq_head_sym(e, "Times")) return NULL;
    size_t n = e->data.function.arg_count;
    if (n < 3) return NULL;

    int idx_2v_fact = -1;       /* Factorial[2 v]            */
    int idx_v_fact_inv = -1;    /* Power[Factorial[v], -1]   */
    int idx_2_pow_neg = -1;     /* Power[2, -v]              */
    Expr* v_expr = NULL;

    for (size_t i = 0; i < n; i++) {
        Expr* a = e->data.function.args[i];
        /* Factorial[Times[2, v]] */
        if (idx_2v_fact == -1 &&
            simp_eq_head_sym(a, "Factorial") &&
            a->data.function.arg_count == 1) {
            Expr* arg = a->data.function.args[0];
            Expr* arg_exp = call_unary_copy("Expand", arg);
            if (arg_exp && simp_eq_head_sym(arg_exp, "Times") &&
                arg_exp->data.function.arg_count == 2) {
                Expr* a0 = arg_exp->data.function.args[0];
                Expr* a1 = arg_exp->data.function.args[1];
                if (a0->type == EXPR_INTEGER && a0->data.integer == 2) {
                    if (!v_expr || expr_eq(v_expr, a1)) {
                        if (!v_expr) v_expr = expr_copy(a1);
                        idx_2v_fact = (int)i;
                        expr_free(arg_exp);
                        continue;
                    }
                }
            }
            if (arg_exp) expr_free(arg_exp);
        }
        /* Power[Factorial[v], -1] */
        if (idx_v_fact_inv == -1 &&
            simp_eq_head_sym(a, "Power") &&
            a->data.function.arg_count == 2 &&
            a->data.function.args[1]->type == EXPR_INTEGER &&
            a->data.function.args[1]->data.integer == -1 &&
            simp_eq_head_sym(a->data.function.args[0], "Factorial") &&
            a->data.function.args[0]->data.function.arg_count == 1) {
            Expr* v = a->data.function.args[0]->data.function.args[0];
            if (!v_expr || expr_eq(v_expr, v)) {
                if (!v_expr) v_expr = expr_copy(v);
                idx_v_fact_inv = (int)i;
                continue;
            }
        }
        /* Power[2, Times[-1, v]] or Power[2, -v] depending on shape. */
        if (idx_2_pow_neg == -1 &&
            simp_eq_head_sym(a, "Power") &&
            a->data.function.arg_count == 2 &&
            a->data.function.args[0]->type == EXPR_INTEGER &&
            a->data.function.args[0]->data.integer == 2) {
            Expr* exp = a->data.function.args[1];
            /* Build Times[-1, v_expr_candidate] and compare. We instead
             * compare via -exp == v: evaluate Times[-1, exp] and check. */
            Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(exp) };
            Expr* neg_call = expr_new_function(expr_new_symbol(SYM_Times),
                                               neg_args, 2);
            Expr* neg_eval = evaluate(neg_call);
            expr_free(neg_call);
            if (neg_eval) {
                if (!v_expr || expr_eq(v_expr, neg_eval)) {
                    if (!v_expr) v_expr = expr_copy(neg_eval);
                    idx_2_pow_neg = (int)i;
                    expr_free(neg_eval);
                    continue;
                }
                expr_free(neg_eval);
            }
        }
    }

    if (idx_2v_fact < 0 || idx_v_fact_inv < 0 || idx_2_pow_neg < 0 || !v_expr) {
        if (v_expr) expr_free(v_expr);
        return NULL;
    }

    /* Build Factorial2[2 v - 1]. */
    Expr* two_v_args[2] = { expr_new_integer(2), expr_copy(v_expr) };
    Expr* two_v = expr_new_function(expr_new_symbol(SYM_Times), two_v_args, 2);
    Expr* arg_args[2] = { two_v, expr_new_integer(-1) };
    Expr* arg_plus = expr_new_function(expr_new_symbol(SYM_Plus), arg_args, 2);
    Expr* fac2_args[1] = { arg_plus };
    Expr* fac2 = expr_new_function(expr_new_symbol("Factorial2"), fac2_args, 1);
    expr_free(v_expr);

    /* If there are residual cofactors, multiply them in. */
    size_t residue = 0;
    for (size_t i = 0; i < n; i++) {
        if ((int)i == idx_2v_fact || (int)i == idx_v_fact_inv ||
            (int)i == idx_2_pow_neg) continue;
        residue++;
    }
    if (residue == 0) return fac2;

    Expr** all = (Expr**)calloc(residue + 1, sizeof(Expr*));
    size_t w = 0;
    all[w++] = fac2;
    for (size_t i = 0; i < n; i++) {
        if ((int)i == idx_2v_fact || (int)i == idx_v_fact_inv ||
            (int)i == idx_2_pow_neg) continue;
        all[w++] = expr_copy(e->data.function.args[i]);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_Times), all, residue + 1);
    free(all);
    return out;
}

/* Walk the tree applying simp_fact_double_factorial at every Times node. */
static Expr* simp_fact_double_factorial_walk(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    Expr* new_head = simp_fact_double_factorial_walk(e->data.function.head);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = simp_fact_double_factorial_walk(e->data.function.args[i]);
    }
    Expr* out = expr_new_function(new_head, new_args, n);
    free(new_args);
    if (simp_eq_head_sym(out, "Times")) {
        Expr* d = simp_fact_double_factorial(out);
        if (d) { expr_free(out); out = d; }
    }
    return out;
}

/* Top-level transform. Returns NULL when the input contains no factorial
 * or when the rewrite produced no change. Otherwise returns the
 * candidate (a freshly-allocated Expr*; caller takes ownership). */
Expr* simp_factorial(const Expr* e) {
    if (!contains_factorial(e)) return NULL;
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    SimpFactGroupSet groups;
    simp_fgs_init(&groups);
    simp_fact_gather(e, &groups);

    /* If every group has exactly one offset (no shifting possible) AND the
     * tree contains only one factorial, the substitution step is a no-op.
     * Still proceed to the double-factorial pattern check below. */
    bool any_shift = false;
    for (size_t i = 0; i < groups.count; i++) {
        if (groups.items[i].max_c > groups.items[i].min_c) {
            any_shift = true;
            break;
        }
    }

    Expr* result = NULL;
    if (any_shift) {
        /* Step C — substitute. The rewrite emits Pochhammer products
         * literally (e.g. Factorial[n] -> Factorial[n-2]*(n-1)*n in a
         * group whose base is n-2) and the evaluator's auto-cancellation
         * through normal Times/Power semantics is sufficient to remove
         * any common Factorial[base] factor between numerator and
         * denominator. We deliberately do NOT run Together here because
         * Together would expand the polynomial denominator (turning
         * 1/(n*(n-1)) into 1/(n^2 - n) and Cancel cannot recover the
         * factored form afterwards because the original Power[Times,
         * -1] structure is gone). The auto-cancel already gives the
         * minimal canonical fraction.
         *
         * Step D' (the lighter version) -- run Cancel on the result so
         * that any residual rational form left by the evaluator is
         * reduced. Skip if the evaluator's own auto-cancel already
         * produced a form without nested Times/Power[Times,-1] noise. */
        Expr* rewritten = simp_fact_rewrite(e, &groups);
        /* Force a full evaluator pass so the auto-cancellation between
         * Factorial[base] in the numerator and the same factor in a
         * Pochhammer-expanded denominator collapses. evaluate() does NOT
         * consume its argument, so we still own `rewritten` and must free
         * it after taking the evaluated form. */
        Expr* evaluated = evaluate(rewritten);
        if (!evaluated) evaluated = expr_copy(rewritten);
        expr_free(rewritten);
        rewritten = NULL;

        /* Expand distributes products over sums, which is what makes
         *   Factorial[n]*(1+n) - n*Factorial[n] -> Factorial[n]
         *   (n+1)*n*Factorial[n-1]            -> n^2 Factorial[n-1] + ...
         * fire. Without it, the evaluator does not distribute (a + b)
         * across c, so Plus's Orderless/FLAT cancellation never sees
         * the like terms. Harmless on inputs that have nothing to
         * distribute (e.g. plain rational forms like 1/(n*(n-1))). */
        Expr* expanded = call_unary_copy("Expand", evaluated);
        if (!expanded) expanded = expr_copy(evaluated);

        /* Two paths feed the re-fold pass:
         *   A. expanded form        -- preserves a factored
         *                              Times[Power[a,-1], Power[b,-1]]
         *                              denominator so (n-2)!/n! lands
         *                              at 1/((n-1)*n) instead of being
         *                              re-expanded to 1/(n^2-n).
         *   B. Together'd form      -- combines additive fractions
         *                              over a common denominator,
         *                              needed for 1/n! - 1/(n+1)! to
         *                              fold to n/(n+1)!.
         * Both paths run Factor (to surface (b+j) linear factors in
         * the cofactor) and combine_inverses (to re-coalesce the
         * Power[a,-1]*Power[b,-1] factors that Mathilda's evaluator
         * leaves separated, which is what makes the combined Times
         * denominator visible to the re-fold walker). The score
         * tiebreak below picks whichever path lands at the lowest
         * SimplifyCount. */
        Expr* together = call_unary_copy("Together", expanded);

        Expr* prep_a = NULL;
        if (!has_non_integer_power(expanded)) {
            factor_memo_push(NULL);
            prep_a = call_unary_copy("Factor", expanded);
            factor_memo_pop();
        }
        if (!prep_a) prep_a = expr_copy(expanded);
        {
            Expr* coal = simp_fact_combine_inverses(prep_a);
            if (coal) { expr_free(prep_a); prep_a = coal; }
        }
        Expr* refold_a = simp_fact_refold(prep_a);

        /* Path B: Together expands the polynomial denominator (e.g.
         *   1/n! - 1/((n+1) n!)  -> n / (n! + n n!)
         * landing in a Plus inside Power[..., -1]). Factor pulls the
         * Plus apart back into a Times of Power[a, -1] factors:
         *   n / (n! + n n!) -> Times[n, Power[n!, -1], Power[1+n, -1]].
         * combine_inverses then coalesces those into a single
         * Power[Times[n!, 1+n], -1], and the re-fold walker descends
         * into that Times to fold Factorial[n] * (n+1) -> Factorial[n+1]. */
        Expr* prep_b = NULL;
        if (together) {
            if (!has_non_integer_power(together)) {
                /* Push a NULL memo to opt out of Factor's inside-Simplify
                 * variable-list narrowing. The narrowing collapses
                 * num/den variable scopes and prevents Factor from
                 * pulling Factorial[n] out of a denominator like
                 * Factorial[n] + n*Factorial[n] -> Factorial[n]*(1+n)
                 * (the polynomial-in-n viewer treats Factorial[n] as a
                 * coefficient that the inside_simplify path then
                 * refuses to factor). With separate variable lists the
                 * factorisation succeeds. We still consult / populate
                 * the outer memo via factor_memo_active() inside the
                 * helper -- only the gating in builtin_factor that
                 * branches on factor_memo_top() != NULL is what we
                 * silence here. */
                factor_memo_push(NULL);
                prep_b = call_unary_copy("Factor", together);
                factor_memo_pop();
            }
            if (!prep_b) prep_b = expr_copy(together);
            Expr* coal = simp_fact_combine_inverses(prep_b);
            if (coal) { expr_free(prep_b); prep_b = coal; }
        }
        Expr* refold_b = prep_b ? simp_fact_refold(prep_b) : NULL;

        Expr* refolded = refold_a;
        if (refold_b && simp_default_complexity(refold_b) <
                        simp_default_complexity(refold_a)) {
            expr_free(refold_a);
            refolded = refold_b;
        } else if (refold_b) {
            expr_free(refold_b);
        }
        if (prep_a) expr_free(prep_a);
        if (prep_b) expr_free(prep_b);
        if (together) expr_free(together);
        expr_free(expanded);

        /* Note on Cancel as an alternate: tempting, but Cancel/Together
         * on the post-refold form expand a Times[Power[a,-1], Power[b,-1]]
         * factored denominator into Power[Plus[a*b...], -1] (e.g.
         * 1/(n*(n-1)) -> 1/(n^2 - n)). The expanded form scores lower
         * under SimplifyCount but is the canonical "wrong" answer for
         * our purposes -- the user's expected forms keep the factored
         * denominator. We rely on the evaluator's own auto-cancel
         * (which produced `evaluated` from `rewritten`) to give the
         * canonical reduction; no further smoothing is wanted. */

        expr_free(evaluated);
        result = refolded;
    } else {
        /* No shift step; still try Factor + fold over the input (covers
         * n*(n-1)! -> n! and (n^2-1)*(n-2)! shapes where there is only
         * one factorial group but the cofactor needs factoring before
         * the (base+j) factors become visible). The NULL memo push
         * mirrors the any_shift path -- see comment there. */
        Expr* factored = NULL;
        if (!has_non_integer_power(e)) {
            factor_memo_push(NULL);
            factored = call_unary_copy("Factor", e);
            factor_memo_pop();
        }
        if (!factored) factored = expr_copy((Expr*)e);
        result = simp_fact_refold(factored);
        expr_free(factored);
    }

    /* Double-factorial pattern recognition. */
    {
        Expr* d = simp_fact_double_factorial_walk(result);
        if (d) { expr_free(result); result = d; }
    }

    /* Final canonicalization: coalesce all `Power[X, -1]` factors in
     * Times nodes into a single `Power[Times[..., -1]]`. Mathilda's
     * evaluator does not auto-perform this combine, which leaves the
     * factored result at a higher SimplifyCount than it deserves
     * (count 12 vs count 9 on `1/(n*(n-1))`). Coalescing brings the
     * score under the original input, so the FactorialRules seed
     * wins the round-loop tiebreak instead of being undone. */
    {
        Expr* c = simp_fact_combine_inverses(result);
        if (c) { expr_free(result); result = c; }
    }

    simp_fgs_free(&groups);
    if (result && expr_eq(result, e)) {
        if (dbg) simp_debug_log("FactorialRules", e, NULL,
                                simp_debug_elapsed_ms(t0));
        expr_free(result);
        return NULL;
    }
    if (dbg) simp_debug_log("FactorialRules", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

