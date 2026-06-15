/* Trig / hyperbolic ratio canonicalization for builtin_times.
 *
 * Each of the six standard trig functions can be expressed as a product of
 * integer powers of Sin and Cos:
 *
 *     Sin[x]  -> ( 1,  0)
 *     Cos[x]  -> ( 0,  1)
 *     Tan[x]  -> ( 1, -1)
 *     Cot[x]  -> (-1,  1)
 *     Sec[x]  -> ( 0, -1)
 *     Csc[x]  -> (-1,  0)
 *
 * After builtin_times has grouped factors by base, we walk the groups,
 * decompose each trig/hyperbolic factor into its (sin_exp, cos_exp) pair,
 * sum the pairs over factors that share the same argument, and re-emit the
 * shortest canonical form for the resulting (a, b). The hyperbolic family
 * {Sinh, Cosh, Tanh, Coth, Sech, Csch} is handled identically but separately.
 *
 * Only integer exponents on candidate factors trigger the rewrite. Anything
 * else (symbolic, real, rational, nested Power) is left as the evaluator
 * delivered it.
 */

#include "trig_canon.h"
#include "sym_names.h"
#include "eval.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* Reentrant suppress counter. See trig_canon.h for the rationale. Mathilda's
 * evaluator is single-threaded, so a plain global is fine. */
static int g_suppress = 0;

void trig_canon_suppress_inc(void) { g_suppress++; }
void trig_canon_suppress_dec(void) {
    if (g_suppress > 0) g_suppress--;
    /* Cache-coherence on the evaluator's last_evaluated_at timestamps:
     * subexpressions reduced while suppression was active have been
     * stamped at g_eval_clock, but their evaluation result is only valid
     * under the suppress flag.  Once the suppress region ends, the outer
     * evaluator's evaluate() short-circuit would otherwise return those
     * cached forms unchanged, even when trig_canon would now collapse
     * them (e.g. Cos[a] Sec[a] -> 1 inside a Plus child of a TrigReduce
     * result -- exactly the failure mode behind the Tan-addition
     * Simplify regression).  Bumping the clock on the outermost dec
     * invalidates the suppressed-mode cache so the next evaluate() pass
     * actually walks the tree and the trig canonicalization applies. */
    if (g_suppress == 0) eval_clock_bump();
}

/* Classify a base expression as a member of the trig or hyperbolic ratio
 * family. Returns 1 on hit, 0 otherwise. On hit, fills *arg_out (a borrowed
 * pointer into the base subtree), the unit Sin/Cos exponent contribution
 * (sin_per, cos_per), and *kind (0 = trig, 1 = hyperbolic). */
static int classify(Expr* base, Expr** arg_out,
                    int* sin_per, int* cos_per, int* kind) {
    if (!base || base->type != EXPR_FUNCTION) return 0;
    Expr* head = base->data.function.head;
    if (!head || head->type != EXPR_SYMBOL) return 0;
    if (base->data.function.arg_count != 1) return 0;
    const char* h = head->data.symbol;

    int sp = 0, cp = 0, k = 0;
    if      (h == SYM_Sin)  { sp =  1; cp =  0; k = 0; }
    else if (h == SYM_Cos)  { sp =  0; cp =  1; k = 0; }
    else if (h == SYM_Tan)  { sp =  1; cp = -1; k = 0; }
    else if (h == SYM_Cot)  { sp = -1; cp =  1; k = 0; }
    else if (h == SYM_Sec)  { sp =  0; cp = -1; k = 0; }
    else if (h == SYM_Csc)  { sp = -1; cp =  0; k = 0; }
    else if (h == SYM_Sinh) { sp =  1; cp =  0; k = 1; }
    else if (h == SYM_Cosh) { sp =  0; cp =  1; k = 1; }
    else if (h == SYM_Tanh) { sp =  1; cp = -1; k = 1; }
    else if (h == SYM_Coth) { sp = -1; cp =  1; k = 1; }
    else if (h == SYM_Sech) { sp =  0; cp = -1; k = 1; }
    else if (h == SYM_Csch) { sp = -1; cp =  0; k = 1; }
    else return 0;

    *arg_out = base->data.function.args[0];
    *sin_per = sp;
    *cos_per = cp;
    *kind = k;
    return 1;
}

static Expr* make_unary(const char* name, Expr* arg /* borrowed */) {
    Expr* args[1];
    args[0] = expr_copy(arg);
    return expr_new_function(expr_new_symbol(name), args, 1);
}

/* Bound for accepted exponents. Anything beyond this magnitude is assumed
 * uninteresting (and risks integer overflow when multiplied by sin_per/cos_per
 * accumulators). 2**20 is far beyond any exponent a human is likely to type. */
#define TRIG_CANON_EXP_BOUND 1048576

/* Map a (sin_exp, cos_exp) pair back to a sequence of (base, exponent) groups
 * using the shortest naming available. Writes at most 2 BasePower entries to
 * out[]. Returns the number written.
 *
 * The cases are exhaustive:
 *
 *   (0, 0)            -> 0 entries (factor of 1)
 *   (a, 0), a != 0    -> Sin^a (a>0) or Csc^|a| (a<0)
 *   (0, b), b != 0    -> Cos^b (b>0) or Sec^|b| (b<0)
 *   (a, b), a*b > 0   -> { Sin^a, Cos^b }   (or { Csc^|a|, Sec^|b| } if both negative)
 *   (a > 0, b < 0):     pull out Tan^min(a,|b|), leftover stays as Sin or Sec
 *   (a < 0, b > 0):     pull out Cot^min(|a|,b), leftover stays as Cos or Csc
 */
static size_t canonical_emit(int a, int b, Expr* arg, int kind, BasePower* out) {
    const char* sin_n = kind ? "Sinh" : "Sin";
    const char* cos_n = kind ? "Cosh" : "Cos";
    const char* tan_n = kind ? "Tanh" : "Tan";
    const char* cot_n = kind ? "Coth" : "Cot";
    const char* sec_n = kind ? "Sech" : "Sec";
    const char* csc_n = kind ? "Csch" : "Csc";

    size_t n = 0;

    if (a == 0 && b == 0) return 0;

    if (a == 0) {
        out[n].base     = make_unary(b > 0 ? cos_n : sec_n, arg);
        out[n].exponent = expr_new_integer(b > 0 ? b : -b);
        n++;
        return n;
    }
    if (b == 0) {
        out[n].base     = make_unary(a > 0 ? sin_n : csc_n, arg);
        out[n].exponent = expr_new_integer(a > 0 ? a : -a);
        n++;
        return n;
    }

    /* Both nonzero. */
    if (a > 0 && b > 0) {
        out[n].base = make_unary(sin_n, arg);
        out[n].exponent = expr_new_integer(a);
        n++;
        out[n].base = make_unary(cos_n, arg);
        out[n].exponent = expr_new_integer(b);
        n++;
        return n;
    }
    if (a < 0 && b < 0) {
        out[n].base = make_unary(csc_n, arg);
        out[n].exponent = expr_new_integer(-a);
        n++;
        out[n].base = make_unary(sec_n, arg);
        out[n].exponent = expr_new_integer(-b);
        n++;
        return n;
    }

    if (a > 0 && b < 0) {
        int s = a + b;  /* leftover: positive => Sin extras, negative => Sec extras */
        if (s > 0) {
            out[n].base = make_unary(sin_n, arg);
            out[n].exponent = expr_new_integer(s);
            n++;
            out[n].base = make_unary(tan_n, arg);
            out[n].exponent = expr_new_integer(-b);
            n++;
        } else if (s == 0) {
            out[n].base = make_unary(tan_n, arg);
            out[n].exponent = expr_new_integer(a);
            n++;
        } else {
            out[n].base = make_unary(tan_n, arg);
            out[n].exponent = expr_new_integer(a);
            n++;
            out[n].base = make_unary(sec_n, arg);
            out[n].exponent = expr_new_integer(-s);
            n++;
        }
        return n;
    }

    /* a < 0 && b > 0 */
    {
        int s = a + b;
        if (s > 0) {
            out[n].base = make_unary(cot_n, arg);
            out[n].exponent = expr_new_integer(-a);
            n++;
            out[n].base = make_unary(cos_n, arg);
            out[n].exponent = expr_new_integer(s);
            n++;
        } else if (s == 0) {
            out[n].base = make_unary(cot_n, arg);
            out[n].exponent = expr_new_integer(b);
            n++;
        } else {
            out[n].base = make_unary(cot_n, arg);
            out[n].exponent = expr_new_integer(b);
            n++;
            out[n].base = make_unary(csc_n, arg);
            out[n].exponent = expr_new_integer(-s);
            n++;
        }
        return n;
    }
}

typedef struct {
    Expr* arg;            /* borrowed: points into one of the contributing groups */
    int kind;             /* 0 = trig, 1 = hyperbolic */
    int a;                /* accumulated Sin exponent */
    int b;                /* accumulated Cos exponent */
    size_t* idxs;         /* indices into groups[] that fed this bucket */
    size_t idx_count;
    size_t idx_cap;
} Bucket;

static void bucket_add_idx(Bucket* b, size_t idx) {
    if (b->idx_count == b->idx_cap) {
        size_t nc = b->idx_cap ? b->idx_cap * 2 : 4;
        b->idxs = realloc(b->idxs, sizeof(size_t) * nc);
        b->idx_cap = nc;
    }
    b->idxs[b->idx_count++] = idx;
}

void trig_canon_groups(BasePower* groups, size_t* group_count) {
    if (g_suppress) return;
    size_t n = *group_count;
    if (n == 0) return;

    Bucket* buckets = NULL;
    size_t bucket_count = 0;
    size_t bucket_cap = 0;

    /* Pass 1: scan groups, classify trig/hyp factors, accumulate into buckets
     * keyed by (kind, arg). */
    for (size_t i = 0; i < n; i++) {
        Expr* arg = NULL;
        int sin_per = 0, cos_per = 0, kind = 0;
        if (!classify(groups[i].base, &arg, &sin_per, &cos_per, &kind)) continue;

        Expr* exp = groups[i].exponent;
        if (!exp || exp->type != EXPR_INTEGER) continue;
        int64_t ev = exp->data.integer;
        if (ev > TRIG_CANON_EXP_BOUND || ev < -TRIG_CANON_EXP_BOUND) continue;

        int da = sin_per * (int)ev;
        int db = cos_per * (int)ev;

        int found = -1;
        for (size_t k = 0; k < bucket_count; k++) {
            if (buckets[k].kind != kind) continue;
            if (!expr_eq(buckets[k].arg, arg)) continue;
            found = (int)k;
            break;
        }

        if (found < 0) {
            if (bucket_count == bucket_cap) {
                size_t nc = bucket_cap ? bucket_cap * 2 : 4;
                buckets = realloc(buckets, sizeof(Bucket) * nc);
                bucket_cap = nc;
            }
            Bucket* nb = &buckets[bucket_count++];
            nb->arg = arg;
            nb->kind = kind;
            nb->a = da;
            nb->b = db;
            nb->idxs = NULL;
            nb->idx_count = 0;
            nb->idx_cap = 0;
            bucket_add_idx(nb, i);
        } else {
            buckets[found].a += da;
            buckets[found].b += db;
            bucket_add_idx(&buckets[found], i);
        }
    }

    if (bucket_count == 0) {
        free(buckets);
        return;
    }

    /* Pass 2: for each bucket, decide whether a rewrite is worth doing.
     *
     * - A solo contributor (idx_count == 1) with a positive integer exponent
     *   is already canonical (Sin^k, Cos^k, Tan^k, ... for k > 0) and we leave
     *   it alone.
     *
     * - For multi-contributor buckets, after computing the canonical_emit
     *   output we compare it (as a set of (head, exp) pairs) against the
     *   contributors. If they match, the rewrite would just rebuild the
     *   same expression in possibly-different array order -- which would
     *   perturb the canonical sort order downstream (Times is Orderless
     *   and the args were sorted before builtin_times saw them; reordering
     *   here means callers see the unsorted order). Skipping pure no-ops
     *   keeps the output stable.
     *
     * - Anything else (real merge, sign flip to Tan/Sec/etc., or full
     *   cancellation) is dropped + re-emitted. */
    bool* drop = calloc(n, sizeof(bool));
    BasePower* fresh = NULL;
    size_t fresh_count = 0;
    size_t fresh_cap = 0;

    for (size_t k = 0; k < bucket_count; k++) {
        Bucket* bk = &buckets[k];

        if (bk->idx_count == 1) {
            int64_t ev = groups[bk->idxs[0]].exponent->data.integer;
            if (ev > 0) continue;  /* already canonical */
        }

        BasePower out[2];
        size_t emitted = canonical_emit(bk->a, bk->b, bk->arg, bk->kind, out);

        /* No-op detection: if the emit produces the same multiset of
         * (head_symbol, exponent) pairs as the contributors, leave them
         * alone to preserve their existing positions in groups[]. */
        if (emitted == bk->idx_count) {
            bool matched_all = true;
            bool used[2] = { false, false };
            for (size_t j = 0; j < bk->idx_count; j++) {
                Expr* cb = groups[bk->idxs[j]].base;
                Expr* ce = groups[bk->idxs[j]].exponent;
                bool found_match = false;
                for (size_t m = 0; m < emitted; m++) {
                    if (used[m]) continue;
                    if (!expr_eq(cb, out[m].base)) continue;
                    if (!expr_eq(ce, out[m].exponent)) continue;
                    used[m] = true;
                    found_match = true;
                    break;
                }
                if (!found_match) { matched_all = false; break; }
            }
            if (matched_all) {
                for (size_t m = 0; m < emitted; m++) {
                    expr_free(out[m].base);
                    expr_free(out[m].exponent);
                }
                continue;
            }
        }

        for (size_t j = 0; j < bk->idx_count; j++) drop[bk->idxs[j]] = true;

        for (size_t j = 0; j < emitted; j++) {
            if (fresh_count == fresh_cap) {
                size_t nc = fresh_cap ? fresh_cap * 2 : 4;
                fresh = realloc(fresh, sizeof(BasePower) * nc);
                fresh_cap = nc;
            }
            fresh[fresh_count++] = out[j];
        }
    }

    /* Pass 3: compact groups[] to remove dropped entries, then append the
     * freshly-emitted ones. The total never exceeds the original count
     * (see the header comment for why), so writing back into groups[] is
     * safe without reallocation. */
    size_t out_idx = 0;
    for (size_t i = 0; i < n; i++) {
        if (drop[i]) {
            expr_free(groups[i].base);
            expr_free(groups[i].exponent);
            continue;
        }
        if (out_idx != i) groups[out_idx] = groups[i];
        out_idx++;
    }
    for (size_t j = 0; j < fresh_count; j++) {
        groups[out_idx++] = fresh[j];
    }
    *group_count = out_idx;

    free(drop);
    free(fresh);
    for (size_t k = 0; k < bucket_count; k++) free(buckets[k].idxs);
    free(buckets);
}

Expr* trig_canon_power(Expr* base, int64_t exp) {
    if (g_suppress) return NULL;
    Expr* arg = NULL;
    int sin_per = 0, cos_per = 0, kind = 0;
    if (!classify(base, &arg, &sin_per, &cos_per, &kind)) return NULL;
    /* Positive and zero exponents are either already canonical (Sin^k, Tan^k,
     * Sec^k, ...) or handled elsewhere (anything^0 -> 1 by the Power
     * evaluator). Only negative integer exponents need this rewrite. */
    if (exp >= 0) return NULL;
    if (exp < -TRIG_CANON_EXP_BOUND) return NULL;

    int a = sin_per * (int)exp;
    int b = cos_per * (int)exp;

    BasePower out[2];
    size_t emitted = canonical_emit(a, b, arg, kind, out);

    /* For a single trig/hyp ratio base, exactly one of sin_per/cos_per is
     * nonzero (Sin/Cos/Sec/Csc) or sin_per+cos_per == 0 (Tan/Cot), so the
     * accumulated (a, b) always falls in a "single emit" branch of
     * canonical_emit. Defensively bail out otherwise rather than building a
     * Times node here. */
    if (emitted != 1) {
        for (size_t i = 0; i < emitted; i++) {
            expr_free(out[i].base);
            expr_free(out[i].exponent);
        }
        return NULL;
    }

    Expr* nb = out[0].base;
    Expr* ne = out[0].exponent;
    if (ne->type == EXPR_INTEGER && ne->data.integer == 1) {
        expr_free(ne);
        return nb;
    }
    Expr* args[2] = { nb, ne };
    return expr_new_function(expr_new_symbol(SYM_Power), args, 2);
}
