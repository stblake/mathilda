/* diractr.c -- Experimental`FeynCalc`DiracTr: a native Dirac trace kernel.
 *
 * DEFINITION.  The builtin returns exactly -- SameQ, term order included --
 * what this recursive definition (Rolf Mertig's benchmark, Woxi issue #2142)
 * returns when g is a plain symbol:
 *
 *     diracTr[] := 4;
 *     diracTr[a_, b_] := 4 g[a, b];
 *     diracTr[a_, rest___] := Module[{list = {rest}, n},
 *       Sum[(-1)^(m + 1) * Expand[g[a, list[[m]]] * diracTr @@ Delete[list, m]],
 *         {m, 1, Length[list]}]];
 *
 * Unrolled, that is the 4-dimensional trace of n gamma matrices with the
 * metric written g:
 *
 *     Tr[a1 ... an] = 4 Sum_{perfect matchings M} sign(M) Prod_{(i<j) in M} g[ai, aj]
 *
 * with (n-1)!! terms (135135 at n = 14), the first index of each g the one
 * that comes first in the argument list, and sign(M) the product of the
 * (-1)^(m+1) the recursion picks up.  An odd n (and n = 1, an empty Sum) is 0;
 * n = 0 is 4.  g is Global`g.
 *
 * WHY NATIVE.  The definition above spends its time in the evaluator: a
 * Module, a Delete and an Expand per recursion node, and Plus/Times
 * re-canonicalising ever larger sums.  Here the evaluator is never called.
 * The answer is built directly in its canonical form and stamped as a fixed
 * point, so evaluating the returned value is O(1) as well.
 *
 * ALGORITHM.
 *   1. Pair table.  The n(n-1)/2 nodes g[ai, aj] are built once and sorted by
 *      expr_compare; each gets a dense id (equal nodes, from repeated
 *      arguments, share one).  Every term references these nodes by refcount
 *      (expr_copy is a refcount bump), so a term costs one Times node.
 *   2. Enumeration.  A depth-n/2 C recursion over the remaining positions
 *      writes each matching as its n/2 pair ids (sorted, so they are the
 *      term's factor order) and its sign.  Pure integer work, no Expr.
 *   3. Order.
 *        - Every pair node distinct (the common case: distinct arguments):
 *          no two terms share a monomial and no factor repeats.  The canonical
 *          order of such Times terms in a Plus is reverse-lexicographic on the
 *          sorted factor ids -- the LAST factor is the most significant key --
 *          so an LSD radix sort over the n/2 positions orders all terms in
 *          O(T n) with no comparisons.  Within a term, ascending id IS the
 *          Times order (g against g compares argument by argument).
 *        - Repeated pair nodes: equal factors fold into Power[g, k], equal
 *          monomials add their coefficients and zeros drop out, exactly as
 *          Expand would.  Factor and term order then come from expr_compare,
 *          the comparator Times and Plus themselves sort by.
 *   4. Build + stamp.  Plus[Times[c, g.., g..], ...] is assembled bottom-up
 *      and each node is stamped with eval_mark_fixed_point (INERT when the
 *      arguments allow), so the evaluator reuses it without a re-walk.
 *
 * WHEN g IS NOT A PLAIN SYMBOL.  If Global`g has values, rules or attributes,
 * an argument carries UpValues, or Plus/Times/Power have user rules, the
 * definition above would let the evaluator rewrite g[..] and its products.
 * Then the unsorted sum of products is returned inside Expand[...], unstamped,
 * and the evaluator does what the .m would have done.  Correct; not fast.
 */
#include "diractr.h"

#include "context.h"
#include "eval.h"
#include "expr.h"
#include "message.h"
#include "sym_names.h"
#include "symtab.h"
#include "attr.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define DIRACTR_HEAD "Experimental`FeynCalc`DiracTr"

/* ------------------------------------------------------------------------
 * Guards: is the evaluator guaranteed to leave our canonical tree alone?
 * ------------------------------------------------------------------------ */

/* The metric head: a symbol with no values, rules, builtin, or attributes
 * other than Protected.  (NValues and FormatValues change neither the value
 * nor the order of g[a, b], so they do not matter.) */
static bool metric_is_plain(const char* name) {
    const SymbolDef* d = symtab_lookup(name);
    if (!d) return true;
    return !d->own_values && !d->down_values && !d->sub_values && !d->up_values &&
           !d->builtin_func && (d->attributes & ~(uint32_t)ATTR_PROTECTED) == 0;
}

/* Plus, Times and Power: their builtins canonicalise exactly as we do unless
 * someone Unprotected them and added rules. */
static bool arith_head_is_builtin_only(const char* name) {
    const SymbolDef* d = symtab_lookup(name);
    return !d || (!d->own_values && !d->down_values && !d->up_values && !d->sub_values);
}

/* An argument whose head symbol carries UpValues could rewrite g[.., a, ..]. */
static bool arg_has_upvalues(const Expr* a) {
    while (a && a->type == EXPR_FUNCTION) a = a->data.function.head;
    if (!a || a->type != EXPR_SYMBOL) return false;
    const SymbolDef* d = symtab_lookup(a->data.symbol.name);
    return d && d->up_values;
}

/* ------------------------------------------------------------------------
 * Enumeration of perfect matchings
 * ------------------------------------------------------------------------ */

typedef struct {
    size_t   n, h;      /* argument count, factors per term (n / 2) */
    const uint16_t* pid;/* pid[i*n + j], i < j: id of g[ai, aj] */
    uint16_t* ids;      /* T * h factor ids, each row sorted ascending */
    int8_t*  sign;      /* T signs (+1 / -1) */
    size_t   t;         /* terms written so far */
    int*     rem;       /* (h + 1) * n scratch: remaining positions per level */
    uint16_t cur[256];  /* factor ids of the matching being built (h <= 256) */
} Enum;

static void enum_rec(Enum* E, size_t level, size_t nrem, int sign) {
    if (nrem == 0) {
        /* Insertion sort: h is small and the row is usually sorted already
         * (it is whenever the arguments are in canonical order). */
        uint16_t* row = E->ids + E->t * E->h;
        for (size_t k = 0; k < E->h; k++) {
            uint16_t v = E->cur[k];
            size_t j = k;
            while (j > 0 && row[j - 1] > v) { row[j] = row[j - 1]; j--; }
            row[j] = v;
        }
        E->sign[E->t++] = (int8_t)sign;
        return;
    }
    const int* rem = E->rem + level * E->n;
    int* next = E->rem + (level + 1) * E->n;
    const int a = rem[0];
    const uint16_t* pa = E->pid + (size_t)a * E->n;
    /* next = rem minus {rem[0], rem[m]}; built incrementally: going from m to
     * m + 1 only swaps which of rem[m], rem[m+1] is left out. */
    for (size_t k = 1; k + 1 < nrem; k++) next[k - 1] = rem[k + 1];   /* m = 1 */
    for (size_t m = 1; m < nrem; m++) {
        if (m > 1) next[m - 2] = rem[m - 1];   /* rem[m-1] back in, rem[m] out */
        E->cur[level] = pa[rem[m]];
        /* (-1)^(m+1) with m counted from 1 */
        enum_rec(E, level + 1, nrem - 2, (m & 1) ? sign : -sign);
    }
}

/* ------------------------------------------------------------------------
 * Sorting helpers
 * ------------------------------------------------------------------------ */

static const Expr* const* g_sort_nodes;     /* qsort context (single-threaded) */
static int cmp_pair_index(const void* x, const void* y) {
    return expr_compare(g_sort_nodes[*(const size_t*)x], g_sort_nodes[*(const size_t*)y]);
}
static int cmp_expr_ptr(const void* x, const void* y) {
    return expr_compare(*(Expr* const*)x, *(Expr* const*)y);
}
static const uint16_t* g_sort_ids;
static size_t g_sort_h;
static int cmp_row(const void* x, const void* y) {
    const uint16_t* a = g_sort_ids + (size_t)*(const uint32_t*)x * g_sort_h;
    const uint16_t* b = g_sort_ids + (size_t)*(const uint32_t*)y * g_sort_h;
    for (size_t k = 0; k < g_sort_h; k++)
        if (a[k] != b[k]) return a[k] < b[k] ? -1 : 1;
    return 0;
}

/* Stable LSD radix sort of the id rows themselves (with their signs), the
 * last factor the most significant key.  The rows MOVE, rather than an index
 * array being permuted: every pass then streams through memory, where
 * gathering rows through an index missed cache on every element (2M rows x
 * 8 passes at n = 16).  Two factor positions form one digit whenever nid^2
 * buckets stay small, halving the passes.  O(T h). */
static bool radix_sort_rows(uint16_t** ids_p, int8_t** sign_p, size_t T,
                            size_t h, size_t nid) {
    const bool pair_digit = nid <= 256;             /* <= 65536 buckets */
    const size_t nb = pair_digit ? nid * nid : nid;
    uint16_t* ids = *ids_p;
    int8_t* sign = *sign_p;
    uint16_t* ids2 = malloc(T * h * sizeof *ids2);
    int8_t* sign2 = malloc(T * sizeof *sign2);
    size_t* cnt = malloc((nb + 1) * sizeof *cnt);
    if (!ids2 || !sign2 || !cnt) { free(ids2); free(sign2); free(cnt); return false; }
    for (size_t pos = 0; pos < h; ) {
        /* digit = row[pos] (low) and row[pos + 1] (high, more significant) */
        const bool two = pair_digit && pos + 1 < h;
        memset(cnt, 0, (nb + 1) * sizeof *cnt);
        const uint16_t* r = ids + pos;
        for (size_t i = 0; i < T; i++, r += h)
            cnt[(two ? (size_t)r[1] * nid + r[0] : r[0]) + 1]++;
        for (size_t k = 1; k <= nb; k++) cnt[k] += cnt[k - 1];
        const uint16_t* row = ids;
        for (size_t i = 0; i < T; i++, row += h) {
            size_t d = two ? (size_t)row[pos + 1] * nid + row[pos] : row[pos];
            size_t at = cnt[d]++;
            memcpy(ids2 + at * h, row, h * sizeof *row);
            sign2[at] = sign[i];
        }
        uint16_t* ti = ids; ids = ids2; ids2 = ti;
        int8_t* ts = sign; sign = sign2; sign2 = ts;
        pos += two ? 2 : 1;
    }
    free(ids2); free(sign2); free(cnt);
    *ids_p = ids;
    *sign_p = sign;
    return true;
}

/* ------------------------------------------------------------------------
 * The builtin
 * ------------------------------------------------------------------------ */

/* Wrap `terms` (count T, ownership taken; the array itself is freed) as the
 * value: 0, the single term, or Plus[...].  `stamp`: mark the Plus a fixed
 * point (its terms already are). */
static Expr* finish_sum(Expr** terms, size_t T, bool stamp) {
    Expr* out;
    if (T == 0) out = expr_new_integer(0);
    else if (T == 1) out = terms[0];
    else {
        out = expr_new_function(expr_new_symbol(SYM_Plus), terms, T);
        if (stamp) eval_mark_fixed_point(out);
    }
    free(terms);
    return out;
}

Expr* builtin_experimental_diractr(Expr* res) {
    const size_t n = res->data.function.arg_count;
    Expr* const* a = res->data.function.args;

    if (n == 0) return expr_new_integer(4);
    if (n & 1) return expr_new_integer(0);

    const size_t h = n / 2;
    const size_t P = n * (n - 1) / 2;

    /* T = (n-1)!!, guarded: the term index is a uint32_t and a pair id a
     * uint16_t.  Past either, the answer would need hundreds of GB anyway. */
    uint64_t T64 = 1;
    for (uint64_t k = n - 1; k > 1; k -= 2) {
        if (T64 > UINT32_MAX / k) {
            mth_message(DIRACTR_HEAD, "toobig",
                        "The trace of %u gamma matrices has more than 2^32 terms.",
                        (unsigned)n);
            return NULL;
        }
        T64 *= k;
    }
    const size_t T = (size_t)T64;

    /* Metric head.  Global`g is stored under the name the resolver gives. */
    char* gname = context_resolve_name("Global`g");
    if (!gname) return NULL;
    bool plain = metric_is_plain(gname) &&
                 arith_head_is_builtin_only(SYM_Plus) &&
                 arith_head_is_builtin_only(SYM_Times) &&
                 arith_head_is_builtin_only(SYM_Power);
    for (size_t i = 0; plain && i < n; i++)
        if (arg_has_upvalues(a[i])) plain = false;

    /* 1. Pair nodes g[ai, aj], i < j, in position order. */
    Expr** pnode = malloc(P * sizeof *pnode);
    uint16_t* pid = malloc(n * n * sizeof *pid);
    size_t* order = malloc(P * sizeof *order);
    Expr** idnode = malloc(P * sizeof *idnode);      /* id -> node (borrowed) */
    if (!pnode || !pid || !order || !idnode) {
        free(pnode); free(pid); free(order); free(idnode); free(gname);
        mth_message(DIRACTR_HEAD, "nomem", "Out of memory building the pair table.");
        return NULL;
    }
    Expr* gsym = expr_new_symbol(gname);
    free(gname);
    {
        size_t p = 0;
        for (size_t i = 0; i < n; i++)
            for (size_t j = i + 1; j < n; j++, p++) {
                Expr* ga[2] = { expr_copy(a[i]), expr_copy(a[j]) };
                pnode[p] = expr_new_function(expr_copy(gsym), ga, 2);
                if (plain) eval_mark_fixed_point(pnode[p]);
            }
    }
    expr_free(gsym);

    /* Ids.  Plain: by canonical order, equal nodes sharing one id.  Otherwise
     * position order -- nothing is sorted, the evaluator will do it. */
    size_t nid = 0;
    for (size_t p = 0; p < P; p++) order[p] = p;
    if (plain) {
        g_sort_nodes = (const Expr* const*)pnode;
        qsort(order, P, sizeof *order, cmp_pair_index);
    }
    {
        uint16_t* id_of = malloc(P * sizeof *id_of);
        if (!id_of) { nid = 0; goto oom_pairs; }
        for (size_t r = 0; r < P; r++) {
            size_t p = order[r];
            if (r == 0 || !plain || !expr_eq(pnode[order[r - 1]], pnode[p]))
                idnode[nid++] = pnode[p];
            id_of[p] = (uint16_t)(nid - 1);
        }
        size_t p = 0;
        for (size_t i = 0; i < n; i++)
            for (size_t j = i + 1; j < n; j++, p++) pid[i * n + j] = id_of[p];
        free(id_of);
    }
    const bool all_distinct = plain && nid == P;

    /* 2. Enumerate the matchings into flat id rows. */
    Enum E;
    E.n = n; E.h = h; E.pid = pid; E.t = 0;
    E.ids = malloc(T * h * sizeof *E.ids);
    E.sign = malloc(T * sizeof *E.sign);
    E.rem = malloc((h + 1) * n * sizeof *E.rem);
    uint32_t* idx = malloc(T * sizeof *idx);
    if (!E.ids || !E.sign || !E.rem || !idx) {
        free(E.ids); free(E.sign); free(E.rem); free(idx);
        goto oom_pairs;
    }
    for (size_t i = 0; i < n; i++) E.rem[i] = (int)i;
    enum_rec(&E, 0, n, 1);
    free(E.rem);

    Expr* c4  = expr_new_integer(4);
    Expr* cm4 = expr_new_integer(-4);
    Expr* times_sym = expr_new_symbol(SYM_Times);
    Expr* targs[257];
    Expr* out = NULL;

    if (all_distinct) {
        /* 3a. Radix order, 4a. build + stamp. */
        if (!radix_sort_rows(&E.ids, &E.sign, T, h, nid)) goto oom_terms;
        Expr** terms = malloc(T * sizeof *terms);
        if (!terms) goto oom_terms;
        /* Every term is Times[+-4, g..] over the same stamped g nodes.  When
         * those all carry one stamp (one clock, one GROUND/INERT flag), the
         * verdict eval_mark_fixed_point reaches is the same for every term:
         * prove it once, copy it T - 1 times. */
        bool uniform = true;
        for (size_t k = 1; k < nid && uniform; k++)
            uniform = idnode[k]->last_evaluated_at == idnode[0]->last_evaluated_at;
        for (size_t k = 0; k < T; k++) {
            const uint16_t* row = E.ids + k * h;
            targs[0] = expr_copy(E.sign[k] > 0 ? c4 : cm4);
            for (size_t f = 0; f < h; f++) targs[f + 1] = expr_copy(idnode[row[f]]);
            terms[k] = expr_new_function(expr_copy(times_sym), targs, h + 1);
            if (uniform && k > 0) eval_mark_fixed_point_like(terms[k], terms[0]);
            else eval_mark_fixed_point(terms[k]);
        }
        out = finish_sum(terms, T, true);
    } else if (plain) {
        /* 3b. Repeated factors: group equal monomials, add coefficients. */
        for (size_t t = 0; t < T; t++) idx[t] = (uint32_t)t;
        g_sort_ids = E.ids; g_sort_h = h;
        qsort(idx, T, sizeof *idx, cmp_row);
        Expr** terms = malloc(T * sizeof *terms);
        if (!terms) goto oom_terms;
        size_t nt = 0;
        Expr* pow_sym = expr_new_symbol(SYM_Power);
        for (size_t k = 0; k < T; ) {
            const uint16_t* row = E.ids + (size_t)idx[k] * h;
            int64_t coef = 0;
            size_t k2 = k;
            while (k2 < T && cmp_row(&idx[k], &idx[k2]) == 0) coef += 4 * E.sign[idx[k2++]];
            k = k2;
            if (coef == 0) continue;
            size_t nf = 0;
            targs[nf++] = expr_new_integer(coef);
            for (size_t f = 0; f < h; ) {
                size_t f2 = f;
                while (f2 < h && row[f2] == row[f]) f2++;
                Expr* base = expr_copy(idnode[row[f]]);
                if (f2 - f == 1) targs[nf++] = base;
                else {
                    Expr* pa[2] = { base, expr_new_integer((int64_t)(f2 - f)) };
                    Expr* pw = expr_new_function(expr_copy(pow_sym), pa, 2);
                    eval_mark_fixed_point(pw);
                    targs[nf++] = pw;
                }
                f = f2;
            }
            /* Times order is expr_compare order; the number stays first. */
            qsort(targs + 1, nf - 1, sizeof *targs, cmp_expr_ptr);
            terms[nt] = expr_new_function(expr_copy(times_sym), targs, nf);
            eval_mark_fixed_point(terms[nt]);
            nt++;
        }
        expr_free(pow_sym);
        qsort(terms, nt, sizeof *terms, cmp_expr_ptr);
        out = finish_sum(terms, nt, true);
    } else {
        /* 3c. g (or an argument, or arithmetic) carries rules: hand the raw
         * sum to Expand and let the evaluator apply them, as the .m would. */
        Expr** terms = malloc(T * sizeof *terms);
        if (!terms) goto oom_terms;
        for (size_t t = 0; t < T; t++) {
            const uint16_t* row = E.ids + t * h;
            targs[0] = expr_copy(E.sign[t] > 0 ? c4 : cm4);
            for (size_t f = 0; f < h; f++) targs[f + 1] = expr_copy(idnode[row[f]]);
            terms[t] = expr_new_function(expr_copy(times_sym), targs, h + 1);
        }
        Expr* sum = finish_sum(terms, T, false);
        Expr* ea[1] = { sum };
        out = expr_new_function(expr_new_symbol(SYM_Expand), ea, 1);
    }

oom_terms:
    if (!out) mth_message(DIRACTR_HEAD, "nomem", "Out of memory building the trace.");
    expr_free(c4); expr_free(cm4); expr_free(times_sym);
    free(E.ids); free(E.sign); free(idx);
    for (size_t p = 0; p < P; p++) expr_free(pnode[p]);
    free(pnode); free(pid); free(order); free(idnode);
    return out;

oom_pairs:
    mth_message(DIRACTR_HEAD, "nomem", "Out of memory building the trace.");
    for (size_t p = 0; p < P; p++) expr_free(pnode[p]);
    free(pnode); free(pid); free(order); free(idnode);
    return NULL;
}

void experimental_feyncalc_init(void) {
    symtab_add_builtin(DIRACTR_HEAD, builtin_experimental_diractr);
    symtab_get_def(DIRACTR_HEAD)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(DIRACTR_HEAD,
        "Experimental`FeynCalc`DiracTr[a1, ..., an] is the trace of n Dirac gamma "
        "matrices in 4 dimensions with metric g (Global`g): "
        "4 Sum over perfect matchings M of sign(M) Product of g[ai, aj], "
        "(n-1)!! terms, each g's first index the earlier argument. "
        "DiracTr[] is 4; an odd n gives 0. It is exactly (SameQ, term order "
        "included) the value of the recursion diracTr[] := 4; "
        "diracTr[a_, b_] := 4 g[a, b]; diracTr[a_, rest___] := "
        "Sum[(-1)^(m+1) Expand[g[a, {rest}[[m]]] diracTr @@ Delete[{rest}, m]], "
        "{m, Length[{rest}]}]. Repeated arguments fold equal factors into powers "
        "and combine equal terms; indices are not contracted. Native C, evaluator-"
        "free: the result is built in canonical order and stamped as evaluated. "
        "When g, an argument or Plus/Times/Power carries rules, the sum is returned "
        "through Expand so the rules apply.");
}
