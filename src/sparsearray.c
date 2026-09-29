/*
 * sparsearray.c -- the dense form of a SparseArray specification (Normal).
 *
 * See sparsearray.h for the accepted forms. The work is one flat buffer of
 * prod(dims) slots in row-major order, filled rule by rule (first writer wins),
 * then topped up with the default and folded into nested Lists. Explicit
 * positions and Bands are written directly; only pattern rules visit every
 * slot, which is inherent: the dense result has prod(dims) entries anyway.
 */

#include "sparsearray.h"
#include "sym_names.h"
#include "match.h"
#include "arithmetic.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* Largest dense result Normal will build. Beyond this the call is left
 * unevaluated rather than attempting a multi-gigabyte allocation. */
#define SA_MAX_ELEMENTS ((int64_t)1 << 27)
#define SA_MAX_RANK 32

static bool head_named(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

bool is_sparse_array(const Expr* e) {
    return head_named(e, SYM_SparseArray);
}

static bool is_rule(const Expr* e) {
    return (head_named(e, SYM_Rule) || head_named(e, SYM_RuleDelayed))
        && e->data.function.arg_count == 2;
}

/* A List whose every element is a machine Integer. */
static bool is_int_list(const Expr* e) {
    if (!head_named(e, SYM_List)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (e->data.function.args[i]->type != EXPR_INTEGER) return false;
    return true;
}

/* A position: an Integer (rank 1) or a non-empty List of Integers. */
static bool is_position(const Expr* e) {
    return e->type == EXPR_INTEGER || (is_int_list(e) && e->data.function.arg_count > 0);
}

/* A List of positions all of one rank: the LHS of {p1, p2} -> {v1, v2}. */
static bool is_position_list(const Expr* e) {
    if (!head_named(e, SYM_List) || e->data.function.arg_count == 0) return false;
    size_t rank = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        const Expr* p = e->data.function.args[i];
        if (!is_int_list(p) || p->data.function.arg_count == 0) return false;
        if (i == 0) rank = p->data.function.arg_count;
        else if (p->data.function.arg_count != rank) return false;
    }
    return true;
}

static size_t position_rank(const Expr* p) {
    return p->type == EXPR_INTEGER ? 1 : p->data.function.arg_count;
}

static int64_t position_coord(const Expr* p, size_t k) {
    return p->type == EXPR_INTEGER ? p->data.integer : p->data.function.args[k]->data.integer;
}

/* ---- the dense buffer ---------------------------------------------------- */

typedef struct {
    size_t   rank;
    int64_t  dims[SA_MAX_RANK];
    int64_t  total;
    Expr**   slots;        /* total entries, NULL = not yet written */
} Dense;

/* Row-major offset of the 1-based coordinates `c`, resolving negative indices
 * from the end; -1 when out of range. */
static int64_t dense_offset(const Dense* d, const int64_t* c) {
    int64_t off = 0;
    for (size_t k = 0; k < d->rank; k++) {
        int64_t i = c[k];
        if (i < 0) i = d->dims[k] + 1 + i;
        if (i < 1 || i > d->dims[k]) return -1;
        off = off * d->dims[k] + (i - 1);
    }
    return off;
}

/* Write `v` (taking ownership) at `off` unless an earlier rule got there. */
static void dense_put(Dense* d, int64_t off, Expr* v) {
    if (off < 0 || d->slots[off]) { expr_free(v); return; }
    d->slots[off] = v;
}

static bool dense_put_position(Dense* d, const Expr* pos, const Expr* val) {
    if (position_rank(pos) != d->rank) return false;
    int64_t c[SA_MAX_RANK];
    for (size_t k = 0; k < d->rank; k++) c[k] = position_coord(pos, k);
    int64_t off = dense_offset(d, c);
    if (off < 0) return false;             /* position outside the array */
    dense_put(d, off, expr_copy((Expr*)val));
    return true;
}

/* Band[start], Band[start, end], Band[start, end, step] -> v. `v` a scalar
 * fills the whole band; a List supplies successive entries and ends it early. */
static bool dense_put_band(Dense* d, const Expr* band, const Expr* val) {
    size_t na = band->data.function.arg_count;
    if (na < 1 || na > 3) return false;
    const Expr* st = band->data.function.args[0];
    if (!is_int_list(st) || st->data.function.arg_count != d->rank) return false;
    int64_t cur[SA_MAX_RANK], end[SA_MAX_RANK], step[SA_MAX_RANK];
    for (size_t k = 0; k < d->rank; k++) {
        cur[k]  = st->data.function.args[k]->data.integer;
        if (cur[k] < 0) cur[k] = d->dims[k] + 1 + cur[k];
        end[k]  = d->dims[k];
        step[k] = 1;
    }
    if (na >= 2) {
        const Expr* en = band->data.function.args[1];
        if (!is_int_list(en) || en->data.function.arg_count != d->rank) return false;
        for (size_t k = 0; k < d->rank; k++) {
            end[k] = en->data.function.args[k]->data.integer;
            if (end[k] < 0) end[k] = d->dims[k] + 1 + end[k];
        }
    }
    if (na == 3) {
        const Expr* sp = band->data.function.args[2];
        if (!is_int_list(sp) || sp->data.function.arg_count != d->rank) return false;
        for (size_t k = 0; k < d->rank; k++) {
            step[k] = sp->data.function.args[k]->data.integer;
            if (step[k] <= 0) return false;
        }
    }
    bool list_val = head_named(val, SYM_List);
    size_t n = list_val ? val->data.function.arg_count : SIZE_MAX;
    for (size_t j = 0; j < n; j++) {
        for (size_t k = 0; k < d->rank; k++)
            if (cur[k] < 1 || cur[k] > end[k] || cur[k] > d->dims[k]) return true;
        const Expr* v = list_val ? val->data.function.args[j] : val;
        dense_put(d, dense_offset(d, cur), expr_copy((Expr*)v));
        for (size_t k = 0; k < d->rank; k++) cur[k] += step[k];
    }
    return true;
}

/* A pattern rule: try it at every still-empty slot, matching the position
 * List {i, j, ...} (and, for a vector, the bare index i as well). */
static void dense_put_pattern(Dense* d, const Expr* lhs, const Expr* rhs) {
    int64_t c[SA_MAX_RANK];
    for (size_t k = 0; k < d->rank; k++) c[k] = 1;
    for (int64_t off = 0; off < d->total; off++) {
        if (!d->slots[off]) {
            Expr** items = malloc(sizeof(Expr*) * d->rank);
            for (size_t k = 0; k < d->rank; k++) items[k] = expr_new_integer(c[k]);
            Expr* pos = expr_new_function(expr_new_symbol(SYM_List), items, d->rank);
            free(items);
            MatchEnv* env = env_new();
            bool ok = match(pos, (Expr*)lhs, env);
            if (!ok && d->rank == 1) {
                env_free(env);
                env = env_new();
                Expr* bare = expr_new_integer(c[0]);
                ok = match(bare, (Expr*)lhs, env);
                expr_free(bare);
            }
            if (ok) d->slots[off] = replace_bindings((Expr*)rhs, env);
            env_free(env);
            expr_free(pos);
        }
        for (size_t k = d->rank; k-- > 0; ) {          /* odometer advance */
            if (++c[k] <= d->dims[k]) break;
            c[k] = 1;
        }
    }
}

/* Fold slots [off, off + prod(dims[level..])) into nested Lists, moving the
 * slot Exprs into the result. */
static Expr* dense_fold(Dense* d, size_t level, int64_t* off) {
    int64_t n = d->dims[level];
    Expr** items = malloc(sizeof(Expr*) * (size_t)(n ? n : 1));
    for (int64_t i = 0; i < n; i++) {
        if (level + 1 == d->rank) {
            items[i] = d->slots[*off];
            d->slots[*off] = NULL;
            (*off)++;
        } else {
            items[i] = dense_fold(d, level + 1, off);
        }
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);
    return out;
}

static void dense_free(Dense* d) {
    if (!d->slots) return;
    for (int64_t i = 0; i < d->total; i++) if (d->slots[i]) expr_free(d->slots[i]);
    free(d->slots);
    d->slots = NULL;
}

static bool dense_alloc(Dense* d) {
    d->total = 1;
    for (size_t k = 0; k < d->rank; k++) {
        if (d->dims[k] < 0) return false;
        if (d->dims[k] > 0 && d->total > SA_MAX_ELEMENTS / d->dims[k]) return false;
        d->total *= d->dims[k];
    }
    d->slots = calloc((size_t)(d->total ? d->total : 1), sizeof(Expr*));
    return d->slots != NULL;
}

/* Fill the unwritten slots with `def`, fold, and release the buffer. */
static Expr* dense_finish(Dense* d, const Expr* def) {
    for (int64_t i = 0; i < d->total; i++)
        if (!d->slots[i]) d->slots[i] = def ? expr_copy((Expr*)def) : expr_new_integer(0);
    int64_t off = 0;
    Expr* out = d->rank ? dense_fold(d, 0, &off) : NULL;
    dense_free(d);
    return out;
}

/* Parse a dims spec (Integer n or a List of non-negative Integers). */
static bool parse_dims(const Expr* spec, Dense* d) {
    if (spec->type == EXPR_INTEGER) {
        d->rank = 1;
        d->dims[0] = spec->data.integer;
        return d->dims[0] >= 0;
    }
    if (!is_int_list(spec) || spec->data.function.arg_count == 0
        || spec->data.function.arg_count > SA_MAX_RANK) return false;
    d->rank = spec->data.function.arg_count;
    for (size_t k = 0; k < d->rank; k++) {
        d->dims[k] = spec->data.function.args[k]->data.integer;
        if (d->dims[k] < 0) return false;
    }
    return true;
}

/* ---- the rule forms ------------------------------------------------------ */

/* Collect the rules of a specification into a flat array (borrowed). A single
 * rule, a List of rules, or the empty List. */
static bool collect_rules(const Expr* spec, const Expr*** rules, size_t* n) {
    if (is_rule(spec)) {
        *rules = malloc(sizeof(Expr*));
        (*rules)[0] = spec;
        *n = 1;
        return true;
    }
    if (!head_named(spec, SYM_List)) return false;
    size_t m = spec->data.function.arg_count;
    for (size_t i = 0; i < m; i++) if (!is_rule(spec->data.function.args[i])) return false;
    *rules = malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) (*rules)[i] = spec->data.function.args[i];
    *n = m;
    return true;
}

/* Infer dims from explicit positions (the per-coordinate maximum). Fails on a
 * pattern or Band rule, which do not determine a size. */
static bool infer_dims(const Expr** rules, size_t n, Dense* d) {
    d->rank = 0;
    for (size_t r = 0; r < n; r++) {
        const Expr* lhs = rules[r]->data.function.args[0];
        const Expr* rhs = rules[r]->data.function.args[1];
        size_t np = 1;
        bool many = is_position_list(lhs) && head_named(rhs, SYM_List);
        if (many) np = lhs->data.function.arg_count;
        else if (!is_position(lhs)) return false;
        for (size_t j = 0; j < np; j++) {
            const Expr* p = many ? lhs->data.function.args[j] : lhs;
            size_t rk = position_rank(p);
            if (rk > SA_MAX_RANK) return false;
            if (d->rank == 0) {
                d->rank = rk;
                for (size_t k = 0; k < rk; k++) d->dims[k] = 0;
            } else if (rk != d->rank) {
                return false;
            }
            for (size_t k = 0; k < rk; k++) {
                int64_t c = position_coord(p, k);
                if (c < 1) return false;            /* negative index needs dims */
                if (c > d->dims[k]) d->dims[k] = c;
            }
        }
    }
    return d->rank > 0;
}

static bool apply_rules(Dense* d, const Expr** rules, size_t n) {
    for (size_t r = 0; r < n; r++) {
        const Expr* lhs = rules[r]->data.function.args[0];
        const Expr* rhs = rules[r]->data.function.args[1];
        if (head_named(lhs, SYM_Band)) {
            if (!dense_put_band(d, lhs, rhs)) return false;
        } else if (is_position_list(lhs) && head_named(rhs, SYM_List)
                   && position_rank(lhs->data.function.args[0]) == d->rank) {
            if (lhs->data.function.arg_count != rhs->data.function.arg_count) return false;
            for (size_t j = 0; j < lhs->data.function.arg_count; j++)
                if (!dense_put_position(d, lhs->data.function.args[j], rhs->data.function.args[j]))
                    return false;
        } else if (is_position(lhs)) {
            if (!dense_put_position(d, lhs, rhs)) return false;
        } else {
            dense_put_pattern(d, lhs, rhs);
        }
    }
    return true;
}

/* ---- the dense-list form ------------------------------------------------- */

/* Dimensions of a rectangular nested List (full depth). */
static void list_dims(const Expr* e, Dense* d) {
    d->rank = 0;
    while (head_named(e, SYM_List) && d->rank < SA_MAX_RANK) {
        size_t n = e->data.function.arg_count;
        d->dims[d->rank++] = (int64_t)n;
        if (n == 0) break;
        const Expr* first = e->data.function.args[0];
        /* every sibling must agree on being a List of the same length */
        for (size_t i = 0; i < n; i++) {
            const Expr* s = e->data.function.args[i];
            if (head_named(first, SYM_List) != head_named(s, SYM_List)
                || (head_named(s, SYM_List)
                    && s->data.function.arg_count != first->data.function.arg_count))
                return;
        }
        e = first;
    }
}

/* Copy the entries of a dense List into `d` (dims possibly larger or smaller:
 * entries outside d are dropped, missing ones take the default). */
static void put_dense_list(Dense* d, const Expr* e, size_t level, int64_t* c) {
    if (level == d->rank) {
        dense_put(d, dense_offset(d, c), expr_copy((Expr*)e));
        return;
    }
    if (!head_named(e, SYM_List)) return;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        c[level] = (int64_t)i + 1;
        put_dense_list(d, e->data.function.args[i], level + 1, c);
    }
}

/* ---- the internal form --------------------------------------------------- */

/* SparseArray[Automatic, dims, default, {1, {rowptr, colidx}, vals}]. */
static Expr* from_internal(const Expr* sa) {
    const Expr* const* A = (const Expr* const*)sa->data.function.args;
    Dense d; memset(&d, 0, sizeof d);
    if (!parse_dims(A[1], &d)) return NULL;
    const Expr* data = A[3];
    if (!head_named(data, SYM_List) || data->data.function.arg_count != 3) return NULL;
    const Expr* idx  = data->data.function.args[1];
    const Expr* vals = data->data.function.args[2];
    if (!head_named(idx, SYM_List) || idx->data.function.arg_count != 2
        || !head_named(vals, SYM_List)) return NULL;
    const Expr* rowptr = idx->data.function.args[0];
    const Expr* colidx = idx->data.function.args[1];
    if (!is_int_list(rowptr) || !head_named(colidx, SYM_List)) return NULL;
    size_t nnz = vals->data.function.arg_count;
    if (colidx->data.function.arg_count != nnz) return NULL;
    size_t nrow = rowptr->data.function.arg_count;
    /* rank 1 has a single pseudo-row; otherwise one row per first-dim index */
    size_t want_rows = d.rank == 1 ? 2 : (size_t)d.dims[0] + 1;
    size_t tail = d.rank == 1 ? 1 : d.rank - 1;
    if (nrow != want_rows) return NULL;
    if (!dense_alloc(&d)) return NULL;
    int64_t c[SA_MAX_RANK];
    for (size_t row = 0; row + 1 < nrow; row++) {
        int64_t lo = rowptr->data.function.args[row]->data.integer;
        int64_t hi = rowptr->data.function.args[row + 1]->data.integer;
        if (lo < 0 || hi < lo || hi > (int64_t)nnz) { dense_free(&d); return NULL; }
        for (int64_t k = lo; k < hi; k++) {
            const Expr* ci = colidx->data.function.args[k];
            if (!is_int_list(ci) || ci->data.function.arg_count != tail) {
                dense_free(&d); return NULL;
            }
            size_t j = 0;
            if (d.rank > 1) c[j++] = (int64_t)row + 1;
            for (size_t t = 0; t < tail; t++) c[j++] = ci->data.function.args[t]->data.integer;
            int64_t off = dense_offset(&d, c);
            if (off < 0) { dense_free(&d); return NULL; }
            dense_put(&d, off, expr_copy(vals->data.function.args[k]));
        }
    }
    return dense_finish(&d, A[2]);
}

/* ---- entry point --------------------------------------------------------- */

Expr* sparse_array_to_dense(const Expr* sa) {
    if (!is_sparse_array(sa)) return NULL;
    size_t na = sa->data.function.arg_count;
    if (na < 1 || na > 4) return NULL;
    const Expr* spec = sa->data.function.args[0];

    if (spec->type == EXPR_SYMBOL && spec->data.symbol.name == SYM_Automatic) {
        return na == 4 ? from_internal(sa) : NULL;
    }
    if (na == 4) return NULL;

    const Expr* dims_e = na >= 2 ? sa->data.function.args[1] : NULL;
    const Expr* def    = na == 3 ? sa->data.function.args[2] : NULL;
    Dense d; memset(&d, 0, sizeof d);

    const Expr** rules = NULL;
    size_t nrules = 0;
    bool rule_form = collect_rules(spec, &rules, &nrules);
    /* {} is both an empty rule list and an empty dense list: with dims it is
     * an all-default array either way. */

    if (rule_form) {
        bool ok = dims_e ? parse_dims(dims_e, &d) : infer_dims(rules, nrules, &d);
        if (ok) ok = dense_alloc(&d);
        if (ok) ok = apply_rules(&d, rules, nrules);
        free(rules);
        if (!ok) { dense_free(&d); return NULL; }
        return dense_finish(&d, def);
    }

    if (!head_named(spec, SYM_List)) return NULL;
    /* A dense array. Without dims it is its own Normal. */
    if (!dims_e) return expr_copy((Expr*)spec);
    Dense src; memset(&src, 0, sizeof src);
    list_dims(spec, &src);
    if (!parse_dims(dims_e, &d) || d.rank > src.rank) return NULL;
    if (!dense_alloc(&d)) return NULL;
    int64_t c[SA_MAX_RANK];
    put_dense_list(&d, spec, 0, c);
    return dense_finish(&d, def);
}
