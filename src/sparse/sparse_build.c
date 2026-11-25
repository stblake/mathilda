/*
 * sparse_build.c -- the SparseArray constructor.
 *
 *   SparseArray[{pos1 -> v1, ...}]                    dims inferred from positions
 *   SparseArray[rules, dims]                          dims: n, {d1, ...} or Automatic
 *   SparseArray[rules, dims, val]                     background val
 *   SparseArray[{pos1, pos2, ...} -> {v1, v2, ...}]   (or -> v for all of them)
 *   SparseArray[patt -> v | patt :> v, dims]          every matching position
 *   SparseArray[Band[start(, end(, step))] -> v, ...] diagonals and subblocks
 *   SparseArray[list(, dims(, val))]                  a dense (packed) array
 *   SparseArray[s(, dims(, val))]                     minimise / resize / rebackground
 *   SparseArray[Automatic, dims, val, {1, {rowptr, colidx}, vals}]   canonical form
 *
 * Every form funnels into the COO builder (sparse_core.c), which sorts into
 * canonical order with first-rule-wins and drops background-valued entries.
 * A malformed specification emits the SparseArray:: message and leaves the call
 * unevaluated. The large-array paths -- a packed {positions} -> {values} rule
 * and a packed dense array -- read the machine buffers directly and never box
 * an element.
 */

#include "sparse.h"
#include "sym_names.h"
#include "eval.h"
#include "match.h"
#include "message.h"
#include "ndarray.h"
#include "print.h"
#include "pack.h"
#include "linkedlist.h"   /* dense rows may be LinkedLists */
#include <stdlib.h>
#include <string.h>

/* ======================================================================== */
/* Helpers                                                                  */
/* ======================================================================== */

static bool head_is_sym(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

static bool is_sym(const Expr* e, const char* sym) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == sym;
}

static bool is_rule(const Expr* e) {
    return (head_is_sym(e, SYM_Rule) || head_is_sym(e, SYM_RuleDelayed)) &&
           e->data.function.arg_count == 2;
}

/* A plain List, a packed list, or a visible NDArray: something with list
 * elements. */
static bool lv_is_list(const Expr* e) {
    return head_is_sym(e, SYM_List) || (e && e->type == EXPR_NDARRAY);
}

/* A value no SparseArray may hold as an element: any array (a List, a buffer,
 * or another SparseArray). */
static bool is_array_value(const Expr* e) {
    return lv_is_list(e) || (e && e->type == EXPR_SPARSE);
}

static int64_t lv_len(const Expr* e) {
    if (e->type == EXPR_NDARRAY) return e->data.ndarray.dims[0];
    return (int64_t)e->data.function.arg_count;
}

/* Element i of a list (fresh reference). For a buffer of rank >= 2 it is the
 * sub-array, as Part gives it. */
static Expr* lv_item(const Expr* e, int64_t i) {
    if (e->type != EXPR_NDARRAY) return expr_copy(e->data.function.args[i]);
    if (e->data.ndarray.rank == 1) return ndarray_element_to_expr(e, (size_t)i);
    Expr* ix = expr_new_integer(i + 1);
    bool degrade = false;
    Expr* r = ndarray_part(e, &ix, 1, &degrade);
    expr_free(ix);
    return r;
}

/* Is `e` an Integer (machine), reading through a buffer element if needed? */
static bool lv_int_at(const Expr* e, int64_t i, int64_t* out) {
    if (e->type == EXPR_NDARRAY) {
        if (e->data.ndarray.rank != 1 || e->data.ndarray.dtype != NDT_INT64) return false;
        *out = ndt_get_i(e->data.ndarray.data, (size_t)i, NDT_INT64);
        return true;
    }
    const Expr* x = e->data.function.args[i];
    if (x->type != EXPR_INTEGER) return false;
    *out = x->data.integer;
    return true;
}

/* A position: an Integer, or a non-empty List of Integers. Reads it into `c`
 * (1-based, possibly negative) and its length into *rank. */
static bool read_position(const Expr* e, int64_t* c, int* rank) {
    if (e->type == EXPR_INTEGER) { c[0] = e->data.integer; *rank = 1; return true; }
    if (!lv_is_list(e)) return false;
    if (e->type == EXPR_NDARRAY && e->data.ndarray.rank != 1) return false;
    int64_t n = lv_len(e);
    if (n < 1 || n > SP_MAX_RANK) return false;
    for (int64_t k = 0; k < n; k++) if (!lv_int_at(e, k, &c[k])) return false;
    *rank = (int)n;
    return true;
}

/* Does `e` contain a pattern object anywhere? */
static bool has_pattern(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* n = h->data.symbol.name;
        if (n == SYM_Blank || n == SYM_BlankSequence || n == SYM_BlankNullSequence ||
            n == SYM_Pattern || n == SYM_Alternatives || n == SYM_Condition ||
            n == SYM_PatternTest || n == SYM_Optional || n == SYM_Except ||
            n == SYM_Repeated || n == SYM_RepeatedNull)
            return true;
    }
    if (has_pattern(h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_pattern(e->data.function.args[i])) return true;
    return false;
}

/* `_` or `x_`. */
static bool is_plain_blank(const Expr* e, const char** name) {
    if (head_is_sym(e, SYM_Blank) && e->data.function.arg_count == 0) { *name = NULL; return true; }
    if (head_is_sym(e, SYM_Pattern) && e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_SYMBOL &&
        head_is_sym(e->data.function.args[1], SYM_Blank) &&
        e->data.function.args[1]->data.function.arg_count == 0) {
        *name = e->data.function.args[0]->data.symbol.name;
        return true;
    }
    return false;
}

/* A pattern that matches every position of a rank-`rank` array: `_`, `x_`, or
 * a List of `rank` plain blanks with distinct names. */
static bool matches_everything(const Expr* p, int rank) {
    const char* nm;
    if (is_plain_blank(p, &nm)) return true;
    if (!head_is_sym(p, SYM_List) || (int)p->data.function.arg_count != rank) return false;
    const char* names[SP_MAX_RANK];
    for (int k = 0; k < rank; k++) {
        if (!is_plain_blank(p->data.function.args[k], &names[k])) return false;
        for (int j = 0; j < k; j++) if (names[k] && names[j] == names[k]) return false;
    }
    return true;
}

/* {i_, i_, ..., i_}: the main diagonal. */
static bool is_diagonal_pattern(const Expr* p, int rank) {
    if (rank < 2 || !head_is_sym(p, SYM_List) || (int)p->data.function.arg_count != rank)
        return false;
    const char* first = NULL;
    for (int k = 0; k < rank; k++) {
        const char* nm;
        if (!is_plain_blank(p->data.function.args[k], &nm) || !nm) return false;
        if (k == 0) first = nm;
        else if (nm != first) return false;
    }
    return true;
}

static char* show(const Expr* e) { return expr_to_string_outputform((Expr*)e); }

/* Emit SparseArray::tag with up to three expression arguments (NULL to omit). */
static void sa_msg(const char* tag, const char* fmt, const Expr* a, const Expr* b, const Expr* c) {
    char* sa = a ? show(a) : NULL;
    char* sb = b ? show(b) : NULL;
    char* sc = c ? show(c) : NULL;
    mth_message("SparseArray", tag, fmt, sa ? sa : "", sb ? sb : "", sc ? sc : "");
    free(sa); free(sb); free(sc);
}

static Expr* dims_expr(int rank, const int64_t* dims) {
    Expr** it = malloc(sizeof(Expr*) * (size_t)(rank ? rank : 1));
    for (int k = 0; k < rank; k++) it[k] = expr_new_integer(dims[k]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), it, (size_t)rank);
    free(it);
    return out;
}

/* The empty array of the given dims (a zero among them) as a nested List:
 * {} for {0} or {0, 3}, {{}, {}, {}} for {3, 0}. */
static Expr* empty_nested(int rank, const int64_t* dims) {
    int64_t n = rank > 1 && dims[0] > 0 ? dims[0] : 0;
    Expr** it = malloc(sizeof(Expr*) * (size_t)(n ? n : 1));
    for (int64_t k = 0; k < n; k++) it[k] = empty_nested(rank - 1, dims + 1);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), it, (size_t)n);
    free(it);
    return out;
}

/* ---- dims specification -------------------------------------------------- */

typedef struct {
    bool    given;          /* a dims argument other than Automatic */
    bool    scalar;         /* dims given as one Integer n (every dimension n) */
    int     rank;           /* 0 = not yet known */
    int64_t dims[SP_MAX_RANK];
} DimSpec;

/* Parse the dims argument; false (after SparseArray::adims) when malformed. */
static bool parse_dims(const Expr* d, DimSpec* ds) {
    memset(ds, 0, sizeof *ds);
    if (!d || is_sym(d, SYM_Automatic)) return true;
    ds->given = true;
    if (d->type == EXPR_INTEGER && d->data.integer >= 0) {
        ds->scalar = true;
        ds->dims[0] = d->data.integer;
        return true;
    }
    if (lv_is_list(d) && (d->type != EXPR_NDARRAY || d->data.ndarray.rank == 1)) {
        int64_t n = lv_len(d);
        bool ok = n >= 1 && n <= SP_MAX_RANK;
        for (int64_t k = 0; ok && k < n; k++)
            ok = lv_int_at(d, k, &ds->dims[k]) && ds->dims[k] >= 0;
        if (ok) { ds->rank = (int)n; return true; }
    }
    sa_msg("adims", "Array dimension specification %s should be Automatic, a non-negative "
           "machine integer, or a list of non-negative machine integers.", d, NULL, NULL);
    return false;
}

/* Fix the rank of a scalar dims spec once it is known. */
static void dims_set_rank(DimSpec* ds, int rank) {
    if (ds->rank) return;
    ds->rank = rank;
    if (ds->scalar) for (int k = 1; k < rank; k++) ds->dims[k] = ds->dims[0];
}

/* ======================================================================== */
/* Dense arrays                                                             */
/* ======================================================================== */

/* A dense row: a List, a buffer, or a LinkedList -- whose elements are read
 * from its slots (dl_at: borrowed for a boxed list, a fresh box of a typed
 * slot, *own set, otherwise). The rule paths below keep lv_is_list. */
static bool dl_is_row(const Expr* e) { return lv_is_list(e) || (e && e->type == EXPR_LINKEDLIST); }
static int64_t dl_len(const Expr* e) {
    return e->type == EXPR_LINKEDLIST ? ll_length(e->data.llist) : lv_len(e);
}
static Expr* dl_at(const Expr* e, int64_t i, bool* own) {
    if (e->type == EXPR_LINKEDLIST) return ll_elem_view(e->data.llist, i, own);
    *own = false;
    return e->data.function.args[i];
}

/* Full-depth dimensions of a rectangular nested list; false when ragged, or
 * when a leaf is itself a List at a shallower level. Packed buffers answer
 * from their shape. */
static bool dense_dims(const Expr* e, int* rank, int64_t* dims) {
    if (e->type == EXPR_NDARRAY) {
        *rank = e->data.ndarray.rank;
        memcpy(dims, e->data.ndarray.dims, sizeof(int64_t) * (size_t)*rank);
        return true;
    }
    /* Walk the first branch for the candidate shape, then verify. */
    int r = 0;
    const Expr* x = e;
    Expr* held = NULL;     /* a boxed typed slot: a leaf, so the walk ends there */
    bool ok = true;
    while (dl_is_row(x)) {
        if (x->type == EXPR_NDARRAY) {
            for (int k = 0; k < x->data.ndarray.rank && r < SP_MAX_RANK; k++)
                dims[r++] = x->data.ndarray.dims[k];
            break;
        }
        if (r >= SP_MAX_RANK) { ok = false; break; }
        int64_t n = dl_len(x);
        dims[r++] = n;
        if (n == 0) break;
        bool own;
        Expr* nx = dl_at(x, 0, &own);
        if (own) { if (held) expr_free(held); held = nx; }
        x = nx;
    }
    if (held) expr_free(held);
    *rank = r;
    return ok;
}

/* Does every element of `e` at level `level` agree with the shape? */
static bool dense_check(const Expr* e, int level, int rank, const int64_t* dims) {
    if (level == rank) return !is_array_value(e) && e->type != EXPR_LINKEDLIST;
    if (!dl_is_row(e) || dl_len(e) != dims[level]) return false;
    if (e->type == EXPR_NDARRAY) {
        if (e->data.ndarray.rank != rank - level) return false;
        for (int k = 0; k < e->data.ndarray.rank; k++)
            if (e->data.ndarray.dims[k] != dims[level + k]) return false;
        return true;
    }
    int64_t n = dl_len(e);
    for (int64_t i = 0; i < n; i++) {
        bool own;
        Expr* x = dl_at(e, i, &own);
        bool ok = dense_check(x, level + 1, rank, dims);
        ll_elem_done(x, own);
        if (!ok) return false;
    }
    return true;
}

/* Are all leaves machine Reals / machine Complexes (the "inexact array" whose
 * zero is 0.)? Are all leaves machine numbers (Integer or Real)? */
typedef struct { bool all_inexact; bool all_machine_real_or_int; bool any; } LeafScan;

static void scan_leaves(const Expr* e, int level, int rank, LeafScan* s) {
    if (e->type == EXPR_NDARRAY) {
        NDType dt = e->data.ndarray.dtype;
        s->any = s->any || ndarray_size(e) > 0;
        if (dt == NDT_INT64 || dt == NDT_BOOL) s->all_inexact = false;
        if (ndt_is_complex(dt) || dt == NDT_BOOL) s->all_machine_real_or_int = false;
        return;
    }
    if (level == rank) {
        s->any = true;
        bool real = e->type == EXPR_REAL;
        bool cplx = head_is_sym(e, SYM_Complex) && e->data.function.arg_count == 2 &&
                    e->data.function.args[0]->type == EXPR_REAL &&
                    e->data.function.args[1]->type == EXPR_REAL;
        if (!real && !cplx) s->all_inexact = false;
        if (!real && e->type != EXPR_INTEGER) s->all_machine_real_or_int = false;
        return;
    }
    int64_t n = dl_len(e);
    for (int64_t i = 0; i < n; i++) {
        bool own;
        Expr* x = dl_at(e, i, &own);
        scan_leaves(x, level + 1, rank, s);
        ll_elem_done(x, own);
    }
}

/* Push the entries of a dense list (with dims `src` cropped to `dims`) in
 * row-major order. `to_real` converts Integer leaves to Reals. */
static bool push_dense(SpCoo* c, const Expr* e, int level, int64_t* coord,
                       const Expr* bg, bool to_real) {
    int rank = c->rank;
    if (e->type == EXPR_NDARRAY && level < rank) {
        /* A (sub-)buffer: walk it flat. */
        int r = e->data.ndarray.rank;
        NDType dt = e->data.ndarray.dtype;
        int64_t total = (int64_t)ndarray_size(e);
        int64_t sub[SP_MAX_RANK];
        memset(sub, 0, sizeof sub);
        double bre = 0;
        bool bg_real = bg->type == EXPR_REAL, bg_int = bg->type == EXPR_INTEGER;
        if (bg_real) bre = bg->data.real;
        if (bg_int) bre = (double)bg->data.integer;
        for (int64_t i = 0; i < total; i++) {
            /* coordinates of flat element i */
            bool inside = true;
            for (int k = 0; k < r; k++) {
                coord[level + k] = sub[k];
                if (sub[k] >= c->dims[level + k]) inside = false;
            }
            if (inside) {
                bool ok = true;
                if (dt == NDT_FLOAT64 || dt == NDT_FLOAT32) {
                    double re, im;
                    ndt_get(e->data.ndarray.data, (size_t)i, dt, &re, &im);
                    if (!((bg_real || bg_int) && re == bre))
                        ok = c->kind == SPV_F64 ? sp_coo_push_f64(c, coord, re)
                                                : sp_coo_push_expr(c, coord, expr_new_real(re));
                } else if (dt == NDT_INT64) {
                    int64_t v = ndt_get_i(e->data.ndarray.data, (size_t)i, dt);
                    bool is_bg = bg_int ? v == bg->data.integer : (bg_real && (double)v == bre);
                    if (!is_bg) {
                        if (c->kind == SPV_INT64) ok = sp_coo_push_i64(c, coord, v);
                        else if (c->kind == SPV_F64) ok = sp_coo_push_f64(c, coord, (double)v);
                        else ok = sp_coo_push_expr(c, coord,
                                                   to_real ? expr_new_real((double)v) : expr_new_integer(v));
                    }
                } else {
                    Expr* v = ndarray_element_to_expr(e, (size_t)i);
                    if (sp_value_is_background(v, bg)) expr_free(v);
                    else ok = sp_coo_push_expr(c, coord, v);
                }
                if (!ok) return false;
            }
            for (int k = r - 1; k >= 0; k--) {       /* odometer */
                if (++sub[k] < e->data.ndarray.dims[k]) break;
                sub[k] = 0;
            }
        }
        return true;
    }
    if (level == rank) {
        Expr* v = (to_real && e->type == EXPR_INTEGER) ? expr_new_real((double)e->data.integer)
                                                     : expr_copy((Expr*)e);
        if (sp_value_is_background(v, bg)) { expr_free(v); return true; }
        return sp_coo_push_expr(c, coord, v);
    }
    int64_t n = dl_len(e);
    if (n > c->dims[level]) n = c->dims[level];
    for (int64_t i = 0; i < n; i++) {
        coord[level] = i;
        bool own;
        Expr* x = dl_at(e, i, &own);
        bool ok = push_dense(c, x, level + 1, coord, bg, to_real);
        ll_elem_done(x, own);
        if (!ok) return false;
    }
    return true;
}

/* A packed Real or Integer array kept at its own dims: two passes over the
 * flat buffer -- count the non-background elements, then emit their composite
 * keys straight from the flat index (row = i / span, inner key = i mod span),
 * already in canonical order. NULL when the shape or dtype is not this case
 * (the caller then takes the general path). `bg` is borrowed. */
static Expr* from_packed_dense(const Expr* nd, const Expr* bg, bool to_real) {
    NDType dt = nd->data.ndarray.dtype;
    if (dt != NDT_FLOAT64 && dt != NDT_INT64) return NULL;
    int rank = nd->data.ndarray.rank;
    const int64_t* dims = nd->data.ndarray.dims;
    int kbits;
    if (rank > SP_MAX_RANK || !sp_coo_composite(rank, dims, &kbits)) return NULL;
    int64_t total = sp_total(rank, dims);
    if (total <= 0) return NULL;
    int64_t span = rank >= 2 ? total / dims[0] : total;
    bool bg_int = bg->type == EXPR_INTEGER, bg_real = bg->type == EXPR_REAL;
    double bd = bg_int ? (double)bg->data.integer : bg_real ? bg->data.real : 0.0;
    int64_t bi = bg_int ? bg->data.integer : 0;
    bool bg_num = bg_int || bg_real;
    SpValKind kind = (dt == NDT_FLOAT64 || to_real) ? SPV_F64 : SPV_INT64;
    int64_t m = 0;
    if (dt == NDT_FLOAT64) {
        const double* b = nd->data.ndarray.data;
        if (bg_num) { for (int64_t i = 0; i < total; i++) m += b[i] != bd; }
        else m = total;
    } else {
        const int64_t* b = nd->data.ndarray.data;
        if (bg_int) { for (int64_t i = 0; i < total; i++) m += b[i] != bi; }
        else if (bg_real) { for (int64_t i = 0; i < total; i++) m += (double)b[i] != bd; }
        else m = total;
    }
    uint64_t* ck = malloc(sizeof(uint64_t) * (size_t)(m ? m : 1));
    void* vals = malloc(sizeof(double) * (size_t)(m ? m : 1));   /* 8 bytes either kind */
    if (!ck || !vals) { free(ck); free(vals); return NULL; }
    int64_t q = 0;
    for (int64_t i = 0; i < total; i++) {
        bool keep;
        double dv = 0; int64_t iv = 0;
        if (dt == NDT_FLOAT64) { dv = ((const double*)nd->data.ndarray.data)[i]; keep = !bg_num || dv != bd; }
        else {
            iv = ((const int64_t*)nd->data.ndarray.data)[i];
            keep = bg_int ? iv != bi : (!bg_real || (double)iv != bd);
        }
        if (!keep) continue;
        uint64_t row = rank >= 2 ? (uint64_t)(i / span) : 0;
        uint64_t key = rank >= 2 ? (uint64_t)(i % span) : (uint64_t)i;
        ck[q] = rank >= 2 ? (row << kbits) | key : key;
        if (kind == SPV_F64) ((double*)vals)[q] = dt == NDT_FLOAT64 ? dv : (double)iv;
        else ((int64_t*)vals)[q] = iv;
        q++;
    }
    SpCoo c;
    sp_coo_init(&c, rank, dims, kind);
    c.ck = ck;
    c.vals = vals;
    c.n = c.cap = m;
    c.ordered = true;
    return sp_to_expr(sp_coo_finish(&c, expr_copy((Expr*)bg), false));
}

/* SparseArray[list(, dims(, val))]. */
static Expr* from_dense(const Expr* call, const Expr* list, const Expr* dims_arg,
                        const Expr* bg_arg) {
    int rank;
    int64_t dims[SP_MAX_RANK];
    if (!dense_dims(list, &rank, dims) || rank == 0 || !dense_check(list, 0, rank, dims) ||
        sp_total(rank, dims) == 0) {
        sa_msg("rect", "List of rules or non-empty rectangular array expected at position 1 in %s.",
               call, NULL, NULL);
        return NULL;
    }
    DimSpec ds;
    if (!parse_dims(dims_arg, &ds)) return NULL;
    if (ds.given) {
        if (ds.scalar) dims_set_rank(&ds, rank);
        if (ds.rank != rank) {
            char* sd = show(dims_arg);
            mth_message("SparseArray", "drnk", "The requested dimensions, %s, have length "
                        "inconsistent with the tensor rank (%d) of the input.", sd, rank);
            free(sd);
            return NULL;
        }
    }
    LeafScan ls = { true, true, false };
    scan_leaves(list, 0, rank, &ls);
    Expr* bg = bg_arg ? expr_copy((Expr*)bg_arg)
                      : (ls.any && ls.all_inexact ? expr_new_real(0.0) : expr_new_integer(0));
    /* A machine-number array against an inexact background is an inexact
     * array: its Integer leaves become Reals. */
    bool to_real = bg_arg && bg->type == EXPR_REAL && ls.all_machine_real_or_int;

    bool own_dims = !ds.given;
    if (ds.given && ds.rank == rank) {
        own_dims = true;
        for (int k = 0; k < rank; k++) own_dims = own_dims && ds.dims[k] == dims[k];
    }
    if (list->type == EXPR_NDARRAY && own_dims) {
        Expr* fast = from_packed_dense(list, bg, to_real);
        if (fast) { expr_free(bg); return fast; }
    }
    SpValKind kind = SPV_EXPR;
    if (list->type == EXPR_NDARRAY) {
        NDType dt = list->data.ndarray.dtype;
        if (dt == NDT_FLOAT64 || dt == NDT_FLOAT32 || (dt == NDT_INT64 && to_real)) kind = SPV_F64;
        else if (dt == NDT_INT64) kind = SPV_INT64;
    }
    SpCoo c;
    sp_coo_init(&c, rank, ds.given ? ds.dims : dims, kind);
    c.ordered = true;
    int64_t coord[SP_MAX_RANK];
    if (!push_dense(&c, list, 0, coord, bg, to_real)) {
        sp_coo_free(&c);
        expr_free(bg);
        return NULL;
    }
    return sp_to_expr(sp_coo_finish(&c, bg, true));
}

/* ======================================================================== */
/* Rules                                                                    */
/* ======================================================================== */

typedef struct {
    const Expr* call;
    const Expr* rules_e;     /* the whole rule specification, for messages */
    DimSpec     ds;
    SpCoo       coo;
    bool        coo_open;
    Expr*       bg;          /* background (owned) */
    bool        bg_explicit;
    bool        done;        /* an all-matching rule set the background */
} RuleCtx;

static void msg_posd(RuleCtx* x, const Expr* rule) {
    Expr* d = dims_expr(x->ds.rank, x->ds.dims);
    sa_msg("posd", "The left-hand side of %s in %s is not a position or a pattern that will match "
           "the position of an element in an array with dimensions %s.", rule, x->rules_e, d);
    expr_free(d);
}

static void msg_valnl(const Expr* rule) {
    sa_msg("valnl", "The value specified by the rule %s should not be a List.", rule, NULL, NULL);
}

/* Resolve 1-based (possibly negative) coordinates against the dims into
 * 0-based ones; false when out of range. */
static bool resolve(const DimSpec* ds, const int64_t* c, int64_t* out) {
    for (int k = 0; k < ds->rank; k++) {
        int64_t i = c[k];
        if (i < 0) i = ds->dims[k] + 1 + i;
        if (i < 1 || i > ds->dims[k]) return false;
        out[k] = i - 1;
    }
    return true;
}

/* ---- pass 1: rank and inferred dims -------------------------------------- */

/* Array depth of a value list for Band (0 scalar, 1 vector, ...). */
static int value_depth(const Expr* v) {
    int d = 0;
    while (lv_is_list(v) && lv_len(v) > 0) {
        if (v->type == EXPR_NDARRAY) return d + v->data.ndarray.rank;
        d++;
        v = v->data.function.args[0];
    }
    return lv_is_list(v) ? d + 1 : d;
}

static bool band_vec(const Expr* e, int64_t* v, int rank) {
    int r;
    return read_position(e, v, &r) && r == rank;
}

/* Contribute a rule's literal extent to the inferred dims (max |coord|). */
static void infer_extent(int64_t* ext, int rank, const int64_t* c) {
    for (int k = 0; k < rank; k++) {
        int64_t a = c[k] < 0 ? -c[k] : c[k];
        if (a > ext[k]) ext[k] = a;
    }
}

/* Rank of the first rule that determines one (0 = none does). */
static int rule_rank(const Expr* rule) {
    const Expr* lhs = rule->data.function.args[0];
    const Expr* rhs = rule->data.function.args[1];
    int64_t c[SP_MAX_RANK];
    int r;
    (void)rhs;
    if (head_is_sym(lhs, SYM_Band) && lhs->data.function.arg_count >= 1)
        return read_position(lhs->data.function.args[0], c, &r) ? r : 0;
    if (read_position(lhs, c, &r)) return r;
    if (lv_is_list(lhs) && lv_len(lhs) > 0) {
        /* {p1, p2, ...}: the rank of the (List) positions inside */
        Expr* p0 = lv_item(lhs, 0);
        bool ok = p0 && p0->type != EXPR_INTEGER && read_position(p0, c, &r);
        if (p0) expr_free(p0);
        if (ok) return r;
    }
    if (head_is_sym(lhs, SYM_List) && has_pattern(lhs)) return (int)lhs->data.function.arg_count;
    return 0;
}

/* Infer dims from the literal positions and Bands. */
static void infer_dims(Expr* const* rules, size_t n, int rank, int64_t* ext) {
    for (int k = 0; k < rank; k++) ext[k] = 0;
    for (size_t i = 0; i < n; i++) {
        const Expr* lhs = rules[i]->data.function.args[0];
        const Expr* rhs = rules[i]->data.function.args[1];
        int64_t c[SP_MAX_RANK];
        int r;
        if (head_is_sym(lhs, SYM_Band)) {
            if (lhs->data.function.arg_count < 1 || !band_vec(lhs->data.function.args[0], c, rank))
                continue;
            int64_t step[SP_MAX_RANK];
            for (int k = 0; k < rank; k++) step[k] = 1;
            if (lhs->data.function.arg_count == 3) band_vec(lhs->data.function.args[2], step, rank);
            int d = value_depth(rhs);
            if (d == 0) {
                infer_extent(ext, rank, c);
            } else if (d == 1) {
                int64_t len = lv_len(rhs), e[SP_MAX_RANK];
                for (int k = 0; k < rank; k++) e[k] = c[k] + (len - 1) * step[k];
                infer_extent(ext, rank, c);
                infer_extent(ext, rank, e);
            } else if (d == rank || d == rank + 1) {
                /* one block, or a run of blocks down the diagonal */
                int64_t at[SP_MAX_RANK];
                memcpy(at, c, sizeof at);
                int64_t nb = d == rank ? 1 : lv_len(rhs);
                for (int64_t b = 0; b < nb; b++) {
                    Expr* blk = d == rank ? expr_copy((Expr*)rhs) : lv_item(rhs, b);
                    int br; int64_t bd[SP_MAX_RANK];
                    if (blk && dense_dims(blk, &br, bd) && br == rank) {
                        int64_t e[SP_MAX_RANK];
                        for (int k = 0; k < rank; k++) { e[k] = at[k] + bd[k] - 1; at[k] += bd[k]; }
                        infer_extent(ext, rank, e);
                    }
                    if (blk) expr_free(blk);
                }
            }
            continue;
        }
        if (read_position(lhs, c, &r)) {
            if (r == rank) infer_extent(ext, rank, c);
            continue;
        }
        if (lhs->type == EXPR_NDARRAY && lhs->data.ndarray.dtype == NDT_INT64 &&
            lhs->data.ndarray.rank == 2 && lhs->data.ndarray.dims[1] == rank) {
            /* a packed position matrix: straight off the buffer */
            const int64_t* P = lhs->data.ndarray.data;
            int64_t m = lhs->data.ndarray.dims[0];
            for (int64_t j = 0; j < m; j++) infer_extent(ext, rank, P + j * rank);
            continue;
        }
        if (lv_is_list(lhs) && lv_len(lhs) > 0) {
            int64_t m = lv_len(lhs);
            for (int64_t j = 0; j < m; j++) {
                Expr* p = lv_item(lhs, j);
                if (p && read_position(p, c, &r) && r == rank) infer_extent(ext, rank, c);
                if (p) expr_free(p);
            }
        }
    }
}

/* ---- pass 2: entries ----------------------------------------------------- */

static bool push_value(RuleCtx* x, const int64_t* c0, Expr* v, const Expr* rule) {
    if (is_array_value(v)) { msg_valnl(rule); expr_free(v); return false; }
    return sp_coo_push_expr(&x->coo, c0, v);
}

static bool put_position(RuleCtx* x, const Expr* pos, const Expr* val, const Expr* rule) {
    int64_t c[SP_MAX_RANK], c0[SP_MAX_RANK];
    int r;
    if (!read_position(pos, c, &r)) return false;
    if (r != x->ds.rank || !resolve(&x->ds, c, c0)) { msg_posd(x, rule); return false; }
    return push_value(x, c0, expr_copy((Expr*)val), rule);
}

static void msg_bndtr(const Expr* band_full, const Expr* vals, const Expr* rule) {
    sa_msg("bndtr", "There are insufficient positions available in the %s to fit all the values %s "
           "so some of the values will be omitted in the specification %s.", band_full, vals, rule);
}

/* Band[start(, end(, step))] -> v. */
static bool put_band(RuleCtx* x, const Expr* band, const Expr* val, const Expr* rule) {
    int rank = x->ds.rank;
    size_t na = band->data.function.arg_count;
    int64_t cur[SP_MAX_RANK], end[SP_MAX_RANK], step[SP_MAX_RANK];
    if (na < 1 || na > 3 || !band_vec(band->data.function.args[0], cur, rank)) {
        msg_posd(x, rule);
        return false;
    }
    for (int k = 0; k < rank; k++) {
        if (cur[k] < 0) cur[k] = x->ds.dims[k] + 1 + cur[k];
        step[k] = 1;
    }
    if (na == 3 && !band_vec(band->data.function.args[2], step, rank)) { msg_posd(x, rule); return false; }
    for (int k = 0; k < rank; k++) end[k] = step[k] > 0 ? x->ds.dims[k] : 1;
    if (na >= 2) {
        if (!band_vec(band->data.function.args[1], end, rank)) { msg_posd(x, rule); return false; }
        for (int k = 0; k < rank; k++) if (end[k] < 0) end[k] = x->ds.dims[k] + 1 + end[k];
    }
    int d = value_depth(val);
    if (d == rank || (d == rank + 1 && rank >= 2)) {
        /* Subblocks: one rank-r block at `start`, or a list of them placed
         * successively down the diagonal. */
        int64_t nb = d == rank ? 1 : lv_len(val);
        for (int64_t b = 0; b < nb; b++) {
            Expr* blk = d == rank ? expr_copy((Expr*)val) : lv_item(val, b);
            int br; int64_t bd[SP_MAX_RANK];
            if (!blk || !dense_dims(blk, &br, bd) || br != rank || !dense_check(blk, 0, rank, bd)) {
                if (blk) expr_free(blk);
                msg_valnl(rule);
                return false;
            }
            int64_t sub[SP_MAX_RANK];
            int64_t tot = sp_total(rank, bd);
            for (int k = 0; k < rank; k++) sub[k] = 0;
            for (int64_t t = 0; t < tot; t++) {
                int64_t c0[SP_MAX_RANK];
                bool inside = true;
                for (int k = 0; k < rank; k++) {
                    int64_t i = cur[k] + sub[k];
                    if (i < 1 || i > x->ds.dims[k]) inside = false;
                    c0[k] = i - 1;
                }
                if (inside) {
                    /* the element at `sub` of the block */
                    Expr* el = expr_copy(blk);
                    for (int k = 0; k < rank && el; k++) {
                        Expr* nx = lv_item(el, sub[k]);
                        expr_free(el);
                        el = nx;
                    }
                    if (!el || !sp_coo_push_expr(&x->coo, c0, el)) { expr_free(blk); return false; }
                }
                for (int k = rank - 1; k >= 0; k--) {
                    if (++sub[k] < bd[k]) break;
                    sub[k] = 0;
                }
            }
            for (int k = 0; k < rank; k++) cur[k] += bd[k];
            expr_free(blk);
        }
        return true;
    }
    if (d > 1) { msg_valnl(rule); return false; }
    bool vec = d == 1;
    int64_t n = vec ? lv_len(val) : INT64_MAX;
    for (int64_t j = 0; j < n; j++) {
        bool inside = true;
        for (int k = 0; k < rank && inside; k++) {
            if (cur[k] < 1 || cur[k] > x->ds.dims[k]) inside = false;
            else if (step[k] > 0 ? cur[k] > end[k] : cur[k] < end[k]) inside = false;
        }
        if (!inside) {
            if (vec) {
                Expr* full_args[3] = { NULL, NULL, NULL };
                int64_t st[SP_MAX_RANK];
                band_vec(band->data.function.args[0], st, rank);
                for (int k = 0; k < rank; k++) if (st[k] < 0) st[k] = x->ds.dims[k] + 1 + st[k];
                full_args[0] = dims_expr(rank, st);
                full_args[1] = dims_expr(rank, end);
                full_args[2] = dims_expr(rank, step);
                Expr* full = expr_new_function(expr_new_symbol(SYM_Band), full_args, 3);
                Expr* rr_args[2] = { expr_copy(full), expr_copy((Expr*)val) };
                Expr* rr = expr_new_function(expr_new_symbol(SYM_Rule), rr_args, 2);
                msg_bndtr(full, val, rr);
                expr_free(full); expr_free(rr);
            }
            break;
        }
        int64_t c0[SP_MAX_RANK];
        for (int k = 0; k < rank; k++) c0[k] = cur[k] - 1;
        Expr* v = vec ? lv_item(val, j) : expr_copy((Expr*)val);
        if (!push_value(x, c0, v, rule)) return false;
        for (int k = 0; k < rank; k++) cur[k] += step[k];
    }
    return true;
}

/* A pattern rule at one position: match {i, j, ...} (and, for a vector, the
 * bare index as well), then the substituted value, evaluated. */
static bool put_pattern_at(RuleCtx* x, const int64_t* c0, const Expr* lhs, const Expr* rhs,
                           const Expr* rule, bool* matched) {
    int rank = x->ds.rank;
    Expr* items[SP_MAX_RANK];
    for (int k = 0; k < rank; k++) items[k] = expr_new_integer(c0[k] + 1);
    Expr* pos = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)rank);
    MatchEnv* env = env_new();
    bool ok = match(pos, (Expr*)lhs, env);
    if (!ok && rank == 1) {
        env_free(env);
        env = env_new();
        Expr* bare = expr_new_integer(c0[0] + 1);
        ok = match(bare, (Expr*)lhs, env);
        expr_free(bare);
    }
    expr_free(pos);
    *matched = ok;
    if (!ok) { env_free(env); return true; }
    Expr* sub = replace_bindings((Expr*)rhs, env);
    env_free(env);
    Expr* v = evaluate(sub);
    expr_free(sub);
    return push_value(x, c0, v, rule);
}

static bool put_pattern(RuleCtx* x, const Expr* lhs, const Expr* rhs, const Expr* rule) {
    int rank = x->ds.rank;
    int64_t c0[SP_MAX_RANK];
    bool matched;
    if (is_diagonal_pattern(lhs, rank)) {
        int64_t m = x->ds.dims[0];
        for (int k = 1; k < rank; k++) if (x->ds.dims[k] < m) m = x->ds.dims[k];
        for (int64_t i = 0; i < m; i++) {
            for (int k = 0; k < rank; k++) c0[k] = i;
            if (!put_pattern_at(x, c0, lhs, rhs, rule, &matched)) return false;
        }
        return true;
    }
    int64_t total = sp_total(rank, x->ds.dims);
    if (total < 0) return false;
    for (int k = 0; k < rank; k++) c0[k] = 0;
    for (int64_t t = 0; t < total; t++) {
        if (!put_pattern_at(x, c0, lhs, rhs, rule, &matched)) return false;
        for (int k = rank - 1; k >= 0; k--) {
            if (++c0[k] < x->ds.dims[k]) break;
            c0[k] = 0;
        }
    }
    return true;
}

/* {p1, p2, ...} -> {v1, v2, ...} with a packed position matrix and a packed
 * value vector: straight from the buffers. Returns 1 handled, 0 not this
 * shape, -1 error (message emitted). */
static int put_packed_multi(RuleCtx* x, const Expr* lhs, const Expr* rhs, const Expr* rule) {
    if (lhs->type != EXPR_NDARRAY || lhs->data.ndarray.dtype != NDT_INT64) return 0;
    if (rhs->type != EXPR_NDARRAY || rhs->data.ndarray.rank != 1) return 0;
    int rank = x->ds.rank;
    int lr = lhs->data.ndarray.rank;
    int64_t n = lhs->data.ndarray.dims[0];
    int w = lr == 1 ? 1 : (int)lhs->data.ndarray.dims[1];
    if (lr > 2 || w != rank || rhs->data.ndarray.dims[0] != n) return 0;
    NDType vt = rhs->data.ndarray.dtype;
    if (vt != NDT_FLOAT64 && vt != NDT_INT64 && vt != NDT_COMPLEX64) return 0;
    if (x->coo.n > 0) return 0;      /* only as the sole (first) rule: kinds can't mix */
    int kbits;
    if (!sp_coo_composite(rank, x->ds.dims, &kbits)) return 0;
    SpValKind kind = vt == NDT_FLOAT64 ? SPV_F64 : vt == NDT_INT64 ? SPV_INT64 : SPV_C128;
    size_t vsz = sp_kind_size(kind);
    uint64_t* ck = malloc(sizeof(uint64_t) * (size_t)(n ? n : 1));
    void* vals = malloc(vsz * (size_t)(n ? n : 1));
    if (!ck || !vals) { free(ck); free(vals); return -1; }
    const int64_t* P = lhs->data.ndarray.data;
    const int64_t* D = x->ds.dims;
    int first = rank >= 2 ? 1 : 0;
    int64_t bad = -1;
    if (rank == 2) {
        /* the matrix case, unrolled: one pass validates and keys */
        int64_t d0 = D[0], d1 = D[1];
        for (int64_t i = 0; i < n; i++) {
            int64_t r = P[2 * i], cc = P[2 * i + 1];
            if (r < 0) r += d0 + 1;
            if (cc < 0) cc += d1 + 1;
            if (r < 1 || r > d0 || cc < 1 || cc > d1) { bad = i; break; }
            ck[i] = ((uint64_t)(r - 1) << kbits) | (uint64_t)(cc - 1);
        }
    } else {
        for (int64_t i = 0; i < n && bad < 0; i++) {
            uint64_t k = 0;
            for (int j = first; j < rank; j++) {
                int64_t q = P[i * w + j];
                if (q < 0) q += D[j] + 1;
                if (q < 1 || q > D[j]) { bad = i; break; }
                k = k * (uint64_t)D[j] + (uint64_t)(q - 1);
            }
            if (bad >= 0) break;
            if (first) {
                int64_t r = P[i * w];
                if (r < 0) r += D[0] + 1;
                if (r < 1 || r > D[0]) { bad = i; break; }
                k |= (uint64_t)(r - 1) << kbits;
            }
            ck[i] = k;
        }
    }
    if (bad >= 0) {
        free(ck); free(vals);
        Expr* ra[2] = { lv_item(lhs, bad), lv_item(rhs, bad) };
        Expr* one = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
        msg_posd(x, one);
        expr_free(one);
        (void)rule;
        return -1;
    }
    memcpy(vals, rhs->data.ndarray.data, vsz * (size_t)n);
    sp_coo_free(&x->coo);
    sp_coo_init(&x->coo, rank, x->ds.dims, kind);
    x->coo.ck = ck;
    x->coo.vals = vals;
    x->coo.n = x->coo.cap = n;
    return 1;
}

static bool apply_rule(RuleCtx* x, const Expr* rule, size_t index, size_t nrules) {
    const Expr* lhs = rule->data.function.args[0];
    const Expr* rhs = rule->data.function.args[1];
    bool delayed = head_is_sym(rule, SYM_RuleDelayed);
    int64_t c[SP_MAX_RANK];
    int r;

    if (head_is_sym(lhs, SYM_Band)) return put_band(x, lhs, rhs, rule);

    if (read_position(lhs, c, &r)) {
        if (r != x->ds.rank) { msg_posd(x, rule); return false; }
        if (delayed) {
            Expr* v = evaluate((Expr*)rhs);
            bool ok = put_position(x, lhs, v, rule);
            expr_free(v);
            return ok;
        }
        return put_position(x, lhs, rhs, rule);
    }

    /* A list of positions. */
    if (lv_is_list(lhs) && lv_len(lhs) > 0 && !has_pattern(lhs)) {
        int handled = (index == 0 && nrules == 1) ? put_packed_multi(x, lhs, rhs, rule) : 0;
        if (handled == 1) return true;
        if (handled < 0) return false;
        int64_t m = lv_len(lhs);
        bool many_vals = lv_is_list(rhs);
        if (many_vals && lv_len(rhs) != m) {
            sa_msg("rule1", "The rule %s is not of the form indices -> values.", rule, NULL, NULL);
            return false;
        }
        for (int64_t j = 0; j < m; j++) {
            Expr* p = lv_item(lhs, j);
            Expr* v = many_vals ? lv_item(rhs, j) : expr_copy((Expr*)rhs);
            bool ok = p && v;
            if (ok && !read_position(p, c, &r)) {
                sa_msg("pos", "The left-hand side of %s in %s is not a position or pattern that will "
                       "match a position.", rule, x->rules_e, NULL);
                ok = false;
            }
            if (ok) {
                int64_t c0[SP_MAX_RANK];
                if (r != x->ds.rank || !resolve(&x->ds, c, c0)) {
                    /* the message names the one offending p -> v */
                    Expr* pr_args[2] = { expr_copy(p), expr_copy(v) };
                    Expr* pr = expr_new_function(expr_new_symbol(SYM_Rule), pr_args, 2);
                    msg_posd(x, pr);
                    expr_free(pr);
                    ok = false;
                } else if (is_array_value(v)) {
                    msg_valnl(rule);
                    ok = false;
                } else {
                    ok = sp_coo_push_expr(&x->coo, c0, v);
                    v = NULL;   /* adopted */
                }
            }
            if (p) expr_free(p);
            if (v) expr_free(v);
            if (!ok) return false;
        }
        return true;
    }

    if (has_pattern(lhs)) {
        if (!delayed && !x->bg_explicit && matches_everything(lhs, x->ds.rank)) {
            /* {_, _} -> v names the background; later rules are shadowed. */
            if (is_array_value(rhs)) { msg_valnl(rule); return false; }
            expr_free(x->bg);
            x->bg = expr_copy((Expr*)rhs);
            x->done = true;
            return true;
        }
        return put_pattern(x, lhs, rhs, rule);
    }

    sa_msg("pos", "The left-hand side of %s in %s is not a position or pattern that will match a "
           "position.", rule, x->rules_e, NULL);
    return false;
}

/* SparseArray[rules(, dims(, val))]. */
static Expr* from_rules(const Expr* call, const Expr* spec, Expr* const* rules, size_t n,
                        const Expr* dims_arg, const Expr* bg_arg) {
    RuleCtx x;
    memset(&x, 0, sizeof x);
    x.call = call;
    x.rules_e = spec;
    if (!parse_dims(dims_arg, &x.ds)) return NULL;

    /* A literal position with a 0 index is no position: say so, naming
     * the rule (it fell through to ::exdims "from the positions {}"). */
    for (size_t i = 0; i < n; i++) {
        const Expr* lhs = rules[i]->data.function.args[0];
        if (lhs->type != EXPR_FUNCTION || lhs->data.function.head->type != EXPR_SYMBOL ||
            lhs->data.function.head->data.symbol.name != SYM_List) continue;
        bool ints = lhs->data.function.arg_count > 0, bad = false;
        for (size_t k = 0; ints && k < lhs->data.function.arg_count; k++) {
            const Expr* c = lhs->data.function.args[k];
            if (c->type != EXPR_INTEGER) ints = false;
            else if (c->data.integer == 0) bad = true;
        }
        if (ints && bad) {
            sa_msg("pos", "The left-hand side of %s in %s is not a position or pattern that will "
                   "match a position.", rules[i], spec, NULL);
            return NULL;
        }
    }

    /* Rank: from the dims, else from the first rule that determines one. */
    int rank = x.ds.rank;
    for (size_t i = 0; i < n && !rank; i++) rank = rule_rank(rules[i]);
    if (!rank && x.ds.scalar) rank = 1;
    if (!rank) {
        Expr* pos = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        sa_msg("exdims", "The dimensions cannot be determined from the positions %s.", pos, NULL, NULL);
        expr_free(pos);
        return NULL;
    }
    /* With the rank inferred, every literal position must have it (against
     * given dims a mismatch is SparseArray::posd, reported per rule). */
    for (size_t i = 0; i < n && !x.ds.given; i++) {
        const Expr* lhs = rules[i]->data.function.args[0];
        int64_t c[SP_MAX_RANK];
        int r;
        if (read_position(lhs, c, &r) && r != rank) {
            Expr* rk = expr_new_integer(rank);
            sa_msg("posr", "The left-hand side of %s in %s is not a position or a pattern that will "
                   "match the position of an element in an array with depth %s.", rules[i], spec, rk);
            expr_free(rk);
            return NULL;
        }
    }
    if (x.ds.scalar || !x.ds.given) {
        if (x.ds.scalar) {
            dims_set_rank(&x.ds, rank);
        } else {
            x.ds.rank = rank;
            infer_dims(rules, n, rank, x.ds.dims);
            /* Patterns alone (or nothing) determine no size. */
            if (sp_total(rank, x.ds.dims) == 0) {
                Expr* pos = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
                sa_msg("exdims", "The dimensions cannot be determined from the positions %s.", pos, NULL, NULL);
                expr_free(pos);
                return NULL;
            }
        }
    }
    if (bg_arg && is_array_value(bg_arg)) {
        sa_msg("ilist", "The implicit value %s in %s cannot be a list.", bg_arg, call, NULL);
        return NULL;
    }
    /* Given dimensions with a zero among them hold no element at all: the
     * answer is that empty array as an ordinary nested List ({} for {0},
     * {{}, {}, {}} for {3, 0}), whatever the rules say. */
    if (sp_total(rank, x.ds.dims) == 0) return empty_nested(rank, x.ds.dims);
    x.bg = bg_arg ? expr_copy((Expr*)bg_arg) : expr_new_integer(0);
    x.bg_explicit = bg_arg != NULL;
    sp_coo_init(&x.coo, rank, x.ds.dims, SPV_EXPR);
    for (size_t i = 0; i < n && !x.done; i++) {
        if (!apply_rule(&x, rules[i], i, n)) {
            sp_coo_free(&x.coo);
            expr_free(x.bg);
            return NULL;
        }
    }
    return sp_to_expr(sp_coo_finish(&x.coo, x.bg, true));
}

/* ======================================================================== */
/* Canonical form and existing arrays                                       */
/* ======================================================================== */

/* SparseArray[Automatic, dims, val, {1, {rowptr, colidx}, vals}]. */
static Expr* from_canonical(const Expr* call) {
    Expr* const* A = call->data.function.args;
    DimSpec ds;
    if (!parse_dims(A[1], &ds) || !ds.given || ds.scalar || ds.rank == 0) return NULL;
    const Expr* data = A[3];
    if (!head_is_sym(data, SYM_List) || data->data.function.arg_count != 3) return NULL;
    const Expr* idx = data->data.function.args[1];
    const Expr* vals = data->data.function.args[2];
    if (!head_is_sym(idx, SYM_List) || idx->data.function.arg_count != 2 || !lv_is_list(vals))
        return NULL;
    const Expr* rowptr = idx->data.function.args[0];
    const Expr* colidx = idx->data.function.args[1];
    if (!lv_is_list(rowptr) || !lv_is_list(colidx)) return NULL;
    int rank = ds.rank;
    int64_t nrows = rank >= 2 ? ds.dims[0] : 1;
    int width = rank >= 2 ? rank - 1 : 1;
    int64_t nnz = lv_len(vals);
    if (lv_len(rowptr) != nrows + 1 || lv_len(colidx) != nnz) return NULL;
    if (is_array_value(A[2])) return NULL;

    SpCoo c;
    sp_coo_init(&c, rank, ds.dims, SPV_EXPR);
    int64_t prev = 0;
    bool ok = lv_int_at(rowptr, 0, &prev) && prev == 0;
    for (int64_t r = 0; ok && r < nrows; r++) {
        int64_t hi;
        ok = lv_int_at(rowptr, r + 1, &hi) && hi >= prev && hi <= nnz;
        for (int64_t k = prev; ok && k < hi; k++) {
            Expr* t = lv_item(colidx, k);
            int64_t tc[SP_MAX_RANK], c0[SP_MAX_RANK];
            int tr;
            ok = t && read_position(t, tc, &tr) && tr == width && t->type != EXPR_INTEGER;
            if (t) expr_free(t);
            if (!ok) break;
            int j = 0;
            if (rank >= 2) c0[j++] = r;
            for (int q = 0; q < width; q++) {
                int64_t dim = ds.dims[(rank >= 2 ? 1 : 0) + q];
                if (tc[q] < 1 || tc[q] > dim) { ok = false; break; }
                c0[j++] = tc[q] - 1;
            }
            if (!ok) break;
            Expr* v = lv_item(vals, k);
            if (!v || is_array_value(v)) { if (v) expr_free(v); ok = false; break; }
            ok = sp_coo_push_expr(&c, c0, v);
        }
        prev = hi;
    }
    if (ok && prev != nnz) ok = false;
    if (!ok) { sp_coo_free(&c); return NULL; }
    /* The canonical form is taken at its word: explicit entries stay explicit
     * even when they equal the background. */
    return sp_to_expr(sp_coo_finish(&c, expr_copy(A[2]), false));
}

/* SparseArray[s(, dims(, val))]. */
static Expr* from_sparse(const Expr* call, const Expr* s, const Expr* dims_arg, const Expr* bg_arg) {
    const SparseArray* sa = sp_of(s);
    DimSpec ds;
    if (!parse_dims(dims_arg, &ds)) return NULL;
    if (ds.given) {
        if (ds.scalar) dims_set_rank(&ds, sa->rank);
        if (ds.rank != sa->rank) {
            char* sd = show(dims_arg);
            mth_message("SparseArray", "drnk", "The requested dimensions, %s, have length "
                        "inconsistent with the tensor rank (%d) of the input.", sd, sa->rank);
            free(sd);
            return NULL;
        }
    } else {
        ds.rank = sa->rank;
        memcpy(ds.dims, sa->dims, sizeof(int64_t) * (size_t)sa->rank);
    }
    if (bg_arg && is_array_value(bg_arg)) {
        sa_msg("ilist", "The implicit value %s in %s cannot be a list.", bg_arg, call, NULL);
        return NULL;
    }
    if (bg_arg && !sp_value_is_background(bg_arg, sa->bg)) {
        /* A new background: every old implicit position becomes explicit, so
         * this is the dense conversion of the (resized) array. */
        Expr* dense = sp_normal(sa);
        if (!dense) return NULL;
        Expr* d = dims_expr(ds.rank, ds.dims);
        Expr* out = from_dense(call, dense, d, bg_arg);
        expr_free(dense);
        expr_free(d);
        return out;
    }
    /* Same background: keep the explicit entries inside the new dims, minus
     * those that equal the background. */
    const SparsePattern* p = sa->pat;
    SpCoo c;
    sp_coo_init(&c, sa->rank, ds.dims, sa->kind);
    c.ordered = true;
    int first = sa->rank >= 2 ? 1 : 0;
    for (int64_t r = 0; r < p->nrows; r++) {
        if (first && r >= ds.dims[0]) break;
        for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) {
            int64_t c0[SP_MAX_RANK];
            bool inside = true;
            if (first) c0[0] = r;
            for (int j = 0; j < p->width; j++) {
                c0[first + j] = sp_col(p, k, j);
                if (c0[first + j] >= ds.dims[first + j]) inside = false;
            }
            if (!inside) continue;
            bool ok;
            switch (sa->kind) {
                case SPV_INT64: ok = sp_coo_push_i64(&c, c0, ((const int64_t*)sa->vals)[k]); break;
                case SPV_F64:   ok = sp_coo_push_f64(&c, c0, ((const double*)sa->vals)[k]); break;
                case SPV_C128:  ok = sp_coo_push_c128(&c, c0, ((const double*)sa->vals)[2 * k],
                                                     ((const double*)sa->vals)[2 * k + 1]); break;
                default:        ok = sp_coo_push_expr(&c, c0, sp_value_expr(sa, k)); break;
            }
            if (!ok) { sp_coo_free(&c); return NULL; }
        }
    }
    return sp_to_expr(sp_coo_finish(&c, expr_copy(bg_arg ? (Expr*)bg_arg : sa->bg), true));
}

/* ======================================================================== */
/* The builtin                                                              */
/* ======================================================================== */

static Expr* sparsearray_impl(Expr* res);

/* A Rule first argument may carry packed sides (the evaluator keeps them for
 * the buffer path, see eval_sparse_rule_arg in eval.c). When the constructor
 * declines, unpack them in place, so the call that comes to rest is an
 * ordinary tree with no packed node inside it. */
Expr* builtin_sparsearray(Expr* res) {
    Expr* out = sparsearray_impl(res);
    if (out || res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return out;
    Expr* r = res->data.function.args[0];
    if (is_rule(r) && (is_packed_list(r->data.function.args[0]) ||
                       is_packed_list(r->data.function.args[1]))) {
        /* pack_unpack is NULL for a side that is not a buffer: keep that
         * side as it is (only one of the two may have been packed). */
        Expr* a[2];
        for (int k = 0; k < 2; k++) {
            Expr* side = r->data.function.args[k];
            a[k] = is_packed_list(side) ? pack_unpack(side) : expr_copy(side);
        }
        Expr* plain = expr_new_function(expr_copy(r->data.function.head), a, 2);
        expr_free(r);
        res->data.function.args[0] = plain;
        expr_invalidate_hash(res);
    }
    return NULL;
}

static Expr* sparsearray_impl(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t na = res->data.function.arg_count;
    if (na < 1 || na > 4) return NULL;
    Expr* const* A = res->data.function.args;
    const Expr* spec = A[0];

    if (is_sym(spec, SYM_Automatic)) return na == 4 ? from_canonical(res) : NULL;
    if (na == 4) return NULL;
    const Expr* dims_arg = na >= 2 ? A[1] : NULL;
    const Expr* bg_arg = na == 3 ? A[2] : NULL;

    if (is_sparse(spec)) return from_sparse(res, spec, dims_arg, bg_arg);

    if (is_rule(spec)) return from_rules(res, spec, (Expr* const*)&A[0], 1, dims_arg, bg_arg);
    if (head_is_sym(spec, SYM_List)) {
        size_t n = spec->data.function.arg_count;
        bool all_rules = true;
        for (size_t i = 0; i < n && all_rules; i++) all_rules = is_rule(spec->data.function.args[i]);
        if (all_rules)
            return from_rules(res, spec, spec->data.function.args, n, dims_arg, bg_arg);
    }
    if (lv_is_list(spec) || spec->type == EXPR_LINKEDLIST) {
        if (bg_arg && is_array_value(bg_arg)) {
            sa_msg("ilist", "The implicit value %s in %s cannot be a list.", bg_arg, res, NULL);
            return NULL;
        }
        return from_dense(res, spec, dims_arg, bg_arg);
    }
    sa_msg("list", "List expected at position 1 in %s.", res, NULL, NULL);
    return NULL;
}
