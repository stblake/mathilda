// Part and Head

#include "eval.h"
#include "ndstruct.h"   /* ndstruct_delist_repack — packed-argument fallback */
#include "symtab.h"
#include "sym_names.h"
#include "assoc.h"
#include "ds_ops.h"
#include "pvec_ops.h"
#include "pvec.h"
#include "ndarray.h"
#include "pack.h"      /* pack_force -- lift a small gather source up to the buffer */
#include "part.h"
#include "sparse/sparse.h"
#include "message.h"
#include "print.h"     /* expr_to_string -- Set::ndarray */
#include "common.h"
#include "linkedlist.h"

static bool is_atomic(Expr* e);

/* True for a 2-argument Rule/RuleDelayed, i.e. a well-formed association entry.
 * is_association() only checks the head, so Association[1, 2] can reach us and
 * must never be dereferenced as though its arguments were rules. */
static bool part_is_entry(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.arg_count == 2;
}

/*
 * Resolve one slot of a Span spec against a container length.
 *
 * `which` names the slot (0 = start, 1 = end, 2 = step) and selects the
 * default that All stands for: the first element, the last element, and a
 * step of 1 respectively.  Negative indices and UpTo[n] are meaningful for
 * the bounds only, never for the step.  Returns false for a slot that is not
 * a valid index spec.
 */
static bool span_slot(const Expr* slot, int64_t len, int which, int64_t* out) {
    if (slot->type == EXPR_INTEGER) {
        int64_t v = slot->data.integer;
        if (which != 2 && v < 0) v = len + v + 1;
        *out = v;
        return true;
    }
    if (slot->type == EXPR_SYMBOL && slot->data.symbol.name == SYM_All) {
        *out = (which == 1) ? len : 1;
        return true;
    }
    if (which != 2 && slot->type == EXPR_FUNCTION &&
        slot->data.function.head->type == EXPR_SYMBOL &&
        slot->data.function.head->data.symbol.name == SYM_UpTo &&
        slot->data.function.arg_count == 1 &&
        slot->data.function.args[0]->type == EXPR_INTEGER) {
        int64_t v = slot->data.function.args[0]->data.integer;
        if (v > len) v = len;
        if (v < 0) v = len + v + 1;
        *out = v;
        return true;
    }
    return false;
}

/*
 * Normalize Span[a, b] / Span[a, b, step] against a container of `len`
 * elements into the 1-based index of the first selected element, the stride
 * between selections, and how many there are -- so a caller iterates
 *
 *     for (i = 0, k = *start; i < *count; i++, k += *step)
 *
 * with every k guaranteed to be in 1..len.
 *
 * Returns false for a malformed spec -- a slot that is not an index, or a
 * zero step.  A span that selects nothing is *not* malformed: it succeeds
 * with *count == 0.  Shared by Part and the position-path walker so a read
 * and a rewrite can never disagree about which elements a span names.
 */
static bool span_resolve(const Expr* span, int64_t len,
                         int64_t* start, int64_t* step, int64_t* count) {
    int64_t s = 1, e = len, st = 1;
    size_t argc = span->data.function.arg_count;

    if (argc >= 1 && !span_slot(span->data.function.args[0], len, 0, &s)) return false;
    if (argc >= 2 && !span_slot(span->data.function.args[1], len, 1, &e)) return false;
    if (argc >= 3) {
        if (!span_slot(span->data.function.args[2], len, 2, &st)) return false;
        if (st == 0) return false;
    }

    int64_t n = 0;
    if (st > 0) {
        if (s <= e && s >= 1 && e <= len) n = (e - s) / st + 1;
    } else {
        if (s >= e && s <= len && e >= 1) n = (s - e) / (-st) + 1;
    }
    if (n < 0) n = 0;

    *start = s; *step = st; *count = n;
    return true;
}
/*
 * Whether Span `span` lies within a container of `len` elements, the rule
 * Part, Take and the LinkedList reader share: both ends inside, or the empty
 * span one short of the start in the step's direction (3;;2 is {}, 3;;1 and
 * 2;;10 of a 5-list are errors).  When it does not, Part::take names the
 * positions and `whole` (audit Part#1) and false is returned; a malformed
 * spec is left to span_resolve.
 */
static bool span_check_take(const Expr* span, int64_t len, const Expr* whole) {
    int64_t s = 1, e = len, st = 1;
    size_t argc = span->data.function.arg_count;
    if (argc >= 1 && !span_slot(span->data.function.args[0], len, 0, &s)) return true;
    if (argc >= 2 && !span_slot(span->data.function.args[1], len, 1, &e)) return true;
    if (argc >= 3 && !span_slot(span->data.function.args[2], len, 2, &st)) return true;
    if (st == 0) return true;
    bool ok = (st > 0) ? (s >= 1 && s <= len + 1 && e >= 0 && e <= len && e >= s - 1)
                       : (s >= 0 && s <= len && e >= 1 && e <= len + 1 && e <= s + 1);
    if (!ok) {
        char* l = expr_to_string((Expr*)whole);
        mth_message("Part", "take", "Cannot take positions %lld through %lld in %s.",
                    (long long)s, (long long)e, l ? l : "");
        free(l);
    }
    return ok;
}

bool part_span_resolve(const Expr* span, int64_t len, int64_t* start, int64_t* step,
                       int64_t* count) {
    return span_resolve(span, len, start, step, count);
}

static Expr* expr_part_assign_rec(Expr* expr, Expr** indices, size_t nindices, Expr* rhs);
static bool ll_part_assign(Expr* llx, Expr** indices, size_t nindices, Expr* rhs);
static void ll_part_msg_partw(const char* head, Expr* ll, Expr* spec);

/* The piece of `rhs` that the k-th of `nsel` positions selected at one level
 * receives. Returns a new reference.
 *
 * Wolfram distributes a List right-hand side LEVEL BY LEVEL. At a subscript
 * that selects several positions (a list of positions, All, a Span), a List rhs
 * whose length equals the number of positions selected there is split across
 * them, and each piece carries on to the deeper subscripts. Any other rhs (an
 * atom, or a List of a different length) is broadcast whole to every position
 * and also carries on. So, on a 3x3 `m`:
 *     m[[{1,2},{1,2}]] = {{a,b},{c,d}}   rows get {a,b} then {c,d}, split again
 *     m[[{1,2},{1,2}]] = {a,b}           row 1 gets a twice, row 2 gets b twice
 *     m[[All, 2]] = {a,b}                every row's column 2 gets {a,b}
 * A single integer or key subscript passes rhs through untouched, which is why
 * `m[[2]] = {1,2,3}` stores the whole list; callers pass nsel == 0 for that.
 *
 * This replaced a single running counter over all the leaves, which distributed
 * {{a,b},{c,d}} as two whole lists into four cells and split a too-short list
 * partway (`v[[2;;4]] = {0,0}` gave {1, 0, 0, {0, 0}, ...}). */
static Expr* part_rhs_piece(Expr* rhs, size_t k, size_t nsel) {
    if (nsel == 0) return expr_copy(rhs);       /* a single-position subscript */
    if (rhs->type == EXPR_FUNCTION && rhs->data.function.head->type == EXPR_SYMBOL &&
        rhs->data.function.head->data.symbol.name == SYM_List &&
        rhs->data.function.arg_count == nsel)
        return expr_copy(rhs->data.function.args[k]);
    if (is_packed_list(rhs) && rhs->data.ndarray.dims[0] == (int64_t)nsel) {
        /* Row k of a packed rhs, kept on the buffer when it is itself a row. */
        Expr* ki = expr_new_integer((int64_t)k + 1);
        bool degrade = false;
        Expr* piece = ndarray_part(rhs, &ki, 1, &degrade);
        expr_free(ki);
        if (piece) return piece;
        Expr* plain = ndarray_to_nested_list(rhs);
        piece = expr_copy(plain->data.function.args[k]);
        expr_free(plain);
        return piece;
    }
    return expr_copy(rhs);
}

/* Assign through the value of association entry `entry`, recursing for any
 * remaining indices, and return the resulting entry. Takes ownership of `entry`
 * (the caller's array slot) and hands back a fresh one -- never an in-place
 * write, because expr_copy is a refcount bump and `entry` is the very node the
 * caller's association holds. Assigning into args[1] here is what made
 * `u2 = u1; u2[["x"]] = 9` corrupt `u1`. Returns `entry` unchanged when it is
 * not a two-argument rule (Association[1, 2] can exist). */
static Expr* assoc_entry_assigned(Expr* entry, Expr** rest, size_t nrest, Expr* rhs) {
    if (!part_is_entry(entry)) return entry;
    Expr* nv = expr_part_assign_rec(entry->data.function.args[1], rest, nrest, rhs);
    Expr* fresh = assoc_entry_with_value(entry, nv);   /* adopts nv */
    expr_free(entry);
    return fresh;
}

/* assoc_entry_assigned for the k-th of nsel selected entries. */
static Expr* assoc_entry_assigned_kth(Expr* entry, Expr** rest, size_t nrest,
                                      Expr* rhs, size_t k, size_t nsel) {
    Expr* piece = part_rhs_piece(rhs, k, nsel);
    Expr* out = assoc_entry_assigned(entry, rest, nrest, piece);
    expr_free(piece);
    return out;
}

/* Replace *slot (owned) by the assignment of the k-th of nsel rhs pieces through
 * the remaining subscripts. */
static void part_assign_slot(Expr** slot, Expr** rest, size_t nrest,
                             Expr* rhs, size_t k, size_t nsel) {
    Expr* piece = part_rhs_piece(rhs, k, nsel);
    Expr* replaced = expr_part_assign_rec(*slot, rest, nrest, piece);
    expr_free(piece);
    expr_free(*slot);
    *slot = replaced;
}


/* Assign into a dense NDArray on the buffer. Value semantics, as everywhere else
 * in Part assignment: a fresh NDArray comes back and `arr` is untouched.
 *
 * The whole Part spec vocabulary works — integer (including negative), All,
 * Span, a list of positions, and fewer subscripts than the rank — because the
 * write reuses ndarray_part_set, which selects positions through the same
 * per-axis selector ndarray_part reads through.
 *
 * A plain numeric List rhs is packed at the array's dtype first, so
 * `m[[2]] = {0, 0, 0, 0}` stays on the buffer (and the result stays packed).
 *
 * Returns NULL -- a DECLINE, never the array unchanged -- when the buffer cannot
 * hold the answer: a symbolic value, a shape that broadcasts or goes ragged, an
 * out-of-range subscript. The caller then answers through the List path. This
 * used to return the input untouched, which made `m[[2]] = {0, 0, 0, 0}` on any
 * packed matrix a silent no-op: ndarray_part_set takes only a scalar or an
 * NDArray rhs, and a literal List is neither. */
static Expr* nd_part_assign(Expr* arr, Expr** indices, size_t nindices, Expr* rhs) {
    NDType dt = arr->data.ndarray.dtype;
    Expr* packed_rhs = NULL;
    if (rhs->type == EXPR_FUNCTION) {
        packed_rhs = pack_force(expr_copy(rhs), true, dt);
        if (!is_ndarray(packed_rhs)) { expr_free(packed_rhs); return NULL; }
        rhs = packed_rhs;
    }
    /* An INDEPENDENT node with its own buffer.  Not expr_copy: that is a
     * refcount bump, so writing through it would reach every other holder of
     * the same array — including the symbol this assignment is meant to leave
     * untouched until the caller stores the result back. */
    size_t n = 1;
    for (int i = 0; i < arr->data.ndarray.rank; i++) n *= (size_t)arr->data.ndarray.dims[i];
    size_t bytes = ndt_elem_size(dt) * n;
    Expr* fresh = NULL;
    void* nb = malloc(bytes ? bytes : 1);
    if (nb) {
        memcpy(nb, arr->data.ndarray.data, bytes);
        fresh = expr_new_ndarray_like(arr, arr->data.ndarray.rank, arr->data.ndarray.dims,
                                      nb, dt);                          /* adopts nb */
        if (!fresh) free(nb);
    }
    if (fresh && !ndarray_part_set(fresh, (Expr* const*)indices, nindices, rhs)) {
        expr_free(fresh);
        fresh = NULL;
    }
    if (packed_rhs) expr_free(packed_rhs);
    return fresh;
}

/* True when every leaf of `rhs` can go into a `dt` buffer with its HEAD intact:
 * an exact Integer for int64, a machine Real for a float dtype.
 *
 * Deliberately stricter than "is this a machine number". Writing the exact
 * Integer 1 into a float64 buffer reads back as 1., so `p[[1]] = 1` would
 * answer {1., 2., 3.} where the same assignment on an ordinary list answers
 * {1, 2., 3.} -- an exact head on an inexact tail, which a uniform buffer
 * cannot express (src/numloop.h states the same rule for compiled loops).
 * Rationals, BigInts, arbitrary-precision reals, Complex and symbols all fail
 * here for the same reason. A packed rhs qualifies when its dtype is `dt`.
 *
 * A List right-hand side distributes its elements across the selected
 * positions, so every leaf has to qualify. */
static bool nd_rhs_fits_dtype(const Expr* rhs, NDType dt) {
    if (!rhs) return false;
    if (is_ndarray(rhs)) return rhs->data.ndarray.dtype == dt;
    if (rhs->type == EXPR_FUNCTION &&
        rhs->data.function.head->type == EXPR_SYMBOL &&
        rhs->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < rhs->data.function.arg_count; i++)
            if (!nd_rhs_fits_dtype(rhs->data.function.args[i], dt)) return false;
        return true;
    }
    return (dt == NDT_INT64) ? (rhs->type == EXPR_INTEGER)
                             : (rhs->type == EXPR_REAL);
}

/* Part assignment into an NDArray `arr` (packed list or visible NDArray).
 *
 * The buffer answers when it can (nd_part_assign). Otherwise the array goes
 * through the List path and is offered back to a buffer:
 *
 *   PACKED LIST   repacked only when the result is uniform with no change of
 *                 head (pack_repack_like), so `p[[1]] = 1` on {1., 2., 3.} keeps
 *                 the exact 1 and a ragged result stays a List -- as in WL.
 *                 A packed list unpacks up front when the value cannot be stored
 *                 with its head intact (nd_rhs_fits_dtype): the buffer path
 *                 would coerce it.
 *   VISIBLE       repacked at its own dtype, coercing -- the user asked for a
 *                 machine buffer and the dtype is part of the value. A result no
 *                 buffer holds (a symbol, a ragged shape) cannot be stored: the
 *                 array is left unchanged, and said so with Set::ndarray rather
 *                 than silently. */
static Expr* nd_part_assign_or_list(Expr* arr, Expr** indices, size_t nindices, Expr* rhs) {
    bool packed = is_packed_list(arr);
    NDType dt = arr->data.ndarray.dtype;
    /* A packed list is a List: l[[All]] = 0 on a matrix replaces each ROW by
     * 0 ({0, 0, 0}), it does not fill the buffer (that is a visible NDArray's
     * machine-buffer semantics). A scalar written above the element level
     * takes the List path. */
    bool above_elements = packed && nindices < (size_t)arr->data.ndarray.rank &&
                          !(rhs->type == EXPR_NDARRAY ||
                            (rhs->type == EXPR_FUNCTION && rhs->data.function.head->type == EXPR_SYMBOL &&
                             rhs->data.function.head->data.symbol.name == SYM_List));
    if (!above_elements && (!packed || nd_rhs_fits_dtype(rhs, dt))) {
        Expr* r = nd_part_assign(arr, indices, nindices, rhs);
        if (r) return r;
    }
    Expr* plain = ndarray_to_nested_list(arr);
    Expr* out = expr_part_assign_rec(plain, indices, nindices, rhs);
    expr_free(plain);
    if (packed) return pack_repack_like(arr, out);
    Expr* back = out ? ndarray_from_nested_list_like(arr, out, dt) : NULL;
    if (back) { expr_free(out); return back; }
    char* s = expr_to_string(rhs);
    mth_message("Set", "ndarray",
                "%s cannot be stored in an NDArray of type \"%s\"; the array is unchanged.",
                s, ndt_to_string(dt));
    free(s);
    if (out) expr_free(out);
    return expr_copy(arr);
}

/* Part assignment below the top level, per payload kind: the bodies of the
 * ds_ops part_assign slot.  Each returns the new value of `expr` (a new
 * reference), or NULL (after a message) to leave the Set unevaluated. */
bool ll_assign_in_place(Expr* llx, Expr** indices, size_t nindices, Expr* rhs) {
    return ll_part_assign(llx, indices, nindices, rhs);
}

Expr* ll_assign_row(Expr* expr, Expr** indices, size_t nindices, Expr* rhs) {
    /* A LinkedList reached below the top ({ll}[[1, 2]] = v): written in
     * place, as everywhere -- the enclosing immutable value is rebuilt around
     * the same object. A failed write leaves it unchanged (after its message). */
        ll_part_assign(expr, indices, nindices, rhs);
        return expr_copy(expr);
}

Expr* nd_assign_row(Expr* expr, Expr** indices, size_t nindices, Expr* rhs) {
    return nd_part_assign_or_list(expr, indices, nindices, rhs);
}

Expr* pvec_assign_row(Expr* expr, Expr** indices, size_t nindices, Expr* rhs) {
    /* A PVector is a value: the write builds a new PVector (path copy) and the
     * old one is untouched.  An integer first subscript is native; any other
     * spec goes through the List and back. */
        PVec* pv = expr->data.pvec;
        Expr* ix = indices[0];
        if (ix->type == EXPR_INTEGER) {
            int64_t len = pv_length(pv), k = ix->data.integer;
            if (k < 0) k = len + k + 1;
            if (k < 1 || k > len) { ll_part_msg_partw("Set", expr, ix); return expr_copy(expr); }   /* unchanged, as the other kinds */
            Expr* ne;
            if (nindices == 1) ne = expr_copy(rhs);
            else {
                Expr* cur = pv_get(pv, k - 1);
                ne = expr_part_assign_rec(cur, indices + 1, nindices - 1, rhs);
                expr_free(cur);
                if (!ne) return expr_copy(expr);
            }
            return expr_new_pvec(pv_set(pv_retain(pv), k - 1, ne));
        }
        Expr* l = pv_to_list(pv);
        Expr* r = expr_part_assign_rec(l, indices, nindices, rhs);
        expr_free(l);
        if (!r) return expr_copy(expr);
        Expr* back = pvec_make(r);
        if (back) { expr_free(r); return back; }
        return r;
}

Expr* sparse_assign_row(Expr* expr, Expr** indices, size_t nindices, Expr* rhs) {
        /* One element: rebuilt sparsely. Anything else (spans, All, a List
         * right-hand side) goes through the dense array and back, keeping the
         * background. A part that does not exist leaves the array as it was. */
        bool err = false;
        Expr* r = sp_part_assign(expr, indices, nindices, rhs, &err);
        if (r) return r;
        if (err) return expr_copy(expr);
        const SparseArray* sa = expr->data.sparse;
        Expr* dense = sp_normal(sa);
        if (!dense) return expr_copy(expr);
        /* A SparseArray right-hand side is spliced as the values it holds:
         * s[[All, 1]] = SparseArray[{p, q, r}] fills a column, it does not put
         * the vector object into every row. */
        Expr* rhs_dense = NULL;
        if (rhs->type == EXPR_SPARSE) {
            bool refused = false;
            rhs_dense = sp_gate_dense(rhs, &refused);
            if (!rhs_dense) { expr_free(dense); return expr_copy(expr); }
        }
        Expr* nd = expr_part_assign_rec(dense, indices, nindices, rhs_dense ? rhs_dense : rhs);
        expr_free(dense);
        if (rhs_dense) expr_free(rhs_dense);
        if (!nd) return expr_copy(expr);
        /* An assignment that breaks the array (a row replaced by an atom,
         * s[[All]] = 0) leaves an ordinary List, as on the dense form. */
        if (!sp_same_shape(nd, sa)) return nd;
        Expr** d = malloc(sizeof(Expr*) * (size_t)sa->rank);
        for (int k = 0; k < sa->rank; k++) d[k] = expr_new_integer(sa->dims[k]);
        Expr* dims = expr_new_function(expr_new_symbol(SYM_List), d, (size_t)sa->rank);
        free(d);
        Expr* args[3] = { nd, dims, expr_copy(sa->bg) };
        Expr* call = expr_new_function(expr_new_symbol(SYM_SparseArray), args, 3);
        Expr* out = evaluate(call);
        expr_free(call);
        return out;
}

static Expr* expr_part_assign_rec(Expr* expr, Expr** indices, size_t nindices, Expr* rhs) {
    if (nindices == 0) return expr_copy(rhs);

    /* A payload kind (LinkedList, NDArray, PVector, SparseArray): its row's
     * part_assign slot writes below this level. */
    {
        bool handled = false;
        Expr* r = ds_part_assign(expr, indices, nindices, rhs, &handled);
        if (handled) return r;
    }

    /* A subscript below an atom: Set::partd, and the value is unchanged
     * (audit Set#1: l[[1, 1]] = 9 on {1, 2, 3} was silent). */
    if (is_atomic(expr)) {
        mth_message("Set", "partd", "Part specification is longer than depth of object.");
        return expr_copy(expr);
    }

    Expr* idx = indices[0];
    Expr** rest = indices + 1;
    size_t nrest = nindices - 1;

    size_t len = expr->data.function.arg_count;
    Expr** new_args = NULL;
    if (len > 0) new_args = malloc(sizeof(Expr*) * len);
    for (size_t i = 0; i < len; i++) {
        new_args[i] = expr_copy(expr->data.function.args[i]);
    }
    Expr* new_head = expr_copy(expr->data.function.head);

    if (is_association(expr)) {
        /* Association assignment always targets entry *values*:
         *   a[[Key[k]]]/a[["s"]] = v  update (or append when absent),
         *   a[[i]] = v                positional value,
         *   a[[All]] = v              every value,
         *   a[[i;;j]] = v             spanned values,
         *   a[[{k1,k2,...}]] = v      the listed keys'/positions' values.
         * Non-key structural indices must NOT be appended as literal keys. */
        if (idx->type == EXPR_SYMBOL && idx->data.symbol.name == SYM_All) {
            for (size_t i = 0; i < len; i++)
                new_args[i] = assoc_entry_assigned_kth(new_args[i], rest, nrest, rhs, i, len);
        } else if (idx->type == EXPR_FUNCTION && idx->data.function.head->type == EXPR_SYMBOL &&
                   idx->data.function.head->data.symbol.name == SYM_Span) {
            int64_t start = 1, end = (int64_t)len, step = 1;
            size_t sa = idx->data.function.arg_count;
            if (sa >= 1 && idx->data.function.args[0]->type == EXPR_INTEGER) {
                start = idx->data.function.args[0]->data.integer; if (start < 0) start = (int64_t)len + start + 1;
            }
            if (sa >= 2 && idx->data.function.args[1]->type == EXPR_INTEGER) {
                end = idx->data.function.args[1]->data.integer; if (end < 0) end = (int64_t)len + end + 1;
            }
            if (sa >= 3 && idx->data.function.args[2]->type == EXPR_INTEGER) {
                step = idx->data.function.args[2]->data.integer;
            }
            if (step > 0) {
                size_t nsel = 0;
                for (int64_t i = start; i <= end && i >= 1 && i <= (int64_t)len; i += step) nsel++;
                size_t k = 0;
                for (int64_t i = start; i <= end && i >= 1 && i <= (int64_t)len; i += step)
                    new_args[i - 1] = assoc_entry_assigned_kth(new_args[i - 1], rest, nrest, rhs, k++, nsel);
            }
        } else if (idx->type == EXPR_FUNCTION && idx->data.function.head->type == EXPR_SYMBOL &&
                   idx->data.function.head->data.symbol.name == SYM_List) {
            size_t nsel = idx->data.function.arg_count;
            for (size_t j = 0; j < nsel; j++) {
                Expr* sub = idx->data.function.args[j];
                if (sub->type == EXPR_INTEGER) {
                    int64_t p = sub->data.integer; if (p < 0) p = (int64_t)len + p + 1;
                    if (p >= 1 && p <= (int64_t)len)
                        new_args[p - 1] = assoc_entry_assigned_kth(new_args[p - 1], rest, nrest, rhs, j, nsel);
                } else {
                    Expr* k = (sub->type == EXPR_FUNCTION && sub->data.function.head->type == EXPR_SYMBOL &&
                               sub->data.function.head->data.symbol.name == SYM_Key && sub->data.function.arg_count == 1)
                              ? sub->data.function.args[0] : sub;
                    for (size_t i = 0; i < len; i++)
                        if (new_args[i]->type == EXPR_FUNCTION && new_args[i]->data.function.arg_count == 2 &&
                            expr_eq(new_args[i]->data.function.args[0], k)) {
                            new_args[i] = assoc_entry_assigned_kth(new_args[i], rest, nrest, rhs, j, nsel); break;
                        }
                }
            }
        } else if (idx->type == EXPR_INTEGER) {
            int64_t pos = idx->data.integer;
            if (pos < 0) pos = (int64_t)len + pos + 1;
            if (pos >= 1 && pos <= (int64_t)len)
                new_args[pos - 1] = assoc_entry_assigned(new_args[pos - 1], rest, nrest, rhs);
        } else {
            /* Single key: Key[k] (unwrapped) or a literal key. */
            Expr* lookup_key = (idx->type == EXPR_FUNCTION && idx->data.function.head->type == EXPR_SYMBOL &&
                                idx->data.function.head->data.symbol.name == SYM_Key && idx->data.function.arg_count == 1)
                               ? idx->data.function.args[0] : idx;
            int64_t found = -1;
            for (size_t i = 0; i < len; i++) {
                Expr* rule = new_args[i];
                if (rule->type == EXPR_FUNCTION && rule->data.function.arg_count == 2 &&
                    expr_eq(rule->data.function.args[0], lookup_key)) { found = (int64_t)i; break; }
            }
            if (found >= 0) {
                new_args[found] = assoc_entry_assigned(new_args[found], rest, nrest, rhs);
            } else {
                /* The key is absent, so CREATE it.
                 *
                 * nrest == 0: this is the last subscript, and the recursive call
                 * returns the RHS value itself.
                 *
                 * nrest > 0: a nested path runs through a key that does not
                 * exist yet, so open an empty Association there and assign on
                 * through it -- `a = <||>; a["k", "s"] = 7` gives
                 * `<|"k" -> <|"s" -> 7|>|>`, as in Wolfram. This branch used to
                 * do nothing at all for nrest > 0, which made a deep write into
                 * a fresh Association a SILENT no-op: the Set was left
                 * unevaluated, `;` discarded it, and the association stayed
                 * empty with no message. Data loss with no diagnostic is worse
                 * than an error, and the memo-table idiom
                 * `$cache[key, "field"] = v` depends on the auto-vivification. */
                Expr* seed = (nrest == 0)
                    ? NULL
                    : expr_new_function(expr_new_symbol(SYM_Association), NULL, 0);
                Expr* nv = expr_part_assign_rec(seed ? seed : new_head, rest, nrest, rhs);
                if (seed) expr_free(seed);
                Expr* krule_args[2] = { expr_copy(lookup_key), nv };
                Expr* nrule = expr_new_function(expr_new_symbol(SYM_Rule), krule_args, 2);
                new_args = realloc(new_args, sizeof(Expr*) * (len + 1));
                new_args[len] = nrule;
                len++;
            }
        }
    } else if (idx->type == EXPR_INTEGER) {
        int64_t k = idx->data.integer;
        if (k == 0) {
            part_assign_slot(&new_head, rest, nrest, rhs, 0, 0);
        } else {
            if (k < 0) k = len + k + 1;
            if (k >= 1 && k <= (int64_t)len)
                part_assign_slot(&new_args[k - 1], rest, nrest, rhs, 0, 0);
            else {
                /* Set::partw, the value unchanged (audit Set#1). */
                char* ks = expr_to_string(idx);
                char* es = expr_to_string(expr);
                mth_message("Set", "partw", "Part %s of %s does not exist.", ks ? ks : "", es ? es : "");
                free(ks); free(es);
            }
        }
    } else if (idx->type == EXPR_SYMBOL && idx->data.symbol.name == SYM_All) {
        for (size_t i = 0; i < len; i++)
            part_assign_slot(&new_args[i], rest, nrest, rhs, i, len);
    } else if (idx->type == EXPR_FUNCTION && idx->data.function.head->data.symbol.name == SYM_Span) {
        int64_t start = 1, end = len, step = 1;
        size_t span_argc = idx->data.function.arg_count;
        if (span_argc >= 1) {
            Expr* a1 = idx->data.function.args[0];
            if (a1->type == EXPR_INTEGER) { start = a1->data.integer; if (start < 0) start = len + start + 1; }
            else if (a1->type == EXPR_SYMBOL && a1->data.symbol.name == SYM_All) start = 1;
            else if (a1->type == EXPR_FUNCTION && a1->data.function.head->data.symbol.name == SYM_UpTo && a1->data.function.arg_count == 1 && a1->data.function.args[0]->type == EXPR_INTEGER) {
                start = a1->data.function.args[0]->data.integer;
                if (start > (int64_t)len) start = len;
                if (start < 0) start = len + start + 1;
            }
        }
        if (span_argc >= 2) {
            Expr* a2 = idx->data.function.args[1];
            if (a2->type == EXPR_INTEGER) { end = a2->data.integer; if (end < 0) end = len + end + 1; }
            else if (a2->type == EXPR_SYMBOL && a2->data.symbol.name == SYM_All) end = len;
            else if (a2->type == EXPR_FUNCTION && a2->data.function.head->data.symbol.name == SYM_UpTo && a2->data.function.arg_count == 1 && a2->data.function.args[0]->type == EXPR_INTEGER) {
                end = a2->data.function.args[0]->data.integer;
                if (end > (int64_t)len) end = len;
                if (end < 0) end = len + end + 1;
            }
        }
        if (span_argc >= 3) {
            Expr* a3 = idx->data.function.args[2];
            if (a3->type == EXPR_INTEGER) step = a3->data.integer;
            else if (a3->type == EXPR_SYMBOL && a3->data.symbol.name == SYM_All) step = 1;
        }

        /* Count first: the rhs splits only when its length equals the count. */
        size_t nsel = 0, k = 0;
        if (step > 0) {
            for (int64_t i = start; i <= end && i >= 1 && i <= (int64_t)len; i += step) nsel++;
            for (int64_t i = start; i <= end && i >= 1 && i <= (int64_t)len; i += step)
                part_assign_slot(&new_args[i - 1], rest, nrest, rhs, k++, nsel);
        } else if (step < 0) {
            for (int64_t i = start; i >= end && i >= 1 && i <= (int64_t)len; i += step) nsel++;
            for (int64_t i = start; i >= end && i >= 1 && i <= (int64_t)len; i += step)
                part_assign_slot(&new_args[i - 1], rest, nrest, rhs, k++, nsel);
        }
    } else if (idx->type == EXPR_FUNCTION && idx->data.function.head->data.symbol.name == SYM_List) {
        size_t nsel = idx->data.function.arg_count;
        for (size_t i = 0; i < nsel; i++) {
            Expr* sub_idx = idx->data.function.args[i];
            if (sub_idx->type == EXPR_INTEGER) {
                int64_t k = sub_idx->data.integer;
                if (k == 0) {
                    part_assign_slot(&new_head, rest, nrest, rhs, i, nsel);
                } else {
                    if (k < 0) k = len + k + 1;
                    if (k >= 1 && k <= (int64_t)len)
                        part_assign_slot(&new_args[k - 1], rest, nrest, rhs, i, nsel);
                }
            }
        }
    }

    Expr* new_func = expr_new_function(new_head, new_args, len);
    if (new_args) free(new_args);
    return new_func;
}

/* Write ONE element of an exclusively owned packed array in place: every index
 * an in-range machine integer (negatives count from the end), exactly `rank`
 * of them, and a scalar of the buffer's own type (an exact Integer into int64, a
 * machine Real into float64) -- the cases that cannot fail halfway, so a refusal
 * leaves the array untouched and the caller falls back to the copying path.
 * This is what makes `Do[p[[i]] = v, {i, n}]` O(n) instead of O(n^2): the copy
 * path duplicated the whole buffer on every assignment. */
static bool nd_set_element_inplace(Expr* arr, Expr* const* idx, size_t n, const Expr* rhs) {
    if (arr->type != EXPR_NDARRAY || n != (size_t)arr->data.ndarray.rank) return false;
    NDType dt = arr->data.ndarray.dtype;
    if (dt == NDT_INT64) { if (rhs->type != EXPR_INTEGER) return false; }
    else if (dt == NDT_FLOAT64) { if (rhs->type != EXPR_REAL) return false; }
    else return false;
    size_t off = 0;
    for (size_t k = 0; k < n; k++) {
        if (idx[k]->type != EXPR_INTEGER) return false;
        int64_t d = arr->data.ndarray.dims[k], i = idx[k]->data.integer;
        if (i < 0) i += d + 1;
        if (i < 1 || i > d) return false;
        off = off * (size_t)d + (size_t)(i - 1);
    }
    if (dt == NDT_INT64) ((int64_t*)arr->data.ndarray.data)[off] = rhs->data.integer;
    else ((double*)arr->data.ndarray.data)[off] = rhs->data.real;
    expr_invalidate_hash(arr);
    return true;
}

Expr* expr_part_assign(Expr* lhs, Expr* rhs) {
    if (lhs->type != EXPR_FUNCTION || lhs->data.function.head->type != EXPR_SYMBOL || lhs->data.function.head->data.symbol.name != SYM_Part) return NULL;
    if (lhs->data.function.arg_count < 2) return NULL;

    Expr* sym = lhs->data.function.args[0];
    /* A nested target x[[i]][[j]] = v is x[[i, j]] = v when every inner index
     * is an integer (then composing the Parts is concatenating the specs; an
     * inner All or Span would not compose that way and stays unsupported).
     * The inner indices were held with the target, so they are evaluated
     * here: v[[k]][[i]] = ... with k = 3 updates v[[3, i]]. */
    if (sym->type == EXPR_FUNCTION && sym->data.function.head->type == EXPR_SYMBOL &&
        sym->data.function.head->data.symbol.name == SYM_Part &&
        sym->data.function.arg_count >= 2) {
        size_t nin = sym->data.function.arg_count - 1;
        size_t nout = lhs->data.function.arg_count - 1;
        Expr** idx = malloc(sizeof(Expr*) * (nin + nout + 1));
        bool ok = true;
        for (size_t i = 0; i < nin; i++) {
            idx[i + 1] = evaluate(sym->data.function.args[i + 1]);
            if (idx[i + 1]->type != EXPR_INTEGER) ok = false;
        }
        if (!ok) {
            for (size_t i = 0; i < nin; i++) expr_free(idx[i + 1]);
            free(idx);
            return NULL;
        }
        for (size_t i = 0; i < nout; i++) idx[nin + 1 + i] = expr_copy(lhs->data.function.args[i + 1]);
        idx[0] = expr_copy(sym->data.function.args[0]);
        Expr* flat = expr_new_function(expr_new_symbol(SYM_Part), idx, nin + nout + 1);
        free(idx);
        Expr* r = expr_part_assign(flat, rhs);
        expr_free(flat);
        return r;
    }
    /* The target is the LinkedList object itself -- what a function body
     * sees after its parameter is substituted (setFirst[l_, v_] :=
     * (l[[1]] = v)): reference semantics make that write the caller's list. */
    if (ds_has_in_place(sym->type))
        return ds_assign_in_place(sym, lhs->data.function.args + 1, lhs->data.function.arg_count - 1, rhs)
               ? expr_copy(sym) : NULL;
    if (sym->type != EXPR_SYMBOL) return NULL;

    Expr* current_val = symtab_get_own_values(sym->data.symbol.name) ? evaluate(sym) : NULL;
    if (!current_val) return NULL;

    /* A variable holding a LinkedList: write the shared object in place; the
     * OwnValue already points at it and is not rebound. */
    if (ds_has_in_place(current_val->type)) {
        bool ok = ds_assign_in_place(current_val, lhs->data.function.args + 1,
                                 lhs->data.function.arg_count - 1, rhs);
        if (!ok) { expr_free(current_val); return NULL; }
        return current_val;
    }

    /* v[[i]] = x on a PVector this variable owns exclusively: pv_set edits the
     * payload in place (O(1) when its nodes are unshared, else one path copy),
     * and the node is re-pointed at the result.  The node being exclusively
     * owned is what makes this invisible to every other holder. */
    if (current_val->type == EXPR_PVEC && lhs->data.function.arg_count == 2 &&
        lhs->data.function.args[1]->type == EXPR_INTEGER &&
        symtab_value_exclusively_owned(sym, current_val)) {
        int64_t len = pv_length(current_val->data.pvec), k = lhs->data.function.args[1]->data.integer;
        if (k < 0) k = len + k + 1;
        if (k >= 1 && k <= len) {
            current_val->data.pvec = pv_set(current_val->data.pvec, k - 1, expr_copy(rhs));
            expr_invalidate_hash(current_val);
            eval_clock_bump_binding(symtab_get_def(sym->data.symbol.name));
            return current_val;
        }
    }

    /* p[[i, ...]] = scalar on a packed array this variable owns exclusively:
     * write the element in place (see nd_set_element_inplace). */
    if (current_val->type == EXPR_NDARRAY && symtab_value_exclusively_owned(sym, current_val) &&
        nd_set_element_inplace(current_val, lhs->data.function.args + 1,
                               lhs->data.function.arg_count - 1, rhs)) {
        eval_clock_bump_binding(symtab_get_def(sym->data.symbol.name));
        return current_val;
    }

    /* a[[Key[k]]] = v / a["k"] = v on an association this variable owns
     * exclusively: upsert in place (O(1) amortised) instead of rebuilding the
     * whole node.  Positional/All/Span specs and nested paths take the general
     * path below.  See assoc_try_inplace_update for the ownership proof. */
    if (is_association(current_val) && lhs->data.function.arg_count == 2) {
        Expr* ix = lhs->data.function.args[1];
        Expr* key = NULL;
        if (ix->type == EXPR_STRING) key = ix;
        else if (ix->type == EXPR_FUNCTION && ix->data.function.head->type == EXPR_SYMBOL &&
                 ix->data.function.head->data.symbol.name == SYM_Key &&
                 ix->data.function.arg_count == 1) key = ix->data.function.args[0];
        if (key) {
            Expr* ra[2] = { expr_copy(key), expr_copy(rhs) };
            Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
            bool done = assoc_try_inplace_update(sym, current_val, rule, false);
            expr_free(rule);
            if (done) return current_val;
        }
    }

    Expr* new_val = expr_part_assign_rec(current_val, lhs->data.function.args + 1,
                                         lhs->data.function.arg_count - 1, rhs);
    expr_free(current_val);

    if (new_val) {
        symtab_add_own_value(sym->data.symbol.name, sym, new_val);
        return new_val;
    }
    return NULL;
}

static bool is_atomic(Expr* e) {
    if (!e) return true;
    if (e->type != EXPR_FUNCTION) return true;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Complex || h == SYM_Rational) return true;
    }
    return false;
}

Expr* expr_head(Expr* e) {
    if (!e) return NULL;
    /* Payload kinds (packed array, SparseArray, LinkedList, ...) name their own
     * Head through the ds_ops registry.  A packed list IS a List and must say so:
     * is_atomic() above deliberately still reports true for it (the 13 sites that
     * consult it go on to read data.function.arg_count, which aliases the ndarray
     * payload's `data` pointer, so flipping that predicate would turn a heap
     * pointer into an element count); presentation is decided in the NDArray row
     * and in the AtomQ builtin, not in the internal walker guard. */
    const DsOps* ops = ds_ops(e->type);
    if (ops) return ops->head(e);
    switch (e->type) {
        case EXPR_INTEGER: return expr_new_symbol(SYM_Integer);
        case EXPR_BIGINT: return expr_new_symbol(SYM_Integer);
        case EXPR_REAL: return expr_new_symbol(SYM_Real);
#ifdef USE_MPFR
        case EXPR_MPFR: return expr_new_symbol(SYM_Real);
#endif
        case EXPR_SYMBOL: return expr_new_symbol(SYM_Symbol);
        case EXPR_STRING: return expr_new_symbol(SYM_String);
        case EXPR_COMPILED: return expr_new_symbol(SYM_CompiledFunction);
        case EXPR_FUNCTION: return expr_copy(e->data.function.head);
        default: return NULL;
    }
}

Expr* builtin_head(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;

    Expr* head = expr_head(res->data.function.args[0]);
    if (!head) return NULL;
    if (argc == 1) return head;

    /* Head[expr, h] -> h[Head[expr]]. The wrapper h itself is left to the
     * evaluator to evaluate further (e.g. Head[{a,b}, f] -> f[List]). */
    Expr** wargs = malloc(sizeof(Expr*));
    wargs[0] = head;
    Expr* wrapper = expr_copy(res->data.function.args[1]);
    Expr* out = expr_new_function(wrapper, wargs, 1);
    free(wargs);
    return out;
}

/* Part of an association at one index level; defined below. */
static Expr* assoc_part(Expr* assoc, Expr* idx, Expr** rest, size_t nrest);

/* ------------------------------------------------------------------------
 * LinkedList Part (LINKED_LIST_DESIGN.md §3.2)
 *
 * The first subscript addresses the list's elements in O(1) each:
 *   - an integer (negative counts from the end) gives that element, and the
 *     remaining subscripts are its own Part;
 *   - All, a Span or a List of integers gives a NEW LinkedList of the
 *     selected elements, each with the remaining subscripts applied -- the
 *     closure rule: what is a List for a List is a LinkedList here.
 * [[0]] (the head) never reaches this; expr_part's integer branch answers
 * LinkedList through expr_head.
 *
 * Errors are Part::partw (an integer out of range) and Part::take (a span
 * running off the end), naming the list as LinkedList[<n>], and the Part
 * stays unevaluated. A deeper subscript that does not exist in an element
 * declines exactly as the element's own Part does.
 * ---------------------------------------------------------------------- */

static void ll_part_msg_partw(const char* head, Expr* ll, Expr* spec) {
    char* s = expr_to_string(spec);
    char* l = expr_to_string(ll);
    mth_message(head, "partw", "Part %s of %s does not exist.", s ? s : "", l ? l : "");
    free(s); free(l);
}

/* Part of one element at the remaining subscripts. The element may itself be
 * a packed array (a pushed row), which expr_part treats as atomic. */
static Expr* ll_elem_part(Expr* e, Expr** rest, size_t nrest) {
    if (nrest == 0) return e;                                  /* adopt */
    Expr* r;
    if (is_ndarray(e) && !(rest[0]->type == EXPR_INTEGER && rest[0]->data.integer == 0)) {
        bool degrade = false;
        r = ndarray_part(e, rest, nrest, &degrade);
        if (!r && degrade) {
            Expr* plain = pack_unpack(e);
            r = plain ? expr_part(plain, rest, nrest) : NULL;
            expr_free(plain);
        }
    } else {
        r = expr_part(e, rest, nrest);
    }
    expr_free(e);
    return r;
}

/* The 0-based positions a selecting subscript names, or NULL (with a message
 * when the spec is well formed but out of range). *out_n receives the count. */
static int64_t* ll_part_positions(const char* head, Expr* llx, Expr* idx, int64_t len,
                                  int64_t* out_n) {
    int64_t* pos = NULL;
    if (idx->type == EXPR_SYMBOL && idx->data.symbol.name == SYM_All) {
        pos = malloc(sizeof(int64_t) * (size_t)(len ? len : 1));
        for (int64_t i = 0; i < len; i++) pos[i] = i;
        *out_n = len;
        return pos;
    }
    if (idx->type != EXPR_FUNCTION || idx->data.function.head->type != EXPR_SYMBOL)
        return NULL;
    const char* h = idx->data.function.head->data.symbol.name;
    if (h == SYM_Span) {
        size_t argc = idx->data.function.arg_count;
        if (argc < 1 || argc > 3) return NULL;
        int64_t s = 1, e = len, st = 1;
        if (!span_slot(idx->data.function.args[0], len, 0, &s)) return NULL;
        if (argc >= 2 && !span_slot(idx->data.function.args[1], len, 1, &e)) return NULL;
        if (argc == 3 && (!span_slot(idx->data.function.args[2], len, 2, &st) || st == 0))
            return NULL;
        /* Valid: both ends inside the list, or the empty span one short of
         * the start in the step's direction (3;;2 is {}, 3;;1 is an error). */
        bool ok;
        int64_t n = 0;
        if (st > 0) {
            ok = s >= 1 && s <= len + 1 && e >= 0 && e <= len && e >= s - 1;
            if (ok && e >= s) n = (e - s) / st + 1;
        } else {
            ok = s >= 0 && s <= len && e >= 1 && e <= len + 1 && e <= s + 1;
            if (ok && e <= s) n = (s - e) / (-st) + 1;
        }
        if (!ok) {
            char* l = expr_to_string(llx);
            mth_message(head, "take", "Cannot take positions %lld through %lld in %s.",
                        (long long)s, (long long)e, l ? l : "");
            free(l);
            return NULL;
        }
        pos = malloc(sizeof(int64_t) * (size_t)(n ? n : 1));
        for (int64_t i = 0, k = s; i < n; i++, k += st) pos[i] = k - 1;
        *out_n = n;
        return pos;
    }
    if (h == SYM_List) {
        size_t n = idx->data.function.arg_count;
        pos = malloc(sizeof(int64_t) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) {
            Expr* q = idx->data.function.args[i];
            if (q->type != EXPR_INTEGER) { free(pos); return NULL; }
            int64_t k = q->data.integer;
            if (k < 0) k = len + k + 1;
            if (k < 1 || k > len) { ll_part_msg_partw(head, llx, q); free(pos); return NULL; }
            pos[i] = k - 1;
        }
        *out_n = (int64_t)n;
        return pos;
    }
    return NULL;
}

Expr* ll_part_expr(Expr* llx, Expr** indices, size_t nindices);
static Expr* ll_part(Expr* llx, Expr** indices, size_t nindices) {
    LinkedList* ll = llx->data.llist;
    int64_t len = ll_length(ll);
    Expr* idx = indices[0];
    Expr** rest = indices + 1;
    size_t nrest = nindices - 1;

    if (idx->type == EXPR_INTEGER) {
        int64_t k = idx->data.integer;
        if (k < 0) k = len + k + 1;
        if (k < 1 || k > len) { ll_part_msg_partw("Part", llx, idx); return NULL; }
        return ll_elem_part(ll_get(ll, k - 1), rest, nrest);
    }

    int64_t n = 0;
    int64_t* pos = ll_part_positions("Part", llx, idx, len, &n);
    if (!pos) return NULL;
    LinkedList* out = ll_new();
    for (int64_t i = 0; i < n; i++) {
        if (nrest == 0) {                         /* typed slots copy unboxed */
            ll_push_back_from(out, ll, pos[i]);
            continue;
        }
        Expr* r = ll_elem_part(ll_get(ll, pos[i]), rest, nrest);
        if (!r) { ll_release(out); free(pos); return NULL; }
        ll_push_back(out, r);
    }
    free(pos);
    return expr_new_linkedlist(out);
}

Expr* ll_part_expr(Expr* llx, Expr** indices, size_t nindices) { return ll_part(llx, indices, nindices); }

/* PVector Part: ll_part's reading rules over a persistent vector.  An integer
 * subscript is one element; All / Span / List select a new PVector (sharing
 * nothing with the source -- a selection is a fresh, short vector). */
Expr* pvec_part(Expr* pvx, Expr** indices, size_t nindices) {
    PVec* pv = pvx->data.pvec;
    int64_t len = pv_length(pv);
    Expr* idx = indices[0];
    Expr** rest = indices + 1;
    size_t nrest = nindices - 1;

    if (idx->type == EXPR_INTEGER) {
        int64_t k = idx->data.integer;
        if (k < 0) k = len + k + 1;
        if (k < 1 || k > len) { ll_part_msg_partw("Part", pvx, idx); return NULL; }
        return ll_elem_part(pv_get(pv, k - 1), rest, nrest);
    }
    int64_t n = 0;
    int64_t* pos = ll_part_positions("Part", pvx, idx, len, &n);
    if (!pos) return NULL;
    PVec* out = pv_new();
    for (int64_t i = 0; i < n; i++) {
        Expr* r = nrest == 0 ? pv_get(pv, pos[i]) : ll_elem_part(pv_get(pv, pos[i]), rest, nrest);
        if (!r) { pv_release(out); free(pos); return NULL; }
        out = pv_push(out, r);
    }
    free(pos);
    return expr_new_pvec(out);
}

/* ------------------------------------------------------------------------
 * LinkedList Part assignment, IN PLACE (LINKED_LIST_DESIGN.md §3.3)
 *
 * ll[[spec...]] = v writes the selected slots of the list itself -- the one
 * mutable object in the system, shared by every variable that holds it -- in
 * O(1) per written slot. The first subscript takes the full spec: an integer
 * (negative from the end), All, a Span or a List of integers; a List or
 * LinkedList right-hand side whose length equals the number of selected
 * positions is distributed over them, anything else is broadcast, exactly as
 * for a List (part_rhs_piece). Further subscripts rewrite the element at each
 * position as a List assignment would (copy-on-write for an immutable
 * element; a LinkedList element is itself written in place).
 *
 * All-or-nothing: every new value is computed and cycle-checked before any
 * slot changes, so Set::partw / Set::take / Set::llcyc leave the list as it
 * was and the Set unevaluated.
 * ---------------------------------------------------------------------- */

static Expr* expr_part_assign_rec(Expr* expr, Expr** indices, size_t nindices, Expr* rhs);

/* The piece of `rhs` the k-th of nsel positions receives: as part_rhs_piece,
 * plus a LinkedList right-hand side of matching length. */
static Expr* ll_rhs_piece(Expr* rhs, size_t k, size_t nsel) {
    if (nsel > 0 && rhs->type == EXPR_LINKEDLIST && ll_length(rhs->data.llist) == (int64_t)nsel)
        return ll_get(rhs->data.llist, (int64_t)k);
    return part_rhs_piece(rhs, k, nsel);
}

/* Raw buffer-to-slot write of a packed rank-1 right-hand side into a typed
 * list of the same kind (ll[[a ;; b]] = packedvec): no Expr per element. */
static bool ll_assign_packed(LinkedList* ll, const int64_t* pos, int64_t n, Expr* rhs) {
    if (!is_packed_list(rhs) || rhs->data.ndarray.rank != 1 || rhs->data.ndarray.dims[0] != n)
        return false;
    NDType dt = rhs->data.ndarray.dtype;
    if (dt == NDT_INT64 && ll_kind(ll) == LL_I64) {
        const int64_t* d = (const int64_t*)rhs->data.ndarray.data;
        for (int64_t j = 0; j < n; j++) ll_set_i64(ll, pos[j], d[j]);
        return true;
    }
    if (dt == NDT_FLOAT64 && ll_kind(ll) == LL_REAL) {
        const double* d = (const double*)rhs->data.ndarray.data;
        for (int64_t j = 0; j < n; j++) ll_set_real(ll, pos[j], d[j]);
        return true;
    }
    return false;
}

static bool ll_part_assign(Expr* llx, Expr** indices, size_t nindices, Expr* rhs) {
    LinkedList* ll = llx->data.llist;
    int64_t len = ll_length(ll);
    Expr* idx = indices[0];
    Expr** rest = indices + 1;
    size_t nrest = nindices - 1;

    int64_t one, n, *pos;
    size_t nsel;
    if (idx->type == EXPR_INTEGER) {
        int64_t k = idx->data.integer;
        if (k < 0) k = len + k + 1;
        if (k < 1 || k > len) { ll_part_msg_partw("Set", llx, idx); return false; }
        one = k - 1;
        pos = &one;
        n = 1;
        nsel = 0;                         /* a single position takes rhs whole */
    } else {
        pos = ll_part_positions("Set", llx, idx, len, &n);
        if (!pos) return false;
        nsel = (size_t)n;
        if (nrest == 0 && ll_assign_packed(ll, pos, n, rhs)) { free(pos); return true; }
    }

    Expr** vals = malloc(sizeof(Expr*) * (size_t)(n ? n : 1));
    for (int64_t j = 0; j < n; j++) {
        Expr* piece = ll_rhs_piece(rhs, (size_t)j, nsel);
        if (nrest == 0) {
            vals[j] = piece;
        } else {
            Expr* cur = ll_get(ll, pos[j]);
            vals[j] = expr_part_assign_rec(cur, rest, nrest, piece);
            expr_free(cur);
            expr_free(piece);
        }
    }
    for (int64_t j = 0; j < n; j++) {
        if (ll_would_cycle(ll, vals[j])) {
            ll_cycle_message("Set", vals[j]);
            for (int64_t t = 0; t < n; t++) expr_free(vals[t]);
            free(vals);
            if (pos != &one) free(pos);
            return false;
        }
    }
    for (int64_t j = 0; j < n; j++) ll_set(ll, pos[j], vals[j]);   /* consumes */
    free(vals);
    if (pos != &one) free(pos);
    return true;
}

Expr* expr_part(Expr* expr, Expr** indices, size_t nindices) {
    if (!expr) return NULL;
    if (nindices == 0) return expr_copy(expr);   /* the empty path: expr itself */
    if (!indices) return NULL;

    Expr* idx = indices[0];
    Expr** rest = indices + 1;
    size_t nrest = nindices - 1;

    /* Association indexing (Wolfram semantics, see assoc_part): a single
     * integer / key / Key[k] selects one VALUE; All, a Span or a List of
     * positions or keys selects a sub-ASSOCIATION. */
    if (is_association(expr)) return assoc_part(expr, idx, rest, nrest);

    /* A payload kind reached below the top level ({c0, sa}[[2, 1]], {ll}[[1, 3]],
     * {pv}[[1, 2]]): the remaining spec is the kind's own Part, through its
     * ds_ops row.  A SparseArray answers [[0]] itself (List); the other kinds
     * decline it so head extraction below handles it.  An NDArray keeps its own
     * path (it is atomic here). */
    if (expr->type != EXPR_NDARRAY) {
        int st = DS_PART_DECLINE;
        Expr* r = ds_part(expr, indices, nindices, &st);
        if (st == DS_PART_DONE) return r;
    }

    // Handle integer index
    if (idx->type == EXPR_INTEGER) {
        int64_t k = idx->data.integer;
        
        if (k == 0) {
            // Head extraction - allowed even for atoms
            Expr* head = expr_head(expr);
            if (!head) return NULL;
            Expr* result = expr_part(head, rest, nrest);
            expr_free(head);
            return result;
        } else {
            if (is_atomic(expr)) return NULL;
            int64_t len = (int64_t)expr->data.function.arg_count;
            if (k < 0) k = len + k + 1;
            if (k < 1 || k > len) return NULL;
            return expr_part(expr->data.function.args[k - 1], rest, nrest);
        }
    }

    // Handle "Span"
    if (idx->type == EXPR_FUNCTION && idx->data.function.head->type == EXPR_SYMBOL &&
        idx->data.function.head->data.symbol.name == SYM_Span) {
        if (is_atomic(expr)) return NULL;
        int64_t len = (int64_t)expr->data.function.arg_count;
        int64_t start, step, count;
        if (!span_resolve(idx, len, &start, &step, &count)) return NULL;
        if (!span_check_take(idx, len, expr)) return NULL;

        Expr** args = NULL;
        if (count > 0) {
            args = malloc(sizeof(Expr*) * count);
            if (!args) return NULL;
        }
        
        int64_t current = start;
        for (int64_t i = 0; i < count; i++) {
            Expr* fake_idx = expr_new_integer(current);
            Expr** new_indices = malloc(sizeof(Expr*) * nindices);
            new_indices[0] = fake_idx;
            for (size_t j = 0; j < nrest; j++) new_indices[j + 1] = rest[j];
            
            args[i] = expr_part(expr, new_indices, nindices);
            
            free(new_indices);
            expr_free(fake_idx);
            
            if (!args[i]) {
                for (int64_t k = 0; k < i; k++) expr_free(args[k]);
                if (args) free(args);
                return NULL;
            }
            current += step;
        }
        
        Expr* result = expr_new_function(expr_copy(expr->data.function.head), args, count);
        if (args) free(args);
        return result;
    }

    // Handle "All"
    if (idx->type == EXPR_SYMBOL && idx->data.symbol.name == SYM_All) {
        if (is_atomic(expr)) return NULL;
        size_t len = expr->data.function.arg_count;
        Expr** args = NULL;
        if (len > 0) {
            args = malloc(sizeof(Expr*) * len);
            if (!args) return NULL;
        }
        
        for (size_t i = 0; i < len; i++) {
            Expr* fake_idx = expr_new_integer((int64_t)i + 1);
            Expr** new_indices = malloc(sizeof(Expr*) * nindices);
            new_indices[0] = fake_idx;
            for (size_t j = 0; j < nrest; j++) new_indices[j + 1] = rest[j];
            
            args[i] = expr_part(expr, new_indices, nindices);
            
            free(new_indices);
            expr_free(fake_idx);
            
            if (!args[i]) {
                for (size_t k = 0; k < i; k++) expr_free(args[k]);
                free(args);
                return NULL;
            }
        }
        Expr* result = expr_new_function(expr_copy(expr->data.function.head), args, len);
        free(args);
        return result;
    }

    // Handle List of indices
    if (idx->type == EXPR_FUNCTION && idx->data.function.head->data.symbol.name == SYM_List) {
        if (is_atomic(expr)) return NULL;
        size_t len = idx->data.function.arg_count;
        Expr** args = NULL;
        if (len > 0) {
            args = malloc(sizeof(Expr*) * len);
            if (!args) return NULL;
        }
        
        for (size_t i = 0; i < len; i++) {
            Expr** new_indices = malloc(sizeof(Expr*) * nindices);
            new_indices[0] = idx->data.function.args[i];
            for (size_t j = 0; j < nrest; j++) new_indices[j + 1] = rest[j];
            
            args[i] = expr_part(expr, new_indices, nindices);
            
            free(new_indices);
            
            if (!args[i]) {
                for (size_t k = 0; k < i; k++) expr_free(args[k]);
                free(args);
                return NULL;
            }
        }
        Expr* result = expr_new_function(expr_copy(expr->data.function.head), args, len);
        free(args);
        return result;
    }

    // Invalid index type
    return NULL;
}

/* How one element of an association Part spec addresses the entries. */
typedef enum { APS_POS, APS_KEY, APS_BAD } AssocSpecKind;

/* Classify an association Part spec element. A machine integer is a position;
 * a string or Key[k] is a key (*key is set to the BORROWED key expression).
 * Anything else -- a bare symbol, a real, ... -- is not a valid part spec in
 * Mathematica (Part::pkspec1), so the whole Part stays unevaluated. */
static AssocSpecKind assoc_spec_kind(const Expr* s, const Expr** key) {
    if (s->type == EXPR_INTEGER) return APS_POS;
    if (s->type == EXPR_STRING) { *key = s; return APS_KEY; }
    if (s->type == EXPR_FUNCTION && s->data.function.head->type == EXPR_SYMBOL &&
        s->data.function.head->data.symbol.name == SYM_Key &&
        s->data.function.arg_count == 1) {
        *key = s->data.function.args[0];
        return APS_KEY;
    }
    return APS_BAD;
}

/* Missing["KeyAbsent", spec] -- spec is the index as written (Key[k] or "k"),
 * exactly as Mathematica reports it. */
static Expr* assoc_key_absent(const Expr* spec) {
    Expr* margs[2] = { expr_new_string("KeyAbsent"), expr_copy((Expr*)spec) };
    return expr_new_function(expr_new_symbol(SYM_Missing), margs, 2);
}

/* The entry `e` with its value replaced by value[[rest]], or NULL when that
 * deeper Part fails. With no deeper indices the entry itself is shared (a
 * refcount bump), which keeps a RuleDelayed entry delayed. */
static Expr* assoc_entry_part(const Expr* e, Expr** rest, size_t nrest) {
    if (nrest == 0) return expr_copy((Expr*)e);
    Expr* v = expr_part(e->data.function.args[1], rest, nrest);
    if (!v) return NULL;
    return assoc_entry_with_value(e, v);
}


/* Build the sub-association from `n` owned entries (NULL entries mean a
 * deeper Part failed: everything is freed and NULL returned). Duplicate keys
 * collapse as in the Association constructor. */
static Expr* assoc_part_collect(Expr** ents, size_t n, bool ok) {
    Expr* out = ok ? assoc_from_rules(ents, n) : NULL;
    for (size_t i = 0; i < n; i++) if (ents[i]) expr_free(ents[i]);
    free(ents);
    return out;
}

/*
 * Part of an association, one index level (`rest` are the deeper indices,
 * applied to the selected values). Mirrors Mathematica 15:
 *
 *   a[[i]] / a[["k"]] / a[[Key[k]]]   the value (then [[rest]]); an absent key
 *                                     gives Missing["KeyAbsent", spec] with
 *                                     the spec as written, e.g. Key[k]
 *   a[[0]]                            the head, Association
 *   a[[All]], a[[m;;n;;s]]            a sub-association of those entries
 *   a[[{i, j}]], a[[{"k", Key[k]}]]   a sub-association in the given order; an
 *                                     absent key maps to its Missing[...]
 *
 * NULL (the Part stays unevaluated, as after Part::partw / pkspec1 / pmix /
 * take / partd) for a position out of range, a Span running off the end, a
 * spec list mixing positions with keys, a non-key atom such as a bare symbol,
 * or a deeper index that fails on any selected value.
 */
static Expr* assoc_part(Expr* assoc, Expr* idx, Expr** rest, size_t nrest) {
    size_t na = assoc->data.function.arg_count;

    const Expr* key = NULL;
    AssocSpecKind kind = assoc_spec_kind(idx, &key);

    /* Single-key and single-position reads touch one entry, so they must not
     * pay an O(n) well-formedness scan (that made a[[Key[k]]] O(n) per read):
     * assoc_lookup_value is index-backed and scan-safe on a malformed node,
     * and the positional read checks just the entry it dereferences. */
    if (kind == APS_KEY) {
        Expr* v = assoc_lookup_value(assoc, key);   /* O(1) via the key index */
        return v ? expr_part(v, rest, nrest) : assoc_key_absent(idx);
    }
    if (kind == APS_POS && idx->data.integer != 0) {
        int64_t pos = idx->data.integer;
        if (pos < 0) pos = (int64_t)na + pos + 1;
        if (pos < 1 || pos > (int64_t)na) return NULL;
        if (!is_rule2(assoc->data.function.args[pos - 1])) return NULL;   /* malformed */
        return expr_part(assoc->data.function.args[pos - 1]->data.function.args[1], rest, nrest);
    }

    for (size_t i = 0; i < na; i++)
        if (!is_rule2(assoc->data.function.args[i])) return NULL;  /* malformed */

    if (kind == APS_POS) {
        int64_t pos = idx->data.integer;
        if (pos == 0) {
            Expr* head = expr_head(assoc);
            if (!head) return NULL;
            Expr* result = expr_part(head, rest, nrest);
            expr_free(head);
            return result;
        }
        if (pos < 0) pos = (int64_t)na + pos + 1;
        if (pos < 1 || pos > (int64_t)na) return NULL;
        return expr_part(assoc->data.function.args[pos - 1]->data.function.args[1], rest, nrest);
    }
    if (kind == APS_KEY) {
        Expr* v = assoc_lookup_value(assoc, key);   /* O(1) via the key index */
        return v ? expr_part(v, rest, nrest) : assoc_key_absent(idx);
    }

    bool is_all  = idx->type == EXPR_SYMBOL && idx->data.symbol.name == SYM_All;
    bool is_span = idx->type == EXPR_FUNCTION && idx->data.function.head->type == EXPR_SYMBOL &&
                   idx->data.function.head->data.symbol.name == SYM_Span;
    bool is_list = idx->type == EXPR_FUNCTION && idx->data.function.head->type == EXPR_SYMBOL &&
                   idx->data.function.head->data.symbol.name == SYM_List;

    if (is_all || is_span) {
        int64_t start = 1, step = 1, count = (int64_t)na;
        if (is_span && (!span_check_take(idx, (int64_t)na, assoc) ||
                        !span_resolve(idx, (int64_t)na, &start, &step, &count)))
            return NULL;
        Expr** ents = malloc(sizeof(Expr*) * (count > 0 ? (size_t)count : 1));
        bool ok = true;
        int64_t k = start;
        for (int64_t i = 0; i < count; i++, k += step) {
            ents[i] = ok ? assoc_entry_part(assoc->data.function.args[k - 1], rest, nrest) : NULL;
            if (!ents[i]) ok = false;
        }
        return assoc_part_collect(ents, (size_t)count, ok);
    }

    if (!is_list) return NULL;                       /* pkspec1 */

    /* A list of positions or of keys -- never a mixture (Part::pmix). */
    size_t m = idx->data.function.arg_count;
    AssocSpecKind lkind = APS_BAD;
    for (size_t j = 0; j < m; j++) {
        const Expr* kj = NULL;
        AssocSpecKind kd = assoc_spec_kind(idx->data.function.args[j], &kj);
        if (kd == APS_BAD || (j > 0 && kd != lkind)) return NULL;
        lkind = kd;
        if (kd == APS_POS) {
            int64_t p = idx->data.function.args[j]->data.integer;
            if (p < 0) p = (int64_t)na + p + 1;
            if (p < 1 || p > (int64_t)na) return NULL;   /* partw */
        }
    }
    Expr** ents = malloc(sizeof(Expr*) * (m ? m : 1));
    bool ok = true;
    for (size_t j = 0; j < m; j++) {
        Expr* sj = idx->data.function.args[j];
        ents[j] = NULL;
        if (!ok) continue;
        if (lkind == APS_POS) {
            int64_t p = sj->data.integer;
            if (p < 0) p = (int64_t)na + p + 1;
            ents[j] = assoc_entry_part(assoc->data.function.args[p - 1], rest, nrest);
        } else {
            const Expr* kj = NULL;
            (void)assoc_spec_kind(sj, &kj);
            Expr* v = assoc_lookup_value(assoc, kj);
            Expr* nv = v ? expr_part(v, rest, nrest) : assoc_key_absent(sj);
            if (nv) {
                Expr* rargs[2] = { expr_copy((Expr*)kj), nv };
                ents[j] = expr_new_function(expr_new_symbol(SYM_Rule), rargs, 2);
            }
        }
        if (!ents[j]) ok = false;
    }
    return assoc_part_collect(ents, m, ok);
}

Expr* builtin_part(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) {
        return NULL;
    }

    Expr* expr = res->data.function.args[0];
    Expr** indices = res->data.function.args + 1;
    size_t nindices = res->data.function.arg_count - 1;

    /* A LinkedList of indices names parts as a List of them does:
     * expr[[..., {i1, i2, ...}, ...]] is the parts expr[[..., ik, ...]], under
     * expr's own head (a LinkedList for a LinkedList). Read from the slots. */
    for (size_t j = 0; j < nindices; j++) {
        if (indices[j]->type != EXPR_LINKEDLIST) continue;
        LinkedList* ix = ll_copy(indices[j]->data.llist);
        int64_t m = ll_length(ix);
        bool to_ll = expr->type == EXPR_LINKEDLIST;
        Expr** parts = malloc(sizeof(Expr*) * (size_t)(m ? m : 1));
        int64_t got = 0;
        for (; got < m; got++) {
            Expr** a = malloc(sizeof(Expr*) * (nindices + 1));
            a[0] = expr_copy(expr);
            for (size_t t = 0; t < nindices; t++)
                a[t + 1] = (t == j) ? ll_get(ix, got) : expr_copy(indices[t]);
            Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Part), a, nindices + 1));
            free(a);
            if (head_is(r, SYM_Part) || eval_is_inflight_unwind(r)) {
                /* an index out of range: the message has been given; the call
                 * stays as it stands */
                bool unwind = eval_is_inflight_unwind(r);
                for (int64_t k = 0; k < got; k++) expr_free(parts[k]);
                free(parts);
                ll_release(ix);
                if (unwind) return r;
                expr_free(r);
                return NULL;
            }
            parts[got] = r;
        }
        ll_release(ix);
        Expr* out;
        if (to_ll) {
            LinkedList* o = ll_new();
            for (int64_t k = 0; k < m; k++) ll_push_back(o, parts[k]);
            out = expr_new_linkedlist(o);
        } else if (expr->type == EXPR_FUNCTION) {
            out = pack_offer(expr_new_function(expr_copy(expr->data.function.head), parts, (size_t)m));
        } else {
            out = pack_offer(expr_new_function(expr_new_symbol(SYM_List), parts, (size_t)m));
        }
        free(parts);
        return out;
    }

    /* PACK THE SMALL SOURCE UP. A gather's cost is set by the number of
     * INDICES, but the packing decision was made about the SOURCE in isolation:
     * a 65-entry lookup table is correctly left unpacked, and a 262144-entry
     * packed index then has to be met on the boxed path -- one Expr per output
     * element, out of a table that is a single strided read: 52.7 ms for
     * 262144 indices into a 65-entry table, against 3.1 ms once it is lifted. Same rule as Dot, Join,
     * Outer and the Listable gate: pack the small operand up, never materialise
     * the large one down. Looking up sphere parameters by hit index is what a
     * ray tracer ends with, and the table is always small.
     *
     * A packed index list is already at or above the packing threshold, so this
     * only ever fires on a genuinely large gather. Anything ndarray_part
     * declines falls straight through to the ordinary path below with the
     * ORIGINAL argument, so lifting cannot change an answer -- including an
     * out-of-range diagnostic. */
    if (!is_ndarray(expr) && nindices == 1 && pack_enabled() &&
        is_packed_list(indices[0]) &&
        indices[0]->data.ndarray.rank == 1 &&
        indices[0]->data.ndarray.dtype == NDT_INT64 &&
        expr->type == EXPR_FUNCTION &&
        expr->data.function.head->type == EXPR_SYMBOL &&
        expr->data.function.head->data.symbol.name == SYM_List) {
        Expr* lifted = pack_force(expr_copy(expr), false, NDT_FLOAT64);
        if (is_ndarray(lifted)) {
            bool ldeg = false;
            Expr* r = ndarray_part(lifted, indices, nindices, &ldeg);
            expr_free(lifted);
            if (r) return r;
        } else {
            expr_free(lifted);
        }
    }

    /* Packed array / SparseArray / LinkedList: element access through the
     * ds_ops row. A NDArray [[0]] (head extraction) declines to the general
     * path; Span/All/List positions on an NDArray degrade to the general List
     * Part via delist-and-reeval; a part that does not exist leaves the Part
     * unevaluated. */
    {
        int st = DS_PART_DECLINE;
        Expr* r = ds_part(expr, indices, nindices, &st);
        if (st == DS_PART_DEGRADE) return ndarray_delist_and_reeval(res);
        if (st == DS_PART_DONE) return r;
    }

    // [[0]] is allowed even for atoms
    if (nindices > 0 && indices[0]->type == EXPR_INTEGER && indices[0]->data.integer == 0) {
        // Allowed
    } else {
        if (is_atomic(expr)) return NULL;
    }

    return expr_part(expr, indices, nindices);
}

/* Part of `expr` at one position path, honouring a BUFFER.
 *
 * expr_part cannot index an NDArray: is_atomic() deliberately reports true for
 * one (see the note beside expr_head), so it returns NULL, builtin_extract
 * returns NULL with it, and the call comes to rest -- where the POST-gate
 * materialises the buffer and the ordinary List path answers. Correct, and
 * 99.7 ms to read three elements out of 10^6 (HPC plan item C.4), because the
 * cost is the million Expr nodes nobody wanted.
 *
 * ndarray_part is the gather Part has used since experiment 12. Extract is the
 * same operation on the same value and was simply never wired to it; [[0]] is
 * excluded here exactly as it is in builtin_part, because head extraction is
 * expr_head's job and an NDArray's head is not one of its elements.
 *
 * NULL preserves the old behaviour exactly rather than being an error: the
 * caller declines, the post-gate materialises, the ordinary path answers. */
static Expr* extract_part(Expr* expr, Expr** indices, size_t nindices) {
    if (nindices > 0 && is_ndarray(expr) &&
        !(indices[0]->type == EXPR_INTEGER && indices[0]->data.integer == 0)) {
        bool degrade = false;
        return ndarray_part(expr, indices, nindices, &degrade);
    }
    return expr_part(expr, indices, nindices);
}

static Expr* extract_single(Expr* expr, Expr* pos, Expr* h) {
    size_t nindices;
    Expr** indices;
    Expr* single_index_buf[1];
    if (pos->type == EXPR_FUNCTION && pos->data.function.head->type == EXPR_SYMBOL &&
        pos->data.function.head->data.symbol.name == SYM_List) {
        nindices = pos->data.function.arg_count;
        indices = pos->data.function.args;
    } else {
        /* Scalar position (e.g. integer) is treated as a one-element path,
         * matching Mathematica: Extract[expr, n] == Extract[expr, {n}]. */
        single_index_buf[0] = pos;
        indices = single_index_buf;
        nindices = 1;
    }
    Expr* part = extract_part(expr, indices, nindices);
    if (!part) {
        return NULL;
    }
    if (h) {
        Expr* args[1] = { part };
        return expr_new_function(expr_copy(h), args, 1);
    }
    return part;
}

Expr* builtin_extract(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 1) { // Operator form
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner_extract = expr_new_function(expr_new_symbol(SYM_Extract), inner_args, 2);
        Expr* func_args[1] = { inner_extract };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }
    if (argc < 2 || argc > 3) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* pos = res->data.function.args[1];
    Expr* h = (argc == 3) ? res->data.function.args[2] : NULL;

    /* A LinkedList position: one path when its elements are indices
     * (Extract[e, LinkedList[{i, j}]] is e[[i, j]]), else a list of positions
     * (a List of the extracts, as for a List of positions). */
    if (pos->type == EXPR_LINKEDLIST) {
        LinkedList* pl = ll_copy(pos->data.llist);
        int64_t m = ll_length(pl);
        bool paths = false;
        for (int64_t i = 0; i < m && !paths; i++) {
            bool own;
            Expr* el = ll_elem_view(pl, i, &own);
            paths = ll_seq_length(el) >= 0;
            ll_elem_done(el, own);
        }
        Expr* out;
        if (!paths) {
            Expr** a = malloc(sizeof(Expr*) * (size_t)(m + 1));
            a[0] = expr_copy(expr);
            for (int64_t i = 0; i < m; i++) a[i + 1] = ll_get(pl, i);
            out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Part), a, (size_t)m + 1));
            free(a);
            if (head_is(out, SYM_Part)) {          /* no such part: as for a List */
                expr_free(out);
                ll_release(pl);
                return NULL;
            }
            if (h) {
                Expr* ha[1] = { out };
                out = expr_new_function(expr_copy(h), ha, 1);
            }
        } else {
            Expr** items = malloc(sizeof(Expr*) * (size_t)(m ? m : 1));
            for (int64_t i = 0; i < m; i++) {
                Expr* a[3] = { expr_copy(expr), ll_get(pl, i), h ? expr_copy(h) : NULL };
                items[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Extract), a, h ? 3 : 2));
            }
            out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)m);
            free(items);
        }
        ll_release(pl);
        return out;
    }

    bool is_list_of_pos = false;
    if (pos->type == EXPR_FUNCTION && pos->data.function.head->type == EXPR_SYMBOL &&
        pos->data.function.head->data.symbol.name == SYM_List) {
        if (pos->data.function.arg_count > 0 && pos->data.function.args[0]->type == EXPR_FUNCTION &&
            pos->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
            pos->data.function.args[0]->data.function.head->data.symbol.name == SYM_List) {
            is_list_of_pos = true;
        }
    }
    /* Otherwise pos is treated as a single position (a scalar like an integer,
     * or a single-position list); extract_single handles both forms. The
     * exception is {}: an empty LIST of positions, so the answer is {} ({{}}
     * is the one position naming the whole expression). */
    if (pos->type == EXPR_FUNCTION && pos->data.function.head->type == EXPR_SYMBOL &&
        pos->data.function.head->data.symbol.name == SYM_List &&
        pos->data.function.arg_count == 0)
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    if (is_list_of_pos) {
        size_t npos = pos->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * npos);
        for (size_t i = 0; i < npos; i++) {
            args[i] = extract_single(expr, pos->data.function.args[i], h);
            if (!args[i]) {
                for (size_t j = 0; j < i; j++) expr_free(args[j]);
                free(args);
                return NULL;
            }
        }
        Expr* r = expr_new_function(expr_new_symbol(SYM_List), args, npos);
        free(args);   /* expr_new_function copies the array; free our copy */
        return r;
    } else {
        return extract_single(expr, pos, h);
    }
}
/* For an association entry (a Rule[k,v]), the "element" First/Last yields is the
 * value v, not the whole rule -- matching Wolfram. Returns the entry unchanged
 * for non-associations. */
static Expr* first_last_element(Expr* container, Expr* entry) {
    if (is_association(container) && entry->type == EXPR_FUNCTION &&
        entry->data.function.arg_count == 2) {
        return entry->data.function.args[1];
    }
    return entry;
}

/* First / Last of a LinkedList: O(1). An empty list gives the default when
 * there is one, else First::nofirst / Last::nolast and stays unevaluated. */
Expr* ll_first_last_expr(Expr* arg, Expr* deflt, bool last);
static Expr* ll_first_last(Expr* arg, Expr* deflt, bool last) {
    LinkedList* ll = arg->data.llist;
    int64_t n = ll_length(ll);
    if (n > 0) return ll_get(ll, last ? n - 1 : 0);
    if (deflt) return expr_copy(deflt);
    char* s = expr_to_string(arg);
    if (last) mth_message("Last", "nolast", "%s has zero length and no last element.", s ? s : "");
    else mth_message("First", "nofirst", "%s has zero length and no first element.", s ? s : "");
    free(s);
    return NULL;
}

Expr* ll_first_last_expr(Expr* arg, Expr* deflt, bool last) { return ll_first_last(arg, deflt, last); }

/* First / Last of an atom or of an expression with no elements and no
 * default: ::normal or ::nofirst / ::nolast, and the call stays unevaluated
 * (NULL), as in the reference (Last[x], Last[{}]). */
static Expr* first_last_empty(const Expr* call, Expr* arg, bool last) {
    const char* head = last ? "Last" : "First";
    if (expr_is_plain_atom(arg)) return builtin_normal_error(head, call, 1);
    /* a payload kind (SparseArray, ...) is left to its gate, silently */
    if (arg->type != EXPR_FUNCTION || arg->data.function.arg_count > 0) return NULL;
    char* s = expr_to_string(arg);
    if (last) mth_message("Last", "nolast", "%s has zero length and no last element.", s ? s : "");
    else mth_message("First", "nofirst", "%s has zero length and no first element.", s ? s : "");
    free(s);
    return NULL;
}

Expr* builtin_first(Expr* res) {
    if (res->type != EXPR_FUNCTION ||
        res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr* deflt = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    /* An NDArray is atomic, so is_atomic below rejects it and the call comes
     * back UNEVALUATED, while the identical List call works.  First is a
     * leading-axis slice of the buffer, so take it directly; delist_repack (which
     * materialises every element) remains the fallback for what it declines. */
    { bool dh; Expr* dr = ds_native(arg, SYM_First, res, &dh); if (dh) return dr; }
    if (is_ndarray(arg)) {
        Expr* fast = ndstruct_head_tail(res, ND_FIRST);
        if (fast) return fast;
        return ndstruct_delist_repack(res, arg);
    }
    /* First[expr, def] returns def when expr has no elements (or is atomic);
     * with no default the call is First::normal (an atom) or First::nofirst
     * (no elements) and stays unevaluated. */
    if (is_atomic(arg) || arg->data.function.arg_count < 1)
        return deflt ? expr_copy(deflt) : first_last_empty(res, arg, false);
    return expr_copy(first_last_element(arg, arg->data.function.args[0]));
}

Expr* builtin_last(Expr* res) {
    if (res->type != EXPR_FUNCTION ||
        res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr* deflt = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    /* An NDArray is atomic, so is_atomic below rejects it and the call comes
     * back UNEVALUATED, while the identical List call works.  Last is a
     * leading-axis slice of the buffer, so take it directly; delist_repack (which
     * materialises every element) remains the fallback for what it declines. */
    { bool dh; Expr* dr = ds_native(arg, SYM_Last, res, &dh); if (dh) return dr; }
    if (is_ndarray(arg)) {
        Expr* fast = ndstruct_head_tail(res, ND_LAST);
        if (fast) return fast;
        return ndstruct_delist_repack(res, arg);
    }
    if (is_atomic(arg) || arg->data.function.arg_count < 1)
        return deflt ? expr_copy(deflt) : first_last_empty(res, arg, true);
    return expr_copy(first_last_element(arg, arg->data.function.args[arg->data.function.arg_count - 1]));
}

Expr* builtin_most(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    { bool dh; Expr* dr = ds_native(arg, SYM_Most, res, &dh); if (dh) return dr; }
    /* An NDArray is atomic, so is_atomic below rejects it and the call comes
     * back UNEVALUATED, while the identical List call works.  Most is a
     * leading-axis slice of the buffer, so take it directly; delist_repack (which
     * materialises every element) remains the fallback for what it declines. */
    if (is_ndarray(arg)) {
        Expr* fast = ndstruct_head_tail(res, ND_MOST);
        if (fast) return fast;
        return ndstruct_delist_repack(res, arg);
    }
    if (is_atomic(arg) || arg->data.function.arg_count < 1) return NULL;
    
    size_t new_count = arg->data.function.arg_count - 1;
    Expr** new_args = NULL;
    if (new_count > 0) {
        new_args = malloc(sizeof(Expr*) * new_count);
        for (size_t i = 0; i < new_count; i++) {
            new_args[i] = expr_copy(arg->data.function.args[i]);
        }
    }
    Expr* result = expr_new_function(expr_copy(arg->data.function.head), new_args, new_count);
    free(new_args);
    return result;
}

Expr* builtin_rest(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    { bool dh; Expr* dr = ds_native(arg, SYM_Rest, res, &dh); if (dh) return dr; }
    /* An NDArray is atomic, so is_atomic below rejects it and the call comes
     * back UNEVALUATED, while the identical List call works.  Rest is a
     * leading-axis slice of the buffer, so take it directly; delist_repack (which
     * materialises every element) remains the fallback for what it declines. */
    if (is_ndarray(arg)) {
        Expr* fast = ndstruct_head_tail(res, ND_REST);
        if (fast) return fast;
        return ndstruct_delist_repack(res, arg);
    }
    if (is_atomic(arg) || arg->data.function.arg_count < 1) return NULL;
    
    size_t new_count = arg->data.function.arg_count - 1;
    Expr** new_args = NULL;
    if (new_count > 0) {
        new_args = malloc(sizeof(Expr*) * new_count);
        for (size_t i = 0; i < new_count; i++) {
            new_args[i] = expr_copy(arg->data.function.args[i + 1]);
        }
    }
    Expr* result = expr_new_function(expr_copy(arg->data.function.head), new_args, new_count);
    free(new_args);
    return result;
}

static Expr* insert_path(Expr* expr, Expr* elem, Expr** path, size_t path_len) {
    if (is_atomic(expr)) return expr_copy(expr);
    if (path_len == 0) return expr_copy(expr);

    int64_t n = 0;
    if (path[0]->type != EXPR_INTEGER) return expr_copy(expr);
    n = path[0]->data.integer;

    if (path_len == 1) {
        size_t len = expr->data.function.arg_count;
        int64_t idx = n;
        if (idx < 0) idx = (int64_t)len + idx + 2;
        if (idx < 1) idx = 1;
        if (idx > (int64_t)len + 1) idx = (int64_t)len + 1;
        
        size_t insert_pos = (size_t)idx - 1;
        size_t new_count = len + 1;
        Expr** new_args = malloc(sizeof(Expr*) * new_count);
        
        size_t j = 0;
        for (size_t i = 0; i < new_count; i++) {
            if (i == insert_pos) {
                new_args[i] = expr_copy(elem);
            } else {
                new_args[i] = expr_copy(expr->data.function.args[j++]);
            }
        }
        Expr* result = expr_new_function(expr_copy(expr->data.function.head), new_args, new_count);
        free(new_args);
        return result;
    } else {
        size_t len = expr->data.function.arg_count;
        int64_t target = n;
        if (target < 0) target = (int64_t)len + target + 1;
        if (target < 1 || target > (int64_t)len) return expr_copy(expr);
        
        size_t target_idx = (size_t)target - 1;
        Expr** new_args = malloc(sizeof(Expr*) * len);
        for (size_t i = 0; i < len; i++) {
            if (i == target_idx) {
                new_args[i] = insert_path(expr->data.function.args[i], elem, path + 1, path_len - 1);
            } else {
                new_args[i] = expr_copy(expr->data.function.args[i]);
            }
        }
        Expr* result = expr_new_function(expr_copy(expr->data.function.head), new_args, len);
        free(new_args);
        return result;
    }
}

/* Descending order on paths, for applying several Insert/Delete positions one
 * after another: a lexicographically later position is applied first, so it
 * never shifts the part an earlier one names. The paths must already be
 * normalised (normalize_path) -- a raw negative index sorts below every
 * positive one although it may name a later part. A non-integer step (Key[k]
 * in an association) orders canonically after the integers. */
static int compare_paths(const void* a, const void* b) {
    Expr* pa = *(Expr**)a;
    Expr* pb = *(Expr**)b;

    size_t la = (pa->type == EXPR_FUNCTION) ? pa->data.function.arg_count : 1;
    size_t lb = (pb->type == EXPR_FUNCTION) ? pb->data.function.arg_count : 1;

    size_t n = (la < lb) ? la : lb;
    for (size_t i = 0; i < n; i++) {
        Expr* ea = (pa->type == EXPR_FUNCTION) ? pa->data.function.args[i] : pa;
        Expr* eb = (pb->type == EXPR_FUNCTION) ? pb->data.function.args[i] : pb;
        bool ia = ea->type == EXPR_INTEGER, ib = eb->type == EXPR_INTEGER;
        if (ia && ib) {
            int64_t va = ea->data.integer, vb = eb->data.integer;
            if (va < vb) return 1; // Reverse sort
            if (va > vb) return -1;
            continue;
        }
        if (ia != ib) return ia ? 1 : -1;
        int c = expr_compare(ea, eb);
        if (c) return c < 0 ? 1 : -1;
    }
    if (la < lb) return 1;
    if (la > lb) return -1;
    return 0;
}

/* One position of a list of positions, with every negative index resolved
 * against the ORIGINAL expression: each position in Insert[expr, e, {p1, p2,
 * ...}] / Delete[expr, {p1, p2, ...}] names a part (or, for Insert's last
 * index, a slot) of expr itself, not of the expression the earlier positions
 * have already changed. Resolving first is what makes the descending
 * sequential application correct: Delete[{a, b, c}, {{-3}, {2}}] is {c}, and
 * the aliases {-1} and {3} of a length-3 list become the same path, which the
 * caller then removes once. Resolution stops (keeping the rest verbatim) at a
 * non-integer step, an association, an atom, or an index out of range -- the
 * walkers handle those as before. A new List. */
static Expr* normalize_path(Expr* expr, Expr* path, bool insert) {
    size_t n = path->type == EXPR_FUNCTION ? path->data.function.arg_count : 1;
    Expr** steps = path->type == EXPR_FUNCTION ? path->data.function.args : &path;
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    Expr* e = expr;
    for (size_t d = 0; d < n; d++) {
        Expr* p = steps[d];
        if (!e || p->type != EXPR_INTEGER || is_atomic(e) || is_association(e)) {
            out[d] = expr_copy(p);
            e = NULL;
            continue;
        }
        int64_t len = (int64_t)e->data.function.arg_count;
        int64_t k = p->data.integer;
        bool last = d + 1 == n;
        if (k < 0) {
            int64_t r = (insert && last) ? len + 2 + k : len + 1 + k;
            if (r >= 1) k = r;
        }
        out[d] = (k == p->data.integer) ? expr_copy(p) : expr_new_integer(k);
        if (k == 0) e = e->data.function.head;
        else if (k >= 1 && k <= len) e = e->data.function.args[k - 1];
        else e = NULL;
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), out, n);
    free(out);
    return r;
}

Expr* expr_insert(Expr* expr, Expr* elem, Expr* pos) {
    if (!expr || !elem || !pos) return NULL;
    
    // Case 3: List of positions
    if (pos->type == EXPR_FUNCTION && pos->data.function.head->data.symbol.name == SYM_List) {
        bool all_lists = true;
        for (size_t i = 0; i < pos->data.function.arg_count; i++) {
            if (pos->data.function.args[i]->type != EXPR_FUNCTION) {
                all_lists = false;
                break;
            }
        }
        
        if (all_lists && pos->data.function.arg_count > 0) {
            // Sort positions in descending order to apply them safely
            size_t np = pos->data.function.arg_count;
            Expr** sorted_pos = malloc(sizeof(Expr*) * np);
            for (size_t i = 0; i < np; i++)
                sorted_pos[i] = normalize_path(expr, pos->data.function.args[i], true);
            qsort(sorted_pos, np, sizeof(Expr*), compare_paths);

            /* A repeated slot inserts once per listing (the reference's
             * Insert[{a, b}, x, {{1}, {1}}] is {x, x, a, b}), so no dedup. */
            Expr* current = expr_copy(expr);
            for (size_t i = 0; i < np; i++) {
                Expr* next = insert_path(current, elem, sorted_pos[i]->data.function.args,
                                         sorted_pos[i]->data.function.arg_count);
                expr_free(current);
                current = next;
            }
            for (size_t i = 0; i < np; i++) expr_free(sorted_pos[i]);
            free(sorted_pos);
            return current;
        }
        
        // Case 2: Single path
        return insert_path(expr, elem, pos->data.function.args, pos->data.function.arg_count);
    }
    
    // Case 1: Single integer
    if (pos->type == EXPR_INTEGER) {
        return insert_path(expr, elem, &pos, 1);
    }
    
    return expr_copy(expr);
}

static int64_t assoc_entry_index(const Expr* assoc, const Expr* idx);

/* The insertion slot (0..n: before entry slot, n = at the end) named by one
 * Insert position on an association of n entries -- an integer (negative
 * counts from the end, -1 being the end), a string key or Key[k] (insert
 * before that entry) -- or -1 when out of range or the key is absent. */
static int64_t assoc_insert_slot(const Expr* assoc, const Expr* p) {
    int64_t n = (int64_t)assoc->data.function.arg_count;
    if (p->type == EXPR_INTEGER) {
        int64_t v = p->data.integer;
        if (v > 0) return v <= n + 1 ? v - 1 : -1;
        if (v < 0) { v = n + v + 2; return v >= 1 ? v - 1 : -1; }
        return -1;
    }
    if (p->type == EXPR_STRING ||
        (p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL &&
         p->data.function.head->data.symbol.name == SYM_Key && p->data.function.arg_count == 1))
        return assoc_entry_index(assoc, p);
    return -1;
}

static int slot_cmp(const void* a, const void* b) {
    int64_t x = *(const int64_t*)a, y = *(const int64_t*)b;
    return (x > y) - (x < y);
}

/*
 * Insert[assoc, rule(s), pos] (Mathematica 15). The rule, or list of rules, is
 * placed at pos -- an integer, a key, Key[k], {pos} or {{pos1}, {pos2}, ...}
 * -- and an existing entry with the same key is removed, so the NEW position
 * wins: Insert[<|a -> 1, b -> 2, c -> 3|>, a -> 9, 3] is <|b -> 2, a -> 9, c -> 3|>.
 * A deeper path {p, i, ...} inserts into the value at p. Anything but rules
 * leaves the association unchanged (Mathematica says Insert::invdt and returns
 * it); an out-of-range position or absent key leaves the call unevaluated.
 */
static Expr* assoc_insert(Expr* assoc, Expr* elem, Expr* pos) {
    size_t n = assoc->data.function.arg_count;
    for (size_t i = 0; i < n; i++)
        if (!is_rule2(assoc->data.function.args[i])) return NULL;
    bool is_list_pos = pos->type == EXPR_FUNCTION && pos->data.function.head->type == EXPR_SYMBOL &&
                       pos->data.function.head->data.symbol.name == SYM_List;
    bool list_of_paths = is_list_pos && pos->data.function.arg_count > 0 &&
                         pos->data.function.args[0]->type == EXPR_FUNCTION &&
                         pos->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
                         pos->data.function.args[0]->data.function.head->data.symbol.name == SYM_List;

    /* A deep path {p, rest...}: insert elem into the value at p. */
    if (is_list_pos && !list_of_paths && pos->data.function.arg_count >= 2) {
        int64_t e = assoc_entry_index(assoc, pos->data.function.args[0]);
        if (e < 0) return NULL;
        Expr* entry = assoc->data.function.args[e];
        Expr* val = entry->data.function.args[1];
        if (is_atomic(val)) return NULL;
        Expr* nv = insert_path(val, elem, pos->data.function.args + 1, pos->data.function.arg_count - 1);
        Expr** ents = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++)
            ents[i] = (i == (size_t)e) ? assoc_entry_with_value(entry, nv)
                                       : expr_copy(assoc->data.function.args[i]);
        Expr* out = expr_new_function(expr_copy(assoc->data.function.head), ents, n);
        free(ents);
        return out;
    }

    /* The rules to insert. */
    Expr** nr; size_t nn;
    if (is_rule2(elem)) { nr = &elem; nn = 1; }
    else if (elem->type == EXPR_FUNCTION && elem->data.function.head->type == EXPR_SYMBOL &&
             elem->data.function.head->data.symbol.name == SYM_List) {
        nr = elem->data.function.args; nn = elem->data.function.arg_count;
        for (size_t i = 0; i < nn; i++) if (!is_rule2(nr[i])) return expr_copy(assoc);
    } else {
        return expr_copy(assoc);                        /* Insert::invdt */
    }

    /* The insertion slots. */
    size_t ns = 1;
    int64_t one;
    int64_t* slots = &one;
    if (is_list_pos) {
        size_t m = pos->data.function.arg_count;
        if (m == 0) return NULL;
        if (!list_of_paths) {
            one = assoc_insert_slot(assoc, pos->data.function.args[0]);   /* {p} */
        } else {
            slots = malloc(sizeof(int64_t) * m);
            ns = m;
            for (size_t i = 0; i < m; i++) {
                Expr* pi = pos->data.function.args[i];
                bool ok = pi->type == EXPR_FUNCTION && pi->data.function.arg_count == 1 &&
                          pi->data.function.head->type == EXPR_SYMBOL &&
                          pi->data.function.head->data.symbol.name == SYM_List;
                slots[i] = ok ? assoc_insert_slot(assoc, pi->data.function.args[0]) : -1;
            }
            qsort(slots, ns, sizeof(int64_t), slot_cmp);
        }
    } else {
        one = assoc_insert_slot(assoc, pos);
    }
    for (size_t i = 0; i < ns; i++)
        if (slots[i] < 0) { if (slots != &one) free(slots); return NULL; }

    /* Lay out old entries (minus those whose key is being inserted) and the new
     * rules at each slot; assoc_from_rules then collapses repeated new keys. */
    Expr* newset = assoc_from_rules(nr, nn);            /* O(1) key membership */
    Expr** out = malloc(sizeof(Expr*) * (n + nn * ns + 1));
    size_t k = 0, s = 0;
    for (size_t i = 0; i <= n; i++) {
        while (s < ns && slots[s] == (int64_t)i) {
            for (size_t j = 0; j < nn; j++) out[k++] = nr[j];
            s++;
        }
        if (i < n) {
            Expr* r = assoc->data.function.args[i];
            if (!assoc_lookup_value(newset, r->data.function.args[0])) out[k++] = r;
        }
    }
    Expr* result = assoc_from_rules(out, k);            /* copies the (borrowed) rules */
    free(out);
    expr_free(newset);
    if (slots != &one) free(slots);
    return result;
}

static Expr* number_atom_delete_insert(Expr* res, Expr* atom, Expr* spec, bool insert);

/* ---- Out-of-range positions (Insert::ins, Delete::partw) -----------------
 *
 * insert_path and delete_path clamp or ignore a position they cannot reach:
 * Insert[{1, 2}, x, 10] gave {1, 2, x} and Delete[{1, 2}, 5] gave {1, 2}. The
 * reference raises a message and leaves the call unevaluated, which is also
 * what Check and Quiet need. So every position is checked BEFORE anything is
 * built; a path through an association, or with a non-integer index (Key[k]),
 * is left to the walkers, which handle those. */

enum { PC_OK, PC_BAD, PC_BAD_ATOM };

/* Can `p` (n indices) be reached in `expr`? For Insert the last index is a
 * slot, 1 .. len + 1 or -1 .. -(len + 1); for Delete an element, or 0 (the
 * head, which Delete replaces by Sequence). PC_BAD_ATOM: the path runs into
 * an atom at depth > 0; *atom / *idx are that atom and the index tried on it. */
static int pos_check(Expr* expr, Expr** p, size_t n, bool insert, Expr** atom, Expr** idx) {
    Expr* e = expr;
    for (size_t d = 0; d < n; d++) {
        if (p[d]->type != EXPR_INTEGER || is_association(e)) return PC_OK;
        if (is_atomic(e)) {
            if (d == 0) return PC_BAD;
            *atom = e; *idx = p[d];
            return PC_BAD_ATOM;
        }
        int64_t len = (int64_t)e->data.function.arg_count;
        int64_t k = p[d]->data.integer;
        if (d + 1 == n) {
            if (insert) return (k != 0 && (k > 0 ? k <= len + 1 : -k <= len + 1)) ? PC_OK : PC_BAD;
            return (k == 0 || (k > 0 ? k <= len : -k <= len)) ? PC_OK : PC_BAD;
        }
        if (k == 0) { e = e->data.function.head; continue; }
        int64_t j = k > 0 ? k : len + k + 1;
        if (j < 1 || j > len) return PC_BAD;
        e = e->data.function.args[j - 1];
    }
    return PC_OK;
}

static Expr* pos_as_list(Expr* path) {
    if (path->type == EXPR_FUNCTION && head_is(path, SYM_List)) return expr_copy(path);
    Expr* one = expr_copy(path);
    return expr_new_function(expr_new_symbol(SYM_List), &one, 1);
}

static void pos_message(bool insert, Expr* expr, Expr* path, int why, Expr* atom, Expr* idx) {
    char* es = expr_to_string(why == PC_BAD_ATOM ? atom : expr);
    Expr* shown = why == PC_BAD_ATOM ? expr_copy(idx) : pos_as_list(path);
    char* ps = expr_to_string(shown);
    if (insert)
        mth_message("Insert", "ins", "Cannot insert at position %s in %s", ps ? ps : "", es ? es : "");
    else
        mth_message("Delete", "partw", "Part %s of %s does not exist.", ps ? ps : "", es ? es : "");
    free(ps);
    free(es);
    expr_free(shown);
}

/* Check every path of `pos` against `expr`, reading the spec exactly as
 * expr_insert / expr_delete do. False after the message. A list of positions
 * names its first bad one for Insert and its last for Delete, which processes
 * the list from the end -- both as the reference reports them. */
static bool positions_reachable(Expr* expr, Expr* pos, bool insert) {
    Expr** paths;
    size_t np;
    bool listed = pos->type == EXPR_FUNCTION && pos->data.function.head->type == EXPR_SYMBOL
               && pos->data.function.head->data.symbol.name == SYM_List;
    if (pos->type == EXPR_INTEGER) { paths = &pos; np = 1; }
    else if (!listed) return true;
    else {
        size_t m = pos->data.function.arg_count;
        if (m == 0) return true;                    /* {}: no position */
        bool several = true;
        for (size_t i = 0; i < m && several; i++) {
            Expr* a = pos->data.function.args[i];
            several = insert ? a->type == EXPR_FUNCTION
                             : (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL
                                && a->data.function.head->data.symbol.name == SYM_List);
        }
        if (several) { paths = pos->data.function.args; np = m; }
        else { paths = &pos; np = 1; }
    }
    for (size_t t = 0; t < np; t++) {
        Expr* path = paths[insert ? t : np - 1 - t];
        Expr** ix;
        size_t n;
        if (path->type == EXPR_INTEGER) { ix = &path; n = 1; }
        else if (path->type == EXPR_FUNCTION) { ix = path->data.function.args; n = path->data.function.arg_count; }
        else continue;
        Expr* atom = NULL;
        Expr* idx = NULL;
        int why = pos_check(expr, ix, n, insert, &atom, &idx);
        /* Insert names the whole position even when it runs into an atom;
         * Delete names the index and the atom ("Part 1 of 1"). */
        if (why == PC_BAD_ATOM && insert) why = PC_BAD;
        if (why != PC_OK) { pos_message(insert, expr, path, why, atom, idx); return false; }
    }
    return true;
}

/* One index of an Insert position path: a machine integer, or a key (a string
 * or Key[k]) for an association level. */
static bool insert_index_ok(const Expr* e) {
    return e->type == EXPR_INTEGER || e->type == EXPR_STRING || head_is((Expr*)e, SYM_Key);
}

/* An Insert position spec: an index, a path {i, j, ...} of indices, or a list
 * of such paths -- anything else (a symbol, 1.5, 10^30, a mix {{1}, 2}) is
 * Insert::psl and the call stays unevaluated (audit Insert#3; it used to come
 * back as the input unchanged). */
static bool insert_pos_spec_ok(const Expr* pos) {
    if (insert_index_ok(pos)) return true;
    if (!head_is((Expr*)pos, SYM_List)) return false;
    size_t m = pos->data.function.arg_count;
    if (m == 0) return true;
    bool paths = head_is(pos->data.function.args[0], SYM_List);
    for (size_t i = 0; i < m; i++) {
        const Expr* a = pos->data.function.args[i];
        if (paths) {
            if (!head_is((Expr*)a, SYM_List)) return false;
            for (size_t j = 0; j < a->data.function.arg_count; j++)
                if (!insert_index_ok(a->data.function.args[j])) return false;
        } else if (!insert_index_ok(a)) {
            return false;
        }
    }
    return true;
}

Expr* builtin_insert(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    if (!insert_pos_spec_ok(res->data.function.args[2])) {
        char* ps = expr_to_string(res->data.function.args[2]);
        char* cs = expr_to_string(res);
        mth_message("Insert", "psl", "Position specification %s in %s is not a machine-sized "
                    "integer or a list of machine-sized integers.", ps ? ps : "", cs ? cs : "");
        free(ps); free(cs);
        return NULL;
    }
    { bool dh; Expr* dr = ds_native(res->data.function.args[0], SYM_Insert, res, &dh); if (dh) return dr; }
    if (expr_is_number_atom(res->data.function.args[0]))
        return number_atom_delete_insert(res, res->data.function.args[0], res->data.function.args[2], true);
    if (is_association(res->data.function.args[0]))
        return assoc_insert(res->data.function.args[0], res->data.function.args[1],
                            res->data.function.args[2]);
    /* An atom at the top keeps the walker's answer (it is not a position error). */
    if (!is_atomic(res->data.function.args[0]) &&
        !positions_reachable(res->data.function.args[0], res->data.function.args[2], true))
        return NULL;
    return expr_insert(res->data.function.args[0], res->data.function.args[1], res->data.function.args[2]);
}

static Expr* delete_path(Expr* expr, Expr** path, size_t path_len) {
    if (path_len == 0) {
        return expr_new_function(expr_new_symbol(SYM_Sequence), NULL, 0);
    }

    if (is_atomic(expr)) return expr_copy(expr);

    /* Association: a non-integer index is a key (Key[k] or a literal key). Find
     * the entry, then either drop it (last index) or recurse into its value. */
    if (is_association(expr) && path[0]->type != EXPR_INTEGER) {
        Expr* key = path[0];
        if (path[0]->type == EXPR_FUNCTION && path[0]->data.function.head->type == EXPR_SYMBOL &&
            path[0]->data.function.head->data.symbol.name == SYM_Key && path[0]->data.function.arg_count == 1)
            key = path[0]->data.function.args[0];
        size_t len = expr->data.function.arg_count;
        int64_t t = -1;
        for (size_t i = 0; i < len; i++) {
            Expr* rule = expr->data.function.args[i];
            if (rule->type == EXPR_FUNCTION && rule->data.function.arg_count == 2 &&
                expr_eq(rule->data.function.args[0], key)) { t = (int64_t)i; break; }
        }
        if (t < 0) return expr_copy(expr);
        if (path_len == 1) {
            Expr** na = (len > 1) ? malloc(sizeof(Expr*) * (len - 1)) : NULL;
            size_t j = 0;
            for (size_t i = 0; i < len; i++)
                if ((int64_t)i != t) na[j++] = expr_copy(expr->data.function.args[i]);
            Expr* result = expr_new_function(expr_copy(expr->data.function.head), na, len - 1);
            if (na) free(na);
            return result;
        }
        Expr** na = malloc(sizeof(Expr*) * len);
        for (size_t i = 0; i < len; i++) na[i] = expr_copy(expr->data.function.args[i]);
        /* na[t] is a refcount-bumped share of the source entry, so writing its
         * args[1] in place would corrupt every other holder (the exact aliasing
         * bug assoc_entry_with_value exists to prevent — cf. path_descend_child).
         * Rebuild the entry around the recursed value instead. */
        Expr* rule = na[t];
        Expr* nv = delete_path(rule->data.function.args[1], path + 1, path_len - 1);
        na[t] = assoc_entry_with_value(rule, nv);   /* fresh entry: copies key+head, adopts nv */
        expr_free(rule);                            /* drop our ref on the shared entry */
        Expr* result = expr_new_function(expr_copy(expr->data.function.head), na, len);
        free(na);
        return result;
    }

    int64_t n = 0;
    if (path[0]->type != EXPR_INTEGER) return expr_copy(expr);
    n = path[0]->data.integer;

    if (path_len == 1) {
        if (n == 0) {
            Expr** new_args = NULL;
            if (expr->data.function.arg_count > 0) {
                new_args = malloc(sizeof(Expr*) * expr->data.function.arg_count);
                for (size_t i = 0; i < expr->data.function.arg_count; i++) {
                    new_args[i] = expr_copy(expr->data.function.args[i]);
                }
            }
            Expr* result = expr_new_function(expr_new_symbol(SYM_Sequence), new_args, expr->data.function.arg_count);
            if (new_args) free(new_args);
            return result;
        }

        size_t len = expr->data.function.arg_count;
        int64_t idx = n;
        if (idx < 0) idx = (int64_t)len + idx + 1;
        if (idx < 1 || idx > (int64_t)len) return expr_copy(expr);

        size_t target_idx = (size_t)idx - 1;
        size_t new_count = len - 1;
        Expr** new_args = NULL;
        if (new_count > 0) {
            new_args = malloc(sizeof(Expr*) * new_count);
            size_t j = 0;
            for (size_t i = 0; i < len; i++) {
                if (i != target_idx) {
                    new_args[j++] = expr_copy(expr->data.function.args[i]);
                }
            }
        }
        Expr* result = expr_new_function(expr_copy(expr->data.function.head), new_args, new_count);
        if (new_args) free(new_args);
        return result;
    } else {
        if (n == 0) {
            Expr* new_head = delete_path(expr->data.function.head, path + 1, path_len - 1);
            Expr** new_args = NULL;
            if (expr->data.function.arg_count > 0) {
                new_args = malloc(sizeof(Expr*) * expr->data.function.arg_count);
                for (size_t i = 0; i < expr->data.function.arg_count; i++) {
                    new_args[i] = expr_copy(expr->data.function.args[i]);
                }
            }
            Expr* result = expr_new_function(new_head, new_args, expr->data.function.arg_count);
            if (new_args) free(new_args);
            return result;
        }

        size_t len = expr->data.function.arg_count;
        int64_t target = n;
        if (target < 0) target = (int64_t)len + target + 1;
        if (target < 1 || target > (int64_t)len) return expr_copy(expr);
        
        size_t target_idx = (size_t)target - 1;
        Expr** new_args = NULL;
        if (len > 0) {
            new_args = malloc(sizeof(Expr*) * len);
            for (size_t i = 0; i < len; i++) {
                if (i == target_idx) {
                    new_args[i] = delete_path(expr->data.function.args[i], path + 1, path_len - 1);
                } else {
                    new_args[i] = expr_copy(expr->data.function.args[i]);
                }
            }
        }
        Expr* result = expr_new_function(expr_copy(expr->data.function.head), new_args, len);
        if (new_args) free(new_args);
        return result;
    }
}

Expr* expr_delete(Expr* expr, Expr* pos) {
    if (!expr || !pos) return NULL;

    // Case 3: List of positions
    if (pos->type == EXPR_FUNCTION && pos->data.function.head->data.symbol.name == SYM_List) {
        /* Multiple-paths iff every element is itself a List (a sub-path).
         * A non-List step such as Key[k] means `pos` is a single path, e.g.
         * {Key["a"], Key["x"]} is one nested path, not two. */
        bool all_lists = pos->data.function.arg_count > 0;
        for (size_t i = 0; i < pos->data.function.arg_count; i++) {
            Expr* e = pos->data.function.args[i];
            if (!(e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
                  e->data.function.head->data.symbol.name == SYM_List)) {
                all_lists = false;
                break;
            }
        }
        
        if (all_lists && pos->data.function.arg_count > 0) {
            size_t np = pos->data.function.arg_count;
            Expr** sorted_pos = malloc(sizeof(Expr*) * np);
            for (size_t i = 0; i < np; i++)
                sorted_pos[i] = normalize_path(expr, pos->data.function.args[i], false);
            qsort(sorted_pos, np, sizeof(Expr*), compare_paths);

            Expr* current = expr_copy(expr);
            for (size_t i = 0; i < np; i++) {
                /* A repeated position names the same original part once:
                 * Delete[Range[12], {{3}, {3}}] removes 3, not 3 and then the
                 * 4 that slid into its place. Normalised, an alias ({-1} and
                 * {3} of a length-3 list) is a repeat too. */
                if (i > 0 && expr_eq(sorted_pos[i - 1], sorted_pos[i])) continue;
                Expr* next = delete_path(current, sorted_pos[i]->data.function.args,
                                         sorted_pos[i]->data.function.arg_count);
                expr_free(current);
                current = next;
            }
            for (size_t i = 0; i < np; i++) expr_free(sorted_pos[i]);
            free(sorted_pos);
            return current;
        }
        
        /* {} is an empty list of positions, not the path to the whole
         * expression ({{}} is that): nothing is deleted.
         * Delete[l, Position[l, 0]] is the idiom it serves. */
        if (pos->data.function.arg_count == 0) return expr_copy(expr);

        // Case 2: Single path
        return delete_path(expr, pos->data.function.args, pos->data.function.arg_count);
    }
    
    // Case 1: Single integer
    if (pos->type == EXPR_INTEGER) {
        return delete_path(expr, &pos, 1);
    }
    
    return expr_copy(expr);
}

/* The integer the reference names in Delete::partw / Insert::ins for a
 * position spec: n, {n, ...} -> n, {{...}, ..., {m, ...}} -> the first index of
 * the LAST position (the reference processes a position list from the end). */
static Expr* atom_spec_first_index(Expr* spec) {
    if (spec->type == EXPR_FUNCTION && head_is(spec, SYM_List) && spec->data.function.arg_count > 0) {
        Expr* last = spec->data.function.args[spec->data.function.arg_count - 1];
        if (last->type == EXPR_FUNCTION && head_is(last, SYM_List)) {
            return last->data.function.arg_count > 0 ? last->data.function.args[0] : NULL;
        }
        return spec->data.function.args[0];
    }
    return spec;
}

/* Delete/Insert on a Rational/Complex number (an atom, no parts): {} positions
 * leave it unchanged; otherwise the reference's message and the call stays
 * unevaluated. Returns the result, or NULL after the message. */
static Expr* number_atom_delete_insert(Expr* res, Expr* atom, Expr* spec, bool insert) {
    if (spec->type == EXPR_FUNCTION && head_is(spec, SYM_List) && spec->data.function.arg_count == 0)
        return expr_copy(atom);
    Expr* idx = atom_spec_first_index(spec);
    char* is = idx ? expr_to_string(idx) : NULL;
    char* as = expr_to_string(atom);
    if (insert)
        mth_message("Insert", "ins", "Cannot insert at position {%s} in %s", is ? is : "", as ? as : "");
    else
        mth_message("Delete", "partw", "Part %s of %s does not exist.", is ? is : "", as ? as : "");
    free(is);
    free(as);
    (void)res;
    return NULL;
}

Expr* builtin_delete(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count == 2 && expr_is_number_atom(res->data.function.args[0]))
        return number_atom_delete_insert(res, res->data.function.args[0], res->data.function.args[1], false);
    if (res->data.function.arg_count == 1) {
        return NULL; // Operator form, stays unevaluated
    } else if (res->data.function.arg_count == 2) {
        { bool dh; Expr* dr = ds_native(res->data.function.args[0], SYM_Delete, res, &dh); if (dh) return dr; }
        /* Delete[assoc, Key[k]] and Delete[assoc, "k"] name one entry by its
         * key, as the one-element path {Key[k]} does (audit Delete#3 /
         * Association#1 / Key#2). */
        {
            Expr* pos = res->data.function.args[1];
            if (is_association(res->data.function.args[0]) &&
                (head_is(pos, SYM_Key) || pos->type == EXPR_STRING)) {
                Expr* path_args[1] = { expr_copy(pos) };
                Expr* path = expr_new_function(expr_new_symbol(SYM_List), path_args, 1);
                Expr* out = NULL;
                if (positions_reachable(res->data.function.args[0], path, false))
                    out = expr_delete(res->data.function.args[0], path);
                expr_free(path);
                return out;
            }
        }
        if (!positions_reachable(res->data.function.args[0], res->data.function.args[1], false))
            return NULL;
        return expr_delete(res->data.function.args[0], res->data.function.args[1]);
    }
    return NULL;
}

/* ---------------- Position-path walker (MapAt / ReplaceAt) ---------------- */

/*
 * Resolve one association index spec to a 0-based entry index: Key[k] and a
 * bare literal key look the key up, an integer is positional over the entries
 * (negatives counting from the end).  Returns -1 when the key is absent or the
 * position is out of range.  All, Span and the head index 0 are handled by the
 * caller, before this triage, because they address the entry *list* rather
 * than naming a single key.
 */
static int64_t assoc_entry_index(const Expr* assoc, const Expr* idx) {
    size_t n = assoc->data.function.arg_count;
    const Expr* key;

    if (idx->type == EXPR_FUNCTION && idx->data.function.head->type == EXPR_SYMBOL &&
        idx->data.function.head->data.symbol.name == SYM_Key &&
        idx->data.function.arg_count == 1) {
        key = idx->data.function.args[0];
    } else if (idx->type == EXPR_INTEGER) {
        int64_t p = idx->data.integer;
        if (p < 0) p = (int64_t)n + p + 1;
        if (p < 1 || p > (int64_t)n) return -1;
        return p - 1;
    } else {
        key = idx;   /* a string or any other literal key */
    }

    for (size_t i = 0; i < n; i++) {
        Expr* rule = assoc->data.function.args[i];
        if (part_is_entry(rule) && expr_eq(rule->data.function.args[0], (Expr*)key))
            return (int64_t)i;
    }
    return -1;
}

/*
 * Walk the remaining path into child `i` of the rebuilt argument array and
 * store the result back.  For an association the child is an entry, and it is
 * the entry's *value* that is walked; a FRESH Rule node is built rather than
 * assigning into the existing one, because expr_copy is a refcount bump --
 * new_args[i] is the very node the caller passed in, and writing through it
 * would mutate the caller's expression.  Returns false if the walk failed.
 */
static bool path_descend_child(Expr** new_args, size_t i, bool assoc,
                               Expr** path, size_t plen, PartLeafFn fn, void* ctx) {
    Expr* child = new_args[i];

    if (!assoc) {
        Expr* r = expr_apply_at_path(child, path, plen, fn, ctx);
        if (!r) return false;
        expr_free(new_args[i]);
        new_args[i] = r;
        return true;
    }

    if (!part_is_entry(child)) return false;
    Expr* nv = expr_apply_at_path(child->data.function.args[1], path, plen, fn, ctx);
    if (!nv) return false;

    new_args[i] = assoc_entry_with_value(child, nv);   /* adopts nv */
    expr_free(child);
    return true;
}

Expr* expr_apply_at_path(Expr* expr, Expr** path, size_t plen,
                         PartLeafFn fn, void* ctx) {
    if (!expr || !fn) return NULL;
    if (plen == 0) return fn(ctx, expr);
    if (!path) return NULL;
    /* A LinkedList: a new LinkedList with the named elements walked (an
     * integer, negative from the end, All or a Span; position 0 -- the head --
     * is not a part of a LinkedList). A walked element is evaluated as it is
     * stored, since the evaluator never looks inside a LinkedList. */
    if (expr->type == EXPR_LINKEDLIST) {
        Expr* idx = path[0];
        int64_t len = ll_length(expr->data.llist), start = 1, step = 1, count = 0;
        if (idx->type == EXPR_INTEGER && idx->data.integer != 0) {
            int64_t k = idx->data.integer;
            if (k < 0) k = len + k + 1;
            if (k < 1 || k > len) return NULL;
            start = k; count = 1;
        } else if (idx->type == EXPR_SYMBOL && idx->data.symbol.name == SYM_All) {
            count = len;
        } else if (idx->type == EXPR_FUNCTION && idx->data.function.head->type == EXPR_SYMBOL &&
                   idx->data.function.head->data.symbol.name == SYM_Span) {
            if (!span_resolve(idx, len, &start, &step, &count)) return NULL;
        } else {
            return NULL;
        }
        LinkedList* out = ll_copy(expr->data.llist);
        int64_t cur = start;
        for (int64_t c = 0; c < count; c++, cur += step) {
            Expr* el = ll_get(out, cur - 1);
            Expr* r = expr_apply_at_path(el, path + 1, plen - 1, fn, ctx);
            expr_free(el);
            if (!r) { ll_release(out); return NULL; }
            r = eval_and_free(r);
            if (eval_is_inflight_unwind(r)) { ll_release(out); return r; }
            ll_set(out, cur - 1, r);
        }
        return expr_new_linkedlist(out);
    }
    /* Cannot descend into an atom. is_atomic (not a bare type check) so that
     * Rational and Complex -- stored as EXPR_FUNCTION nodes -- stay leaves, as
     * they do for Part and MapIndexed: MapAt[f, 1/2, 1] must not manufacture
     * Rational[f[1], 2]. */
    if (is_atomic(expr)) return NULL;

    Expr*  idx   = path[0];
    Expr** rest  = path + 1;
    size_t nrest = plen - 1;
    size_t len   = expr->data.function.arg_count;
    bool   assoc = is_association(expr);

    Expr** new_args = NULL;
    if (len > 0) {
        new_args = malloc(sizeof(Expr*) * len);
        if (!new_args) return NULL;
        for (size_t i = 0; i < len; i++)
            new_args[i] = expr_copy(expr->data.function.args[i]);
    }
    Expr* new_head = expr_copy(expr->data.function.head);
    bool ok = true;

    if (idx->type == EXPR_INTEGER && idx->data.integer == 0) {
        /* Position 0 is the head, for associations too. */
        Expr* r = expr_apply_at_path(expr->data.function.head, rest, nrest, fn, ctx);
        if (r) { expr_free(new_head); new_head = r; }
        else ok = false;
    } else if (idx->type == EXPR_SYMBOL && idx->data.symbol.name == SYM_All) {
        for (size_t i = 0; i < len && ok; i++)
            ok = path_descend_child(new_args, i, assoc, rest, nrest, fn, ctx);
    } else if (idx->type == EXPR_FUNCTION && idx->data.function.head->type == EXPR_SYMBOL &&
               idx->data.function.head->data.symbol.name == SYM_Span) {
        int64_t start = 1, step = 1, count = 0;
        if (!span_resolve(idx, (int64_t)len, &start, &step, &count)) {
            ok = false;
        } else {
            int64_t cur = start;
            for (int64_t c = 0; c < count && ok; c++, cur += step)
                ok = path_descend_child(new_args, (size_t)(cur - 1), assoc,
                                        rest, nrest, fn, ctx);
        }
    } else if (assoc) {
        int64_t target = assoc_entry_index(expr, idx);
        if (target < 0) ok = false;
        else ok = path_descend_child(new_args, (size_t)target, true, rest, nrest, fn, ctx);
    } else if (idx->type == EXPR_INTEGER) {
        int64_t k = idx->data.integer;
        if (k < 0) k = (int64_t)len + k + 1;
        if (k < 1 || k > (int64_t)len) ok = false;
        else ok = path_descend_child(new_args, (size_t)(k - 1), false, rest, nrest, fn, ctx);
    } else {
        ok = false;   /* Key[...] on a non-association, a real, a string, ... */
    }

    if (!ok) {
        for (size_t i = 0; i < len; i++) expr_free(new_args[i]);
        free(new_args);
        expr_free(new_head);
        return NULL;
    }

    Expr* result = expr_new_function(new_head, new_args, len);
    free(new_args);
    return result;
}

Expr* expr_apply_at_positions(Expr* expr, Expr* pos, PartLeafFn fn, void* ctx) {
    if (!expr || !pos || !fn) return NULL;

    bool pos_is_list = head_is(pos, SYM_List);

    /* {} is an empty list of *positions*, so nothing is mapped.  ({{}} is one
     * position -- the empty path -- and addresses the whole expression.) */
    if (pos_is_list && pos->data.function.arg_count == 0) return expr_copy(expr);

    /* A List whose first element is itself a List is a list of paths; anything
     * else is a single path, possibly spelled as a List of indices. */
    if (pos_is_list && head_is(pos->data.function.args[0], SYM_List)) {
        Expr* current = expr_copy(expr);
        for (size_t i = 0; i < pos->data.function.arg_count; i++) {
            Expr* sub = pos->data.function.args[i];
            Expr* next = head_is(sub, SYM_List)
                ? expr_apply_at_path(current, sub->data.function.args,
                                     sub->data.function.arg_count, fn, ctx)
                : expr_apply_at_path(current, &sub, 1, fn, ctx);
            expr_free(current);
            if (!next) return NULL;
            current = next;
        }
        return current;
    }

    if (pos_is_list)
        return expr_apply_at_path(expr, pos->data.function.args,
                                  pos->data.function.arg_count, fn, ctx);
    return expr_apply_at_path(expr, &pos, 1, fn, ctx);
}
