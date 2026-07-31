// Part and Head

#include "eval.h"
#include "ndstruct.h"   /* ndstruct_delist_repack — packed-argument fallback */
#include "symtab.h"
#include "sym_names.h"
#include "assoc.h"
#include "ndarray.h"
#include "part.h"
#include "common.h"

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
static Expr* expr_part_assign_rec(Expr* expr, Expr** indices, size_t nindices, Expr* rhs, size_t* rhs_idx, bool is_rhs_list);
/* Single-index association Part: resolve one key/Key[k]/positional index into
 * its value (recursing for any remaining indices), or Missing["KeyAbsent", k].
 * Returns NULL only for an out-of-range positional index. */
static Expr* assoc_part_single(Expr* assoc, Expr* idx, Expr** rest, size_t nrest);

/* Assign through the value of association entry `entry`, recursing for any
 * remaining indices, and return the resulting entry. Takes ownership of `entry`
 * (the caller's array slot) and hands back a fresh one -- never an in-place
 * write, because expr_copy is a refcount bump and `entry` is the very node the
 * caller's association holds. Assigning into args[1] here is what made
 * `u2 = u1; u2[["x"]] = 9` corrupt `u1`. Returns `entry` unchanged when it is
 * not a two-argument rule (Association[1, 2] can exist). */
static Expr* assoc_entry_assigned(Expr* entry, Expr** rest, size_t nrest,
                                  Expr* rhs, size_t* rhs_idx, bool is_rhs_list) {
    if (!part_is_entry(entry)) return entry;
    Expr* nv = expr_part_assign_rec(entry->data.function.args[1], rest, nrest,
                                    rhs, rhs_idx, is_rhs_list);
    Expr* fresh = assoc_entry_with_value(entry, nv);   /* adopts nv */
    expr_free(entry);
    return fresh;
}


/* Assign into a dense NDArray. Value semantics, as everywhere else in Part
 * assignment: a fresh NDArray comes back and `arr` is untouched.
 *
 * The whole Part spec vocabulary works — integer (including negative), All,
 * Span, a list of positions, and fewer subscripts than the rank — because the
 * write reuses ndarray_part_set, which selects positions through the same
 * per-axis selector ndarray_part reads through. Before that shared selector
 * existed this handled only a full run of integer subscripts and returned the
 * array UNCHANGED for anything else, so `a[[2 ;; 4]] = 0` on an NDArray silently
 * did nothing while the same code on a List worked.
 *
 * Returns the array unchanged (the old behaviour) for a spec or right-hand side
 * outside what the buffer can represent — a symbolic value, a complex value into
 * a real buffer, an out-of-range subscript. */
static Expr* nd_part_assign(Expr* arr, Expr** indices, size_t nindices, Expr* rhs) {
    /* An INDEPENDENT node with its own buffer.  Not expr_copy: that is a
     * refcount bump, so writing through it would reach every other holder of
     * the same array — including the symbol this assignment is meant to leave
     * untouched until the caller stores the result back. */
    size_t n = 1;
    for (int i = 0; i < arr->data.ndarray.rank; i++) n *= (size_t)arr->data.ndarray.dims[i];
    size_t bytes = ndt_elem_size(arr->data.ndarray.dtype) * n;
    void* nb = malloc(bytes ? bytes : 1);
    if (!nb) return expr_copy(arr);
    memcpy(nb, arr->data.ndarray.data, bytes);
    Expr* fresh = expr_new_ndarray_like(arr, arr->data.ndarray.rank, arr->data.ndarray.dims,
                                   nb, arr->data.ndarray.dtype);   /* adopts nb */
    if (!fresh) { free(nb); return expr_copy(arr); }
    if (ndarray_part_set(fresh, (Expr* const*)indices, nindices, rhs)) return fresh;
    expr_free(fresh);
    return expr_copy(arr);
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
 * here for the same reason.
 *
 * A List right-hand side distributes its elements across the selected
 * positions, so every leaf has to qualify. */
static bool nd_rhs_fits_dtype(const Expr* rhs, NDType dt) {
    if (!rhs) return false;
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

static Expr* expr_part_assign_rec(Expr* expr, Expr** indices, size_t nindices, Expr* rhs, size_t* rhs_idx, bool is_rhs_list) {
    if (nindices == 0) {
        if (is_rhs_list) {
            if (*rhs_idx < rhs->data.function.arg_count) {
                return expr_copy(rhs->data.function.args[(*rhs_idx)++]);
            } else {
                return expr_copy(rhs);
            }
        }
        return expr_copy(rhs);
    }

    /* Dense NDArray element assignment (a[[i]] = v, a[[i,j]] = v).
     *
     * A PACKED LIST unpacks first when the value being written cannot be stored
     * with its head intact. The buffer path coerces -- which is right for a
     * visible NDArray[...], where the user asked for a machine buffer and the
     * dtype is part of the value -- but on a packed list that would be a silent
     * answer change on ordinary list code: `p[[1]] = 1` must give {1, 2., 3.},
     * not {1., 2., 3.}, and `p[[1]] = 1/3` must give the exact Rational.
     *
     * The result is left unpacked. That is not laziness: every case that lands
     * here produces a list the buffer could not have held anyway (mixed heads,
     * a Rational, a symbol), so there is nothing to repack. */
    if (expr->type == EXPR_NDARRAY) {
        if (is_packed_list(expr) &&
            !nd_rhs_fits_dtype(rhs, expr->data.ndarray.dtype)) {
            Expr* plain = ndarray_to_nested_list(expr);
            Expr* out = expr_part_assign_rec(plain, indices, nindices,
                                             rhs, rhs_idx, is_rhs_list);
            expr_free(plain);
            return out;
        }
        return nd_part_assign(expr, indices, nindices, rhs);
    }

    if (is_atomic(expr)) return expr_copy(expr);

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
                new_args[i] = assoc_entry_assigned(new_args[i], rest, nrest, rhs, rhs_idx, is_rhs_list);
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
            if (step > 0)
                for (int64_t i = start; i <= end && i >= 1 && i <= (int64_t)len; i += step)
                    new_args[i - 1] = assoc_entry_assigned(new_args[i - 1], rest, nrest, rhs, rhs_idx, is_rhs_list);
        } else if (idx->type == EXPR_FUNCTION && idx->data.function.head->type == EXPR_SYMBOL &&
                   idx->data.function.head->data.symbol.name == SYM_List) {
            for (size_t j = 0; j < idx->data.function.arg_count; j++) {
                Expr* sub = idx->data.function.args[j];
                if (sub->type == EXPR_INTEGER) {
                    int64_t p = sub->data.integer; if (p < 0) p = (int64_t)len + p + 1;
                    if (p >= 1 && p <= (int64_t)len)
                        new_args[p - 1] = assoc_entry_assigned(new_args[p - 1], rest, nrest, rhs, rhs_idx, is_rhs_list);
                } else {
                    Expr* k = (sub->type == EXPR_FUNCTION && sub->data.function.head->type == EXPR_SYMBOL &&
                               sub->data.function.head->data.symbol.name == SYM_Key && sub->data.function.arg_count == 1)
                              ? sub->data.function.args[0] : sub;
                    for (size_t i = 0; i < len; i++)
                        if (new_args[i]->type == EXPR_FUNCTION && new_args[i]->data.function.arg_count == 2 &&
                            expr_eq(new_args[i]->data.function.args[0], k)) {
                            new_args[i] = assoc_entry_assigned(new_args[i], rest, nrest, rhs, rhs_idx, is_rhs_list); break;
                        }
                }
            }
        } else if (idx->type == EXPR_INTEGER) {
            int64_t pos = idx->data.integer;
            if (pos < 0) pos = (int64_t)len + pos + 1;
            if (pos >= 1 && pos <= (int64_t)len)
                new_args[pos - 1] = assoc_entry_assigned(new_args[pos - 1], rest, nrest, rhs, rhs_idx, is_rhs_list);
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
                new_args[found] = assoc_entry_assigned(new_args[found], rest, nrest, rhs, rhs_idx, is_rhs_list);
            } else if (nrest == 0) {
                /* nrest == 0 -> the recursive call returns the RHS value. */
                Expr* nv = expr_part_assign_rec(new_head, rest, nrest, rhs, rhs_idx, is_rhs_list);
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
            Expr* replaced = expr_part_assign_rec(new_head, rest, nrest, rhs, rhs_idx, is_rhs_list);
            expr_free(new_head);
            new_head = replaced;
        } else {
            if (k < 0) k = len + k + 1;
            if (k >= 1 && k <= (int64_t)len) {
                Expr* replaced = expr_part_assign_rec(new_args[k - 1], rest, nrest, rhs, rhs_idx, is_rhs_list);
                expr_free(new_args[k - 1]);
                new_args[k - 1] = replaced;
            }
        }
    } else if (idx->type == EXPR_SYMBOL && idx->data.symbol.name == SYM_All) {
        for (size_t i = 0; i < len; i++) {
            Expr* replaced = expr_part_assign_rec(new_args[i], rest, nrest, rhs, rhs_idx, is_rhs_list);
            expr_free(new_args[i]);
            new_args[i] = replaced;
        }
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

        if (step > 0) {
            for (int64_t i = start; i <= end && i >= 1 && i <= (int64_t)len; i += step) {
                Expr* replaced = expr_part_assign_rec(new_args[i - 1], rest, nrest, rhs, rhs_idx, is_rhs_list);
                expr_free(new_args[i - 1]);
                new_args[i - 1] = replaced;
            }
        } else if (step < 0) {
            for (int64_t i = start; i >= end && i >= 1 && i <= (int64_t)len; i += step) {
                Expr* replaced = expr_part_assign_rec(new_args[i - 1], rest, nrest, rhs, rhs_idx, is_rhs_list);
                expr_free(new_args[i - 1]);
                new_args[i - 1] = replaced;
            }
        }
    } else if (idx->type == EXPR_FUNCTION && idx->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < idx->data.function.arg_count; i++) {
            Expr* sub_idx = idx->data.function.args[i];
            if (sub_idx->type == EXPR_INTEGER) {
                int64_t k = sub_idx->data.integer;
                if (k == 0) {
                    Expr* replaced = expr_part_assign_rec(new_head, rest, nrest, rhs, rhs_idx, is_rhs_list);
                    expr_free(new_head);
                    new_head = replaced;
                } else {
                    if (k < 0) k = len + k + 1;
                    if (k >= 1 && k <= (int64_t)len) {
                        Expr* replaced = expr_part_assign_rec(new_args[k - 1], rest, nrest, rhs, rhs_idx, is_rhs_list);
                        expr_free(new_args[k - 1]);
                        new_args[k - 1] = replaced;
                    }
                }
            }
        }
    }

    Expr* new_func = expr_new_function(new_head, new_args, len);
    if (new_args) free(new_args);
    return new_func;
}

Expr* expr_part_assign(Expr* lhs, Expr* rhs) {
    if (lhs->type != EXPR_FUNCTION || lhs->data.function.head->type != EXPR_SYMBOL || lhs->data.function.head->data.symbol.name != SYM_Part) return NULL;
    if (lhs->data.function.arg_count < 2) return NULL;
    
    Expr* sym = lhs->data.function.args[0];
    if (sym->type != EXPR_SYMBOL) return NULL;
    
    Expr* current_val = symtab_get_own_values(sym->data.symbol.name) ? evaluate(sym) : NULL;
    if (!current_val) return NULL;
    
    bool is_rhs_list = (rhs->type == EXPR_FUNCTION && rhs->data.function.head->type == EXPR_SYMBOL && rhs->data.function.head->data.symbol.name == SYM_List);
    
    size_t rhs_idx = 0;
    Expr* new_val = expr_part_assign_rec(current_val, lhs->data.function.args + 1, lhs->data.function.arg_count - 1, rhs, &rhs_idx, is_rhs_list);
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
    switch (e->type) {
        case EXPR_INTEGER: return expr_new_symbol(SYM_Integer);
        case EXPR_BIGINT: return expr_new_symbol(SYM_Integer);
        case EXPR_REAL: return expr_new_symbol(SYM_Real);
#ifdef USE_MPFR
        case EXPR_MPFR: return expr_new_symbol(SYM_Real);
#endif
        /* A packed list IS a List and must say so. Note is_atomic() above
         * deliberately still reports true for it: the 13 sites that consult it
         * all go on to read data.function.arg_count, which aliases the ndarray
         * payload's `data` pointer, so flipping that predicate would turn a
         * heap pointer into an element count. Presentation is decided here and
         * in the AtomQ builtin, not in the internal walker guard. */
        case EXPR_NDARRAY:
            return expr_new_symbol(is_packed_list(e) ? SYM_List : SYM_NDArray);
        case EXPR_SYMBOL: return expr_new_symbol(SYM_Symbol);
        case EXPR_STRING: return expr_new_symbol(SYM_String);
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

Expr* expr_part(Expr* expr, Expr** indices, size_t nindices) {
    if (!expr || !indices) return NULL;
    if (nindices == 0) return expr_copy(expr);

    Expr* idx = indices[0];
    Expr** rest = indices + 1;
    size_t nrest = nindices - 1;

    /* Association indexing: assoc[[Key[k]]] / assoc[["strkey"]] look a key up
     * by value; assoc[[i]] with a positive/negative integer is positional
     * over the values (Wolfram semantics); assoc[[{k1,...}]] maps over the
     * sub-indices, giving {assoc[[k1]], ...}. A missing key yields
     * Missing["KeyAbsent", key]. */
    if (is_association(expr)) {
        if (idx->type == EXPR_FUNCTION && idx->data.function.head->type == EXPR_SYMBOL &&
            idx->data.function.head->data.symbol.name == SYM_List) {
            size_t m = idx->data.function.arg_count;
            Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
            for (size_t j = 0; j < m; j++) {
                Expr* sub = idx->data.function.args[j];
                Expr* v = assoc_part_single(expr, sub, rest, nrest);
                if (!v) {  /* e.g. positional out of range: report as missing */
                    Expr* margs[2] = { expr_new_string("KeyAbsent"), expr_copy(sub) };
                    v = expr_new_function(expr_new_symbol(SYM_Missing), margs, 2);
                }
                out[j] = v;
            }
            Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, m);
            free(out);
            return result;
        }
        return assoc_part_single(expr, idx, rest, nrest);
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

static Expr* assoc_part_single(Expr* assoc, Expr* idx, Expr** rest, size_t nrest) {
    size_t na = assoc->data.function.arg_count;
    Expr* lookup_key = NULL;   /* borrowed */
    bool positional = false;
    int64_t pos = 0;

    if (idx->type == EXPR_FUNCTION && idx->data.function.head->type == EXPR_SYMBOL &&
        idx->data.function.head->data.symbol.name == SYM_Key && idx->data.function.arg_count == 1) {
        lookup_key = idx->data.function.args[0];
    } else if (idx->type == EXPR_INTEGER) {
        positional = true;
        pos = idx->data.integer;
    } else {
        lookup_key = idx;  /* strings and other literal keys */
    }

    if (positional) {
        if (pos == 0) {
            /* assoc[[0]] gives the head Association (as for any expression);
             * an integer *key* must be requested with Key[k]. */
            Expr* head = expr_head(assoc);
            if (!head) return NULL;
            Expr* result = expr_part(head, rest, nrest);
            expr_free(head);
            return result;
        }
        if (pos < 0) pos = (int64_t)na + pos + 1;
        if (pos < 1 || pos > (int64_t)na) return NULL;
        Expr* rule = assoc->data.function.args[pos - 1];
        return expr_part(rule->data.function.args[1], rest, nrest);
    }
    for (size_t i = 0; i < na; i++) {
        Expr* rule = assoc->data.function.args[i];
        if (expr_eq(rule->data.function.args[0], lookup_key))
            return expr_part(rule->data.function.args[1], rest, nrest);
    }
    /* Key absent. */
    Expr* margs[2] = { expr_new_string("KeyAbsent"), expr_copy(lookup_key) };
    return expr_new_function(expr_new_symbol(SYM_Missing), margs, 2);
}

Expr* builtin_part(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) {
        return NULL;
    }

    Expr* expr = res->data.function.args[0];
    Expr** indices = res->data.function.args + 1;
    size_t nindices = res->data.function.arg_count - 1;

    /* NDArray: index the flat buffer directly. [[0]] (head extraction) falls
     * through to expr_part/expr_head, which reports NDArray as the head. Plain
     * integer subscripts index natively; Span/All/List positions degrade to the
     * general List Part via delist-and-reeval. */
    if (is_ndarray(expr) &&
        !(nindices > 0 && indices[0]->type == EXPR_INTEGER &&
          indices[0]->data.integer == 0)) {
        bool degrade = false;
        Expr* r = ndarray_part(expr, indices, nindices, &degrade);
        if (r) return r;
        if (degrade) return ndarray_delist_and_reeval(res);
        return NULL;  /* out-of-range / too many subscripts: leave unevaluated */
    }

    // Mathematica allows [[0]] even for atoms
    if (nindices > 0 && indices[0]->type == EXPR_INTEGER && indices[0]->data.integer == 0) {
        // Allowed
    } else {
        if (is_atomic(expr)) return NULL;
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
    Expr* part = expr_part(expr, indices, nindices);
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
     * or a single-position list); extract_single handles both forms. */

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

Expr* builtin_first(Expr* res) {
    if (res->type != EXPR_FUNCTION ||
        res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr* deflt = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    /* An NDArray is atomic, so is_atomic below rejects it and the call comes
     * back UNEVALUATED, while the identical List call works.  First is a
     * leading-axis slice of the buffer, so take it directly; delist_repack (which
     * materialises every element) remains the fallback for what it declines. */
    if (is_ndarray(arg)) {
        Expr* fast = ndstruct_head_tail(res, ND_FIRST);
        if (fast) return fast;
        return ndstruct_delist_repack(res, arg);
    }
    /* First[expr, def] returns def when expr has no elements (or is atomic);
     * with no default the 1-arg form is left unevaluated, as before. */
    if (is_atomic(arg) || arg->data.function.arg_count < 1)
        return deflt ? expr_copy(deflt) : NULL;
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
    if (is_ndarray(arg)) {
        Expr* fast = ndstruct_head_tail(res, ND_LAST);
        if (fast) return fast;
        return ndstruct_delist_repack(res, arg);
    }
    if (is_atomic(arg) || arg->data.function.arg_count < 1)
        return deflt ? expr_copy(deflt) : NULL;
    return expr_copy(first_last_element(arg, arg->data.function.args[arg->data.function.arg_count - 1]));
}

Expr* builtin_most(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
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

// Simple path comparison for sorting
static int compare_paths(const void* a, const void* b) {
    Expr* pa = *(Expr**)a;
    Expr* pb = *(Expr**)b;
    
    // Both are Lists of Integers (or single Integers converted to lists)
    size_t la = (pa->type == EXPR_FUNCTION) ? pa->data.function.arg_count : 1;
    size_t lb = (pb->type == EXPR_FUNCTION) ? pb->data.function.arg_count : 1;
    
    size_t n = (la < lb) ? la : lb;
    for (size_t i = 0; i < n; i++) {
        int64_t va = (pa->type == EXPR_FUNCTION) ? pa->data.function.args[i]->data.integer : pa->data.integer;
        int64_t vb = (pb->type == EXPR_FUNCTION) ? pb->data.function.args[i]->data.integer : pb->data.integer;
        if (va < vb) return 1; // Reverse sort
        if (va > vb) return -1;
    }
    if (la < lb) return 1;
    if (la > lb) return -1;
    return 0;
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
            Expr** sorted_pos = malloc(sizeof(Expr*) * pos->data.function.arg_count);
            memcpy(sorted_pos, pos->data.function.args, sizeof(Expr*) * pos->data.function.arg_count);
            qsort(sorted_pos, pos->data.function.arg_count, sizeof(Expr*), compare_paths);
            
            Expr* current = expr_copy(expr);
            for (size_t i = 0; i < pos->data.function.arg_count; i++) {
                Expr* next = NULL;
                if (sorted_pos[i]->type == EXPR_FUNCTION) {
                    next = insert_path(current, elem, sorted_pos[i]->data.function.args, sorted_pos[i]->data.function.arg_count);
                } else {
                    Expr* p = sorted_pos[i];
                    next = insert_path(current, elem, &p, 1);
                }
                expr_free(current);
                current = next;
            }
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

Expr* builtin_insert(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
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
        Expr* rule = na[t];
        Expr* nv = delete_path(rule->data.function.args[1], path + 1, path_len - 1);
        expr_free(rule->data.function.args[1]);
        rule->data.function.args[1] = nv;
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
            Expr** sorted_pos = malloc(sizeof(Expr*) * pos->data.function.arg_count);
            memcpy(sorted_pos, pos->data.function.args, sizeof(Expr*) * pos->data.function.arg_count);
            qsort(sorted_pos, pos->data.function.arg_count, sizeof(Expr*), compare_paths);
            
            Expr* current = expr_copy(expr);
            for (size_t i = 0; i < pos->data.function.arg_count; i++) {
                Expr* next = NULL;
                if (sorted_pos[i]->type == EXPR_FUNCTION) {
                    next = delete_path(current, sorted_pos[i]->data.function.args, sorted_pos[i]->data.function.arg_count);
                } else {
                    Expr* p = sorted_pos[i];
                    next = delete_path(current, &p, 1);
                }
                expr_free(current);
                current = next;
            }
            free(sorted_pos);
            return current;
        }
        
        // Case 2: Single path
        return delete_path(expr, pos->data.function.args, pos->data.function.arg_count);
    }
    
    // Case 1: Single integer
    if (pos->type == EXPR_INTEGER) {
        return delete_path(expr, &pos, 1);
    }
    
    return expr_copy(expr);
}

Expr* builtin_delete(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count == 1) {
        return NULL; // Operator form, stays unevaluated
    } else if (res->data.function.arg_count == 2) {
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
