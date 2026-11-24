#include "list_common.h"
#include "transpose.h"
#include "ndarray.h"
#include "ndstruct.h"
#include "message.h"
#include "linkedlist.h"

static int get_array_dimensions(Expr* e, int64_t* dims, const char* head_name) {
    if (!head_is(e, intern_symbol(head_name))) {
        return 0;
    }
    dims[0] = (int64_t)e->data.function.arg_count;
    if (dims[0] == 0) return 1;
    int64_t sub_dims[64];
    int depth = get_array_dimensions(e->data.function.args[0], sub_dims, head_name);
    for (size_t i = 1; i < e->data.function.arg_count; i++) {
        int64_t cur_dims[64];
        if (get_array_dimensions(e->data.function.args[i], cur_dims, head_name) != depth) return 1;
        for (int j = 0; j < depth; j++) if (cur_dims[j] != sub_dims[j]) return 1;
    }
    for (int i = 0; i < depth; i++) dims[i + 1] = sub_dims[i];
    return depth + 1;
}

/* For the 1-arg Transpose the top two levels must be rectangular -- a List of
 * rows, each a List of the same length -- while entries at the third level and
 * below are opaque leaves. Returns 1 and sets nrows, ncols on success, else 0.
 * (A8: Mathematica's 1-arg Transpose swaps only the first two levels, so a list
 * of equal-length rows whose entries have differing deeper shapes is valid.) */
static int get_top2_dims(Expr* e, const char* head_name,
                         int64_t* nrows, int64_t* ncols) {
    const char* h = intern_symbol(head_name);
    if (!head_is(e, h)) return 0;
    int64_t r = (int64_t)e->data.function.arg_count;
    if (r == 0) return 0;
    Expr* row0 = e->data.function.args[0];
    if (!head_is(row0, h)) return 0;              /* not a list of lists */
    int64_t c = (int64_t)row0->data.function.arg_count;
    for (int64_t i = 1; i < r; i++) {
        Expr* row = e->data.function.args[i];
        if (!head_is(row, h)) return 0;
        if ((int64_t)row->data.function.arg_count != c) return 0;  /* ragged top level */
    }
    *nrows = r; *ncols = c;
    return 1;
}

/* Transpose only the top two levels: result[j][i] = list[i][j], with each
 * list[i][j] kept whole (opaque). */
static Expr* transpose_top2(Expr* list, const char* head, int64_t nr, int64_t nc) {
    Expr** cols = malloc(sizeof(Expr*) * (size_t)(nc ? nc : 1));
    for (int64_t j = 0; j < nc; j++) {
        Expr** rowj = malloc(sizeof(Expr*) * (size_t)nr);
        for (int64_t i = 0; i < nr; i++)
            rowj[i] = expr_copy(list->data.function.args[i]->data.function.args[j]);
        cols[j] = expr_new_function(expr_new_symbol(head), rowj, (size_t)nr);
        free(rowj);
    }
    Expr* result = expr_new_function(expr_new_symbol(head), cols, (size_t)nc);
    free(cols);
    return result;
}

static Expr* get_element_at(Expr* e, int64_t* indices, size_t depth) {
    Expr* curr = e;
    for (size_t i = 0; i < depth; i++) {
        curr = curr->data.function.args[indices[i]];
    }
    return curr;
}

static Expr* build_transposed(const char* head, int64_t* out_dims, size_t out_depth, int64_t* out_indices_base, int64_t* current_out_indices, 
                             int64_t* in_indices, int64_t* perm, size_t in_depth, Expr* original) {
    if (out_depth == 0) {
        // level k in list is Subscript[n, k]-th level in result.
        // So in_indices[k] = out_indices_base[perm[k] - 1]
        for (size_t k = 0; k < in_depth; k++) {
            in_indices[k] = out_indices_base[perm[k] - 1];
        }
        return expr_copy(get_element_at(original, in_indices, in_depth));
    }

    size_t len = (size_t)out_dims[0];
    Expr** args = malloc(sizeof(Expr*) * len);
    for (size_t i = 0; i < len; i++) {
        current_out_indices[0] = (int64_t)i;
        args[i] = build_transposed(head, out_dims + 1, out_depth - 1, out_indices_base, current_out_indices + 1, in_indices, perm, in_depth, original);
    }
    Expr* res = expr_new_function(expr_new_symbol(head), args, len);
    free(args);
    return res;
}

static void tr_msg(const char* tag, const char* fmt, const Expr* a, const Expr* b) {
    char* as = a ? expr_to_string_outputform((Expr*)a) : NULL;
    char* bs = b ? expr_to_string_outputform((Expr*)b) : NULL;
    mth_message("Transpose", tag, fmt, as ? as : "", bs ? bs : "");
    free(as);
    free(bs);
}

/* Is no element of `list` a List (a vector, whose transpose is itself)? */
static bool is_flat_vector(const Expr* list) {
    for (size_t i = 0; i < list->data.function.arg_count; i++) {
        const Expr* x = list->data.function.args[i];
        if (x->type == EXPR_NDARRAY || head_is(x, SYM_List)) return false;
    }
    return true;
}

/* Read and check the level specification of Transpose[list, spec] against an
 * array of depth `in_depth` (dims `in_dims`): perm[k] (1-based) is the result
 * level of input level k. Returns the spec's length, 0 when the call must stay
 * unevaluated (after the message that says why), or -1 for a spec longer than
 * the array that is the identity (the array comes back unchanged). */
static int transpose_perm_spec(const Expr* spec, int in_depth, const int64_t* in_dims,
                               const Expr* call, int64_t* perm) {
    if (!head_is(spec, SYM_List)) {
        if (head_is(spec, SYM_TwoWayRule)) return 0;
        tr_msg("permspec", "Invalid permutation specification found at position 2 in %s.", call, NULL);
        return 0;
    }
    size_t L = spec->data.function.arg_count;
    if (L == 0 || L > 64) return 0;
    int64_t sp[64];
    for (size_t i = 0; i < L; i++) {
        const Expr* x = spec->data.function.args[i];
        if (x->type != EXPR_INTEGER || x->data.integer < 1) {
            tr_msg("perm1", "Entry %s in permutation %s is not a positive machine integer.", x, spec);
            return 0;
        }
        sp[i] = x->data.integer;
    }
    for (size_t i = 0; i < L; i++)
        if (sp[i] > (int64_t)L) {
            char* xs = expr_to_string_outputform(spec->data.function.args[i]);
            char* ss = expr_to_string_outputform((Expr*)spec);
            mth_message("Transpose", "perm2",
                        "Entry %s in %s is out of bounds for a permutation of length %d.",
                        xs ? xs : "", ss ? ss : "", (int)L);
            free(xs);
            free(ss);
            return 0;
        }
    int64_t mx = 0;
    for (size_t i = 0; i < L; i++) if (sp[i] > mx) mx = sp[i];
    for (int64_t lvl = 1; lvl <= mx; lvl++) {
        bool hit = false;
        for (size_t i = 0; i < L && !hit; i++) hit = sp[i] == lvl;
        if (!hit) {
            char* ss = expr_to_string_outputform((Expr*)spec);
            mth_message("Transpose", "newdims",
                        "Level rearrangement %s does not specify destination for level %d.",
                        ss ? ss : "", (int)lvl);
            free(ss);
            return 0;
        }
    }
    if ((int)L > in_depth) {
        bool ident = true;
        for (size_t i = 0; i < L; i++) if (sp[i] != (int64_t)i + 1) ident = false;
        if (ident) return -1;
        Expr** d = malloc(sizeof(Expr*) * (size_t)(in_depth ? in_depth : 1));
        for (int i = 0; i < in_depth; i++) d[i] = expr_new_integer(in_dims[i]);
        Expr* dl = expr_new_function(expr_new_symbol(SYM_List), d, (size_t)in_depth);
        free(d);
        tr_msg("tperm", "Permutation %s is longer than the dimensions %s of the expression.", spec, dl);
        expr_free(dl);
        return 0;
    }
    /* levels past the spec stay where they are */
    for (int i = 0; i < in_depth; i++) perm[i] = i < (int)L ? sp[i] : i + 1;
    /* levels sent to one result level form a diagonal: equal lengths only */
    for (int i = 0; i < in_depth; i++)
        for (int j = i + 1; j < in_depth; j++)
            if (perm[i] == perm[j] && in_dims[i] != in_dims[j]) {
                tr_msg("diagnl", "Level rearrangement %s would require collapsing dimensions of "
                       "unequal lengths.", spec, NULL);
                return 0;
            }
    return (int)L;
}

/* ---- a nest with a LinkedList level --------------------------------------
 * A LinkedList is the List of its elements, so these are the List walks above
 * over element views. The answer keeps a LinkedList at every level the
 * argument had one (along the first-element spine). */

static int64_t seq_n(const Expr* e) {
    return e->type == EXPR_LINKEDLIST ? ll_length(e->data.llist)
                                      : (int64_t)e->data.function.arg_count;
}

/* get_array_dimensions over Lists and LinkedLists: the depth to which `e` is
 * rectangular. */
static int ll_array_dimensions(const Expr* e, int64_t* dims) {
    if (!ll_is_seq(e)) return 0;
    dims[0] = seq_n(e);
    if (dims[0] == 0) return 1;
    if (e->type == EXPR_LINKEDLIST && ll_kind(e->data.llist) != LL_BOXED) return 1;
    int64_t sub[64], cur[64];
    bool own;
    Expr* el = ll_seq_view(e, 0, &own);
    int depth = ll_array_dimensions(el, sub);
    ll_elem_done(el, own);
    for (int64_t i = 1; i < dims[0]; i++) {
        el = ll_seq_view(e, i, &own);
        int d = ll_array_dimensions(el, cur);
        ll_elem_done(el, own);
        if (d != depth) return 1;
        for (int j = 0; j < depth; j++) if (cur[j] != sub[j]) return 1;
    }
    for (int i = 0; i < depth; i++) dims[i + 1] = sub[i];
    return depth + 1;
}

static Expr* transpose_ll(Expr* res) {
    Expr* list = res->data.function.args[0];
    size_t argc = res->data.function.arg_count;
    int64_t in_dims[64];
    int in_depth = ll_array_dimensions(list, in_dims);
    bool vector = in_depth == 1;
    if (vector) {
        for (int64_t i = 0, n = seq_n(list); i < n && vector; i++) {
            bool own;
            Expr* el = ll_seq_view(list, i, &own);
            if (ll_is_seq(el) || el->type == EXPR_NDARRAY) vector = false;
            ll_elem_done(el, own);
        }
    }
    if (vector && list->type == EXPR_LINKEDLIST) {
        /* a vector is its own transpose: a copy, as the copying heads give */
        const Expr* sp = argc == 2 ? res->data.function.args[1] : NULL;
        if (!sp || (head_is(sp, SYM_List) && sp->data.function.arg_count == 1 &&
                    sp->data.function.args[0]->type == EXPR_INTEGER &&
                    sp->data.function.args[0]->data.integer == 1))
            return expr_new_linkedlist(ll_copy(list->data.llist));
    }
    if (in_depth < 2 && argc == 1) {
        if (vector) return expr_copy(list);
        /* the top two levels, deeper entries opaque (A8) */
        int64_t nr = seq_n(list), nc = -1;
        bool ok = nr > 0;
        for (int64_t i = 0; i < nr && ok; i++) {
            bool own;
            Expr* row = ll_seq_view(list, i, &own);
            if (!ll_is_seq(row)) ok = false;
            else if (nc < 0) nc = seq_n(row);
            else if (seq_n(row) != nc) ok = false;
            ll_elem_done(row, own);
        }
        if (!ok) {
            tr_msg("nmtx", "The first two levels of %s cannot be transposed.", list, NULL);
            return NULL;
        }
        in_depth = 2;
        in_dims[0] = nr;
        in_dims[1] = nc;
    }
    if (in_depth < 2 && !vector) return NULL;

    int64_t perm[64];
    if (argc == 1) {
        perm[0] = 2; perm[1] = 1;
        for (int i = 2; i < in_depth; i++) perm[i] = i + 1;
    } else {
        int L = transpose_perm_spec(res->data.function.args[1], in_depth, in_dims, res, perm);
        if (L <= 0) return L < 0 ? (list->type == EXPR_LINKEDLIST
                                     ? expr_new_linkedlist(ll_copy(list->data.llist))
                                     : expr_copy(list)) : NULL;
    }
    int out_depth = 0;
    for (int i = 0; i < in_depth; i++) if (perm[i] > out_depth) out_depth = (int)perm[i];
    int64_t out_dims[64];
    for (int i = 0; i < out_depth; i++) out_dims[i] = -1;
    for (int i = 0; i < in_depth; i++) {
        int t = (int)perm[i] - 1;
        if (out_dims[t] == -1 || in_dims[i] < out_dims[t]) out_dims[t] = in_dims[i];
    }
    int64_t nin = 1, nout = 1;
    for (int i = 0; i < in_depth; i++) nin *= in_dims[i];
    for (int i = 0; i < out_depth; i++) nout *= out_dims[i];
    /* each input cell read on demand (ll_leaf), never all of them at once */
    LLLeaves L;
    if (!ll_leaves_nest(&L, list, in_depth, in_dims)) return NULL;
    (void)nin;
    Expr** out = malloc(sizeof(Expr*) * (size_t)(nout ? nout : 1));
    int64_t stride[64];
    { int64_t s = 1; for (int i = in_depth - 1; i >= 0; i--) { stride[i] = s; s *= in_dims[i]; } }
    int64_t oi[64];
    for (int64_t o = 0; o < nout; o++) {
        int64_t rem = o;
        for (int i = out_depth - 1; i >= 0; i--) { oi[i] = rem % out_dims[i]; rem /= out_dims[i]; }
        int64_t src = 0;
        for (int k = 0; k < in_depth; k++) src += oi[perm[k] - 1] * stride[k];
        bool own;
        Expr* e = ll_leaf(&L, src, &own);
        out[o] = own ? e : expr_copy(e);
    }
    bool lvl[64];
    ll_tensor_levels(list, out_depth < in_depth ? out_depth : in_depth, lvl);
    for (int i = in_depth; i < out_depth; i++) lvl[i] = false;
    Expr* r = ll_tensor_build(out_dims, out_depth, out, lvl);
    free(out);
    return r;
}

Expr* builtin_transpose(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* list = res->data.function.args[0];
    size_t argc = res->data.function.arg_count;

    /* Transpose[list, m <-> n] swaps levels m and n: the permutation that
     * exchanges them, as a list spec (audit Transpose#1). */
    if (argc == 2 && head_is(res->data.function.args[1], SYM_TwoWayRule) &&
        res->data.function.args[1]->data.function.arg_count == 2) {
        const Expr* a = res->data.function.args[1]->data.function.args[0];
        const Expr* b = res->data.function.args[1]->data.function.args[1];
        if (a->type != EXPR_INTEGER || b->type != EXPR_INTEGER ||
            a->data.integer < 1 || b->data.integer < 1 ||
            a->data.integer > 64 || b->data.integer > 64) {
            tr_msg("permspec", "Invalid permutation specification found at position 2 in %s.", res, NULL);
            return NULL;
        }
        int64_t m = a->data.integer, n = b->data.integer, L = m > n ? m : n;
        Expr** pa = malloc(sizeof(Expr*) * (size_t)L);
        for (int64_t k = 1; k <= L; k++)
            pa[k - 1] = expr_new_integer(k == m ? n : k == n ? m : k);
        Expr* spec = expr_new_function(expr_new_symbol(SYM_List), pa, (size_t)L);
        free(pa);
        Expr* ta[2] = { expr_copy(list), spec };
        Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Transpose), ta, 2));
        if (head_is(r, SYM_Transpose)) { expr_free(r); return NULL; }
        return r;
    }

    /* The transpose of a vector is the vector itself. */
    if (is_ndarray(list) && list->data.ndarray.rank == 1) {
        const Expr* sp = argc == 2 ? res->data.function.args[1] : NULL;
        if (!sp || (head_is(sp, SYM_List) && sp->data.function.arg_count == 1 &&
                    sp->data.function.args[0]->type == EXPR_INTEGER &&
                    sp->data.function.args[0]->data.integer == 1))
            return expr_copy(list);
    }

    /* Transpose[ndarray] (rank 2) swaps axes on the buffer (see ndstruct.c). */
    if (is_ndarray(list)) return ndstruct_transpose(res);

    /* A LinkedList level: the same walk over element views (above). */
    if (ll_has_ll_level1(list)) return transpose_ll(res);

    if (list->type != EXPR_FUNCTION || list->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* head = list->data.function.head->data.symbol.name;

    int64_t in_dims[64];
    int in_depth = get_array_dimensions(list, in_dims, head);
    if (in_depth < 1) return NULL;
    bool vector = in_depth == 1 && head == SYM_List && is_flat_vector(list);
    if (in_depth < 2 && argc == 1) {
        /* A vector is its own transpose. Ragged deeper entries: the 1-arg form
         * still transposes the top two levels when they are rectangular,
         * treating level-2 entries as opaque (A8). */
        if (vector) return expr_copy(list);
        int64_t nr, nc;
        if (get_top2_dims(list, head, &nr, &nc)) return transpose_top2(list, head, nr, nc);
        if (head == SYM_List)
            tr_msg("nmtx", "The first two levels of %s cannot be transposed.", list, NULL);
        return NULL;
    }
    /* The n-arg permutation form genuinely needs full rectangularity. */
    if (in_depth < 2 && !vector) return NULL;

    int64_t* perm = malloc(sizeof(int64_t) * in_depth);
    if (argc == 1) {
        perm[0] = 2; perm[1] = 1;
        for (int i = 2; i < in_depth; i++) perm[i] = i + 1;
    } else {
        int L = transpose_perm_spec(res->data.function.args[1], in_depth, in_dims, res, perm);
        if (L <= 0) {
            free(perm);
            return L < 0 ? expr_copy(list) : NULL;
        }
    }

    int out_depth = 0;
    for (int i = 0; i < in_depth; i++) if (perm[i] > out_depth) out_depth = (int)perm[i];

    int64_t* out_dims = malloc(sizeof(int64_t) * out_depth);
    for (int i = 0; i < out_depth; i++) out_dims[i] = -1;

    for (int i = 0; i < in_depth; i++) {
        int target_idx = (int)perm[i] - 1;
        if (out_dims[target_idx] == -1 || in_dims[i] < out_dims[target_idx]) {
            out_dims[target_idx] = in_dims[i];
        }
    }

    int64_t* out_indices_base = calloc(out_depth, sizeof(int64_t));
    int64_t* in_indices = malloc(sizeof(int64_t) * in_depth);

    Expr* result = build_transposed(head, out_dims, (size_t)out_depth, out_indices_base, out_indices_base, in_indices, perm, (size_t)in_depth, list);

    free(perm); free(out_dims); free(out_indices_base); free(in_indices);
    return result;
}

/* ConjugateTranspose[m] is equivalent to Conjugate[Transpose[m]];
 * ConjugateTranspose[m, spec] is equivalent to Conjugate[Transpose[m, spec]].
 * For a 1-D vector input (1-arg form) the result keeps the vector shape and
 * the entries are conjugated, matching Mathematica's behaviour. Conjugate is
 * Listable, so wrapping the (possibly nested) transposed array in Conjugate
 * lets the evaluator thread conjugation down to every leaf at fixed point. */
Expr* builtin_conjugate_transpose(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    Expr* m = res->data.function.args[0];

    /* ConjugateTranspose[ndarray, ...] = Conjugate[Transpose[ndarray, ...]]:
     * Transpose has a native NDArray path (ndstruct_transpose) and Conjugate is
     * an NDArray element-wise kernel, so build the composite and let both fire.
     *
     * RANK 1 IS NOT A TRANSPOSE. The one-argument form on a vector keeps the
     * shape and only conjugates (the branch below does this for a List), and
     * Transpose declines a rank-1 buffer -- so the composite came back as an
     * unevaluated Conjugate[Transpose[v]]. That is why this head sat in
     * tools/check_packed_aware.py's EXEMPT table rather than on pack.c's AWARE
     * list, and why a packed 1000x1000 ConjugateTranspose cost 430 ms against
     * NumPy's 2.20 ms: the gate materialised 10^6 elements first. Handling the
     * vector here is what lets the head opt in. */
    if (is_ndarray(m)) {
        if (argc == 1 && m->data.ndarray.rank < 2) {
            if (m->data.ndarray.rank != 1) return NULL;
            Expr* ca[1] = { expr_copy(m) };
            return eval_and_free(expr_new_function(expr_new_symbol(SYM_Conjugate), ca, 1));
        }
        Expr** tr = malloc(sizeof(Expr*) * argc);
        for (size_t i = 0; i < argc; i++) tr[i] = expr_copy(res->data.function.args[i]);
        Expr* tcall = eval_and_free(expr_new_function(expr_new_symbol(SYM_Transpose), tr, argc));
        free(tr);
        /* Transpose could not reduce (a bad spec): surface the unevaluated
         * ConjugateTranspose, as the List path below does, rather than a
         * Conjugate[Transpose[...]] wrapper. */
        if (head_is(tcall, SYM_Transpose)) { expr_free(tcall); return NULL; }
        /* Conjugation of REAL data is the identity, and running it anyway is a
         * second full pass over the buffer -- 6 ms of the 7.8 ms a 1000x1000
         * took, for a kernel that copies each element to itself. */
        if (is_ndarray(tcall) && !ndt_is_complex(tcall->data.ndarray.dtype))
            return tcall;
        Expr* ca[1] = { tcall };
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Conjugate), ca, 1));
    }

    int64_t dims[64];
    int depth = 0;
    if (ll_has_ll_level1(m)) {
        depth = ll_array_dimensions(m, dims);
    } else if (m->type == EXPR_FUNCTION && m->data.function.head->type == EXPR_SYMBOL &&
        m->data.function.head->data.symbol.name == SYM_List) {
        depth = get_array_dimensions(m, dims, "List");
    }

    /* Not a (rectangular) List: leave the call unevaluated so that symbolic
     * matrices like ConjugateTranspose[A] stay intact. */
    if (depth == 0) return NULL;

    /* 1-arg form on a ragged list (depth 1, but with List entries): the top
     * two levels are transposed when they are rectangular, else
     * ConjugateTranspose::nmtx (audit ConjugateTranspose#1: {{1, I}, {3}} was
     * conjugated as if it were a vector). */
    if (argc == 1 && depth == 1 && m->type == EXPR_FUNCTION && !is_flat_vector(m)) {
        int64_t nr, nc;
        if (!get_top2_dims(m, SYM_List, &nr, &nc)) {
            char* ms = expr_to_string_outputform(m);
            mth_message("ConjugateTranspose", "nmtx",
                        "The first two levels of %s cannot be transposed.", ms ? ms : "");
            free(ms);
            return NULL;
        }
        depth = 2;
    }

    /* 1-arg form on a 1-D vector: conjugate elementwise, keep the shape. */
    if (argc == 1 && depth == 1) {
        Expr** conj_args = malloc(sizeof(Expr*) * 1);
        conj_args[0] = expr_copy(m);
        Expr* conj = expr_new_function(expr_new_symbol(SYM_Conjugate), conj_args, 1);
        free(conj_args);
        return eval_and_free(conj);
    }

    /* 1-arg form on lower depth (depth < 2) we cannot transpose: leave it. */
    if (argc == 1 && depth < 2) return NULL;

    /* Build Transpose[m] or Transpose[m, spec]. */
    Expr** tr_args = malloc(sizeof(Expr*) * argc);
    tr_args[0] = expr_copy(m);
    if (argc == 2) tr_args[1] = expr_copy(res->data.function.args[1]);
    Expr* transposed_call = expr_new_function(expr_new_symbol(SYM_Transpose), tr_args, argc);
    free(tr_args);

    Expr* transposed = eval_and_free(transposed_call);

    /* If Transpose could not reduce (e.g. invalid spec), surface the
     * unevaluated ConjugateTranspose rather than a spurious
     * Conjugate[Transpose[...]] wrapper. */
    if (head_is(transposed, SYM_Transpose)) {
        expr_free(transposed);
        return NULL;
    }

    Expr** conj_args = malloc(sizeof(Expr*) * 1);
    conj_args[0] = transposed;
    Expr* conj = expr_new_function(expr_new_symbol(SYM_Conjugate), conj_args, 1);
    free(conj_args);
    return eval_and_free(conj);
}
