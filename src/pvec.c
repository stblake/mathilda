/* pvec.c -- persistent vector.  See pvec.h for the contract.
 *
 * Layout: PVec { count, shift, root, tail }.  `root` is a trie of branch nodes
 * (PV_WIDTH children) over full leaves of PV_WIDTH elements; `tail` is the
 * rightmost, possibly partial leaf, kept outside the trie so appends touch one
 * small array.  Everything except the root/tail pointers lives in refcounted
 * PNodes, which is what makes structural sharing safe: a node with refcount 1 is
 * private to one path and may be edited in place, a node with refcount > 1 is
 * copied before it is edited (path copying).
 *
 * All slots are one 8-byte word (int64_t, double, Expr*, or a child pointer),
 * so one node type serves branches and leaves of every kind.
 */
#include "pvec.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

typedef union PWord {
    int64_t i;
    double  d;
    Expr*   e;
    struct PNode* c;
} PWord;

typedef struct PNode {
    uint32_t rc;
    uint8_t  leaf;       /* 1 = leaf of elements, 0 = branch of children */
    uint8_t  n;          /* leaf: live elements; branch: unused */
    PWord    w[PV_WIDTH];
} PNode;

struct PVec {
    uint32_t rc;
    PVKind   kind;
    int64_t  count;
    int      shift;      /* trie height * PV_BITS; PV_BITS for a one-level trie */
    PNode*   root;       /* branch node, never NULL once count > PV_WIDTH */
    PNode*   tail;       /* leaf, NULL while empty */
    uint64_t id;
};

static uint64_t g_next_id = 1;

/* ---- nodes -------------------------------------------------------------- */
static PNode* node_new(bool leaf) {
    PNode* n = (PNode*)calloc(1, sizeof(PNode));
    n->rc = 1;
    n->leaf = leaf ? 1 : 0;
    return n;
}
static PNode* node_retain(PNode* n) { if (n) n->rc++; return n; }

static void node_release(PNode* n, PVKind kind) {
    if (!n || --n->rc) return;
    if (n->leaf) {
        if (kind == PV_BOXED)
            for (int i = 0; i < n->n; i++) expr_free(n->w[i].e);
    } else {
        for (int i = 0; i < PV_WIDTH; i++)
            if (n->w[i].c) node_release(n->w[i].c, kind);
    }
    free(n);
}

/* An exclusive version of n: n itself when held once, else a copy that retains
 * what it points to (children, or boxed elements).  Consumes the caller's
 * reference to n. */
static PNode* node_own(PNode* n, PVKind kind) {
    if (n->rc == 1) return n;
    PNode* c = node_new(n->leaf);
    c->n = n->n;
    if (n->leaf) {
        memcpy(c->w, n->w, sizeof c->w);
        if (kind == PV_BOXED)
            for (int i = 0; i < n->n; i++) c->w[i].e = expr_copy(n->w[i].e);
    } else {
        for (int i = 0; i < PV_WIDTH; i++) c->w[i].c = node_retain(n->w[i].c);
    }
    n->rc--;                      /* n is shared, so this cannot reach zero */
    return c;
}

/* ---- vector shell ------------------------------------------------------- */
PVec* pv_new(void) {
    PVec* v = (PVec*)calloc(1, sizeof(PVec));
    v->rc = 1;
    v->shift = PV_BITS;
    v->id = g_next_id++;
    return v;
}
PVec* pv_retain(PVec* v) { v->rc++; return v; }
bool pv_exclusive(const PVec* v) { return v->rc == 1; }

void pv_release(PVec* v) {
    if (!v || --v->rc) return;
    node_release(v->root, v->kind);
    node_release(v->tail, v->kind);
    free(v);
}

/* An exclusive PVec: v itself when held once, else a fresh shell sharing v's
 * nodes.  Consumes the caller's reference. */
static PVec* pv_own(PVec* v) {
    if (v->rc == 1) return v;
    PVec* w = (PVec*)malloc(sizeof(PVec));
    *w = *v;
    w->rc = 1;
    w->id = g_next_id++;
    node_retain(w->root);
    node_retain(w->tail);
    v->rc--;
    return w;
}

int64_t pv_length(const PVec* v) { return v->count; }
PVKind  pv_kind(const PVec* v) { return v->kind; }
uint64_t pv_id(const PVec* v) { return v->id; }
bool pv_same(const PVec* a, const PVec* b) {
    return a == b || (a->count == b->count && a->root == b->root && a->tail == b->tail &&
                      a->kind == b->kind);
}

static int64_t tail_off(const PVec* v) {
    return v->count < PV_WIDTH ? 0 : ((v->count - 1) >> PV_BITS) << PV_BITS;
}

/* A node is shared when it, or any ancestor, is held by more than one
 * reference: a whole subtree hangs off a shared branch even though its own nodes
 * have refcount 1. */
static void stats_rec(const PNode* n, bool under_shared, int64_t* nodes, int64_t* shared) {
    if (!n) return;
    bool sh = under_shared || n->rc > 1;
    (*nodes)++;
    if (sh) (*shared)++;
    if (!n->leaf)
        for (int i = 0; i < PV_WIDTH; i++) stats_rec(n->w[i].c, sh, nodes, shared);
}
void pv_node_stats(const PVec* v, int64_t* nodes, int64_t* shared) {
    *nodes = 0; *shared = 0;
    stats_rec(v->root, false, nodes, shared);
    stats_rec(v->tail, false, nodes, shared);
}

/* ---- element <-> word --------------------------------------------------- */
static PVKind classify(const Expr* e) {
    if (e->type == EXPR_INTEGER) return PV_I64;
    if (e->type == EXPR_REAL) return PV_REAL;
    return PV_BOXED;
}
static Expr* word_to_expr(PWord w, PVKind k) {
    switch (k) {
        case PV_I64:  return expr_new_integer(w.i);
        case PV_REAL: return expr_new_real(w.d);
        default:      return expr_copy(w.e);
    }
}
/* Consumes e. */
static PWord expr_to_word(Expr* e, PVKind k) {
    PWord w;
    switch (k) {
        case PV_I64:  w.i = e->data.integer; expr_free(e); break;
        case PV_REAL: w.d = e->data.real;    expr_free(e); break;
        default:      w.e = e; break;
    }
    return w;
}

/* ---- navigation ---------------------------------------------------------- */
static const PNode* leaf_for(const PVec* v, int64_t i) {
    if (i >= tail_off(v)) return v->tail;
    const PNode* n = v->root;
    for (int level = v->shift; level > 0; level -= PV_BITS)
        n = n->w[(i >> level) & (PV_WIDTH - 1)].c;
    return n;
}

Expr* pv_get(const PVec* v, int64_t i) {
    const PNode* l = leaf_for(v, i);
    return word_to_expr(l->w[i & (PV_WIDTH - 1)], v->kind);
}
bool pv_get_i64(const PVec* v, int64_t i, int64_t* out) {
    if (v->kind != PV_I64) return false;
    *out = leaf_for(v, i)->w[i & (PV_WIDTH - 1)].i;
    return true;
}
bool pv_get_real(const PVec* v, int64_t i, double* out) {
    if (v->kind != PV_REAL) return false;
    *out = leaf_for(v, i)->w[i & (PV_WIDTH - 1)].d;
    return true;
}
int64_t pv_span_at(const PVec* v, int64_t i, const void** p) {
    const PNode* l = leaf_for(v, i);
    int off = (int)(i & (PV_WIDTH - 1));
    *p = &l->w[off];
    return (int64_t)l->n - off;
}

/* ---- kind conversion (typed -> boxed, once) ------------------------------ */
/* Rewrite every node of an exclusive v as boxed.  Shared nodes must be copied
 * first, so the tree is rebuilt: O(n), and it happens at most once per vector
 * lineage (boxed never goes back). */
static PNode* box_copy(const PNode* n, PVKind from) {
    PNode* c = node_new(n->leaf);
    c->n = n->n;
    if (n->leaf) {
        for (int i = 0; i < n->n; i++)
            c->w[i].e = from == PV_I64 ? expr_new_integer(n->w[i].i) : expr_new_real(n->w[i].d);
    } else {
        for (int i = 0; i < PV_WIDTH; i++)
            if (n->w[i].c) c->w[i].c = box_copy(n->w[i].c, from);
    }
    return c;
}
static PVec* pv_to_boxed(PVec* v) {
    if (v->kind == PV_BOXED || v->kind == PV_EMPTY) { v->kind = PV_BOXED; return v; }
    PVKind from = v->kind;
    PNode* r = v->root ? box_copy(v->root, from) : NULL;
    PNode* t = v->tail ? box_copy(v->tail, from) : NULL;
    node_release(v->root, from);
    node_release(v->tail, from);
    v->root = r; v->tail = t; v->kind = PV_BOXED;
    return v;
}

/* ---- set ----------------------------------------------------------------- */
static PNode* assoc_rec(PNode* n, int level, int64_t i, PWord val, PVKind kind) {
    n = node_own(n, kind);
    if (level == 0) {
        int j = (int)(i & (PV_WIDTH - 1));
        if (kind == PV_BOXED) expr_free(n->w[j].e);
        n->w[j] = val;
    } else {
        int j = (int)((i >> level) & (PV_WIDTH - 1));
        n->w[j].c = assoc_rec(n->w[j].c, level - PV_BITS, i, val, kind);
    }
    return n;
}

PVec* pv_set(PVec* v, int64_t i, Expr* e) {
    v = pv_own(v);
    PVKind k = classify(e);
    if (v->kind != k && v->kind != PV_BOXED) v = pv_to_boxed(v);
    PWord w = expr_to_word(e, v->kind);
    if (i >= tail_off(v)) {
        v->tail = node_own(v->tail, v->kind);
        int j = (int)(i & (PV_WIDTH - 1));
        if (v->kind == PV_BOXED) expr_free(v->tail->w[j].e);
        v->tail->w[j] = w;
    } else {
        v->root = assoc_rec(v->root, v->shift, i, w, v->kind);
    }
    v->id = g_next_id++;
    return v;
}

/* ---- push ---------------------------------------------------------------- */
/* Hang `leaf` (a full leaf, ownership transferred) at position `i` of the trie
 * under `n` at `level`, creating branch nodes on the way. */
static PNode* push_tail(PNode* n, int level, int64_t i, PNode* leaf, PVKind kind) {
    n = node_own(n, kind);
    int j = (int)((i >> level) & (PV_WIDTH - 1));
    if (level == PV_BITS) {
        n->w[j].c = leaf;
    } else if (n->w[j].c) {
        n->w[j].c = push_tail(n->w[j].c, level - PV_BITS, i, leaf, kind);
    } else {
        /* a fresh path of single-child branches down to the leaf */
        PNode* child = leaf;
        for (int l = PV_BITS; l < level; l += PV_BITS) {
            PNode* b = node_new(false);
            b->w[0].c = child;                 /* fresh subtree: first leaf is index 0 */
            (void)l;
            child = b;
        }
        n->w[j].c = child;
    }
    return n;
}

PVec* pv_push(PVec* v, Expr* e) {
    v = pv_own(v);
    PVKind k = classify(e);
    if (v->kind == PV_EMPTY) v->kind = k;
    else if (v->kind != k && v->kind != PV_BOXED) v = pv_to_boxed(v);
    PWord w = expr_to_word(e, v->kind);

    if (!v->tail) v->tail = node_new(true);
    int64_t toff = tail_off(v);
    int tl = (int)(v->count - toff);

    if (tl < PV_WIDTH) {
        v->tail = node_own(v->tail, v->kind);
        v->tail->w[tl] = w;
        v->tail->n = (uint8_t)(tl + 1);
    } else {
        /* tail is full: push it into the trie, start a new tail */
        PNode* full = v->tail;                       /* ownership moves to the trie */
        PNode* nt = node_new(true);
        nt->w[0] = w;
        nt->n = 1;
        if (!v->root) {
            v->root = node_new(false);
            v->root->w[0].c = full;
            v->shift = PV_BITS;
        } else if ((v->count >> PV_BITS) > ((int64_t)1 << v->shift)) {
            /* root overflow: grow one level */
            PNode* nr = node_new(false);
            nr->w[0].c = v->root;
            PNode* child = full;
            for (int l = 0; l < v->shift; l += PV_BITS) {   /* fresh single-child path */
                PNode* b = node_new(false);
                b->w[0].c = child;
                child = b;
            }
            nr->w[1].c = child;
            v->root = nr;
            v->shift += PV_BITS;
        } else {
            v->root = push_tail(v->root, v->shift, v->count - PV_WIDTH, full, v->kind);
        }
        v->tail = nt;
    }
    v->count++;
    v->id = g_next_id++;
    return v;
}

PVec* pv_push_i64(PVec* v, int64_t x)  { return pv_push(v, expr_new_integer(x)); }
PVec* pv_push_real(PVec* v, double x)  { return pv_push(v, expr_new_real(x)); }

/* ---- pop ----------------------------------------------------------------- */
/* Remove the rightmost leaf under n at `level`; returns the new subtree (NULL
 * when it becomes empty).  Consumes n. */
static PNode* pop_tail(PNode* n, int level, int64_t last, PVKind kind) {
    int j = (int)((last >> level) & (PV_WIDTH - 1));
    if (level > PV_BITS) {
        PNode* c = pop_tail(n->w[j].c ? node_retain(n->w[j].c) : NULL, level - PV_BITS, last, kind);
        if (!c && j == 0) { node_release(n, kind); return NULL; }
        n = node_own(n, kind);
        node_release(n->w[j].c, kind);
        n->w[j].c = c;
        return n;
    }
    if (j == 0) { node_release(n, kind); return NULL; }
    n = node_own(n, kind);
    node_release(n->w[j].c, kind);
    n->w[j].c = NULL;
    return n;
}

PVec* pv_pop(PVec* v, Expr** out) {
    if (out) *out = NULL;
    if (v->count == 0) return v;
    v = pv_own(v);
    int64_t toff = tail_off(v);
    int tl = (int)(v->count - toff);
    if (out) *out = word_to_expr(v->tail->w[tl - 1], v->kind);

    if (v->count == 1) {
        node_release(v->tail, v->kind); v->tail = NULL;
        node_release(v->root, v->kind); v->root = NULL;
        v->count = 0; v->shift = PV_BITS; v->kind = PV_EMPTY;
    } else if (tl > 1) {
        v->tail = node_own(v->tail, v->kind);
        if (v->kind == PV_BOXED) expr_free(v->tail->w[tl - 1].e);
        v->tail->w[tl - 1].i = 0;
        v->tail->n = (uint8_t)(tl - 1);
        v->count--;
    } else {
        /* the tail holds one element: the last leaf of the trie becomes the tail */
        int64_t last = v->count - 2;                 /* index inside the previous leaf */
        PNode* nt = node_retain((PNode*)leaf_for(v, last));
        PNode* r = pop_tail(v->root, v->shift, v->count - 2, v->kind);
        v->root = NULL;
        node_release(v->tail, v->kind);
        v->tail = nt;
        if (r && v->shift > PV_BITS && !r->w[1].c) {
            PNode* only = node_retain(r->w[0].c);
            node_release(r, v->kind);
            r = only;
            v->shift -= PV_BITS;
        }
        v->root = r;
        v->count--;
    }
    v->id = g_next_id++;
    return v;
}

/* ---- conversion ---------------------------------------------------------- */
PVec* pv_from_elems(Expr* const* elems, size_t n) {
    PVec* v = pv_new();
    for (size_t i = 0; i < n; i++) v = pv_push(v, expr_copy(elems[i]));
    return v;
}

Expr* pv_to_list(const PVec* v) {
    size_t n = (size_t)v->count;
    Expr** items = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    size_t k = 0;
    int64_t i = 0;
    while (i < v->count) {
        const void* p;
        int64_t span = pv_span_at(v, i, &p);
        const PWord* w = (const PWord*)p;
        for (int64_t j = 0; j < span; j++) items[k++] = word_to_expr(w[j], v->kind);
        i += span;
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return r;
}

/* ---- equality, hash, order ------------------------------------------------ */
bool pv_equal(const PVec* a, const PVec* b) {
    if (a == b || pv_same(a, b)) return true;
    if (a->count != b->count) return false;
    if (a->kind == b->kind && a->kind != PV_BOXED) {
        for (int64_t i = 0; i < a->count; i++) {
            const PWord wa = leaf_for(a, i)->w[i & (PV_WIDTH - 1)];
            const PWord wb = leaf_for(b, i)->w[i & (PV_WIDTH - 1)];
            if (a->kind == PV_I64 ? wa.i != wb.i : wa.d != wb.d) return false;
        }
        return true;
    }
    for (int64_t i = 0; i < a->count; i++) {
        Expr* x = pv_get(a, i);
        Expr* y = pv_get(b, i);
        bool eq = expr_eq(x, y);
        expr_free(x); expr_free(y);
        if (!eq) return false;
    }
    return true;
}

/* Equal vectors hash equal whatever their storage kind: a sample of at most 16
 * elements (spread over the vector) hashed as Exprs, mixed with the length. */
uint64_t pv_hash(const PVec* v) {
    uint64_t h = 0x9E3779B97F4A7C15ULL ^ (uint64_t)v->count;
    int64_t n = v->count;
    int64_t step = n > 16 ? n / 16 : 1;
    for (int64_t i = 0; i < n; i += step) {
        Expr* x = pv_get(v, i);
        h = (h ^ expr_hash(x)) * 1099511628211ULL;
        expr_free(x);
    }
    return h;
}

int pv_compare(const PVec* a, const PVec* b) {
    if (a->count != b->count) return a->count < b->count ? -1 : 1;
    for (int64_t i = 0; i < a->count; i++) {
        Expr* x = pv_get(a, i);
        Expr* y = pv_get(b, i);
        int c = expr_compare(x, y);
        expr_free(x); expr_free(y);
        if (c) return c < 0 ? -1 : 1;
    }
    return 0;
}
