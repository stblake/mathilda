/* Mathilda -- sparse linear algebra: exact rational elimination.
 *
 * LinearSolve, Det and MatrixRank of an exact (integer / rational) SparseArray
 * without densifying: a right-looking LU over FLINT's fmpq with a Markowitz-
 * style pivot choice -- the shortest active row (kept in count buckets), and
 * in it the entry whose column has the fewest active entries. Fill-free
 * systems (banded, block-diagonal, ...) then cost O(nnz) rational operations;
 * the dense fraction-free path (FLINT fmpq_mat) costs O(n^3) whatever the
 * sparsity, and is refused outright past the 2^27-element densify cap.
 *
 * LinearSolve eliminates the right-hand sides alongside and back-substitutes
 * in pivot order; a singular system takes the free (non-pivot) unknowns as 0
 * and is consistent iff every unpivoted row's right-hand side is 0 -- the
 * reference's particular-solution semantics for exact input (LinearSolve::nosol
 * otherwise). Det is the signed product of the pivots, MatrixRank their count.
 *
 * Used for exact sparse matrices of order >= 64 that are sparse enough
 * (density < 10%) or that do not fit the densify cap at all; anything else
 * (symbolic entries, a nonzero background, small or dense matrices) keeps the
 * dense exact path.
 */
#include "sparse_la.h"
#include "sparse_internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "eval.h"
#include "message.h"
#include <stdlib.h>
#include <string.h>

#ifdef USE_FLINT
#include <flint/fmpq.h>
#include "flint_bridge_internal.h"

#define EXACT_MIN_N 64

typedef struct { int64_t col; fmpq v; } Ent;
typedef struct { Ent* e; int64_t n, cap; } Row;

typedef struct {
    int64_t n;           /* rows = cols */
    Row* rows;
    int64_t* colcnt;     /* active entries per column (approximate upper bound) */
    /* row-count buckets of the active rows: doubly linked lists */
    int64_t* head;       /* head[c]: first active row with c entries, -1 none */
    int64_t *next, *prev;
    int64_t* cnt;        /* current bucket of each row (-1: retired) */
} Elim;

static void row_free(Row* r) {
    for (int64_t k = 0; k < r->n; k++) fmpq_clear(&r->e[k].v);
    free(r->e);
    r->e = NULL;
    r->n = r->cap = 0;
}

static void bucket_remove(Elim* E, int64_t i) {
    int64_t c = E->cnt[i];
    if (c < 0) return;
    if (E->prev[i] >= 0) E->next[E->prev[i]] = E->next[i];
    else E->head[c] = E->next[i];
    if (E->next[i] >= 0) E->prev[E->next[i]] = E->prev[i];
    E->cnt[i] = -1;
}
static void bucket_insert(Elim* E, int64_t i) {
    int64_t c = E->rows[i].n;
    if (c > E->n) c = E->n;
    E->cnt[i] = c;
    E->prev[i] = -1;
    E->next[i] = E->head[c];
    if (E->head[c] >= 0) E->prev[E->head[c]] = i;
    E->head[c] = i;
}

static void elim_free(Elim* E) {
    if (E->rows) for (int64_t i = 0; i < E->n; i++) row_free(&E->rows[i]);
    free(E->rows); free(E->colcnt); free(E->head); free(E->next); free(E->prev); free(E->cnt);
    memset(E, 0, sizeof *E);
}

/* Load an exact rational square SparseArray (zero background). */
static bool elim_load(const SparseArray* sa, Elim* E) {
    memset(E, 0, sizeof *E);
    if (sa->rank != 2 || sa->dims[0] != sa->dims[1] || !sp_is_zero(sa->bg)) return false;
    if (sa->kind != SPV_INT64 && sa->kind != SPV_EXPR) return false;
    int64_t n = sa->dims[0];
    E->n = n;
    E->rows = calloc((size_t)n, sizeof(Row));
    E->colcnt = calloc((size_t)n, sizeof(int64_t));
    E->head = malloc(sizeof(int64_t) * (size_t)(n + 1));
    E->next = malloc(sizeof(int64_t) * (size_t)n);
    E->prev = malloc(sizeof(int64_t) * (size_t)n);
    E->cnt = malloc(sizeof(int64_t) * (size_t)n);
    if (!E->rows || !E->colcnt || !E->head || !E->next || !E->prev || !E->cnt) { elim_free(E); return false; }
    for (int64_t c = 0; c <= n; c++) E->head[c] = -1;
    const SparsePattern* p = sa->pat;
    for (int64_t i = 0; i < n; i++) {
        int64_t r0 = p->rowptr[i], r1 = p->rowptr[i + 1];
        Row* R = &E->rows[i];
        R->cap = r1 - r0 > 0 ? r1 - r0 : 1;
        R->e = malloc(sizeof(Ent) * (size_t)R->cap);
        if (!R->e) { elim_free(E); return false; }
        for (int64_t k = r0; k < r1; k++) {
            Ent* t = &R->e[R->n];
            fmpq_init(&t->v);
            t->col = sp_col(p, k, 0);
            R->n++;
            if (sa->kind == SPV_INT64) {
                fmpq_set_si(&t->v, ((const int64_t*)sa->vals)[k], 1);
            } else if (!fb_expr_to_fmpq(((Expr* const*)sa->vals)[k], &t->v)) {
                elim_free(E);
                return false;                  /* not a rational: symbolic */
            }
            if (fmpq_is_zero(&t->v)) {
                /* an explicit stored zero (sparse arithmetic does not
                 * re-minimise): never a pivot, never an entry */
                fmpq_clear(&t->v);
                R->n--;
                continue;
            }
            E->colcnt[t->col]++;
        }
        bucket_insert(E, i);
    }
    return true;
}

/* dst = dst - m * src over the columns not yet pivoted (both sorted). */
static bool row_axpy(Row* dst, const fmpq_t m, const Row* src, int64_t* colcnt, const char* pivoted) {
    Ent* out = malloc(sizeof(Ent) * (size_t)(dst->n + src->n + 1));
    if (!out) return false;
    int64_t a = 0, b = 0, o = 0;
    fmpq_t t;
    fmpq_init(t);
    while (a < dst->n || b < src->n) {
        int64_t ca = a < dst->n ? dst->e[a].col : INT64_MAX;
        int64_t cb = b < src->n ? src->e[b].col : INT64_MAX;
        if (cb < INT64_MAX && pivoted[cb]) { b++; continue; }
        if (ca < cb) {
            out[o++] = dst->e[a++];              /* move (no copy of the fmpq) */
        } else if (cb < ca) {
            fmpq_init(&out[o].v);
            fmpq_mul(&out[o].v, m, &src->e[b].v);
            fmpq_neg(&out[o].v, &out[o].v);
            out[o].col = cb;
            colcnt[cb]++;
            o++;
            b++;
        } else {
            fmpq_mul(t, m, &src->e[b].v);
            fmpq_sub(&dst->e[a].v, &dst->e[a].v, t);
            if (fmpq_is_zero(&dst->e[a].v)) {
                fmpq_clear(&dst->e[a].v);
                colcnt[ca]--;
            } else {
                out[o++] = dst->e[a];
            }
            a++;
            b++;
        }
    }
    fmpq_clear(t);
    free(dst->e);
    dst->e = out;
    dst->n = o;
    dst->cap = dst->n + src->n + 1;
    return true;
}

/* The elimination. piv_row[k], piv_col[k] for k < *rank; rhs (n x nrhs,
 * row-major fmpq) is updated alongside when non-NULL. Column lists are
 * recomputed lazily by scanning the active rows that contain the pivot
 * column -- found through a per-column row index rebuilt per step from the
 * candidate rows (rows are scanned only in the pivot column's support). */
typedef struct { int64_t* rows; int64_t n, cap; } IdxList;

static bool elim_run(Elim* E, fmpq* rhs, int64_t nrhs, int64_t* piv_row, int64_t* piv_col,
                     int64_t* rank, int* sign) {
    int64_t n = E->n;
    char* pivoted = calloc((size_t)n, 1);
    char* rowdone = calloc((size_t)n, 1);
    /* column -> rows containing it (may hold stale entries) */
    IdxList* col = calloc((size_t)n, sizeof(IdxList));
    if (!pivoted || !rowdone || !col) { free(pivoted); free(rowdone); free(col); return false; }
    bool ok = true;
    for (int64_t i = 0; i < n && ok; i++)
        for (int64_t k = 0; k < E->rows[i].n && ok; k++) {
            IdxList* L = &col[E->rows[i].e[k].col];
            if (L->n == L->cap) {
                int64_t nc = L->cap ? 2 * L->cap : 4;
                int64_t* nr = realloc(L->rows, sizeof(int64_t) * (size_t)nc);
                if (!nr) { ok = false; break; }
                L->rows = nr; L->cap = nc;
            }
            L->rows[L->n++] = i;
        }
    fmpq_t m;
    fmpq_init(m);
    int64_t r = 0;
    *sign = 1;
    /* permutation parity: count inversions of the pivot sequences at the end */
    while (ok) {
        /* the shortest nonempty active row */
        int64_t pr = -1;
        for (int64_t c = 1; c <= n && pr < 0; c++) if (E->head[c] >= 0) pr = E->head[c];
        if (pr < 0) break;                      /* no active entries left */
        Row* P = &E->rows[pr];
        /* in it, the entry whose column has the fewest active entries */
        int64_t best = 0;
        for (int64_t k = 1; k < P->n; k++)
            if (E->colcnt[P->e[k].col] < E->colcnt[P->e[best].col]) best = k;
        int64_t pc = P->e[best].col;
        bucket_remove(E, pr);
        rowdone[pr] = 1;
        pivoted[pc] = 1;
        piv_row[r] = pr;
        piv_col[r] = pc;
        r++;
        for (int64_t k = 0; k < P->n; k++) E->colcnt[P->e[k].col]--;
        /* eliminate column pc from every other active row that holds it */
        IdxList* L = &col[pc];
        for (int64_t t = 0; t < L->n && ok; t++) {
            int64_t i = L->rows[t];
            if (rowdone[i]) continue;
            Row* R = &E->rows[i];
            int64_t at = -1;
            for (int64_t k = 0; k < R->n; k++) if (R->e[k].col == pc) { at = k; break; }
            if (at < 0) continue;               /* stale */
            fmpq_div(m, &R->e[at].v, &P->e[best].v);
            int64_t before = R->n;
            /* remember the columns that are new to row i, for the column lists */
            ok = row_axpy(R, m, P, E->colcnt, pivoted);
            if (!ok) break;
            /* drop the eliminated pc entry (row_axpy skipped pc from P since it
             * is pivoted; R's own pc entry is now dst - m*piv = 0) */
            for (int64_t k = 0; k < R->n; k++)
                if (R->e[k].col == pc) {
                    fmpq_clear(&R->e[k].v);
                    memmove(&R->e[k], &R->e[k + 1], sizeof(Ent) * (size_t)(R->n - k - 1));
                    R->n--;
                    break;
                }
            /* fill: index the row's columns (a duplicate or stale index entry is
             * harmless: a visit re-checks that the row still holds the column) */
            if (R->n >= before)
                for (int64_t k = 0; k < R->n && ok; k++) {
                    IdxList* C = &col[R->e[k].col];
                    bool have = false;
                    for (int64_t q = C->n; q-- > 0 && q >= C->n - 8;) if (C->rows[q] == i) { have = true; break; }
                    if (have) continue;
                    if (C->n == C->cap) {
                        int64_t nc = C->cap ? 2 * C->cap : 4;
                        int64_t* nr = realloc(C->rows, sizeof(int64_t) * (size_t)nc);
                        if (!nr) { ok = false; break; }
                        C->rows = nr; C->cap = nc;
                    }
                    C->rows[C->n++] = i;
                }
            if (rhs)
                for (int64_t q = 0; q < nrhs; q++) {
                    fmpq_t t2;
                    fmpq_init(t2);
                    fmpq_mul(t2, m, rhs + pr * nrhs + q);
                    fmpq_sub(rhs + i * nrhs + q, rhs + i * nrhs + q, t2);
                    fmpq_clear(t2);
                }
            bucket_remove(E, i);
            bucket_insert(E, i);
        }
        /* the pivot's own pc entry stays in P for back-substitution; P is
         * frozen from here on */
    }
    fmpq_clear(m);
    for (int64_t j = 0; j < n; j++) free(col[j].rows);
    free(col);
    free(pivoted);
    free(rowdone);
    *rank = r;
    if (ok && r == n) {
        /* sign of the row and column permutations (selection order) */
        for (int pass = 0; pass < 2; pass++) {
            const int64_t* perm = pass ? piv_col : piv_row;
            /* parity by cycle decomposition */
            char* seen = calloc((size_t)n, 1);
            if (!seen) return false;
            for (int64_t s = 0; s < n; s++) {
                if (seen[s]) continue;
                int64_t len = 0;
                for (int64_t t = s; !seen[t]; t = perm[t]) { seen[t] = 1; len++; }
                if (len % 2 == 0) *sign = -*sign;
            }
            free(seen);
        }
    }
    return ok;
}

/* The pivot value of row pr at column pc (it is still in the frozen row). */
static const fmpq* pivot_val(const Elim* E, int64_t pr, int64_t pc) {
    const Row* R = &E->rows[pr];
    for (int64_t k = 0; k < R->n; k++) if (R->e[k].col == pc) return &R->e[k].v;
    return NULL;
}

/* Should this exact SparseArray take the sparse elimination? */
static bool want_exact_sparse(const SparseArray* sa) {
    if (sa->rank != 2 || sa->dims[0] != sa->dims[1] || sa->dims[0] < EXACT_MIN_N) return false;
    if (sa->kind != SPV_INT64 && sa->kind != SPV_EXPR) return false;
    if (!sp_is_zero(sa->bg)) return false;
    int64_t n = sa->dims[0];
    bool fits = n <= ((int64_t)1 << 27) / n;
    double density = (double)sa->pat->nnz / ((double)n * (double)n);
    return !fits || density < 0.1;
}

Expr* spla_exact_det(const Expr* A) {
    const SparseArray* sa = sp_of(A);
    if (!sa || !want_exact_sparse(sa)) return NULL;
    Elim E;
    if (!elim_load(sa, &E)) return NULL;
    int64_t n = E.n, rank = 0;
    int sign = 1;
    int64_t* pr = malloc(sizeof(int64_t) * (size_t)n);
    int64_t* pc = malloc(sizeof(int64_t) * (size_t)n);
    Expr* out = NULL;
    if (pr && pc && elim_run(&E, NULL, 0, pr, pc, &rank, &sign)) {
        fmpq_t d;
        fmpq_init(d);
        if (rank == n) {
            fmpq_set_si(d, sign, 1);
            for (int64_t k = 0; k < n; k++) fmpq_mul(d, d, pivot_val(&E, pr[k], pc[k]));
        }
        out = fb_fmpq_to_expr(d);
        fmpq_clear(d);
    }
    free(pr); free(pc);
    elim_free(&E);
    return out;
}

Expr* spla_exact_rank(const Expr* A) {
    const SparseArray* sa = sp_of(A);
    if (!sa || !want_exact_sparse(sa)) return NULL;
    Elim E;
    if (!elim_load(sa, &E)) return NULL;
    int64_t n = E.n, rank = 0;
    int sign = 1;
    int64_t* pr = malloc(sizeof(int64_t) * (size_t)n);
    int64_t* pc = malloc(sizeof(int64_t) * (size_t)n);
    Expr* out = NULL;
    if (pr && pc && elim_run(&E, NULL, 0, pr, pc, &rank, &sign)) out = expr_new_integer(rank);
    free(pr); free(pc);
    elim_free(&E);
    return out;
}

/* LinearSolve[A, b] for exact A (sparse) and exact b (vector or matrix);
 * NULL to decline; expr_copy(call) after ::nosol. */
Expr* spla_exact_solve(const Expr* call, const Expr* A, const Expr* B) {
    const SparseArray* sa = sp_of(A);
    if (!sa || !want_exact_sparse(sa)) return NULL;
    int64_t n = sa->dims[0];
    /* the right-hand side as rows of fmpq */
    Expr* bd = NULL;
    const Expr* b = B;
    if (B->type == EXPR_SPARSE) { bd = sp_normal(sp_of(B)); if (!bd) return NULL; b = bd; }
    if (b->type != EXPR_FUNCTION || b->data.function.head->type != EXPR_SYMBOL ||
        b->data.function.head->data.symbol.name != SYM_List || (int64_t)b->data.function.arg_count != n) {
        expr_free(bd);
        return NULL;
    }
    bool matrix = false;
    int64_t nrhs = 1;
    const Expr* b0 = b->data.function.args[0];
    if (b0->type == EXPR_FUNCTION && b0->data.function.head->type == EXPR_SYMBOL &&
        b0->data.function.head->data.symbol.name == SYM_List) {
        matrix = true;
        nrhs = (int64_t)b0->data.function.arg_count;
        if (nrhs == 0) { expr_free(bd); return NULL; }
    }
    fmpq* rhs = malloc(sizeof(fmpq) * (size_t)(n * nrhs));
    if (!rhs) { expr_free(bd); return NULL; }
    for (int64_t k = 0; k < n * nrhs; k++) fmpq_init(rhs + k);
    bool ok = true;
    for (int64_t i = 0; i < n && ok; i++) {
        const Expr* bi = b->data.function.args[i];
        if (!matrix) { ok = fb_expr_to_fmpq(bi, rhs + i) != 0; continue; }
        if (bi->type != EXPR_FUNCTION || bi->data.function.head->type != EXPR_SYMBOL ||
            bi->data.function.head->data.symbol.name != SYM_List || (int64_t)bi->data.function.arg_count != nrhs) {
            ok = false;
            break;
        }
        for (int64_t q = 0; q < nrhs && ok; q++) ok = fb_expr_to_fmpq(bi->data.function.args[q], rhs + i * nrhs + q) != 0;
    }
    expr_free(bd);
    Elim E;
    if (!ok || !elim_load(sa, &E)) {
        for (int64_t k = 0; k < n * nrhs; k++) fmpq_clear(rhs + k);
        free(rhs);
        return NULL;
    }
    int64_t rank = 0;
    int sign = 1;
    int64_t* pr = malloc(sizeof(int64_t) * (size_t)n);
    int64_t* pc = malloc(sizeof(int64_t) * (size_t)n);
    Expr* out = NULL;
    if (pr && pc && elim_run(&E, rhs, nrhs, pr, pc, &rank, &sign)) {
        /* consistency of the unpivoted rows */
        char* used = calloc((size_t)n, 1);
        bool consistent = used != NULL;
        if (used) {
            for (int64_t k = 0; k < rank; k++) used[pr[k]] = 1;
            for (int64_t i = 0; i < n && consistent; i++)
                if (!used[i])
                    for (int64_t q = 0; q < nrhs; q++)
                        if (!fmpq_is_zero(rhs + i * nrhs + q)) { consistent = false; break; }
        }
        free(used);
        if (!consistent) {
            mth_message("LinearSolve", "nosol", "Linear equation encountered that has no solution.");
            out = expr_copy((Expr*)call);
        } else {
            fmpq* x = malloc(sizeof(fmpq) * (size_t)(n * nrhs));
            if (x) {
                for (int64_t k = 0; k < n * nrhs; k++) fmpq_init(x + k);    /* free unknowns: 0 */
                fmpq_t t;
                fmpq_init(t);
                for (int64_t k = rank; k-- > 0;) {
                    const Row* R = &E.rows[pr[k]];
                    const fmpq* pv = pivot_val(&E, pr[k], pc[k]);
                    for (int64_t q = 0; q < nrhs; q++) {
                        fmpq_set(t, rhs + pr[k] * nrhs + q);
                        for (int64_t e = 0; e < R->n; e++) {
                            int64_t c = R->e[e].col;
                            if (c == pc[k]) continue;
                            fmpq_submul(t, &R->e[e].v, x + c * nrhs + q);
                        }
                        fmpq_div(x + pc[k] * nrhs + q, t, pv);
                    }
                }
                fmpq_clear(t);
                Expr** items = malloc(sizeof(Expr*) * (size_t)n);
                if (items) {
                    for (int64_t i = 0; i < n; i++) {
                        if (!matrix) { items[i] = fb_fmpq_to_expr(x + i); continue; }
                        Expr** row = malloc(sizeof(Expr*) * (size_t)nrhs);
                        for (int64_t q = 0; q < nrhs; q++) row[q] = fb_fmpq_to_expr(x + i * nrhs + q);
                        items[i] = expr_new_function(expr_new_symbol(SYM_List), row, (size_t)nrhs);
                        free(row);
                    }
                    out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
                    free(items);
                }
                for (int64_t k = 0; k < n * nrhs; k++) fmpq_clear(x + k);
                free(x);
            }
        }
    }
    free(pr); free(pc);
    elim_free(&E);
    for (int64_t k = 0; k < n * nrhs; k++) fmpq_clear(rhs + k);
    free(rhs);
    return out;
}

#else  /* no FLINT: the dense exact path answers */

Expr* spla_exact_det(const Expr* A) { (void)A; return NULL; }
Expr* spla_exact_rank(const Expr* A) { (void)A; return NULL; }
Expr* spla_exact_solve(const Expr* call, const Expr* A, const Expr* B) { (void)call; (void)A; (void)B; return NULL; }

#endif
