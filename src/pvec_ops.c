/* pvec_ops.c -- PVector at the expression level.  See pvec_ops.h. */
#include "pvec_ops.h"
#include "pvec.h"
#include "eval.h"
#include "attr.h"
#include "symtab.h"
#include "sym_names.h"
#include "ndarray.h"
#include <stdlib.h>
#include <string.h>

/* ---- construction --------------------------------------------------------- */
static bool is_plain_list(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

Expr* pvec_make(const Expr* e) {
    if (is_plain_list(e)) {
        PVec* v = pv_new();
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            v = pv_push(v, expr_copy(e->data.function.args[i]));
        return expr_new_pvec(v);
    }
    if (is_ndarray(e)) {
        /* A rank-1 int64 / float64 buffer copies straight into typed leaves. */
        if (e->data.ndarray.rank == 1 && e->data.ndarray.dtype == NDT_INT64) {
            const int64_t* d = (const int64_t*)e->data.ndarray.data;
            PVec* v = pv_new();
            for (int64_t i = 0; i < e->data.ndarray.dims[0]; i++) v = pv_push_i64(v, d[i]);
            return expr_new_pvec(v);
        }
        if (e->data.ndarray.rank == 1 && e->data.ndarray.dtype == NDT_FLOAT64) {
            const double* d = (const double*)e->data.ndarray.data;
            PVec* v = pv_new();
            for (int64_t i = 0; i < e->data.ndarray.dims[0]; i++) v = pv_push_real(v, d[i]);
            return expr_new_pvec(v);
        }
        if (is_packed_list(e)) {                         /* any other packed List */
            Expr* l = ndarray_to_nested_list(e);
            Expr* r = l ? pvec_make(l) : NULL;
            if (l) expr_free(l);
            return r;
        }
    }
    return NULL;
}

static Expr* builtin_pvector(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return expr_new_pvec(pv_new());
    if (n != 1) return NULL;
    Expr* a = res->data.function.args[0];
    if (a->type == EXPR_PVEC) return expr_copy(a);        /* immutable: itself */
    return pvec_make(a);                                  /* NULL leaves it unevaluated */
}

static Expr* builtin_pvectorq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return expr_new_symbol(res->data.function.args[0]->type == EXPR_PVEC ? SYM_True : SYM_False);
}

/* ---- the gate ------------------------------------------------------------- */
typedef enum { PG_PLAIN, PG_WRAP } PGAct;
static const struct { const char* name; PGAct act; } g_heads[] = {
    /* structure-preserving: the answer is a PVector */
    {"Map", PG_WRAP}, {"Select", PG_WRAP}, {"Take", PG_WRAP}, {"Drop", PG_WRAP},
    {"Rest", PG_WRAP}, {"Most", PG_WRAP}, {"Join", PG_WRAP}, {"Reverse", PG_WRAP},
    {"Sort", PG_WRAP}, {"Delete", PG_WRAP}, {"Insert", PG_WRAP}, {"ReplacePart", PG_WRAP},
    {"Riffle", PG_WRAP}, {"DeleteCases", PG_WRAP}, {"Accumulate", PG_WRAP},
    {"FoldList", PG_WRAP}, {"Differences", PG_WRAP}, {"Union", PG_WRAP},
    {"DeleteDuplicates", PG_WRAP}, {"RotateLeft", PG_WRAP}, {"RotateRight", PG_WRAP},
    {"PadLeft", PG_WRAP}, {"PadRight", PG_WRAP}, {"Prepend", PG_WRAP}, {"Append", PG_WRAP},
    /* reductions and queries: answer exactly as for the List */
    {"Total", PG_PLAIN}, {"Mean", PG_PLAIN}, {"Median", PG_PLAIN}, {"Max", PG_PLAIN},
    {"Min", PG_PLAIN}, {"Apply", PG_PLAIN}, {"Fold", PG_PLAIN}, {"Scan", PG_PLAIN},
    {"Count", PG_PLAIN}, {"Position", PG_PLAIN}, {"MemberQ", PG_PLAIN}, {"FreeQ", PG_PLAIN},
    {"Dot", PG_PLAIN}, {"Norm", PG_PLAIN}, {"Variance", PG_PLAIN},
    {"StandardDeviation", PG_PLAIN}, {"Tally", PG_PLAIN}, {"Counts", PG_PLAIN},
    {"Flatten", PG_PLAIN}, {"Ordering", PG_PLAIN}, {"Cases", PG_PLAIN}, {"First", PG_PLAIN},
    {"Last", PG_PLAIN}
};

static bool any_pvec_arg(const Expr* res) {
    for (size_t i = 0; i < res->data.function.arg_count; i++)
        if (res->data.function.args[i]->type == EXPR_PVEC) return true;
    return false;
}

/* The call with every PVector argument replaced by its List. */
static Expr* delisted_call(const Expr* res) {
    size_t n = res->data.function.arg_count;
    Expr** a = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        const Expr* x = res->data.function.args[i];
        a[i] = x->type == EXPR_PVEC ? pv_to_list(x->data.pvec) : expr_copy((Expr*)x);
    }
    Expr* c = expr_new_function(expr_copy(res->data.function.head), a, n);
    free(a);
    return c;
}

/* Wrap a List answer back into a PVector; anything else is returned as is. */
static Expr* wrap_result(Expr* r) {
    if (!r) return NULL;
    Expr* w = pvec_make(r);
    if (!w) return r;
    expr_free(r);
    return w;
}

Expr* pvec_gate(Expr* res, const char* head, uint32_t attrs) {
    if (res->type != EXPR_FUNCTION || !any_pvec_arg(res)) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr* a0 = argc ? res->data.function.args[0] : NULL;

    /* Native operations: no materialisation. */
    if (a0 && a0->type == EXPR_PVEC) {
        PVec* pv = a0->data.pvec;
        int64_t len = pv_length(pv);
        if (head == SYM_Append && argc == 2)
            return expr_new_pvec(pv_push(pv_retain(pv), expr_copy(res->data.function.args[1])));
        if (head == SYM_Prepend && argc == 2) {
            PVec* v = pv_push(pv_new(), expr_copy(res->data.function.args[1]));
            for (int64_t i = 0; i < len; i++) v = pv_push(v, pv_get(pv, i));
            return expr_new_pvec(v);
        }
        if (argc == 1 && len > 0) {
            if (head == SYM_First) return pv_get(pv, 0);
            if (head == SYM_Last)  return pv_get(pv, len - 1);
            if (head == SYM_Most)  return expr_new_pvec(pv_pop(pv_retain(pv), NULL));
            if (head == SYM_Rest) {
                PVec* v = pv_new();
                for (int64_t i = 1; i < len; i++) v = pv_push(v, pv_get(pv, i));
                return expr_new_pvec(v);
            }
        }
    }
    /* Append / Prepend with the vector not first: the vector is just an element. */
    if (head == SYM_Append || head == SYM_Prepend) return NULL;

    PGAct act = PG_PLAIN;
    bool listed = false;
    for (size_t i = 0; i < sizeof g_heads / sizeof *g_heads; i++)
        if (strcmp(head, g_heads[i].name) == 0) { act = g_heads[i].act; listed = true; break; }
    if (!listed) {
        if (!(attrs & ATTR_LISTABLE)) return NULL;        /* not ours: the vector is data */
        act = PG_WRAP;                                    /* listable arithmetic: pv + 1 */
    }
    Expr* r = eval_and_free(delisted_call(res));
    return act == PG_WRAP ? wrap_result(r) : r;
}

/* ---- registration ---------------------------------------------------------- */
void pvec_init(void) {
    symtab_add_builtin(SYM_PVector, builtin_pvector);
    symtab_get_def(SYM_PVector)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PVector",
        "PVector[]\n\tAn empty persistent vector. PVector[list] holds the elements of a List\n"
        "\tor packed array. A persistent vector is an immutable VALUE: Append, Prepend,\n"
        "\tMost, Rest, Part assignment v[[i]] = x and the rest return a new vector that\n"
        "\tshares all untouched storage with the old one, so keeping the old value is\n"
        "\tfree. v[[i]] and Append[v, x] are O(log32 n) (effectively constant), as\n"
        "\tis v[[i]] = x; AppendTo[v, x] and Part assignment on a vector no one else\n"
        "\tholds update in place. Length, Dimensions, Part, First, Last, Most, Rest,\n"
        "\tNormal and pattern _PVector work on it; Map, Select, Take, Drop, Join,\n"
        "\tReverse, Sort and Listable functions (v + 1) return a PVector, and Total,\n"
        "\tMax, Count and the like answer as for the List. Storage is unboxed while\n"
        "\tevery element is a machine integer or real. Prints as PVector[<n>].");
    symtab_add_builtin(SYM_PVectorQ, builtin_pvectorq);
    symtab_get_def(SYM_PVectorQ)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PVectorQ", "PVectorQ[e]\n\tGives True if e is a persistent vector, and False otherwise.");
}
