/*
 * sparse_gate.c -- the SparseArray transparency gate.
 *
 * A SparseArray is a raw EXPR_SPARSE node, so a builtin that reads its
 * argument's data.function.args sees an atom: before this gate Total[s],
 * Max[s] and Tr[s] all answered `s` itself, and Dot / Transpose / Det stayed
 * unevaluated. The packed surface met the same problem with its transparency
 * gate (pack.c), which MATERIALISES a buffer for every head that has not opted
 * in. That default cannot be copied here: the language treats a SparseArray as
 * an atom for the pattern heads (Count[s, 0] is 0, Cases/Position/MemberQ
 * see nothing inside), and a container -- List[s], Rule[a, s], a user's
 * f[s_] -- must keep the object. (Sort, Fold, Union, Split, Gather, Ordering
 * reject an array with ::normal in the reference; here they answer as on
 * Normal[s], a deliberate extension -- see data-structures.md.)
 *
 * So the gate is OPT-IN per head, by SymbolDef.sparse_mode, and the modes
 * record which of three behaviours the head has on an array argument
 * (tools/sparse_stress.py re-checks the invariant over every listed head):
 *
 *   SPG_DENSE     the head answers as it would on Normal[s] (Total, Mean,
 *                 Inverse, Accumulate, Equal, s . v): densify, evaluate.
 *   SPG_RESPARSE  as SPG_DENSE, but an array answer comes back sparse, over the
 *                 first sparse argument's background (Reverse, Take, Flatten,
 *                 Join, Transpose, KroneckerProduct, ...).
 *   SPG_ROWS      the head walks level 1, and level-1 parts of a SparseArray
 *                 are themselves SparseArrays (Apply, MapIndexed, Plus @@ s):
 *                 replace the array by the List of its rows.
 *   SPG_MAP       Map: rows for a matrix; for a vector the elements, and the
 *                 answer is sparse again (f applied to the background too).
 *
 * Two flags qualify a mode:
 *
 *   SPG_NATIVE    the builtin has its own sparse path (sparse_dot.c, ...). The
 *                 gate does not run before it; it runs only if the builtin
 *                 declines, as the fallback that keeps the head correct.
 *   SPG_DEEP      also look one and two List levels down (ArrayFlatten's
 *                 block matrix, MapThread's list of arrays).
 *   SPG_ANYPOS    RESPARSE when ANY argument is an array, not only the first
 *                 (Outer[f, s, v], ArrayFlatten[{{s, s}}]).
 *   SPG_ALLARR    RESPARSE only when every array argument is sparse
 *                 (Join[s, m] with a dense m is dense).
 *
 * A densification that would exceed SP_GATE_MAX elements is refused with
 * SparseArray::dense and the call is left unevaluated: an implicit 80 GB Normal
 * of a 10^5 x 10^5 matrix is never the right answer to a Total.
 *
 * If the densified call does not evaluate (the head declined it), the ORIGINAL
 * call is left at rest, so the user sees their SparseArray, not its Normal.
 */

#include "sparse_internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "eval.h"
#include "message.h"
#include "print.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Largest dense element count the gate will materialise implicitly. */
#define SP_GATE_MAX ((int64_t)1 << 27)

/* ---------------------------------------------------------------------- */
/* Densification                                                           */
/* ---------------------------------------------------------------------- */

/* The dense form of one argument for the gate, or NULL with *refused set when
 * it is too large (after SparseArray::dense). */
static Expr* gate_dense(const Expr* s, bool* refused) {
    const SparseArray* sa = sp_of(s);
    int64_t total = sp_total(sa->rank, sa->dims);
    Expr* d = (total >= 0 && total <= SP_GATE_MAX) ? sp_normal(sa) : NULL;
    if (!d) {
        char* shown = expr_to_string_outputform((Expr*)s);
        mth_message("SparseArray", "dense",
                    "%s is too large to convert to a dense array for this operation.",
                    shown ? shown : "");
        free(shown);
        *refused = true;
    }
    return d;
}

/* The List of level-1 parts of `s`: sparse rows for rank >= 2, elements for
 * rank 1. NULL on allocation failure. */
static Expr* gate_rows(const Expr* s) {
    const SparseArray* sa = sp_of(s);
    int64_t n = sa->dims[0];
    Expr** items = malloc(sizeof(Expr*) * (size_t)(n ? n : 1));
    if (!items) return NULL;
    for (int64_t i = 0; i < n; i++) {
        if (sa->rank == 1) {
            int64_t c = i;
            items[i] = sp_get(sa, &c);
        } else {
            Expr* idx = expr_new_integer(i + 1);
            bool err = false;
            items[i] = sp_part(s, &idx, 1, &err);
            expr_free(idx);
        }
        if (!items[i]) {
            for (int64_t j = 0; j < i; j++) expr_free(items[j]);
            free(items);
            return NULL;
        }
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);
    return out;
}

static bool is_list(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* Rewrite one argument (or, under SPG_DEEP, the elements of a List argument,
 * two levels down). Returns a fresh Expr when something changed, NULL when the
 * argument holds no SparseArray. *refused and *fail report a refusal or an
 * allocation failure. */
static Expr* gate_arg(const Expr* a, int mode, int depth, bool* refused, bool* fail) {
    if (a->type == EXPR_SPARSE) {
        Expr* r = (mode == SPG_ROWS) ? gate_rows(a) : gate_dense(a, refused);
        if (!r && !*refused) *fail = true;
        return r;
    }
    if (depth == 0 || !is_list(a)) return NULL;
    size_t n = a->data.function.arg_count;
    Expr** items = NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* r = gate_arg(a->data.function.args[i], mode, depth - 1, refused, fail);
        if (*refused || *fail) { expr_free(r); break; }
        if (!r) continue;
        if (!items) {
            items = calloc(n ? n : 1, sizeof(Expr*));
            if (!items) { expr_free(r); *fail = true; break; }
        }
        items[i] = r;
    }
    if (!items) return NULL;
    if (*refused || *fail) {
        for (size_t i = 0; i < n; i++) if (items[i]) expr_free(items[i]);
        free(items);
        return NULL;
    }
    for (size_t i = 0; i < n; i++)
        if (!items[i]) items[i] = expr_copy(a->data.function.args[i]);
    Expr* out = expr_new_function(expr_copy(a->data.function.head), items, n);
    free(items);
    return out;
}

/* The first SparseArray in `a`, looking `depth` List levels down. */
static const Expr* first_sparse(const Expr* a, int depth) {
    if (a->type == EXPR_SPARSE) return a;
    if (depth == 0 || !is_list(a)) return NULL;
    for (size_t i = 0; i < a->data.function.arg_count; i++) {
        const Expr* f = first_sparse(a->data.function.args[i], depth - 1);
        if (f) return f;
    }
    return NULL;
}

/* Did evaluation give back the densified call unchanged -- the head declined
 * it? (A same-head answer that DID change, such as Max[0, x] from Max over an
 * array holding x, is an answer, not a decline.) */
static bool declined(const Expr* r, const Expr* rewritten) {
    return r->type == EXPR_FUNCTION && expr_eq(r, rewritten);
}

/* The most common element of the full array `sa` (its Normal flattened to
 * the array level), as the background of a re-sparsified answer whose own
 * background never occurs. NULL if it cannot be found. */
static Expr* commonest_element(const Expr* s) {
    const SparseArray* sa = sp_of(s);
    Expr* nrm = sp_normal(sa);
    if (!nrm) return NULL;
    Expr* fa[2] = { nrm, expr_new_integer(sa->rank - 1) };
    Expr* flat = expr_new_function(expr_new_symbol("Flatten"), fa, 2);
    Expr* ca[2] = { flat, expr_new_integer(1) };
    Expr* cc = expr_new_function(expr_new_symbol("Commonest"), ca, 2);
    Expr* r = evaluate(cc);
    expr_free(cc);
    Expr* out = NULL;
    if (r->type == EXPR_FUNCTION && is_list(r) && r->data.function.arg_count >= 1)
        out = expr_copy(r->data.function.args[0]);
    expr_free(r);
    return out;
}

/* An array answer, made sparse over `bg`. Consumes `r`.
 *
 * `bg` is the leading array's background, which is the answer's background
 * whenever the head maps a background element to itself (Reverse, Take,
 * Transpose, ...). A head that changes elements may not: Clip[v - 1, {0, 10}]
 * maps the background -1 to 0. Kept anyway, every element would be explicit
 * over a background that occurs nowhere -- and the next step would carry that
 * stale background along, so two arrays with equal Normal stop being SameQ.
 * So when `bg` does not occur in the answer at all, the answer's commonest
 * element becomes the background instead. */
static Expr* resparse(Expr* r, const Expr* bg) {
    if (!(r->type == EXPR_NDARRAY || is_list(r))) return r;
    Expr* args[3] = { r, expr_new_symbol(SYM_Automatic), expr_copy((Expr*)bg) };
    Expr* call = expr_new_function(expr_new_symbol(SYM_SparseArray), args, 3);
    Expr* out = evaluate(call);
    if (out->type == EXPR_SPARSE) {
        const SparseArray* sa = sp_of(out);
        int64_t total = sp_total(sa->rank, sa->dims);
        if (total > 0 && sp_nnz(sa) == total) {
            Expr* nbg = commonest_element(out);
            if (nbg && !expr_eq(nbg, sa->bg)) {
                Expr* a2[3] = { expr_copy(call->data.function.args[0]),
                                expr_new_symbol(SYM_Automatic), nbg };
                Expr* c2 = expr_new_function(expr_new_symbol(SYM_SparseArray), a2, 3);
                Expr* again = evaluate(c2);
                expr_free(c2);
                if (again->type == EXPR_SPARSE) { expr_free(out); out = again; }
                else expr_free(again);
            } else {
                expr_free(nbg);
            }
        }
        expr_free(call);
        return out;
    }
    expr_free(call);
    /* Not a full array (ragged, or a List of non-array values): keep the
     * dense answer the densified call produced. */
    Expr* dense = expr_copy(out->data.function.arg_count ? out->data.function.args[0] : out);
    expr_free(out);
    return dense;
}

/* Public forms of the helpers above, for the native kernels (sparse.h). */
Expr* sp_gate_dense(const Expr* s, bool* refused) {
    *refused = false;
    return gate_dense(s, refused);
}

Expr* sp_gate_rows(const Expr* s, bool* refused) {
    *refused = false;
    const SparseArray* sa = sp_of(s);
    if (sa->rank == 1 && sa->dims[0] > SP_GATE_MAX) {
        Expr* none = gate_dense(s, refused);   /* the SparseArray::dense refusal */
        expr_free(none);
        return NULL;
    }
    return gate_rows(s);
}

Expr* sp_gate_resparse(Expr* r, const Expr* bg) { return resparse(r, bg); }

/* ---------------------------------------------------------------------- */
/* Diagnostics: MATHILDA_SPARSE_DIAG=gate                                  */
/* ---------------------------------------------------------------------- */

#define SP_DIAG_MAX 256
static struct { const char* name; uint64_t calls; } g_diag[SP_DIAG_MAX];
static int  g_diag_n = 0;
static int  g_diag_on = -1;

static void diag_dump(void) {
    if (!g_diag_n) {
        fprintf(stderr, "\n[sparse] gate: no SparseArray argument was ever converted.\n");
        return;
    }
    fprintf(stderr, "\n[sparse] gate converted a SparseArray argument for these heads\n"
                    "         (each one is a native sparse path not taken):\n\n");
    for (int i = 0; i < g_diag_n; i++)
        fprintf(stderr, "         %-28s %10llu\n", g_diag[i].name,
                (unsigned long long)g_diag[i].calls);
}

static void diag_note(const Expr* head) {
    if (g_diag_on < 0) {
        const char* v = getenv("MATHILDA_SPARSE_DIAG");
        g_diag_on = (v && strcmp(v, "gate") == 0);
        if (g_diag_on) atexit(diag_dump);
    }
    if (!g_diag_on || head->type != EXPR_SYMBOL) return;
    const char* nm = head->data.symbol.name;
    for (int i = 0; i < g_diag_n; i++)
        if (g_diag[i].name == nm) { g_diag[i].calls++; return; }
    if (g_diag_n < SP_DIAG_MAX) { g_diag[g_diag_n].name = nm; g_diag[g_diag_n].calls = 1; g_diag_n++; }
}

/* ---------------------------------------------------------------------- */
/* Entry point                                                             */
/* ---------------------------------------------------------------------- */

Expr* sp_gate(const Expr* call, int mode_bits, bool* rest) {
    *rest = false;
    int mode = mode_bits & SPG_MODE_MASK;
    int depth = (mode_bits & SPG_DEEP) ? 2 : 0;
    size_t n = call->data.function.arg_count;
    /* Map over a vector: densify, map, and re-sparsify; over a matrix: rows. */
    const Expr* map_vec = NULL;
    if (mode == SPG_MAP) {
        const Expr* tgt = n >= 2 ? call->data.function.args[1] : NULL;
        bool vec = tgt && tgt->type == EXPR_SPARSE && sp_of(tgt)->rank == 1;
        mode = vec ? SPG_RESPARSE : SPG_ROWS;
        if (vec) map_vec = tgt;
    }
    bool refused = false, fail = false;
    Expr** args = NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* r = gate_arg(call->data.function.args[i], mode, depth, &refused, &fail);
        if (refused || fail) { expr_free(r); break; }
        if (!r) continue;
        if (!args) {
            args = calloc(n, sizeof(Expr*));
            if (!args) { expr_free(r); fail = true; break; }
        }
        args[i] = r;
    }
    if (refused || fail || !args) {
        if (args) {
            for (size_t i = 0; i < n; i++) if (args[i]) expr_free(args[i]);
            free(args);
        }
        *rest = refused;
        return NULL;
    }
    for (size_t i = 0; i < n; i++)
        if (!args[i]) args[i] = expr_copy(call->data.function.args[i]);
    Expr* rewritten = expr_new_function(expr_copy(call->data.function.head), args, n);
    free(args);
    diag_note(call->data.function.head);

    /* Level-1 rewriting changes only the argument's spelling: the head walks
     * the List of rows on the next step exactly as it would the array. */
    if (mode == SPG_ROWS) return rewritten;

    Expr* r = evaluate(rewritten);
    bool no_answer = declined(r, rewritten);
    expr_free(rewritten);
    if (no_answer) {
        expr_free(r);
        *rest = true;
        return NULL;
    }
    /* The answer is sparse only when the array leads: f[s, v] is sparse,
     * f[v, s] and f[m, s] are dense (measured across the RESPARSE heads).
     * SPG_ANYPOS heads take their arrays after a function or inside a block
     * List (Outer[f, s, v], ArrayFlatten[{{s, s}}]): any array decides. */
    if (map_vec) {
        /* the background maps too: Map[f, s] has background f[bg] */
        Expr* fa[1] = { expr_copy(sp_of(map_vec)->bg) };
        Expr* fc = expr_new_function(expr_copy(call->data.function.args[0]), fa, 1);
        Expr* fbg = evaluate(fc);
        expr_free(fc);
        r = resparse(r, fbg);
        expr_free(fbg);
    } else if (mode == SPG_RESPARSE) {
        const Expr* lead = NULL;
        if (mode_bits & SPG_ALLARR) {
            for (size_t i = 0; i < n; i++) {
                const Expr* a = call->data.function.args[i];
                if (a->type == EXPR_SPARSE) { if (!lead) lead = a; }
                else if (is_list(a) || a->type == EXPR_NDARRAY) { lead = NULL; break; }
            }
        } else if (mode_bits & SPG_ANYPOS) {
            for (size_t i = 0; i < n && !lead; i++)
                lead = first_sparse(call->data.function.args[i], depth);
        } else if (n > 0 && call->data.function.args[0]->type == EXPR_SPARSE) {
            lead = call->data.function.args[0];
        }
        if (lead) r = resparse(r, sp_of(lead)->bg);
    }
    return r;
}

/* ---------------------------------------------------------------------- */
/* Registration                                                            */
/* ---------------------------------------------------------------------- */

/* Each list was settled call shape by call shape. A head absent from
 * every list keeps treating a SparseArray as an atom -- which is right for the
 * pattern heads (Count, Cases, Position, MemberQ, ReplaceAll, Level) and for
 * every container. */

/* Answer as on Normal[s]. */
static const char* const GATE_DENSE[] = {
    /* reductions and statistics */
    "Total", "Max", "Min", "MinMax", "Mean", "Median", "Variance",
    "StandardDeviation", "Moment", "CentralMoment", "Kurtosis", "Skewness",
    "MeanDeviation", "MedianDeviation", "InterquartileRange", "Quartiles",
    "Covariance", "Correlation", "Commonest", "Tr", "Norm",
    /* linear algebra with a dense answer */
    "Det", "Inverse", "LinearSolve", "Eigenvalues", "Eigenvectors",
    "Eigensystem", "MatrixRank", "NullSpace", "RowReduce", "LUDecomposition",
    "QRDecomposition", "SingularValueDecomposition", "SingularValueList",
    "CholeskyDecomposition", "SchurDecomposition", "JordanDecomposition",
    "HermiteDecomposition", "CharacteristicPolynomial", "LeastSquares",
    "PseudoInverse", "MatrixExp", "MatrixLog", "MatrixFunction",
    "HessenbergDecomposition", "Minors", "Adjugate", "Orthogonalize",
    /* sequence operations with a dense answer */
    "Accumulate", "Differences", "Ratios", "ListConvolve", "ListCorrelate",
    "Fourier", "InverseFourier", "FourierDCT", "FourierDST",
    "FromDigits", "Depth", "Equal", "Unequal", "Inner",
    /* Dot: s . v is dense; s . s is sparse (a native path follows) */
    "Dot",
    /* distances and predicates over the elements */
    "EuclideanDistance", "SquaredEuclideanDistance", "ManhattanDistance",
    "ChessboardDistance", "CosineDistance", "HammingDistance", "EditDistance",
    "SymmetricMatrixQ", "HermitianMatrixQ", "DiagonalMatrixQ",
    "PositiveDefiniteMatrixQ", "LowerTriangularMatrixQ", "UpperTriangularMatrixQ",
    "IntersectingQ", "DuplicateFreeQ", "SubsetQ", "OrderedQ", "SquareMatrixQ",
    /* selection, ordering and combinatorics over the level-1 parts */
    "Tally", "Pick", "Riffle", "TakeLargest", "TakeSmallest",
    "RankedMax", "RankedMin", "Quantile", "RootMeanSquare", "Standardize",
    "FlattenAt", "HankelMatrix",
    "ToeplitzMatrix", "PrincipalComponents",
    /* answered on Normal[s] (the reference refuses them with ::normal) */
    "Union", "Split", "SplitBy", "Gather", "GatherBy", "Ordering",
    "Counts", "CountsBy", "Nearest",
    /* functional heads whose answer is a List in the reference */
    "MapAt", "DeleteDuplicates", "DeleteDuplicatesBy", "AnyTrue", "AllTrue",
    "NoneTrue", "Subsets", "Complement", "Intersection", "MapAll",
    /* an array is not a polynomial or an objective: answer (or decline) as
     * on Normal[s], never treat the atom as a variable */
    "CoefficientRules", "CoefficientArrays", "CoefficientList", "MonomialList",
    "Decompose", "FactorTermsList", "Minimize", "Maximize",
    "NMinimize", "NMaximize", "FindMinimum", "FindMaximum",
    NULL
};

/* As SPG_DENSE, with an array answer re-sparsified over the leading array's
 * background. */
static const char* const GATE_RESPARSE[] = {
    "Transpose", "ConjugateTranspose", "ArrayReshape", "Reverse",
    "Take", "Drop", "Rest", "Most", "First", "Last", "Append",
    "Prepend", "Partition", "RotateLeft", "RotateRight", "ArrayPad",
    "Diagonal", "DiagonalMatrix", "KroneckerProduct", "MatrixPower",
    "Extract", "ReplacePart", "Insert", "Delete", "Chop", "Clip",
    "Expand", "ExpandAll", "Together", "Simplify", "FullSimplify", "Factor", "D", "Dt",
    "MovingAverage", "ExponentialMovingAverage", "LowerTriangularize",
    "UpperTriangularize", "Permute", "MovingMedian", "PadLeft", "PadRight",
    "Normalize",
    /* order-only rearrangements: answered on Normal[s], kept sparse */
    "Sort", "SortBy", "ReverseSort",
    /* an array of polynomials is an array of answers */
    "HornerForm",
    NULL
};

/* Walk level 1: the parts are SparseArrays themselves. */
static const char* const GATE_ROWS[] = {
    "Apply", "MapIndexed", "Scan", "Select",
    NULL
};

void sparse_gate_init(void) {
    for (const char* const* p = GATE_DENSE; *p; p++) symtab_set_sparse_mode(*p, SPG_DENSE);
    for (const char* const* p = GATE_RESPARSE; *p; p++) symtab_set_sparse_mode(*p, SPG_RESPARSE);
    for (const char* const* p = GATE_ROWS; *p; p++) symtab_set_sparse_mode(*p, SPG_ROWS);
    /* Arrays inside List arguments, two levels down. */
    symtab_set_sparse_mode("ArrayFlatten", SPG_RESPARSE | SPG_DEEP | SPG_ANYPOS);
    /* Outer is dense unless f is Times over a leading sparse array -- the
     * native kernel (sparse_func.c) decides; this is its fallback. */
    symtab_set_sparse_mode("Outer", SPG_DENSE | SPG_DEEP);
    /* Arrays nested in a List argument: Flatten[{s, s}]. (Thread reads the
     * arrays inside its f[...] itself, funcprog.c.) */
    symtab_set_sparse_mode("Flatten", SPG_RESPARSE | SPG_DEEP);
    symtab_set_sparse_mode("Join", SPG_RESPARSE | SPG_ALLARR);
    symtab_set_sparse_mode("Map", SPG_MAP);
    symtab_set_sparse_mode("MapThread", SPG_ROWS | SPG_DEEP);
}
