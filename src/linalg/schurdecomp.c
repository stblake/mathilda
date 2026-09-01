/*
 * schurdecomp.c -- SchurDecomposition dispatcher.
 *
 * Parses options, distinguishes the standard form SchurDecomposition[m] from
 * the generalized form SchurDecomposition[{m, a}], classifies the numeric
 * precision of the input, and routes to the right kernel:
 *
 *   - non-numeric (symbolic) matrix        -> NULL (call left unevaluated;
 *                                              a generic matrix has no
 *                                              closed-form Schur decomposition)
 *   - standard, real, arbitrary precision  -> schur_mpfr_standard_real
 *     (RealBlockDiagonalForm -> True, no Pivoting; also the no-LAPACK fallback)
 *   - everything else numeric              -> schur_machine_standard /
 *                                              schur_machine_generalized (LAPACK)
 *
 * An NDArray / packed matrix is handled transparently: the machine kernel loads
 * it via numarray.c's na_load_matrix (which accepts both an NDArray and a boxed
 * List-of-Lists), and schur_matrix_order reads its rank-2 shape directly.  So no
 * separate ndla_* guard is needed -- but the head IS on src/pack.c's AWARE list,
 * so the packing gate hands a packed argument straight through as an NDArray
 * rather than materialising it (see docs/design/packed_arrays.md).
 *
 * Memory contract: standard builtin ownership (SPEC.md §4).  Never frees `res`.
 */

#include "schurdecomp.h"
#include "schurdecomp_internal.h"
#include "linalg.h"
#include "common.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "print.h"
#include "pack.h"
#include "ndarray.h"
#include "numarray.h"
#include "lapack.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------------ *
 *  Option parsing (mirrors qrdecomp.c: parse_bool_value / parse_target...). *
 * ------------------------------------------------------------------------ */

static bool schur_parse_bool(const Expr* rhs, bool* out) {
    if (rhs->type != EXPR_SYMBOL) return false;
    const char* s = rhs->data.symbol.name;
    if (strcmp(s, "True")  == 0) { *out = true;  return true; }
    if (strcmp(s, "False") == 0) { *out = false; return true; }
    return false;
}

static bool schur_parse_targetstructure(const Expr* rhs, bool* structured) {
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) {
        *structured = false;
        return true;
    }
    if (rhs->type != EXPR_STRING) return false;
    if (strcmp(rhs->data.string, "Dense")      == 0) { *structured = false; return true; }
    if (strcmp(rhs->data.string, "Structured") == 0) { *structured = true;  return true; }
    return false;
}

bool schur_parse_options(const Expr* res, SchurOpts* opts) {
    opts->pivoting = false;
    opts->real_block_diagonal_form = true;
    opts->structured = false;

    size_t argc = res->data.function.arg_count;
    for (size_t i = 1; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type != EXPR_FUNCTION
            || opt->data.function.head->type != EXPR_SYMBOL
            || opt->data.function.arg_count != 2) return false;
        const char* hd = opt->data.function.head->data.symbol.name;
        if (hd != SYM_Rule && hd != SYM_RuleDelayed) return false;
        Expr* lhs = opt->data.function.args[0];
        Expr* rhs = opt->data.function.args[1];
        if (lhs->type != EXPR_SYMBOL) return false;
        const char* ln = lhs->data.symbol.name;

        if (ln == SYM_Pivoting) {
            if (!schur_parse_bool(rhs, &opts->pivoting)) return false;
        } else if (ln == SYM_RealBlockDiagonalForm) {
            if (!schur_parse_bool(rhs, &opts->real_block_diagonal_form)) return false;
        } else if (ln == SYM_TargetStructure) {
            if (!schur_parse_targetstructure(rhs, &opts->structured)) return false;
        } else {
            return false;  /* unknown option -> leave the call unevaluated */
        }
    }
    return true;
}

/* ------------------------------------------------------------------------ *
 *  Shape classification.                                                    *
 * ------------------------------------------------------------------------ */

int schur_matrix_order(const Expr* e) {
    if (is_ndarray(e)) {
        const NDArrayData* nd = &e->data.ndarray;
        if (nd->rank == 2 && nd->dims[0] > 0 && nd->dims[0] == nd->dims[1])
            return (int)nd->dims[0];
        return -1;
    }
    int64_t dims[64];
    int rank = get_tensor_dims((Expr*)e, dims);
    if (rank == 2 && dims[0] > 0 && dims[1] > 0 && dims[0] == dims[1])
        return (int)dims[0];
    return -1;
}

/* {m, a} generalized pair: a length-2 List whose FIRST element is itself a
 * square matrix (which distinguishes it from a 2xN ordinary matrix, whose first
 * element is a rank-1 row).  Both matrices must be square and the same order. */
static bool schur_is_pair(const Expr* arg0, const Expr** out_m,
                          const Expr** out_a, int* out_n) {
    if (arg0->type != EXPR_FUNCTION) return false;
    Expr* head = arg0->data.function.head;
    if (head->type != EXPR_SYMBOL || head->data.symbol.name != SYM_List) return false;
    if (arg0->data.function.arg_count != 2) return false;
    const Expr* m = arg0->data.function.args[0];
    const Expr* a = arg0->data.function.args[1];
    int nm = schur_matrix_order(m);
    int na = schur_matrix_order(a);
    if (nm < 1 || nm != na) return false;
    *out_m = m;
    *out_a = a;
    *out_n = nm;
    return true;
}

/* True when every leaf of the matrix is a real number (no Complex[...] leaf).
 * A non-numeric leaf reports "not complex" here; the numeric load downstream is
 * what ultimately rejects it (returning NULL), so the final behaviour on a
 * symbolic matrix is an unevaluated call either way. */
MATHILDA_MAYBE_UNUSED static bool schur_all_real(const Expr* m, int n) {
    if (is_ndarray(m)) return !ndt_is_complex(m->data.ndarray.dtype);
    if (m->type != EXPR_FUNCTION) return false;
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION) return false;
        int cols = (int)row->data.function.arg_count;
        for (int j = 0; j < cols; j++) {
            double re = 0.0, im = 0.0;
            if (na_read_scalar(row->data.function.args[j], &re, &im) && im != 0.0)
                return false;
        }
    }
    return true;
}

/* ------------------------------------------------------------------------ *
 *  Builtin entry point.                                                     *
 * ------------------------------------------------------------------------ */

Expr* builtin_schurdecomposition(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t argc = res->data.function.arg_count;
    if (argc < 1) {
        fprintf(stderr,
                "SchurDecomposition::argx: SchurDecomposition called with 0 "
                "arguments; 1 argument is expected.\n");
        return NULL;
    }

    SchurOpts opts;
    if (!schur_parse_options(res, &opts)) return NULL;

    const Expr* arg0 = res->data.function.args[0];
    const Expr* m = NULL;
    const Expr* a = NULL;
    int n = 0;
    bool generalized = schur_is_pair(arg0, &m, &a, &n);

    if (!generalized) {
        n = schur_matrix_order(arg0);
        if (n < 1) {
            char* s = expr_to_string_fullform((Expr*)arg0);
            fprintf(stderr,
                    "SchurDecomposition::sqma: Argument %s at position 1 is not "
                    "a non-empty square matrix or a pair of square matrices.\n",
                    s ? s : "?");
            free(s);
            return NULL;
        }
        m = arg0;
    }

    if (generalized) {
        /* Generalized (QZ) Schur -- machine LAPACK only; an arbitrary-precision
         * pair falls back to machine precision (MPFR QZ is a documented gap). */
        return schur_machine_generalized(m, a, n, &opts);
    }

    /* Standard form.  Classify precision: an NDArray is machine by dtype; a
     * boxed matrix is machine iff every inexact leaf is <= 53 bits. */
    bool machine = true;
    long bits = 53;
    if (!is_ndarray(m)) {
        CommonInexactInfo info = common_scan_inexact(m);
        if (info.has_inexact && info.min_bits > 53) {
            machine = false;
            bits = info.min_bits;
        }
    }

#ifdef USE_MPFR
    /* The MPFR real-Schur kernel handles the standard, real, block-diagonal,
     * un-pivoted case.  Use it when the input is arbitrary precision (to keep
     * the precision) OR when LAPACK is unavailable (it is the only path). */
    if (opts.real_block_diagonal_form && !opts.pivoting && schur_all_real(m, n)
        && (!machine || !mathilda_lapack_probe())) {
        mpfr_prec_t pb = (mpfr_prec_t)(bits > 53 ? bits : 53);
        Expr* r = schur_mpfr_standard_real(m, n, pb, &opts);
        if (r) return r;
        /* fall through to the machine kernel on MPFR non-convergence */
    }
#endif

    return schur_machine_standard(m, n, &opts);
}

void schurdecomp_init(void) {
    symtab_add_builtin("SchurDecomposition", builtin_schurdecomposition);
    symtab_get_def("SchurDecomposition")->attributes |= ATTR_PROTECTED;
}
