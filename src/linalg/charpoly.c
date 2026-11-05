/* CharacteristicPolynomial -- the characteristic polynomial of a square matrix.
 *
 *   CharacteristicPolynomial[m, x]      == Det[m - x I]
 *   CharacteristicPolynomial[{m, a}, x] == Det[m - x a]   (generalized)
 *
 * This is exactly the polynomial Eigenvalues solves for, so it reuses the eigen
 * module's char-poly machinery (declared in eigen_internal.h):
 *
 *   - Ordinary case: Faddeev-Leverrier-Souriau (O(n^4)), which returns
 *     det(lambda I - m).  Wolfram's CharacteristicPolynomial is
 *     Det[m - x I] = (-1)^n det(x I - m), so the odd-n result is negated.
 *     This keeps the 100x100 machine-matrix case sub-second: the naive
 *     Expand[Det[m - x I]] would face an O(n!) Laplace expansion of a
 *     symbolic-in-x matrix.
 *   - Generalized case: build (m - lambda a) and take its Laplace determinant,
 *     which is det(m - lambda a) with the correct sign directly.  Generalized
 *     inputs are small in practice (<= 3x3), matching the existing generalized
 *     Eigenvalues path.
 *
 * The polynomial is built in a private internal lambda symbol and the user's
 * variable is substituted in at the end, so the second argument may be a
 * symbol, a number, or any expression (the char-poly coefficient identity makes
 * det(lambda I - m)|_{lambda -> x} valid even when m itself contains x).  The
 * result is Expand-ed (Det / Laplace do not multiply out) to reach the
 * collected polynomial form.
 *
 * This head returns a symbolic Plus expression, not a machine buffer, so it is
 * exempt from the packed/NDArray/Compile surfaces.  It consumes packed input
 * safely: the NDArray guard delists a visible NDArray, faddeev pack_unpacks its
 * own input, and the generalized path materialises m / a before indexing them.
 */

#include "eigen.h"
#include "eigen_internal.h"
#include "linalg.h"
#include "ndlinalg.h"
#include "eval.h"
#include "common.h"
#include "pack.h"
#include "flint_mat_bridge.h"
#include "sym_names.h"
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

/* The all-rational ordinary case through FLINT: det(x I - m) from
 * fmpq_mat_charpoly, signed to Det[m - x I] = (-1)^n det(x I - m) and summed
 * as c_0 + c_1 x + ... + c_n x^n.  Faddeev-Leverrier over Expr arithmetic was
 * O(n^4) bignum Expr operations (1.2 s for a 40 x 40 integer matrix, ~6x
 * sympy); FLINT is multimodular.  Returns NULL (caller keeps the generic path)
 * for a matrix with any non-rational entry or without FLINT. */
static Expr* charpoly_flint(const Expr* m, int64_t n, const Expr* var) {
    if (n < 1 || n > 4096) return NULL;
    Expr* mm = pack_unpack(m);
    const Expr* rows = mm ? mm : m;
    if (rows->type != EXPR_FUNCTION || rows->data.function.arg_count != (size_t)n) {
        expr_free(mm); return NULL;
    }
    size_t nn = (size_t)n * (size_t)n;
    Expr** flat = (Expr**)malloc(sizeof(Expr*) * nn);
    if (!flat) { expr_free(mm); return NULL; }
    for (int64_t i = 0; i < n; i++) {
        const Expr* row = rows->data.function.args[i];
        if (row->type != EXPR_FUNCTION || row->data.function.arg_count != (size_t)n) {
            free(flat); expr_free(mm); return NULL;
        }
        for (int64_t j = 0; j < n; j++) flat[i * n + j] = row->data.function.args[j];
    }
    Expr** c = flint_mat_charpoly(flat, (int)n);
    free(flat);
    expr_free(mm);
    if (!c) return NULL;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)(n + 1));
    if (!terms) {
        for (int64_t k = 0; k <= n; k++) expr_free(c[k]);
        free(c);
        return NULL;
    }
    bool neg = (n & 1) != 0;
    for (int64_t k = 0; k <= n; k++) {
        Expr* coef = c[k];
        if (neg) coef = expr_new_function(expr_new_symbol(SYM_Times),
                                          (Expr*[]){ expr_new_integer(-1), coef }, 2);
        Expr* pw = (k == 0) ? NULL
                 : (k == 1) ? expr_copy((Expr*)var)
                 : expr_new_function(expr_new_symbol(SYM_Power),
                                     (Expr*[]){ expr_copy((Expr*)var), expr_new_integer(k) }, 2);
        terms[k] = pw ? expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ coef, pw }, 2)
                      : coef;
    }
    free(c);
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)(n + 1));
    free(terms);
    /* A non-symbol variable (x + 1, a number) is multiplied out exactly as the
     * generic path does it, with Expand. */
    if (var->type != EXPR_SYMBOL)
        sum = expr_new_function(expr_new_symbol(SYM_Expand), (Expr*[]){ sum }, 1);
    return eval_and_free(sum);
}

Expr* builtin_characteristicpolynomial(Expr* res) {
    size_t argc = (res->type == EXPR_FUNCTION) ? res->data.function.arg_count : 0;
    if (argc != 2)
        return builtin_arg_error("CharacteristicPolynomial", argc, 2, 2);

    /* A visible NDArray first argument is delisted and re-evaluated, mirroring
     * builtin_eigenvalues; the boxed-List form then flows through the path
     * below.  (NDArrays nested inside a {m, a} List are handled per-matrix in
     * the generalized branch.) */
    if (linalg_call_has_ndarray(res)) return linalg_delist_and_reeval(res);

    Expr* arg0 = res->data.function.args[0];
    Expr* var  = res->data.function.args[1];

    Expr *m, *a; int64_t n;
    if (!eigen_extract_matrix_pair(arg0, &m, &a, &n))
        return NULL;   /* not a square matrix / matrix pair: leave unevaluated */

    const char* lam = eigen_lambda_name();

    Expr* poly;
    if (a == NULL) {
        Expr* fast = charpoly_flint(m, n, var);
        if (fast) return fast;
        /* Ordinary: Faddeev-Leverrier gives det(lambda I - m). */
        poly = eigen_char_poly_faddeev(m, lam, (int)n);
        if (!poly) return NULL;
    } else {
        /* Generalized: det(m - lambda a) via Laplace expansion (correct sign).
         * m / a may have arrived as raw NDArrays nested in the {m, a} List
         * (past the top-level guard); materialise them, as faddeev does for the
         * ordinary path.  pack_unpack returns NULL for a plain List, in which
         * case a cheap refcount-bump copy is used. */
        Expr* mm = pack_unpack(m); if (!mm) mm = expr_copy(m);
        Expr* aa = pack_unpack(a); if (!aa) aa = expr_copy(a);
        Expr* M = eigen_build_lambda_matrix(mm, aa, lam, n);
        expr_free(mm); expr_free(aa);
        poly = eigen_compute_det(M, (int)n);
        expr_free(M);
        if (!poly) return NULL;
    }

    /* Sign fix for the ordinary path: Det[m - x I] = (-1)^n det(x I - m). */
    Expr* signed_poly = poly;
    if (a == NULL && (n & 1)) {
        signed_poly = expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_integer(-1), poly }, 2);
    }

    /* Expand[ReplaceAll[poly, lambda -> var]] -- substitute the user's variable
     * for the internal lambda, then multiply out.  Built as one tree and
     * evaluated once. */
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol(lam), expr_copy(var) }, 2);
    Expr* subst = expr_new_function(expr_new_symbol(SYM_ReplaceAll),
        (Expr*[]){ signed_poly, rule }, 2);
    Expr* out = expr_new_function(expr_new_symbol(SYM_Expand),
        (Expr*[]){ subst }, 1);
    return eval_and_free(out);
}
