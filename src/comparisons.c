/*
 * Mathilda Comparisons
 *
 * Implements structural comparison functions similar to the Wolfram Language.
 */

#include "comparisons.h"
#include "symtab.h"
#include "sym_names.h"
#include "attr.h"
#include "expr.h"
#include "arithmetic.h"
#include "eval.h"
#include "numeric.h"
#include <gmp.h>
#include <stdbool.h>
#include <string.h>

static bool get_numeric_value(Expr* e, double* val, int64_t* num, int64_t* den, bool* is_exact) {
    if (e->type == EXPR_INTEGER) {
        *val = (double)e->data.integer;
        *num = e->data.integer;
        *den = 1;
        *is_exact = true;
        return true;
    } else if (e->type == EXPR_REAL) {
        *val = e->data.real;
        *is_exact = false;
        return true;
    } else if (is_rational(e, num, den)) {
        *val = (double)(*num) / (*den);
        *is_exact = true;
        return true;
    }

    /* Fallback: drive the expression through N[] (machine precision).
     * Captures NumericQ-style values like Sqrt[3], Pi, or 1 + Sqrt[2]/3
     * that aren't raw numeric literals but have a definite real value.
     * If numericalize can't reduce the expression to a real number, we
     * leave it as unevaluated for the caller. The result is treated as
     * inexact so comparisons use the floating-point tolerance branch. */
    Expr* approx = numericalize(e, numeric_machine_spec());
    if (!approx) return false;
    bool ok = false;
    if (approx->type == EXPR_INTEGER) {
        *val = (double)approx->data.integer;
        *num = approx->data.integer;
        *den = 1;
        *is_exact = false;
        ok = true;
    } else if (approx->type == EXPR_REAL) {
        *val = approx->data.real;
        *is_exact = false;
        ok = true;
    } else if (is_rational(approx, num, den)) {
        *val = (double)(*num) / (*den);
        *is_exact = false;
        ok = true;
    }
    expr_free(approx);
    return ok;
}

static int compare_numeric(Expr* a, Expr* b, bool* can_compare) {
    double va, vb;
    int64_t na, da, nb, db;
    bool exact_a, exact_b;

    *can_compare = false;

    /* Exact integer-like comparison via GMP. Avoids precision loss when
     * either operand is a BigInt whose magnitude exceeds 2^53 (the limit
     * for distinguishing adjacent values through double). Without this,
     * Less[10^30, 10^30 + 1] coerces both sides to the same double and
     * answers False. */
    if (expr_is_integer_like(a) && expr_is_integer_like(b)) {
        mpz_t ma, mb;
        mpz_init(ma); mpz_init(mb);
        expr_to_mpz(a, ma);
        expr_to_mpz(b, mb);
        int cmp = mpz_cmp(ma, mb);
        mpz_clear(ma); mpz_clear(mb);
        *can_compare = true;
        if (cmp < 0) return -1;
        if (cmp > 0) return 1;
        return 0;
    }

    if (!get_numeric_value(a, &va, &na, &da, &exact_a)) return 0;
    if (!get_numeric_value(b, &vb, &nb, &db, &exact_b)) return 0;

    *can_compare = true;
    if (exact_a && exact_b) {
        long double val_a = (long double)na / da;
        long double val_b = (long double)nb / db;
        if (val_a < val_b) return -1;
        if (val_a > val_b) return 1;
        return 0;
    } else {
        double diff = va - vb;
        if (diff < 0) diff = -diff;
        double max_val = (va < 0 ? -va : va);
        if ((vb < 0 ? -vb : vb) > max_val) max_val = (vb < 0 ? -vb : vb);
        
        // 2^(-46) = 1.4210854715202004e-14
        if (diff <= max_val * 1.4210854715202004e-14 || diff == 0.0) {
            return 0;
        }
        if (va < vb) return -1;
        return 1;
    }
}

static bool is_raw_data(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_STRING) return true;
    int64_t n, d;
    if (is_rational(e, &n, &d)) return true;
    return false;
}

/*
 * builtin_sameq: Implements SameQ[x, y, ...].
 * SameQ returns True if all its arguments are identical, and False otherwise.
 * Identical means they have the same structure and atoms.
 */
Expr* builtin_sameq(Expr* res) {
    if (res->type != EXPR_FUNCTION) {
        return NULL;
    }
    
    // SameQ[] and SameQ[x] return True by convention.
    if (res->data.function.arg_count < 2) {
        return expr_new_symbol(SYM_True);
    }
    
    Expr* first = res->data.function.args[0];
    for (size_t i = 1; i < res->data.function.arg_count; i++) {
        if (!expr_eq(first, res->data.function.args[i])) {
            return expr_new_symbol(SYM_False);
        }
    }
    
    return expr_new_symbol(SYM_True);
}

/*
 * builtin_unsameq: Implements UnsameQ[e1, e2, ...].
 * UnsameQ returns True if no two of the ei are identical, and False otherwise.
 */
Expr* builtin_unsameq(Expr* res) {
    if (res->type != EXPR_FUNCTION) {
        return NULL;
    }
    
    // UnsameQ[] and UnsameQ[x] return True by convention.
    if (res->data.function.arg_count < 2) {
        return expr_new_symbol(SYM_True);
    }
    
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        for (size_t j = i + 1; j < res->data.function.arg_count; j++) {
            if (expr_eq(res->data.function.args[i], res->data.function.args[j])) {
                return expr_new_symbol(SYM_False);
            }
        }
    }
    
    return expr_new_symbol(SYM_True);
}

/*
 * builtin_equal: Implements Equal[lhs, rhs, ...].
 * Equal returns True if its arguments are identical (SameQ) or numerically equal (2 == 2.0).
 * Otherwise it returns NULL (unevaluated).
 */
Expr* builtin_equal(Expr* res) {
    if (res->type != EXPR_FUNCTION) {
        return NULL;
    }
    
    // Equal[] and Equal[x] return True.
    if (res->data.function.arg_count < 2) {
        return expr_new_symbol(SYM_True);
    }
    
    bool all_equal = true;
    for (size_t i = 0; i < res->data.function.arg_count - 1; i++) {
        Expr* a = res->data.function.args[i];
        Expr* b = res->data.function.args[i+1];

        bool equal = false;
        bool definitely_unequal = false;

        if (expr_eq(a, b)) {
            equal = true;
        } else {
            bool can_compare = false;
            int cmp = compare_numeric(a, b, &can_compare);
            if (can_compare) {
                if (cmp == 0) equal = true;
                else definitely_unequal = true;
            }
        }

        if (!equal) {
            if (definitely_unequal || (is_raw_data(a) && is_raw_data(b))) {
                return expr_new_symbol(SYM_False);
            }
            all_equal = false;
        }
    }

    if (all_equal) {
        return expr_new_symbol(SYM_True);
    }

    return NULL;
}

/*
 * builtin_unequal: Implements Unequal[e1, e2, ...].
 * Unequal returns False if any two arguments are numerically equal,
 * and True if all pairs are determined to be unequal.
 */
Expr* builtin_unequal(Expr* res) {
    if (res->type != EXPR_FUNCTION) {
        return NULL;
    }
    
    // Unequal[] and Unequal[x] return True.
    if (res->data.function.arg_count < 2) {
        return expr_new_symbol(SYM_True);
    }
    
    bool all_definitely_unequal = true;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        for (size_t j = i + 1; j < res->data.function.arg_count; j++) {
            Expr* a = res->data.function.args[i];
            Expr* b = res->data.function.args[j];
            
            bool equal = false;
            bool definitely_unequal = false;
            
            // Structural identity
            if (expr_eq(a, b)) {
                equal = true;
            } else {
                bool can_compare;
                int cmp = compare_numeric(a, b, &can_compare);
                if (can_compare) {
                    if (cmp == 0) {
                        equal = true;
                    } else {
                        definitely_unequal = true;
                    }
                } else if (is_raw_data(a) && is_raw_data(b)) {
                    // Different raw data types or values (and not equal numerically)
                    definitely_unequal = true;
                }
            }
            
            if (equal) {
                return expr_new_symbol(SYM_False);
            }
            if (!definitely_unequal) {
                all_definitely_unequal = false;
            }
        }
    }
    
    if (all_definitely_unequal) {
        return expr_new_symbol(SYM_True);
    }
    
    return NULL;
}

/* Inequality helpers */
static Expr* evaluate_inequality(Expr* res, int expected_cmp_1, int expected_cmp_2) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 2) return expr_new_symbol(SYM_True);

    for (size_t i = 0; i < res->data.function.arg_count - 1; i++) {
        Expr* a = res->data.function.args[i];
        Expr* b = res->data.function.args[i+1];
        
        bool can_compare;
        int cmp = compare_numeric(a, b, &can_compare);
        
        if (!can_compare) {
            return NULL; // Unevaluated
        }
        
        if (cmp != expected_cmp_1 && cmp != expected_cmp_2) {
            return expr_new_symbol(SYM_False);
        }
    }
    return expr_new_symbol(SYM_True);
}

Expr* builtin_less(Expr* res) {
    return evaluate_inequality(res, -1, -1);
}

Expr* builtin_greater(Expr* res) {
    return evaluate_inequality(res, 1, 1);
}

Expr* builtin_lessequal(Expr* res) {
    return evaluate_inequality(res, -1, 0);
}

Expr* builtin_greaterequal(Expr* res) {
    return evaluate_inequality(res, 1, 0);
}

/* Map an operator-symbol head (Less/LessEqual/Greater/GreaterEqual/Equal)
 * to a pairwise decision. Returns:
 *    1  if `a OP b` is definitely True,
 *    0  if it is definitely False,
 *   -1  if the comparison is not decidable (operand non-numeric or the
 *       op is not a recognised chain head).
 * Equal/Unequal share a small fast path so that, e.g., x == x is True
 * even when x has no numeric value. */
static int decide_pair(const char* op, Expr* a, Expr* b) {
    if (op == SYM_Equal) {
        if (expr_eq(a, b)) return 1;
        bool can; int cmp = compare_numeric(a, b, &can);
        if (!can) return -1;
        return cmp == 0 ? 1 : 0;
    }
    if (op == SYM_Less || op == SYM_LessEqual
     || op == SYM_Greater || op == SYM_GreaterEqual) {
        bool can; int cmp = compare_numeric(a, b, &can);
        if (!can) return -1;
        if (op == SYM_Less)         return cmp <  0 ? 1 : 0;
        if (op == SYM_LessEqual)    return cmp <= 0 ? 1 : 0;
        if (op == SYM_Greater)      return cmp >  0 ? 1 : 0;
        /* GreaterEqual */          return cmp >= 0 ? 1 : 0;
    }
    return -1;
}

/* Inequality[v0, op0, v1, op1, v2, ...]: Mathematica's variadic chained
 * comparison head, produced by the parser for `a < b <= c == d > e` and
 * similar chains. Semantics:
 *   - True iff every adjacent pair `v_i op_i v_{i+1}` is True.
 *   - False if any adjacent pair is decidably False.
 *   - Otherwise, drop pairs that we can prove True, and either collapse
 *     to a single binary head (when one pair remains) or rebuild a
 *     smaller Inequality with the surviving pairs.
 *
 * The argument layout is strict: 2n+1 args, even slots are values, odd
 * slots are bare operator symbols. Anything else returns NULL so the
 * caller's input is preserved unevaluated. */
Expr* builtin_inequality(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return expr_new_symbol(SYM_True);
    if (n == 1) return expr_new_symbol(SYM_True);  /* one value, no pair */
    if ((n & 1u) == 0) return NULL;              /* must be 2k+1 */

    /* Validate operator slots upfront — every odd index must be one of the
     * five chain heads. */
    for (size_t i = 1; i < n; i += 2) {
        Expr* op = res->data.function.args[i];
        if (op->type != EXPR_SYMBOL) return NULL;
        const char* s = op->data.symbol.name;
        if (s != SYM_Less && s != SYM_LessEqual
         && s != SYM_Greater && s != SYM_GreaterEqual
         && s != SYM_Equal) return NULL;
    }

    /* Decide each pair; any False short-circuits the whole chain. We also
     * remember which pairs were "undecided" so we can rebuild a residual
     * Inequality if needed. */
    size_t npairs = (n - 1) / 2;
    bool* undecided = (bool*)calloc(npairs, sizeof(bool));
    if (!undecided) return NULL;
    bool any_undecided = false;
    for (size_t k = 0; k < npairs; k++) {
        Expr* a  = res->data.function.args[2*k];
        Expr* op = res->data.function.args[2*k + 1];
        Expr* b  = res->data.function.args[2*k + 2];
        int d = decide_pair(op->data.symbol.name, a, b);
        if (d == 0) {                    /* definitely False */
            free(undecided);
            return expr_new_symbol(SYM_False);
        }
        if (d < 0) {
            undecided[k] = true;
            any_undecided = true;
        }
    }
    if (!any_undecided) {
        free(undecided);
        return expr_new_symbol(SYM_True);
    }

    /* Build a residual Inequality from the undecided pairs. Each surviving
     * pair contributes (lhs, op, rhs). When the previous emitted RHS is
     * structurally equal to the next pair's LHS (the common case when
     * undecided pairs are adjacent, but also when a dropped pair was
     * `b == b`), we omit the duplicate value so the chain stays well-
     * formed (2n+1 args). */
    size_t cap = 2 * npairs + 1;
    Expr** out_args = (Expr**)malloc(sizeof(Expr*) * cap);
    if (!out_args) { free(undecided); return NULL; }
    size_t out_n = 0;
    Expr* last_emitted_val = NULL;
    for (size_t k = 0; k < npairs; k++) {
        if (!undecided[k]) continue;
        Expr* a  = res->data.function.args[2*k];
        Expr* op = res->data.function.args[2*k + 1];
        Expr* b  = res->data.function.args[2*k + 2];
        if (!last_emitted_val || !expr_eq(last_emitted_val, a)) {
            out_args[out_n++] = expr_copy(a);
        }
        out_args[out_n++] = expr_copy(op);
        out_args[out_n++] = expr_copy(b);
        last_emitted_val = b;
    }
    free(undecided);

    if (out_n == 3) {
        /* Exactly one undecided pair survived — return it as the binary
         * comparison so downstream consumers see the familiar shape. */
        Expr* binary_args[2] = { out_args[0], out_args[2] };
        const char* head_sym = out_args[1]->data.symbol.name;
        expr_free(out_args[1]);
        Expr* result = expr_new_function(expr_new_symbol(head_sym), binary_args, 2);
        free(out_args);
        return result;
    }
    Expr* residual = expr_new_function(expr_new_symbol(SYM_Inequality), out_args, out_n);
    free(out_args);
    return residual;
}

/*
 * comparisons_init: Registers comparison-related builtins in the symbol table.
 */
void comparisons_init(void) {
    symtab_add_builtin("SameQ", builtin_sameq);
    symtab_add_builtin("UnsameQ", builtin_unsameq);
    symtab_add_builtin("Equal", builtin_equal);
    symtab_add_builtin("Unequal", builtin_unequal);
    symtab_add_builtin("Less", builtin_less);
    symtab_add_builtin("Greater", builtin_greater);
    symtab_add_builtin("LessEqual", builtin_lessequal);
    symtab_add_builtin("GreaterEqual", builtin_greaterequal);
    symtab_add_builtin("Inequality", builtin_inequality);
    set_attributes("Inequality", ATTR_PROTECTED);
}
