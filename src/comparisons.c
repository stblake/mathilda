/*
 * Mathilda Comparisons
 *
 * Implements structural comparison functions similar to the Wolfram Language.
 */

#include "comparisons.h"
#include "symtab.h"
#include "expr.h"
#include "arithmetic.h"
#include "eval.h"
#include "numeric.h"
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
        return expr_new_symbol("True");
    }
    
    Expr* first = res->data.function.args[0];
    for (size_t i = 1; i < res->data.function.arg_count; i++) {
        if (!expr_eq(first, res->data.function.args[i])) {
            return expr_new_symbol("False");
        }
    }
    
    return expr_new_symbol("True");
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
        return expr_new_symbol("True");
    }
    
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        for (size_t j = i + 1; j < res->data.function.arg_count; j++) {
            if (expr_eq(res->data.function.args[i], res->data.function.args[j])) {
                return expr_new_symbol("False");
            }
        }
    }
    
    return expr_new_symbol("True");
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
        return expr_new_symbol("True");
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
                return expr_new_symbol("False");
            }
            all_equal = false;
        }
    }

    if (all_equal) {
        return expr_new_symbol("True");
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
        return expr_new_symbol("True");
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
                return expr_new_symbol("False");
            }
            if (!definitely_unequal) {
                all_definitely_unequal = false;
            }
        }
    }
    
    if (all_definitely_unequal) {
        return expr_new_symbol("True");
    }
    
    return NULL;
}

/* Inequality helpers */
static Expr* evaluate_inequality(Expr* res, int expected_cmp_1, int expected_cmp_2) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 2) return expr_new_symbol("True");

    for (size_t i = 0; i < res->data.function.arg_count - 1; i++) {
        Expr* a = res->data.function.args[i];
        Expr* b = res->data.function.args[i+1];
        
        bool can_compare;
        int cmp = compare_numeric(a, b, &can_compare);
        
        if (!can_compare) {
            return NULL; // Unevaluated
        }
        
        if (cmp != expected_cmp_1 && cmp != expected_cmp_2) {
            return expr_new_symbol("False");
        }
    }
    return expr_new_symbol("True");
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
}
