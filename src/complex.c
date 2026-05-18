#include "complex.h"
#include "symtab.h"
#include "eval.h"
#include "arithmetic.h"
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

void complex_init(void) {
    symtab_add_builtin("Re", builtin_re);
    symtab_add_builtin("Im", builtin_im);
    symtab_add_builtin("ReIm", builtin_reim);
    symtab_add_builtin("Abs", builtin_abs);
    symtab_add_builtin("Conjugate", builtin_conjugate);
    symtab_add_builtin("Arg", builtin_arg);

    symtab_get_def("Re")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Im")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("ReIm")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Abs")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Conjugate")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Arg")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
}

Expr* builtin_re(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        return expr_copy(re);
    }
    int64_t n, d;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || is_rational(arg, &n, &d)) {
        return expr_copy(arg);
    }
    return NULL;
}

Expr* builtin_im(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        return expr_copy(im);
    }
    int64_t n, d;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || is_rational(arg, &n, &d)) {
        return expr_new_integer(0);
    }
    return NULL;
}

Expr* builtin_reim(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        Expr** results = malloc(sizeof(Expr*) * 2);
        results[0] = expr_copy(re);
        results[1] = expr_copy(im);
        return expr_new_function(expr_new_symbol("List"), results, 2);
    }
    int64_t n, d;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || is_rational(arg, &n, &d)) {
        Expr** results = malloc(sizeof(Expr*) * 2);
        results[0] = expr_copy(arg);
        results[1] = expr_new_integer(0);
        return expr_new_function(expr_new_symbol("List"), results, 2);
    }
    return NULL;
}

Expr* builtin_abs(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        Expr* two1 = expr_new_integer(2);
        Expr* two2 = expr_new_integer(2);
        
        Expr** args_re = malloc(sizeof(Expr*) * 2);
        args_re[0] = expr_copy(re);
        args_re[1] = two1;
        Expr* pow_re = expr_new_function(expr_new_symbol("Power"), args_re, 2);
        free(args_re);
        
        Expr** args_im = malloc(sizeof(Expr*) * 2);
        args_im[0] = expr_copy(im);
        args_im[1] = two2;
        Expr* pow_im = expr_new_function(expr_new_symbol("Power"), args_im, 2);
        free(args_im);
        
        Expr** args_plus = malloc(sizeof(Expr*) * 2);
        args_plus[0] = pow_re;
        args_plus[1] = pow_im;
        Expr* sum = expr_new_function(expr_new_symbol("Plus"), args_plus, 2);
        free(args_plus);
        
        Expr* rat_args[2] = { expr_new_integer(1), expr_new_integer(2) };
        Expr* half = expr_new_function(expr_new_symbol("Rational"), rat_args, 2);
        Expr* pow_args[2] = { sum, half };
        return expr_new_function(expr_new_symbol("Power"), pow_args, 2);
    }
    int64_t n, d;
    if (arg->type == EXPR_BIGINT) {
        mpz_t r;
        mpz_init(r);
        mpz_abs(r, arg->data.bigint);
        Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
        mpz_clear(r);
        return result;
    }
    if (arg->type == EXPR_INTEGER) {
        int64_t val = arg->data.integer;
        return expr_new_integer(val < 0 ? -val : val);
    }
    if (arg->type == EXPR_REAL) {
        double val = arg->data.real;
        return expr_new_real(val < 0.0 ? -val : val);
    }
    if (is_rational(arg, &n, &d)) {
        return make_rational(n < 0 ? -n : n, d);
    }
    return NULL;
}

Expr* builtin_conjugate(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        Expr** args_times = malloc(sizeof(Expr*) * 2);
        args_times[0] = expr_new_integer(-1);
        args_times[1] = expr_copy(im);
        Expr* minus_im = expr_new_function(expr_new_symbol("Times"), args_times, 2);
        /* expr_new_function copies the pointers out of args_times into its
         * own internal storage; the caller still owns the temporary array. */
        free(args_times);
        return make_complex(expr_copy(re), minus_im);
    }
    int64_t n, d;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || is_rational(arg, &n, &d)) {
        return expr_copy(arg);
    }
    return NULL;
}

Expr* builtin_arg(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    
    Expr* re = NULL;
    Expr* im = NULL;
    int64_t n, d;
    
    if (is_complex(arg, &re, &im)) {
        // re and im are assigned
    } else if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || is_rational(arg, &n, &d)) {
        re = arg;
    } else {
        return NULL;
    }
    
    double re_val = 0.0;
    double im_val = 0.0;
    bool inexact = false;
    
    if (re) {
        if (re->type == EXPR_INTEGER) re_val = (double)re->data.integer;
        else if (re->type == EXPR_REAL) { re_val = re->data.real; inexact = true; }
        else if (is_rational(re, &n, &d)) re_val = (double)n / d;
        else return NULL;
    }
    if (im) {
        if (im->type == EXPR_INTEGER) im_val = (double)im->data.integer;
        else if (im->type == EXPR_REAL) { im_val = im->data.real; inexact = true; }
        else if (is_rational(im, &n, &d)) im_val = (double)n / d;
        else return NULL;
    }
    
    if (re_val == 0.0 && im_val == 0.0) {
        return expr_new_integer(0);
    }
    
    if (!inexact) {
        if (im_val == 0.0) {
            if (re_val > 0.0) return expr_new_integer(0);
            else return expr_new_symbol("Pi");
        }
        if (re_val == 0.0) {
            Expr* args[2];
            args[0] = make_rational(im_val > 0.0 ? 1 : -1, 2);
            args[1] = expr_new_symbol("Pi");
            return expr_new_function(expr_new_symbol("Times"), args, 2);
        }
        if (re_val == im_val) {
            Expr* args[2];
            args[0] = make_rational(re_val > 0.0 ? 1 : -3, 4);
            args[1] = expr_new_symbol("Pi");
            return expr_new_function(expr_new_symbol("Times"), args, 2);
        }
        if (re_val == -im_val) {
            Expr* args[2];
            args[0] = make_rational(re_val > 0.0 ? -1 : 3, 4);
            args[1] = expr_new_symbol("Pi");
            return expr_new_function(expr_new_symbol("Times"), args, 2);
        }

        Expr* args[2];
        args[0] = expr_copy(re);
        args[1] = im ? expr_copy(im) : expr_new_integer(0);
        return expr_new_function(expr_new_symbol("ArcTan"), args, 2);
    }
    
    return expr_new_real(atan2(im_val, re_val));
}
