#include "list_common.h"

bool is_overflow(Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == SYM_Overflow;
}

bool is_listq(Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == SYM_List;
}

bool is_infinity(Expr* e) {
    return e->type == EXPR_SYMBOL && e->data.symbol == SYM_Infinity;
}

bool is_minus_infinity(Expr* e) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Times &&
        e->data.function.arg_count == 2) {
        Expr* a1 = e->data.function.args[0];
        Expr* a2 = e->data.function.args[1];
        if (a1->type == EXPR_INTEGER && a1->data.integer == -1 && is_infinity(a2)) return true;
        if (a2->type == EXPR_INTEGER && a2->data.integer == -1 && is_infinity(a1)) return true;
    }
    return false;
}

Expr* make_minus_infinity(void) {
    Expr* args[2] = { expr_new_integer(-1), expr_new_symbol(SYM_Infinity) };
    return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
}

bool is_real_numeric(Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
    if (is_rational(e, NULL, NULL)) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    Expr* re, *im;
    if (is_complex(e, &re, &im)) {
        if (im->type == EXPR_INTEGER && im->data.integer == 0) return true;
        if (im->type == EXPR_REAL && im->data.real == 0.0) return true;
        if (im->type == EXPR_BIGINT && mpz_sgn(im->data.bigint) == 0) return true;
#ifdef USE_MPFR
        if (im->type == EXPR_MPFR && mpfr_zero_p(im->data.mpfr)) return true;
#endif
    }
    return false;
}
