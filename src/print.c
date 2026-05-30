#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

// Print

#include "print.h"
#include "expr.h"
#include "arithmetic.h"
#include "context.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include <inttypes.h>

static void print_standard(Expr* e, int parent_prec);
static void print_series_data(Expr* e, int parent_prec);
static void print_tex(Expr* e, int parent_prec);

/* Returns a fresh expression holding -b, where b is a bigint node. The result
 * is normalised back to a machine int when it fits. Used by the Plus printers
 * to pull the sign of a negative bigint term onto the operator. */
static Expr* negate_bigint(Expr* b) {
    mpz_t n;
    mpz_init(n);
    mpz_neg(n, b->data.bigint);
    Expr* r = expr_new_bigint_from_mpz(n);
    mpz_clear(n);
    return expr_bigint_normalize(r);
}

/*
 * When non-zero, print_standard renders heads like SeriesData using their
 * literal head[args...] form rather than the pretty mathematical form.
 * Toggled while printing inside an InputForm[...] wrapper.
 */
static int g_inputform_depth = 0;

static int get_expr_prec(Expr* e) {
    if (e->type != EXPR_FUNCTION) return 1000;
    if (e->data.function.head->type != EXPR_SYMBOL) return 1000;
    
    const char* head = e->data.function.head->data.symbol;

    if (head == SYM_Set || head == SYM_SetDelayed) return 40;
    if (head == SYM_Rule || head == SYM_RuleDelayed) return 120;
    if (head == SYM_Condition) return 130;
    if (head == SYM_Alternatives) return 160;
    if (head == SYM_Repeated || head == SYM_RepeatedNull) return 170;
    if (head == SYM_And) return 215;
    if (head == SYM_Or) return 215;
    if (head == SYM_Equal || head == SYM_Unequal || head == SYM_Less || head == SYM_Greater || head == SYM_LessEqual || head == SYM_GreaterEqual || head == SYM_SameQ || head == SYM_UnsameQ || head == SYM_Inequality) return 290;
    if (head == SYM_Plus) return 310;

    if (head == SYM_Times) return 400;
    if (head == SYM_Divide) return 470;
    if (head == SYM_Power) {
        if (e->data.function.arg_count == 2) {
            Expr* exp = e->data.function.args[1];
            if (exp->type == EXPR_INTEGER && exp->data.integer < 0) return 470;
            int64_t n, d;
            if (is_rational(exp, &n, &d) && n < 0) return 470;
            if (exp->type == EXPR_REAL && exp->data.real < 0.0) return 470;
        }
        return 590;
    }
    if (head == SYM_Rational) return 470;
    if (head == SYM_Complex) {
        /* Only the bare imaginary unit `I` (Complex[0, 1]) is an atom.
         * `-I` (Complex[0, -1]) prints as a unary-minus form and so
         * takes Plus-level precedence — otherwise Power[-I, n] would
         * print as "-I^n" (i.e. -(I^n)) instead of "(-I)^n". General
         * Complex[a, b] also prints at Plus level. */
        if (e->data.function.arg_count == 2) {
            Expr* re = e->data.function.args[0];
            Expr* im = e->data.function.args[1];
            bool re_zero = (re->type == EXPR_INTEGER && re->data.integer == 0)
                        || (re->type == EXPR_REAL && re->data.real == 0.0);
            if (re_zero && im->type == EXPR_INTEGER && im->data.integer == 1) {
                return 1000;
            }
        }
        return 310;
    }
    /* Function[body] is the slot-based pure function rendered "body &".
     * Mathematica gives "&" precedence 90 (just above Postfix). */
    if (head == SYM_Function && e->data.function.arg_count == 1) return 90;
    return 1000;
}

static void print_function_fullform(Expr* e) {
    expr_print_fullform(e->data.function.head);
    printf("[");
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (i > 0) printf(", ");
        expr_print_fullform(e->data.function.args[i]);
    }    
    printf("]");
}

void expr_print_fullform(Expr* e) {
    if (!e) { printf("Null"); return; }
    switch (e->type) {
        case EXPR_INTEGER: printf("%" PRId64, e->data.integer); break;
        case EXPR_REAL: {
            char buf[64];
            snprintf(buf, sizeof(buf), "%g", e->data.real);
            if (strchr(buf, '.') == NULL && strchr(buf, 'e') == NULL && strchr(buf, 'E') == NULL) printf("%s.0", buf);
            else printf("%s", buf);
            break;
        }
        case EXPR_SYMBOL: printf("%s", context_display_name(e->data.symbol)); break;
        case EXPR_STRING: printf("\"%s\"", e->data.string); break;
        case EXPR_FUNCTION: print_function_fullform(e); break;
        case EXPR_BIGINT: {
            char* str = mpz_get_str(NULL, 10, e->data.bigint);
            printf("%s", str);
            free(str);
            break;
        }
#ifdef USE_MPFR
        case EXPR_MPFR: {
            /* Print with decimal-digit precision matching the value's
             * stored binary precision. `%.*Rg` uses the `g` style (drops
             * trailing zeros) so most common values print cleanly; we
             * append ".0" for integer-valued outputs for consistency
             * with EXPR_REAL formatting (e.g. 3. → 3.0). */
            long digits = (long)ceil((double)mpfr_get_prec(e->data.mpfr)
                                     / 3.3219280948873626);
            if (digits < 1) digits = 1;
            char* buf = NULL;
            int len = mpfr_asprintf(&buf, "%.*Rg", (int)digits, e->data.mpfr);
            if (len >= 0 && buf) {
                printf("%s", buf);
                if (strchr(buf, '.') == NULL && strchr(buf, 'e') == NULL
                    && strchr(buf, 'E') == NULL && strchr(buf, 'n') == NULL
                    && strchr(buf, 'i') == NULL /* nan / inf */) {
                    printf(".0");
                }
                mpfr_free_str(buf);
            } else {
                printf("<mpfr print error>");
            }
            break;
        }
#endif
    }
}

static void print_standard(Expr* e, int parent_prec) {
    if (!e) { printf("Null"); return; }
    if (e->type != EXPR_FUNCTION) {
        expr_print_fullform(e);
        return;
    }

    int my_prec = get_expr_prec(e);
    bool need_parens = (my_prec < parent_prec);

    if (need_parens) printf("(");

    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* head = e->data.function.head->data.symbol;
        
        if (head == SYM_FullForm && e->data.function.arg_count == 1) {
            expr_print_fullform(e->data.function.args[0]);
        }
        else if (head == SYM_TeXForm && e->data.function.arg_count == 1) {
            print_tex(e->data.function.args[0], 0);
        }
        else if (head == SYM_InputForm && e->data.function.arg_count == 1) {
            g_inputform_depth++;
            print_standard(e->data.function.args[0], parent_prec);
            g_inputform_depth--;
        }
        else if (head == SYM_SeriesData && e->data.function.arg_count == 6 && g_inputform_depth == 0) {
            print_series_data(e, parent_prec);
        }
        else if (head == SYM_HoldForm && e->data.function.arg_count == 1) {
            print_standard(e->data.function.args[0], parent_prec);
        }
        else if (head == SYM_Rational && e->data.function.arg_count == 2) {
            print_standard(e->data.function.args[0], 470);
            printf("/");
            print_standard(e->data.function.args[1], 470);
        }
        else if (head == SYM_Complex && e->data.function.arg_count == 2) {
            Expr* re = e->data.function.args[0];
            Expr* im = e->data.function.args[1];
            bool re_zero = (re->type == EXPR_INTEGER && re->data.integer == 0);
            bool im_one = (im->type == EXPR_INTEGER && im->data.integer == 1);
            bool im_minus_one = (im->type == EXPR_INTEGER && im->data.integer == -1);

            bool im_neg = im_minus_one;
            Expr* im_abs = NULL;
            int64_t rn, rd;
            if (!im_neg) {
                if (im->type == EXPR_INTEGER && im->data.integer < 0) {
                    im_neg = true;
                    im_abs = expr_new_integer(-im->data.integer);
                } else if (im->type == EXPR_REAL && im->data.real < 0.0) {
                    im_neg = true;
                    im_abs = expr_new_real(-im->data.real);
                } else if (is_rational(im, &rn, &rd) && rn < 0) {
                    im_neg = true;
                    im_abs = make_rational(-rn, rd);
                }
            }

            if (!re_zero) {
                print_standard(re, 310);
                printf(im_neg ? " - " : " + ");
            } else if (im_neg) {
                printf("-");
            }
            if (im_one || im_minus_one) printf("I");
            else if (im_abs) { print_standard(im_abs, 400); printf("*I"); expr_free(im_abs); }
            else { print_standard(im, 400); printf("*I"); }
        }
        else if (head == SYM_Power && e->data.function.arg_count == 2) {
            Expr* exp = e->data.function.args[1];
            int64_t n, d;
            bool is_negative = false;
            Expr* pos_exp = NULL;
            if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
                is_negative = true;
                pos_exp = expr_new_integer(-exp->data.integer);
            } else if (is_rational(exp, &n, &d) && n < 0) {
                is_negative = true;
                pos_exp = make_rational(-n, d);
            } else if (exp->type == EXPR_REAL && exp->data.real < 0.0) {
                is_negative = true;
                pos_exp = expr_new_real(-exp->data.real);
            }

            if (is_negative) {
                printf("1/");
                if (pos_exp->type == EXPR_INTEGER && pos_exp->data.integer == 1) {
                    print_standard(e->data.function.args[0], 470);
                } else {
                    Expr* args[2] = { expr_copy(e->data.function.args[0]), expr_copy(pos_exp) };
                    Expr* tmp = expr_new_function(expr_new_symbol("Power"), args, 2);
                    print_standard(tmp, 470);
                    expr_free(tmp);
                }
                if (pos_exp) expr_free(pos_exp);
            } else if (is_rational(exp, &n, &d) && n == 1 && d == 2) {
                printf("Sqrt[");
                print_standard(e->data.function.args[0], 0);
                printf("]");
            } else {
                /* Negative numeric atoms parse as unary-minus, which has lower
                 * precedence than Power. So `(-1)^(1/3)` round-trips as
                 * `-1^(1/3)` -> `-(1^(1/3))`. Wrap explicitly. */
                Expr* base = e->data.function.args[0];
                bool base_paren = false;
                int64_t bn, bd;
                if (base->type == EXPR_INTEGER && base->data.integer < 0) base_paren = true;
                else if (base->type == EXPR_REAL && base->data.real < 0.0) base_paren = true;
                else if (base->type == EXPR_BIGINT && mpz_sgn(base->data.bigint) < 0) base_paren = true;
                else if (is_rational(base, &bn, &bd) && bn < 0) base_paren = true;
                if (base_paren) {
                    printf("(");
                    print_standard(base, 0);
                    printf(")");
                } else {
                    print_standard(base, 590);
                }
                printf("^");
                print_standard(e->data.function.args[1], 590);
            }
        }
        else if (head == SYM_Times) {
            size_t count = e->data.function.arg_count;

            /* Pull a leading negative numeric factor out as unary minus so
             * Times[-1, Infinity] prints "-Infinity", Times[-1, a, b] prints
             * "-a b", and Times[-3, x] prints "-3 x" (we keep the magnitude
             * as a factor when it isn't 1). Without this the numerator loop
             * below would emit "-1 Infinity". */
            int64_t lead_start = 0;
            bool flipped_sign = false;
            Expr* flipped_head = NULL; /* only used when lead arg != -1 */
            if (count > 0) {
                Expr* a0 = e->data.function.args[0];
                if (a0->type == EXPR_INTEGER && a0->data.integer < 0) {
                    flipped_sign = true;
                    if (a0->data.integer == -1 && count > 1) {
                        lead_start = 1; /* drop the -1, emit only "-" */
                    } else {
                        flipped_head = expr_new_integer(-a0->data.integer);
                    }
                } else if (a0->type == EXPR_REAL && a0->data.real < 0.0) {
                    flipped_sign = true;
                    if (a0->data.real == -1.0 && count > 1) {
                        lead_start = 1;
                    } else {
                        flipped_head = expr_new_real(-a0->data.real);
                    }
                } else if (a0->type == EXPR_FUNCTION &&
                           a0->data.function.head->type == EXPR_SYMBOL &&
                           a0->data.function.head->data.symbol == SYM_Rational &&
                           a0->data.function.arg_count == 2 &&
                           a0->data.function.args[0]->type == EXPR_INTEGER &&
                           a0->data.function.args[0]->data.integer < 0) {
                    flipped_sign = true;
                    Expr* num = a0->data.function.args[0];
                    Expr* den = a0->data.function.args[1];
                    Expr* rargs[2] = { expr_new_integer(-num->data.integer), expr_copy(den) };
                    flipped_head = expr_new_function(expr_new_symbol("Rational"), rargs, 2);
                } else if (a0->type == EXPR_FUNCTION &&
                           a0->data.function.head->type == EXPR_SYMBOL &&
                           a0->data.function.head->data.symbol == SYM_Complex &&
                           a0->data.function.arg_count == 2 &&
                           a0->data.function.args[0]->type == EXPR_INTEGER &&
                           a0->data.function.args[0]->data.integer == 0 &&
                           a0->data.function.args[1]->type == EXPR_INTEGER &&
                           a0->data.function.args[1]->data.integer < 0) {
                    /* Pure-imaginary-negative lead factor: Times[Complex[0, -n], ...]
                     * pulls the sign out so it prints as "-n I ..." instead of
                     * "(-n I) ..." or, for n == 1, "-I ..." instead of "(-I) ...". */
                    flipped_sign = true;
                    Expr* im = a0->data.function.args[1];
                    if (im->data.integer == -1 && count > 1) {
                        /* Drop the imaginary unit's -1; the loop will print
                         * Complex[0, 1] = "I" as the new lead factor. */
                        Expr* cargs[2] = { expr_new_integer(0), expr_new_integer(1) };
                        flipped_head = expr_new_function(expr_new_symbol("Complex"), cargs, 2);
                    } else {
                        Expr* cargs[2] = { expr_new_integer(0), expr_new_integer(-im->data.integer) };
                        flipped_head = expr_new_function(expr_new_symbol("Complex"), cargs, 2);
                    }
                }
            }
            if (flipped_sign) printf("-");

            Expr** num_args = malloc(sizeof(Expr*) * count);
            Expr** den_args = malloc(sizeof(Expr*) * count);
            size_t num_count = 0, den_count = 0;

            for (size_t i = lead_start; i < count; i++) {
                Expr* arg = (i == 0 && flipped_head) ? flipped_head : e->data.function.args[i];
                if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && arg->data.function.head->data.symbol == SYM_Power && arg->data.function.arg_count == 2) {
                    Expr* exp = arg->data.function.args[1];
                    bool is_neg = false;
                    int64_t n, d;
                    Expr* pos_exp = NULL;
                    if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
                        is_neg = true;
                        pos_exp = expr_new_integer(-exp->data.integer);
                    } else if (is_rational(exp, &n, &d) && n < 0) {
                        is_neg = true;
                        pos_exp = make_rational(-n, d);
                    } else if (exp->type == EXPR_REAL && exp->data.real < 0.0) {
                        is_neg = true;
                        pos_exp = expr_new_real(-exp->data.real);
                    }
                    if (is_neg) {
                        if (pos_exp->type == EXPR_INTEGER && pos_exp->data.integer == 1) {
                            den_args[den_count++] = expr_copy(arg->data.function.args[0]);
                        } else {
                            Expr* p_args[2] = { expr_copy(arg->data.function.args[0]), pos_exp };
                            den_args[den_count++] = expr_new_function(expr_new_symbol("Power"), p_args, 2);
                            pos_exp = NULL; 
                        }
                        if (pos_exp) expr_free(pos_exp);
                    } else {
                        num_args[num_count++] = expr_copy(arg);
                    }
                } else {
                    num_args[num_count++] = expr_copy(arg);
                }
            }
            
            if (den_count == 0) {
                for (size_t i = 0; i < num_count; i++) {
                    if (i > 0) printf(" ");
                    print_standard(num_args[i], 400);
                }
            } else {
                Expr* num = NULL;
                if (num_count == 0) num = expr_new_integer(1);
                else if (num_count == 1) num = expr_copy(num_args[0]);
                else {
                    Expr** nc = malloc(sizeof(Expr*) * num_count);
                    for(size_t i=0; i<num_count; i++) nc[i] = expr_copy(num_args[i]);
                    num = expr_new_function(expr_new_symbol("Times"), nc, num_count);
                    free(nc);
                }
                
                Expr* den = NULL;
                if (den_count == 1) den = expr_copy(den_args[0]);
                else {
                    Expr** dc = malloc(sizeof(Expr*) * den_count);
                    for(size_t i=0; i<den_count; i++) dc[i] = expr_copy(den_args[i]);
                    den = expr_new_function(expr_new_symbol("Times"), dc, den_count);
                    free(dc);
                }
                
                print_standard(num, 470);
                printf("/");
                print_standard(den, 470);
                
                expr_free(num);
                expr_free(den);
            }
            
            for (size_t i = 0; i < num_count; i++) expr_free(num_args[i]);
            free(num_args);
            for (size_t i = 0; i < den_count; i++) expr_free(den_args[i]);
            free(den_args);
            if (flipped_head) expr_free(flipped_head);
        }

        else if (head == SYM_Inequality
                 && e->data.function.arg_count >= 3
                 && (e->data.function.arg_count & 1u) == 1) {
            /* Inequality[v0, op0, v1, op1, v2, ...]. Walk the alternating
             * (value, op-symbol) sequence, printing each binary operator
             * using the same glyphs as the corresponding binary heads. */
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                Expr* a = e->data.function.args[i];
                if ((i & 1u) == 0) {
                    print_standard(a, my_prec);
                } else {
                    const char* op = " ?? ";
                    if (a->type == EXPR_SYMBOL) {
                        const char* s = a->data.symbol;
                        if      (s == SYM_Equal)        op = " == ";
                        else if (s == SYM_Less)         op = " < ";
                        else if (s == SYM_Greater)      op = " > ";
                        else if (s == SYM_LessEqual)    op = " <= ";
                        else if (s == SYM_GreaterEqual) op = " >= ";
                    }
                    printf("%s", op);
                }
            }
        }
        else if ((head == SYM_Equal || head == SYM_Unequal || head == SYM_Less || head == SYM_Greater || head == SYM_LessEqual || head == SYM_GreaterEqual || head == SYM_SameQ || head == SYM_UnsameQ || head == SYM_Set || head == SYM_SetDelayed || head == SYM_Rule || head == SYM_RuleDelayed || head == SYM_Condition || head == SYM_And || head == SYM_Or || head == SYM_Alternatives) && e->data.function.arg_count >= 2) {
            const char* op = "";
            if (head == SYM_Equal) op = " == ";
            else if (head == SYM_Unequal) op = " != ";
            else if (head == SYM_Less) op = " < ";
            else if (head == SYM_Greater) op = " > ";
            else if (head == SYM_LessEqual) op = " <= ";
            else if (head == SYM_GreaterEqual) op = " >= ";
            else if (head == SYM_SameQ) op = " === ";
            else if (head == SYM_UnsameQ) op = " =!= ";
            else if (head == SYM_Set) op = " = ";
            else if (head == SYM_SetDelayed) op = " := ";
            else if (head == SYM_Rule) op = " -> ";
            else if (head == SYM_RuleDelayed) op = " :> ";
            else if (head == SYM_Condition) op = " /; ";
            else if (head == SYM_And) op = " && ";
            else if (head == SYM_Or) op = " || ";
            else if (head == SYM_Alternatives) op = " | ";

            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (i > 0) printf("%s", op);
                print_standard(e->data.function.args[i], my_prec);
            }
        }
        else if (head == SYM_Plus) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                Expr* arg = e->data.function.args[i];
                bool is_negative = false;
                Expr* to_print = arg;
                Expr* t_copy = NULL;

                if (arg->type == EXPR_INTEGER && arg->data.integer < 0) {
                    is_negative = true;
                    t_copy = expr_new_integer(-arg->data.integer);
                    to_print = t_copy;
                } else if (arg->type == EXPR_BIGINT && mpz_sgn(arg->data.bigint) < 0) {
                    is_negative = true;
                    t_copy = negate_bigint(arg);
                    to_print = t_copy;
                } else if (arg->type == EXPR_REAL && arg->data.real < 0) {
                    is_negative = true;
                    t_copy = expr_new_real(-arg->data.real);
                    to_print = t_copy;
                } else if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL) {
                    const char* h = arg->data.function.head->data.symbol;
                    if (h == SYM_Times && arg->data.function.arg_count > 0) {
                        Expr* f_arg = arg->data.function.args[0];
                        if (f_arg->type == EXPR_INTEGER && f_arg->data.integer < 0) {
                            is_negative = true;
                            t_copy = expr_unshare(expr_copy(arg));
                            expr_free(t_copy->data.function.args[0]);
                            t_copy->data.function.args[0] = expr_new_integer(-f_arg->data.integer);

                            if (t_copy->data.function.args[0]->data.integer == 1 && t_copy->data.function.arg_count > 1) {
                                Expr** new_args = malloc(sizeof(Expr*) * (t_copy->data.function.arg_count - 1));
                                for (size_t k = 1; k < t_copy->data.function.arg_count; k++) {
                                    new_args[k-1] = expr_copy(t_copy->data.function.args[k]);
                                }
                                Expr* new_times = expr_new_function(expr_new_symbol("Times"), new_args, t_copy->data.function.arg_count - 1);
                                expr_free(t_copy);
                                t_copy = new_times;
                                free(new_args);
                            }
                            to_print = t_copy;
                        } else if (f_arg->type == EXPR_BIGINT && mpz_sgn(f_arg->data.bigint) < 0) {
                            is_negative = true;
                            t_copy = expr_unshare(expr_copy(arg));
                            expr_free(t_copy->data.function.args[0]);
                            t_copy->data.function.args[0] = negate_bigint(f_arg);

                            if (t_copy->data.function.args[0]->type == EXPR_INTEGER &&
                                t_copy->data.function.args[0]->data.integer == 1 &&
                                t_copy->data.function.arg_count > 1) {
                                Expr** new_args = malloc(sizeof(Expr*) * (t_copy->data.function.arg_count - 1));
                                for (size_t k = 1; k < t_copy->data.function.arg_count; k++) {
                                    new_args[k-1] = expr_copy(t_copy->data.function.args[k]);
                                }
                                Expr* new_times = expr_new_function(expr_new_symbol("Times"), new_args, t_copy->data.function.arg_count - 1);
                                expr_free(t_copy);
                                t_copy = new_times;
                                free(new_args);
                            }
                            to_print = t_copy;
                        } else if (f_arg->type == EXPR_REAL && f_arg->data.real < 0) {
                            is_negative = true;
                            t_copy = expr_unshare(expr_copy(arg));
                            expr_free(t_copy->data.function.args[0]);
                            t_copy->data.function.args[0] = expr_new_real(-f_arg->data.real);
                            to_print = t_copy;
                        } else if (f_arg->type == EXPR_FUNCTION && f_arg->data.function.head->type == EXPR_SYMBOL && f_arg->data.function.head->data.symbol == SYM_Rational) {
                            Expr* num = f_arg->data.function.args[0];
                            if (num->type == EXPR_INTEGER && num->data.integer < 0) {
                                is_negative = true;
                                t_copy = expr_unshare(expr_copy(arg));
                                /* Need to mutate the rational nested inside
                                 * t_copy.args[0]: unshare that level too. */
                                Expr* rat = expr_unshare(t_copy->data.function.args[0]);
                                t_copy->data.function.args[0] = rat;
                                expr_free(rat->data.function.args[0]);
                                rat->data.function.args[0] = expr_new_integer(-num->data.integer);
                                to_print = t_copy;
                            }
                        } else if (f_arg->type == EXPR_FUNCTION &&
                                   f_arg->data.function.head->type == EXPR_SYMBOL &&
                                   f_arg->data.function.head->data.symbol == SYM_Complex &&
                                   f_arg->data.function.arg_count == 2 &&
                                   f_arg->data.function.args[0]->type == EXPR_INTEGER &&
                                   f_arg->data.function.args[0]->data.integer == 0 &&
                                   f_arg->data.function.args[1]->type == EXPR_INTEGER &&
                                   f_arg->data.function.args[1]->data.integer < 0) {
                            /* Plus[..., Times[Complex[0, -n], rest...]]:
                             * print "- Times[Complex[0, n], rest...]" so the
                             * sign attaches to the operator instead of the
                             * factor. */
                            is_negative = true;
                            t_copy = expr_unshare(expr_copy(arg));
                            Expr* cplx = expr_unshare(t_copy->data.function.args[0]);
                            t_copy->data.function.args[0] = cplx;
                            int64_t im_val = cplx->data.function.args[1]->data.integer;
                            expr_free(cplx->data.function.args[1]);
                            cplx->data.function.args[1] = expr_new_integer(-im_val);
                            to_print = t_copy;
                        }
                    } else if (h == SYM_Rational && arg->data.function.arg_count == 2) {
                        Expr* num = arg->data.function.args[0];
                        if (num->type == EXPR_INTEGER && num->data.integer < 0) {
                            is_negative = true;
                            t_copy = expr_unshare(expr_copy(arg));
                            expr_free(t_copy->data.function.args[0]);
                            t_copy->data.function.args[0] = expr_new_integer(-num->data.integer);
                            to_print = t_copy;
                        }
                    } else if (h == SYM_Complex && arg->data.function.arg_count == 2 &&
                               arg->data.function.args[0]->type == EXPR_INTEGER &&
                               arg->data.function.args[0]->data.integer == 0 &&
                               arg->data.function.args[1]->type == EXPR_INTEGER &&
                               arg->data.function.args[1]->data.integer < 0) {
                        /* Plus[..., Complex[0, -n]] -> "- Complex[0, n]" so the
                         * sign sits with the operator. */
                        is_negative = true;
                        int64_t im_val = arg->data.function.args[1]->data.integer;
                        Expr* cargs[2] = { expr_new_integer(0), expr_new_integer(-im_val) };
                        t_copy = expr_new_function(expr_new_symbol("Complex"), cargs, 2);
                        to_print = t_copy;
                    }
                }

                if (i > 0) {
                    if (is_negative) {
                        printf(" - ");
                    } else {
                        printf(" + ");
                    }
                } else {
                    if (is_negative) {
                        printf("-");
                    }
                }

                print_standard(to_print, 310);
                if (t_copy) expr_free(t_copy);
            }
        }
        else if (head == SYM_List) {
            printf("{");
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (i > 0) printf(", ");
                print_standard(e->data.function.args[i], 0);
            }
            printf("}");
        }
        else if (head == SYM_Blank && e->data.function.arg_count <= 1) {
            printf("_");
            if (e->data.function.arg_count == 1) print_standard(e->data.function.args[0], 1000);
        }
        else if (head == SYM_BlankSequence && e->data.function.arg_count <= 1) {
            printf("__");
            if (e->data.function.arg_count == 1) print_standard(e->data.function.args[0], 1000);
        }
        else if (head == SYM_BlankNullSequence && e->data.function.arg_count <= 1) {
            printf("___");
            if (e->data.function.arg_count == 1) print_standard(e->data.function.args[0], 1000);
        }
        else if (head == SYM_Slot && e->data.function.arg_count == 1) {
            /* Slot[n] -> #n (e.g. #1). Atom-level precedence already returned
             * by get_expr_prec via the fall-through, so no extra parens. */
            Expr* idx = e->data.function.args[0];
            printf("#");
            if (idx->type == EXPR_INTEGER) {
                printf("%" PRId64, idx->data.integer);
            } else if (idx->type == EXPR_STRING) {
                printf("\"%s\"", idx->data.string);
            } else {
                print_standard(idx, 1000);
            }
        }
        else if (head == SYM_SlotSequence && e->data.function.arg_count == 1) {
            /* SlotSequence[n] -> ##n; SlotSequence[1] is conventionally "##". */
            Expr* idx = e->data.function.args[0];
            printf("##");
            if (idx->type == EXPR_INTEGER) {
                if (idx->data.integer != 1) printf("%" PRId64, idx->data.integer);
            } else {
                print_standard(idx, 1000);
            }
        }
        else if (head == SYM_Function && e->data.function.arg_count == 1) {
            /* Function[body] -> body &. The body normally has high precedence
             * (Plus/Times/Power etc.) so it prints without extra parens; if
             * something genuinely lower-precedence shows up, print_standard
             * inserts them at parent_prec=90. */
            print_standard(e->data.function.args[0], 90);
            printf(" &");
        }
        else if (head == SYM_Pattern && e->data.function.arg_count == 2 &&
                 e->data.function.args[0]->type == EXPR_SYMBOL &&
                 e->data.function.args[1]->type == EXPR_FUNCTION &&
                 e->data.function.args[1]->data.function.head->type == EXPR_SYMBOL &&
                 (e->data.function.args[1]->data.function.head->data.symbol == SYM_Blank ||
                  e->data.function.args[1]->data.function.head->data.symbol == SYM_BlankSequence ||
                  e->data.function.args[1]->data.function.head->data.symbol == SYM_BlankNullSequence) &&
                 e->data.function.args[1]->data.function.arg_count <= 1) {
            /* Pattern[name, Blank[...]] -> name_ / name_h / name__ / name___ */
            printf("%s", e->data.function.args[0]->data.symbol);
            print_standard(e->data.function.args[1], 1000);
        }
        else {
            print_standard(e->data.function.head, 1000);
            printf("[");
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (i > 0) printf(", ");
                print_standard(e->data.function.args[i], 0);
            }
            printf("]");
        }
    } else {
        print_standard(e->data.function.head, 1000);
        printf("[");
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (i > 0) printf(", ");
            print_standard(e->data.function.args[i], 0);
        }    
        printf("]");
    }

    if (need_parens) printf(")");
}

void expr_print(Expr* e) {
    print_standard(e, 0);
    printf("\n");
}

char* expr_to_string(Expr* e) {
    // A bit of a hack: print to a memory buffer
    char* buffer = NULL;
    size_t len;
    FILE* stream = open_memstream(&buffer, &len);
    if (!stream) return NULL;

    // Temporarily redirect stdout
    FILE* old_stdout = stdout;
    stdout = stream;
    
    print_standard(e, 0);
    
    // Restore stdout and close stream
    stdout = old_stdout;
    fclose(stream);
    
    return buffer;
}

char* expr_to_string_fullform(Expr* e) {
    char* buffer = NULL;
    size_t len;
    FILE* stream = open_memstream(&buffer, &len);
    if (!stream) return NULL;

    FILE* old_stdout = stdout;
    stdout = stream;
    
    expr_print_fullform(e);
    
    stdout = old_stdout;
    fclose(stream);
    
    return buffer;
}

#include "symtab.h"

Expr* builtin_print(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        print_standard(res->data.function.args[i], 0);
    }
    printf("\n");
    return expr_new_symbol("Null");
}

Expr* builtin_fullform(Expr* res) {
    (void)res;
    return NULL; // Unevaluated wrapper
}

Expr* builtin_inputform(Expr* res) {
    (void)res;
    return NULL; // Unevaluated wrapper
}

Expr* builtin_texform(Expr* res) {
    (void)res;
    return NULL; // Unevaluated wrapper; rendering happens in print_standard
}

/*
 * Pretty-prints SeriesData[x, x0, {a0, ..., a_{k-1}}, nmin, nmax, den] as
 *     a0 + a1 (x-x0) + ... + a_{k-1} (x-x0)^((nmin+k-1)/den) + O[x-x0]^(nmax/den).
 * The rendering is built by constructing an equivalent Plus-expression of
 * Times/Power nodes and delegating to print_standard so that sign handling,
 * precedence, and fractional-exponent formatting are shared with ordinary
 * expressions. If the arguments are not in a shape we can render (for
 * example non-literal nmin/nmax/den or a non-List coefficient container),
 * we fall back to the generic head[args...] form.
 */

static int64_t abs_i64(int64_t v) { return v < 0 ? -v : v; }

/* Builds (x - x0), or just x when x0 is a literal zero. The result is a
 * fresh expression owned by the caller. Common negative-number cases for
 * x0 are normalised (e.g. x - (-2) becomes x + 2) so that the Plus printer
 * produces a clean subtraction. */
static Expr* series_build_base(Expr* x, Expr* x0) {
    if (x0->type == EXPR_INTEGER && x0->data.integer == 0) return expr_copy(x);
    if (x0->type == EXPR_REAL && x0->data.real == 0.0) return expr_copy(x);

    if (x0->type == EXPR_INTEGER && x0->data.integer < 0) {
        Expr* args[2] = { expr_copy(x), expr_new_integer(-x0->data.integer) };
        return expr_new_function(expr_new_symbol("Plus"), args, 2);
    }
    if (x0->type == EXPR_REAL && x0->data.real < 0.0) {
        Expr* args[2] = { expr_copy(x), expr_new_real(-x0->data.real) };
        return expr_new_function(expr_new_symbol("Plus"), args, 2);
    }

    Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(x0) };
    Expr* neg = expr_new_function(expr_new_symbol("Times"), neg_args, 2);
    Expr* plus_args[2] = { expr_copy(x), neg };
    return expr_new_function(expr_new_symbol("Plus"), plus_args, 2);
}

/* Builds base^(num/den) as an Expr, reducing the exponent to lowest terms.
 * Returns NULL to signal "exponent is zero" (in which case the caller
 * should use the coefficient alone) and a fresh base-copy when the
 * reduced exponent is 1. */
static Expr* series_build_power(Expr* base, int64_t num, int64_t den) {
    if (num == 0) return NULL;
    int64_t g = gcd(abs_i64(num), abs_i64(den));
    if (g == 0) g = 1;
    int64_t n = num / g;
    int64_t d = den / g;
    if (d < 0) { d = -d; n = -n; }

    if (d == 1 && n == 1) return expr_copy(base);

    Expr* exp;
    if (d == 1) exp = expr_new_integer(n);
    else exp = make_rational(n, d);

    Expr* pargs[2] = { expr_copy(base), exp };
    return expr_new_function(expr_new_symbol("Power"), pargs, 2);
}

/* Combines a coefficient with an optional power term into a single summand.
 * Returns NULL when the coefficient is a literal zero (i.e. the term is
 * omitted from the displayed sum). A NULL power corresponds to (x-x0)^0
 * and yields just the coefficient. */
static Expr* series_build_term(Expr* coef, Expr* power /* borrowed, may be NULL */) {
    if (coef->type == EXPR_INTEGER && coef->data.integer == 0) return NULL;
    if (coef->type == EXPR_REAL && coef->data.real == 0.0) return NULL;

    if (power == NULL) return expr_copy(coef);

    bool coef_is_one = (coef->type == EXPR_INTEGER && coef->data.integer == 1);
    if (coef_is_one) {
        /* A bare Plus would merge visually with the surrounding series
         * Plus (equal precedence -> no parens), so wrap it in a 1-arg
         * Times to force the Times printer to add parentheses. */
        if (power->type == EXPR_FUNCTION &&
            power->data.function.head->type == EXPR_SYMBOL &&
            power->data.function.head->data.symbol == SYM_Plus) {
            Expr* targs[1] = { expr_copy(power) };
            return expr_new_function(expr_new_symbol("Times"), targs, 1);
        }
        return expr_copy(power);
    }

    Expr* args[2] = { expr_copy(coef), expr_copy(power) };
    return expr_new_function(expr_new_symbol("Times"), args, 2);
}

/* Renders e as the generic head[arg1, arg2, ...] form. Used as a fallback
 * when SeriesData has unexpected argument shapes. */
static void print_series_generic(Expr* e, int parent_prec) {
    (void)parent_prec;
    print_standard(e->data.function.head, 1000);
    printf("[");
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (i > 0) printf(", ");
        print_standard(e->data.function.args[i], 0);
    }
    printf("]");
}

static void print_series_data(Expr* e, int parent_prec) {
    Expr* x        = e->data.function.args[0];
    Expr* x0       = e->data.function.args[1];
    Expr* coefs    = e->data.function.args[2];
    Expr* nmin_e   = e->data.function.args[3];
    Expr* nmax_e   = e->data.function.args[4];
    Expr* den_e    = e->data.function.args[5];

    bool coefs_ok = (coefs->type == EXPR_FUNCTION &&
                     coefs->data.function.head->type == EXPR_SYMBOL &&
                     coefs->data.function.head->data.symbol == SYM_List);
    bool ints_ok  = (nmin_e->type == EXPR_INTEGER &&
                     nmax_e->type == EXPR_INTEGER &&
                     den_e->type  == EXPR_INTEGER &&
                     den_e->data.integer != 0);
    if (!coefs_ok || !ints_ok) {
        print_series_generic(e, parent_prec);
        return;
    }

    int64_t nmin = nmin_e->data.integer;
    int64_t nmax = nmax_e->data.integer;
    int64_t den  = den_e->data.integer;
    size_t  k    = coefs->data.function.arg_count;

    Expr* base = series_build_base(x, x0);

    /* Collect non-zero summands plus the trailing O-term. */
    Expr** terms = malloc(sizeof(Expr*) * (k + 1));
    size_t tcount = 0;

    for (size_t i = 0; i < k; i++) {
        int64_t num = nmin + (int64_t)i;
        Expr* power = series_build_power(base, num, den);
        Expr* term  = series_build_term(coefs->data.function.args[i], power);
        if (power) expr_free(power);
        if (term) terms[tcount++] = term;
    }

    /* Trailing O[base]^(nmax/den) term. Always emitted, even when all
     * coefficients are zero, because it conveys the truncation order. */
    {
        Expr* oargs[1] = { expr_copy(base) };
        Expr* o_call   = expr_new_function(expr_new_symbol("O"), oargs, 1);
        int64_t g  = gcd(abs_i64(nmax), abs_i64(den));
        if (g == 0) g = 1;
        int64_t n  = nmax / g;
        int64_t d  = den  / g;
        if (d < 0) { d = -d; n = -n; }
        Expr* exp_expr = (d == 1) ? expr_new_integer(n) : make_rational(n, d);
        Expr* pargs[2] = { o_call, exp_expr };
        terms[tcount++] = expr_new_function(expr_new_symbol("Power"), pargs, 2);
    }

    Expr* sum;
    if (tcount == 1) {
        sum = terms[0];
    } else {
        sum = expr_new_function(expr_new_symbol("Plus"), terms, tcount);
    }
    free(terms);

    print_standard(sum, parent_prec);
    expr_free(sum);
    expr_free(base);
}

/* ===================== TeXForm renderer =====================
 *
 * Renders an Expr as AMS-LaTeX-compatible TeX. TeXForm[expr] in
 * print_standard dispatches here. Precedence uses the same numeric
 * levels as get_expr_prec so parens are inserted consistently.
 * Single-character symbol names are emitted bare (italicised by TeX
 * default math mode); multi-character names are wrapped in \text{}.
 */

/* Emit a number in TeX form (integers, reals, bigints, mpfr). */
static void print_tex_number(Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER: printf("%" PRId64, e->data.integer); break;
        case EXPR_REAL: {
            char buf[64];
            snprintf(buf, sizeof(buf), "%g", e->data.real);
            if (strchr(buf, '.') == NULL && strchr(buf, 'e') == NULL && strchr(buf, 'E') == NULL) printf("%s.", buf);
            else printf("%s", buf);
            break;
        }
        case EXPR_BIGINT: {
            char* s = mpz_get_str(NULL, 10, e->data.bigint);
            printf("%s", s);
            free(s);
            break;
        }
#ifdef USE_MPFR
        case EXPR_MPFR: {
            long digits = (long)ceil((double)mpfr_get_prec(e->data.mpfr) / 3.3219280948873626);
            if (digits < 1) digits = 1;
            char* buf = NULL;
            int len = mpfr_asprintf(&buf, "%.*Rg", (int)digits, e->data.mpfr);
            if (len >= 0 && buf) { printf("%s", buf); mpfr_free_str(buf); }
            break;
        }
#endif
        default: expr_print_fullform(e); break;
    }
}

static void print_tex_symbol(const char* raw_name) {
    const char* name = context_display_name(raw_name);
    if (strcmp(name, "Pi") == 0)              { printf("\\pi "); return; }
    if (strcmp(name, "E") == 0)               { printf("e"); return; }
    if (strcmp(name, "I") == 0)               { printf("i"); return; }
    if (strcmp(name, "Infinity") == 0)        { printf("\\infty "); return; }
    if (strcmp(name, "ComplexInfinity") == 0) { printf("\\tilde{\\infty }"); return; }
    if (strcmp(name, "Indeterminate") == 0)   { printf("\\text{Indeterminate}"); return; }
    if (strcmp(name, "EulerGamma") == 0)      { printf("\\gamma "); return; }
    if (strcmp(name, "GoldenRatio") == 0)     { printf("\\phi "); return; }
    if (strcmp(name, "Catalan") == 0)         { printf("C"); return; }
    if (strcmp(name, "Degree") == 0)          { printf("{}^{\\circ }"); return; }
    if (strcmp(name, "True") == 0)            { printf("\\text{True}"); return; }
    if (strcmp(name, "False") == 0)           { printf("\\text{False}"); return; }
    if (strcmp(name, "Null") == 0)            { printf("\\text{Null}"); return; }
    size_t len = strlen(name);
    if (len == 1) printf("%s", name);
    else printf("\\text{%s}", name);
}

/* Looks up a trig/hyperbolic head and returns its TeX command and
 * whether it's an inverse variant. Returns false if the head is not a
 * standard trig/hyperbolic function. */
static bool tex_trig_lookup(const char* h, const char** tex, bool* inv) {
    static const struct { const char* head; const char* tex; bool inv; } table[] = {
        {"Sin",  "\\sin",  false}, {"Cos",  "\\cos",  false},
        {"Tan",  "\\tan",  false}, {"Cot",  "\\cot",  false},
        {"Sec",  "\\sec",  false}, {"Csc",  "\\csc",  false},
        {"Sinh", "\\sinh", false}, {"Cosh", "\\cosh", false},
        {"Tanh", "\\tanh", false}, {"Coth", "\\coth", false},
        {"Sech", "\\text{sech}", false}, {"Csch", "\\text{csch}", false},
        {"ArcSin", "\\sin", true}, {"ArcCos", "\\cos", true},
        {"ArcTan", "\\tan", true}, {"ArcCot", "\\cot", true},
        {"ArcSec", "\\sec", true}, {"ArcCsc", "\\csc", true},
        {"ArcSinh","\\sinh", true}, {"ArcCosh","\\cosh", true},
        {"ArcTanh","\\tanh", true}, {"ArcCoth","\\coth", true},
        {"ArcSech","\\text{sech}", true}, {"ArcCsch","\\text{csch}", true},
    };
    size_t n = sizeof(table) / sizeof(table[0]);
    for (size_t i = 0; i < n; i++) {
        if (strcmp(h, table[i].head) == 0) {
            *tex = table[i].tex;
            *inv = table[i].inv;
            return true;
        }
    }
    return false;
}

/* Maps a named head to a TeX symbol when the function is to be rendered
 * like f(x). Returns NULL if the head should use its default name. */
static const char* tex_function_name(const char* h) {
    if (strcmp(h, "Log") == 0)      return "\\log ";
    if (strcmp(h, "Exp") == 0)      return "\\exp ";
    if (strcmp(h, "Sqrt") == 0)     return NULL; /* handled via \sqrt{} */
    if (strcmp(h, "Abs") == 0)      return NULL; /* handled via |...|  */
    if (strcmp(h, "Floor") == 0)    return NULL; /* handled via \lfloor */
    if (strcmp(h, "Ceiling") == 0)  return NULL; /* handled via \lceil  */
    if (strcmp(h, "Gamma") == 0)    return "\\Gamma ";
    if (strcmp(h, "Re") == 0)       return "\\Re ";
    if (strcmp(h, "Im") == 0)       return "\\Im ";
    if (strcmp(h, "Conjugate") == 0)return NULL; /* handled via \overline */
    return NULL;
}

/* Emits a comma-separated argument list wrapped in (...). */
static void print_tex_args_parens(Expr* e, size_t start) {
    printf("(");
    for (size_t i = start; i < e->data.function.arg_count; i++) {
        if (i > start) printf(",");
        print_tex(e->data.function.args[i], 0);
    }
    printf(")");
}

/* Detects whether `im` is a negative numeric, and if so produces a
 * freshly-owned Expr holding its absolute value (caller must free). */
static bool tex_extract_negative(Expr* im, Expr** abs_out) {
    int64_t rn, rd;
    if (im->type == EXPR_INTEGER && im->data.integer < 0) {
        *abs_out = expr_new_integer(-im->data.integer);
        return true;
    }
    if (im->type == EXPR_REAL && im->data.real < 0.0) {
        *abs_out = expr_new_real(-im->data.real);
        return true;
    }
    if (is_rational(im, &rn, &rd) && rn < 0) {
        *abs_out = make_rational(-rn, rd);
        return true;
    }
    return false;
}

static void print_tex(Expr* e, int parent_prec) {
    if (!e) { printf("\\text{Null}"); return; }

    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT
#ifdef USE_MPFR
        || e->type == EXPR_MPFR
#endif
    ) {
        print_tex_number(e);
        return;
    }
    if (e->type == EXPR_SYMBOL) { print_tex_symbol(e->data.symbol); return; }
    if (e->type == EXPR_STRING) { printf("\\text{\"%s\"}", e->data.string); return; }
    if (e->type != EXPR_FUNCTION) return;

    int my_prec = get_expr_prec(e);
    bool need_parens = (my_prec < parent_prec);
    if (need_parens) printf("\\left(");

    Expr* head_expr = e->data.function.head;
    if (head_expr->type == EXPR_SYMBOL) {
        const char* head = head_expr->data.symbol;
        size_t argc = e->data.function.arg_count;

        if (head == SYM_Rational && argc == 2) {
            printf("\\frac{");
            print_tex(e->data.function.args[0], 0);
            printf("}{");
            print_tex(e->data.function.args[1], 0);
            printf("}");
        }
        else if (head == SYM_Complex && argc == 2) {
            Expr* re = e->data.function.args[0];
            Expr* im = e->data.function.args[1];
            bool re_zero = (re->type == EXPR_INTEGER && re->data.integer == 0);
            bool im_one = (im->type == EXPR_INTEGER && im->data.integer == 1);
            bool im_minus_one = (im->type == EXPR_INTEGER && im->data.integer == -1);
            Expr* im_abs = NULL;
            bool im_neg = im_minus_one;
            if (!im_neg) im_neg = tex_extract_negative(im, &im_abs);
            if (!re_zero) {
                print_tex(re, 310);
                printf(im_neg ? "-" : "+");
            } else if (im_neg) {
                printf("-");
            }
            if (im_one || im_minus_one) printf("i");
            else if (im_abs) { print_tex(im_abs, 400); printf(" i"); expr_free(im_abs); }
            else { print_tex(im, 400); printf(" i"); }
        }
        else if (head == SYM_Sqrt && argc == 1) {
            printf("\\sqrt{");
            print_tex(e->data.function.args[0], 0);
            printf("}");
        }
        else if (head == SYM_Abs && argc == 1) {
            printf("\\left| ");
            print_tex(e->data.function.args[0], 0);
            printf("\\right| ");
        }
        else if (head == SYM_Floor && argc == 1) {
            printf("\\left\\lfloor ");
            print_tex(e->data.function.args[0], 0);
            printf("\\right\\rfloor ");
        }
        else if (head == SYM_Ceiling && argc == 1) {
            printf("\\left\\lceil ");
            print_tex(e->data.function.args[0], 0);
            printf("\\right\\rceil ");
        }
        else if (head == SYM_Conjugate && argc == 1) {
            printf("\\overline{");
            print_tex(e->data.function.args[0], 0);
            printf("}");
        }
        else if (head == SYM_Power && argc == 2) {
            Expr* base = e->data.function.args[0];
            Expr* exp = e->data.function.args[1];
            int64_t n, d;
            /* Negative exponent -> 1/base^|exp| with \frac */
            bool exp_neg = false;
            Expr* pos_exp = NULL;
            if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
                exp_neg = true; pos_exp = expr_new_integer(-exp->data.integer);
            } else if (is_rational(exp, &n, &d) && n < 0) {
                exp_neg = true; pos_exp = make_rational(-n, d);
            } else if (exp->type == EXPR_REAL && exp->data.real < 0.0) {
                exp_neg = true; pos_exp = expr_new_real(-exp->data.real);
            }
            if (exp_neg) {
                printf("\\frac{1}{");
                if (pos_exp->type == EXPR_INTEGER && pos_exp->data.integer == 1) {
                    print_tex(base, 0);
                } else {
                    Expr* pargs[2] = { expr_copy(base), expr_copy(pos_exp) };
                    Expr* tmp = expr_new_function(expr_new_symbol("Power"), pargs, 2);
                    print_tex(tmp, 0);
                    expr_free(tmp);
                }
                printf("}");
                expr_free(pos_exp);
            } else if (is_rational(exp, &n, &d) && n == 1 && d == 2) {
                printf("\\sqrt{");
                print_tex(base, 0);
                printf("}");
            } else if (is_rational(exp, &n, &d) && n == 1 && d > 2) {
                printf("\\sqrt[%" PRId64 "]{", d);
                print_tex(base, 0);
                printf("}");
            } else {
                /* Negative numeric atoms parse as unary-minus, so render
                 * `-1^{1/3}` as `(-1)^{1/3}` to keep the meaning unambiguous. */
                bool base_paren = false;
                int64_t bn, bd;
                if (base->type == EXPR_INTEGER && base->data.integer < 0) base_paren = true;
                else if (base->type == EXPR_REAL && base->data.real < 0.0) base_paren = true;
                else if (base->type == EXPR_BIGINT && mpz_sgn(base->data.bigint) < 0) base_paren = true;
                else if (is_rational(base, &bn, &bd) && bn < 0) base_paren = true;
                if (base_paren) {
                    printf("\\left(");
                    print_tex(base, 0);
                    printf("\\right)");
                } else {
                    print_tex(base, 600); /* force parens for compound bases */
                }
                printf("^{");
                print_tex(exp, 0);
                printf("}");
            }
        }
        else if (head == SYM_Times) {
            size_t count = argc;
            int64_t lead_start = 0;
            bool flipped_sign = false;
            Expr* flipped_head = NULL;
            if (count > 0) {
                Expr* a0 = e->data.function.args[0];
                if (a0->type == EXPR_INTEGER && a0->data.integer < 0) {
                    flipped_sign = true;
                    if (a0->data.integer == -1 && count > 1) lead_start = 1;
                    else flipped_head = expr_new_integer(-a0->data.integer);
                } else if (a0->type == EXPR_REAL && a0->data.real < 0.0) {
                    flipped_sign = true;
                    if (a0->data.real == -1.0 && count > 1) lead_start = 1;
                    else flipped_head = expr_new_real(-a0->data.real);
                } else {
                    int64_t rn, rd;
                    if (is_rational(a0, &rn, &rd) && rn < 0) {
                        flipped_sign = true;
                        flipped_head = make_rational(-rn, rd);
                    }
                }
            }
            if (flipped_sign) printf("-");

            /* Collect numerator / denominator via negative-exponent Powers,
             * matching print_standard so a/b renders as \frac{a}{b}. */
            Expr** num_args = malloc(sizeof(Expr*) * count);
            Expr** den_args = malloc(sizeof(Expr*) * count);
            size_t num_count = 0, den_count = 0;
            for (size_t i = (size_t)lead_start; i < count; i++) {
                Expr* arg = (i == 0 && flipped_head) ? flipped_head : e->data.function.args[i];
                bool den_ok = false;
                if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL
                    && arg->data.function.head->data.symbol == SYM_Power
                    && arg->data.function.arg_count == 2) {
                    Expr* exp = arg->data.function.args[1];
                    int64_t n, d;
                    Expr* pos_exp = NULL;
                    bool neg = false;
                    if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
                        neg = true; pos_exp = expr_new_integer(-exp->data.integer);
                    } else if (is_rational(exp, &n, &d) && n < 0) {
                        neg = true; pos_exp = make_rational(-n, d);
                    } else if (exp->type == EXPR_REAL && exp->data.real < 0.0) {
                        neg = true; pos_exp = expr_new_real(-exp->data.real);
                    }
                    if (neg) {
                        if (pos_exp->type == EXPR_INTEGER && pos_exp->data.integer == 1) {
                            den_args[den_count++] = expr_copy(arg->data.function.args[0]);
                            expr_free(pos_exp);
                        } else {
                            Expr* p_args[2] = { expr_copy(arg->data.function.args[0]), pos_exp };
                            den_args[den_count++] = expr_new_function(expr_new_symbol("Power"), p_args, 2);
                        }
                        den_ok = true;
                    }
                }
                if (!den_ok) num_args[num_count++] = expr_copy(arg);
            }

            if (den_count == 0) {
                for (size_t i = 0; i < num_count; i++) {
                    if (i > 0) printf(" ");
                    print_tex(num_args[i], 400);
                }
            } else {
                printf("\\frac{");
                if (num_count == 0) { printf("1"); }
                else {
                    for (size_t i = 0; i < num_count; i++) {
                        if (i > 0) printf(" ");
                        print_tex(num_args[i], 0);
                    }
                }
                printf("}{");
                for (size_t i = 0; i < den_count; i++) {
                    if (i > 0) printf(" ");
                    print_tex(den_args[i], 0);
                }
                printf("}");
            }
            for (size_t i = 0; i < num_count; i++) expr_free(num_args[i]);
            for (size_t i = 0; i < den_count; i++) expr_free(den_args[i]);
            free(num_args);
            free(den_args);
            if (flipped_head) expr_free(flipped_head);
        }
        else if (head == SYM_Plus) {
            for (size_t i = 0; i < argc; i++) {
                Expr* arg = e->data.function.args[i];
                bool is_neg = false;
                Expr* to_print = arg;
                Expr* owned = NULL;

                if (arg->type == EXPR_INTEGER && arg->data.integer < 0) {
                    is_neg = true; owned = expr_new_integer(-arg->data.integer); to_print = owned;
                } else if (arg->type == EXPR_BIGINT && mpz_sgn(arg->data.bigint) < 0) {
                    is_neg = true; owned = negate_bigint(arg); to_print = owned;
                } else if (arg->type == EXPR_REAL && arg->data.real < 0.0) {
                    is_neg = true; owned = expr_new_real(-arg->data.real); to_print = owned;
                } else if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL
                           && arg->data.function.head->data.symbol == SYM_Times
                           && arg->data.function.arg_count > 0) {
                    Expr* f0 = arg->data.function.args[0];
                    if (f0->type == EXPR_INTEGER && f0->data.integer < 0) {
                        is_neg = true;
                        owned = expr_unshare(expr_copy(arg));
                        expr_free(owned->data.function.args[0]);
                        owned->data.function.args[0] = expr_new_integer(-f0->data.integer);
                        if (owned->data.function.args[0]->data.integer == 1 && owned->data.function.arg_count > 1) {
                            Expr** na = malloc(sizeof(Expr*) * (owned->data.function.arg_count - 1));
                            for (size_t k = 1; k < owned->data.function.arg_count; k++) na[k-1] = expr_copy(owned->data.function.args[k]);
                            Expr* t = expr_new_function(expr_new_symbol("Times"), na, owned->data.function.arg_count - 1);
                            expr_free(owned); owned = t; free(na);
                        }
                        to_print = owned;
                    } else if (f0->type == EXPR_BIGINT && mpz_sgn(f0->data.bigint) < 0) {
                        is_neg = true;
                        owned = expr_unshare(expr_copy(arg));
                        expr_free(owned->data.function.args[0]);
                        owned->data.function.args[0] = negate_bigint(f0);
                        if (owned->data.function.args[0]->type == EXPR_INTEGER &&
                            owned->data.function.args[0]->data.integer == 1 &&
                            owned->data.function.arg_count > 1) {
                            Expr** na = malloc(sizeof(Expr*) * (owned->data.function.arg_count - 1));
                            for (size_t k = 1; k < owned->data.function.arg_count; k++) na[k-1] = expr_copy(owned->data.function.args[k]);
                            Expr* t = expr_new_function(expr_new_symbol("Times"), na, owned->data.function.arg_count - 1);
                            expr_free(owned); owned = t; free(na);
                        }
                        to_print = owned;
                    } else if (f0->type == EXPR_REAL && f0->data.real < 0.0) {
                        is_neg = true;
                        owned = expr_unshare(expr_copy(arg));
                        expr_free(owned->data.function.args[0]);
                        owned->data.function.args[0] = expr_new_real(-f0->data.real);
                        to_print = owned;
                    } else {
                        int64_t rn, rd;
                        if (is_rational(f0, &rn, &rd) && rn < 0) {
                            is_neg = true;
                            owned = expr_unshare(expr_copy(arg));
                            Expr* new_rat = make_rational(-rn, rd);
                            expr_free(owned->data.function.args[0]);
                            owned->data.function.args[0] = new_rat;
                            to_print = owned;
                        }
                    }
                } else {
                    int64_t rn, rd;
                    if (is_rational(arg, &rn, &rd) && rn < 0) {
                        is_neg = true; owned = make_rational(-rn, rd); to_print = owned;
                    }
                }

                if (i > 0) printf(is_neg ? "-" : "+");
                else if (is_neg) printf("-");

                print_tex(to_print, 310);
                if (owned) expr_free(owned);
            }
        }
        else if (head == SYM_List) {
            printf("\\{");
            for (size_t i = 0; i < argc; i++) {
                if (i > 0) printf(",");
                print_tex(e->data.function.args[i], 0);
            }
            printf("\\}");
        }
        else if (head == SYM_Inequality && argc >= 3 && (argc & 1u) == 1) {
            for (size_t i = 0; i < argc; i++) {
                Expr* a = e->data.function.args[i];
                if ((i & 1u) == 0) {
                    print_tex(a, my_prec);
                } else {
                    const char* op = "";
                    if (a->type == EXPR_SYMBOL) {
                        const char* s = a->data.symbol;
                        if      (s == SYM_Equal)        op = "=";
                        else if (s == SYM_Less)         op = "<";
                        else if (s == SYM_Greater)      op = ">";
                        else if (s == SYM_LessEqual)    op = "\\leq ";
                        else if (s == SYM_GreaterEqual) op = "\\geq ";
                    }
                    printf("%s", op);
                }
            }
        }
        else if ((head == SYM_Equal || head == SYM_Unequal
                  || head == SYM_Less || head == SYM_Greater
                  || head == SYM_LessEqual || head == SYM_GreaterEqual
                  || head == SYM_SameQ || head == SYM_UnsameQ
                  || head == SYM_Rule || head == SYM_RuleDelayed
                  || head == SYM_Set || head == SYM_SetDelayed
                  || head == SYM_And || head == SYM_Or) && argc >= 2) {
            const char* op = "";
            if (head == SYM_Equal) op = "=";
            else if (head == SYM_Unequal) op = "\\neq ";
            else if (head == SYM_Less) op = "<";
            else if (head == SYM_Greater) op = ">";
            else if (head == SYM_LessEqual) op = "\\leq ";
            else if (head == SYM_GreaterEqual) op = "\\geq ";
            else if (head == SYM_SameQ) op = "\\equiv ";
            else if (head == SYM_UnsameQ) op = "\\not\\equiv ";
            else if (head == SYM_Rule || head == SYM_RuleDelayed) op = "\\to ";
            else if (head == SYM_Set) op = "=";
            else if (head == SYM_SetDelayed) op = ":=";
            else if (head == SYM_And) op = "\\land ";
            else if (head == SYM_Or)  op = "\\lor ";
            for (size_t i = 0; i < argc; i++) {
                if (i > 0) printf("%s", op);
                print_tex(e->data.function.args[i], my_prec);
            }
        }
        else if (head == SYM_Not && argc == 1) {
            printf("\\neg ");
            print_tex(e->data.function.args[0], 1000);
        }
        else {
            const char* tex_fn = NULL;
            bool is_inv = false;
            if (tex_trig_lookup(head, &tex_fn, &is_inv) && argc == 1) {
                /* \sin(x) or \sin ^{-1}(x) */
                printf("%s", tex_fn);
                if (is_inv) printf(" ^{-1}");
                print_tex_args_parens(e, 0);
            } else {
                const char* named = tex_function_name(head);
                if (named) {
                    printf("%s", named);
                    print_tex_args_parens(e, 0);
                } else {
                    /* Generic: f(a, b). Head gets symbol treatment so single-
                     * char heads render italic, multi-char heads use \text{}. */
                    print_tex(head_expr, 1000);
                    print_tex_args_parens(e, 0);
                }
            }
        }
    } else {
        /* Non-symbol head: render as (head)(args...). */
        print_tex(head_expr, 1000);
        print_tex_args_parens(e, 0);
    }

    if (need_parens) printf("\\right)");
}
