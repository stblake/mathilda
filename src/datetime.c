#include "datetime.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

Expr* builtin_timing(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL;
    }
    
    Expr* arg = res->data.function.args[0];
    
    clock_t start = clock();
    Expr* evaluated = evaluate(arg);
    clock_t end = clock();
    
    double time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    Expr** results = malloc(sizeof(Expr*) * 2);
    results[0] = expr_new_real(time_used);
    results[1] = evaluated;
    
    Expr* final_res = expr_new_function(expr_new_symbol("List"), results, 2);
    free(results);
    return final_res;
}

static int compare_doubles(const void* a, const void* b) {
    double arg1 = *(const double*)a;
    double arg2 = *(const double*)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

Expr* builtin_repeated_timing(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)) {
        return NULL;
    }
    
    Expr* arg = res->data.function.args[0];
    double t_target = 1.0; // Default 1 second
    
    if (res->data.function.arg_count == 2) {
        Expr* t_expr = res->data.function.args[1];
        if (t_expr->type == EXPR_INTEGER) {
            t_target = (double)t_expr->data.integer;
        } else if (t_expr->type == EXPR_REAL) {
            t_target = t_expr->data.real;
        } else {
            Expr* t_eval = evaluate(t_expr);
            if (t_eval->type == EXPR_INTEGER) {
                t_target = (double)t_eval->data.integer;
            } else if (t_eval->type == EXPR_REAL) {
                t_target = t_eval->data.real;
            } else {
                expr_free(t_eval);
                return NULL;
            }
            expr_free(t_eval);
        }
    }
    
    Expr* first_evaluated = NULL;
    
    size_t cap = 16;
    size_t count = 0;
    double* timings = malloc(sizeof(double) * cap);
    
    clock_t total_start = clock();
    double elapsed_total = 0.0;
    
    // Minimum 4 evaluations, or until t_target seconds is reached
    while (count < 4 || elapsed_total < t_target) {
        clock_t start = clock();
        Expr* evaluated = evaluate(arg);
        clock_t end = clock();
        
        if (count == 0) {
            first_evaluated = evaluated;
        } else {
            expr_free(evaluated);
        }
        
        if (count == cap) {
            cap *= 2;
            timings = realloc(timings, sizeof(double) * cap);
        }
        
        double time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        timings[count++] = time_used;
        
        clock_t current_time = clock();
        elapsed_total = ((double) (current_time - total_start)) / CLOCKS_PER_SEC;
    }
    
    qsort(timings, count, sizeof(double), compare_doubles);
    
    // Trimmed mean (discard lower and upper quartiles)
    size_t q_size = count / 4;
    size_t start_idx = q_size;
    size_t end_idx = count - q_size;
    
    double sum = 0.0;
    for (size_t i = start_idx; i < end_idx; i++) {
        sum += timings[i];
    }
    double average_time = sum / (double)(end_idx - start_idx);
    
    free(timings);
    
    Expr** results = malloc(sizeof(Expr*) * 2);
    results[0] = expr_new_real(average_time);
    results[1] = first_evaluated;
    
    Expr* final_res = expr_new_function(expr_new_symbol("List"), results, 2);
    free(results);
    return final_res;
}

/*
 * Days from the proleptic Gregorian epoch (1900-01-01) for the given
 * (y, m, d).  The Fliegel & Van Flandern Julian-day-number formula is
 * used directly; it tolerates out-of-range month and day inputs and
 * normalises them implicitly (e.g. {2022,2,31} = {2022,3,3}).  Month
 * normalisation for m <= 0 or m > 12 is handled explicitly first so
 * the JDN expression remains within the algorithm's documented range.
 */
static int64_t days_since_1900(int64_t y, int64_t m, int64_t d) {
    if (m > 12) {
        int64_t k = (m - 1) / 12;
        y += k;
        m -= 12 * k;
    } else if (m < 1) {
        int64_t k = (12 - m) / 12;
        y -= k;
        m += 12 * k;
    }
    int64_t a = (14 - m) / 12;
    int64_t yy = y + 4800 - a;
    int64_t mm = m + 12 * a - 3;
    int64_t jdn = d + (153 * mm + 2) / 5
                + 365 * yy + yy / 4 - yy / 100 + yy / 400
                - 32045;
    return jdn - 2415021;
}

static int expr_to_double_strict(const Expr* e, double* out) {
    if (!e) return 0;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer; return 1;
        case EXPR_REAL:    *out = e->data.real;            return 1;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return 1;
        default: return 0;
    }
}

Expr* builtin_absolute_time(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 0) {
        time_t now = time(NULL);
        if (now == (time_t)-1) return NULL;
        struct tm lt;
        struct tm* lp = localtime(&now);
        if (!lp) return NULL;
        lt = *lp;

        int64_t y   = (int64_t)lt.tm_year + 1900;
        int64_t mo  = (int64_t)lt.tm_mon + 1;
        int64_t d   = (int64_t)lt.tm_mday;
        int64_t h   = (int64_t)lt.tm_hour;
        int64_t mi  = (int64_t)lt.tm_min;
        int64_t s   = (int64_t)lt.tm_sec;

        int64_t days = days_since_1900(y, mo, d);
        double total = (double)days * 86400.0
                     + (double)h * 3600.0
                     + (double)mi * 60.0
                     + (double)s;
        return expr_new_real(total);
    }

    if (argc != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* AbsoluteTime[t] for numeric t passes through as an existing spec. */
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || arg->type == EXPR_BIGINT) {
        return expr_copy(arg);
    }

    /* AbsoluteTime[{y, m, d, h, m, s}] -- list with possible elision. */
    if (arg->type == EXPR_FUNCTION
        && arg->data.function.head->type == EXPR_SYMBOL
        && strcmp(arg->data.function.head->data.symbol, "List") == 0) {

        size_t n = arg->data.function.arg_count;
        if (n < 1 || n > 6) return NULL;

        /* Defaults match Mathematica: {y, 1, 1, 0, 0, 0}. */
        double parts[6] = {0.0, 1.0, 1.0, 0.0, 0.0, 0.0};
        for (size_t i = 0; i < n; i++) {
            if (!expr_to_double_strict(arg->data.function.args[i], &parts[i])) {
                return NULL;
            }
        }

        /* Year and month must be integer-valued; day/h/m/s may be fractional. */
        if (parts[0] != floor(parts[0])) return NULL;
        if (parts[1] != floor(parts[1])) return NULL;

        int64_t y  = (int64_t)parts[0];
        int64_t mo = (int64_t)parts[1];

        /* Split day into integer + fractional parts so the JDN formula sees an integer. */
        double d_floor = floor(parts[2]);
        int64_t d_int  = (int64_t)d_floor;
        double d_frac  = parts[2] - d_floor;

        int64_t days = days_since_1900(y, mo, d_int);
        double total = (double)days * 86400.0
                     + d_frac * 86400.0
                     + parts[3] * 3600.0
                     + parts[4] * 60.0
                     + parts[5];

        /* Return an exact integer when all inputs and the total are integral. */
        int all_int = (d_frac == 0.0)
                   && (parts[3] == floor(parts[3]))
                   && (parts[4] == floor(parts[4]))
                   && (parts[5] == floor(parts[5]));
        if (all_int) {
            double rounded = floor(total + 0.5);
            if (rounded == total && rounded >= (double)INT64_MIN && rounded <= (double)INT64_MAX) {
                return expr_new_integer((int64_t)rounded);
            }
        }
        return expr_new_real(total);
    }

    return NULL;
}

void datetime_init(void) {
    symtab_add_builtin("Timing", builtin_timing);
    symtab_add_builtin("RepeatedTiming", builtin_repeated_timing);
    symtab_add_builtin("AbsoluteTime", builtin_absolute_time);
    symtab_get_def("AbsoluteTime")->attributes |= ATTR_PROTECTED;
}