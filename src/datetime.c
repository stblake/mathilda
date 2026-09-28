/* AbsoluteTiming uses clock_gettime(CLOCK_MONOTONIC), which is POSIX, not C99.
 * glibc hides it under -std=c99 while Darwin exposes it implicitly, so the
 * feature-test macro must come before the first #include -- below it the header
 * has already been parsed with the wrong namespace.  Matches src/core.c. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "datetime.h"
#include "common.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "numeric.h"
#include "message.h"
#include "print.h"
#include <time.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include <gmp.h>

/* Emit a Wolfram-style diagnostic ("Head::tag: text") unless Quiet is active.
 * The firing is always noted so Check[] sees it. Mirrors ops_msg in
 * src/assoc_ops.c. */
static void dt_msg(const char* fmt, ...) {
    mth_msg_note_fired();
    if (mth_msg_suppressed()) return;
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/*
 * Elapsed wall-clock seconds from a monotonic source.
 *
 * WHY THIS IS NOT clock().  clock() reports CPU time summed across threads, so
 * any operation that uses nd_parallel_for / nd_parallel_reduce or the platform
 * BLAS reports roughly cores x its true duration -- Timing[Total[bigArray]] on
 * an 8-core host reads ~8x the time the user actually waited.  That makes
 * Timing[] unusable for benchmarking the threaded NDArray paths.
 *
 * CLOCK_MONOTONIC rather than CLOCK_REALTIME so an NTP step or a manual clock
 * change during a long evaluation cannot produce a negative interval.
 */
static double dt_wall_seconds(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0)
        return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
    /* No monotonic clock: fall back to CPU time, which is at least a duration.
     * Over-reports threaded work, exactly as Timing[] does. */
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

/* Wall-clock reading captured at datetime_init(), i.e. kernel start-up. It is
 * the zero point SessionTime[] measures from -- the same monotonic source as
 * AbsoluteTiming, so time spent inside Pause[] is counted here too. */
static double g_session_start = 0.0;

Expr* builtin_absolute_timing(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    double start = dt_wall_seconds();
    Expr* evaluated = evaluate(arg);
    double elapsed = dt_wall_seconds() - start;

    Expr** results = malloc(sizeof(Expr*) * 2);
    if (!results) { expr_free(evaluated); return NULL; }
    results[0] = expr_new_real(elapsed);
    results[1] = evaluated;

    Expr* final_res = expr_new_function(expr_new_symbol(SYM_List), results, 2);
    free(results);
    return final_res;
}

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
    
    Expr* final_res = expr_new_function(expr_new_symbol(SYM_List), results, 2);
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
    
    Expr* final_res = expr_new_function(expr_new_symbol(SYM_List), results, 2);
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

/*
 * Convert a {y, m, d, h, mi, s} parts array to absolute seconds since
 * 1900-01-01 (local, no TZ/DST/leap-second correction). Year and month must be
 * integer-valued -- the lengths of years and months vary, so those two fields
 * cannot be fractional; returns 0 (leaving *total untouched) if either is not.
 * Day / hour / minute / second may be fractional and out of range; they are
 * normalised implicitly by the JDN + seconds arithmetic (e.g. {2022,2,31} =
 * {2022,3,3}). Shared by AbsoluteTime and DateList so the two heads agree
 * bit-for-bit on the epoch and the accumulation order.
 */
static int datelist_parts_to_abstime(const double parts[6], double* total) {
    if (parts[0] != floor(parts[0])) return 0;   /* year  */
    if (parts[1] != floor(parts[1])) return 0;   /* month */

    int64_t y  = (int64_t)parts[0];
    int64_t mo = (int64_t)parts[1];

    /* Split the day into integer + fractional parts so the JDN formula sees an
     * integer; the fractional day is carried as seconds. */
    double d_floor = floor(parts[2]);
    int64_t d_int  = (int64_t)d_floor;
    double d_frac  = parts[2] - d_floor;

    int64_t days = days_since_1900(y, mo, d_int);
    *total = (double)days * 86400.0
           + d_frac * 86400.0
           + parts[3] * 3600.0
           + parts[4] * 60.0
           + parts[5];
    return 1;
}

/*
 * Inverse of datelist_parts_to_abstime: split absolute seconds since 1900-01-01
 * into a normalised broken-down local date/time. The date part inverts
 * days_since_1900 with the standard inverse Fliegel & Van Flandern formula; the
 * intraday part is split from the same Real that AbsoluteTime would return, so
 * DateList[AbsoluteTime[spec]] and DateList[spec] agree, and any fractional
 * round-trip carries whatever floating residue that single Real held. Year..
 * minute come back as exact integers; seconds as a (possibly fractional) double.
 */
static void gregorian_from_abstime(double total,
                                   int64_t* y, int64_t* mo, int64_t* d,
                                   int64_t* h, int64_t* mi, double* s) {
    /* Whole days since epoch (floor, so a negative time lands on the earlier
     * day) and the in-day remainder in [0, 86400). The guard also absorbs a
     * division that rounded floor(total/86400) one day too high. */
    double day_f = floor(total / 86400.0);
    double rem   = total - day_f * 86400.0;
    if (rem < 0.0)       { rem += 86400.0; day_f -= 1.0; }
    if (rem >= 86400.0)  { rem -= 86400.0; day_f += 1.0; }

    int64_t jdn = (int64_t)day_f + 2415021;   /* days_since_1900 + JDN(1900-01-01) */

    /* Inverse Fliegel & Van Flandern (pure integer arithmetic). */
    int64_t a  = jdn + 32044;
    int64_t b  = (4 * a + 3) / 146097;
    int64_t c  = a - (146097 * b) / 4;
    int64_t dd = (4 * c + 3) / 1461;
    int64_t e  = c - (1461 * dd) / 4;
    int64_t m  = (5 * e + 2) / 153;

    *d  = e - (153 * m + 2) / 5 + 1;
    *mo = m + 3 - 12 * (m / 10);
    *y  = 100 * b + dd - 4800 + (m / 10);

    /* Intraday split, same accumulation order as the forward direction. */
    double hh = floor(rem / 3600.0);
    rem -= hh * 3600.0;
    double mm = floor(rem / 60.0);
    rem -= mm * 60.0;
    *h  = (int64_t)hh;
    *mi = (int64_t)mm;
    *s  = rem;
}

/* Build the {y, m, d, h, mi, s} result list: five integers plus a Real second.
 * Adopts nothing from the caller; returns a freshly-owned List. */
static Expr* datelist_make_result(int64_t y, int64_t mo, int64_t d,
                                  int64_t h, int64_t mi, double s) {
    Expr** parts = malloc(sizeof(Expr*) * 6);
    parts[0] = expr_new_integer(y);
    parts[1] = expr_new_integer(mo);
    parts[2] = expr_new_integer(d);
    parts[3] = expr_new_integer(h);
    parts[4] = expr_new_integer(mi);
    parts[5] = expr_new_real(s);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), parts, 6);
    free(parts);
    return out;
}

/*
 * Read a {y, m, d, h, mi, s} list (as built by datelist_make_result, or any
 * DateList result) back into absolute seconds since 1900-01-01 -- the exact
 * inverse of datelist_make_result for normalised fields. Returns 0 (leaving
 * *total untouched) unless the argument is a List of six numeric elements. Lets
 * UnixTime reuse the DateList string/format parsers verbatim without re-deriving
 * the calendar arithmetic.
 */
static int datelist_result_to_abstime(const Expr* list, double* total) {
    if (!head_is(list, SYM_List) || list->data.function.arg_count != 6) return 0;
    double parts[6];
    for (int i = 0; i < 6; i++) {
        if (!expr_to_double_strict(list->data.function.args[i], &parts[i])) return 0;
    }
    return datelist_parts_to_abstime(parts, total);
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
    if (head_is(arg, SYM_List)) {

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
        double total;
        if (!datelist_parts_to_abstime(parts, &total)) return NULL;

        /* Return an exact integer when all inputs and the total are integral. */
        int all_int = (parts[2] == floor(parts[2]))
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

/* ======================================================================
 * DateList -- broken-down date/time {y, m, d, h, mi, s}.
 *
 * Every form (current time, an absolute-time number, a {y,m,...} spec with
 * elision, a date string, or a {"string", {elements}} format spec) is reduced
 * to a raw parts[6] array and then run through the same
 * datelist_parts_to_abstime -> gregorian_from_abstime pipeline, so out-of-range
 * fields normalise identically to AbsoluteTime and to Mathematica. Uses local
 * time with no time-zone / DST / leap-second correction, exactly like
 * AbsoluteTime.
 * ====================================================================== */

/* Format element codes for the {"string", {e1, e2, ...}} form. */
enum {
    EL_NONE = 0, EL_YEAR, EL_YEARSHORT, EL_QUARTER, EL_MONTH, EL_MONTHNAME,
    EL_DAY, EL_DAYNAME, EL_HOUR, EL_HOUR12, EL_AMPM, EL_MINUTE, EL_SECOND,
    EL_MILLISECOND
};

static int element_code(const char* s) {
    if (!strcmp(s, "Year"))        return EL_YEAR;
    if (!strcmp(s, "YearShort"))   return EL_YEARSHORT;
    if (!strcmp(s, "Quarter"))     return EL_QUARTER;
    if (!strcmp(s, "Month"))       return EL_MONTH;
    if (!strcmp(s, "MonthName"))   return EL_MONTHNAME;
    if (!strcmp(s, "Day"))         return EL_DAY;
    if (!strcmp(s, "DayName"))     return EL_DAYNAME;
    if (!strcmp(s, "Hour"))        return EL_HOUR;
    if (!strcmp(s, "Hour12"))      return EL_HOUR12;
    if (!strcmp(s, "AMPM"))        return EL_AMPM;
    if (!strcmp(s, "Minute"))      return EL_MINUTE;
    if (!strcmp(s, "Second"))      return EL_SECOND;
    if (!strcmp(s, "Millisecond")) return EL_MILLISECOND;
    return EL_NONE;
}

/* Month 1..12 for a full or abbreviated (>=3 letter) English month name,
 * case-insensitive prefix match on the canonical name; 0 if not a month. */
static int month_from_name(const char* tok) {
    static const char* names[12] = {
        "january", "february", "march", "april", "may", "june",
        "july", "august", "september", "october", "november", "december"
    };
    size_t len = strlen(tok);
    if (len < 3) return 0;
    for (int m = 0; m < 12; m++) {
        size_t nlen = strlen(names[m]);
        if (len > nlen) continue;
        int ok = 1;
        for (size_t i = 0; i < len; i++) {
            if (tolower((unsigned char)tok[i]) != names[m][i]) { ok = 0; break; }
        }
        if (ok) return m + 1;
    }
    return 0;
}

/* True for a full or abbreviated (>=3 letter) English weekday name. */
static int is_day_name(const char* tok) {
    static const char* names[7] = {
        "monday", "tuesday", "wednesday", "thursday", "friday", "saturday", "sunday"
    };
    size_t len = strlen(tok);
    if (len < 3) return 0;
    for (int i = 0; i < 7; i++) {
        size_t nlen = strlen(names[i]);
        if (len > nlen) continue;
        int ok = 1;
        for (size_t j = 0; j < len; j++) {
            if (tolower((unsigned char)tok[j]) != names[i][j]) { ok = 0; break; }
        }
        if (ok) return 1;
    }
    return 0;
}

/* 1 for "pm", 0 for "am" (case-insensitive), -1 otherwise. */
static int ampm_of(const char* t) {
    if ((t[0] == 'a' || t[0] == 'A') && (t[1] == 'm' || t[1] == 'M') && t[2] == '\0') return 0;
    if ((t[0] == 'p' || t[0] == 'P') && (t[1] == 'm' || t[1] == 'M') && t[2] == '\0') return 1;
    return -1;
}

static int tok_is_numeric(const char* t) {
    if (!*t) return 0;
    for (const char* p = t; *p; p++) if (!isdigit((unsigned char)*p)) return 0;
    return 1;
}

/* Two- or one-digit years map to 2000+yy ("26" -> 2026, "1" -> 2001), matching
 * Mathematica's default; longer strings are taken verbatim. */
static double norm_year(const char* t) {
    long v = atol(t);
    if (strlen(t) <= 2) return (double)(2000 + v);
    return (double)v;
}

#define DT_MAX_TOKENS 16
#define DT_TOKEN_LEN  64

/* Split s into maximal alphanumeric runs (any non-alphanumeric run is a
 * separator). Fills tok[][], returns the count (capped at DT_MAX_TOKENS). */
static int dt_tokenize(const char* s, char tok[DT_MAX_TOKENS][DT_TOKEN_LEN]) {
    int count = 0;
    const char* p = s;
    while (*p && count < DT_MAX_TOKENS) {
        while (*p && !isalnum((unsigned char)*p)) p++;
        if (!*p) break;
        int len = 0;
        while (*p && isalnum((unsigned char)*p)) {
            if (len < DT_TOKEN_LEN - 1) tok[count][len++] = *p;
            p++;
        }
        tok[count][len] = '\0';
        count++;
    }
    return count;
}

/* Apply one token to the parts array for a given format element. Returns 0 only
 * on a genuine parse error (e.g. an unrecognised month name). */
static int apply_element(int code, const char* t, double parts[6],
                         int* pm, int* have_hour12) {
    switch (code) {
        case EL_YEAR:        parts[0] = (double)atol(t); break;
        case EL_YEARSHORT:   parts[0] = (double)(2000 + atol(t)); break;
        case EL_QUARTER: {
            long q = atol(t);
            if (q < 1) q = 1;
            parts[1] = (double)((q - 1) * 3 + 1);
            break;
        }
        case EL_MONTH:       parts[1] = (double)atol(t); break;
        case EL_MONTHNAME: {
            int m = month_from_name(t);
            if (!m) return 0;
            parts[1] = (double)m;
            break;
        }
        case EL_DAY:         parts[2] = strtod(t, NULL); break;
        case EL_DAYNAME:     break;  /* recognised, contributes nothing */
        case EL_HOUR:        parts[3] = strtod(t, NULL); break;
        case EL_HOUR12:      parts[3] = strtod(t, NULL); *have_hour12 = 1; break;
        case EL_AMPM:        *pm = (tolower((unsigned char)t[0]) == 'p') ? 1 : 0; break;
        case EL_MINUTE:      parts[4] = strtod(t, NULL); break;
        case EL_SECOND:      parts[5] = strtod(t, NULL); break;
        case EL_MILLISECOND: parts[5] += strtod(t, NULL) / 1000.0; break;
        default: return 0;
    }
    return 1;
}

/* Current local calendar year, for filling year-less string specs. */
static int64_t current_year(void) {
    time_t now = time(NULL);
    if (now == (time_t)-1) return 1900;
    struct tm* lp = localtime(&now);
    if (!lp) return 1900;
    return (int64_t)lp->tm_year + 1900;
}

/* DateList[] -- current local date and time, with sub-second precision from the
 * realtime clock (localtime only resolves to the second). */
static Expr* datelist_now(void) {
    time_t secs;
    long   nsec = 0;
    struct timespec ts;
    if (clock_gettime(CLOCK_REALTIME, &ts) == 0) {
        secs = ts.tv_sec;
        nsec = ts.tv_nsec;
    } else {
        secs = time(NULL);
        if (secs == (time_t)-1) return NULL;
    }
    struct tm lt;
    struct tm* lp = localtime(&secs);
    if (!lp) return NULL;
    lt = *lp;
    return datelist_make_result(
        (int64_t)lt.tm_year + 1900,
        (int64_t)lt.tm_mon + 1,
        (int64_t)lt.tm_mday,
        (int64_t)lt.tm_hour,
        (int64_t)lt.tm_min,
        (double)lt.tm_sec + (double)nsec * 1e-9);
}

/* Finish a string/format spec: apply 12-hour clock, normalise, build result. */
static Expr* datelist_finish(double parts[6], int pm, int have_hour12) {
    if (have_hour12 && pm >= 0) {
        if (pm == 1 && parts[3] < 12.0) parts[3] += 12.0;
        if (pm == 0 && parts[3] == 12.0) parts[3] = 0.0;
    }
    double total;
    if (!datelist_parts_to_abstime(parts, &total)) return NULL;
    int64_t y, mo, d, h, mi; double s;
    gregorian_from_abstime(total, &y, &mo, &d, &h, &mi, &s);
    return datelist_make_result(y, mo, d, h, mi, s);
}

/* DateList["string"] -- best-effort parse of a free-form date string. */
static Expr* datelist_from_string(const char* s) {
    char tok[DT_MAX_TOKENS][DT_TOKEN_LEN];
    int nt = dt_tokenize(s, tok);
    if (nt == 0) return NULL;

    double parts[6] = { (double)current_year(), 1.0, 1.0, 0.0, 0.0, 0.0 };
    int pm = -1;
    int month = 0;
    int num_idx[DT_MAX_TOKENS], nn = 0;

    /* Classify word tokens; collect numeric-token indices in order. */
    for (int i = 0; i < nt; i++) {
        if (tok_is_numeric(tok[i])) { num_idx[nn++] = i; continue; }
        int m = month_from_name(tok[i]);
        if (m) { if (month) return NULL; month = m; continue; }
        if (is_day_name(tok[i])) continue;
        int ap = ampm_of(tok[i]);
        if (ap >= 0) { pm = ap; continue; }
        return NULL;   /* an unrecognised word -> not a date we can read */
    }

    int ambiguous = 0;
    if (month) {
        /* Month name present: numerics are day, 4-digit year, then h, m, s. */
        parts[1] = (double)month;
        int day_set = 0, year_set = 0, time_i = 3;
        for (int k = 0; k < nn; k++) {
            const char* t = tok[num_idx[k]];
            long v = atol(t);
            if (!year_set && strlen(t) == 4) { parts[0] = (double)v; year_set = 1; }
            else if (!day_set && v >= 1 && v <= 31 && strlen(t) <= 2) { parts[2] = (double)v; day_set = 1; }
            else if (time_i <= 5) { parts[time_i++] = (double)v; }
        }
        if (!day_set) return NULL;
    } else {
        /* All numeric. A leading 4-digit token reads as ISO Y/M/D; otherwise the
         * order is genuinely ambiguous -- take the US default M/D/Y and warn. */
        if (nn == 0) return NULL;
        int base = 3;
        if (strlen(tok[num_idx[0]]) == 4) {
            parts[0] = (double)atol(tok[num_idx[0]]);
            if (nn > 1) parts[1] = (double)atol(tok[num_idx[1]]);
            if (nn > 2) parts[2] = (double)atol(tok[num_idx[2]]);
        } else {
            ambiguous = 1;
            parts[1] = (double)atol(tok[num_idx[0]]);
            if (nn > 1) parts[2] = (double)atol(tok[num_idx[1]]);
            if (nn > 2) parts[0] = norm_year(tok[num_idx[2]]);
        }
        for (int k = 3; k < nn && (base + k - 3) <= 5; k++) {
            parts[base + k - 3] = (double)atol(tok[num_idx[k]]);
        }
    }

    if (ambiguous) {
        dt_msg("DateList::ambig: Warning: the interpretation of the string %s "
               "as a date is ambiguous.", s);
    }
    return datelist_finish(parts, pm, 0);
}

/* DateList[{"string", {e1, e2, ...}}] -- format-element-driven parse. Separator
 * strings in the element list are ignored: the string is tokenised on any
 * non-alphanumeric run and tokens are assigned to the element entries in order. */
static Expr* datelist_from_format(Expr* arg) {
    size_t n = arg->data.function.arg_count;
    const char* s = arg->data.function.args[0]->data.string;

    if (n < 2 || !head_is(arg->data.function.args[1], SYM_List)) {
        return datelist_from_string(s);   /* {"string"} alone */
    }
    Expr* fmt = arg->data.function.args[1];
    size_t nf = fmt->data.function.arg_count;

    char tok[DT_MAX_TOKENS][DT_TOKEN_LEN];
    int nt = dt_tokenize(s, tok);

    double parts[6] = { (double)current_year(), 1.0, 1.0, 0.0, 0.0, 0.0 };
    int pm = -1, have_hour12 = 0, ti = 0;

    for (size_t fi = 0; fi < nf; fi++) {
        Expr* fe = fmt->data.function.args[fi];
        if (fe->type != EXPR_STRING) continue;
        int code = element_code(fe->data.string);
        if (code == EL_NONE) continue;      /* explicit separator: skip */
        if (ti >= nt) continue;             /* ran out of tokens: keep default */
        if (!apply_element(code, tok[ti++], parts, &pm, &have_hour12)) return NULL;
    }
    return datelist_finish(parts, pm, have_hour12);
}

Expr* builtin_date_list(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 0) return datelist_now();
    if (argc != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* DateList[time] -- an AbsoluteTime specification (seconds since 1900). */
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || arg->type == EXPR_BIGINT) {
        double total;
        expr_to_double_strict(arg, &total);
        int64_t y, mo, d, h, mi; double s;
        gregorian_from_abstime(total, &y, &mo, &d, &h, &mi, &s);
        return datelist_make_result(y, mo, d, h, mi, s);
    }

    /* DateList["string"]. */
    if (arg->type == EXPR_STRING) {
        return datelist_from_string(arg->data.string);
    }

    if (head_is(arg, SYM_List)) {
        size_t n = arg->data.function.arg_count;

        /* {"string", ...} -- a DateString / format specification. */
        if (n >= 1 && arg->data.function.args[0]->type == EXPR_STRING) {
            return datelist_from_format(arg);
        }

        /* {y, m, d, h, mi, s} numeric spec, elided from the right. */
        if (n < 1 || n > 6) return NULL;
        double parts[6] = {0.0, 1.0, 1.0, 0.0, 0.0, 0.0};
        for (size_t i = 0; i < n; i++) {
            if (!expr_to_double_strict(arg->data.function.args[i], &parts[i])) {
                return NULL;   /* non-numeric element: leave unevaluated */
            }
        }
        double total;
        if (!datelist_parts_to_abstime(parts, &total)) {
            /* Year or month not integer-valued: the lengths of years and months
             * vary, so those fields cannot be fractional. */
            char* str = expr_to_string(arg);
            dt_msg("DateList::arg: Argument %s cannot be interpreted as a date "
                   "or time input.", str ? str : "?");
            free(str);
            return NULL;
        }
        int64_t y, mo, d, h, mi; double s;
        gregorian_from_abstime(total, &y, &mo, &d, &h, &mi, &s);
        return datelist_make_result(y, mo, d, h, mi, s);
    }

    return NULL;
}

/*
 * Seconds from 1900-01-01 to 1970-01-01: days_since_1900(1970, 1, 1) * 86400 =
 * 25567 * 86400. The fixed offset between an AbsoluteTime (seconds since 1900)
 * and a Unix time (seconds since 1970). No timezone / DST / leap-second
 * correction is applied, matching the rest of this module.
 */
#define DT_UNIX_EPOCH_OFFSET 2208988800.0

/*
 * Turn absolute seconds-since-1900 into a Unix time: subtract the epoch offset
 * and round to the nearest whole second. UnixTime always yields an Integer -- the
 * nearest whole second, as Mathematica does -- so a fractional-second spec rounds
 * rather than returning a Real (this is a deliberate difference from
 * AbsoluteTime). A Real is returned only in the corner case where the value does
 * not fit an int64.
 */
static Expr* unixtime_from_abstime(double total_1900) {
    double unix_secs = total_1900 - DT_UNIX_EPOCH_OFFSET;
    double rounded = floor(unix_secs + 0.5);
    if (rounded >= (double)INT64_MIN && rounded <= (double)INT64_MAX) {
        return expr_new_integer((int64_t)rounded);
    }
    return expr_new_real(unix_secs);
}

/*
 * UnixTime[] / UnixTime[date] -- seconds since 1970-01-01 00:00:00 GMT.
 *
 *   UnixTime[]                     current time, the true POSIX epoch second.
 *   UnixTime[t]                    t taken as an AbsoluteTime (seconds since
 *                                  1900), shifted to the Unix epoch.
 *   UnixTime[{y, m, d, h, mi, s}]  a DateList spec (elided from the right,
 *                                  fields normalised), shifted to the epoch.
 *   UnixTime["string"]             a DateString spec.
 *   UnixTime[{"string", {e, ...}}] a date string with explicit format elements.
 *
 * The date-list and string forms share the DateList backend, so UnixTime and
 * DateList agree on how any spec is interpreted; UnixTime just reports the shifted
 * epoch second instead of the broken-down list.
 */
Expr* builtin_unix_time(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    /* UnixTime[] -- the true POSIX epoch second (GMT), not the local clock. */
    if (argc == 0) {
        time_t now = time(NULL);
        if (now == (time_t)-1) return NULL;
        return expr_new_integer((int64_t)now);
    }

    if (argc != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* UnixTime[t] -- a number is an AbsoluteTime spec (seconds since 1900). */
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || arg->type == EXPR_BIGINT) {
        double total;
        expr_to_double_strict(arg, &total);
        return unixtime_from_abstime(total);
    }

    /* UnixTime["string"] -- reuse the DateList free-form parser. */
    if (arg->type == EXPR_STRING) {
        Expr* dl = datelist_from_string(arg->data.string);
        if (!dl) return NULL;
        double total;
        int ok = datelist_result_to_abstime(dl, &total);
        expr_free(dl);
        return ok ? unixtime_from_abstime(total) : NULL;
    }

    if (head_is(arg, SYM_List)) {
        size_t n = arg->data.function.arg_count;

        /* {"string", ...} -- reuse the DateList format-element parser. */
        if (n >= 1 && arg->data.function.args[0]->type == EXPR_STRING) {
            Expr* dl = datelist_from_format(arg);
            if (!dl) return NULL;
            double total;
            int ok = datelist_result_to_abstime(dl, &total);
            expr_free(dl);
            return ok ? unixtime_from_abstime(total) : NULL;
        }

        /* {y, m, d, h, mi, s} numeric spec, elided from the right. */
        if (n < 1 || n > 6) return NULL;
        double parts[6] = {0.0, 1.0, 1.0, 0.0, 0.0, 0.0};
        for (size_t i = 0; i < n; i++) {
            if (!expr_to_double_strict(arg->data.function.args[i], &parts[i])) {
                return NULL;   /* non-numeric element: leave unevaluated */
            }
        }
        double total;
        if (!datelist_parts_to_abstime(parts, &total)) {
            /* Year or month not integer-valued: the lengths of years and months
             * vary, so those fields cannot be fractional. */
            char* str = expr_to_string(arg);
            dt_msg("UnixTime::arg: Argument %s cannot be interpreted as a date "
                   "or time input.", str ? str : "?");
            free(str);
            return NULL;
        }
        return unixtime_from_abstime(total);
    }

    return NULL;
}

/*
 * Coerce a Pause argument to a machine double. Accepts integers, reals,
 * bignums, rationals, MPFR reals, and any NumericQ symbolic form (Pi, Sqrt[2],
 * ...) via numericalize -- so Pause[1/4] and Pause[Pi] behave like Mathematica.
 * Returns 0 for a genuinely non-numeric argument, leaving Pause[x] unevaluated.
 * Mirrors clip_to_double_value in core.c. */
static int pause_seconds(const Expr* e, double* out) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer;  return 1; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;             return 1; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return 1; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return 1; }
#endif
    int64_t n, d;
    if (is_rational(e, &n, &d) && d != 0) { *out = (double)n / (double)d; return 1; }

    /* Symbolic numeric forms (Pi, E, Sqrt[2], ...) -- ask numericalize. */
    Expr* approx = numericalize(e, numeric_machine_spec());
    if (!approx) return 0;
    int ok = 0;
    if (approx->type == EXPR_INTEGER)      { *out = (double)approx->data.integer; ok = 1; }
    else if (approx->type == EXPR_REAL)    { *out = approx->data.real; ok = isfinite(*out); }
    else if (approx->type == EXPR_BIGINT)  { *out = mpz_get_d(approx->data.bigint); ok = 1; }
#ifdef USE_MPFR
    else if (approx->type == EXPR_MPFR)    { *out = mpfr_get_d(approx->data.mpfr, MPFR_RNDN); ok = isfinite(*out); }
#endif
    else if (is_rational(approx, &n, &d) && d != 0) { *out = (double)n / (double)d; ok = 1; }
    expr_free(approx);
    return ok;
}

/*
 * Pause[n] blocks for at least n seconds of wall-clock time, then returns Null.
 *
 * It sleeps with nanosleep, which consumes no CPU, so the elapsed time is
 * counted by the wall clocks (AbsoluteTiming, SessionTime) but is invisible to
 * the CPU clocks (Timing, TimeUsed) -- no special-casing is needed; the two
 * clock sources do it for free.  nanosleep can return early on a signal, so we
 * resume with the reported remainder until the full duration has elapsed,
 * honouring the "at least n seconds" contract.
 */
Expr* builtin_pause(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL;
    }

    double seconds;
    if (!pause_seconds(res->data.function.args[0], &seconds)) {
        /* Symbolic / non-numeric argument: leave Pause[x] unevaluated. */
        return NULL;
    }

    if (seconds > 0.0 && isfinite(seconds)) {
        double whole = floor(seconds);
        long   nsec  = (long)((seconds - whole) * 1e9);
        if (nsec < 0)         nsec = 0;
        if (nsec > 999999999L) nsec = 999999999L;

        struct timespec req = { (time_t)whole, nsec };
        struct timespec rem;
        while (nanosleep(&req, &rem) != 0 && errno == EINTR) {
            req = rem;
        }
    }

    /* seconds <= 0 (or non-finite) means no wait, matching Pause[0]. */
    return expr_new_symbol(SYM_Null);
}

/* SessionTime[] -- wall-clock seconds since kernel start-up (see g_session_start). */
Expr* builtin_session_time(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 0) {
        return NULL;
    }
    return expr_new_real(dt_wall_seconds() - g_session_start);
}

/* TimeUsed[] -- total CPU seconds used so far this session, via clock(). Does
 * not advance during Pause[], exactly as Timing does not. */
Expr* builtin_time_used(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 0) {
        return NULL;
    }
    return expr_new_real((double)clock() / (double)CLOCKS_PER_SEC);
}

void datetime_init(void) {
    g_session_start = dt_wall_seconds();

    symtab_add_builtin("Timing", builtin_timing);
    symtab_add_builtin("AbsoluteTiming", builtin_absolute_timing);
    symtab_add_builtin("RepeatedTiming", builtin_repeated_timing);
    symtab_add_builtin("AbsoluteTime", builtin_absolute_time);
    symtab_add_builtin("DateList", builtin_date_list);
    symtab_add_builtin("UnixTime", builtin_unix_time);
    symtab_add_builtin("Pause", builtin_pause);
    symtab_add_builtin("SessionTime", builtin_session_time);
    symtab_add_builtin("TimeUsed", builtin_time_used);

    symtab_get_def("AbsoluteTime")->attributes |= ATTR_PROTECTED;
    symtab_get_def("DateList")->attributes      |= ATTR_PROTECTED;
    symtab_get_def("UnixTime")->attributes      |= ATTR_PROTECTED;
    symtab_get_def("Pause")->attributes         |= ATTR_PROTECTED;
    symtab_get_def("SessionTime")->attributes   |= ATTR_PROTECTED;
    symtab_get_def("TimeUsed")->attributes      |= ATTR_PROTECTED;
}