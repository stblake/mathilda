/* clock_gettime(CLOCK_MONOTONIC) is POSIX, not C99; glibc hides it under
 * -std=c99 unless this macro precedes every #include. Matches src/datetime.c. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"

/* Wall-clock seconds from a monotonic source, for verifying Pause actually
 * waited and that SessionTime tracks elapsed time. */
static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Evaluate `input`, assert the result is the symbol Null, and free it. */
static void assert_evals_to_null(const char* input) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_SYMBOL);
    assert(strcmp(e->data.symbol.name, "Null") == 0);
    expr_free(e);
}

void run_test(const char* input) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("Failed to parse: %s\n", input);
        assert(false);
    }
    Expr* res = evaluate(e);
    if (!res) {
        printf("Failed to evaluate: %s\n", input);
        expr_free(e);
        assert(false);
    }
    
    char* s = expr_to_string_fullform(res);
    printf("EVAL: %s -> %s\n", input, s);
    free(s);
    expr_free(e);
    expr_free(res);
}

void test_timing() {
    // Timing should return a list: {time, result}
    // E.g. Timing[Plus[2, 3]] -> {time, 5}
    Expr* p = parse_expression("Timing[2 + 3]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_FUNCTION);
    assert(strcmp(e->data.function.head->data.symbol.name, "List") == 0);
    assert(e->data.function.arg_count == 2);
    assert(e->data.function.args[0]->type == EXPR_REAL); // time
    assert(e->data.function.args[1]->type == EXPR_INTEGER);
    assert(e->data.function.args[1]->data.integer == 5);
    expr_free(e);
}

void test_repeated_timing() {
    Expr* p = parse_expression("RepeatedTiming[2 + 3]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_FUNCTION);
    assert(strcmp(e->data.function.head->data.symbol.name, "List") == 0);
    assert(e->data.function.arg_count == 2);
    assert(e->data.function.args[0]->type == EXPR_REAL); // time
    assert(e->data.function.args[1]->type == EXPR_INTEGER);
    assert(e->data.function.args[1]->data.integer == 5);
    expr_free(e);
}

static int64_t eval_to_int(const char* input) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_INTEGER);
    int64_t v = e->data.integer;
    expr_free(e);
    return v;
}

static double eval_to_real(const char* input) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    expr_free(p);
    double v;
    if (e->type == EXPR_REAL)         v = e->data.real;
    else if (e->type == EXPR_INTEGER) v = (double)e->data.integer;
    else { assert(0); v = 0.0; }
    expr_free(e);
    return v;
}

void test_absolute_time_datelist_int() {
    /* Exact Mathematica reference values. */
    assert(eval_to_int("AbsoluteTime[{2022,1,1,0,0,0}]") == 3849984000LL);
    assert(eval_to_int("AbsoluteTime[{2022,1}]")          == 3849984000LL);
    assert(eval_to_int("AbsoluteTime[{2022}]")            == 3849984000LL);
}

void test_absolute_time_normalization() {
    /* {2022, 2, 31} normalizes to {2022, 3, 3}. */
    int64_t a = eval_to_int("AbsoluteTime[{2022,2,31}]");
    int64_t b = eval_to_int("AbsoluteTime[{2022,3,3}]");
    assert(a == 3855254400LL);
    assert(a == b);
}

void test_absolute_time_fractional_day() {
    /* {2022, 3, 15.5} = {2022, 3, 15} + half a day (43200 s). */
    double v = eval_to_real("AbsoluteTime[{2022,3,15.5}]");
    int64_t base = eval_to_int("AbsoluteTime[{2022,3,15,0,0,0}]");
    assert(v == (double)base + 43200.0);
}

void test_absolute_time_fractional_hour() {
    /* {2022, 3, 15, 12.3} = {2022, 3, 15} + 12.3 hours. */
    double v = eval_to_real("AbsoluteTime[{2022,3,15,12.3}]");
    int64_t base = eval_to_int("AbsoluteTime[{2022,3,15,0,0,0}]");
    double expected = (double)base + 12.3 * 3600.0;
    assert(v == expected);
}

void test_absolute_time_passthrough() {
    /* AbsoluteTime[t] for numeric t returns t unchanged. */
    assert(eval_to_int("AbsoluteTime[3849984000]") == 3849984000LL);
}

void test_absolute_time_now() {
    /* AbsoluteTime[] returns a real and stays in a plausible 21st-century range. */
    Expr* p = parse_expression("AbsoluteTime[]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_REAL);
    /* 2020-01-01 .. 2200-01-01 in AbsoluteTime seconds. */
    assert(e->data.real > 3.78e9);
    assert(e->data.real < 9.46e9);
    expr_free(e);
}

void test_absolute_time_attributes() {
    /* Must be Protected per the spec. */
    Expr* p = parse_expression("MemberQ[Attributes[AbsoluteTime], Protected]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_SYMBOL);
    assert(strcmp(e->data.symbol.name, "True") == 0);
    expr_free(e);
}

/* Evaluate a DateList expression and assert the result is {y,m,d,h,mi,s} with
 * the given integer y..mi and second within tol of `sec`. Structural, so it does
 * not couple to the printer's Real formatting. */
static void assert_datelist(const char* input,
                            int64_t y, int64_t mo, int64_t d,
                            int64_t h, int64_t mi, double sec, double tol) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_FUNCTION);
    assert(strcmp(e->data.function.head->data.symbol.name, "List") == 0);
    assert(e->data.function.arg_count == 6);
    Expr** a = e->data.function.args;
    assert(a[0]->type == EXPR_INTEGER && a[0]->data.integer == y);
    assert(a[1]->type == EXPR_INTEGER && a[1]->data.integer == mo);
    assert(a[2]->type == EXPR_INTEGER && a[2]->data.integer == d);
    assert(a[3]->type == EXPR_INTEGER && a[3]->data.integer == h);
    assert(a[4]->type == EXPR_INTEGER && a[4]->data.integer == mi);
    assert(a[5]->type == EXPR_REAL);
    assert(a[5]->data.real >= sec - tol && a[5]->data.real <= sec + tol);
    expr_free(e);
}

void test_datelist_absolute_time() {
    /* Invert an AbsoluteTime number; round-trips the existing AbsoluteTime anchor
     * (AbsoluteTime[{2022,1,1,0,0,0}] == 3849984000). */
    assert_datelist("DateList[3849984000]", 2022, 1, 1, 0, 0, 0.0, 0.0);
    assert_datelist("DateList[3957775000]", 2025, 6, 1, 13, 56, 40.0, 1e-6);
}

void test_datelist_elision() {
    /* Trailing fields default to {_,1,1,0,0,0}; a Real day still yields 0 h/m/s. */
    assert_datelist("DateList[{2026,9}]",    2026, 9, 1, 0, 0, 0.0, 0.0);
    assert_datelist("DateList[{2026,9,28.}]", 2026, 9, 28, 0, 0, 0.0, 1e-6);
}

void test_datelist_normalization() {
    /* Out-of-range and fractional fields reduce, exactly like AbsoluteTime. */
    assert_datelist("DateList[{2022,3,15.5}]", 2022, 3, 15, 12, 0, 0.0, 1e-6);
    assert_datelist("DateList[{2022,0}]",      2021, 12, 1, 0, 0, 0.0, 0.0);
    assert_datelist("DateList[{2022,1,0}]",    2021, 12, 31, 0, 0, 0.0, 0.0);
}

void test_datelist_fractional_hour() {
    /* {2026,9,28,8.1}: date exact, and the intraday split reconstructs to
     * 8.1 hours (29160 s) within a loose tolerance regardless of the last-ULP
     * residue (Mathematica shows 8h 6m 4.77e-7 s; we land on 8h 6m 0s). */
    Expr* p = parse_expression("DateList[{2026,9,28,8.1}]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_FUNCTION && e->data.function.arg_count == 6);
    Expr** a = e->data.function.args;
    assert(a[0]->data.integer == 2026 && a[1]->data.integer == 9 && a[2]->data.integer == 28);
    double intraday = (double)a[3]->data.integer * 3600.0
                    + (double)a[4]->data.integer * 60.0
                    + a[5]->data.real;
    assert(intraday >= 29160.0 - 1e-3 && intraday <= 29160.0 + 1e-3);
    expr_free(e);
}

void test_datelist_roundtrip() {
    /* DateList[] is the inverse of AbsoluteTime for an exact spec. */
    assert(eval_to_int("AbsoluteTime[DateList[3849984000]]") == 3849984000LL);
}

void test_datelist_strings() {
    assert_datelist("DateList[\"28 Sep, 2026\"]", 2026, 9, 28, 0, 0, 0.0, 0.0);
    assert_datelist("DateList[\"30 Oct 2026\"]",  2026, 10, 30, 0, 0, 0.0, 0.0);
    /* Ambiguous numeric string: US M/D/Y default, YearShort 1 -> 2001. */
    assert_datelist("DateList[\"05/10/1\"]",      2001, 5, 10, 0, 0, 0.0, 0.0);
}

void test_datelist_format() {
    assert_datelist("DateList[{\"09/28/26\",{\"Day\",\"Month\",\"YearShort\"}}]",
                    2028, 4, 9, 0, 0, 0.0, 0.0);
    assert_datelist("DateList[{\"9/28/2026\",{\"Month\",\"Day\",\"Year\"}}]",
                    2026, 9, 28, 0, 0, 0.0, 0.0);
    /* Explicit separators between elements. */
    assert_datelist("DateList[{\"9/28/2026\",{\"Month\",\"/\",\"Day\",\"/\",\"Year\"}}]",
                    2026, 9, 28, 0, 0, 0.0, 0.0);
    /* YearShort 05 -> 2005. */
    assert_datelist("DateList[{\"05/10/1\",{\"YearShort\",\"Day\",\"Month\"}}]",
                    2005, 1, 10, 0, 0, 0.0, 0.0);
}

void test_datelist_current_year_fill() {
    /* A string spec with no year fills the current calendar year. */
    time_t now = time(NULL);
    struct tm* lp = localtime(&now);
    int64_t cur_year = (int64_t)lp->tm_year + 1900;
    assert_datelist("DateList[{\"2/15\",{\"Month\",\"Day\"}}]",
                    cur_year, 2, 15, 0, 0, 0.0, 0.0);
}

void test_datelist_bad_month_unevaluated() {
    /* Non-integer month cannot be a date (DateList::arg): stays unevaluated. */
    Expr* p = parse_expression("DateList[{2022,3.5}]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_FUNCTION);
    assert(strcmp(e->data.function.head->data.symbol.name, "DateList") == 0);
    assert(e->data.function.arg_count == 1);
    expr_free(e);
}

void test_datelist_now_shape() {
    /* DateList[] is a 6-element list: five integers plus a Real second. */
    Expr* p = parse_expression("DateList[]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_FUNCTION);
    assert(strcmp(e->data.function.head->data.symbol.name, "List") == 0);
    assert(e->data.function.arg_count == 6);
    for (int i = 0; i < 5; i++) assert(e->data.function.args[i]->type == EXPR_INTEGER);
    assert(e->data.function.args[5]->type == EXPR_REAL);
    /* Plausible 21st-century year. */
    assert(e->data.function.args[0]->data.integer >= 2020);
    assert(e->data.function.args[0]->data.integer <= 2200);
    expr_free(e);
}

void test_datelist_attributes() {
    Expr* p = parse_expression("MemberQ[Attributes[DateList], Protected]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_SYMBOL);
    assert(strcmp(e->data.symbol.name, "True") == 0);
    expr_free(e);
}

void test_pause_null() {
    /* Pause returns Null; zero and negative durations return immediately. */
    assert_evals_to_null("Pause[0]");
    assert_evals_to_null("Pause[-1]");
    assert_evals_to_null("Pause[0.05]");
}

void test_pause_waits() {
    /* Pause[t] blocks for at least t seconds of wall-clock time. */
    double t0 = now_seconds();
    assert_evals_to_null("Pause[0.05]");
    double elapsed = now_seconds() - t0;
    assert(elapsed >= 0.045);   /* nanosleep guarantees >= t; small slack */
}

void test_pause_symbolic_unevaluated() {
    /* Pause[x] with a non-numeric argument stays unevaluated. */
    Expr* p = parse_expression("Pause[x]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_FUNCTION);
    assert(strcmp(e->data.function.head->data.symbol.name, "Pause") == 0);
    assert(e->data.function.arg_count == 1);
    expr_free(e);
}

void test_pause_attributes() {
    /* Must be Protected per the spec. */
    Expr* p = parse_expression("MemberQ[Attributes[Pause], Protected]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_SYMBOL);
    assert(strcmp(e->data.symbol.name, "True") == 0);
    expr_free(e);
}

void test_pause_timing_cpu_zero() {
    /* Timing (CPU clock) must NOT count Pause: near-zero CPU, result Null. */
    Expr* p = parse_expression("Timing[Pause[0.05]]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_FUNCTION);
    assert(strcmp(e->data.function.head->data.symbol.name, "List") == 0);
    assert(e->data.function.arg_count == 2);
    assert(e->data.function.args[0]->type == EXPR_REAL);
    assert(e->data.function.args[0]->data.real < 0.03);   /* no CPU burned */
    assert(e->data.function.args[1]->type == EXPR_SYMBOL);
    assert(strcmp(e->data.function.args[1]->data.symbol.name, "Null") == 0);
    expr_free(e);
}

void test_pause_absolute_timing_wall() {
    /* AbsoluteTiming (wall clock) MUST count Pause: >= requested, result Null. */
    Expr* p = parse_expression("AbsoluteTiming[Pause[0.05]]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_FUNCTION);
    assert(e->data.function.arg_count == 2);
    assert(e->data.function.args[0]->type == EXPR_REAL);
    assert(e->data.function.args[0]->data.real >= 0.045);
    assert(e->data.function.args[1]->type == EXPR_SYMBOL);
    assert(strcmp(e->data.function.args[1]->data.symbol.name, "Null") == 0);
    expr_free(e);
}

void test_time_unit() {
    double u = eval_to_real("$TimeUnit");
    assert(u > 0.0);
    assert(u <= 0.01);   /* a small fraction of a second */
}

void test_session_time_counts_pause() {
    /* SessionTime is wall-clock: it advances across a Pause. */
    double before = eval_to_real("SessionTime[]");
    assert(before >= 0.0);
    assert_evals_to_null("Pause[0.05]");
    double after = eval_to_real("SessionTime[]");
    assert(after - before >= 0.045);
}

void test_time_used_excludes_pause() {
    /* TimeUsed is CPU time: it barely moves across a Pause. */
    double before = eval_to_real("TimeUsed[]");
    assert(before >= 0.0);
    assert_evals_to_null("Pause[0.05]");
    double after = eval_to_real("TimeUsed[]");
    assert(after - before < 0.03);
}

int main() {
    symtab_init();
    core_init();

    printf("Running datetime tests...\n");
    test_timing();
    test_repeated_timing();
    test_absolute_time_datelist_int();
    test_absolute_time_normalization();
    test_absolute_time_fractional_day();
    test_absolute_time_fractional_hour();
    test_absolute_time_passthrough();
    test_absolute_time_now();
    test_absolute_time_attributes();
    test_datelist_absolute_time();
    test_datelist_elision();
    test_datelist_normalization();
    test_datelist_fractional_hour();
    test_datelist_roundtrip();
    test_datelist_strings();
    test_datelist_format();
    test_datelist_current_year_fill();
    test_datelist_bad_month_unevaluated();
    test_datelist_now_shape();
    test_datelist_attributes();
    test_pause_null();
    test_pause_waits();
    test_pause_symbolic_unevaluated();
    test_pause_attributes();
    test_pause_timing_cpu_zero();
    test_pause_absolute_timing_wall();
    test_time_unit();
    test_session_time_counts_pause();
    test_time_used_excludes_pause();
    printf("All datetime tests passed!\n");
    symtab_clear();
    return 0;
}
