/* Builtin-audit regressions: time and date (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* DateList#1 (high): "DateList returns negative day numbers for years before
 * about -4713".  The Fliegel & Van Flandern formulas use floor division. */
static void test_datelist_far_past(void) {
    audit_chk("DateList#1", "DateList[{-5000, 3, 15}]", "{-5000, 3, 15, 0, 0, 0.0}");
    audit_chk("DateList#1", "DateList[{-100000, 1, 1}]", "{-100000, 1, 1, 0, 0, 0.0}");
    audit_chk("DateList#1", "DateList[{-4713, 11, 24}]", "{-4713, 11, 24, 0, 0, 0.0}");
    audit_chk("DateList#1", "DateList[{-4801, 2, 29}]", "{-4801, 3, 1, 0, 0, 0.0}");
    audit_chk("DateList#1", "DateList[{-4800, 2, 29}]", "{-4800, 2, 29, 0, 0, 0.0}");
    /* every day of a far-past 400-year cycle round-trips */
    audit_true("DateList#1",
               "Module[{t0 = AbsoluteTime[{-6000, 1, 1}]}, "
               "And @@ Table[AbsoluteTime[DateList[t0 + 86400 k]] == t0 + 86400 k, {k, 0, 146097, 997}]]");
    audit_true("DateList#1", "First[DateList[-10^12]] < -29000 && 1 <= DateList[-10^12][[3]] <= 31");
    /* ordinary dates unchanged */
    audit_chk("DateList#1", "DateList[{2024, 2, 29}]", "{2024, 2, 29, 0, 0, 0.0}");
    audit_chk("DateList#1", "DateList[0]", "{1900, 1, 1, 0, 0, 0.0}");
}

/* UnixTime#1 (medium): "UnixTime overflows int64 cast at 2^63 (undefined
 * behaviour, off-by-one here)". */
static void test_unixtime_int64_edge(void) {
    audit_true("UnixTime#1", "UnixTime[2.^63 + 2208988800.] === 2^63");
    audit_true("UnixTime#1", "UnixTime[2.^64 + 2208988800.] === 2^64");
    audit_true("UnixTime#1", "UnixTime[-2.^63 + 2208988800.] === -2^63");
    audit_true("UnixTime#1", "UnixTime[2.^62 + 2208988800.] === 2^62");
    audit_true("UnixTime#1", "UnixTime[{1970, 1, 2}] === 86400");
}

/* UnixTime#4 (low): an exact AbsoluteTime is shifted exactly. */
static void test_s14_unixtime_exact(void) {
    audit_true("UnixTime#4", "UnixTime[10^30] === 10^30 - 2208988800");
    audit_true("UnixTime#4", "UnixTime[2208988800] === 0");
    audit_true("UnixTime#4", "UnixTime[0] === -2208988800");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Timing#1", "Quiet[Check[Timing[1, 2], fired, Timing::argx]]", "fired");
    audit_chk("Timing#1", "Quiet[Timing[1, 2]]", "Timing[1, 2]");
    audit_chk("AbsoluteTiming#1", "Quiet[Check[AbsoluteTiming[1, 2], fired, AbsoluteTiming::argx]]", "fired");
    audit_chk("AbsoluteTiming#1", "Quiet[AbsoluteTiming[1, 2]]", "AbsoluteTiming[1, 2]");
    audit_chk("Pause#1", "Quiet[Check[Pause[1, 2], fired, Pause::argx]]", "fired");
    audit_chk("Pause#1", "Quiet[Pause[1, 2]]", "Pause[1, 2]");
    audit_chk("RepeatedTiming#2", "Quiet[Check[RepeatedTiming[1, 2, 3], fired, RepeatedTiming::argt]]", "fired");
    audit_chk("RepeatedTiming#2", "Quiet[RepeatedTiming[1, 2, 3]]", "RepeatedTiming[1, 2, 3]");
    audit_chk("SessionTime#2", "Quiet[Check[SessionTime[1], fired, SessionTime::argr]]", "fired");
    audit_chk("SessionTime#2", "Quiet[SessionTime[1]]", "SessionTime[1]");
    audit_chk("TimeConstrained#3", "Quiet[Check[TimeConstrained[1, 2, 3, 4], fired, TimeConstrained::argt]]", "fired");
    audit_chk("TimeConstrained#3", "Quiet[TimeConstrained[1, 2, 3, 4]]", "TimeConstrained[1, 2, 3, 4]");
    audit_chk("TimeUsed#2", "Quiet[Check[TimeUsed[1], fired, TimeUsed::argr]]", "fired");
    audit_chk("TimeUsed#2", "Quiet[TimeUsed[1]]", "TimeUsed[1]");
}

static void test_wave4_dates(void) {
    audit_chk("AbsoluteTime#3", "AbsoluteTime[\"2022-03-03\"]", "3855254400");
    audit_chk("AbsoluteTime#3", "AbsoluteTime[{\"2022-03-03\", {\"Year\",\"-\",\"Month\",\"-\",\"Day\"}}]", "3855254400");
    audit_chk("DateList#3", "DateList[\"Sep 28 2026 3 PM\"]", "{2026, 9, 28, 15, 0, 0.0}");
    audit_chk("DateList#4", "Take[DateList[\"3:15 PM\"], {4, 6}]", "{15, 15, 0.0}");
    audit_chk("DateList#4", "Take[DateList[\"3:15 PM\"], 3] == Take[DateList[], 3]", "True");
    audit_chk("DateList#5", "DateList[\"2026-09-28T12:30:00\"]", "{2026, 9, 28, 12, 30, 0.0}");
    audit_chk("DateList#6", "DateList[\"31/12/2022\"]", "{2022, 12, 31, 0, 0, 0.0}");
    audit_chk("DateList#7", "DateList[1/2]", "{1900, 1, 1, 0, 0, 0.5}");
    audit_chk("DateList#7", "DateList[{2022, 1, 1, 0, 0, 1/2}]", "{2022, 1, 1, 0, 0, 0.5}");
    audit_chk("UnixTime#3", "UnixTime[2208988800 + 3/2]", "2");
    audit_chk("UnixTime#3", "UnixTime[{2022, 1, 1, 0, 0, 1/2}]", "1640995200");
}

/* Wave 5 G7 (arith-misc): low-severity findings. Message counts are read
 * off $MessageList inside a Block (each diagnostic once per evaluation). */
static void test_wave5_g7_time_and_date(void) {
    audit_chk("AbsoluteTime#4", "AbsoluteTime[{2022, 1, 1, 0, 0, 10^20}]", "100000000003849984000");
    audit_true("DateList#8", "Abs[DateList[{2022, 5, 6, 7, 8, 9.123456789}][[6]] - 9.123456789] < 10^-12");
    audit_chk("DateList#8", "Most[DateList[{2022, 5, 6, 7, 8, 9.123456789}]]", "{2022, 5, 6, 7, 8}");
    audit_chk("RepeatedTiming#1", "RepeatedTiming[1 + 1, 1/10][[2]]", "2");
    audit_chk("RepeatedTiming#1", "RepeatedTiming[1 + 1, Pi/100][[2]]", "2");
    audit_chk("TimeConstrained#2", "Quiet[Check[TimeConstrained[1 + 1, x], fired, TimeConstrained::timc]]", "fired");
    audit_chk("TimeConstrained#2", "Quiet[Check[TimeConstrained[1 + 1, -Infinity], fired, TimeConstrained::timc]]", "fired");
    audit_chk("TimeConstrained#2", "Quiet[TimeConstrained[1 + 1, 1 + I]]", "TimeConstrained[1 + 1, 1 + I]");
}

/* Wave 6 G1: a BigInt year converts exactly, both ways. */
static void test_wave6_g1_time_and_date(void) {
    audit_chk("AbsoluteTime#2", "AbsoluteTime[{10^30}]", "31556951999999999999999999940041769600");
    audit_chk("AbsoluteTime#2", "AbsoluteTime[{2022, 1, 1, 0, 0, 10^20}]", "100000000003849984000");
    audit_chk("DateList#2", "DateList[{10^20}]", "{100000000000000000000, 1, 1, 0, 0, 0.0}");
    audit_chk("DateList#2", "DateList[AbsoluteTime[{10^30, 7, 15, 3, 4, 5}]]",
              "{1000000000000000000000000000000, 7, 15, 3, 4, 5.0}");
    audit_chk("DateList#2", "With[{d = DateList[10^30]}, {1 <= d[[2]] <= 12, 1 <= d[[3]] <= 31, AbsoluteTime[d] == 10^30}]",
              "{True, True, True}");
    audit_chk("DateList#2", "DateList[{-10^20, 3}]", "{-100000000000000000000, 3, 1, 0, 0, 0.0}");
    audit_chk("UnixTime#2", "UnixTime[{10^20}]", "3155695199999999937832780800");
    audit_chk("UnixTime#2", "UnixTime[{10^20}] == AbsoluteTime[{10^20}] - 2208988800", "True");
}

int main(void) {
    audit_init();
    TEST(test_datelist_far_past);
    TEST(test_unixtime_int64_edge);
    TEST(test_s14_unixtime_exact);
    TEST(test_wave3_arity_options);
    TEST(test_wave4_dates);
    TEST(test_wave5_g7_time_and_date);
    TEST(test_wave6_g1_time_and_date);
    return audit_finish("audit_time_and_date");
}
