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

int main(void) {
    audit_init();
    TEST(test_datelist_far_past);
    TEST(test_unixtime_int64_edge);
    return audit_finish("audit_time_and_date");
}
