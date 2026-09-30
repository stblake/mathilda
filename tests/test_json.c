/* End-to-end tests for JSON: ImportString / ExportString with "RawJSON" and
 * "JSON", Import / Export of .json files, and round trips (src/json.c).
 *
 * Expected values were produced by Mathematica 15 (wolframscript) on the same
 * input.  Where Mathematica's own error text or message position differs, the
 * test checks only the $Failed result. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "json.h"
#include <stdio.h>

static int g_checks = 0;
#define Q(in, out) do { assert_eval_eq((in), (out), 0); g_checks++; } while (0)

/* Quiet[...] keeps the expected error messages out of the test log. */
#define QF(in) Q("Quiet[" in "]", "$Failed")

/* ---- import: values ---- */

static void test_import_values(void) {
    Q("ImportString[\"{\\\"a\\\": 1, \\\"b\\\": [1, 2.5, true, false, null], \\\"c\\\": {\\\"d\\\": \\\"x\\\"}}\", \"RawJSON\"]",
      "<|\"a\" -> 1, \"b\" -> {1, 2.5, True, False, Null}, \"c\" -> <|\"d\" -> \"x\"|>|>");
    Q("ImportString[\"[1,2,3]\", \"RawJSON\"]", "{1, 2, 3}");
    Q("ImportString[\"42\", \"RawJSON\"]", "42");
    Q("ImportString[\"-7\", \"RawJSON\"]", "-7");
    Q("ImportString[\"\\\"hi\\\"\", \"RawJSON\"]", "\"hi\"");
    Q("ImportString[\"null\", \"RawJSON\"]", "Null");
    Q("ImportString[\"true\", \"RawJSON\"]", "True");
    Q("ImportString[\"{}\", \"RawJSON\"]", "<||>");
    Q("ImportString[\"[]\", \"RawJSON\"]", "{}");
    Q("ImportString[\"  [ 1 , { \\\"b\\\" : [ ] } ]\\n \", \"RawJSON\"]", "{1, <|\"b\" -> {}|>}");
    /* duplicate keys: last value wins, first position kept */
    Q("ImportString[\"{\\\"a\\\":1,\\\"b\\\":0,\\\"a\\\":2}\", \"RawJSON\"]", "<|\"a\" -> 2, \"b\" -> 0|>");
    /* numbers */
    Q("ImportString[\"1.0\", \"RawJSON\"]", "1.0");
    Q("ImportString[\"0.1\", \"RawJSON\"]", "0.1");
    Q("ImportString[\"1.5e3\", \"RawJSON\"] === 1500.", "True");
    Q("ImportString[\"-3e2\", \"RawJSON\"]", "-300");        /* exponent, no fraction: exact */
    Q("ImportString[\"1E-2\", \"RawJSON\"]", "1/100");
    Q("ImportString[\"2.5E+1\", \"RawJSON\"] === 25.", "True");
    Q("ImportString[\"-0\", \"RawJSON\"]", "0");
    Q("ImportString[\"01\", \"RawJSON\"]", "1");             /* leading zero accepted, as Mathematica */
    Q("ImportString[\"123456789012345678901234567890\", \"RawJSON\"]",
      "123456789012345678901234567890");
    Q("ImportString[\"9223372036854775808\", \"RawJSON\"]", "9223372036854775808");
    Q("ImportString[\"1e30\", \"RawJSON\"]", "1000000000000000000000000000000");
    Q("NumberQ[ImportString[\"1.5e400\", \"RawJSON\"]]", "True");
    /* strings and escapes */
    Q("ImportString[\"\\\"a\\\\n\\\\t\\\\\\\"\\\\\\\\\\\\/b\\\"\", \"RawJSON\"] === \"a\\n\\t\\\"\\\\/b\"", "True");
    Q("ImportString[\"\\\"\\\\u00e9\\\"\", \"RawJSON\"] === \"\xc3\xa9\"", "True");
    Q("ImportString[\"\\\"\\\\ud83d\\\\ude00\\\"\", \"RawJSON\"] === \"\xf0\x9f\x98\x80\"", "True");
    Q("ImportString[\"\\\"\\\\u20AC\\\"\", \"RawJSON\"] === \"\xe2\x82\xac\"", "True");
    Q("ImportString[\"\\\"\xc3\xa9\\\"\", \"RawJSON\"] === \"\xc3\xa9\"", "True");   /* raw UTF-8 */
    Q("StringLength[ImportString[\"\\\"\\\\u0041\\\\u0042\\\"\", \"RawJSON\"]]", "2");
}

/* ---- import: the "JSON" (rule list) form and plain text ---- */

static void test_import_json_rules(void) {
    Q("ImportString[\"{\\\"a\\\": 1, \\\"b\\\": [1, true, null], \\\"c\\\": {\\\"d\\\": \\\"x\\\"}}\", \"JSON\"]",
      "{\"a\" -> 1, \"b\" -> {1, True, Null}, \"c\" -> {\"d\" -> \"x\"}}");
    Q("ImportString[\"[{\\\"a\\\":1}, 2]\", \"JSON\"]", "{{\"a\" -> 1}, 2}");
    Q("ImportString[\"{}\", \"JSON\"]", "{}");
    Q("ImportString[\"x\"]", "\"x\"");
}

/* ---- import: errors ---- */

static void test_import_errors(void) {
    QF("ImportString[\"[1,]\", \"RawJSON\"]");
    QF("ImportString[\"{\\\"a\\\" 1}\", \"RawJSON\"]");
    QF("ImportString[\"tru\", \"RawJSON\"]");
    QF("ImportString[\"\", \"RawJSON\"]");
    QF("ImportString[\"   \", \"RawJSON\"]");
    QF("ImportString[\"[1] x\", \"RawJSON\"]");
    QF("ImportString[\"[NaN]\", \"RawJSON\"]");
    QF("ImportString[\"\\\"\\\\x\\\"\", \"RawJSON\"]");
    QF("ImportString[\"{\\\"a\\\":1,}\", \"RawJSON\"]");
    QF("ImportString[\"{1:2}\", \"RawJSON\"]");
    QF("ImportString[\"\\\"abc\", \"RawJSON\"]");
    QF("ImportString[\"\\\"\\\\ud83d\\\"\", \"RawJSON\"]");       /* stray high surrogate */
    QF("ImportString[\"\\\"\\\\ude00\\\"\", \"RawJSON\"]");       /* stray low surrogate */
    QF("ImportString[\"\\\"\\\\u12\\\"\", \"RawJSON\"]");
    QF("ImportString[\"[1 2]\", \"RawJSON\"]");
    QF("ImportString[\"[.5]\", \"RawJSON\"]");
    QF("ImportString[\"[+1]\", \"RawJSON\"]");
    QF("ImportString[\"[1e]\", \"RawJSON\"]");
    QF("ImportString[\"[1,]\", \"JSON\"]");
    /* the error is a message: Check sees it */
    Q("Check[ImportString[\"[\", \"RawJSON\"], caught]", "caught");
}

/* ---- export ---- */

/* ExportString results are compared with === against a Mathilda string
 * literal (the test printer does not escape quotes inside strings). */
#define QS(expr, lit) Q(expr " === " lit, "True")

static void test_export(void) {
    QS("ExportString[<|\"a\" -> 1, \"b\" -> {1, 2.5}, \"c\" -> <|\"x\" -> False|>|>, \"RawJSON\", \"Compact\" -> True]",
       "\"{\\\"a\\\":1,\\\"b\\\":[1,2.5],\\\"c\\\":{\\\"x\\\":false}}\"");
    QS("ExportString[<|\"a\" -> 1, \"b\" -> {1, 2}|>, \"RawJSON\"]",
       "\"{\\n\\t\\\"a\\\":1,\\n\\t\\\"b\\\":[\\n\\t\\t1,\\n\\t\\t2\\n\\t]\\n}\"");
    QS("ExportString[<|\"a\" -> 1, \"b\" -> {1, 2}|>, \"RawJSON\", \"Compact\" -> False]",
       "\"{\\n\\t\\\"a\\\":1,\\n\\t\\\"b\\\":[\\n\\t\\t1,\\n\\t\\t2\\n\\t]\\n}\"");
    QS("ExportString[{<|\"a\" -> 1|>, {}}, \"RawJSON\"]",
       "\"[\\n\\t{\\n\\t\\t\\\"a\\\":1\\n\\t},\\n\\t[]\\n]\"");
    QS("ExportString[<|\"a\" -> {}, \"b\" -> <||>|>, \"RawJSON\"]",
       "\"{\\n\\t\\\"a\\\":[],\\n\\t\\\"b\\\":{}\\n}\"");
    QS("ExportString[{}, \"RawJSON\"]", "\"[]\"");
    QS("ExportString[<||>, \"RawJSON\"]", "\"{}\"");
    QS("ExportString[{True, False, Null}, \"RawJSON\", \"Compact\" -> True]", "\"[true,false,null]\"");
    QS("ExportString[-5, \"RawJSON\"]", "\"-5\"");
    QS("ExportString[10^30, \"RawJSON\"]", "\"1000000000000000000000000000000\"");
    QS("ExportString[1/2, \"RawJSON\"]", "\"0.5\"");
    QS("ExportString[1/3, \"RawJSON\"]", "\"0.3333333333333333\"");
    /* reals: Mathematica's JSON exponent style */
    QS("ExportString[1.0, \"RawJSON\"]", "\"1.0\"");
    QS("ExportString[0., \"RawJSON\"]", "\"0.0\"");
    QS("ExportString[0.12, \"RawJSON\"]", "\"0.12\"");
    QS("ExportString[-0.5, \"RawJSON\"]", "\"-0.5\"");
    QS("ExportString[25., \"RawJSON\"]", "\"2.5e1\"");
    QS("ExportString[0.05, \"RawJSON\"]", "\"5.0e-2\"");
    QS("ExportString[100000., \"RawJSON\"]", "\"1.0e5\"");
    QS("ExportString[123456789.123, \"RawJSON\"]", "\"1.23456789123e8\"");
    QS("ExportString[1.*^-7, \"RawJSON\"]", "\"1.0e-7\"");
    QS("ExportString[2.^53, \"RawJSON\"]", "\"9.007199254740992e15\"");
    QS("ExportString[0.1 + 0.2, \"RawJSON\"]", "\"0.30000000000000004\"");
    /* strings: escapes, / as \/, control bytes as \u00XX, UTF-8 passes through */
    QS("ExportString[\"x\\ny\\\"z\", \"RawJSON\"]", "\"\\\"x\\\\ny\\\\\\\"z\\\"\"");
    QS("ExportString[\"a/b\", \"RawJSON\"]", "\"\\\"a\\\\/b\\\"\"");
    QS("ExportString[ImportString[\"\\\"\\\\u0001\\\"\", \"RawJSON\"], \"RawJSON\"]", "\"\\\"\\\\u0001\\\"\"");
    QS("ExportString[\"\xc3\xa9\", \"RawJSON\"]", "\"\\\"\xc3\xa9\\\"\"");
    /* "JSON" writes rule lists as objects */
    QS("ExportString[{\"a\" -> 1, \"b\" -> {1, 2}}, \"JSON\", \"Compact\" -> True]",
       "\"{\\\"a\\\":1,\\\"b\\\":[1,2]}\"");
    QS("ExportString[{{\"a\" -> 1}, {\"b\" -> 2}}, \"JSON\", \"Compact\" -> True]",
       "\"[{\\\"a\\\":1},{\\\"b\\\":2}]\"");
    QS("ExportString[<|\"a\" -> <|\"b\" -> 1|>|>, \"JSON\", \"Compact\" -> True]",
       "\"{\\\"a\\\":{\\\"b\\\":1}}\"");
}

static void test_export_errors(void) {
    QF("ExportString[x, \"RawJSON\"]");
    QF("ExportString[<|1 -> 2|>, \"RawJSON\"]");
    QF("ExportString[{f[1]}, \"RawJSON\"]");
    QF("ExportString[Pi, \"RawJSON\"]");
    QF("ExportString[{\"a\" -> 1}, \"RawJSON\"]");      /* rule lists are not RawJSON objects */
    QF("ExportString[{\"a\" -> 1, 2}, \"JSON\"]");
    QF("ExportString[{x -> 1}, \"JSON\"]");
    QF("ExportString[Missing[], \"RawJSON\"]");
}

/* ---- round trips ---- */

static void test_roundtrip(void) {
    Q("jr = <|\"n\" -> {1, -2, 3.25, 10^20, 0.1, 1.*^-300, 6.02*^23}, \"s\" -> \"q\\\"\\\\\\n\\t/\xc3\xbc\", "
      "\"b\" -> {True, False, Null}, \"o\" -> <|\"x\" -> <|\"y\" -> {}|>, \"e\" -> <||>|>|>; "
      "ImportString[ExportString[jr, \"RawJSON\"], \"RawJSON\"] === jr", "True");
    Q("ImportString[ExportString[jr, \"RawJSON\", \"Compact\" -> True], \"RawJSON\"] === jr", "True");
    Q("ImportString[ExportString[{\"a\" -> {1, {\"b\" -> 2}}}, \"JSON\"], \"JSON\"]",
      "{\"a\" -> {1, {\"b\" -> 2}}}");
    /* Every machine real survives Export -> Import exactly. */
    Q("SeedRandom[7]; rr = RandomReal[{-10^6, 10^6}, 200]; "
      "ImportString[ExportString[rr, \"RawJSON\"], \"RawJSON\"] === rr", "True");
    Q("rs = Table[10.^k * 1.2345678901234567, {k, -300, 300, 37}]; "
      "ImportString[ExportString[rs, \"RawJSON\"], \"RawJSON\"] === rs", "True");
    /* \u0000 is kept (as modified UTF-8) and written back as \u0000. */
    Q("ExportString[ImportString[\"\\\"a\\\\u0000b\\\"\", \"RawJSON\"], \"RawJSON\"] === \"\\\"a\\\\u0000b\\\"\"",
      "True");
    Q("ImportString[ExportString[Range[1000], \"RawJSON\"], \"RawJSON\"] === Range[1000]", "True");
}

/* ---- Import / Export files ---- */

static void test_files(void) {
    Q("Export[\"/tmp/mathilda_test_json.json\", <|\"a\" -> {1, 2}, \"b\" -> \"c\"|>]",
      "\"/tmp/mathilda_test_json.json\"");
    Q("Import[\"/tmp/mathilda_test_json.json\"]", "{\"a\" -> {1, 2}, \"b\" -> \"c\"}");
    Q("Import[\"/tmp/mathilda_test_json.json\", \"RawJSON\"]", "<|\"a\" -> {1, 2}, \"b\" -> \"c\"|>");
    Q("Export[\"/tmp/mathilda_test_json.txt\", {1, 2}, \"RawJSON\", \"Compact\" -> True]",
      "\"/tmp/mathilda_test_json.txt\"");
    Q("Import[\"/tmp/mathilda_test_json.txt\", \"RawJSON\"]", "{1, 2}");
    QF("Import[\"/tmp/mathilda_no_such_file.json\"]");
}

/* ---- the C API ---- */

static void test_c_api(void) {
    Expr* e = json_parse("{\"k\": [1, 2]}", 13, true, "Import");
    char* s;
    ASSERT(e != NULL);
    s = json_serialize(e, true, true);
    ASSERT_STR_EQ(s, "{\"k\":[1,2]}");
    free(s);
    expr_free(e);
    g_checks += 2;
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_import_values);
    TEST(test_import_json_rules);
    TEST(test_import_errors);
    TEST(test_export);
    TEST(test_export_errors);
    TEST(test_roundtrip);
    TEST(test_files);
    TEST(test_c_api);

    printf("All JSON tests passed (%d checks).\n", g_checks);
    return 0;
}
