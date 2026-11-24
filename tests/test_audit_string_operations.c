/* Builtin-audit regressions: string operations (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* StringRiffle#1 (medium): "StringRiffle segfaults on a deeply nested list
 * (unbounded recursion)".  Both walks are iterative now. */
static void test_stringriffle_deep_nesting(void) {
    audit_run("auditDeep = Nest[List, \"a\", 100000]");
    audit_chk("StringRiffle#1", "StringLength[StringRiffle[auditDeep]]", "1");
    audit_true("StringRiffle#1", "StringRiffle[auditDeep] === \"a\"");
    audit_run("Clear[auditDeep]");
    /* Deep and wide: one separator per level, the default scheme. */
    audit_true("StringRiffle#1",
               "StringRiffle[Nest[{#, #} &, \"a\", 3]] === "
               "\"a a\\na a\\n\\na a\\na a\"");
    /* Every documented form is unchanged. */
    audit_true("StringRiffle#1", "StringRiffle[{\"a\", \"b\", \"c\"}] === \"a b c\"");
    audit_true("StringRiffle#1", "StringRiffle[{{1, 2}, {3, 4}}] === \"1 2\\n3 4\"");
    audit_true("StringRiffle#1",
               "StringRiffle[{\"a\", \"b\"}, {\"(\", \", \", \")\"}] === \"(a, b)\"");
    audit_true("StringRiffle#1",
               "StringRiffle[{{a, b}, {c}}, {\"[\", \";\", \"]\"}, \"-\"] === \"[a-b;c]\"");
    audit_true("StringRiffle#1", "StringRiffle[{}] === \"\"");
    audit_true("StringRiffle#1", "StringRiffle[{{}, {x}}] === \"\\nx\"");
}

/* FromCharacterCode#1 (high): "Code point 0 (NUL) truncates the string:
 * FromCharacterCode[0] is the empty string and {65,0,66} gives A".  U+0000 is
 * held as C0 80 internally (src/strings/utf8.h). */
static void test_fromcharactercode_nul(void) {
    audit_chk("FromCharacterCode#1", "StringLength[FromCharacterCode[{65, 0, 66}]]", "3");
    audit_chk("FromCharacterCode#1", "ToCharacterCode[FromCharacterCode[{65, 0, 66}]]", "{65, 0, 66}");
    audit_chk("FromCharacterCode#1", "StringLength[FromCharacterCode[0]]", "1");
    audit_chk("FromCharacterCode#1", "ToCharacterCode[FromCharacterCode[0]]", "{0}");
    /* the external UTF-8 byte form carries a 0 byte, and round-trips */
    audit_chk("FromCharacterCode#1", "ToCharacterCode[FromCharacterCode[{65, 0, 66}], \"UTF8\"]", "{65, 0, 66}");
    audit_true("FromCharacterCode#1",
               "FromCharacterCode[{65, 0, 66}, \"UTF8\"] === FromCharacterCode[{65, 0, 66}]");
    audit_true("FromCharacterCode#1",
               "StringReverse[FromCharacterCode[{65, 0, 66}]] === FromCharacterCode[{66, 0, 65}]");
    audit_true("FromCharacterCode#1",
               "StringTake[FromCharacterCode[{65, 0, 66}], {2}] === FromCharacterCode[0]");
    audit_chk("FromCharacterCode#1", "StringPosition[FromCharacterCode[{65, 0, 66}], \"B\"]", "{{3, 3}}");
    audit_chk("FromCharacterCode#1", "ToCharacterCode[FromCharacterCode[{0, 233, 0}]]", "{0, 233, 0}");
    /* ordinary text unchanged */
    audit_chk("FromCharacterCode#1", "ToCharacterCode[FromCharacterCode[{104, 233, 108}]]", "{104, 233, 108}");
    audit_chk("FromCharacterCode#1", "ToCharacterCode[FromCharacterCode[{104, 233}], \"UTF8\"]", "{104, 195, 169}");
}

/* StringReplaceList#1 (high): "StringReplaceList with n = Infinity returns {}
 * instead of every single replacement". */
static void test_stringreplacelist_infinity(void) {
    audit_true("StringReplaceList#1",
               "StringReplaceList[\"aaa\", \"a\" -> \"b\", Infinity] === {\"baa\", \"aba\", \"aab\"}");
    audit_true("StringReplaceList#1",
               "StringReplaceList[\"aaa\", \"a\" -> \"b\", Infinity] === StringReplaceList[\"aaa\", \"a\" -> \"b\"]");
    audit_true("StringReplaceList#1", "StringReplaceList[\"aaa\", \"a\" -> \"b\", 0] === {}");
    audit_true("StringReplaceList#1", "StringReplaceList[\"aaa\", \"a\" -> \"b\", 2] === {\"baa\", \"aba\"}");
    audit_true("StringReplaceList#1", "StringReplace[\"aaa\", \"a\" -> \"b\", Infinity] === \"bbb\"");
    /* n = 0 replaces nothing in StringReplace either */
    audit_true("StringReplaceList#1", "StringReplace[\"aaa\", \"a\" -> \"b\", 0] === \"aaa\"");
    audit_true("StringReplaceList#1", "StringReplace[\"aaa\", \"a\" -> \"b\", 2] === \"bba\"");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Decapitalize#1", "Quiet[Check[Decapitalize[\"A\", \"B\"], fired, Decapitalize::argx]]", "fired");
    audit_chk("Decapitalize#1", "Quiet[Decapitalize[\"A\", \"B\"]]", "Decapitalize[\"A\", \"B\"]");
    audit_chk("LetterQ#1", "Quiet[Check[LetterQ[\"a\", \"b\"], fired, LetterQ::argx]]", "fired");
    audit_chk("LetterQ#1", "Quiet[LetterQ[\"a\", \"b\"]]", "LetterQ[\"a\", \"b\"]");
    audit_chk("StringContainsQ#2", "Quiet[Check[StringContainsQ[\"abc\", \"b\", Foo -> 1], fired, StringContainsQ::optx]]", "fired");
    audit_chk("StringContainsQ#2", "Quiet[StringContainsQ[\"abc\", \"b\", Foo -> 1]]", "StringContainsQ[\"abc\", \"b\", Foo -> 1]");
    audit_chk("StringCount#2", "Quiet[Check[StringCount[\"abc\", \"a\", Foo -> 1], fired, StringCount::optx]]", "fired");
    audit_chk("StringCount#2", "Quiet[StringCount[\"abc\", \"a\", Foo -> 1]]", "StringCount[\"abc\", \"a\", Foo -> 1]");
    audit_chk("StringDelete#1", "Quiet[Check[StringDelete[], fired, StringDelete::argt]]", "fired");
    audit_chk("StringDelete#1", "Quiet[StringDelete[]]", "StringDelete[]");
    audit_chk("StringLength#2", "Quiet[Check[StringLength[\"a\", \"b\"], fired, StringLength::argx]]", "fired");
    audit_chk("StringLength#2", "Quiet[StringLength[\"a\", \"b\"]]", "StringLength[\"a\", \"b\"]");
    audit_chk("StringMatchQ#3", "Quiet[Check[StringMatchQ[\"abc\", \"abc\", Foo -> 1], fired, StringMatchQ::optx]]", "fired");
    audit_chk("StringMatchQ#3", "Quiet[StringMatchQ[\"abc\", \"abc\", Foo -> 1]]", "StringMatchQ[\"abc\", \"abc\", Foo -> 1]");
    audit_chk("StringContainsQ#2", "StringContainsQ[\"abc\", \"B\", IgnoreCase -> True]", "True");
}

int main(void) {
    audit_init();
    TEST(test_stringriffle_deep_nesting);
    TEST(test_fromcharactercode_nul);
    TEST(test_stringreplacelist_infinity);
    TEST(test_wave3_arity_options);
    return audit_finish("audit_string_operations");
}
