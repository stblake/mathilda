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

static void test_wave4_strings(void) {
    audit_chk("StringSplit#1", "StringSplit[\"abc\", \"\", 2]", "{\"a\", \"bc\"}");
    audit_chk("StringPosition#1", "StringPosition[\"aaa\",\"a\",0]", "{}");
    audit_chk("StringMatchQ#1", "StringMatchQ[\"abc\", \"abc\", SpellingCorrection -> False]", "True");
    audit_chk("StringMatchQ#1", "Options[StringMatchQ]", "{IgnoreCase -> False, SpellingCorrection -> False}");
    audit_chk("Alphabet#2", "Alphabet[\"English\"] === Alphabet[]", "True");
    audit_chk("Alphabet#2", "Length[Alphabet[\"Greek\"]]", "24");
    audit_chk("Alphabet#2", "Length[Alphabet[\"Russian\"]]", "33");
    audit_chk("Alphabet#2", "ToCharacterCode[StringJoin[Alphabet[\"Greek\"][[{1, 18, 24}]]]]", "{945, 963, 969}");
    audit_chk("Capitalize#1", "Capitalize[\"the quick brown fox\",\"AllWords\"]", "\"The Quick Brown Fox\"");
    audit_chk("Capitalize#1", "Capitalize[\"the lord of the rings\",\"TitleCase\"]", "\"The Lord of the Rings\"");
    audit_chk("StringExpression#1", "\"a\" ~~ \"b\"", "\"ab\"");
    audit_chk("StringExpression#1", "\"a\" ~~ \"b\" ~~ x ~~ \"c\" ~~ \"d\"", "\"ab\"~~x~~\"cd\"");
    audit_chk("StringExpression#1", "StringExpression[\"a\"]", "\"a\"");
    audit_chk("StringExpression#1", "StringExpression[\"a\", \"b\"] // StringLength", "2");
    audit_chk("StringExtract#1", "StringExtract[\"a b c\", 2;;5]", "{\"b\", \"c\", Missing[\"PartAbsent\", 4], Missing[\"PartAbsent\", 5]}");
    audit_chk("RegularExpression#1", "Check[StringContainsQ[StringJoin[Table[\"a\", {30}]] <> \"b\", RegularExpression[\"(a+)+$\"]], \"fired\", RegularExpression::limit]", "\"fired\"");
}

/* Wave 5 G7 (arith-misc): low-severity findings. Message counts are read
 * off $MessageList inside a Block (each diagnostic once per evaluation). */
static void test_wave5_g7_string_operations(void) {
    audit_chk("StringJoin#1", "Block[{$MessageList = {}}, StringJoin[\"a\", StringJoin[\"b\", x]]; Length[$MessageList]]", "2");
    audit_chk("StringJoin#2", "Quiet[Check[StringJoin[\"a\", x], fired, StringJoin::string]]", "fired");
    audit_chk("StringSplit#2", "Quiet[Check[StringSplit[x], fired, StringSplit::strse]]", "fired");
    audit_chk("StringSplit#2", "Quiet[Check[StringSplit[123, \",\"], fired, StringSplit::strse]]", "fired");
    audit_chk("StringSplit#2", "Quiet[StringSplit[{\"a,b\", 5}, \",\"]]", "StringSplit[{\"a,b\", 5}, \",\"]");
    audit_chk("Capitalize#2", "Quiet[Check[Capitalize[], fired, Capitalize::argt]]", "fired");
    audit_chk("Capitalize#2", "Quiet[Check[Capitalize[\"a\", \"b\"], fired, Capitalize::cmode]]", "fired");
    audit_chk("CharacterRange#3", "Quiet[Check[CharacterRange[], fired, CharacterRange::argrx]]", "fired");
    audit_chk("CharacterRange#3", "Quiet[Check[CharacterRange[\"a\"], fired, CharacterRange::argr]]", "fired");
    audit_chk("DigitQ#1", "Quiet[Check[DigitQ[], fired, DigitQ::argx]]", "fired");
    audit_chk("DigitQ#1", "Quiet[Check[DigitQ[\"1\", \"2\"], fired, DigitQ::argx]]", "fired");
    audit_chk("LowerCaseQ#1", "Quiet[Check[LowerCaseQ[], fired, LowerCaseQ::argx]]", "fired");
    audit_chk("LowerCaseQ#1", "Quiet[Check[LowerCaseQ[\"a\", \"b\"], fired, LowerCaseQ::argx]]", "fired");
    audit_chk("ToLowerCase#1", "Quiet[Check[ToLowerCase[], fired, ToLowerCase::argx]]", "fired");
    audit_chk("ToLowerCase#1", "Quiet[Check[ToLowerCase[\"A\", \"B\"], fired, ToLowerCase::argx]]", "fired");
    audit_chk("ToCharacterCode#3", "Quiet[Check[ToCharacterCode[], fired, ToCharacterCode::argt]]", "fired");
    audit_chk("ToCharacterCode#3", "Quiet[Check[ToCharacterCode[\"a\", \"b\", \"c\"], fired, ToCharacterCode::argt]]", "fired");
    /* One message for a bad regex, RegularExpression::regex at the offset in the
     * user's pattern; the consumer does not repeat it (Wave 5 G8, FileNames#2). */
    audit_chk("RegularExpression#2", "Block[{$MessageList = {}}, StringMatchQ[\"a\", RegularExpression[\"(\"]]; $MessageList]", "{RegularExpression::regex}");
    audit_chk("RegularExpression#2", "Quiet[Check[StringMatchQ[\"a\", RegularExpression[\"(\"]], fired, RegularExpression::regex]]", "fired");
    audit_chk("StringMatchQ#4", "Block[{$MessageList = {}}, StringMatchQ[\"abc\", RegularExpression[\"a(b\"]]; $MessageList]", "{RegularExpression::regex}");
    audit_chk("StringEndsQ#2", "Quiet[Check[StringEndsQ[\"abc\", RegularExpression[\"[\"]], fired, StringEndsQ::regex]]", "fired");
    audit_chk("StringStartsQ#2", "Quiet[Check[StringStartsQ[\"abc\", RegularExpression[\"(\"]], fired, StringStartsQ::regex]]", "fired");
    audit_chk("StringCases#1", "Quiet[Check[StringCases[1, \"b\"], fired, StringCases::strse]]", "fired");
    audit_chk("StringContainsQ#1", "Quiet[Check[StringContainsQ[5, \"a\"], fired, StringContainsQ::strse]]", "fired");
    audit_chk("StringEndsQ#1", "Quiet[Check[StringEndsQ[5, \"a\"], fired, StringEndsQ::strse]]", "fired");
    audit_chk("StringEndsQ#1", "Quiet[StringEndsQ[{\"a\", 5}, \"a\"]]", "StringEndsQ[{\"a\", 5}, \"a\"]");
    audit_chk("StringFreeQ#1", "Quiet[Check[StringFreeQ[x, \"a\"], fired, StringFreeQ::strse]]", "fired");
    audit_chk("StringStartsQ#1", "Quiet[Check[StringStartsQ[x, \"a\"], fired, StringStartsQ::strse]]", "fired");
    audit_chk("StringStartsQ#1", "StringStartsQ[{\"ab\", \"ba\"}, \"a\"]", "{True, False}");
    audit_chk("StringMatchQ#2", "Quiet[Check[StringMatchQ[5, \"a\"], fired, StringMatchQ::strse]]", "fired");
    audit_chk("StringMatchQ#2", "Quiet[StringMatchQ[{\"a\", 5}, \"a\"]]", "StringMatchQ[{\"a\", 5}, \"a\"]");
    audit_chk("StringMatchQ#2", "Quiet[StringMatchQ[{{\"a\"}}, \"a\"]]", "StringMatchQ[{{\"a\"}}, \"a\"]");
    audit_chk("StringCount#1", "StringCount[\"abc\", {}]", "0");
    audit_chk("StringCount#1", "StringCount[{\"abc\", \"d\"}, {}]", "{0, 0}");
    audit_chk("StringDrop#1", "Block[{$MessageList = {}}, StringDrop[{\"abcd\", \"x\"}, 2]; Length[$MessageList]]", "1");
    audit_chk("StringInsert#1", "Block[{$MessageList = {}}, StringInsert[{\"abc\", \"d\"}, \"X\", 3]; Length[$MessageList]]", "1");
    audit_chk("StringPart#1", "Block[{$MessageList = {}}, StringPart[{\"abc\", \"de\"}, 3]; Length[$MessageList]]", "1");
    audit_chk("StringExtract#2", "StringExtract[\"a b c\", {1, x}]", "StringExtract[\"a b c\", {1, x}]");
    audit_chk("StringExtract#2", "StringExtract[\"a b c\", {1, 3}]", "{\"a\", \"c\"}");
    audit_chk("StringPosition#2", "Quiet[Check[StringPosition[\"aaa\", \"a\", -1], fired, StringPosition::innf]]", "fired");
    audit_chk("StringPosition#2", "Quiet[StringPosition[\"aaa\", \"a\", -1]]", "StringPosition[\"aaa\", \"a\", -1]");
    audit_chk("StringRepeat#1", "StringRepeat[\"ab\", 10^30, 5]", "\"ababa\"");
    audit_chk("StringRepeat#1", "StringRepeat[\"ab\", 3, 10^30]", "\"ababab\"");
    audit_chk("StringReplaceList#2", "StringReplaceList[\"aaa\", \"a\" -> \"b\", 10^30]", "{\"baa\", \"aba\", \"aab\"}");
    audit_chk("StringReplacePart#1", "StringFreeQ[Information[\"StringReplacePart\"], \"{s1, s2, ...}\"]", "True");
    audit_chk("StringReplacePart#2", "StringReplacePart[\"abc\", \"X\"]", "StringReplacePart[\"abc\", \"X\"]");
    audit_chk("StringReplacePart#2", "Head[StringReplacePart[\"X\", {1, 1}]]", "StringReplacePart");
    audit_chk("StringReplacePart#2", "StringReplacePart[\"X\", {1, 1}][\"abc\"]", "\"Xbc\"");
    audit_chk("StringRotateLeft#1", "StringRotateLeft[\"abc\", 10^20]", "\"bca\"");
    audit_chk("StringRotateRight#1", "StringRotateRight[\"abcdef\", 10^20]", "\"cdefab\"");
    audit_chk("StringRotateRight#2", "StringRotateRight[\"abcdef\", -9223372036854775808]", "\"cdefab\"");
    audit_chk("StringTrim#1", "Quiet[Check[StringTrim[{\" a \", 1}], fired, StringTrim::strse]]", "fired");
}

/* Wave 6 G1: surrogate codes are characters; the StringJoin and StringTrim
 * messages the spec pages show. */
static void test_wave6_g1_strings(void) {
    audit_chk("ToCharacterCode#2", "ToCharacterCode[\"\\:d83d\"]", "{55357}");
    audit_chk("ToCharacterCode#2", "StringLength[\"\\:d83d\\:de00\"]", "2");
    audit_chk("CharacterRange#2", "Take[ToCharacterCode[CharacterRange[55295, 57344]], 3]", "{{55295}, {55296}, {55297}}");
    audit_chk("CharacterRange#2", "Count[CharacterRange[55295, 57344], \"?\"]", "0");
    audit_chk("CharacterRange#2", "Length[CharacterRange[55295, 57344]]", "2050");
    audit_chk("CharacterRange#2", "ToCharacterCode[CharacterRange[55295, 55297]]", "{{55295}, {55296}, {55297}}");
    audit_chk("CharacterRange#2", "ToCharacterCode[FromCharacterCode[55357]]", "{55357}");
    audit_chk("StringJoin#2", "Check[StringJoin[\"a\", x], fired, StringJoin::string]", "fired");
    audit_chk("StringJoin#2", "Quiet[StringJoin[\"a\", x]]", "StringJoin[\"a\", x]");
    audit_chk("StringTrim#1", "Check[StringTrim[{\" a \", 1}], fired, StringTrim::strse]", "fired");
    audit_chk("StringTrim#1", "Quiet[StringTrim[{\" a \", 1}]]", "StringTrim[{\" a \", 1}]");
}

int main(void) {
    audit_init();
    TEST(test_stringriffle_deep_nesting);
    TEST(test_fromcharactercode_nul);
    TEST(test_stringreplacelist_infinity);
    TEST(test_wave3_arity_options);
    TEST(test_wave4_strings);
    TEST(test_wave5_g7_string_operations);
    TEST(test_wave6_g1_strings);
    return audit_finish("audit_string_operations");
}
