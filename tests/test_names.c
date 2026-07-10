/*
 * test_names.c - unit tests for the Names builtin.
 *
 * Names enumerates symbol-table names by string pattern (with the * and @
 * metacharacters), a RegularExpression[...], or a list of such patterns; the
 * result is a canonically sorted List of Strings.  Assertions favour stable
 * shapes (===, MemberQ, Length) over pinning exact name lists, since the set
 * of registered builtins grows over time.
 *
 * assert_eval_eq's assert() is compiled out under -DNDEBUG, but it still
 * prints "FAIL:" to stderr -- CI greps for that.
 */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

/* ---- exact literal string: matches just that symbol ---- */
void test_names_exact_literal() {
    assert_eval_eq("Names[\"Sin\"]", "{\"Sin\"}", 0);
}

/* ---- no match yields an empty list ---- */
void test_names_no_match() {
    assert_eval_eq("Names[\"ZzzzNoSuchSymbol*\"]", "{}", 0);
}

/* ---- '*' matches zero or more characters ---- */
void test_names_star_prefix() {
    /* "List" itself (zero chars) and longer names are all included. */
    assert_eval_eq("MemberQ[Names[\"List*\"], \"List\"]", "True", 0);
    assert_eval_eq("MemberQ[Names[\"List*\"], \"ListQ\"]", "True", 0);
    /* A non-List name is excluded. */
    assert_eval_eq("MemberQ[Names[\"List*\"], \"Sin\"]", "False", 0);
}

/* ---- '*' in the middle and a bare "*" matching everything ---- */
void test_names_star_forms() {
    assert_eval_eq("MemberQ[Names[\"*Plot\"], \"ListPlot\"]", "True", 0);
    assert_eval_eq("Names[\"*\"] === Names[]", "True", 0);
}

/* ---- '@' matches one or more NON-uppercase characters ---- */
void test_names_at_metacharacter() {
    /* "Ar" + one-or-more lowercase: Arg matches (g), ... */
    assert_eval_eq("MemberQ[Names[\"Ar@\"], \"Arg\"]", "True", 0);
    /* ...but ArcSin does NOT: '@' cannot cross the uppercase S. */
    assert_eval_eq("MemberQ[Names[\"Ar@\"], \"ArcSin\"]", "False", 0);
    /* '@' requires at least one character: "Sin@" matches nothing named "Sin". */
    assert_eval_eq("MemberQ[Names[\"Sin@\"], \"Sin\"]", "False", 0);
}

/* ---- a List argument is a set of alternative patterns ---- */
void test_names_list_of_patterns() {
    assert_eval_eq("Names[{\"Sin\", \"Cos\"}]", "{\"Cos\", \"Sin\"}", 0);
    /* Overlapping alternatives do not duplicate a name. */
    assert_eval_eq("Names[{\"Sin\", \"Sin\"}]", "{\"Sin\"}", 0);
    /* An empty list matches nothing. */
    assert_eval_eq("Names[{}]", "{}", 0);
}

/* ---- output is canonically sorted (=== Sort of itself) ---- */
void test_names_sorted() {
    assert_eval_eq("Names[\"S*\"] === Sort[Names[\"S*\"]]", "True", 0);
    assert_eval_eq("Names[] === Sort[Names[]]", "True", 0);
}

/* ---- Names[] lists the whole namespace ---- */
void test_names_all() {
    assert_eval_eq("Length[Names[]] > 0", "True", 0);
    assert_eval_eq("MemberQ[Names[], \"Plus\"]", "True", 0);
    assert_eval_eq("MemberQ[Names[], \"Names\"]", "True", 0);
}

/* ---- every element of the result is a String ---- */
void test_names_returns_strings() {
    assert_eval_eq("Head[Names[\"Sin\"][[1]]]", "String", 0);
    /* Every element is a String -> the set of heads is exactly {String}. */
    assert_eval_eq("Union[Head /@ Names[\"List*\"]]", "{String}", 0);
    assert_eval_eq("Union[Head /@ Names[]]", "{String}", 0);
}

/* ---- RegularExpression pattern, matched against the WHOLE name (anchored) ---- */
void test_names_regular_expression() {
    /* "Si." is exactly three chars -> {Sin} (anchored, not a substring match). */
    assert_eval_eq("Names[RegularExpression[\"Si.\"]]", "{\"Sin\"}", 0);
    /* Anchoring: a bare "in" must NOT match "Sin". */
    assert_eval_eq("MemberQ[Names[RegularExpression[\"in\"]], \"Sin\"]", "False", 0);
    /* Alternation inside the regex still works. */
    assert_eval_eq("Names[RegularExpression[\"Sin|Cos\"]]", "{\"Cos\", \"Sin\"}", 0);
}

/* ---- arity / bad argument: stays unevaluated ---- */
void test_names_unevaluated() {
    assert_eval_eq("Names[5]", "Names[5]", 0);
    assert_eval_eq("Names[\"Sin\", \"Cos\"]", "Names[\"Sin\", \"Cos\"]", 0);
}

/* ---- Names is Protected ---- */
void test_names_protected() {
    /* Redefining a Protected symbol is refused; the builtin still works. */
    assert_eval_eq("Names[\"Sin\"]", "{\"Sin\"}", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_names_exact_literal);
    TEST(test_names_no_match);
    TEST(test_names_star_prefix);
    TEST(test_names_star_forms);
    TEST(test_names_at_metacharacter);
    TEST(test_names_list_of_patterns);
    TEST(test_names_sorted);
    TEST(test_names_all);
    TEST(test_names_returns_strings);
    TEST(test_names_regular_expression);
    TEST(test_names_unevaluated);
    TEST(test_names_protected);

    printf("All Names tests passed.\n");
    return 0;
}
