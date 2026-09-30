/* End-to-end tests for Query and Dataset (src/assoc_query.c).
 *
 * Every expected value below was produced by Mathematica 15 (wolframscript)
 * on the same input, except where a comment notes a documented difference.
 *
 * Named slots: the record predicates are written #["a"] because `#a` is
 * still parsed as #1*a (the parser fix is a separate stream).  Once `#a`
 * parses as a named slot, the tests marked [#a] should be switched to that
 * spelling -- the expected outputs do not change. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "assoc_query.h"
#include <stdio.h>

static int g_checks = 0;
#define Q(in, out) do { assert_eval_eq((in), (out), 0); g_checks++; } while (0)

static void run(const char* src) {
    Expr* e = parse_expression(src);
    Expr* r;
    assert(e != NULL);
    r = evaluate(e);
    expr_free(e);
    if (r) expr_free(r);
}

#define REC1 "<|\"a\" -> 1, \"b\" -> \"x\", \"c\" -> 3.5|>"
#define REC2 "<|\"a\" -> 2, \"b\" -> \"y\", \"c\" -> 1.5|>"
#define REC3 "<|\"a\" -> 3, \"b\" -> \"z\", \"c\" -> 2.5|>"

static void setup(void) {
    run("qdata = {<|\"a\"->1,\"b\"->\"x\",\"c\"->3.5|>, <|\"a\"->2,\"b\"->\"y\",\"c\"->1.5|>, "
        "<|\"a\"->3,\"b\"->\"z\",\"c\"->2.5|>}");
    run("qmd = {<|\"a\"->1|>, <|\"b\"->2|>, <|\"a\"->5|>}");
    run("qaa = <|\"x\" -> <|\"b\" -> 5|>, \"y\" -> <|\"b\" -> 6|>|>");
    run("qds = Dataset[qdata]");
}

/* ---- part specifications (descending, level-consuming or mapping) ---- */

static void test_query_parts(void) {
    Q("Query[\"a\"][<|\"a\" -> 1, \"b\" -> 2|>]", "1");
    Q("Query[Key[\"a\"]][<|\"a\" -> 1|>]", "1");
    Q("Query[All, \"a\"][qdata]", "{1, 2, 3}");
    Q("Query[2, \"b\"][qdata]", "\"y\"");
    Q("Query[-1][qdata]", REC3);
    Q("Query[{1, 3}][qdata]", "{" REC1 ", " REC3 "}");
    Q("Query[{1, 3}, \"a\"][qdata]", "{1, 3}");
    Q("Query[1 ;; 2, \"a\"][qdata]", "{1, 2}");
    Q("Query[2 ;;][qdata]", "{" REC2 ", " REC3 "}");
    Q("Query[;; ;; 2, \"a\"][qdata]", "{1, 3}");
    Q("Query[All, {\"a\", \"b\"}][qdata]",
      "{<|\"a\" -> 1, \"b\" -> \"x\"|>, <|\"a\" -> 2, \"b\" -> \"y\"|>, <|\"a\" -> 3, \"b\" -> \"z\"|>}");
    Q("Query[{\"a\", \"b\"}][<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>]", "<|\"a\" -> 1, \"b\" -> 2|>");
    Q("Query[\"a\", \"b\"][<|\"a\" -> <|\"b\" -> 5|>|>]", "5");
    Q("Query[\"a\", 1][<|\"a\" -> {5, 6}|>]", "5");
    Q("Query[All, \"b\"][qaa]", "<|\"x\" -> 5, \"y\" -> 6|>");
    Q("Query[1][<|\"x\" -> 5, \"y\" -> 6|>]", "5");
    Q("Query[{1}][<|\"x\" -> 5, \"y\" -> 6|>]", "<|\"x\" -> 5|>");
    Q("Query[All][qdata] === qdata", "True");
    Q("Query[][qdata] === qdata", "True");
    Q("Query[All, All, f][{<|\"a\" -> 1, \"b\" -> 2|>}]", "{<|\"a\" -> f[1], \"b\" -> f[2]|>}");
}

/* ---- Missing propagation ---- */

static void test_query_missing(void) {
    Q("Query[\"z\"][<|\"a\" -> 1|>]", "Missing[\"KeyAbsent\", \"z\"]");
    Q("Query[5][qdata]", "Missing[\"PartAbsent\", 5]");
    Q("Query[All, \"a\"][{1, 2}]",
      "{Missing[\"PartInvalid\", \"a\"], Missing[\"PartInvalid\", \"a\"]}");
    Q("Query[\"a\"][{<|\"a\" -> 1|>}]", "Missing[\"PartInvalid\", \"a\"]");
    Q("Query[\"z\", \"b\"][<|\"a\" -> <|\"b\" -> 5|>|>]", "Missing[\"KeyAbsent\", \"z\"]");
    Q("Query[\"a\", \"z\"][<|\"a\" -> <|\"b\" -> 5|>|>]", "Missing[\"KeyAbsent\", \"z\"]");
    Q("Query[{1, 5}][{1, 2}]", "{1, Missing[\"PartAbsent\", 5]}");
    Q("Query[{\"a\", \"z\"}][<|\"a\" -> 1|>]", "<|\"a\" -> 1, \"z\" -> Missing[\"KeyAbsent\", \"z\"]|>");
    Q("Query[All, \"z\", f][{<|\"a\" -> 1|>}]", "{f[Missing[\"KeyAbsent\", \"z\"]]}");
    Q("Query[All, \"a\"][<|\"x\" -> 1|>]", "<|\"x\" -> Missing[\"PartInvalid\", \"a\"]|>");
    /* Aggregations drop Missing (MissingBehavior -> Automatic) ... */
    Q("Query[Total, \"a\"][qmd]", "6");
    Q("Query[Mean, \"a\"][qmd]", "3");
    Q("Query[Max, \"a\"][qmd]", "5");
    Q("Query[Min, \"a\"][qmd]", "1");
    Q("Query[Total, \"a\"][{<|\"b\" -> 1|>}]", "0");
    /* ... but Length, Counts and arbitrary functions see them. */
    Q("Query[Length, \"a\"][qmd]", "3");
    Q("Query[f, \"a\"][qmd]", "f[{1, Missing[\"KeyAbsent\", \"a\"], 5}]");
    Q("Query[Counts, \"a\"][qmd]", "<|1 -> 1, Missing[\"KeyAbsent\", \"a\"] -> 1, 5 -> 1|>");
}

/* ---- descending operators ---- */

static void test_query_descending(void) {
    Q("Query[Select[#[\"a\"] > 1 &]][qdata]", "{" REC2 ", " REC3 "}");              /* [#a] */
    Q("Query[Select[#[\"a\"] > 1 &], \"b\"][qdata]", "{\"y\", \"z\"}");                /* [#a] */
    Q("Query[SortBy[#[\"c\"] &]][qdata]", "{" REC2 ", " REC3 ", " REC1 "}");           /* [#a] */
    Q("Query[SortBy[#[\"c\"] &], \"a\"][qdata]", "{2, 3, 1}");                          /* [#a] */
    Q("Query[SortBy[\"c\"], \"a\"][qdata]", "{2, 3, 1}");
    Q("Query[MaximalBy[\"c\"], \"a\"][qdata]", "{1}");
    Q("Query[Select[EvenQ]][{1, 2, 3, 4}]", "{2, 4}");
    Q("Query[Select[EvenQ], #^2 &][{1, 2, 3, 4}]", "{4, 16}");
    Q("Query[Reverse, \"a\"][qdata]", "{3, 2, 1}");
    Q("Query[Select[#[\"b\"] > 5 &]][qaa]", "<|\"y\" -> <|\"b\" -> 6|>|>");               /* [#a] */
    Q("Query[GroupBy[#[\"a\"] > 1 &]][qdata]",                                        /* [#a] */
      "<|False -> {" REC1 "}, True -> {" REC2 ", " REC3 "}|>");
    Q("Query[GroupBy[#[\"a\"] > 1 &], Total, \"a\"][qdata]", "<|False -> 1, True -> 5|>"); /* [#a] */
    Q("Query[GroupBy[#[\"a\"] > 1 &], Length][qdata]", "<|False -> 1, True -> 2|>");   /* [#a] */
    Q("Query[GroupBy[\"b\"], All, \"a\"][{<|\"a\" -> 1, \"b\" -> 2|>, <|\"a\" -> 3, \"b\" -> 2|>}]",
      "<|2 -> {1, 3}|>");
    Q("Query[Values][<|\"a\" -> 1, \"b\" -> 2|>]", "{1, 2}");
    Q("Query[Values, \"a\"][<|\"p\" -> <|\"a\" -> 1|>, \"q\" -> <|\"a\" -> 2|>|>]", "{1, 2}");
    Q("Query[Keys][qdata]", "{{\"a\", \"b\", \"c\"}, {\"a\", \"b\", \"c\"}, {\"a\", \"b\", \"c\"}}");
    /* Descending Select, then an ascending Total at level 2: Total of each
     * record's (Missing) parts, exactly as Mathematica gives {0, 0}. */
    Q("Query[Select[#[\"a\"] > 1 &], Total, \"a\"][qdata]", "{0, 0}");                /* [#a] */
    Q("Query[Select[#[\"a\"] > 5 &], Total, \"a\"][qdata]", "{}");                    /* [#a] */
}

/* ---- ascending operators ---- */

static void test_query_ascending(void) {
    Q("Query[Total, \"a\"][qdata]", "6");
    Q("Query[Total][{1, 2, 3}]", "6");
    Q("Query[Total][<|\"x\" -> 1, \"y\" -> 2|>]", "3");
    Q("Query[Total, \"b\"][qaa]", "11");
    Q("Query[Max, \"c\"][qdata]", "3.5");
    Q("Query[Mean, \"c\"][qdata]", "2.5");
    Q("Query[Length][qdata]", "3");
    Q("Query[Counts, \"b\"][qdata]", "<|\"x\" -> 1, \"y\" -> 1, \"z\" -> 1|>");
    Q("Query[Last, \"b\"][qdata]", "\"z\"");
    Q("Query[First][qdata]", REC1);
    Q("Query[Total, #^2 &][{1, 2, 3, 4}]", "30");
    Q("Query[All, Total][{{1, 2}, {3, 4}}]", "{3, 7}");
    Q("Query[Total, All][{{1, 2}, {3, 4}}]", "{4, 6}");
    Q("Query[f, g][{{1, 2}, {3, 4}}]", "f[{g[{1, 2}], g[{3, 4}]}]");
    Q("Query[f][x]", "f[x]");
    Q("Query[All, f][x]", "x");
    Q("Query[All, f][{}]", "{}");
    Q("Query[All, {\"a\", \"c\"}, Total][{<|\"a\" -> 1, \"c\" -> 2|>}]", "{<|\"a\" -> 1, \"c\" -> 2|>}");
    Q("Query[Transpose][{<|\"a\" -> 1, \"b\" -> 2|>, <|\"a\" -> 3, \"b\" -> 4|>}]",
      "<|\"a\" -> {1, 3}, \"b\" -> {2, 4}|>");
    Q("Query[RightComposition[Select[#[\"a\"] > 1 &], Length]][qdata]", "2");        /* [#a] */
    Q("Query[Query[All, \"a\"]][qdata]", "{1, 2, 3}");
    Q("Query[All, Query[\"a\"]][qdata]", "{1, 2, 3}");
}

/* ---- structural operators: {f, g}, <|...|>, {part -> f} ---- */

static void test_query_structural(void) {
    Q("Query[All, {f, g}][{1, 2}]", "{{f[1], g[1]}, {f[2], g[2]}}");
    Q("Query[All, <|\"x\" -> \"a\", \"y\" -> Query[\"c\"]|>][qdata]",
      "{<|\"x\" -> 1, \"y\" -> 3.5|>, <|\"x\" -> 2, \"y\" -> 1.5|>, <|\"x\" -> 3, \"y\" -> 2.5|>}");
    Q("Query[{\"a\" -> f}][<|\"a\" -> 1, \"b\" -> 2|>]", "<|\"a\" -> f[1], \"b\" -> 2|>");
    Q("Query[All, {\"a\" -> f, \"b\" -> g}][{<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>}]",
      "{<|\"a\" -> f[1], \"b\" -> g[2], \"c\" -> 3|>}");
    Q("Query[Nothing, \"a\"][<|\"a\" -> 7|>]", "7");
}

/* ---- Normal[Query[...]] and the printed form ---- */

static void test_query_normal(void) {
    Q("Normal[Query[All]]", "Identity");
    Q("Normal[Query[\"a\"]]", "GeneralUtilities`Slice[\"a\"]");
    Q("Normal[Query[All, \"a\"]]", "GeneralUtilities`Slice[All, \"a\"]");
    Q("Normal[Query[Select[s], h]]", "RightComposition[Select[s], Map[h]]");
    Q("Normal[Query[Total, h]]", "RightComposition[Map[h], Total]");
    Q("Normal[Query[\"a\", h]]", "RightComposition[GeneralUtilities`Slice[\"a\"], h]");
    Q("Normal[Query[f, g, h]]", "RightComposition[Map[RightComposition[Map[h], g]], f]");
    Q("Normal[Query[Select[s], Total, \"a\"]]",
      "RightComposition[Select[s], Map[RightComposition[GeneralUtilities`Slice[All, \"a\"], Total]]]");
    Q("Normal[Query[GroupBy[s], h]]", "RightComposition[GroupBy[s], Map[h]]");
    Q("Normal[Query[Counts, h]]", "RightComposition[Map[h], Counts]");
    Q("Query[All, \"a\"]", "Query[All, \"a\"]");
    Q("Head[Query[All, \"a\"]]", "Query");
    Q("MemberQ[Attributes[Query], Protected]", "True");
}

/* ---- Dataset ---- */

static void test_dataset(void) {
    Q("Head[qds]", "Dataset");
    Q("Normal[qds] === qdata", "True");
    Q("Head[qds[All, \"a\"]]", "Dataset");
    Q("Normal[qds[All, \"a\"]]", "{1, 2, 3}");
    Q("qds[Total, \"a\"]", "6");
    Q("qds[Length]", "3");
    Q("Head[qds[1]]", "Dataset");
    Q("qds[1, \"a\"]", "1");
    Q("Normal[qds[Select[#[\"a\"] > 1 &]]]", "{" REC2 ", " REC3 "}");                 /* [#a] */
    Q("Normal[qds[SortBy[\"c\"], \"b\"]]", "{\"y\", \"z\", \"x\"}");
    Q("Normal[qds[GroupBy[\"b\"], Length]]", "<|\"x\" -> 1, \"y\" -> 1, \"z\" -> 1|>");
    Q("Head[qds[Counts, \"b\"]]", "Dataset");
    Q("Length[qds]", "3");
    Q("Dimensions[qds]", "{3, 3}");
    Q("Dimensions[Dataset[{1, 2, 3}]]", "{3}");
    Q("Dimensions[Dataset[{{1, 2}, {3, 4}}]]", "{2, 2}");
    Q("Dimensions[Dataset[<|\"x\" -> <|\"a\" -> 1, \"b\" -> 2|>, \"y\" -> <|\"a\" -> 1, \"b\" -> 2|>|>]]",
      "{2, 2}");
    Q("Normal[Keys[qds]]", "{{\"a\", \"b\", \"c\"}, {\"a\", \"b\", \"c\"}, {\"a\", \"b\", \"c\"}}");
    Q("Normal[Values[qds]]", "{{1, \"x\", 3.5}, {2, \"y\", 1.5}, {3, \"z\", 2.5}}");
    Q("Normal[Keys[Dataset[<|\"x\" -> 1, \"y\" -> 2|>]]]", "{\"x\", \"y\"}");
    Q("Head[Select[qds, #[\"a\"] > 1 &]]", "Dataset");                                /* [#a] */
    Q("Normal[SortBy[qds, #[\"c\"] &]]", "{" REC2 ", " REC3 ", " REC1 "}");           /* [#a] */
    Q("Normal[Map[f, Dataset[{1, 2}]]]", "{f[1], f[2]}");
    Q("Head[Map[f, qds]]", "Dataset");
    Q("Normal[First[qds]]", REC1);
    Q("Normal[Last[qds]]", REC3);
    Q("Normal[Take[qds, 2]]", "{" REC1 ", " REC2 "}");
    Q("Normal[qds[[All, \"a\"]]]", "{1, 2, 3}");
    Q("qds[[2, \"b\"]]", "\"y\"");
    Q("Head[qds[[2]]]", "Dataset");
    Q("Normal[Query[All, \"a\"][qds]]", "{1, 2, 3}");
    Q("Head[Query[All, \"a\"][qds]]", "Dataset");
    Q("Dataset[{1, 2, 3}][Total]", "6");
    Q("Dataset[{1, 2, 3}][Max]", "3");
    Q("Normal[Dataset[{1, 2, 3}][All, #^2 &]]", "{1, 4, 9}");
    Q("Dataset[<|\"x\" -> 1|>][\"x\"]", "1");
    Q("Dataset[{Missing[]}][1]", "Missing[]");
    Q("Normal[Dataset[5]]", "5");
    Q("Dataset[{}][Length]", "0");
    Q("Dataset[{1, 2}] === Dataset[{1, 2}]", "True");
    Q("Normal[Dataset[Dataset[{1}]]]", "{1}");
    Q("MemberQ[Attributes[Dataset], Protected]", "True");
    /* Normal/Length on ordinary data are unaffected by the Dataset wrappers. */
    Q("Length[{1, 2, 3}]", "3");
    Q("Map[f, {1, 2}]", "{f[1], f[2]}");
}

/* ---- the text table the REPL prints ---- */

static void check_format(const char* src, const char* expected) {
    Expr* e = parse_expression(src);
    Expr* r = evaluate(e);
    char* s = dataset_format(r);
    expr_free(e);
    if (!s || strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: format of %s\n  expected:\n%s\n  actual:\n%s\n", src, expected,
                s ? s : "(null)");
        exit(1);
    }
    free(s);
    expr_free(r);
    g_checks++;
}

static void test_dataset_format(void) {
    check_format("qds",
                 "Dataset <3 x 3>\n"
                 "a | b | c\n"
                 "--+---+----\n"
                 "1 | x | 3.5\n"
                 "2 | y | 1.5\n"
                 "3 | z | 2.5");
    check_format("Dataset[<|\"x\" -> <|\"a\" -> 1, \"b\" -> 2|>, \"y\" -> <|\"a\" -> 3|>|>]",
                 "Dataset <2>\n"
                 "  | a | b\n"
                 "--+---+--\n"
                 "x | 1 | 2\n"
                 "y | 3 | -");
    check_format("Dataset[<|\"x\" -> 1, \"y\" -> \"two\"|>]",
                 "Dataset <2>\n"
                 "x | 1\n"
                 "y | two");
    check_format("Dataset[{{1, 2}, {3, 4}}]",
                 "Dataset <2 x 2>\n"
                 "1 | 2\n"
                 "3 | 4");
    check_format("Dataset[5]", "Dataset[5]");
    check_format("Dataset[Range[22]]",
                 "Dataset <22>\n1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n11\n12\n13\n14\n15\n16\n17\n18\n19\n20\n"
                 "rows 1-20 of 22");
}

int main(void) {
    symtab_init();
    core_init();
    setup();

    TEST(test_query_parts);
    TEST(test_query_missing);
    TEST(test_query_descending);
    TEST(test_query_ascending);
    TEST(test_query_structural);
    TEST(test_query_normal);
    TEST(test_dataset);
    TEST(test_dataset_format);

    printf("All Query/Dataset tests passed (%d checks).\n", g_checks);
    return 0;
}
