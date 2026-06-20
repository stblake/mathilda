#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

/* Shared driver: parse, evaluate, compare the printed form. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("PadLeft test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- PadLeft[list, n]: zero-pad on the left to length n ---- */
void test_padleft_basic() {
    check("PadLeft[{a, b, c}, 10]", "{0, 0, 0, 0, 0, 0, 0, a, b, c}");
    check("PadLeft[{1, 2, 3}, 5]", "{0, 0, 1, 2, 3}");
    check("PadLeft[{a, b, c}, 3]", "{a, b, c}");
    check("PadLeft[{a, b, c}, 0]", "{}");
}

/* ---- PadLeft truncates from the front when n < Length[list] ---- */
void test_padleft_truncate() {
    check("PadLeft[{a, b, c, d, e}, 2]", "{d, e}");
    check("PadLeft[{a, b, c}, 1]", "{c}");
}

/* ---- PadLeft[list, n, x]: pad with a repeated element ---- */
void test_padleft_element() {
    check("PadLeft[{a, b, c}, 10, x]", "{x, x, x, x, x, x, x, a, b, c}");
    check("PadLeft[{a, b, c}, 6, 7]", "{7, 7, 7, a, b, c}");
}

/* ---- PadLeft[list, n, {x...}]: cyclic padding (doc examples) ---- */
void test_padleft_cyclic() {
    check("PadLeft[{a, b, c}, 10, {x, y, z}]",
          "{z, x, y, z, x, y, z, a, b, c}");
    check("PadLeft[{a, b, c}, 7, {u, v}]", "{v, u, v, u, a, b, c}");
}

/* ---- A single-element padding list whose element is itself a list ---- */
void test_padleft_list_element() {
    check("PadLeft[{a, b, c}, 7, {{u, v}}]",
          "{{u, v}, {u, v}, {u, v}, {u, v}, a, b, c}");
}

/* ---- PadLeft[list, n, padding, m]: margin m on the RIGHT ---- */
void test_padleft_margin() {
    check("PadLeft[{a, b, c}, 10, x, 2]", "{x, x, x, x, x, a, b, c, x, x}");
    /* Centre a 1 in a length-19 list of zeros (margin 9). */
    check("PadLeft[{1}, 19, 0, 9]",
          "{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}");
}

/* ---- Negative margin truncates trailing elements ---- */
void test_padleft_negative_margin() {
    check("PadLeft[{a, b, c}, 7, {u, v}, -2]", "{v, u, v, u, v, u, a}");
}

/* ---- Negative length pads on the RIGHT instead ---- */
void test_padleft_negative_length() {
    check("PadLeft[{a, b, c}, -10]", "{a, b, c, 0, 0, 0, 0, 0, 0, 0}");
    check("PadLeft[{a, b, c}, -5, x]", "{a, b, c, x, x}");
}

/* ---- The head need not be List; it is preserved ---- */
void test_padleft_head() {
    check("PadLeft[f[a, b, c], 8, x]", "f[x, x, x, x, x, a, b, c]");
}

/* ---- Empty list with cyclic padding repeats the sequence ---- */
void test_padleft_empty() {
    check("PadLeft[{}, 10, {x, y, z}]",
          "{x, y, z, x, y, z, x, y, z, x}");
    check("PadLeft[{}, 4]", "{0, 0, 0, 0}");
}

/* ---- Multidimensional: PadLeft[list, {n1, n2}] (doc example) ---- */
void test_padleft_multidim() {
    check("PadLeft[{{a, b}, {c}}, {3, 5}]",
          "{{0, 0, 0, 0, 0}, {0, 0, 0, a, b}, {0, 0, 0, 0, c}}");
    /* Atom where a row is expected becomes a 1-element row. */
    check("PadLeft[{a, {b, c}}, {2, 3}]", "{{0, 0, a}, {0, b, c}}");
}

/* ---- Multidimensional with a repeated padding block (doc example) ---- */
void test_padleft_block() {
    check("PadLeft[{{aa, bb}, {cc}}, {4, 4}, {{x, y}, {z}}]",
          "{{x, y, x, y}, {z, z, z, z}, {x, y, aa, bb}, {z, z, z, cc}}");
}

/* ---- Multidimensional with per-level margins (doc example) ---- */
void test_padleft_block_margins() {
    check("PadLeft[{{aa, bb}, {cc}}, {5, 5}, {{x, y}, {z}}, {1, 2}]",
          "{{y, x, y, x, y}, {z, z, z, z, z}, {y, aa, bb, x, y}, "
          "{z, z, cc, z, z}, {y, x, y, x, y}}");
}

/* ---- PadLeft[list]: pad a ragged array with zeros to full (doc example) ---- */
void test_padleft_full() {
    check("PadLeft[{{a, b, c}, {d, e}, {f}}]",
          "{{a, b, c}, {0, d, e}, {0, 0, f}}");
    check("PadLeft[{a, b, c}]", "{a, b, c}");
    check("PadLeft[{{1}, {2, 3}, {4, 5, 6}}]",
          "{{0, 0, 1}, {0, 2, 3}, {4, 5, 6}}");
    /* Pad digit lists to a common length (doc example). */
    check("PadLeft[Table[IntegerDigits[i^2, 2], {i, 5}]]",
          "{{0, 0, 0, 0, 1}, {0, 0, 1, 0, 0}, {0, 1, 0, 0, 1}, "
          "{1, 0, 0, 0, 0}, {1, 1, 0, 0, 1}}");
}

/* ---- PadLeft[list, Automatic, x]: pad to full with x ---- */
void test_padleft_automatic() {
    check("PadLeft[{{a, b, c}, {d, e}, {f}}, Automatic, q]",
          "{{a, b, c}, {q, d, e}, {q, q, f}}");
}

/* ---- Padding values are evaluated/normalised through the engine ---- */
void test_padleft_evaluation() {
    check("PadLeft[{1, 2}, 4, 1 + 1]", "{2, 2, 1, 2}");
    check("PadLeft[{x}, 3, 2 + 3]", "{5, 5, x}");
}

/* ---- Wrong argument counts leave the call unevaluated (after argb) ---- */
void test_padleft_arg_errors() {
    check("PadLeft[]", "PadLeft[]");
    check("PadLeft[1, 2, 3, 4, 5, 6, 6]", "PadLeft[1, 2, 3, 4, 5, 6, 6]");
}

/* ---- A non-array first argument is left unevaluated ---- */
void test_padleft_non_array() {
    check("PadLeft[x, 5]", "PadLeft[x, 5]");
    check("PadLeft[5, 3, 0]", "PadLeft[5, 3, 0]");
}

/* ---- Attributes: Protected ---- */
void test_padleft_attributes() {
    check("MemberQ[Attributes[PadLeft], Protected]", "True");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_padleft_basic);
    TEST(test_padleft_truncate);
    TEST(test_padleft_element);
    TEST(test_padleft_cyclic);
    TEST(test_padleft_list_element);
    TEST(test_padleft_margin);
    TEST(test_padleft_negative_margin);
    TEST(test_padleft_negative_length);
    TEST(test_padleft_head);
    TEST(test_padleft_empty);
    TEST(test_padleft_multidim);
    TEST(test_padleft_block);
    TEST(test_padleft_block_margins);
    TEST(test_padleft_full);
    TEST(test_padleft_automatic);
    TEST(test_padleft_evaluation);
    TEST(test_padleft_arg_errors);
    TEST(test_padleft_non_array);
    TEST(test_padleft_attributes);

    printf("All PadLeft tests passed!\n");
    return 0;
}
