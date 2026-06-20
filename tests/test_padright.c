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
        printf("PadRight test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- PadRight[list, n]: zero-pad on the right to length n ---- */
void test_padright_basic() {
    check("PadRight[{a, b, c}, 10]", "{a, b, c, 0, 0, 0, 0, 0, 0, 0}");
    check("PadRight[{1, 2, 3}, 5]", "{1, 2, 3, 0, 0}");
    /* n equal to the length leaves the list unchanged. */
    check("PadRight[{a, b, c}, 3]", "{a, b, c}");
    /* n == 0 gives the empty list. */
    check("PadRight[{a, b, c}, 0]", "{}");
}

/* ---- PadRight[list, n] truncates when n < Length[list] ---- */
void test_padright_truncate() {
    check("PadRight[{a, b, c, d, e}, 2]", "{a, b}");
    check("PadRight[{a, b, c}, 1]", "{a}");
}

/* ---- PadRight[list, n, x]: pad with a repeated element ---- */
void test_padright_element() {
    check("PadRight[{a, b, c}, 10, x]", "{a, b, c, x, x, x, x, x, x, x}");
    check("PadRight[{a, b, c}, 6, 7]", "{a, b, c, 7, 7, 7}");
}

/* ---- PadRight[list, n, {x...}]: cyclic padding ---- */
void test_padright_cyclic() {
    check("PadRight[{a, b, c}, 10, {x, y, z}]",
          "{a, b, c, x, y, z, x, y, z, x}");
    /* Two-element cycle continues the phase past the data. */
    check("PadRight[{a, b, c}, 7, {u, v}]", "{a, b, c, v, u, v, u}");
}

/* ---- A single-element padding list whose element is itself a list ---- */
void test_padright_list_element() {
    check("PadRight[{a, b, c}, 7, {{u, v}}]",
          "{a, b, c, {u, v}, {u, v}, {u, v}, {u, v}}");
}

/* ---- PadRight[list, n, padding, m]: margin m on the LEFT ---- */
void test_padright_margin() {
    check("PadRight[{a, b, c}, 10, x, 2]", "{x, x, a, b, c, x, x, x, x, x}");
    /* Centre a 1 in a length-19 list of zeros (margin 9). */
    check("PadRight[{1}, 19, 0, 9]",
          "{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}");
}

/* ---- Negative margin truncates leading elements ---- */
void test_padright_negative_margin() {
    check("PadRight[{a, b, c}, 7, {u, v}, -2]", "{c, v, u, v, u, v, u}");
}

/* ---- Negative length pads on the LEFT instead ---- */
void test_padright_negative_length() {
    check("PadRight[{a, b, c}, -10]", "{0, 0, 0, 0, 0, 0, 0, a, b, c}");
    check("PadRight[{a, b, c}, -5, x]", "{x, x, a, b, c}");
}

/* ---- The head need not be List; it is preserved ---- */
void test_padright_head() {
    check("PadRight[f[a, b, c], 8, x]", "f[a, b, c, x, x, x, x, x]");
}

/* ---- Empty list with cyclic padding repeats the sequence ---- */
void test_padright_empty() {
    check("PadRight[{}, 10, {x, y, z}]",
          "{x, y, z, x, y, z, x, y, z, x}");
    check("PadRight[{}, 4]", "{0, 0, 0, 0}");
}

/* ---- Multidimensional: PadRight[list, {n1, n2}] ---- */
void test_padright_multidim() {
    check("PadRight[{{a, b}, {c}}, {3, 5}]",
          "{{a, b, 0, 0, 0}, {c, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}");
    /* Atom where a row is expected becomes a 1-element row. */
    check("PadRight[{a, {b, c}}, {2, 3}]", "{{a, 0, 0}, {b, c, 0}}");
}

/* ---- Multidimensional with a repeated padding block ---- */
void test_padright_block() {
    check("PadRight[{{aa, bb}, {cc}}, {4, 4}, {{x, y}, {z}}]",
          "{{aa, bb, x, y}, {cc, z, z, z}, {x, y, x, y}, {z, z, z, z}}");
}

/* ---- Multidimensional with per-level margins ---- */
void test_padright_block_margins() {
    check("PadRight[{{aa, bb}, {cc}}, {5, 5}, {{x, y}, {z}}, {1, 2}]",
          "{{z, z, z, z, z}, {x, y, aa, bb, x}, {z, z, cc, z, z}, "
          "{x, y, x, y, x}, {z, z, z, z, z}}");
}

/* ---- PadRight[list]: pad a ragged array with zeros to full ---- */
void test_padright_full() {
    check("PadRight[{{a, b, c}, {d, e}, {f}}]",
          "{{a, b, c}, {d, e, 0}, {f, 0, 0}}");
    /* A flat list of atoms is already full. */
    check("PadRight[{a, b, c}]", "{a, b, c}");
    check("PadRight[{{1}, {2, 3}, {4, 5, 6}}]",
          "{{1, 0, 0}, {2, 3, 0}, {4, 5, 6}}");
}

/* ---- PadRight[list, Automatic, x]: pad to full with x ---- */
void test_padright_automatic() {
    check("PadRight[{{a, b, c}, {d, e}, {f}}, Automatic, q]",
          "{{a, b, c}, {d, e, q}, {f, q, q}}");
}

/* ---- Padding values are evaluated/normalised through the engine ---- */
void test_padright_evaluation() {
    check("PadRight[{1, 2}, 4, 1 + 1]", "{1, 2, 2, 2}");
    check("PadRight[{x}, 3, 2 + 3]", "{x, 5, 5}");
}

/* ---- Wrong argument counts leave the call unevaluated (after argb) ---- */
void test_padright_arg_errors() {
    check("PadRight[]", "PadRight[]");
    check("PadRight[1, 2, 3, 4, 5, 6, 6]", "PadRight[1, 2, 3, 4, 5, 6, 6]");
}

/* ---- A non-array first argument is left unevaluated ---- */
void test_padright_non_array() {
    check("PadRight[x, 5]", "PadRight[x, 5]");
    check("PadRight[5, 3, 0]", "PadRight[5, 3, 0]");
}

/* ---- Attributes: Protected ---- */
void test_padright_attributes() {
    check("MemberQ[Attributes[PadRight], Protected]", "True");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_padright_basic);
    TEST(test_padright_truncate);
    TEST(test_padright_element);
    TEST(test_padright_cyclic);
    TEST(test_padright_list_element);
    TEST(test_padright_margin);
    TEST(test_padright_negative_margin);
    TEST(test_padright_negative_length);
    TEST(test_padright_head);
    TEST(test_padright_empty);
    TEST(test_padright_multidim);
    TEST(test_padright_block);
    TEST(test_padright_block_margins);
    TEST(test_padright_full);
    TEST(test_padright_automatic);
    TEST(test_padright_evaluation);
    TEST(test_padright_arg_errors);
    TEST(test_padright_non_array);
    TEST(test_padright_attributes);

    printf("All PadRight tests passed!\n");
    return 0;
}
