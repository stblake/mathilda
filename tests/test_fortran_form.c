/*
 * test_fortran_form.c -- FortranForm / CForm text (src/print_fortran.c), the
 * FormatType/PageWidth options of ToString and OpenWrite, and the
 * MantissaExponent machine path FeynCalc's Write2 relies on.
 *
 * Every expected string was produced by the reference kernel on the same
 * input (ToString[FortranForm[e]], ToString[CForm[e]], ToString[e,
 * FormatType -> FortranForm, PageWidth -> n]).
 */
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Evaluate input; it must give a String equal to expected. */
static void check_str(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    if (!parsed) { fprintf(stderr, "FAIL: does not parse: %s\n", input); abort(); }
    Expr* v = evaluate(parsed);
    expr_free(parsed);
    if (!v || v->type != EXPR_STRING || strcmp(v->data.string, expected) != 0) {
        char* got = v ? expr_to_string(v) : NULL;
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected,
                got ? got : "(null)");
        free(got);
        abort();
    }
    expr_free(v);
}

static void ff(const char* e, const char* fortran, const char* c) {
    char buf[1024];
    snprintf(buf, sizeof buf, "ToString[FortranForm[%s]]", e);
    check_str(buf, fortran);
    snprintf(buf, sizeof buf, "ToString[CForm[%s]]", e);
    check_str(buf, c);
}

static void test_atoms(void) {
    ff("100.", "100.", "100.");
    ff("1.*^6", "1.e6", "1.e6");
    ff("1234567.", "1.234567e6", "1.234567e6");
    ff("0.00001", "0.00001", "0.00001");
    ff("1.5*^-7", "1.5e-7", "1.5e-7");
    ff("2/3", "0.6666666666666666", "0.6666666666666666");
    ff("12345678901234567890", "12345678901234567890", "12345678901234567890");
    ff("2 + 3 I", "(2,3)", "Complex(2,3)");
    ff("Complex[1.5, -2.]", "(1.5,-2.)", "Complex(1.5,-2.)");
    ff("\"str\"", "\"str\"", "\"str\"");
    ff("Infinity", "DirectedInfinity(1)", "DirectedInfinity(1)");
    ff("-Infinity", "DirectedInfinity(-1)", "DirectedInfinity(-1)");
}

static void test_products_and_powers(void) {
    ff("a/100 + 5 x", "a/100. + 5*x", "a/100. + 5*x");
    ff("-2/3 x", "(-2*x)/3.", "(-2*x)/3.");
    ff("-x/(2 y)", "-0.5*x/y", "-0.5*x/y");
    ff("1/(2 x)", "1/(2.*x)", "1/(2.*x)");
    ff("a b c^-1 d^-2", "(a*b)/(c*d**2)", "(a*b)/(c*Power(d,2))");
    ff("-x y", "-(x*y)", "-(x*y)");
    ff("-a/b", "-(a/b)", "-(a/b)");
    ff("-x^-1", "-(1/x)", "-(1/x)");
    ff("-x^-2", "-x**(-2)", "-Power(x,-2)");
    ff("x^-2", "x**(-2)", "Power(x,-2)");
    ff("y x^(-3/2)", "y/x**1.5", "y/Power(x,1.5)");
    ff("1/x^(2/3)", "x**(-0.6666666666666666)", "Power(x,-0.6666666666666666)");
    ff("1/Sqrt[x]", "1/Sqrt(x)", "1/Sqrt(x)");
    ff("x^y^z", "x**y**z", "Power(x,Power(y,z))");
    ff("(x^y)^z", "(x**y)**z", "Power(Power(x,y),z)");
    ff("(-2)^x", "(-2)**x", "Power(-2,x)");
    ff("(a + b)^n", "(a + b)**n", "Power(a + b,n)");
    ff("E^x", "E**x", "Power(E,x)");
    ff("x y^(1/2) z^(-1/2)", "(x*Sqrt(y))/Sqrt(z)", "(x*Sqrt(y))/Sqrt(z)");
}

static void test_sums_and_relations(void) {
    ff("a - b c - d", "a - b*c - d", "a - b*c - d");
    ff("a - I b", "a - (0,1)*b", "a - Complex(0,1)*b");
    ff("a - b/2", "a - b/2.", "a - b/2.");
    ff("Cos[y] Sin[x]", "Cos(y)*Sin(x)", "Cos(y)*Sin(x)");
    ff("{1, 2, x}", "List(1,2,x)", "List(1,2,x)");
    ff("x == y", "x.eq.y", "x == y");
    ff("a != b", "a.ne.b", "a != b");
    ff("a >= b", "a.ge.b", "a >= b");
    ff("x && y", "x.and.y", "x && y");
    ff("!a", ".not.a", "!a");
    ff("Mod[a, b]", "Mod(a,b)", "a % b");
    ff("h[x][y]", "h(x)(y)", "h(x)(y)");
    ff("HoldForm[a + b]", "a + b", "a + b");
}

static void test_format_values(void) {
    /* Format[..., FortranForm] applies; a plain Format rule does not. */
    check_str("fh /: Format[fh[r_Real], FortranForm] := SequenceForm[r, D, 0]; "
              "ToString[FortranForm[x^(-fh[0.6]) + a fh[2.]]]",
              "x**(-0.6D0) + a*2.D0");
    check_str("Format[gg[x_]] := zz[x]; ToString[FortranForm[gg[1]]]", "gg(1)");
}

static void test_line_breaking(void) {
    /* fixed-form continuation: 8-space first line, "     -" + depth + 1 */
    check_str("ToString[Sum[ToExpression[\"v\" <> ToString[k]], {k, 7}], "
              "FormatType -> FortranForm, PageWidth -> 30]",
              "        v1 + v2 + v3 + v4 + \n     -  v5 + v6 + v7");
    check_str("ToString[Sum[ToExpression[\"v\" <> ToString[k]], {k, 6}], "
              "FormatType -> FortranForm, PageWidth -> 30]",
              "v1 + v2 + v3 + v4 + v5 + v6");
    /* an over-long atom is split with a backslash */
    check_str("ToString[aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa + b, "
              "FormatType -> FortranForm, PageWidth -> 30]",
              "        aaaaaaaaaaaaaaaaaaaaa\\\n     -   aaaaaaaaaaaaaaa + b");
}

static void test_stream_formattype(void) {
    char path[256], prog[1024];
    snprintf(path, sizeof path, "/tmp/Mathilda_fortran_%d.f", (int)getpid());
    unlink(path);
    snprintf(prog, sizeof prog,
             "p = \"%s\"; s = OpenWrite[p, FormatType -> FortranForm, PageWidth -> 30]; "
             "WriteString[s, \"        \", FortranForm[t], \" = \"]; "
             "Write[s, Sum[ToExpression[\"v\" <> ToString[k]], {k, 7}]]; Close[s]; "
             "ReadString[p]", path);
    check_str(prog, "        t =         v1 + v2 + v3 + v4 + \n     -  v5 + v6 + v7\n");
    unlink(path);
}

static void test_mantissa_exponent(void) {
    /* the exponent is decided on the double's exact value */
    assert_eval_eq("MantissaExponent[0.01] === {0.1, -1}", "True", 0);
    assert_eval_eq("MantissaExponent[-0.01] === {-0.1, -1}", "True", 0);
    assert_eval_eq("MantissaExponent[0.001] === {0.1, -2}", "True", 0);
    assert_eval_eq("MantissaExponent[100.] === {0.1, 3}", "True", 0);
    assert_eval_eq("MantissaExponent[999.9999999999999] === {0.9999999999999999, 3}", "True", 0);
    assert_eval_eq("MantissaExponent[1.*^-300] === {0.1, -299}", "True", 0);
    assert_eval_eq("MantissaExponent[1024., 2] === {0.5, 11}", "True", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_atoms);
    TEST(test_products_and_powers);
    TEST(test_sums_and_relations);
    TEST(test_format_values);
    TEST(test_line_breaking);
    TEST(test_stream_formattype);
    TEST(test_mantissa_exponent);
    printf("All FortranForm tests passed.\n");
    return 0;
}
