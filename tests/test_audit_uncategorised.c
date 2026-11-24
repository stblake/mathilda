/* Builtin-audit regressions: uncategorised heads (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"

/* SlotSequence#1 (high): "## with no remaining arguments is left
 * unsubstituted instead of becoming an empty sequence". */
static void test_slotsequence_empty(void) {
    audit_run("ClearAll[f, a, b, c]");
    audit_chk("SlotSequence#1", "{##} &[]", "{}");
    audit_chk("SlotSequence#1", "{a, ##} &[]", "{a}");
    audit_chk("SlotSequence#1", "{##3} &[a, b]", "{}");
    audit_chk("SlotSequence#1", "f[## &[]]", "f[]");
    /* sanity */
    audit_chk("SlotSequence#1", "{##2} &[a, b, c]", "{b, c}");
    audit_chk("SlotSequence#1", "f[##] &[a, b]", "f[a, b]");
}

/* Decrement#1 (high): "Decrement evaluates a Part l-value's index twice:
 * returns one element and decrements another".  Shared by every
 * read-modify-write operator (increment_core, AppendTo/PrependTo). */
static void test_decrement_part_once(void) {
    audit_run("ClearAll[di, dl]; di = 1; dl = {10, 20, 30}");
    audit_chk("Decrement#1", "{dl[[di++]]--, di, dl}", "{10, 2, {9, 20, 30}}");
    audit_run("di = 1; dl = {10, 20, 30}");
    audit_chk("Decrement#1", "{dl[[di++]]++, di, dl}", "{10, 2, {11, 20, 30}}");
    audit_run("di = 1; dl = {10, 20, 30}");
    audit_chk("Decrement#1", "{++dl[[di++]], di, dl}", "{11, 2, {11, 20, 30}}");
    audit_run("di = 1; dl = {10, 20, 30}");
    audit_chk("Decrement#1", "{--dl[[di++]], di, dl}", "{9, 2, {9, 20, 30}}");
    audit_run("di = 1; dl = {{1}, {2}}");
    audit_chk("Decrement#1", "AppendTo[dl[[di++]], 9]; {dl, di}", "{{{1, 9}, {2}}, 2}");
    audit_run("di = 1; dl = {{1}, {2}}");
    audit_chk("Decrement#1", "PrependTo[dl[[di++]], 9]; {dl, di}", "{{{9, 1}, {2}}, 2}");
    audit_run("di = 1; dl = {{10, 20}, {30, 40}}");
    audit_chk("Decrement#1", "dl[[di++, di++]] -= 1; {dl, di}", "{{{10, 19}, {30, 40}}, 3}");
    /* sanity */
    audit_run("di = 5");
    audit_chk("Decrement#1", "{di--, di}", "{5, 4}");
    audit_run("dl = {1, 2, 3}");
    audit_chk("Decrement#1", "dl[[2]] += 10; dl", "{1, 12, 3}");
    audit_chk("Decrement#1", "Catch[dl[[Throw[7]]]++]", "7");
}

/* Divide#2 (high): "Divide[0., 0.] (and 0/0., 0./0) gives ComplexInfinity
 * instead of Indeterminate". */
static void test_divide_zero_zero_real(void) {
    audit_chk("Divide#2", "Quiet[Divide[0., 0.]]", "Indeterminate");
    audit_chk("Divide#2", "Quiet[Divide[0, 0.]]", "Indeterminate");
    audit_chk("Divide#2", "Quiet[Divide[0., 0]]", "Indeterminate");
    audit_chk("Divide#2", "Check[Divide[0., 0.], fired, Infinity::indet]", "fired");
    /* sanity: a nonzero dividend is still ComplexInfinity */
    audit_chk("Divide#2", "Quiet[Divide[1., 0.]]", "ComplexInfinity");
    audit_chk("Divide#2", "Quiet[Divide[0, 0]]", "Indeterminate");
    audit_chk("Divide#2", "Divide[1., 4]", "0.25");
}

/* TeXForm#1 (high): "2 3^x typesets as 23^x".  Math mode ignores spaces,
 * so a factor that starts with a digit gets an explicit "\ ". */
static void test_texform_digit_factor(void) {
    /* "\\ " in the language is one backslash and a space */
    audit_true("TeXForm#1", "ToString[TeXForm[2*3^x]] === \"2\\\\ 3^{x}\"");
    audit_true("TeXForm#1", "ToString[TeXForm[2*3^x/y]] === \"\\\\frac{2\\\\ 3^{x}}{y}\"");
    /* sanity: a symbol factor needs no explicit space */
    audit_true("TeXForm#1", "ToString[TeXForm[2 x]] === \"2 x\"");
}

/* ListConvolve#1 (medium): the level argument lev below the list's rank. */
static void test_s14_listconvolve_level(void) {
    audit_run("ClearAll[a, b, c, d, e, f, x, y]");
    audit_true("ListConvolve#1",
               "Expand[ListConvolve[{x, y}, {{a, b}, {c, d}, {e, f}}, {-1, 1}, 0, Times, Plus, 1]] === "
               "Expand[{{c x + a y, d x + b y}, {e x + c y, f x + d y}}]");
    audit_true("ListConvolve#1", "ListConvolve[{1, 2}, {3, 4, 5}, {-1, 1}, 0, Times, Plus, 1] === "
                                 "ListConvolve[{1, 2}, {3, 4, 5}]");
}

int main(void) {
    audit_init();
    TEST(test_slotsequence_empty);
    TEST(test_decrement_part_once);
    TEST(test_divide_zero_zero_real);
    TEST(test_texform_digit_factor);
    TEST(test_s14_listconvolve_level);
    return audit_finish("audit_uncategorised");
}
