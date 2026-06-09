/* Tests for StieltjesGamma[n], the inert Stieltjes constants.
 *
 * Covers the StieltjesGamma[0] -> EulerGamma reduction, inert (symbolic)
 * behaviour for higher indices, Listable threading, and attributes. */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>

void test_stieltjes_zero_is_eulergamma() {
    assert_eval_eq("StieltjesGamma[0]", "EulerGamma", 0);
}

void test_stieltjes_inert() {
    /* Higher indices have no elementary closed form -> stay symbolic. */
    assert_eval_eq("StieltjesGamma[1]", "StieltjesGamma[1]", 0);
    assert_eval_eq("StieltjesGamma[2]", "StieltjesGamma[2]", 0);
    assert_eval_eq("StieltjesGamma[n]", "StieltjesGamma[n]", 0);
}

void test_stieltjes_listable() {
    assert_eval_eq("StieltjesGamma[{0, 1, 2}]",
                   "{EulerGamma, StieltjesGamma[1], StieltjesGamma[2]}", 0);
}

void test_stieltjes_attributes() {
    SymbolDef* d = symtab_get_def("StieltjesGamma");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "StieltjesGamma must be Listable");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "StieltjesGamma must be Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_stieltjes_zero_is_eulergamma);
    TEST(test_stieltjes_inert);
    TEST(test_stieltjes_listable);
    TEST(test_stieltjes_attributes);

    printf("All StieltjesGamma tests passed.\n");
    return 0;
}
