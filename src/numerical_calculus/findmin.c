/* findmin.c — registration hub (findmin_init).
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void findmin_init(void) {
    symtab_add_builtin("FindMinimum", builtin_findminimum);
    symtab_get_def("FindMinimum")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_add_builtin("FindMaximum", builtin_findmaximum);
    symtab_get_def("FindMaximum")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    /* NMinimize/NMaximize are Protected but NOT HoldAll (matching Mathematica,
     * whose Attributes[NMinimize] is {Protected}). Their variables are ordinary
     * unbound symbols that evaluate to themselves, and the objective is
     * re-evaluated per trial point under a Block-style binding of those symbols,
     * so holding the arguments is unnecessary — and holding them is what made a
     * Method sub-option value such as "RandomSeed" -> s (s a Do/Table iterator or
     * any expression) arrive unevaluated and get silently dropped. FindMinimum
     * stays HoldAll: its {x, x0} specs pair a variable with an initial value that
     * must not evaluate. */
    symtab_add_builtin("NMinimize", builtin_nminimize);
    symtab_get_def("NMinimize")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("NMaximize", builtin_nmaximize);
    symtab_get_def("NMaximize")->attributes |= ATTR_PROTECTED;
}
