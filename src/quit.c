/*
 * quit.c - Quit[] / Quit[n] and Exit[] / Exit[n].
 *
 * See quit.h for the design: the head stays inert (it is the sentinel the
 * control-flow constructs already stop at) and evaluation records a quit
 * request that the top-level drivers act on.  Quit[n] with a machine
 * integer sets the process exit status; any other argument emits Quit::intm
 * and does not quit.  Exit is a synonym, rewritten to Quit so the sentinel
 * sites only ever have to recognise one head.
 */

#include "quit.h"
#include "symtab.h"
#include "sym_names.h"
#include "attr.h"
#include "message.h"
#include "print.h"

#include <stdlib.h>

static int g_quit_requested = 0;
static int g_quit_status = 0;

int mth_quit_requested(void) { return g_quit_requested; }
int mth_quit_status(void)    { return g_quit_status; }
void mth_quit_cancel(void)   { g_quit_requested = 0; g_quit_status = 0; }

Expr* builtin_quit(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) {
        g_quit_requested = 1;
        g_quit_status = 0;
    } else if (argc == 1 && res->data.function.args[0]->type == EXPR_INTEGER) {
        g_quit_requested = 1;
        g_quit_status = (int)res->data.function.args[0]->data.integer;
    } else if (argc == 1) {
        char* txt = expr_to_string_outputform(res);
        mth_message("Quit", "intm",
                    "Machine-sized integer expected at position 1 in %s.",
                    txt ? txt : "Quit[...]");
        free(txt);
    }
    /* Always stay unevaluated: Quit[...] is the control-flow sentinel. */
    return NULL;
}

Expr* builtin_exit(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    /* Steal the arguments into Quit[...]; the evaluator frees the shell. */
    Expr** moved = argc ? malloc(argc * sizeof(Expr*)) : NULL;
    if (argc && !moved) return NULL;
    for (size_t i = 0; i < argc; i++) { moved[i] = args[i]; args[i] = NULL; }
    Expr* out = expr_new_function(expr_new_symbol(SYM_Quit), moved, argc);
    free(moved);
    return out;
}

void quit_init(void) {
    symtab_add_builtin("Quit", builtin_quit);
    symtab_get_def("Quit")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Quit",
        "Quit[] terminates the Mathilda session, running $Epilog first. "
        "Quit[n] terminates with process exit status n. Works at the REPL, "
        "in -file scripts and inside files read with Get.");

    symtab_add_builtin("Exit", builtin_exit);
    symtab_get_def("Exit")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Exit",
        "Exit[] terminates the Mathilda session; Exit[n] terminates with exit "
        "status n. A synonym for Quit.");
}
