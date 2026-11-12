/*
 * frontend.h -- builtins for the interactive front end.
 *
 *   SlashCommand["name", f, "summary"]   adds the REPL command /name, which
 *                                        evaluates f["argument text"]
 *   Monitor[expr, mon]                   evaluates expr, showing the current
 *                                        value of mon in the REPL's status
 *                                        line meanwhile; returns expr's value
 *
 * Both are harmless outside the interactive REPL: SlashCommand only records
 * the command (nothing runs it), Monitor just evaluates expr.
 */
#ifndef MATHILDA_FRONTEND_H
#define MATHILDA_FRONTEND_H

#include "expr.h"

Expr* builtin_slashcommand(Expr* res);
Expr* builtin_monitor(Expr* res);

/* The mod being loaded (its file name), recorded as the owner of the
 * commands it defines; NULL outside a mod load. The string is copied. */
void  frontend_set_mod_owner(const char* name);

/* mth_monitor_tick (declared in term/term_progress.h): the evaluator's
 * periodic check or a loop step found g_mth_monitor_due set; it evaluates
 * the innermost Monitor's `mon` and publishes its text to the status line. */

void  frontend_init(void);

#endif /* MATHILDA_FRONTEND_H */
