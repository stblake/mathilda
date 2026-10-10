#ifndef PRINT_H
#define PRINT_H

#include <stdio.h>
#include "expr.h"

/* Print-subsystem output sink.  Printing code (print.c and numberform.c) writes
 * to mth_out() rather than to stdout directly, so the string-capture helpers
 * (expr_to_string & friends) can divert all printed output into an in-memory
 * stream.  They used to do that by reassigning `stdout`, but on musl libc
 * `stdout` is a `FILE *const` and cannot be assigned — it builds on glibc and
 * macOS yet breaks the aarch64-musl package build.  Routing through a sink the
 * subsystem owns is portable.  mth_out() returns the real stdout when no
 * capture is active.  mth_out_push() installs `f` and returns the previous
 * sink; mth_out_pop() restores it. */
FILE* mth_out(void);
FILE* mth_out_push(FILE* f);
void  mth_out_pop(FILE* prev);

void expr_print(Expr* e);
void expr_print_fullform(Expr* e);
char* expr_to_string(Expr* e);
char* expr_to_string_fullform(Expr* e);

/* Builds the pretty Plus/Times/Power/O equivalent of a 6-argument
 * SeriesData[...] node, or NULL if it is not in a renderable shape. Shared by
 * the plain and LaTeX printers. Caller owns and must free the result. */
Expr* series_data_to_display_expr(Expr* e);

Expr* builtin_print(Expr* res);
Expr* builtin_fullform(Expr* res);
Expr* builtin_inputform(Expr* res);
Expr* builtin_texform(Expr* res);


#endif // PRINT_H
