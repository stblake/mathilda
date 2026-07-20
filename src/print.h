#ifndef PRINT_H
#define PRINT_H

#include "expr.h"

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
