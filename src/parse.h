#ifndef PARSE_H
#define PARSE_H

#include "expr.h"

/*
 * Parses Mathematica FullForm expressions into Expr trees
 * Supported syntax:
 * - Symbols:        x, `name`, $var
 * - Numbers:        123, -45, 3.14, 1e5
 * - Strings:        "text"
 * - Functions:      f[x,y], Derivative[1][f][x]
 * - Lists:          {1,2,3}
 * 
 * Returns: New Expr tree on success, NULL on failure
 * Memory: Caller must free result with expr_free()
 */
Expr* parse_expression(const char* input);

/*
 * Parses the next top-level STATEMENT from the input string pointer,
 * advancing the pointer past the parsed statement and its trailing ';'
 * separator (if any). Unlike parse_expression, a ';'-chain is NOT folded into
 * a single CompoundExpression: each ';'-terminated statement is returned
 * separately, so a caller that evaluates between calls (the file loader) lets
 * context-changing prologues such as BeginPackage[]/Begin[] take effect on the
 * symbols parsed in later statements. Leading/empty ';' separators are skipped.
 * Returns NULL when no more statements can be parsed.
 */
Expr* parse_next_expression(const char** input_ptr);

#endif // PARSE_H
