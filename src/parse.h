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
 * Source spans for structural (bottom-up) selection in the notebook.
 *
 * mth_parse_spans parses `input` WITHOUT evaluating and returns, for every
 * subexpression the parser builds, its half-open [start, end) byte range into
 * `input` — including each precedence stage (so `a + b*c` yields `a`, `b`, `c`,
 * `b*c`, `a + b*c`). This is what the editor needs for Mathematica-style
 * balanced selection and cannot get from the Expr tree, which carries no
 * positions. It is tolerant: incomplete/!syntactic input yields the spans of the
 * parts that did parse. The caller frees the result with mth_span_sink_free.
 */
typedef struct { int start; int end; } MthSpan;
typedef struct { MthSpan* v; size_t n, cap; } MthSpanSink;
MthSpanSink* mth_parse_spans(const char* input);
void mth_span_sink_free(MthSpanSink* sink);

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
