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

/*
 * The same, as a FILE is read (Get, Read, ReadList): a top-level line break
 * after a complete expression ends the statement whatever the next line
 * starts with, so the two lines "5" and "-3" are two statements, not 5 - 3.
 * An incomplete expression (an open bracket, a trailing operator) still
 * continues onto the next line.
 */
Expr* parse_next_line_statement(const char** input_ptr);

/*
 * Parses `input` as ToExpression reads a string: a sequence of statements,
 * returned as List[s1, s2, ...] (NULL on a syntax error, with *err_pos -- if
 * err_pos is non-NULL -- set near the offending text).
 *
 * With line_statements set (InputForm), a top-level line break ends a
 * complete statement whatever follows it (`a\n+b` is {a, +b}, `a\n// f` is a
 * syntax error), a `;` ending a line takes an implicit Null (`a;\nb` is
 * {a;, b}), and each empty, blank or comment-only line between statements is
 * a Null statement (`a\n\nb` is {a, Null, b}); LF and CR each terminate a
 * line. A `;` chain on one line stays a single CompoundExpression. Without it
 * (StandardForm), line breaks are whitespace and the whole input must be one
 * expression (`a\nb` is {a b}). Input with no statement at all is {Null}.
 */
Expr* parse_statement_sequence(const char* input, int line_statements,
                               const char** err_pos);

/* The same reader, one statement at a time: each statement is handed to `fn`
 * (which takes ownership) as soon as it is parsed and before the next is read;
 * `fn` returns 0 to stop.  Returns 1 when the input was read to the end, 0 on a
 * syntax error (with *err_pos set, the earlier statements already handed over),
 * 2 when `fn` stopped it. */
typedef int (*ParseStatementFn)(Expr* stmt, void* ctx);
int parse_statement_each(const char* input, int line_statements, const char** err_pos,
                         ParseStatementFn fn, void* ctx);

/* Suppress (quiet = 1) or restore the parser's own syntax-error printing;
 * returns the previous setting.  A failed parse still returns NULL. */
int parse_set_quiet(int quiet);

/* Make every symbol in a parsed user statement an existing symbol (outside a
 * kernel module load); see parse.c. */
void parse_materialize_user_symbols(const Expr* e);

#endif // PARSE_H
