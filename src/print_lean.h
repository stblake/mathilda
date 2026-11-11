/*
 * print_lean.h -- LeanForm: Lean 4 text of an expression.
 *
 * Two renderers share one recursive printer (print_lean.c):
 *
 *   - lean_form_string: the DISPLAY form behind LeanForm[e] and
 *     ToString[e, LeanForm].  Permissive: every expression prints (an
 *     unknown head f[x, y] as the application `f x y`, a machine real as a
 *     decimal literal), and Format[f[...], LeanForm] := rhs display rules
 *     apply.  Never NULL.
 *
 *   - lean_form_string_ex: the STRICT translator used by the certificate
 *     generators (Lean`Internal`ToLean).  It refuses -- returns NULL and names the
 *     reason and the offending subexpression -- every node whose Lean text
 *     would not mean the same thing as the Mathilda expression in the chosen
 *     number type: inexact numbers, division over Int (Int `/` is floor
 *     division), transcendental heads and constants when restricted to core
 *     Lean, infinities, and any head that is not in the generator's
 *     pass-through set.  Format rules do NOT apply here, so a user's display
 *     rule can never change the statement handed to Lean.
 *
 * Precedences are Lean 4's: <-> 20, -> 25 (right), \/ 30, /\ 35, not 40,
 * relations 50, + - 65, * / 70, ^ 75 (right), application 1024.  Unary
 * negation always prints -(...) unless its operand is an atom, binders are
 * parenthesised unless they are the whole expression, and a negative literal
 * in operand position is parenthesised: (-3).  There is no line breaking.
 */
#ifndef MATHILDA_PRINT_LEAN_H
#define MATHILDA_PRINT_LEAN_H

#include <stddef.h>
#include "expr.h"

typedef enum { LEAN_NUM_NONE = 0, LEAN_NUM_INT, LEAN_NUM_RAT, LEAN_NUM_REAL } LeanNumType;

typedef struct {
    LeanNumType num;        /* ascription type for rational literals / closed relations and
                             * the type of ForAll / Exists binders (NONE: untyped) */
    int core_only;          /* 1: reject transcendental / Real-only heads (Tier 0) */
    int allow_division;     /* 0 forbids `/`; forced to 0 under LEAN_NUM_INT */
    const char* const* pass_heads;  /* Lean names allowed as applied heads (strict form) */
    size_t n_pass_heads;
    int closed;             /* 1: the text is a whole theorem statement -- every symbol must be
                             * bound by an enclosing ForAll / Exists / Function or be a
                             * pass-through name, so no free identifier can resolve to a Lean
                             * global (default, id, ...) */
} LeanFormOpts;

/* Display form (LeanForm[e], ToString[e, LeanForm]): permissive, never NULL, caller frees. */
char* lean_form_string(const Expr* e);

/* Strict translator: NULL on an untranslatable node; *err (static string) names the
 * reason and *bad (borrowed from e) the offending subexpression.  err / bad may be
 * NULL.  Caller frees the result. */
char* lean_form_string_ex(const Expr* e, const LeanFormOpts* o, const char** err,
                          const Expr** bad);

/* Lean`Internal`ToLean[expr, "Int"|"Rat"|"Real"|None, coreOnly, {passHeads...}(, closed)]:
 * the strict translation as a String, or $Failed with LeanForm::untrans.  closed -> True
 * sets LeanFormOpts.closed. */
Expr* builtin_lean_tolean(Expr* res);

#endif /* MATHILDA_PRINT_LEAN_H */
