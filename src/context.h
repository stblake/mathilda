#ifndef CONTEXT_H
#define CONTEXT_H

#include "expr.h"
#include "symtab.h"
#include <stdbool.h>
#include <stddef.h>

/* Context system for Mathilda.
 *
 * A context is a namespace identified by a string ending in a backtick, e.g.
 * "Global`" or "MyPackage`Private`". Every symbol conceptually lives in a
 * context, addressed as "context`name".
 *
 * Implementation notes:
 * - $Context holds the current context string; new unqualified symbols are
 *   created in it.
 * - $ContextPath is an ordered list of contexts consulted when resolving a
 *   bare identifier that has not yet been seen.
 * - For backward compatibility with existing code that stores builtins under
 *   bare names (e.g. "Plus"), bare names that already resolve to an existing
 *   builtin are left bare; they are treated as implicitly belonging to
 *   "System`".
 */

/* Initialize/shut down the context subsystem. context_init() must be called
 * from core_init() before any .m file is loaded. */
void context_init(void);
void context_shutdown(void);

/* Current $Context value (read-only borrow; never NULL after context_init). */
const char* context_current(void);

/* Path access (read-only). */
size_t context_path_count(void);
const char* context_path_at(size_t i);

/* Resolve a parsed identifier to its canonical symbol-table name.
 *
 * Rules applied in order:
 *  1. Names starting with '$' are returned unchanged (system variables).
 *  2. A leading backtick is treated as relative to the current context:
 *     "`Private`f"  in context "foo`"  becomes "foo`Private`f".
 *  3. Names containing a backtick elsewhere are absolute. "System`X" is
 *     canonicalized to bare "X" (since builtins are stored under bare names).
 *  4. Otherwise (bare identifier):
 *     a. Search each prefix in { $Context } U $ContextPath. If prefix+name
 *        already exists in the symbol table, return that full name.
 *     b. If a bare symbol with a builtin function is already registered,
 *        return bare (implicit "System`").
 *     c. If $Context == "Global`", return bare (legacy behavior).
 *     d. Otherwise, qualify with $Context.
 *
 * Returned string is malloc'd; caller must free(). */
char* context_resolve_name(const char* raw);

/* Display helper: return a pointer to the short name for printing if the
 * symbol's context prefix is currently visible (equal to $Context or on
 * $ContextPath). Otherwise return `name` unchanged. The pointer is either
 * `name` or points inside it; never NULL; never allocated. */
const char* context_display_name(const char* name);

/* Built-in operations (also exposed for use outside the Context[] builtin). */
bool context_begin(const char* ctx);               /* Begin["ctx`"] */
bool context_begin_package(const char* ctx,
                           const char** needs,
                           size_t needs_count);    /* BeginPackage */
/* On success, *out_ended_ctx is set to a newly malloc'd string containing
 * the just-closed context (caller frees). */
bool context_end(char** out_ended_ctx);
bool context_end_package(char** out_ended_ctx);

/* Produce Expr* copies of current state, for use as $Context/$ContextPath
 * OwnValues. Caller owns the returned expressions. */
Expr* context_as_string(void);
Expr* context_path_as_list(void);

#endif /* CONTEXT_H */
