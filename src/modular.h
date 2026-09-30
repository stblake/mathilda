#ifndef MODULAR_H
#define MODULAR_H

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

Expr* builtin_module(Expr* res);
Expr* builtin_block(Expr* res);
Expr* builtin_with(Expr* res);
Expr* builtin_unique(Expr* res);

/* Capture-avoiding pattern substitution support, used by replace_bindings.
 * expr_is_binding_scope: is `e` a scoping construct (Module/Block/With/Function/
 * Table) that actually binds names? scoping_capture_avoid: if substituting a
 * value carrying any name in danger[] into `e`'s body could capture one of e's
 * bound locals, return a fresh copy with the colliding locals alpha-renamed,
 * else NULL. See MATHILDA_DIVERGENCES.md A11. */
bool expr_is_binding_scope(Expr* e);
Expr* scoping_capture_avoid(Expr* e, const char** danger, size_t ndanger);

/* Block's dynamic-scope frames, exposed for the TimeConstrained unwind only.
 * Block restores its locals' OwnValues, DownValues and attributes on the way
 * out, but a timeout siglongjmps past every C frame in between, so
 * tc_run_guarded records the depth before the body and drains back to it after
 * a jump -- otherwise a Block-installed temporary rule (a rewrite hook, say)
 * would stay installed for the rest of the session. Same pattern as the
 * async-defer count and the message-suppression depth. */
int  mth_block_depth_save(void);
void mth_block_depth_unwind(int depth);

void modular_init(void);

#endif // MODULAR_H
