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

void modular_init(void);

#endif // MODULAR_H
