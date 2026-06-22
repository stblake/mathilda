#ifndef MATHILDA_GRAPHICS_SHOW_H
#define MATHILDA_GRAPHICS_SHOW_H

#include "expr.h"

/* Show[graphics, opts...] -- merges trailing options into `graphics`'s
 * option list (rightmost/later wins), renders it, and returns the
 * (possibly re-styled) Graphics[...] expression. */
Expr* builtin_show(Expr* res);

/* Render `graphics_expr` (a Graphics[...] expression) in a blocking
 * Raylib window; control returns to the caller once the window is
 * closed. When USE_GRAPHICS is not compiled in, prints a one-line
 * graceful-degrade message instead and returns immediately. Does not
 * take ownership of `graphics_expr`. */
void graphics_show(const Expr* graphics_expr);

#endif /* MATHILDA_GRAPHICS_SHOW_H */
