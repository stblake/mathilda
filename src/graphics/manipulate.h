/* manipulate.h — Manipulate[expr, {u, umin, umax}, ...] */
#ifndef MATHILDA_GRAPHICS_MANIPULATE_H
#define MATHILDA_GRAPHICS_MANIPULATE_H

#include "expr.h"

/* Builtin for Manipulate[expr, {u, umin, umax}, ...].
 * HoldAll: expr is held unevaluated; it is re-evaluated each frame with
 * every control variable bound to its current value.  Opens a Raylib
 * window with one row per control (continuous slider or discrete button
 * set), each independently user-driven, and renders until the window is
 * closed, then returns Null. */
Expr* builtin_manipulate(Expr* res);

#endif /* MATHILDA_GRAPHICS_MANIPULATE_H */
