/* animate.h — Animate[expr, {t, tmin, tmax}, opts...] */
#ifndef MATHILDA_GRAPHICS_ANIMATE_H
#define MATHILDA_GRAPHICS_ANIMATE_H

#include "expr.h"

/* Builtin for Animate[expr, {t, tmin, tmax}, opts...].
 * HoldAll: expr is held unevaluated; it is re-evaluated each frame with
 * t bound to the current parameter value.  Opens a Raylib window with a
 * playback control bar, renders frames until the window is closed, then
 * returns Null. */
Expr* builtin_animate(Expr* res);

#endif /* MATHILDA_GRAPHICS_ANIMATE_H */
