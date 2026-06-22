/* primitives.h — small shared value types for the graphics engine.
 *
 * Deliberately free of any Expr/symtab dependency so sampling.c stays a
 * pure-numeric module that can be unit tested without the rest of
 * Mathilda's evaluator, and so render.c has simple POD types to draw. */
#ifndef MATHILDA_GRAPHICS_PRIMITIVES_H
#define MATHILDA_GRAPHICS_PRIMITIVES_H

typedef struct {
    double x, y;
} Point2D;

typedef struct {
    double xmin, xmax, ymin, ymax;
} PlotRange2D;

typedef struct {
    unsigned char r, g, b, a;
} RGBA8;

#endif /* MATHILDA_GRAPHICS_PRIMITIVES_H */
