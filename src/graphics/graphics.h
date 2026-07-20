/* graphics.h — umbrella header for the 2D graphics engine (src/graphics/).
 *
 * External callers continue to `#include "graphics.h"` unchanged; this
 * file re-exports every per-module public header plus the graphics_init()
 * entry point that core_init() calls. */
#ifndef MATHILDA_GRAPHICS_H
#define MATHILDA_GRAPHICS_H

#include "primitives.h"
#include "sampling.h"
#include "plot.h"
#include "plot3d.h"
#include "parametricplot.h"
#include "parametricplot3d.h"
#include "listplot.h"
#include "polarplot.h"
#include "streamplot.h"
#include "contourplot.h"
#include "densityplot.h"
#include "barchart.h"
#include "vectorplot.h"
#include "show.h"
#include "animate.h"
#include "complexplot.h"

void graphics_init(void);

#endif /* MATHILDA_GRAPHICS_H */
