/* image.h -- raster images: the representation and its accessors.
 *
 * The foundation of the image subsystem, and deliberately only that. Every filter needs a
 * validated image with known dimensions, a channel count and a pixel type, so those come
 * first; ImageConvolve and the rest are built on this and land separately.
 */
#ifndef MATHILDA_IMAGE_H
#define MATHILDA_IMAGE_H

#include <stdbool.h>
#include <stddef.h>
#include "expr.h"

void image_init(void);
void imagefilter_init(void);
void imagegeom_init(void);

/* Pixel types, in Wolfram's names. The type fixes the RANGE of a stored value, which is what
 * makes ImageData's scaling well defined: "Bit" is {0, 1}, "Byte" is 0..255, "Real" is already
 * the unit interval. */
typedef enum { IMG_BIT, IMG_BYTE, IMG_REAL } ImgType;

/* Read the canonical form Image[data, type] -- dimensions, channels and type, without copying
 * the pixels. `channels` is 1 for a rank-2 (grey) image, otherwise the last dimension.
 * False for anything that is not a canonical image. */
bool image_info(const Expr* e, size_t* width, size_t* height,
                size_t* channels, ImgType* type);

/* Pixels as a flat height*width*channels buffer of UNIT-INTERVAL reals, row-major with
 * channels innermost -- the same interleaved order the nested form uses. The caller frees.
 *
 * Every filter needs this, and needs it in unit scale rather than raw stored values: a filter's
 * arithmetic is defined on brightnesses, so a "Byte" image must become 0..1 before a kernel
 * touches it or the kernel's own scale would silently depend on the input type. */
bool image_load(const Expr* img, size_t* width, size_t* height, size_t* channels,
                double** buf);

/* Build a canonical Image[data, "Real"] from such a buffer.
 *
 * The result is always "Real", never the input's type. A filter produces values that are not
 * generally representable in the input type -- a Gaussian of bytes is not a byte -- and
 * rounding back would throw away precision the caller did not ask to lose. */
Expr* image_build_real(const double* buf, size_t width, size_t height, size_t channels);

#endif /* MATHILDA_IMAGE_H */
