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

/* Pixel types, in Wolfram's names. The type fixes the RANGE of a stored value, which is what
 * makes ImageData's scaling well defined: "Bit" is {0, 1}, "Byte" is 0..255, "Real" is already
 * the unit interval. */
typedef enum { IMG_BIT, IMG_BYTE, IMG_REAL } ImgType;

/* Read the canonical form Image[data, type] -- dimensions, channels and type, without copying
 * the pixels. `channels` is 1 for a rank-2 (grey) image, otherwise the last dimension.
 * False for anything that is not a canonical image. */
bool image_info(const Expr* e, size_t* width, size_t* height,
                size_t* channels, ImgType* type);

#endif /* MATHILDA_IMAGE_H */
