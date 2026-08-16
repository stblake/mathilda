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

/* Read the canonical Image3D[data, type]. `depth` is the number of SLICES, and the storage order is
 * depth x height x width (slices outermost, indexed data[[z, y, x]]) while ImageDimensions reports
 * {width, height, depth} -- fully reversed, which is Mathematica's convention and the 3-D version
 * of the same trap the 2-D accessors carry. */
bool image3d_info(const Expr* e, size_t* width, size_t* height, size_t* depth,
                  size_t* channels, ImgType* type);

/* Voxels as a flat depth*height*width*channels buffer of unit-interval reals. Caller frees. */
bool image3d_load(const Expr* img, size_t* width, size_t* height, size_t* depth,
                  size_t* channels, double** buf);

/* Pixels as a flat height*width*channels buffer of UNIT-INTERVAL reals, row-major with channels
 * innermost -- the same interleaved order the nested form uses. The caller frees.
 *
 * Every filter needs this, and needs it in unit scale rather than raw stored values: a filter's
 * arithmetic is defined on brightnesses, so a "Byte" image must become 0..1 before a kernel touches
 * it or the kernel's own scale would silently depend on the input type. */
bool image_load(const Expr* img, size_t* width, size_t* height, size_t* channels,
                double** buf);

/* Build a canonical Image[data, "Real"] from such a buffer.
 *
 * The result is always "Real", never the input's type. A filter produces values that are not
 * generally representable in the input type -- a Gaussian of bytes is not a byte -- and rounding
 * back would throw away precision the caller did not ask to lose. */
Expr* image_build_real(const double* buf, size_t width, size_t height, size_t channels);

/* Build a canonical Image3D[data, "Real"] from a depth*height*width*channels buffer.
 *
 * Removed once for having no caller and restored with the 3-D convolution that uses it -- untested
 * unused code is worse than absent code. */
Expr* image3d_build_real(const double* buf, size_t width, size_t height, size_t depth,
                         size_t channels);

/* Build a canonical Image[data, "Bit"] from a 0/1 mask, packed.
 *
 * "Bit" rather than "Real" because these are DECISIONS, not brightnesses: a thinning or a
 * binarisation answers yes-or-no per pixel, and ImageData reports stored values, so a real-typed 1
 * would print as `1.` where Mathematica prints `1`. Defined in imagefilter.c, where the packed
 * builder and its nested-List fallback already live. */
Expr* image_build_bit(const unsigned char* mask, size_t width, size_t height);

/* Serialise an Image or Image3D for the notebook front end.
 *
 * Returns a malloc'd JSON object -- {"w":W,"h":H,"data":"<base64 RGBA>"} plus, for a volume, the depth
 * and which slice was sent -- or NULL when `e` is not an image. The caller frees it.
 *
 * RGBA base64 rather than a JSON array of numbers: a 512x512 colour image is 786432 samples, which as
 * decimal text is about 3 MB per evaluation, where base64 of the RGBA bytes is 1.4 MB and is exactly
 * what a canvas ImageData wants -- no per-pixel work in JavaScript at all.
 *
 * A VOLUME SENDS ITS SIX BOUNDARY FACES plus the middle slice, under "faces" and "data". Three faces
 * are visible from any viewpoint, so six is everything an opaque rotatable block can ever show, and
 * for a 256^3 volume that is 6 * 65536 samples rather than 16.7 million. The middle slice stays as
 * well, for a consumer that can only draw one plane. Internal structure is NOT sent: an opaque box
 * cannot display it, and a cutting plane or a translucent rendering would need the volume itself. */
char* image_to_json(const Expr* e);

#endif /* MATHILDA_IMAGE_H */
