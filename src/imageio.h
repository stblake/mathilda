/* imageio.h -- reading and writing raster image FILES.
 *
 * Kept apart from image.c, which owns the representation and knows nothing about formats. The
 * split matters because the decoders are vendored third-party code: everything that includes
 * them lives in exactly one translation unit, so the rest of the subsystem stays free of them.
 */
#ifndef MATHILDA_IMAGEIO_H
#define MATHILDA_IMAGEIO_H

void imageio_init(void);

#endif /* MATHILDA_IMAGEIO_H */
