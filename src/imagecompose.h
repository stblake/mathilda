/* imagecompose.h -- the alpha channel, and putting images together.
 *
 * Kept apart from imagefilter.c because none of these is a neighbourhood operation: they combine
 * or reshape whole images rather than reading a window around each pixel, and they share one
 * channel-promotion rule that the filters have no use for.
 */
#ifndef MATHILDA_IMAGECOMPOSE_H
#define MATHILDA_IMAGECOMPOSE_H

void imagecompose_init(void);

#endif /* MATHILDA_IMAGECOMPOSE_H */
