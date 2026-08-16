/* imagecolor.h -- ColorReplace, ColorQuantize, HistogramTransform.
 *
 * Separate from imagefilter.c because each needs a GLOBAL pass before it can decide anything per
 * pixel: a replacement needs a distance rule, a quantisation needs a palette drawn from every
 * pixel, and an equalisation needs the whole distribution. Gather, then write -- which is not the
 * shape of a neighbourhood filter.
 */
#ifndef MATHILDA_IMAGECOLOR_H
#define MATHILDA_IMAGECOLOR_H

void imagecolor_init(void);

#endif /* MATHILDA_IMAGECOLOR_H */
