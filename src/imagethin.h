/* imagethin.h -- Thinning and Pruning.
 *
 * Separate from imagefilter.c because neither is a neighbourhood filter: both delete pixels in
 * repeated passes and stop when a pass changes nothing, so the result of one pass is the input to
 * the next. A kernel that reads a window and writes one value cannot express that.
 */
#ifndef MATHILDA_IMAGETHIN_H
#define MATHILDA_IMAGETHIN_H

void imagethin_init(void);

#endif /* MATHILDA_IMAGETHIN_H */
