/* hershey_font.h — single-stroke vector font for graphics text.
 *
 * This is a faithful transcription of the classic Hershey *Roman Simplex*
 * occidental font (drawn as line strokes, not filled outlines, exactly like
 * the historical Hershey vector fonts). Glyph data is generated from
 * tools/hershey.dat by tools/gen_hershey.py into hershey_glyphs.inc and covers
 * the full printable ASCII range (digits, upper- and lowercase letters,
 * punctuation and symbols) with proportional advance widths -- used for axis
 * tick labels, AxesLabel/PlotLabel, and Text[]. Only built when USE_GRAPHICS
 * is on. */
#ifndef MATHILDA_GRAPHICS_HERSHEY_FONT_H
#define MATHILDA_GRAPHICS_HERSHEY_FONT_H

#include "raylib.h"

/* Draws `text` with its baseline starting at (x, y), each glyph's native
 * 4x7 (cap-height) grid scaled by `scale` and the whole string rotated
 * `rotation_deg` degrees about (x, y). Unrecognized characters render as
 * blank space (just advance the cursor) rather than erroring. */
void hershey_draw_text(const char* text, float x, float y, float scale,
                        float rotation_deg, Color color);

/* As hershey_draw_text, but draws each stroke `thickness` pixels wide (with
 * rounded joins so thick glyphs stay clean). thickness <= ~1 falls back to the
 * hairline path. hershey_draw_text() is exactly this with thickness 1. */
void hershey_draw_text_ex(const char* text, float x, float y, float scale,
                          float rotation_deg, Color color, float thickness);

/* Total advance width of `text` at the given scale, with no rotation
 * applied -- use this to center/right-align a label before drawing. */
float hershey_text_width(const char* text, float scale);

#endif /* MATHILDA_GRAPHICS_HERSHEY_FONT_H */
