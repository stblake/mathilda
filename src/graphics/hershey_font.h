/* hershey_font.h — a minimal single-stroke vector font for graphics text.
 *
 * NOTE: despite the "Hershey-style" name (drawn as line strokes, not
 * filled glyph outlines, exactly like the classic Hershey vector fonts),
 * this is a small hand-authored stick alphabet, not a transcription of
 * the historical Hershey datasets. It covers digits, uppercase letters
 * (lowercase input folds to the uppercase glyph), and basic punctuation
 * -- enough for axis tick labels, AxesLabel/PlotLabel, and Text[]. This
 * is a deliberate placeholder pending proper font rendering (tracked in
 * docs/spec/builtins/graphics.md). Only built when USE_GRAPHICS is on. */
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
