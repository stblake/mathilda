/* label_font.h — filled sans-serif text (Helvetica/Arial-like) for graphics
 * text, replacing the retired Hershey single-stroke vector font.
 *
 * Discovers a system TrueType font (Arial on macOS, Liberation Sans/DejaVu
 * Sans on Linux, Arial on Windows) at label_font_load() time and falls back
 * to Raylib's embedded GetFontDefault() bitmap font when none is found --
 * guaranteed on a minimal CI/Docker Linux box with no fonts installed at
 * all. Only built when USE_GRAPHICS is on. */
#ifndef MATHILDA_GRAPHICS_LABEL_FONT_H
#define MATHILDA_GRAPHICS_LABEL_FONT_H

#include "raylib.h"

/* The old Hershey vector font used a 7-unit cap height, with every caller's
 * `scale` meaning "world/screen units per cap height". label_font_draw*
 * keeps that exact convention (deriving an equivalent pixel fontSize from
 * the loaded font's own cap-height metric) so none of render.c's centering/
 * offset math -- built around this constant -- needs to change. */
#define LABEL_FONT_CAP_HEIGHT 7.0f

/* Loads the font pair (regular + bold) for the Raylib GPU context that is
 * currently open. A Font's texture doesn't survive CloseWindow, and this
 * codebase opens four independent, non-nested windows (graphics_show,
 * graphics3d_show, Animate, Manipulate), so every one of them must call
 * this once right after InitWindow and label_font_unload() once right
 * before its own CloseWindow. Idempotent no-op if already loaded. */
void label_font_load(void);

/* Releases the loaded font's GPU texture. No-op if nothing is loaded, and
 * safe to call even when discovery fell back to GetFontDefault() (which
 * Raylib owns and must not be unloaded here). */
void label_font_unload(void);

/* Total advance width of `text` at the given scale (same cap-height
 * convention as label_font_draw), with no rotation applied -- use this to
 * center/right-align a label before drawing. Always measured against the
 * regular weight, matching the old Hershey behaviour where advance widths
 * didn't depend on stroke weight either. */
float label_font_width(const char* text, float scale);

/* Draws `text` with its baseline starting at (x, y), sized so its cap
 * height equals scale * LABEL_FONT_CAP_HEIGHT pixels, and the whole string
 * rotated `rotation_deg` about (x, y). */
void label_font_draw(const char* text, float x, float y, float scale,
                      float rotation_deg, Color color);

/* As label_font_draw, but uses the bold face when weight > 1.0 (regular
 * otherwise). Replaces the old Hershey hairline "thickness" parameter --
 * filled glyphs don't need stroke-width emulation the way 1px hairline
 * strokes did, so weight is a simple regular/bold switch instead. */
void label_font_draw_ex(const char* text, float x, float y, float scale,
                         float rotation_deg, Color color, float weight);

/* Drop-in replacements for Raylib's own MeasureText/DrawText (regular
 * weight, top-left anchored, no rotation, fontSize in literal pixels) --
 * used by render3d.c, render.c's own HUD/toolbar chrome, and Animate/
 * Manipulate's slider and button UI, none of which use the cap-height
 * convention above. */
int  label_font_measure_px(const char* text, int fontSize);
void label_font_draw_px(const char* text, int x, int y, int fontSize, Color color);

#endif /* MATHILDA_GRAPHICS_LABEL_FONT_H */
