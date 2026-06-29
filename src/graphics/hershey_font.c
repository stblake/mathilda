/* hershey_font.c — single-stroke vector text using the classic Hershey
 * Roman Simplex font.
 *
 * Glyph data is transcribed from the historical Hershey occidental dataset
 * (tools/hershey.dat) by tools/gen_hershey.py into hershey_glyphs.inc. Each
 * glyph is a list of stroke vertices in raw Hershey integer coordinates
 * (value = char - 'R'), with y increasing downward, baseline at y=9 and cap
 * top at y=-12 (design cap height 21). A vertex with x == HF_PENUP lifts the
 * pen: the next vertex begins a new disconnected stroke. Per-glyph side
 * bearings give proportional advance widths.
 *
 * At draw time we normalise into the renderer's local frame: a 7-unit cap
 * height (HERSHEY_CAP_HEIGHT in render.c) with baseline at y=0 and +y up, so
 * every caller's `scale` and the renderer's positioning math are unchanged.
 * Only built when USE_GRAPHICS is on. */

#include "hershey_font.h"
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "hershey_glyphs.inc"

/* Raw-Hershey -> local-frame conversion. The raw design cap height is 21
 * (baseline 9, cap top -12); we map it onto the 7-unit cap height the renderer
 * assumes, so 'A' (raw advance 18) lands at local advance 6. */
#define HF_BASELINE   9.0f
#define HF_DESIGN_CAP 21.0f
#define HF_NORM       (7.0f / HF_DESIGN_CAP)

/* Global stroke-weight multiplier applied to every glyph stroke. The base
 * (hairline) weight is 1px, so 2.0 draws all text strokes twice as thick. */
#define HERSHEY_STROKE_SCALE 2.0f

static const hf_glyph* find_glyph(unsigned char c) {
    if (c < 32 || c > 126) return NULL;
    return &hf_ascii[c - 32];
}

/* Local advance width of one character (0 for unrenderable codes). */
static float glyph_advance(unsigned char c) {
    const hf_glyph* g = find_glyph(c);
    if (!g) return 0.0f;
    return (float)(g->right - g->left) * HF_NORM;
}

float hershey_text_width(const char* text, float scale) {
    float w = 0.0f;
    for (const unsigned char* p = (const unsigned char*)text; *p; p++) {
        w += glyph_advance(*p);
    }
    return w * scale;
}

void hershey_draw_text_ex(const char* text, float x, float y, float scale,
                          float rotation_deg, Color color, float thickness) {
    float theta = rotation_deg * (float)M_PI / 180.0f;
    float ct = cosf(theta), st = sinf(theta);
    float cursor = 0.0f;
    thickness *= HERSHEY_STROKE_SCALE;
    /* Above a hairline, draw quads with a disc at each vertex so the stroke
     * joins stay closed (DrawLineEx butt-caps would otherwise notch them). */
    bool thick = thickness > 1.05f;
    float cap_r = thickness * 0.5f;

    for (const unsigned char* p = (const unsigned char*)text; *p; p++) {
        const hf_glyph* g = find_glyph(*p);
        if (*p == '.') {
            /* A period is a dot, not a stroke -- a near-zero-length segment is
             * invisible (or a stray sliver) at small scale, so draw it as an
             * actual filled circle on the baseline instead. */
            float adv = glyph_advance('.');
            float lx = cursor + 0.5f * adv;
            float ly = 0.3f; /* just above the baseline, local units */
            Vector2 center = { x + (lx * ct - ly * st) * scale,
                               y - (lx * st + ly * ct) * scale };
            DrawCircleV(center, 0.45f * scale, color);
            cursor += adv;
            continue;
        }
        if (g && g->verts) {
            bool pen_down = false;
            Vector2 prev = { 0, 0 };
            for (int i = 0; i < g->nverts; i++) {
                hf_vtx v = g->verts[i];
                if (v.x == HF_PENUP) { pen_down = false; continue; }
                float lx = cursor + (float)(v.x - g->left) * HF_NORM;
                float ly = (HF_BASELINE - (float)v.y) * HF_NORM;
                Vector2 cur = { x + (lx * ct - ly * st) * scale,
                                y - (lx * st + ly * ct) * scale };
                if (pen_down) {
                    if (thick) {
                        DrawLineEx(prev, cur, thickness, color);
                        DrawCircleV(prev, cap_r, color);
                        DrawCircleV(cur, cap_r, color);
                    } else {
                        DrawLineV(prev, cur, color);
                    }
                }
                prev = cur;
                pen_down = true;
            }
        }
        cursor += glyph_advance(*p);
    }
}

void hershey_draw_text(const char* text, float x, float y, float scale,
                        float rotation_deg, Color color) {
    hershey_draw_text_ex(text, x, y, scale, rotation_deg, color, 1.0f);
}
