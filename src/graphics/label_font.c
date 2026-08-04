/* label_font.c — filled sans-serif text via a discovered system TrueType
 * font, replacing hershey_font.c.
 *
 * Font discovery tries an ordered, platform-specific list of {regular,
 * bold} TTF path pairs and keeps the first pair whose regular file loads.
 * Nothing is bundled: this project vendors no font assets (see
 * src/external/), matching how raylib/GMP-ECM/LAPACK/FLINT/PCRE2/FFTW are
 * all autodetected system libraries rather than vendored, and just like
 * those, the feature degrades gracefully -- here to Raylib's embedded
 * GetFontDefault() bitmap font -- rather than failing when nothing is
 * found, which is the guaranteed case on a minimal CI/Docker Linux image. */

#include "label_font.h"
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>

/* Loaded at a generous pixel size so drawing at any fontSize this renderer
 * actually asks for (tick labels through PlotLabel to Manipulate's large
 * data-driven Text[] scales) stays a downscale of the atlas, not a blurry
 * upscale. */
#define LABEL_FONT_BASE_SIZE 96

typedef struct {
    const char* regular;
    const char* bold;
} FontCandidate;

#if defined(__APPLE__)
static const FontCandidate CANDIDATES[] = {
    { "/System/Library/Fonts/Supplemental/Arial.ttf",
      "/System/Library/Fonts/Supplemental/Arial Bold.ttf" },
};
#elif defined(_WIN32)
static const FontCandidate CANDIDATES[] = {
    { "C:\\Windows\\Fonts\\arial.ttf", "C:\\Windows\\Fonts\\arialbd.ttf" },
};
#else
/* Linux: Liberation Sans is metric-compatible with Arial/Helvetica and is
 * the common default (fonts-liberation); DejaVu Sans is the near-universal
 * fallback pulled in by X11/fontconfig on most distros. Neither is
 * guaranteed -- see the GetFontDefault() fallback below. */
static const FontCandidate CANDIDATES[] = {
    { "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
      "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf" },
    { "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
      "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" },
    { "/usr/share/fonts/TTF/DejaVuSans.ttf",
      "/usr/share/fonts/TTF/DejaVuSans-Bold.ttf" },
};
#endif

static Font s_regular;
static Font s_bold;
static bool s_loaded = false;
static bool s_have_custom = false; /* false => both weights are GetFontDefault() */
static float s_cap_k = 1.0f;       /* baseSize / (regular 'H' ink height in px) */
static float s_baseline_at_base = 0.0f; /* top-of-line -> baseline, px, at baseSize */

static bool file_exists(const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) return false;
    fclose(f);
    return true;
}

static bool font_ok(Font f) {
    return f.texture.id != 0 && f.glyphCount > 0;
}

void label_font_load(void) {
    if (s_loaded) return;

    s_have_custom = false;
    size_t n = sizeof(CANDIDATES) / sizeof(CANDIDATES[0]);
    for (size_t i = 0; i < n; i++) {
        if (!file_exists(CANDIDATES[i].regular)) continue;
        Font reg = LoadFontEx(CANDIDATES[i].regular, LABEL_FONT_BASE_SIZE, NULL, 0);
        if (!font_ok(reg)) continue;

        Font bold = reg;
        if (file_exists(CANDIDATES[i].bold)) {
            Font b = LoadFontEx(CANDIDATES[i].bold, LABEL_FONT_BASE_SIZE, NULL, 0);
            if (font_ok(b)) bold = b;
        }
        s_regular = reg;
        s_bold = bold;
        s_have_custom = true;
        break;
    }
    if (!s_have_custom) {
        s_regular = GetFontDefault();
        s_bold = GetFontDefault();
    }

    /* Derive the cap-height-to-fontSize ratio and the line-top-to-baseline
     * offset from a reference glyph ('H': flat-bottomed, no descender), so
     * label_font_draw's `scale` produces the exact pixel cap height the old
     * Hershey convention did, and its baseline lands where callers expect. */
    int idx = GetGlyphIndex(s_regular, 'H');
    GlyphInfo gi = s_regular.glyphs[idx];
    float cap_height_px = (float)gi.image.height;
    if (cap_height_px < 1.0f) cap_height_px = (float)s_regular.baseSize * 0.7f;
    s_cap_k = (float)s_regular.baseSize / cap_height_px;
    s_baseline_at_base = (float)(gi.offsetY + gi.image.height);

    s_loaded = true;
}

void label_font_unload(void) {
    if (!s_loaded) return;
    if (s_have_custom) {
        UnloadFont(s_regular);
        if (s_bold.texture.id != s_regular.texture.id) UnloadFont(s_bold);
    }
    s_loaded = false;
}

static float font_size_for_scale(float scale) {
    return scale * LABEL_FONT_CAP_HEIGHT * s_cap_k;
}

float label_font_width(const char* text, float scale) {
    if (!s_loaded) label_font_load();
    float font_size = font_size_for_scale(scale);
    return MeasureTextEx(s_regular, text, font_size, 0.0f).x;
}

void label_font_draw_ex(const char* text, float x, float y, float scale,
                         float rotation_deg, Color color, float weight) {
    if (!s_loaded) label_font_load();
    Font font = (weight > 1.0f) ? s_bold : s_regular;
    float font_size = font_size_for_scale(scale);
    float baseline = s_baseline_at_base * (font_size / (float)s_regular.baseSize);
    Vector2 origin = { 0.0f, baseline };
    Vector2 position = { x, y };
    DrawTextPro(font, text, position, origin, rotation_deg, font_size, 0.0f, color);
}

void label_font_draw(const char* text, float x, float y, float scale,
                      float rotation_deg, Color color) {
    label_font_draw_ex(text, x, y, scale, rotation_deg, color, 1.0f);
}

int label_font_measure_px(const char* text, int fontSize) {
    if (!s_loaded) label_font_load();
    return (int)MeasureTextEx(s_regular, text, (float)fontSize, 0.0f).x;
}

void label_font_draw_px(const char* text, int x, int y, int fontSize, Color color) {
    if (!s_loaded) label_font_load();
    DrawTextEx(s_regular, text, (Vector2){ (float)x, (float)y }, (float)fontSize, 0.0f, color);
}
