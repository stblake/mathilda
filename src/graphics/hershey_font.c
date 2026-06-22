/* hershey_font.c — see hershey_font.h for the honesty note on naming.
 *
 * Each glyph is a flat array of points on a 4-wide x 7-tall (cap height)
 * grid, baseline at y=0. A point at the PU sentinel means "lift the pen"
 * -- the next point starts a new disconnected stroke instead of
 * continuing the previous line. */

#include "hershey_font.h"
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { float x, y; } HPt;
#define PU {-1000.0f, -1000.0f}
static int is_pu(HPt p) { return p.x <= -999.0f; }

#define GLYPH_WIDTH   4.0f
#define GLYPH_HEIGHT  7.0f
#define GLYPH_ADVANCE 6.0f /* includes inter-glyph spacing */

/* ---- Digits ---- */
static const HPt g_0[] = { {1,0},{3,0},{4,1},{4,6},{3,7},{1,7},{0,6},{0,1},{1,0} };
static const HPt g_1[] = { {1,5.5f},{2,7}, PU, {2,7},{2,0}, PU, {1,0},{3,0} };
static const HPt g_2[] = { {0,5},{0,6},{1,7},{3,7},{4,6},{4,4.5f},{0,1},{0,0},{4,0} };
static const HPt g_3[] = { {0,6},{1,7},{3,7},{4,6},{4,4},{2,3.5f},{4,3},{4,1},{3,0},{1,0},{0,1} };
static const HPt g_4[] = { {3,0},{3,7},{0,2.5f},{4,2.5f} };
static const HPt g_5[] = { {4,7},{0,7},{0,4},{2,4},{4,3},{4,1},{3,0},{1,0},{0,1} };
static const HPt g_6[] = { {3,7},{1,7},{0,5},{0,1},{1,0},{3,0},{4,1},{4,3},{3,4},{0,3.5f} };
static const HPt g_7[] = { {0,7},{4,7},{1,0} };
static const HPt g_8[] = { {1,0},{3,0},{4,1},{4,3},{3,3.5f},{4,4},{4,6},{3,7},{1,7},{0,6},{0,4},{1,3.5f},{0,3},{0,1},{1,0} };
static const HPt g_9[] = { {3,3.5f},{1,3.5f},{0,4.5f},{0,6},{1,7},{3,7},{4,6},{4,1},{3,0},{1,0} };

/* ---- Uppercase letters ---- */
static const HPt g_A[] = { {0,0},{2,7},{4,0}, PU, {1,2.5f},{3,2.5f} };
static const HPt g_B[] = { {0,0},{0,7},{3,7},{4,6},{3,3.5f},{0,3.5f}, PU, {3,3.5f},{4,1},{3,0},{0,0} };
static const HPt g_C[] = { {4,6},{3,7},{1,7},{0,6},{0,1},{1,0},{3,0},{4,1} };
static const HPt g_D[] = { {0,0},{0,7},{2,7},{4,5},{4,2},{2,0},{0,0} };
static const HPt g_E[] = { {4,0},{0,0},{0,7},{4,7}, PU, {0,3.5f},{3,3.5f} };
static const HPt g_F[] = { {0,0},{0,7},{4,7}, PU, {0,3.5f},{3,3.5f} };
static const HPt g_G[] = { {4,6},{3,7},{1,7},{0,6},{0,1},{1,0},{3,0},{4,1},{4,3.5f},{2.5f,3.5f} };
static const HPt g_H[] = { {0,0},{0,7}, PU, {4,0},{4,7}, PU, {0,3.5f},{4,3.5f} };
static const HPt g_I[] = { {0,7},{4,7}, PU, {2,7},{2,0}, PU, {0,0},{4,0} };
static const HPt g_J[] = { {3,7},{3,1},{2,0},{1,0},{0,1} };
static const HPt g_K[] = { {0,0},{0,7}, PU, {4,7},{0,3.5f},{4,0} };
static const HPt g_L[] = { {0,7},{0,0},{4,0} };
static const HPt g_M[] = { {0,0},{0,7},{2,3},{4,7},{4,0} };
static const HPt g_N[] = { {0,0},{0,7},{4,0},{4,7} };
static const HPt g_O[] = { {1,0},{3,0},{4,1},{4,6},{3,7},{1,7},{0,6},{0,1},{1,0} };
static const HPt g_P[] = { {0,0},{0,7},{3,7},{4,6},{3,3.5f},{0,3.5f} };
static const HPt g_Q[] = { {1,0},{3,0},{4,1},{4,6},{3,7},{1,7},{0,6},{0,1},{1,0}, PU, {2,1.5f},{4,-0.3f} };
static const HPt g_R[] = { {0,0},{0,7},{3,7},{4,6},{3,3.5f},{0,3.5f}, PU, {2,3.5f},{4,0} };
static const HPt g_S[] = { {4,6},{3,7},{1,7},{0,6},{0,4.2f},{4,2.8f},{4,1},{3,0},{1,0},{0,1} };
static const HPt g_T[] = { {0,7},{4,7}, PU, {2,7},{2,0} };
static const HPt g_U[] = { {0,7},{0,1},{1,0},{3,0},{4,1},{4,7} };
static const HPt g_V[] = { {0,7},{2,0},{4,7} };
static const HPt g_W[] = { {0,7},{1,0},{2,4},{3,0},{4,7} };
static const HPt g_X[] = { {0,0},{4,7}, PU, {0,7},{4,0} };
static const HPt g_Y[] = { {0,7},{2,3.5f},{4,7}, PU, {2,3.5f},{2,0} };
static const HPt g_Z[] = { {0,7},{4,7},{0,0},{4,0} };

/* ---- Punctuation ---- */
static const HPt g_period[]  = { {2,0},{2,0.3f} };
static const HPt g_comma[]   = { {2,1},{1.3f,-0.6f} };
static const HPt g_minus[]   = { {0.5f,3.5f},{3.5f,3.5f} };
static const HPt g_plus[]    = { {2,1.5f},{2,5.5f}, PU, {0,3.5f},{4,3.5f} };
static const HPt g_colon[]   = { {2,5},{2,5.3f}, PU, {2,1.5f},{2,1.8f} };
static const HPt g_lparen[]  = { {3,7},{1,6},{1,1},{3,0} };
static const HPt g_rparen[]  = { {1,7},{3,6},{3,1},{1,0} };
static const HPt g_slash[]   = { {0,0},{4,7} };
static const HPt g_under[]   = { {0,-0.5f},{4,-0.5f} };
static const HPt g_apos[]    = { {2,7},{2.5f,5.5f} };
static const HPt g_bang[]    = { {2,7},{2,2.5f}, PU, {2,1},{2,0.6f} };
static const HPt g_qmark[]   = { {0.3f,6},{1,7},{3,7},{3.7f,6},{3.7f,4.5f},{2,3.5f},{2,2.5f}, PU, {2,1},{2,0.6f} };
static const HPt g_star[]    = { {2,2},{2,5}, PU, {0.5f,2.5f},{3.5f,4.5f}, PU, {0.5f,4.5f},{3.5f,2.5f} };
static const HPt g_equals[]  = { {0.5f,4.3f},{3.5f,4.3f}, PU, {0.5f,2.7f},{3.5f,2.7f} };

typedef struct {
    char       c;
    const HPt* pts;
    int        count;
} HEntry;

#define G(ch, arr) { ch, arr, (int)(sizeof(arr) / sizeof(arr[0])) }

static const HEntry g_table[] = {
    G('0', g_0), G('1', g_1), G('2', g_2), G('3', g_3), G('4', g_4),
    G('5', g_5), G('6', g_6), G('7', g_7), G('8', g_8), G('9', g_9),
    G('A', g_A), G('B', g_B), G('C', g_C), G('D', g_D), G('E', g_E),
    G('F', g_F), G('G', g_G), G('H', g_H), G('I', g_I), G('J', g_J),
    G('K', g_K), G('L', g_L), G('M', g_M), G('N', g_N), G('O', g_O),
    G('P', g_P), G('Q', g_Q), G('R', g_R), G('S', g_S), G('T', g_T),
    G('U', g_U), G('V', g_V), G('W', g_W), G('X', g_X), G('Y', g_Y),
    G('Z', g_Z),
    G('.', g_period), G(',', g_comma), G('-', g_minus), G('+', g_plus),
    G(':', g_colon), G('(', g_lparen), G(')', g_rparen), G('/', g_slash),
    G('_', g_under), G('\'', g_apos), G('!', g_bang), G('?', g_qmark),
    G('*', g_star), G('=', g_equals),
};
#define HTABLE_LEN (sizeof(g_table) / sizeof(g_table[0]))

static const HEntry* find_glyph(char c) {
    if (c >= 'a' && c <= 'z') c = (char)(c - 'a' + 'A');
    for (size_t i = 0; i < HTABLE_LEN; i++) {
        if (g_table[i].c == c) return &g_table[i];
    }
    return NULL;
}

float hershey_text_width(const char* text, float scale) {
    float w = 0.0f;
    for (const char* p = text; *p; p++) w += GLYPH_ADVANCE * scale;
    return w;
}

void hershey_draw_text(const char* text, float x, float y, float scale,
                        float rotation_deg, Color color) {
    float theta = rotation_deg * (float)M_PI / 180.0f;
    float ct = cosf(theta), st = sinf(theta);
    float cursor = 0.0f;

    for (const char* p = text; *p; p++) {
        if (*p == '.') {
            /* A period is a dot, not a stroke -- a near-zero-length line
             * segment is invisible (or a stray sliver) at small scale, so
             * draw it as an actual filled circle instead. */
            float lx = cursor + 2.0f * scale;
            float ly = 0.15f * scale;
            Vector2 center = { x + lx * ct - ly * st, y - (lx * st + ly * ct) };
            DrawCircleV(center, 0.45f * scale, color);
            cursor += GLYPH_ADVANCE * scale;
            continue;
        }
        const HEntry* g = (*p == ' ') ? NULL : find_glyph(*p);
        if (g) {
            bool pen_down = false;
            Vector2 prev = { 0, 0 };
            for (int i = 0; i < g->count; i++) {
                HPt pt = g->pts[i];
                if (is_pu(pt)) { pen_down = false; continue; }
                float lx = cursor + pt.x * scale;
                float ly = pt.y * scale;
                Vector2 cur = { x + lx * ct - ly * st, y - (lx * st + ly * ct) };
                if (pen_down) DrawLineV(prev, cur, color);
                prev = cur;
                pen_down = true;
            }
        }
        cursor += GLYPH_ADVANCE * scale;
    }
}
