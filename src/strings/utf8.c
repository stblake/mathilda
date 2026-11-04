/*
 * utf8.c - character (code point) view of UTF-8 strings; see utf8.h.
 */

#include "utf8.h"

#include <stdlib.h>
#include <string.h>

#define HI8 UINT64_C(0x8080808080808080)

int u8_is_ascii(const char* s, size_t n) {
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {             /* 8 bytes per step */
        uint64_t w;
        memcpy(&w, s + i, 8);
        if (w & HI8) return 0;
    }
    for (; i < n; i++)
        if ((unsigned char)s[i] & 0x80) return 0;
    return 1;
}

static int is_cont(unsigned char c) { return (c & 0xC0) == 0x80; }

size_t u8_char_len(const char* str, size_t n, size_t i) {
    const unsigned char* s = (const unsigned char*)str;
    unsigned char b0 = s[i];
    size_t rem = n - i;
    if (b0 < 0x80) return 1;
    if (b0 >= 0xC2 && b0 <= 0xDF)
        return (rem >= 2 && is_cont(s[i + 1])) ? 2 : 1;
    if (b0 >= 0xE0 && b0 <= 0xEF) {
        if (rem < 3) return 1;
        unsigned char b1 = s[i + 1];
        unsigned char lo = (b0 == 0xE0) ? 0xA0 : 0x80;   /* no overlong */
        unsigned char hi = (b0 == 0xED) ? 0x9F : 0xBF;   /* no surrogates */
        return (b1 >= lo && b1 <= hi && is_cont(s[i + 2])) ? 3 : 1;
    }
    if (b0 >= 0xF0 && b0 <= 0xF4) {
        if (rem < 4) return 1;
        unsigned char b1 = s[i + 1];
        unsigned char lo = (b0 == 0xF0) ? 0x90 : 0x80;
        unsigned char hi = (b0 == 0xF4) ? 0x8F : 0xBF;   /* <= U+10FFFF */
        return (b1 >= lo && b1 <= hi && is_cont(s[i + 2]) && is_cont(s[i + 3])) ? 4 : 1;
    }
    return 1;                                /* stray / invalid byte */
}

size_t u8_length(const char* s, size_t n) {
    size_t count = 0, i = 0;
    while (i < n) {
        /* Skip an ASCII run 8 bytes at a time. */
        while (i + 8 <= n) {
            uint64_t w;
            memcpy(&w, s + i, 8);
            if (w & HI8) break;
            i += 8;
            count += 8;
        }
        if (i >= n) break;
        if (!((unsigned char)s[i] & 0x80)) { i++; count++; continue; }
        i += u8_char_len(s, n, i);
        count++;
    }
    return count;
}

size_t u8_offset(const char* s, size_t n, size_t k) {
    size_t i = 0;
    while (k > 0 && i < n) {
        i += u8_char_len(s, n, i);
        k--;
    }
    return i;
}

int u8_view_init(U8View* v, const char* s, size_t n) {
    v->s = s;
    v->nbytes = n;
    v->off = NULL;
    if (u8_is_ascii(s, n)) { v->nchars = n; return 0; }
    size_t* off = malloc((n + 1) * sizeof(size_t));   /* <= one char per byte */
    if (!off) return -1;
    size_t k = 0, i = 0;
    while (i < n) {
        off[k++] = i;
        i += u8_char_len(s, n, i);
    }
    off[k] = n;
    v->nchars = k;
    v->off = off;
    return 0;
}

void u8_view_free(U8View* v) {
    free(v->off);
    v->off = NULL;
}

uint32_t u8_decode(const char* str, size_t n, size_t* i) {
    const unsigned char* s = (const unsigned char*)str;
    size_t k = u8_char_len(str, n, *i);
    const unsigned char* p = s + *i;
    *i += k;
    switch (k) {
        case 2:  return ((uint32_t)(p[0] & 0x1F) << 6) | (p[1] & 0x3F);
        case 3:  return ((uint32_t)(p[0] & 0x0F) << 12) | ((uint32_t)(p[1] & 0x3F) << 6) |
                        (p[2] & 0x3F);
        case 4:  return ((uint32_t)(p[0] & 0x07) << 18) | ((uint32_t)(p[1] & 0x3F) << 12) |
                        ((uint32_t)(p[2] & 0x3F) << 6) | (p[3] & 0x3F);
        default: return p[0];                 /* ASCII, or a malformed byte */
    }
}

size_t u8_encode(uint32_t cp, char* out) {
    unsigned char* o = (unsigned char*)out;
    if (cp < 0x80) { o[0] = (unsigned char)cp; return 1; }
    if (cp < 0x800) {
        o[0] = (unsigned char)(0xC0 | (cp >> 6));
        o[1] = (unsigned char)(0x80 | (cp & 0x3F));
        return 2;
    }
    if (cp >= 0xD800 && cp <= 0xDFFF) return 0;
    if (cp < 0x10000) {
        o[0] = (unsigned char)(0xE0 | (cp >> 12));
        o[1] = (unsigned char)(0x80 | ((cp >> 6) & 0x3F));
        o[2] = (unsigned char)(0x80 | (cp & 0x3F));
        return 3;
    }
    if (cp <= 0x10FFFF) {
        o[0] = (unsigned char)(0xF0 | (cp >> 18));
        o[1] = (unsigned char)(0x80 | ((cp >> 12) & 0x3F));
        o[2] = (unsigned char)(0x80 | ((cp >> 6) & 0x3F));
        o[3] = (unsigned char)(0x80 | (cp & 0x3F));
        return 4;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Named characters \[Name]                                           */
/* ------------------------------------------------------------------ */

typedef struct { const char* name; uint32_t cp; } U8Named;

/* Sorted by name (strcmp order); generated from Mathematica's
 * UnicodeCharacters table by tools/gen_named_chars.py. */
static const U8Named k_named[] = {
#include "named_chars.inc"
};

long u8_named_char(const char* name, size_t len) {
    size_t lo = 0, hi = sizeof k_named / sizeof k_named[0];
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        const char* m = k_named[mid].name;
        int c = strncmp(m, name, len);
        if (c == 0 && m[len] != '\0') c = 1;       /* m is longer: m > name */
        if (c == 0) return (long)k_named[mid].cp;
        if (c < 0) lo = mid + 1;
        else       hi = mid;
    }
    return -1;
}

/* ------------------------------------------------------------------ */
/* Source-code classes of non-ASCII characters                        */
/* ------------------------------------------------------------------ */

/* Sorted by code point; only the exceptions to "a letter" are listed. */
static const U8CharInfo k_char_class[] = {
#include "named_char_class.inc"
};

void u8_for_each_class_head(void (*fn)(const char* head)) {
    for (size_t i = 0; i < sizeof k_char_class / sizeof k_char_class[0]; i++)
        if (k_char_class[i].head) fn(k_char_class[i].head);
}

const U8CharInfo* u8_char_class(uint32_t cp) {
    size_t lo = 0, hi = sizeof k_char_class / sizeof k_char_class[0];
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (k_char_class[mid].cp == cp) return &k_char_class[mid];
        if (k_char_class[mid].cp < cp) lo = mid + 1;
        else                            hi = mid;
    }
    return NULL;
}
