/*
 * charcode.c - character codes and character ranges.
 *
 *   ToCharacterCode[s]           list of the Unicode code points of s
 *   ToCharacterCode[{s1, ...}]   one list per string
 *   FromCharacterCode[n]         the one-character string for code point n
 *   FromCharacterCode[{n1, ...}] the string of those code points
 *   FromCharacterCode[{{..}, ..}] one string per sublist
 *   ToCharacterCode[s, enc] / FromCharacterCode[codes, enc]
 *                                the same in encoding enc: "UTF8" (bytes),
 *                                "ASCII", "ISOLatin1", "Unicode"
 *   CharacterRange["a", "e"]     {"a", "b", "c", "d", "e"}
 *   CharacterRange[n1, n2]       the characters with codes n1 .. n2
 *   Alphabet[]                   the 26 lowercase English letters
 *
 * Array surfaces. A string's codes are a machine-integer vector, so
 * ToCharacterCode writes them straight into a packed int64 buffer (above the
 * packing threshold) instead of building one Integer node per character, and
 * FromCharacterCode -- on the packed-aware list in src/pack.c -- reads a packed
 * or NDArray int64 buffer directly instead of having the gate box every code.
 * Neither has a Compile[] lowering: one side of each is a String, which no
 * compiled register holds (recorded in tools/check_compile_coverage.py).
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "utf8.h"
#include "sym_names.h"
#include "message.h"
#include "print.h"
#include "pack.h"
#include "ndarray.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ------------------------------------------------------------------ */
/* Character encodings                                                */
/* ------------------------------------------------------------------ */

/* The encodings the two-argument forms accept. ENC_UNICODE is the default
 * (code points); ENC_UTF8 is the byte sequence of the UTF-8 form, which is
 * also Mathilda's internal string representation; ASCII and ISOLatin1 are the
 * single-byte encodings, where a character outside the encoding's range has
 * no code (None from ToCharacterCode). */
typedef enum { ENC_UNICODE, ENC_UTF8, ENC_ASCII, ENC_LATIN1, ENC_UNKNOWN } CharEnc;

static CharEnc enc_of(const Expr* e) {
    if (e->type != EXPR_STRING) return ENC_UNKNOWN;
    const char* s = e->data.string;
    if (!strcmp(s, "UTF8") || !strcmp(s, "UTF-8")) return ENC_UTF8;
    if (!strcmp(s, "Unicode")) return ENC_UNICODE;
    if (!strcmp(s, "ASCII") || !strcmp(s, "PrintableASCII")) return ENC_ASCII;
    if (!strcmp(s, "ISOLatin1") || !strcmp(s, "ISO8859-1")) return ENC_LATIN1;
    return ENC_UNKNOWN;
}

/* The largest code point a single-byte encoding represents. */
static uint32_t enc_max(CharEnc enc) {
    return enc == ENC_ASCII ? 0x7F : 0xFF;
}

/* ------------------------------------------------------------------ */
/* ToCharacterCode                                                    */
/* ------------------------------------------------------------------ */

/* Codes of one string in a non-default encoding. UTF8 is the string's own
 * bytes (a packed int64 vector when packing applies); ASCII / ISOLatin1 give
 * a List with None for each character the encoding cannot represent. */
static Expr* codes_of_enc(const char* s, CharEnc enc) {
    size_t len = strlen(s);
    if (enc == ENC_UTF8) {
        /* The string's own bytes, except U+0000: held internally as C0 80
         * (utf8.h), its UTF-8 byte is 0. */
        size_t m = 0;
        for (size_t i = 0; i < len; i++, m++)
            if ((unsigned char)s[i] == 0xC0 && i + 1 < len && (unsigned char)s[i + 1] == 0x80) i++;
        int64_t* buf = NULL;
        Expr* packed = ndbuild_open_i64((int64_t)m, &buf);
        Expr* out = packed ? NULL : expr_new_function(expr_new_symbol(SYM_List), NULL, m);
        if (!packed && !out) return NULL;
        for (size_t i = 0, k = 0; i < len; i++, k++) {
            int64_t b = (int64_t)(unsigned char)s[i];
            if (b == 0xC0 && i + 1 < len && (unsigned char)s[i + 1] == 0x80) { b = 0; i++; }
            if (packed) buf[k] = b;
            else out->data.function.args[k] = expr_new_integer(b);
        }
        return packed ? packed : out;
    }
    size_t n = u8_length(s, len), k = 0;
    uint32_t hi = enc_max(enc);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, n);
    if (!out) return NULL;
    for (size_t i = 0; i < len; ) {
        uint32_t cp = u8_decode(s, len, &i);
        out->data.function.args[k++] = cp <= hi ? expr_new_integer((int64_t)cp)
                                                : expr_new_symbol(SYM_None);
    }
    return out;
}

/* Codes of one string: a packed int64 vector when packing applies, else a
 * List of Integers. */
static Expr* codes_of(const char* s) {
    size_t len = strlen(s);
    size_t n = u8_length(s, len);
    int64_t* buf = NULL;
    Expr* packed = ndbuild_open_i64((int64_t)n, &buf);
    if (packed) {
        size_t k = 0;
        for (size_t i = 0; i < len; ) buf[k++] = (int64_t)u8_decode(s, len, &i);
        return packed;
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, n);
    if (!out) return NULL;
    size_t k = 0;
    for (size_t i = 0; i < len; )
        out->data.function.args[k++] = expr_new_integer((int64_t)u8_decode(s, len, &i));
    return out;
}

static Expr* codes_of_any(const char* s, CharEnc enc) {
    return enc == ENC_UNICODE ? codes_of(s) : codes_of_enc(s, enc);
}

Expr* builtin_tocharactercode(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    CharEnc enc = argc == 2 ? enc_of(res->data.function.args[1]) : ENC_UNICODE;
    if (enc == ENC_UNKNOWN) return NULL;            /* unsupported encoding */
    Expr* a = res->data.function.args[0];
    if (a->type == EXPR_STRING) return codes_of_any(a->data.string, enc);

    /* A List or a LinkedList of strings; the answer in its image. */
    if (!strings_is_seq(a) || !strings_seq_all_strings(a))
        return strings_arg_message("ToCharacterCode", "strse",
                                   "String or list of strings", 1, res);

    LLReader sv;
    ll_rd_open(a, &sv);
    size_t m = sv.n;
    Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
    size_t built = 0;
    for (; out && built < m; built++) {
        out[built] = codes_of_any(ll_rd_at(&sv, built)->data.string, enc);
        if (!out[built]) break;
    }
    ll_rd_close(&sv);
    Expr* r = NULL;
    if (out && built == m) r = ll_seq_like(a, out, m);
    else for (size_t k = 0; out && k < built; k++) expr_free(out[k]);
    free(out);
    return r;
}

/* ------------------------------------------------------------------ */
/* FromCharacterCode                                                  */
/* ------------------------------------------------------------------ */

/* A code is valid in 0 .. 0x10FFFF, surrogates excluded. */
static int code_ok(int64_t c) {
    return c >= 0 && c <= 0x10FFFF && !(c >= 0xD800 && c <= 0xDFFF);
}

/* Report an invalid code: a negative or non-integer one is ::intnm on the
 * call; an out-of-range one is ::notunicode, naming the code list. */
static Expr* fcc_bad(Expr* res, const Expr* code, size_t pos, const Expr* list) {
    if (code->type == EXPR_INTEGER && code->data.integer >= 0) {
        char* shown = expr_to_string_outputform((Expr*)(list ? list : code));
        mth_message("FromCharacterCode", "notunicode",
                    "A character code, which should be a non-negative integer "
                    "less than 1114112, is expected at position %zu in %s%s%s.",
                    pos, list ? "" : "{", shown ? shown : "", list ? "" : "}");
        free(shown);
        return NULL;
    }
    return strings_arg_message("FromCharacterCode", "intnm",
                               "Non-negative machine-sized integer", pos, res);
}

/* Length of the UTF-8 sequence led by byte b (0 for a byte that cannot lead),
 * and the minimum code point that length may encode (overlong rejection). */
static int u8_lead_len(int64_t b, uint32_t* min_cp) {
    if (b < 0x80) { *min_cp = 0;       return 1; }
    if (b < 0xC0) return 0;                          /* continuation byte */
    if (b < 0xE0) { *min_cp = 0x80;    return 2; }
    if (b < 0xF0) { *min_cp = 0x800;   return 3; }
    if (b < 0xF8) { *min_cp = 0x10000; return 4; }
    return 0;
}

/* Decode n UTF-8 bytes into code points (cps has room for n). A byte that does
 * not start a well-formed sequence is reported ($CharacterEncoding::utf8, as a
 * sequence of one byte) and kept as the code point of the same value, so the
 * decode always completes. Returns the code-point count, or (size_t)-1 when a
 * value is not a byte at all (*bad is its index). */
static size_t utf8_decode_codes(const int64_t* b, size_t n, int64_t* cps, size_t* bad) {
    size_t m = 0;
    for (size_t i = 0; i < n; ) {
        if (b[i] < 0 || b[i] > 0xFF) { *bad = i; return (size_t)-1; }
        uint32_t min_cp = 0;
        int len = u8_lead_len(b[i], &min_cp);
        int ok = len > 0 && i + (size_t)len <= n;
        uint32_t cp = 0;
        if (ok) {
            cp = len == 1 ? (uint32_t)b[i]
                          : (uint32_t)b[i] & (0x7Fu >> (unsigned)len);
            for (int j = 1; ok && j < len; j++) {
                int64_t c = b[i + (size_t)j];
                if (c < 0x80 || c > 0xBF) ok = 0;
                else cp = (cp << 6) | ((uint32_t)c & 0x3F);
            }
            if (ok && (cp < min_cp || cp > 0x10FFFF || (cp >= 0xD800 && cp <= 0xDFFF)))
                ok = 0;
        }
        if (!ok) {
            mth_message("$CharacterEncoding", "utf8",
                        "The byte sequence {%lld} could not be interpreted as a "
                        "character in the UTF-8 character encoding.",
                        (long long)b[i]);
            cps[m++] = b[i];
            i++;
            continue;
        }
        cps[m++] = (int64_t)cp;
        i += (size_t)len;
    }
    return m;
}

static Expr* string_of_codes(const int64_t* codes, size_t n, Expr* res, const Expr* list);

/* The string of n codes in encoding `enc`: UTF8 bytes are decoded to code
 * points first; ISOLatin1 / ASCII codes are code points already, restricted to
 * the encoding's range. NULL when a code is invalid (message emitted) or, for
 * ASCII, outside 0 .. 127 (left unevaluated: no character to name). */
static Expr* string_of_codes_enc(const int64_t* codes, size_t n, CharEnc enc,
                                 Expr* res, const Expr* list) {
    if (enc == ENC_UNICODE) return string_of_codes(codes, n, res, list);
    if (enc == ENC_UTF8) {
        int64_t* cps = malloc((n ? n : 1) * sizeof(int64_t));
        if (!cps) return NULL;
        size_t bad = 0;
        size_t m = utf8_decode_codes(codes, n, cps, &bad);
        Expr* s = NULL;
        if (m == (size_t)-1) {
            Expr* c = expr_new_integer(codes[bad]);
            fcc_bad(res, c, bad + 1, list);
            expr_free(c);
        } else {
            s = string_of_codes(cps, m, res, list);
        }
        free(cps);
        return s;
    }
    for (size_t k = 0; k < n; k++) {
        if (codes[k] < 0) {
            Expr* c = expr_new_integer(codes[k]);
            fcc_bad(res, c, k + 1, list);
            expr_free(c);
            return NULL;
        }
        if (codes[k] > (int64_t)enc_max(enc)) return NULL;
    }
    return string_of_codes(codes, n, res, list);
}

/* The string of the n codes read by get(k); NULL (message emitted) if one is
 * invalid. `res`/`list` are only used for diagnostics. */
static Expr* string_of_codes(const int64_t* codes, size_t n, Expr* res, const Expr* list) {
    char* buf = malloc(4 * n + 1);
    if (!buf) return NULL;
    size_t o = 0;
    for (size_t k = 0; k < n; k++) {
        if (!code_ok(codes[k])) {
            free(buf);
            Expr* c = expr_new_integer(codes[k]);
            fcc_bad(res, c, k + 1, list);
            expr_free(c);
            return NULL;
        }
        o += u8_encode((uint32_t)codes[k], buf + o);
    }
    buf[o] = '\0';
    return expr_new_string_adopt(buf);
}

/* A List -- or LinkedList -- of Integer codes -> String; NULL if not all are
 * valid codes. An int64 LinkedList's slots are the codes already (gathered raw,
 * as the packed buffer is read in place). */
static Expr* string_of_list(const Expr* lst, CharEnc enc, Expr* res) {
    const LinkedList* ll = expr_linkedlist(lst);
    if (ll && ll_kind(ll) == LL_I64) {
        int64_t n = ll_length(ll);
        int64_t* codes = malloc((size_t)(n ? n : 1) * sizeof(int64_t));
        if (!codes) return NULL;
        ll_gather(ll, 0, n, codes);
        Expr* s = string_of_codes_enc(codes, (size_t)n, enc, res, lst);
        free(codes);
        return s;
    }
    LLReader e;
    ll_rd_open(lst, &e);
    size_t n = e.n;
    int64_t* codes = malloc((n ? n : 1) * sizeof(int64_t));
    if (!codes) { ll_rd_close(&e); return NULL; }
    for (size_t k = 0; k < n; k++) {
        const Expr* c = ll_rd_at(&e, k);
        if (c->type != EXPR_INTEGER || !code_ok(c->data.integer)) {
            free(codes);
            Expr* r = fcc_bad(res, c, k + 1, lst);
            ll_rd_close(&e);
            return r;
        }
        codes[k] = c->data.integer;
    }
    ll_rd_close(&e);
    Expr* s = string_of_codes_enc(codes, n, enc, res, lst);
    free(codes);
    return s;
}

/* A List or a LinkedList of codes (or of code lists). */
static int is_list(const Expr* e) {
    return e->type == EXPR_LINKEDLIST ||
           (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
            e->data.function.head->data.symbol.name == SYM_List);
}

/* FromCharacterCode of the value `a`; `res` is the call, for diagnostics. */
static Expr* fcc_value(Expr* a, CharEnc enc, Expr* res) {
    if (a->type == EXPR_INTEGER && enc != ENC_UNICODE) {   /* one encoded code */
        if (!code_ok(a->data.integer)) return fcc_bad(res, a, 1, NULL);
        int64_t c = a->data.integer;
        return string_of_codes_enc(&c, 1, enc, res, NULL);
    }
    if (a->type == EXPR_INTEGER) {
        if (!code_ok(a->data.integer)) return fcc_bad(res, a, 1, NULL);
        char buf[5];
        size_t n = u8_encode((uint32_t)a->data.integer, buf);
        return strings_char_expr(buf, n);
    }

    /* An int64 buffer (packed List or visible NDArray): read it directly. */
    if (a->type == EXPR_NDARRAY && a->data.ndarray.dtype == NDT_INT64 &&
        (a->data.ndarray.rank == 1 || a->data.ndarray.rank == 2)) {
        const NDArrayData* nd = &a->data.ndarray;
        const int64_t* codes = (const int64_t*)nd->data;
        if (nd->rank == 1)
            return string_of_codes_enc(codes, (size_t)nd->dims[0], enc, res, a);
        size_t rows = (size_t)nd->dims[0], cols = (size_t)nd->dims[1];
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, rows);
        for (size_t r = 0; out && r < rows; r++) {
            out->data.function.args[r] = string_of_codes_enc(codes + r * cols, cols, enc, res, a);
            if (!out->data.function.args[r]) { expr_free(out); out = NULL; }
        }
        return out;
    }
    if (a->type == EXPR_NDARRAY) {                /* other dtypes: as a List */
        Expr* plain = pack_unpack(a);
        if (!plain) return NULL;
        Expr* r = fcc_value(plain, enc, res);
        expr_free(plain);
        return r;
    }

    if (!is_list(a)) return fcc_bad(res, a, 1, NULL);

    /* {{codes}, {codes}, ...}: one string per sublist, in the image of the
     * outer list. A typed LinkedList holds codes, never sublists. */
    const LinkedList* al = expr_linkedlist(a);
    if (al && ll_kind(al) != LL_BOXED) return string_of_list(a, enc, res);
    LLReader e;
    ll_rd_open(a, &e);
    size_t m = e.n;
    if (m > 0 && is_list(ll_rd_at(&e, 0))) {
        Expr** out = malloc(sizeof(Expr*) * m);
        size_t built = 0;
        for (; out && built < m; built++) {
            const Expr* sub = ll_rd_at(&e, built);
            out[built] = is_list(sub) ? string_of_list(sub, enc, res)
                                      : fcc_bad(res, sub, built + 1, NULL);
            if (!out[built]) break;
        }
        ll_rd_close(&e);
        Expr* r = NULL;
        if (out && built == m) r = ll_seq_like(a, out, m);
        else for (size_t k = 0; out && k < built; k++) expr_free(out[k]);
        free(out);
        return r;
    }
    ll_rd_close(&e);
    return string_of_list(a, enc, res);
}

Expr* builtin_fromcharactercode(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    CharEnc enc = argc == 2 ? enc_of(res->data.function.args[1]) : ENC_UNICODE;
    if (enc == ENC_UNKNOWN) return NULL;            /* unsupported encoding */
    return fcc_value(res->data.function.args[0], enc, res);
}

/* ------------------------------------------------------------------ */
/* CharacterRange / Alphabet                                          */
/* ------------------------------------------------------------------ */

/* The code point of a one-character string or a non-negative integer, or -1. */
static long range_end(const Expr* e) {
    if (e->type == EXPR_INTEGER && code_ok(e->data.integer)) return (long)e->data.integer;
    if (e->type == EXPR_STRING) {
        const char* s = e->data.string;
        size_t len = strlen(s), i = 0;
        if (len == 0) return -1;
        uint32_t cp = u8_decode(s, len, &i);
        return i == len ? (long)cp : -1;
    }
    return -1;
}

static Expr* char_list(uint32_t lo, uint32_t hi) {
    size_t n = hi >= lo ? (size_t)(hi - lo + 1) : 0;
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, n);
    for (size_t k = 0; out && k < n; k++) {
        char buf[5];
        uint32_t cp = lo + (uint32_t)k;
        size_t len = u8_encode(cp, buf);
        if (len == 0) { buf[0] = '?'; len = 1; }          /* surrogate gap */
        out->data.function.args[k] = strings_char_expr(buf, len);
    }
    return out;
}

Expr* builtin_characterrange(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    long lo = range_end(a), hi = range_end(b);
    /* Both ends must be the same kind: two characters or two codes. */
    if (lo < 0 || hi < 0 || a->type != b->type) {
        char* sa = expr_to_string_outputform(a);
        char* sb = expr_to_string_outputform(b);
        mth_message("CharacterRange", "argtype",
                    "Arguments %s and %s of CharacterRange should be either "
                    "non-negative integers or one-character strings.",
                    sa ? sa : "", sb ? sb : "");
        free(sa);
        free(sb);
        return NULL;
    }
    return char_list((uint32_t)lo, (uint32_t)hi);
}

/* The lowercase letters of a named alphabet, as code-point runs [lo, hi]:
 * "English" a-z; "Greek" alpha-omega without the final sigma (24 letters);
 * "Russian" the 32 letters a-ya with yo after ye (33). */
typedef struct { const char* name; uint32_t runs[3][2]; int nruns; } AlphabetDef;
static const AlphabetDef k_alphabets[] = {
    { "English", { { 'a', 'z' } }, 1 },
    { "Greek",   { { 0x3B1, 0x3C1 }, { 0x3C3, 0x3C9 } }, 2 },
    { "Russian", { { 0x430, 0x435 }, { 0x451, 0x451 }, { 0x436, 0x44F } }, 3 },
};

Expr* builtin_alphabet(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return char_list('a', 'z');
    if (argc != 1 || res->data.function.args[0]->type != EXPR_STRING) return NULL;
    const char* lang = res->data.function.args[0]->data.string;
    for (size_t i = 0; i < sizeof k_alphabets / sizeof k_alphabets[0]; i++) {
        const AlphabetDef* d = &k_alphabets[i];
        if (strcmp(lang, d->name) != 0) continue;
        size_t n = 0;
        for (int r = 0; r < d->nruns; r++) n += d->runs[r][1] - d->runs[r][0] + 1;
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, n);
        size_t k = 0;
        for (int r = 0; out && r < d->nruns; r++)
            for (uint32_t cp = d->runs[r][0]; cp <= d->runs[r][1]; cp++) {
                char buf[5];
                size_t len = u8_encode(cp, buf);
                out->data.function.args[k++] = strings_char_expr(buf, len);
            }
        return out;
    }
    return NULL;                       /* an alphabet not tabulated here */
}
