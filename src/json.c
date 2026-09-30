/* ===========================================================================
 * JSON reader and writer (see json.h)
 *
 * Reader.  A single forward pass over the bytes -- O(n) time, no
 * backtracking.  Values are built on one shared value stack: an array or
 * object pushes its children and, at its closing bracket, turns the top of
 * the stack into one List / Association node, so the parse makes one
 * allocation per value plus amortised stack growth, whatever the nesting.
 * Strings are decoded into a shared scratch buffer (escapes, \uXXXX, UTF-16
 * surrogate pairs -> UTF-8); \u0000 is stored as the two-byte "modified
 * UTF-8" sequence C0 80 because Mathilda strings are NUL-terminated, and the
 * writer turns it back into \u0000, so it round-trips.
 *
 * Numbers follow Mathematica's mapping: an integer (any size) is an Integer;
 * a number with a fraction part is a machine Real (an arbitrary-precision
 * real when it over- or underflows a double); a number with an exponent and
 * no fraction ("1e-2", "3E5") is the EXACT Integer / Rational it denotes, as
 * Mathematica reads it.  Leading zeros ("01") are accepted, as Mathematica
 * does.
 *
 * Errors print Mathematica's messages (Import::jsoninvalidtoken, ...,
 * followed by Import::jsonhintposandchar with the line:column of the
 * offending character) and the import gives $Failed.
 *
 * Writer.  Mathematica's layout: one element per line, tab indentation,
 * "key":value with no space, empty containers inline ([] / {}); "Compact" ->
 * True drops all whitespace.  Machine reals print in shortest round-trip form
 * with Mathematica's JSON exponent style (0.5, 1.0, 2.5e1, 1.0e-2: fixed
 * notation only for decimal exponents 0 and -1).  Strings escape ", \, /, the
 * named control escapes and other control bytes as \u00XX; non-ASCII UTF-8
 * passes through unchanged.
 * ======================================================================== */

#include "json.h"
#include "arithmetic.h"   /* make_rational_mpz */
#include "assoc.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "message.h"
#include "ndarray.h"
#include "print.h"
#include "sym_names.h"
#include "symtab.h"
#include "attr.h"

#include <ctype.h>
#include <gmp.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ----------------------------------------------------------- messages */

/* Print `Head::tag: text` unless messages are Quiet-ed; always noted so an
 * enclosing Check sees it. */
static void json_msg(const char* head, const char* tag, const char* text) {
    mth_msg_note_fired();
    if (!mth_msg_suppressed()) fprintf(stderr, "%s::%s: %s\n", head, tag, text);
}

static Expr* failed(void) { return expr_new_symbol(SYM_DollarFailed); }

/* Shared, refcounted leaf/head nodes: JSON produces the same few symbols
 * millions of times, and sharing one node per symbol costs a refcount bump
 * instead of an intern lookup + allocation. */
static Expr *g_list, *g_assoc, *g_rule, *g_true, *g_false, *g_null;

static void json_shared_init(void) {
    if (g_list) return;
    g_list  = expr_new_symbol(SYM_List);
    g_assoc = expr_new_symbol(SYM_Association);
    g_rule  = expr_new_symbol(SYM_Rule);
    g_true  = expr_new_symbol(SYM_True);
    g_false = expr_new_symbol(SYM_False);
    g_null  = expr_new_symbol(SYM_Null);
}

/* ================================================================ reader */

#define JSON_MAX_DEPTH 2000

typedef struct {
    const char* s;
    size_t n, pos;
    bool raw;
    int depth;
    /* value stack */
    Expr** stk;
    size_t sp, scap;
    /* string scratch */
    char* buf;
    size_t blen, bcap;
    /* error */
    const char* etag;
    char emsg[192];
    size_t epos;
} JP;

static void jp_error(JP* p, const char* tag, const char* msg, size_t at) {
    if (p->etag) return;                      /* keep the first error */
    p->etag = tag;
    snprintf(p->emsg, sizeof p->emsg, "%s", msg);
    p->epos = at;
}

static bool jp_push(JP* p, Expr* e) {
    if (!e) { jp_error(p, "jsonmem", "Out of memory.", p->pos); return false; }
    if (p->sp == p->scap) {
        size_t nc = p->scap ? p->scap * 2 : 256;
        Expr** ns = realloc(p->stk, sizeof(Expr*) * nc);
        if (!ns) { expr_free(e); jp_error(p, "jsonmem", "Out of memory.", p->pos); return false; }
        p->stk = ns;
        p->scap = nc;
    }
    p->stk[p->sp++] = e;
    return true;
}

static bool buf_put(JP* p, const char* s, size_t k) {
    if (p->blen + k + 1 > p->bcap) {
        size_t nc = p->bcap ? p->bcap * 2 : 256;
        char* nb;
        while (nc < p->blen + k + 1) nc *= 2;
        nb = realloc(p->buf, nc);
        if (!nb) return false;
        p->buf = nb;
        p->bcap = nc;
    }
    memcpy(p->buf + p->blen, s, k);
    p->blen += k;
    return true;
}

static void skip_ws(JP* p) {
    const char* s = p->s;
    size_t i = p->pos, n = p->n;
    while (i < n && (s[i] == ' ' || s[i] == '\n' || s[i] == '\r' || s[i] == '\t')) i++;
    p->pos = i;
}

static int hexval(int c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

/* Read 4 hex digits at p->pos (after "\u"); -1 on failure. */
static long read_hex4(JP* p) {
    long v = 0;
    int i;
    if (p->pos + 4 > p->n) return -1;
    for (i = 0; i < 4; i++) {
        int h = hexval((unsigned char)p->s[p->pos + i]);
        if (h < 0) return -1;
        v = v * 16 + h;
    }
    p->pos += 4;
    return v;
}

static bool put_utf8(JP* p, unsigned long cp) {
    char b[4];
    size_t k;
    if (cp == 0) { b[0] = (char)0xC0; b[1] = (char)0x80; k = 2; }   /* modified UTF-8 */
    else if (cp < 0x80) { b[0] = (char)cp; k = 1; }
    else if (cp < 0x800) { b[0] = (char)(0xC0 | (cp >> 6)); b[1] = (char)(0x80 | (cp & 0x3F)); k = 2; }
    else if (cp < 0x10000) {
        b[0] = (char)(0xE0 | (cp >> 12)); b[1] = (char)(0x80 | ((cp >> 6) & 0x3F));
        b[2] = (char)(0x80 | (cp & 0x3F)); k = 3;
    } else {
        b[0] = (char)(0xF0 | (cp >> 18)); b[1] = (char)(0x80 | ((cp >> 12) & 0x3F));
        b[2] = (char)(0x80 | ((cp >> 6) & 0x3F)); b[3] = (char)(0x80 | (cp & 0x3F)); k = 4;
    }
    return buf_put(p, b, k);
}

/* Decode a string whose opening quote is at p->pos into p->buf (NUL-
 * terminated).  Returns false on error. */
static bool parse_string_raw(JP* p) {
    const char* s = p->s;
    size_t n = p->n;
    p->pos++;                                 /* opening quote */
    p->blen = 0;
    for (;;) {
        size_t start = p->pos, i = start;
        /* Fast run of plain bytes. */
        while (i < n && s[i] != '"' && s[i] != '\\') i++;
        if (i > start && !buf_put(p, s + start, i - start)) {
            jp_error(p, "jsonmem", "Out of memory.", i);
            return false;
        }
        p->pos = i;
        if (i >= n) {
            jp_error(p, "jsonfoundendofinput",
                     "Expecting JSON token '\"'. Found end of input.", n);
            return false;
        }
        if (s[i] == '"') { p->pos = i + 1; break; }
        /* Escape. */
        if (i + 1 >= n) {
            jp_error(p, "jsonfoundendofinput",
                     "Expecting JSON token '\"'. Found end of input.", n);
            return false;
        }
        {
            char c = s[i + 1], out = 0;
            p->pos = i + 2;
            switch (c) {
            case '"': out = '"'; break;
            case '\\': out = '\\'; break;
            case '/': out = '/'; break;
            case 'b': out = '\b'; break;
            case 'f': out = '\f'; break;
            case 'n': out = '\n'; break;
            case 'r': out = '\r'; break;
            case 't': out = '\t'; break;
            case 'u': {
                long cp = read_hex4(p), lo;
                if (cp < 0) {
                    jp_error(p, "jsoninvcodepoint",
                             "Expecting 4 hexadecimal digits after token \\u.", i);
                    return false;
                }
                if (cp >= 0xDC00 && cp <= 0xDFFF) {
                    jp_error(p, "jsonsurrogate",
                             "Invalid UTF-16 surrogate pair or stray surrogate code unit encountered.", i);
                    return false;
                }
                if (cp >= 0xD800 && cp <= 0xDBFF) {
                    if (p->pos + 1 < n && s[p->pos] == '\\' && s[p->pos + 1] == 'u') {
                        p->pos += 2;
                        lo = read_hex4(p);
                    } else {
                        lo = -1;
                    }
                    if (lo < 0xDC00 || lo > 0xDFFF) {
                        jp_error(p, "jsonsurrogate",
                                 "Invalid UTF-16 surrogate pair or stray surrogate code unit encountered.", i);
                        return false;
                    }
                    cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
                }
                if (!put_utf8(p, (unsigned long)cp)) {
                    jp_error(p, "jsonmem", "Out of memory.", i);
                    return false;
                }
                continue;
            }
            default: {
                char m[64];
                snprintf(m, sizeof m, "Unexpected escaped character '%c'.", c);
                jp_error(p, "jsoninvescchar", m, i + 1);
                return false;
            }
            }
            if (!buf_put(p, &out, 1)) { jp_error(p, "jsonmem", "Out of memory.", i); return false; }
        }
    }
    if (!buf_put(p, "", 1)) { jp_error(p, "jsonmem", "Out of memory.", p->pos); return false; }
    p->blen--;                                /* keep the NUL out of the length */
    return true;
}

static bool parse_value(JP* p);

/* 10^k as an mpz. */
static void mpz_pow10(mpz_t out, unsigned long k) { mpz_ui_pow_ui(out, 10, k); }

static bool parse_number(JP* p) {
    const char* s = p->s;
    size_t n = p->n, start = p->pos, i = start, int_start, int_end;
    bool neg = false, frac = false, expo = false;
    long ex = 0;
    if (s[i] == '-') { neg = true; i++; }
    int_start = i;
    while (i < n && s[i] >= '0' && s[i] <= '9') i++;
    int_end = i;
    if (int_end == int_start) { jp_error(p, "jsoninvalidnum", "Invalid number.", i); return false; }
    if (i < n && s[i] == '.') {
        frac = true;
        i++;
        while (i < n && s[i] >= '0' && s[i] <= '9') i++;
    }
    if (i < n && (s[i] == 'e' || s[i] == 'E')) {
        size_t es;
        bool eneg = false;
        expo = true;
        i++;
        if (i < n && (s[i] == '+' || s[i] == '-')) { eneg = s[i] == '-'; i++; }
        es = i;
        while (i < n && s[i] >= '0' && s[i] <= '9') {
            if (ex < 100000000L) ex = ex * 10 + (s[i] - '0');
            i++;
        }
        if (i == es) { jp_error(p, "jsoninvalidnum", "Invalid number.", i); return false; }
        if (eneg) ex = -ex;
    }
    p->pos = i;

    if (!frac && !expo) {
        size_t nd = int_end - int_start;
        if (nd <= 18) {
            int64_t v = 0;
            size_t k;
            for (k = int_start; k < int_end; k++) v = v * 10 + (s[k] - '0');
            return jp_push(p, expr_new_integer(neg ? -v : v));
        } else {
            char* tmp = malloc(nd + 2);
            Expr* e;
            if (!tmp) return jp_push(p, NULL);
            tmp[0] = '-';
            memcpy(tmp + 1, s + int_start, nd);
            tmp[nd + 1] = '\0';
            e = expr_new_bigint_from_str(neg ? tmp : tmp + 1);
            free(tmp);
            return jp_push(p, e ? expr_bigint_normalize(e) : NULL);
        }
    }
    if (!frac && expo && ex > -100000 && ex < 100000) {
        /* Exact: mantissa * 10^ex. */
        size_t nd = int_end - int_start;
        char* tmp = malloc(nd + 1);
        mpz_t m, t;
        Expr* e;
        if (!tmp) return jp_push(p, NULL);
        memcpy(tmp, s + int_start, nd);
        tmp[nd] = '\0';
        mpz_init_set_str(m, tmp, 10);
        free(tmp);
        if (neg) mpz_neg(m, m);
        mpz_init(t);
        if (ex >= 0) {
            mpz_pow10(t, (unsigned long)ex);
            mpz_mul(m, m, t);
            e = expr_bigint_normalize(expr_new_bigint_from_mpz(m));
        } else {
            mpz_pow10(t, (unsigned long)(-ex));
            e = make_rational_mpz(m, t);
        }
        mpz_clear(m);
        mpz_clear(t);
        return jp_push(p, e);
    }
    /* Inexact. */
    {
        size_t len = i - start;
        char small[64];
        char* tmp = len < sizeof small ? small : malloc(len + 1);
        double d;
        Expr* e = NULL;
        bool nonzero = false;
        size_t k;
        if (!tmp) return jp_push(p, NULL);
        memcpy(tmp, s + start, len);
        tmp[len] = '\0';
        d = strtod(tmp, NULL);
        for (k = int_start; k < i && s[k] != 'e' && s[k] != 'E'; k++)
            if (s[k] >= '1' && s[k] <= '9') { nonzero = true; break; }
        if (isinf(d) || (d == 0.0 && nonzero)) {
#ifdef USE_MPFR
            e = expr_new_mpfr_from_str(tmp, 53);
#endif
            if (!e) e = expr_new_real(d);
        } else {
            e = expr_new_real(d);
        }
        if (tmp != small) free(tmp);
        return jp_push(p, e);
    }
}

/* Match a literal (true / false / null) at p->pos. */
static bool parse_literal(JP* p, const char* word, Expr* node) {
    size_t k = strlen(word), j;
    for (j = 0; j < k; j++) {
        if (p->pos + j >= p->n) {
            char m[96];
            snprintf(m, sizeof m, "Expecting JSON token '%s'. Found end of input.", word);
            jp_error(p, "jsonfoundendofinput", m, p->n);
            return false;
        }
        if (p->s[p->pos + j] != word[j]) {
            jp_error(p, "jsoninvalidtoken", "Invalid token found.", p->pos);
            return false;
        }
    }
    p->pos += k;
    return jp_push(p, expr_copy(node));
}

/* Pop the top `cnt` stack entries into a fresh node with the given head. */
static Expr* pop_node(JP* p, Expr* head, size_t cnt) {
    Expr* e = expr_new_function(expr_copy(head), cnt ? p->stk + (p->sp - cnt) : NULL, cnt);
    p->sp -= cnt;
    return e;
}

/* True when the `cnt` Rule nodes at the top of the stack have distinct
 * (string) keys.  Small objects compare pairwise; larger ones are left to
 * assoc_from_rules' hash index. */
static bool rules_distinct_small(JP* p, size_t cnt) {
    size_t i, j, base = p->sp - cnt;
    for (i = 0; i < cnt; i++) {
        const char* ki = p->stk[base + i]->data.function.args[0]->data.string;
        for (j = i + 1; j < cnt; j++)
            if (strcmp(ki, p->stk[base + j]->data.function.args[0]->data.string) == 0)
                return false;
    }
    return true;
}

static bool parse_object(JP* p) {
    size_t base = p->sp;
    p->pos++;                                 /* '{' */
    skip_ws(p);
    if (p->pos < p->n && p->s[p->pos] == '}') {
        p->pos++;
        return jp_push(p, pop_node(p, p->raw ? g_assoc : g_list, 0));
    }
    for (;;) {
        Expr* key;
        Expr* rule;
        skip_ws(p);
        if (p->pos >= p->n) {
            jp_error(p, "jsonfoundendofinput", "Expecting JSON token '\"'. Found end of input.", p->n);
            return false;
        }
        if (p->s[p->pos] != '"') {
            jp_error(p, "jsonkeynstr", "Object keys must be strings.", p->pos);
            return false;
        }
        if (!parse_string_raw(p)) return false;
        key = expr_new_string(p->buf);
        skip_ws(p);
        if (p->pos >= p->n || p->s[p->pos] != ':') {
            expr_free(key);
            if (p->pos >= p->n)
                jp_error(p, "jsonfoundendofinput", "Expecting JSON token ':'. Found end of input.", p->n);
            else
                jp_error(p, "jsonkvsep", "Missing key value separator ':'.", p->pos);
            return false;
        }
        p->pos++;
        if (!parse_value(p)) { expr_free(key); return false; }
        {
            Expr* kv[2] = { key, p->stk[p->sp - 1] };
            rule = expr_new_function(expr_copy(g_rule), kv, 2);
            p->stk[p->sp - 1] = rule;
        }
        skip_ws(p);
        if (p->pos >= p->n) {
            jp_error(p, "jsonfoundendofinput", "Expecting JSON token '}'. Found end of input.", p->n);
            return false;
        }
        if (p->s[p->pos] == ',') { p->pos++; continue; }
        if (p->s[p->pos] == '}') { p->pos++; break; }
        jp_error(p, "jsonobjectmissingsep", "Expecting end of object or a value separator.", p->pos);
        return false;
    }
    {
        size_t cnt = p->sp - base;
        Expr* obj;
        if (!p->raw) {
            obj = pop_node(p, g_list, cnt);
        } else if (cnt <= 16 && rules_distinct_small(p, cnt)) {
            obj = pop_node(p, g_assoc, cnt);
        } else {
            /* Duplicate keys (or a large object): last value wins, first
             * position kept -- assoc_from_rules copies, so release ours. */
            size_t i;
            obj = assoc_from_rules(p->stk + base, cnt);
            for (i = base; i < p->sp; i++) expr_free(p->stk[i]);
            p->sp = base;
        }
        return jp_push(p, obj);
    }
}

static bool parse_array(JP* p) {
    size_t base = p->sp;
    p->pos++;                                 /* '[' */
    skip_ws(p);
    if (p->pos < p->n && p->s[p->pos] == ']') {
        p->pos++;
        return jp_push(p, pop_node(p, g_list, 0));
    }
    for (;;) {
        if (!parse_value(p)) return false;
        skip_ws(p);
        if (p->pos >= p->n) {
            jp_error(p, "jsonfoundendofinput", "Expecting JSON token ']'. Found end of input.", p->n);
            return false;
        }
        if (p->s[p->pos] == ',') { p->pos++; continue; }
        if (p->s[p->pos] == ']') { p->pos++; break; }
        jp_error(p, "jsonarraymissingsep", "Expecting end of array or a value separator.", p->pos);
        return false;
    }
    return jp_push(p, pop_node(p, g_list, p->sp - base));
}

static bool parse_value(JP* p) {
    char c;
    bool ok;
    skip_ws(p);
    if (p->pos >= p->n) {
        jp_error(p, "jsonfoundendofinput", "Expecting a JSON value. Found end of input.", p->n);
        return false;
    }
    if (++p->depth > JSON_MAX_DEPTH) {
        jp_error(p, "jsondepth", "Maximum nesting depth exceeded.", p->pos);
        return false;
    }
    c = p->s[p->pos];
    switch (c) {
    case '{': ok = parse_object(p); break;
    case '[': ok = parse_array(p); break;
    case '"':
        ok = parse_string_raw(p) && jp_push(p, expr_new_string(p->buf));
        break;
    case 't': ok = parse_literal(p, "true", g_true); break;
    case 'f': ok = parse_literal(p, "false", g_false); break;
    case 'n': ok = parse_literal(p, "null", g_null); break;
    default:
        if (c == '-' || (c >= '0' && c <= '9')) ok = parse_number(p);
        else { jp_error(p, "jsoninvalidtoken", "Invalid token found.", p->pos); ok = false; }
    }
    p->depth--;
    return ok;
}

Expr* json_parse(const char* text, size_t len, bool raw, const char* msghead) {
    JP p;
    Expr* out = NULL;
    size_t i;
    json_shared_init();
    memset(&p, 0, sizeof p);
    p.s = text;
    p.n = len;
    p.raw = raw;
    skip_ws(&p);
    if (p.pos >= p.n) {
        json_msg(msghead, "jsonnullinput", "Data in input is null.");
        return NULL;
    }
    if (parse_value(&p)) {
        skip_ws(&p);
        if (p.pos < p.n)
            jp_error(&p, "jsonexpendofinput",
                     "Unexpected character found while looking for the end of input.", p.pos);
    }
    if (!p.etag && p.sp == 1) {
        out = p.stk[0];
        p.sp = 0;
    } else {
        /* Report: the error, then where it happened. */
        size_t line = 1, col = 1, k;
        char hint[160], near[16];
        if (!p.etag) jp_error(&p, "jsoninvalidtoken", "Invalid token found.", p.pos);
        for (k = 0; k < p.epos && k < p.n; k++) {
            if (text[k] == '\n') { line++; col = 1; } else col++;
        }
        if (p.epos >= p.n) snprintf(near, sizeof near, "EOF");
        else snprintf(near, sizeof near, "%c", text[p.epos]);
        json_msg(msghead, p.etag, p.emsg);
        snprintf(hint, sizeof hint, "An error occurred near character '%s', at line %lu:%lu.",
                 near, (unsigned long)line, (unsigned long)col);
        json_msg(msghead, "jsonhintposandchar", hint);
    }
    for (i = 0; i < p.sp; i++) expr_free(p.stk[i]);
    free(p.stk);
    free(p.buf);
    return out;
}

/* ================================================================ writer */

typedef struct {
    char* s;
    size_t len, cap;
    bool raw, compact, failed;
} JW;

static void w_put(JW* w, const char* s, size_t k) {
    if (w->failed) return;
    if (w->len + k + 1 > w->cap) {
        size_t nc = w->cap ? w->cap * 2 : 1024;
        char* ns;
        while (nc < w->len + k + 1) nc *= 2;
        ns = realloc(w->s, nc);
        if (!ns) { w->failed = true; return; }
        w->s = ns;
        w->cap = nc;
    }
    memcpy(w->s + w->len, s, k);
    w->len += k;
    w->s[w->len] = '\0';
}
static void w_puts(JW* w, const char* s) { w_put(w, s, strlen(s)); }

static void w_indent(JW* w, int depth) {
    int i;
    if (w->compact) return;
    w_put(w, "\n", 1);
    for (i = 0; i < depth; i++) w_put(w, "\t", 1);
}

/* Machine real in Mathematica's JSON form (see file header). */
static void w_real(JW* w, double v) {
    char buf[64], digits[40], out[80];
    int prec, e10;
    size_t nd = 0, k;
    const char* q;
    char* o = out;
    if (!isfinite(v)) { w->failed = true; return; }
    if (v == 0.0) { w_puts(w, "0.0"); return; }
    /* Shortest round-trip digits in at most three tries: any decimal of <= 15
     * significant digits survives rounding to 15 (DBL_DIG), so when 15 digits
     * round-trip, stripping the trailing zeros below already gives the
     * shortest form; otherwise 16 or 17 digits are needed. */
    for (prec = 15; prec <= 17; prec++) {
        snprintf(buf, sizeof buf, "%.*e", prec - 1, v);
        if (prec == 17 || strtod(buf, NULL) == v) break;
    }
    /* buf = [-]d[.ddd]e[+-]XX */
    q = buf;
    if (*q == '-') { *o++ = '-'; q++; }
    while (*q && *q != 'e') { if (*q != '.') digits[nd++] = *q; q++; }
    digits[nd] = '\0';
    e10 = atoi(q + 1);
    while (nd > 1 && digits[nd - 1] == '0') digits[--nd] = '\0';
    if (e10 == -1) {
        memcpy(o, "0.", 2); o += 2;
        for (k = 0; k < nd; k++) *o++ = digits[k];
        *o = '\0';
    } else {
        *o++ = digits[0];
        *o++ = '.';
        if (nd == 1) *o++ = '0';
        else for (k = 1; k < nd; k++) *o++ = digits[k];
        *o = '\0';
        if (e10 != 0) snprintf(o, sizeof out - (size_t)(o - out), "e%d", e10);
    }
    w_puts(w, out);
}

static void w_string(JW* w, const char* s) {
    const unsigned char* u = (const unsigned char*)s;
    size_t run = 0, i;
    w_put(w, "\"", 1);
    for (i = 0; u[i]; i++) {
        const char* esc = NULL;
        char tmp[8];
        unsigned char c = u[i];
        switch (c) {
        case '"': esc = "\\\""; break;
        case '\\': esc = "\\\\"; break;
        case '/': esc = "\\/"; break;
        case '\b': esc = "\\b"; break;
        case '\f': esc = "\\f"; break;
        case '\n': esc = "\\n"; break;
        case '\r': esc = "\\r"; break;
        case '\t': esc = "\\t"; break;
        default:
            if (c < 0x20) {
                snprintf(tmp, sizeof tmp, "\\u%04x", c);
                esc = tmp;
            } else if (c == 0xC0 && u[i + 1] == 0x80) {
                esc = "\\u0000";               /* modified-UTF-8 NUL (see reader) */
            }
        }
        if (!esc) { run++; continue; }
        if (run) w_put(w, s + i - run, run);
        run = 0;
        w_puts(w, esc);
        if (c == 0xC0) i++;
    }
    if (run) w_put(w, s + i - run, run);
    w_put(w, "\"", 1);
}

static void w_fail(JW* w, const char* tag, const char* msg) {
    if (!w->failed) json_msg("Export", tag, msg);
    w->failed = true;
}

static const char* hname(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL) return NULL;
    return e->data.function.head->data.symbol.name;
}

static bool is_rule(const Expr* e) {
    const char* h = hname(e);
    return (h == SYM_Rule || h == SYM_RuleDelayed) && e->data.function.arg_count == 2;
}

static void w_value(JW* w, const Expr* e, int depth);

/* An object from `n` Rule nodes (an Association's or a rule list's). */
static void w_object(JW* w, Expr* const* rules, size_t n, int depth, bool from_assoc) {
    size_t i;
    if (n == 0) { w_puts(w, "{}"); return; }
    w_put(w, "{", 1);
    for (i = 0; i < n && !w->failed; i++) {
        const Expr* r = rules[i];
        const Expr* k = r->data.function.args[0];
        if (i) w_put(w, ",", 1);
        w_indent(w, depth + 1);
        if (k->type != EXPR_STRING) {
            if (from_assoc) {
                w_fail(w, "jsonassockeynstr", "Association contains a non-string key.");
            } else {
                char* ks = expr_to_string((Expr*)k);
                char m[256];
                snprintf(m, sizeof m, "Invalid non-string key %s.", ks ? ks : "?");
                free(ks);
                w_fail(w, "jsonrulelistkeynstr", m);
            }
            return;
        }
        w_string(w, k->data.string);
        w_put(w, ":", 1);
        w_value(w, r->data.function.args[1], depth + 1);
    }
    w_indent(w, depth);
    w_put(w, "}", 1);
}

static void w_array(JW* w, Expr* const* items, size_t n, int depth) {
    size_t i;
    if (n == 0) { w_puts(w, "[]"); return; }
    w_put(w, "[", 1);
    for (i = 0; i < n && !w->failed; i++) {
        if (i) w_put(w, ",", 1);
        w_indent(w, depth + 1);
        w_value(w, items[i], depth + 1);
    }
    w_indent(w, depth);
    w_put(w, "]", 1);
}

static void w_unencodable(JW* w, const Expr* e) {
    char m[256];
    const char* what = NULL;
    char* printed = NULL;
    if (e->type == EXPR_SYMBOL) what = e->data.symbol.name;
    else if (hname(e)) what = hname(e);
    else what = printed = expr_to_string((Expr*)e);
    snprintf(m, sizeof m, "Expression %s cannot be exported as JSON.", what ? what : "?");
    free(printed);
    w_fail(w, "jsonstrictencoding", m);
}

static void w_value(JW* w, const Expr* e, int depth) {
    char num[32];
    const char* h;
    if (w->failed) return;
    if (depth > JSON_MAX_DEPTH) { w_fail(w, "jsondepth", "Maximum nesting depth exceeded."); return; }
    switch (e->type) {
    case EXPR_INTEGER:
        snprintf(num, sizeof num, "%" PRId64, e->data.integer);
        w_puts(w, num);
        return;
    case EXPR_BIGINT: {
        char* s = mpz_get_str(NULL, 10, e->data.bigint);
        w_puts(w, s);
        free(s);
        return;
    }
    case EXPR_REAL:
        if (!isfinite(e->data.real)) { w_unencodable(w, e); return; }
        w_real(w, e->data.real);
        return;
#ifdef USE_MPFR
    case EXPR_MPFR:
        w_real(w, mpfr_get_d(e->data.mpfr, MPFR_RNDN));
        return;
#endif
    case EXPR_STRING:
        w_string(w, e->data.string);
        return;
    case EXPR_SYMBOL: {
        const char* s = e->data.symbol.name;
        if (s == SYM_True) w_puts(w, "true");
        else if (s == SYM_False) w_puts(w, "false");
        else if (s == SYM_Null) w_puts(w, "null");
        else w_unencodable(w, e);
        return;
    }
    case EXPR_NDARRAY:
        if (is_packed_list(e)) {
            Expr* l = ndarray_to_nested_list(e);
            w_value(w, l, depth);
            expr_free(l);
            return;
        }
        w_unencodable(w, e);
        return;
    default:
        break;
    }
    h = hname(e);
    if (h == SYM_Rational && e->data.function.arg_count == 2) {
        mpq_t q;
        mpz_t a, b;
        mpz_init(a); mpz_init(b);
        expr_to_mpz(e->data.function.args[0], a);
        expr_to_mpz(e->data.function.args[1], b);
        mpq_init(q);
        mpq_set_num(q, a);
        mpq_set_den(q, b);
        w_real(w, mpq_get_d(q));
        mpq_clear(q); mpz_clear(a); mpz_clear(b);
        return;
    }
    if (h == SYM_Association) {
        w_object(w, e->data.function.args, e->data.function.arg_count, depth, true);
        return;
    }
    if (h == SYM_List) {
        size_t n = e->data.function.arg_count, i, nr = 0;
        if (!w->raw) {
            for (i = 0; i < n; i++) if (is_rule(e->data.function.args[i])) nr++;
            if (nr == n && n > 0) {
                w_object(w, e->data.function.args, n, depth, false);
                return;
            }
            if (nr > 0) {
                char* ls = expr_to_string((Expr*)e);
                char m[320];
                snprintf(m, sizeof m, "%.200s contains a mixture of rule and non-rule expressions.",
                         ls ? ls : "?");
                free(ls);
                w_fail(w, "jsonrulelistnonrule", m);
                return;
            }
        }
        w_array(w, e->data.function.args, n, depth);
        return;
    }
    w_unencodable(w, e);
}

char* json_serialize(const Expr* e, bool raw, bool compact) {
    JW w;
    memset(&w, 0, sizeof w);
    w.raw = raw;
    w.compact = compact;
    w_value(&w, e, 0);
    if (w.failed) { free(w.s); return NULL; }
    if (!w.s) return mathilda_strdup("");
    return w.s;
}

/* ============================================================ builtins */

typedef enum { FMT_NONE, FMT_RAWJSON, FMT_JSON, FMT_OTHER } JFmt;

static JFmt fmt_of(const Expr* e) {
    if (!e || e->type != EXPR_STRING) return FMT_OTHER;
    if (strcmp(e->data.string, "RawJSON") == 0) return FMT_RAWJSON;
    if (strcmp(e->data.string, "JSON") == 0) return FMT_JSON;
    return FMT_OTHER;
}

/* Scan trailing options from args[from..): "Compact" -> True/False (a
 * string or symbol name).  Returns false on a non-option argument. */
static bool read_opts(const Expr* res, size_t from, bool* compact) {
    size_t i;
    for (i = from; i < res->data.function.arg_count; i++) {
        const Expr* o = res->data.function.args[i];
        const Expr* k;
        if (!is_rule(o)) return false;
        k = o->data.function.args[0];
        if ((k->type == EXPR_STRING && strcmp(k->data.string, "Compact") == 0) ||
            (k->type == EXPR_SYMBOL && strcmp(k->data.symbol.name, "Compact") == 0)) {
            const Expr* v = o->data.function.args[1];
            *compact = v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True;
        }
        /* Other options are accepted and ignored. */
    }
    return true;
}

/* ImportString["str", "RawJSON" | "JSON"] ; ImportString["str"] is the text. */
static Expr* builtin_importstring(Expr* res) {
    const Expr* s;
    JFmt f;
    bool compact = false;
    Expr* out;
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    s = res->data.function.args[0];
    if (s->type != EXPR_STRING) return NULL;
    if (res->data.function.arg_count == 1) return expr_copy((Expr*)s);
    f = fmt_of(res->data.function.args[1]);
    if (f != FMT_RAWJSON && f != FMT_JSON) return NULL;
    if (!read_opts(res, 2, &compact)) return NULL;
    out = json_parse(s->data.string, strlen(s->data.string), f == FMT_RAWJSON, "Import");
    return out ? out : failed();
}

/* ExportString[expr, "RawJSON" | "JSON", "Compact" -> True|False] */
static Expr* builtin_exportstring(Expr* res) {
    JFmt f;
    bool compact = false;
    char* txt;
    Expr* out;
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    f = fmt_of(res->data.function.args[1]);
    if (f != FMT_RAWJSON && f != FMT_JSON) return NULL;
    if (!read_opts(res, 2, &compact)) return NULL;
    txt = json_serialize(res->data.function.args[0], f == FMT_RAWJSON, compact);
    if (!txt) return failed();
    out = expr_new_string(txt);
    free(txt);
    return out;
}

/* ------------------------------------------------------ Import / Export */

static BuiltinFunc g_prev_import, g_prev_export;

static bool ends_with_json(const char* path) {
    size_t n = strlen(path);
    return n >= 5 && path[n - 5] == '.' &&
           tolower((unsigned char)path[n - 4]) == 'j' &&
           tolower((unsigned char)path[n - 3]) == 's' &&
           tolower((unsigned char)path[n - 2]) == 'o' &&
           tolower((unsigned char)path[n - 1]) == 'n';
}

/* Import["f.json"] (JSON rule form), Import[f, "JSON" | "RawJSON"]. */
static Expr* json_import(Expr* res) {
    const Expr* p;
    JFmt f = FMT_NONE;
    bool compact = false;
    FILE* fp;
    char* data;
    long sz;
    size_t got;
    Expr* out;
    if (res->type == EXPR_FUNCTION && res->data.function.arg_count >= 1 &&
        res->data.function.args[0]->type == EXPR_STRING) {
        p = res->data.function.args[0];
        if (res->data.function.arg_count >= 2) f = fmt_of(res->data.function.args[1]);
        else if (ends_with_json(p->data.string)) f = FMT_JSON;
        if ((f == FMT_JSON || f == FMT_RAWJSON) && read_opts(res, 2, &compact)) {
            fp = fopen(p->data.string, "rb");
            if (!fp) {
                json_msg("Import", "nffil", "File not found during Import.");
                return failed();
            }
            if (fseek(fp, 0, SEEK_END) != 0 || (sz = ftell(fp)) < 0 ||
                fseek(fp, 0, SEEK_SET) != 0) {
                fclose(fp);
                return failed();
            }
            data = malloc((size_t)sz + 1);
            if (!data) { fclose(fp); return failed(); }
            got = fread(data, 1, (size_t)sz, fp);
            fclose(fp);
            data[got] = '\0';
            out = json_parse(data, got, f == FMT_RAWJSON, "Import");
            free(data);
            return out ? out : failed();
        }
    }
    return g_prev_import ? g_prev_import(res) : NULL;
}

/* Export["f.json", expr] (JSON form: associations and rule lists are
 * objects), Export[f, expr, "JSON" | "RawJSON", opts]. */
static Expr* json_export(Expr* res) {
    JFmt f = FMT_NONE;
    bool compact = false;
    size_t optfrom = 2;
    if (res->type == EXPR_FUNCTION && res->data.function.arg_count >= 2 &&
        res->data.function.args[0]->type == EXPR_STRING) {
        const char* path = res->data.function.args[0]->data.string;
        if (res->data.function.arg_count >= 3 &&
            res->data.function.args[2]->type == EXPR_STRING) {
            f = fmt_of(res->data.function.args[2]);
            optfrom = 3;
        } else if (ends_with_json(path)) {
            f = FMT_JSON;
        }
        if ((f == FMT_JSON || f == FMT_RAWJSON) && read_opts(res, optfrom, &compact)) {
            char* txt = json_serialize(res->data.function.args[1], f == FMT_RAWJSON, compact);
            FILE* fp;
            size_t n;
            if (!txt) return failed();
            fp = fopen(path, "wb");
            if (!fp) { free(txt); return failed(); }
            n = strlen(txt);
            if (fwrite(txt, 1, n, fp) != n) { fclose(fp); free(txt); return failed(); }
            fclose(fp);
            free(txt);
            return expr_new_string(path);
        }
    }
    return g_prev_export ? g_prev_export(res) : NULL;
}

/* Extend an existing docstring (Import / Export are registered by imageio.c). */
static void append_doc(const char* name, const char* extra) {
    const char* old = symtab_get_docstring(name);
    size_t n = (old ? strlen(old) : 0) + strlen(extra) + 1;
    char* s = malloc(n);
    if (!s) return;
    snprintf(s, n, "%s%s", old ? old : "", extra);
    symtab_set_docstring(name, s);
    free(s);
}

void json_init(void) {
    SymbolDef* def;
    json_shared_init();

    symtab_add_builtin("ImportString", builtin_importstring);
    symtab_get_def("ImportString")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ImportString",
        "ImportString[\"data\", \"RawJSON\"] parses JSON text: objects become associations, "
        "arrays lists, strings strings, true/false/null True/False/Null; integers of any size "
        "are Integers, numbers with a fraction are Reals, and numbers with an exponent but no "
        "fraction (1e-2) are the exact Integer or Rational they denote. Escapes, including "
        "\\uXXXX and UTF-16 surrogate pairs, are decoded to UTF-8. ImportString[\"data\", "
        "\"JSON\"] gives objects as lists of rules instead. Invalid JSON gives $Failed with an "
        "Import::json... message naming the error and its line:column. ImportString[\"data\"] "
        "gives the text itself.");

    symtab_add_builtin("ExportString", builtin_exportstring);
    symtab_get_def("ExportString")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ExportString",
        "ExportString[expr, \"RawJSON\"] writes expr as JSON text: associations as objects "
        "(keys must be strings), lists as arrays, True/False/Null as true/false/null, integers "
        "exactly, reals and rationals as shortest round-trip reals. The layout matches "
        "Mathematica's (tab indentation, \"key\":value); \"Compact\" -> True drops all "
        "whitespace. ExportString[expr, \"JSON\"] also writes lists of rules as objects. An "
        "expression with no JSON form gives $Failed with an Export::json... message.");

    /* Import / Export: claim the JSON formats, delegate everything else to
     * the image/graphics implementation registered by imageio_init. */
    def = symtab_get_def("Import");
    if (def->builtin_func != json_import) {
        g_prev_import = def->builtin_func;
        def->builtin_func = json_import;
        append_doc("Import",
            " Import[\"file.json\"] reads JSON with objects as lists of rules (the \"JSON\" "
            "format); Import[\"file\", \"RawJSON\"] gives objects as associations.");
    }
    def = symtab_get_def("Export");
    if (def->builtin_func != json_export) {
        g_prev_export = def->builtin_func;
        def->builtin_func = json_export;
        append_doc("Export",
            " Export[\"file.json\", expr] writes JSON (associations and lists of rules as "
            "objects); Export[\"file\", expr, \"RawJSON\"] accepts only associations as objects; "
            "\"Compact\" -> True drops whitespace.");
    }
}
