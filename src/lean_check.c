/*
 * lean_check.c -- the structural source check of LeanVerify and Certify.
 *
 *   Lean`Internal`CheckSource[src, mathlib, certify, prelude | None]
 *       -> <|"OK" -> True, "Source" -> s, "Theorems" -> {names...}|>
 *        | <|"OK" -> False, "Kind" -> k, "What" -> w, "Line" -> n|>
 *
 * Why it exists.  A Lean source can run arbitrary code while it is being
 * elaborated (a tactic elaborator attached with `attribute [tactic ...]`, a
 * macro, `#eval`, `initialize`, ...).  Such code can forge every output of the
 * elaborating process, including the compiled module that the kernel replay
 * then reads.  So LeanVerify accepts only sources that cannot contain
 * user-defined code at all, and it decides that WITHOUT reproducing Lean's
 * lexer: the accepted lexical surface is so small that there is nothing for
 * two lexers to disagree about.
 *
 * Rule A, the minimal lexical surface (kind "lex").  Refused anywhere:
 *   - the double quote (no strings, raw strings, interpolated strings,
 *     s! / m! / f! literals);
 *   - "/-", "-/" and "--" (no comments or doc comments of any kind);
 *   - French quotes, the backquote, '#', '$', '%', '\', "@[";
 *   - tabs, carriage returns, NUL and every other control character;
 *   - a quote ' that is not an identifier prime: it must follow an ASCII
 *     letter, digit, '_' or ', and the run of [A-Za-z0-9_'] it ends must
 *     begin with a letter or '_' (so no character literal can occur);
 *   - invalid or overlong UTF-8, and every non-ASCII codepoint outside a
 *     fixed allowlist of mathematical symbols and Greek letters (below).
 * With no string, comment or character literal there is no span whose
 * content Lean skips: every byte belongs to a token the checker sees.
 *
 * Rule B, finer-than-Lean words (kind "token").  A word is a maximal run of
 * [A-Za-z0-9_]; every other byte is a separator.  A Lean identifier never
 * starts inside such a run (it absorbs [A-Za-z0-9_'!?] and letter-like
 * Unicode), except after a numeric literal: in a run that starts with a
 * digit, every suffix that starts with a letter or '_' is a word too.  So
 * every keyword Lean can see is a word here (at worst the checker refuses
 * more).  Refused anywhere: every command-leading word of Lean's parser
 * (lean_check_tokens.h, derived from Lean's parser tables and tested against
 * them) except an allowed command starter in its starting position, and the
 * escape words below (by_elab, native_decide, sorry, eval, IO, ...).
 *
 * Command grammar (kinds "cmd", "import").  Commands start at column 0;
 * every continuation line is indented.  Each command starts with
 *
 *     import Mathlib | import Mathlib.<Path>    only with `mathlib`, only first
 *     theorem NAME ... | lemma NAME ...          (lemma only with `mathlib`)
 *     example ...                                (renamed, see below)
 *     open A B ...                               namespace names only (no
 *                                                `in`, no `scoped`)
 *     namespace A.B | section [A] | noncomputable section [A] | end [A]
 *     variable ... | universe u ...
 *     set_option maxHeartbeats|maxRecDepth|synthInstance.maxHeartbeats <n> [in]
 *     set_option autoImplicit|relaxedAutoImplicit|linter.<x> true|false [in]
 *     def mathilda_cert_<k>_<x> ...              only with `certify` (generator
 *                                                auxiliaries of Lean.m)
 *
 * and nothing else (no #print, no def / instance / attribute / macro ...).
 *
 * Trusted text.  `prelude` (the text of src/internal/lean/Prelude.lean) may
 * occur once, byte-identical, at the start of a line; it is skipped (its
 * definitions are Mathilda's own, reviewed code).  It must itself contain no
 * string, comment, character literal or French quote, and end with a newline,
 * so Lean's lexer is at the top level where it ends.  With `certify`, the
 * generator's footer lines `#print axioms mathilda_cert_<k>_footer` (a whole
 * line, at column 0) are skipped too.
 *
 * Output.  "Source" is src with every `example` keyword replaced by
 * `noncomputable def _root_.mathilda_example_<k>` (same line, so the line
 * numbers do not move) and "Theorems" the fully qualified names of every
 * theorem / lemma / renamed example, with the namespace / section stack
 * tracked through namespace / section / end.
 */
#include "lean_run.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "assoc.h"
#include "attr.h"
#include "lean_check_tokens.h"
#include "sym_names.h"
#include "symtab.h"

typedef struct {
    const char* kind;    /* "lex", "token", "cmd", "import" */
    char* what;          /* owned */
    int line;
} CheckErr;

static char* copy_n(const char* s, size_t n) {
    char* r = malloc(n + 1);
    if (!r) abort();
    memcpy(r, s, n);
    r[n] = '\0';
    return r;
}

/* 1-based line of byte offset pos. */
static int line_of(const char* s, size_t pos) {
    int line = 1;
    for (size_t i = 0; i < pos; i++)
        if (s[i] == '\n') line++;
    return line;
}

static void set_err(CheckErr* e, const char* kind, int line, const char* what, size_t n) {
    if (e->kind) return;            /* the first problem is reported */
    e->kind = kind;
    e->line = line;
    e->what = copy_n(what, n);
}

static void set_err_s(CheckErr* e, const char* kind, int line, const char* what) {
    set_err(e, kind, line, what, strlen(what));
}

static bool is_alpha(unsigned char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_';
}
static bool is_digit(unsigned char c) { return c >= '0' && c <= '9'; }
static bool is_word(unsigned char c) { return is_alpha(c) || is_digit(c); }

/* ------------------------------------------------------------------ */
/* Word lists                                                          */
/* ------------------------------------------------------------------ */

/* Words that escape to code or bypass the kernel inside an accepted command
 * (tactics, term elaborators and library names; the command keywords are in
 * lean_check_tokens.h).  Hand-maintained: this is the part of the check that
 * is a deny-list, and the documented remaining assumption is that built-in
 * elaborators reachable through ordinary syntax run no user code. */
static const char* const ESCAPE_WORDS[] = {
    "by_elab", "run_tac", "run_elab", "run_meta", "run_cmd", "native_decide", "native",
    "implemented_by", "extern", "unsafe", "sorry", "admit", "ofReduceBool", "ofReduceNat",
    "reduceBool", "trustCompiler", "attribute", "debug", "elab", "macro", "eval",
    "skipKernelTC", "addDecl", "addDeclCore", "addAndCompile", "addDeclWithoutChecking",
    "setEnv", "modifyEnv", "withOptions", "evalExpr", "evalConst", "unsafeBaseIO",
    "unsafeIO", "unsafeEIO", "unsafePerformIO", "unsafeCast", "IO", "EIO", "BaseIO", "Elab",
    "MetaM", "TermElabM", "TacticM", "CommandElabM", "CoreM", "MacroM", "command_elab",
    "term_elab", "tactic_elab", "plausible", "slim_check", "polyrith", "bv_decide",
    "lemma", "module", "export",
    NULL
};

static bool in_list(const char* const* l, const char* s, size_t n) {
    for (size_t k = 0; l[k]; k++)
        if (strlen(l[k]) == n && memcmp(l[k], s, n) == 0) return true;
    return false;
}

static bool refused_word(const char* s, size_t n) {
    return in_list(ESCAPE_WORDS, s, n) || in_list(LEAN_COMMAND_WORDS_CORE, s, n) ||
           in_list(LEAN_COMMAND_WORDS_MATHLIB, s, n);
}

/* ------------------------------------------------------------------ */
/* Rule A: the lexical surface                                         */
/* ------------------------------------------------------------------ */

/* Decode one UTF-8 sequence at s[i..n) strictly (no overlong forms, no
 * surrogates, at most U+10FFFF).  Returns its length, 0 when invalid. */
static size_t utf8_decode(const unsigned char* s, size_t i, size_t n, unsigned long* cp) {
    unsigned char c = s[i];
    size_t len;
    unsigned long v, min;
    if (c < 0x80) { *cp = c; return 1; }
    if (c >= 0xC2 && c <= 0xDF) { len = 2; v = c & 0x1F; min = 0x80; }
    else if (c >= 0xE0 && c <= 0xEF) { len = 3; v = c & 0x0F; min = 0x800; }
    else if (c >= 0xF0 && c <= 0xF4) { len = 4; v = c & 0x07; min = 0x10000; }
    else return 0;
    if (i + len > n) return 0;
    for (size_t k = 1; k < len; k++) {
        if ((s[i + k] & 0xC0) != 0x80) return 0;
        v = (v << 6) | (s[i + k] & 0x3F);
    }
    if (v < min || v > 0x10FFFF || (v >= 0xD800 && v <= 0xDFFF)) return 0;
    *cp = v;
    return len;
}

/* The non-ASCII codepoints a LeanVerify source may contain: the symbols of
 * ordinary theorem statements and proofs, and the Greek letters. */
static bool allowed_codepoint(unsigned long cp) {
    static const unsigned long ok[] = {
        0x2200 /* forall */, 0x2203 /* exists */, 0x2192 /* -> */, 0x2190 /* <- */,
        0x2194 /* <-> */, 0x2227 /* and */, 0x2228 /* or */, 0x00AC /* not */,
        0x2264 /* <= */, 0x2265 /* >= */, 0x2260 /* != */, 0x27E8, 0x27E9 /* angle brackets */,
        0x21A6 /* mapsto */, 0x00B7 /* middle dot */, 0x2218 /* compose */, 0x2208 /* in */,
        0x2209 /* notin */, 0x2223 /* divides */, 0x00D7 /* times */, 0x207B /* sup minus */,
        0x00B9, 0x00B2, 0x00B3 /* sup 1 2 3 */, 0x2115, 0x2124, 0x211A, 0x211D, 0x2102 /* NZQRC */,
        0x22A2 /* turnstile */, 0x2261 /* equiv */, 0x2211 /* sum */, 0x220F /* prod */,
        0x221E /* infinity */, 0x230A, 0x230B, 0x2308, 0x2309 /* floor, ceiling */,
        0x2016 /* norm */, 0
    };
    for (size_t k = 0; ok[k]; k++)
        if (ok[k] == cp) return true;
    return (cp >= 0x0391 && cp <= 0x03A9 && cp != 0x03A2) || (cp >= 0x03B1 && cp <= 0x03C9);
}

/* Is the quote at s[i] an identifier prime? */
static bool prime_ok(const char* s, size_t a, size_t i) {
    if (i == a) return false;
    unsigned char p = (unsigned char)s[i - 1];
    if (!(is_word(p) || p == '\'')) return false;
    size_t j = i;
    while (j > a && (is_word((unsigned char)s[j - 1]) || s[j - 1] == '\'')) j--;
    return is_alpha((unsigned char)s[j]);
}

/* Rule A on s[a, b).  false (with err) on the first refused construct. */
static bool surface_ok(const char* s, size_t a, size_t b, CheckErr* err) {
    const unsigned char* u = (const unsigned char*)s;
    for (size_t i = a; i < b; ) {
        unsigned char c = u[i];
        unsigned char nx = i + 1 < b ? u[i + 1] : 0;
        const char* what = NULL;
        size_t wlen = 0;
        if (c == '\n' || c == ' ') { i++; continue; }
        if (c < 0x20 || c == 0x7F)
            what = c == '\t' ? "a tab" : c == '\r' ? "a carriage return" : "a control character";
        else if (c == '"') what = "a string literal (\")";
        else if (c == '-' && nx == '-') what = "a comment (--)";
        else if (c == '/' && nx == '-') what = "a comment (/-)";
        else if (c == '-' && nx == '/') what = "a comment end (-/)";
        else if (c == '`') what = "a backquote (name literal or syntax quotation)";
        else if (c == '#') {
            size_t j = i + 1;
            while (j < b && is_word(u[j])) j++;
            what = s + i;
            wlen = j > i + 1 ? j - i : 1;
        } else if (c == '%') {
            size_t j = i;
            while (j > a && is_word(u[j - 1])) j--;
            what = j < i ? s + j : "%";
            wlen = j < i ? i + 1 - j : 1;
        } else if (c == '$') what = "$";
        else if (c == '\\') what = "a backslash (escape or \\fun)";
        else if (c == '@' && nx == '[') what = "@[ (an attribute)";
        else if (c == '\'' && !prime_ok(s, a, i))
            what = "a character literal or a quote that is not an identifier prime";
        else if (c >= 0x80) {
            unsigned long cp = 0;
            size_t len = utf8_decode(u, i, b, &cp);
            if (len == 0) what = "invalid UTF-8";
            else if (cp == 0xAB || cp == 0xBB) what = "a French quote (escaped identifier)";
            else if (!allowed_codepoint(cp)) {
                static char buf[48];
                snprintf(buf, sizeof buf, "the character U+%04lX", cp);
                what = buf;
            } else {
                i += len;
                continue;
            }
        }
        if (what) {
            set_err(err, "lex", line_of(s, i), what, wlen ? wlen : strlen(what));
            return false;
        }
        i++;
    }
    return true;
}

/* The prelude must not open any span Lean would skip, and must end at a line
 * end: then Lean's lexer is at the top level where the prelude ends. */
static bool prelude_span_free(const char* p) {
    size_t n = strlen(p);
    if (n == 0 || p[n - 1] != '\n') return false;
    for (size_t i = 0; i < n; i++) {
        unsigned char c = (unsigned char)p[i], nx = (unsigned char)p[i + 1];
        if (c == '"' || c == '`' || c == '\'' || c == '\\' || c == '\t' || c == '\r') return false;
        if ((c == '-' && (nx == '-' || nx == '/')) || (c == '/' && nx == '-')) return false;
        if (c == 0xC2 && (nx == 0xAB || nx == 0xBB)) return false;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/* Rule B: words                                                       */
/* ------------------------------------------------------------------ */

/* Every word of s[a, b) (and, in a digit-led run, every suffix starting with a
 * letter or '_') must not be a refused word.  `skip` words starting at the
 * given offsets (command starters) are exempt. */
static bool words_ok(const char* s, size_t a, size_t b, const size_t* skip, size_t nskip,
                     CheckErr* err) {
    const unsigned char* u = (const unsigned char*)s;
    for (size_t i = a; i < b; ) {
        if (!is_word(u[i])) { i++; continue; }
        size_t j = i;
        while (j < b && is_word(u[j])) j++;
        bool exempt = false;
        for (size_t k = 0; k < nskip; k++)
            if (skip[k] == i) exempt = true;
        if (!exempt && refused_word(s + i, j - i)) {
            set_err(err, "token", line_of(s, i), s + i, j - i);
            return false;
        }
        if (is_digit(u[i]))
            for (size_t k = i + 1; k < j; k++)
                if (is_alpha(u[k]) && refused_word(s + k, j - k)) {
                    set_err(err, "token", line_of(s, k), s + k, j - k);
                    return false;
                }
        /* A glued token (Mathlib's `->o`, `[MOD`, ...) ending inside this word:
         * Lean's next token starts where it ends. */
        for (size_t g = 0; LEAN_GLUED_TOKENS_MATHLIB[g]; g++) {
            const char* t = LEAN_GLUED_TOKENS_MATHLIB[g];
            size_t tl = strlen(t), tail = 0;
            while (tail < tl && is_word((unsigned char)t[tl - 1 - tail])) tail++;
            size_t pre = tl - tail, k = i + tail;
            if (tail == 0 || i < a + pre || k >= j || memcmp(s + i - pre, t, tl) != 0) continue;
            size_t e = j;
            if (refused_word(s + k, e - k) ||
                (is_digit(u[k]) && !words_ok(s, k, e, NULL, 0, err))) {
                set_err(err, "token", line_of(s, k), s + k, e - k);
                return false;
            }
        }
        i = j;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/* Command grammar (on bytes; Rule A already holds)                    */
/* ------------------------------------------------------------------ */

static size_t skip_ws(const char* s, size_t i, size_t b) {
    while (i < b && (s[i] == ' ' || s[i] == '\n')) i++;
    return i;
}

static size_t word_end(const char* s, size_t i, size_t b) {
    while (i < b && is_word((unsigned char)s[i])) i++;
    return i;
}

/* A dotted name [A-Za-z_][A-Za-z0-9_']*(.[A-Za-z_][A-Za-z0-9_']*)* at i;
 * returns its end (i when there is none). */
static size_t name_end(const char* s, size_t i, size_t b, bool primes) {
    size_t j = i;
    for (;;) {
        if (j >= b || !is_alpha((unsigned char)s[j])) return j == i ? i : j - 1;
        j++;
        while (j < b && (is_word((unsigned char)s[j]) || (primes && s[j] == '\''))) j++;
        if (j + 1 < b && s[j] == '.' && is_alpha((unsigned char)s[j + 1])) { j++; continue; }
        return j;
    }
}

/* May a name end at e (so that Lean's identifier ends there too)? */
static bool name_stops(const char* s, size_t e, size_t b, const char* also) {
    if (e >= b || s[e] == ' ' || s[e] == '\n') return true;
    return strchr(also, s[e]) != NULL && s[e] != '\0';
}

/* Only spaces / newlines in s[i, b)? */
static bool blank(const char* s, size_t i, size_t b) { return skip_ws(s, i, b) == b; }

/* The namespace / section stack: entries are dotted paths, "" for a section. */
typedef struct { char** v; size_t n, cap; } StrStack;

static void ss_push(StrStack* st, char* p) {
    if (st->n == st->cap) {
        size_t nc = st->cap ? st->cap * 2 : 8;
        char** q = realloc(st->v, nc * sizeof *q);
        if (!q) abort();
        st->v = q;
        st->cap = nc;
    }
    st->v[st->n++] = p;
}

static char* qualified(const StrStack* st, const char* name) {
    if (strncmp(name, "_root_.", 7) == 0) return copy_n(name + 7, strlen(name) - 7);
    size_t len = strlen(name) + 1;
    for (size_t k = 0; k < st->n; k++) len += strlen(st->v[k]) + 1;
    char* r = malloc(len);
    if (!r) abort();
    r[0] = '\0';
    for (size_t k = 0; k < st->n; k++)
        if (st->v[k][0]) { strcat(r, st->v[k]); strcat(r, "."); }
    strcat(r, name);
    return r;
}

static bool is_text(const char* s, size_t a, size_t e, const char* w) {
    size_t n = strlen(w);
    return e - a == n && memcmp(s + a, w, n) == 0;
}

/* set_option NAME VALUE [in], starting after the keyword at i. */
static bool option_ok(const char* s, size_t i, size_t b) {
    size_t a = skip_ws(s, i, b), e = name_end(s, a, b, false);
    if (e == a || !name_stops(s, e, b, "")) return false;
    char* nm = copy_n(s + a, e - a);
    bool numeric = strcmp(nm, "maxHeartbeats") == 0 || strcmp(nm, "maxRecDepth") == 0 ||
                   strcmp(nm, "synthInstance.maxHeartbeats") == 0;
    bool boolean = strcmp(nm, "autoImplicit") == 0 || strcmp(nm, "relaxedAutoImplicit") == 0 ||
                   (strncmp(nm, "linter.", 7) == 0 && strlen(nm) > 7);
    free(nm);
    size_t v = skip_ws(s, e, b), ve = word_end(s, v, b);
    bool ok = false;
    if (numeric && ve > v) {
        ok = true;
        for (size_t k = v; k < ve; k++)
            if (!is_digit((unsigned char)s[k])) ok = false;
    }
    if (boolean) ok = is_text(s, v, ve, "true") || is_text(s, v, ve, "false");
    if (!ok || !name_stops(s, ve, b, "")) return false;
    size_t r = skip_ws(s, ve, b);
    size_t re = word_end(s, r, b);
    if (is_text(s, r, re, "in")) r = re;
    return blank(s, r, b);
}

/* Is `nm` a generator auxiliary name mathilda_cert_<digits>_<word>? */
static bool cert_aux_name(const char* nm) {
    const char* p = "mathilda_cert_";
    if (strncmp(nm, p, strlen(p)) != 0) return false;
    nm += strlen(p);
    if (!is_digit((unsigned char)*nm)) return false;
    while (is_digit((unsigned char)*nm)) nm++;
    if (*nm != '_' || !nm[1]) return false;
    for (nm++; *nm; nm++)
        if (!is_word((unsigned char)*nm)) return false;
    return true;
}

typedef struct {
    bool mathlib, certify, past_imports;
    StrStack st;
    StrStack thms;
    size_t* examples;
    int nexamples;
} CmdState;

static bool cmd_err(CheckErr* err, const char* s, size_t at, const char* what) {
    set_err_s(err, "cmd", line_of(s, at), what);
    return false;
}

/* Check one command s[a, b) (a at column 0). */
static bool check_command(const char* s, size_t a, size_t b, CmdState* cs, CheckErr* err) {
    size_t we = word_end(s, a, b);
    size_t skip[2] = { a, (size_t)-1 };
    if (we == a || !is_alpha((unsigned char)s[a])) {
        unsigned long cp = 0;
        size_t len = utf8_decode((const unsigned char*)s, a, b, &cp);
        set_err(err, "cmd", line_of(s, a), s + a, len ? len : 1);
        return false;
    }
    /* A starter is followed by a space, a newline or the end: `variable?`,
     * `theorem'` are other tokens. */
    if (we < b && s[we] != ' ' && s[we] != '\n') {
        size_t e = we;
        while (e < b && s[e] != ' ' && s[e] != '\n') e++;
        set_err(err, "cmd", line_of(s, a), s + a, e - a);
        return false;
    }
    if (is_text(s, a, we, "import")) {
        size_t m = skip_ws(s, we, b), me = name_end(s, m, b, false);
        bool ok = cs->mathlib && !cs->past_imports && me > m && blank(s, me, b) &&
                  (is_text(s, m, me, "Mathlib") ||
                   (me - m > 8 && memcmp(s + m, "Mathlib.", 8) == 0));
        if (!ok) {
            size_t e = a;
            while (e < b && s[e] != '\n') e++;
            while (e > a && s[e - 1] == ' ') e--;
            set_err(err, "import", line_of(s, a), s + a, e - a);
        }
        return ok;    /* module names are not scanned for words */
    }
    cs->past_imports = true;
    bool ok = true;
    if (is_text(s, a, we, "theorem") || (cs->mathlib && is_text(s, a, we, "lemma"))) {
        size_t m = skip_ws(s, we, b), me = name_end(s, m, b, true);
        if (me == m || !name_stops(s, me, b, ":({["))
            return cmd_err(err, s, a, "theorem without a plain name (ASCII letters, digits, _, ', .)");
        char* nm = copy_n(s + m, me - m);
        ss_push(&cs->thms, qualified(&cs->st, nm));
        free(nm);
    } else if (is_text(s, a, we, "example")) {
        char nm[64];
        cs->nexamples++;
        snprintf(nm, sizeof nm, "mathilda_example_%d", cs->nexamples);
        ss_push(&cs->thms, copy_n(nm, strlen(nm)));
        cs->examples[cs->nexamples - 1] = a;
    } else if (cs->certify && is_text(s, a, we, "def")) {
        size_t m = skip_ws(s, we, b), me = name_end(s, m, b, false);
        char* nm = copy_n(s + m, me - m);
        ok = name_stops(s, me, b, ":({[") && cert_aux_name(nm);
        free(nm);
        if (!ok) return cmd_err(err, s, a, "def");
    } else if (is_text(s, a, we, "namespace")) {
        size_t m = skip_ws(s, we, b), me = name_end(s, m, b, false);
        if (me == m || !blank(s, me, b)) return cmd_err(err, s, a, "namespace");
        ss_push(&cs->st, copy_n(s + m, me - m));
    } else if (is_text(s, a, we, "section") || is_text(s, a, we, "noncomputable")) {
        size_t m = skip_ws(s, we, b);
        if (is_text(s, a, we, "noncomputable")) {
            size_t se = word_end(s, m, b);
            if (!is_text(s, m, se, "section") || !name_stops(s, se, b, ""))
                return cmd_err(err, s, a, "noncomputable (only noncomputable section)");
            skip[1] = m;
            m = skip_ws(s, se, b);
        }
        size_t me = name_end(s, m, b, false);
        if (!blank(s, me, b)) return cmd_err(err, s, a, "section");
        ss_push(&cs->st, copy_n("", 0));
    } else if (is_text(s, a, we, "end")) {
        size_t m = skip_ws(s, we, b), me = name_end(s, m, b, false);
        if (!blank(s, me, b)) return cmd_err(err, s, a, "end");
        if (cs->st.n > 0) free(cs->st.v[--cs->st.n]);
    } else if (is_text(s, a, we, "variable") || is_text(s, a, we, "universe")) {
        /* the body is scanned below */
    } else if (is_text(s, a, we, "open")) {
        size_t m = skip_ws(s, we, b);
        if (m == b) return cmd_err(err, s, a, "open (a namespace name must follow)");
        while (m < b) {
            size_t me = name_end(s, m, b, false);
            if (me == m || !name_stops(s, me, b, "") || is_text(s, m, me, "in") ||
                is_text(s, m, me, "scoped") || is_text(s, m, me, "hiding") ||
                is_text(s, m, me, "renaming"))
                return cmd_err(err, s, a, "open (only namespace names may follow: no in, "
                                          "scoped, hiding, renaming or (...))");
            m = skip_ws(s, me, b);
        }
    } else if (is_text(s, a, we, "set_option")) {
        if (!option_ok(s, we, b)) {
            size_t e = a;
            while (e < b && s[e] != '\n') e++;
            set_err(err, "cmd", line_of(s, a), s + a, e - a);
            return false;
        }
    } else {
        set_err(err, "cmd", line_of(s, a), s + a, we - a);
        return false;
    }
    return words_ok(s, a, b, skip, skip[1] == (size_t)-1 ? 1 : 2, err);
}

/* Is s[i, e) the generator footer line `#print axioms mathilda_cert_<k>_footer`? */
static bool footer_line(const char* s, size_t i, size_t e) {
    const char* p = "#print axioms mathilda_cert_";
    size_t n = strlen(p);
    if (e - i <= n || memcmp(s + i, p, n) != 0) return false;
    size_t k = i + n;
    if (k >= e || !is_digit((unsigned char)s[k])) return false;
    while (k < e && is_digit((unsigned char)s[k])) k++;
    return is_text(s, k, e, "_footer");
}

/* Check the segment s[a, b) (a at a line start; no trusted text inside). */
static bool check_segment(const char* s, size_t a, size_t b, CmdState* cs, CheckErr* err) {
    if (!surface_ok(s, a, b, err)) return false;
    size_t i = a;
    /* Content before the first column-0 line would continue the previous
     * command (or be an indented first command). */
    while (i < b) {
        size_t e = i;
        while (e < b && s[e] != '\n') e++;
        if (i < e && s[i] != ' ') break;
        if (!blank(s, i, e)) {
            set_err_s(err, "cmd", line_of(s, i), "an indented first command (commands start "
                                                 "at column 0)");
            return false;
        }
        i = e < b ? e + 1 : e;
    }
    while (i < b) {
        size_t e = i;
        /* The command runs to the next line that starts with a non-space. */
        for (;;) {
            while (e < b && s[e] != '\n') e++;
            if (e >= b) break;
            e++;
            if (e < b && s[e] != ' ' && s[e] != '\n') break;
        }
        if (!check_command(s, i, e, cs, err)) return false;
        i = e;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/* The builtin                                                         */
/* ------------------------------------------------------------------ */

static Expr* rule(const char* key, Expr* v) {
    Expr* a[2] = { expr_new_string(key), v };
    return expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
}

static Expr* assoc_of(Expr** rules, size_t n) {
    Expr* r = assoc_from_rules(rules, n);
    for (size_t i = 0; i < n; i++) expr_free(rules[i]);
    return r;
}

static bool is_symbol(const Expr* e, const char* name) {
    return e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

/* Check the whole source: split it at the trusted text (the prelude once,
 * Certify's footer lines) and check every segment in between. */
static bool check_source(const char* s, size_t n, const char* prelude, CmdState* cs,
                         CheckErr* err) {
    size_t plen = prelude ? strlen(prelude) : 0;
    bool pre_seen = false;
    if (plen > 0 && !prelude_span_free(prelude)) {
        set_err_s(err, "lex", 1, "the trusted prelude (it must contain no string, comment, "
                                 "character literal or French quote and end with a newline)");
        return false;
    }
    /* Import lines first, so that the reason reported for a source that
     * imports a meta-level library is the import (not whatever else it
     * contains).  The same rule is enforced again by the command grammar. */
    for (size_t i = 0; i < n; ) {
        size_t e = i;
        while (e < n && s[e] != '\n') e++;
        if (e - i >= 7 && memcmp(s + i, "import ", 7) == 0) {
            size_t m = skip_ws(s, i + 6, e), me = name_end(s, m, e, false);
            bool ok = cs->mathlib && me > m && blank(s, me, e) &&
                      (is_text(s, m, me, "Mathlib") ||
                       (me - m > 8 && memcmp(s + m, "Mathlib.", 8) == 0));
            if (!ok) {
                size_t len = e - i;
                while (len > 0 && s[i + len - 1] == ' ') len--;
                set_err(err, "import", line_of(s, i), s + i, len);
                return false;
            }
        }
        i = e < n ? e + 1 : e;
    }
    size_t seg = 0;
    for (size_t i = 0; i < n; ) {
        /* i is a line start */
        size_t e = i;
        while (e < n && s[e] != '\n') e++;
        size_t skip_to = 0;
        if (plen > 0 && !pre_seen && n - i >= plen && memcmp(s + i, prelude, plen) == 0) {
            pre_seen = true;
            skip_to = i + plen;
        } else if (cs->certify && footer_line(s, i, e)) {
            skip_to = e < n ? e + 1 : e;
        }
        if (skip_to) {
            if (!check_segment(s, seg, i, cs, err)) return false;
            seg = i = skip_to;
            continue;
        }
        i = e < n ? e + 1 : e;
    }
    return check_segment(s, seg, n, cs, err);
}

static Expr* builtin_lean_checksource(Expr* res) {
    if (res->data.function.arg_count != 4) return NULL;
    Expr* const* a = res->data.function.args;
    if (a[0]->type != EXPR_STRING) return NULL;
    if (!is_symbol(a[1], SYM_True) && !is_symbol(a[1], SYM_False)) return NULL;
    if (!is_symbol(a[2], SYM_True) && !is_symbol(a[2], SYM_False)) return NULL;
    if (!(a[3]->type == EXPR_STRING || is_symbol(a[3], SYM_None))) return NULL;
    const char* s = a[0]->data.string;
    size_t n = strlen(s);
    const char* prelude = a[3]->type == EXPR_STRING ? a[3]->data.string : NULL;

    CheckErr err = { NULL, NULL, 0 };
    CmdState cs;
    memset(&cs, 0, sizeof cs);
    cs.mathlib = is_symbol(a[1], SYM_True);
    cs.certify = is_symbol(a[2], SYM_True);
    /* At most one example per line, so n / 8 + 1 bounds their number. */
    cs.examples = malloc((n / 8 + 2) * sizeof *cs.examples);
    if (!cs.examples) abort();
    bool ok = check_source(s, n, prelude, &cs, &err);
    Expr* out;
    if (ok) {
        /* Rename the examples (same line; the source grows). */
        const char* rep = "noncomputable def _root_.mathilda_example_";
        size_t extra = (size_t)cs.nexamples * (strlen(rep) + 24);
        char* buf = malloc(n + extra + 1);
        if (!buf) abort();
        size_t w = 0, r = 0;
        for (int k = 0; k < cs.nexamples; k++) {
            memcpy(buf + w, s + r, cs.examples[k] - r);
            w += cs.examples[k] - r;
            w += (size_t)sprintf(buf + w, "%s%d", rep, k + 1);
            r = cs.examples[k] + strlen("example");
        }
        memcpy(buf + w, s + r, n - r);
        w += n - r;
        buf[w] = '\0';
        Expr** names = malloc((cs.thms.n + 1) * sizeof *names);
        if (!names) abort();
        for (size_t k = 0; k < cs.thms.n; k++) names[k] = expr_new_string(cs.thms.v[k]);
        Expr* rules[3] = {
            rule("OK", expr_new_symbol(SYM_True)),
            rule("Source", expr_new_string(buf)),
            rule("Theorems", expr_new_function(expr_new_symbol(SYM_List), names, cs.thms.n)),
        };
        free(names);
        free(buf);
        out = assoc_of(rules, 3);
    } else {
        Expr* rules[4] = {
            rule("OK", expr_new_symbol(SYM_False)),
            rule("Kind", expr_new_string(err.kind ? err.kind : "lex")),
            rule("What", expr_new_string(err.what ? err.what : "")),
            rule("Line", expr_new_integer(err.line)),
        };
        out = assoc_of(rules, 4);
    }
    free(err.what);
    for (size_t k = 0; k < cs.st.n; k++) free(cs.st.v[k]);
    free(cs.st.v);
    for (size_t k = 0; k < cs.thms.n; k++) free(cs.thms.v[k]);
    free(cs.thms.v);
    free(cs.examples);
    return out;
}

void lean_check_init(void) {
    symtab_add_builtin("Lean`Internal`CheckSource", builtin_lean_checksource);
    symtab_set_docstring("Lean`Internal`CheckSource",
        "Lean`Internal`CheckSource[src, mathlib, certify, prelude]\n"
        "\tchecks the Lean 4 source src against LeanVerify's minimal surface (no\n"
        "\tstrings, comments, character literals, #, $, %, backquotes, \\, @[, French\n"
        "\tquotes or unlisted Unicode; no command keyword or escape word out of place)\n"
        "\tand allowlist grammar (theorem / example / open / namespace / section /\n"
        "\tvariable / universe / numeric set_option at column 0; imports of Mathlib\n"
        "\tonly when mathlib is True) and gives <|\"OK\" -> True, \"Source\",\n"
        "\t\"Theorems\"|> or <|\"OK\" -> False, \"Kind\", \"What\", \"Line\"|>. certify ->\n"
        "\tTrue also accepts Lean.m's generator auxiliaries (def mathilda_cert_k_x)\n"
        "\tand footer lines; prelude (or None) is a trusted text skipped once when it\n"
        "\toccurs byte-identical at a line start.");
    symtab_get_def("Lean`Internal`CheckSource")->attributes |= ATTR_PROTECTED;
}
