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
 * then reads.  No deny-list of words can stop that, so LeanVerify accepts only
 * sources that cannot contain user-defined code at all: an ALLOWLIST grammar.
 *
 * Accepted grammar.  Commands start at column 0; every continuation line is
 * indented (a token at column 0 starts a new command).  Tabs, carriage
 * returns and other control characters are refused.  Each command starts with
 *
 *     import Mathlib | import Mathlib.<Path>    only with `mathlib`, only first
 *     theorem NAME ... | lemma NAME ...          (lemma only with `mathlib`)
 *     example ...                                (renamed, see below)
 *     open [scoped] A B ... [in]                 namespace names only
 *     namespace A.B | section [A] | noncomputable section [A] | end [A]
 *     variable ... | universe u ...
 *     set_option maxHeartbeats|maxRecDepth|synthInstance.maxHeartbeats <n> [in]
 *     set_option autoImplicit|relaxedAutoImplicit|linter.<x> true|false [in]
 *     #print axioms NAME
 *     def mathilda_cert_<k>_<x> ...              only with `certify` (generator
 *                                                auxiliaries of Lean.m)
 *
 * and comments / docstrings anywhere.  Everything else (def, abbrev,
 * instance, structure, inductive, class, axiom, opaque, macro, syntax,
 * notation, elab, attribute, @[...], initialize, #eval, #check, run_cmd,
 * simproc, mutual, module, ...) is refused ("cmd").  Inside an accepted
 * command, the words that escape to code or bypass the kernel (by_elab,
 * run_tac, native_decide, native, sorry, admit, unsafe, IO, ...), every
 * command keyword (so that a command hidden on a continuation line, which
 * Lean would parse as a new command, is refused too), `@[`, `#word`,
 * `word%`, backquotes (name literals and syntax quotations) and interpolated
 * strings are refused ("token").  Words are scanned as maximal runs of
 * [A-Za-z_][A-Za-z0-9_'] outside comments and string / character literals,
 * which splits at least as finely as Lean's own identifiers, so a keyword
 * Lean would see is always seen here too.  An escaped identifier <<...>>
 * (French quotes) is accepted only on one line and without a refused word
 * inside.  A source whose lexical structure is ambiguous to this scanner (an
 * unterminated comment, string or escaped identifier) is refused ("lex").
 *
 * The trusted prelude.  `prelude` (the text of src/internal/lean/Prelude.lean)
 * may occur once, byte-identical, at the start of a line at the top level; it
 * is skipped (its definitions are Mathilda's own, reviewed code).
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
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ */
/* Tokens                                                              */
/* ------------------------------------------------------------------ */

enum { TK_WORD, TK_NUM, TK_STR, TK_CHAR, TK_SYM, TK_QID };

typedef struct {
    size_t start, end;   /* byte range in the source */
    int line;            /* 1-based */
    bool col0;           /* first byte of its line */
    int kind;
} Tok;

typedef struct {
    Tok* v;
    size_t n, cap;
} TokList;

typedef struct {
    const char* kind;    /* "lex", "token", "cmd", "import" */
    char* what;          /* owned */
    int line;
} CheckErr;

static void tok_push(TokList* l, Tok t) {
    if (l->n == l->cap) {
        size_t nc = l->cap ? l->cap * 2 : 256;
        Tok* p = realloc(l->v, nc * sizeof *p);
        if (!p) abort();
        l->v = p;
        l->cap = nc;
    }
    l->v[l->n++] = t;
}

static char* copy_n(const char* s, size_t n) {
    char* r = malloc(n + 1);
    if (!r) abort();
    memcpy(r, s, n);
    r[n] = '\0';
    return r;
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

static bool word_start(unsigned char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_';
}

static bool word_char(unsigned char c) {
    return word_start(c) || (c >= '0' && c <= '9') || c == '\'';
}

/* Length of the UTF-8 sequence starting with c (1 for a stray byte). */
static size_t utf8_len(unsigned char c) {
    if (c < 0x80) return 1;
    if ((c & 0xE0) == 0xC0) return 2;
    if ((c & 0xF0) == 0xE0) return 3;
    if ((c & 0xF8) == 0xF0) return 4;
    return 1;
}

/* Words refused inside every command: meta-level escapes, kernel bypasses,
 * code-running tactics, and every command keyword (Lean does not need a new
 * line to start a command, so a command keyword anywhere would start one). */
static const char* const REFUSED_WORDS[] = {
    /* escapes to code / kernel bypasses */
    "by_elab", "run_tac", "run_elab", "run_meta", "run_cmd", "native_decide", "native",
    "implemented_by", "extern", "unsafe", "sorry", "admit", "ofReduceBool", "ofReduceNat",
    "reduceBool", "trustCompiler", "debug", "skipKernelTC", "addDecl", "addDeclCore",
    "addAndCompile", "addDeclWithoutChecking", "setEnv", "modifyEnv", "withOptions",
    "evalExpr", "evalConst", "unsafeBaseIO", "unsafeIO", "unsafeEIO", "unsafePerformIO",
    "unsafeCast", "IO", "EIO", "BaseIO", "Elab", "MetaM", "TermElabM", "TacticM",
    "CommandElabM", "CoreM", "MacroM", "command_elab", "term_elab", "tactic_elab",
    "plausible", "slim_check", "polyrith",
    /* commands and declaration modifiers */
    "attribute", "set_option", "open", "export", "elab", "elab_rules", "macro",
    "macro_rules", "syntax", "declare_syntax_cat", "initialize", "builtin_initialize",
    "axiom", "theorem", "lemma", "def", "abbrev", "instance", "example", "structure",
    "class", "inductive", "coinductive", "opaque", "notation", "notation3", "infix",
    "infixl", "infixr", "prefix", "postfix", "namespace", "section", "end", "universe",
    "variable", "import", "mutual", "deriving", "private", "protected", "noncomputable",
    "partial", "nonrec", "meta", "public", "module", "omit", "include", "scoped", "local",
    "simproc", "dsimproc", "simproc_decl", "dsimproc_decl", "builtin_simproc",
    "add_decl_doc", "alias", "irreducible_def", "register_option", "declare_config_elab",
    "grind_pattern", "recall", "library_note", "assert_not_exists",
    "initialize_simps_projections", "seal", "unseal", "compile_inductive",
    "suppress_compilation", "binder_predicate",
    NULL
};

#define REFUSED_WORDS_PTR REFUSED_WORDS

static bool tok_is_raw(const char* s, const Tok* t, const char* w) {
    size_t n = strlen(w);
    return t->end - t->start == n && memcmp(s + t->start, w, n) == 0;
}

/* Lex src into tokens; comments are dropped, the prelude (if it occurs at a
 * line start at the top level) is skipped once.  false on a lexical refusal. */
static bool lex(const char* s, size_t n, const char* prelude, TokList* out, CheckErr* err,
                size_t* pre_at, size_t* pre_len) {
    size_t i = 0, ls = 0;        /* position, start of the current line */
    int line = 1;
    size_t plen = prelude ? strlen(prelude) : 0;
    *pre_at = (size_t)-1;
    *pre_len = 0;
    while (i < n) {
        unsigned char c = (unsigned char)s[i];
        if (c == '\n') { line++; i++; ls = i; continue; }
        if (c == ' ') { i++; continue; }
        if (c < 0x20 || c == 0x7F) {
            set_err_s(err, "lex", line, c == '\t' ? "a tab" : c == '\r' ? "a carriage return"
                                                            : "a control character");
            return false;
        }
        if (i == ls && plen > 0 && *pre_at == (size_t)-1 && n - i >= plen &&
            memcmp(s + i, prelude, plen) == 0) {
            *pre_at = i;
            *pre_len = plen;
            for (size_t k = 0; k < plen; k++)
                if (prelude[k] == '\n') { line++; ls = i + k + 1; }
            i += plen;
            continue;
        }
        if (c == '-' && i + 1 < n && s[i + 1] == '-') {          /* line comment */
            while (i < n && s[i] != '\n') i++;
            continue;
        }
        if (c == '/' && i + 1 < n && s[i + 1] == '-') {          /* nested block comment */
            int depth = 1, l0 = line;
            i += 2;
            while (i < n && depth > 0) {
                if (s[i] == '/' && i + 1 < n && s[i + 1] == '-') { depth++; i += 2; }
                else if (s[i] == '-' && i + 1 < n && s[i + 1] == '/') { depth--; i += 2; }
                else { if (s[i] == '\n') { line++; ls = i + 1; } i++; }
            }
            if (depth > 0) { set_err_s(err, "lex", l0, "an unterminated comment"); return false; }
            continue;
        }
        Tok t = { i, i, line, i == ls, TK_SYM };
        if (c == '"') {                                          /* string literal */
            if (i > 0 && s[i - 1] == '!') {
                set_err_s(err, "token", line, "an interpolated string");
                return false;
            }
            size_t j = i + 1;
            int l0 = line;
            while (j < n && s[j] != '"') {
                if (s[j] == '\\' && j + 1 < n) j++;
                if (s[j] == '\n') { line++; ls = j + 1; }
                j++;
            }
            if (j >= n) { set_err_s(err, "lex", l0, "an unterminated string"); return false; }
            t.end = j + 1;
            t.kind = TK_STR;
            tok_push(out, t);
            i = j + 1;
            continue;
        }
        if (c == '\'') {                                         /* char literal or prime */
            bool prime = i > 0 && (word_char((unsigned char)s[i - 1]) ||
                                   (unsigned char)s[i - 1] >= 0x80);
            if (!prime && i + 2 < n && s[i + 1] != '\'' && s[i + 1] != '\n') {
                size_t j = i + 1;
                if (s[j] == '\\') {
                    j++;
                    if (j < n && s[j] == 'x') j += 3;
                    else if (j < n && s[j] == 'u') { while (j < n && s[j] != '}' && s[j] != '\n') j++; j++; }
                    else j++;
                } else {
                    j += utf8_len((unsigned char)s[j]);
                }
                if (j < n && s[j] == '\'') {
                    t.end = j + 1;
                    t.kind = TK_CHAR;
                    tok_push(out, t);
                    i = j + 1;
                    continue;
                }
            }
            t.end = i + 1;
            tok_push(out, t);
            i++;
            continue;
        }
        if (c == '`') {
            set_err_s(err, "token", line, "a backquote (name literal or syntax quotation)");
            return false;
        }
        if (c == 0xC2 && i + 1 < n && (unsigned char)s[i + 1] == 0xAB) {
            /* An escaped identifier <<...>> (Mathilda writes a Greek variable
             * as one).  It is always an identifier, never a keyword, so it
             * cannot start a command or an escape; it is accepted only when
             * it ends on the same line and contains no space, backquote,
             * French quote or refused word (so that the scanner and Lean
             * agree on where it ends). */
            size_t j = i + 2;
            bool bad = false;
            while (j + 1 < n && !((unsigned char)s[j] == 0xC2 && (unsigned char)s[j + 1] == 0xBB)) {
                unsigned char d = (unsigned char)s[j];
                if (d == '\n' || d == ' ' || d == '`' || d < 0x20 ||
                    (d == 0xC2 && (unsigned char)s[j + 1] == 0xAB))
                    bad = true;
                j++;
            }
            if (bad || j + 1 >= n) {
                set_err_s(err, "lex", line, "an unterminated or malformed \xC2\xAB\xC2\xBB identifier");
                return false;
            }
            for (size_t k = i + 2; k < j; ) {
                if (word_start((unsigned char)s[k]) && (k == i + 2 || !word_char((unsigned char)s[k - 1]))) {
                    size_t m = k + 1;
                    while (m < j && word_char((unsigned char)s[m])) m++;
                    Tok w = { k, m, line, false, TK_WORD };
                    for (int q = 0; REFUSED_WORDS_PTR[q]; q++)
                        if (tok_is_raw(s, &w, REFUSED_WORDS_PTR[q])) {
                            set_err(err, "token", line, s + i, j + 2 - i);
                            return false;
                        }
                    k = m;
                } else {
                    k++;
                }
            }
            t.end = j + 2;
            t.kind = TK_QID;
            tok_push(out, t);
            i = t.end;
            continue;
        }
        if (c == 0xC2 && i + 1 < n && (unsigned char)s[i + 1] == 0xBB) {
            set_err_s(err, "lex", line, "a stray \xC2\xBB");
            return false;
        }
        if (c == '@' && i + 1 < n && s[i + 1] == '[') {
            set_err_s(err, "token", line, "@[ (an attribute)");
            return false;
        }
        if (word_start(c)) {
            size_t j = i + 1;
            while (j < n && word_char((unsigned char)s[j])) j++;
            /* raw string r"..." / r#"..."# */
            if (j == i + 1 && c == 'r' && j < n && (s[j] == '"' || s[j] == '#')) {
                size_t h = 0, k = j;
                while (k < n && s[k] == '#') { h++; k++; }
                if (k < n && s[k] == '"') {
                    int l0 = line;
                    size_t m = k + 1;
                    for (;;) {
                        if (m >= n) { set_err_s(err, "lex", l0, "an unterminated raw string"); return false; }
                        if (s[m] == '\n') { line++; ls = m + 1; }
                        if (s[m] == '"') {
                            size_t q = 0;
                            while (q < h && m + 1 + q < n && s[m + 1 + q] == '#') q++;
                            if (q == h) break;
                        }
                        m++;
                    }
                    t.end = m + 1 + h;
                    t.kind = TK_STR;
                    tok_push(out, t);
                    i = t.end;
                    continue;
                }
            }
            t.end = j;
            t.kind = TK_WORD;
            tok_push(out, t);
            i = j;
            continue;
        }
        if (c >= '0' && c <= '9') {
            size_t j = i + 1;
            while (j < n && (word_char((unsigned char)s[j]) ||
                             (s[j] == '.' && j + 1 < n && s[j + 1] >= '0' && s[j + 1] <= '9')))
                j++;
            t.end = j;
            t.kind = TK_NUM;
            tok_push(out, t);
            i = j;
            continue;
        }
        t.end = i + utf8_len(c);
        if (t.end > n) t.end = n;
        tok_push(out, t);
        i = t.end;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/* Command-level grammar                                               */
/* ------------------------------------------------------------------ */

static bool tok_is(const char* s, const Tok* t, const char* w) { return tok_is_raw(s, t, w); }

static bool refused_word(const char* s, const Tok* t) {
    for (int k = 0; REFUSED_WORDS[k]; k++)
        if (tok_is(s, t, REFUSED_WORDS[k])) return true;
    return false;
}

/* Scan tokens [a, b) of one command's body. */
static bool scan_body(const char* s, const TokList* tl, size_t a, size_t b, CheckErr* err) {
    for (size_t k = a; k < b; k++) {
        const Tok* t = &tl->v[k];
        const Tok* nx = k + 1 < b ? &tl->v[k + 1] : NULL;
        if (t->kind == TK_WORD && refused_word(s, t)) {
            set_err(err, "token", t->line, s + t->start, t->end - t->start);
            return false;
        }
        if (t->kind == TK_WORD && nx && nx->start == t->end && tok_is(s, nx, "%")) {
            set_err(err, "token", t->line, s + t->start, nx->end - t->start);
            return false;
        }
        if (t->kind == TK_SYM && tok_is(s, t, "#") && nx && nx->start == t->end &&
            nx->kind == TK_WORD) {
            set_err(err, "token", t->line, s + t->start, nx->end - t->start);
            return false;
        }
    }
    return true;
}

/* A dotted name W(.W)* starting at token k with no gaps; returns the index
 * after it (k when there is none). */
static size_t name_at(const char* s, const TokList* tl, size_t k, size_t b) {
    if (k >= b || tl->v[k].kind != TK_WORD) return k;
    size_t j = k + 1;
    while (j + 1 < b && tok_is(s, &tl->v[j], ".") && tl->v[j].start == tl->v[j - 1].end &&
           tl->v[j + 1].kind == TK_WORD && tl->v[j + 1].start == tl->v[j].end)
        j += 2;
    return j;
}

static char* tok_text(const char* s, const TokList* tl, size_t a, size_t b) {
    return copy_n(s + tl->v[a].start, tl->v[b - 1].end - tl->v[a].start);
}

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

static bool all_digits(const char* s, const Tok* t) {
    if (t->kind != TK_NUM) return false;
    for (size_t k = t->start; k < t->end; k++)
        if (s[k] < '0' || s[k] > '9') return false;
    return true;
}

/* set_option NAME VALUE [in], from token a (the keyword) to b. */
static bool option_ok(const char* s, const TokList* tl, size_t a, size_t b) {
    size_t e = name_at(s, tl, a + 1, b);
    if (e == a + 1 || e >= b) return false;
    char* nm = tok_text(s, tl, a + 1, e);
    const Tok* v = &tl->v[e];
    bool numeric = strcmp(nm, "maxHeartbeats") == 0 || strcmp(nm, "maxRecDepth") == 0 ||
                   strcmp(nm, "synthInstance.maxHeartbeats") == 0;
    bool boolean = strcmp(nm, "autoImplicit") == 0 || strcmp(nm, "relaxedAutoImplicit") == 0 ||
                   (strncmp(nm, "linter.", 7) == 0 && strlen(nm) > 7);
    free(nm);
    bool ok = (numeric && all_digits(s, v)) ||
              (boolean && (tok_is(s, v, "true") || tok_is(s, v, "false")));
    size_t rest = e + 1;
    if (rest < b && tok_is(s, &tl->v[rest], "in")) rest++;
    return ok && rest == b;
}

/* Is `nm` a generator auxiliary name mathilda_cert_<digits>_<word>? */
static bool cert_aux_name(const char* nm) {
    const char* p = "mathilda_cert_";
    if (strncmp(nm, p, strlen(p)) != 0) return false;
    nm += strlen(p);
    if (*nm < '0' || *nm > '9') return false;
    while (*nm >= '0' && *nm <= '9') nm++;
    if (*nm != '_' || !nm[1]) return false;
    for (nm++; *nm; nm++)
        if (!word_char((unsigned char)*nm) || *nm == '\'') return false;
    return true;
}

typedef StrStack Names;

static void names_push(Names* l, char* s) { ss_push(l, s); }

/* Check one command, tokens [a, b). */
static bool check_command(const char* s, const TokList* tl, size_t a, size_t b, bool mathlib,
                          bool certify, bool* past_imports, StrStack* st, Names* thms,
                          size_t* examples, int* nexamples, CheckErr* err) {
    const Tok* t0 = &tl->v[a];
    int line = t0->line;
    if (t0->kind == TK_WORD && tok_is(s, t0, "import")) {
        size_t e = name_at(s, tl, a + 1, b);
        char* mod = e > a + 1 ? tok_text(s, tl, a + 1, e) : copy_n("", 0);
        bool ok = mathlib && !*past_imports && e == b && e > a + 1 &&
                  (strcmp(mod, "Mathlib") == 0 || strncmp(mod, "Mathlib.", 8) == 0);
        if (!ok) {
            char* what = tok_text(s, tl, a, b);
            set_err_s(err, "import", line, what);
            free(what);
        }
        free(mod);
        return ok;
    }
    *past_imports = true;
    if (t0->kind == TK_SYM && tok_is(s, t0, "#")) {
        if (a + 2 < b && tok_is(s, &tl->v[a + 1], "print") && tl->v[a + 1].start == t0->end &&
            tok_is(s, &tl->v[a + 2], "axioms") && name_at(s, tl, a + 3, b) == b && a + 3 < b)
            return true;
        size_t e = a + 1 < b ? a + 2 : a + 1;
        char* what = tok_text(s, tl, a, e);
        set_err_s(err, "cmd", line, what);
        free(what);
        return false;
    }
    if (t0->kind != TK_WORD) {
        char* what = tok_text(s, tl, a, a + 1);
        set_err_s(err, "cmd", line, what);
        free(what);
        return false;
    }
    if (tok_is(s, t0, "theorem") || (mathlib && tok_is(s, t0, "lemma"))) {
        size_t e = name_at(s, tl, a + 1, b);
        if (e == a + 1) { set_err_s(err, "cmd", line, "theorem without a name"); return false; }
        char* nm = tok_text(s, tl, a + 1, e);
        names_push(thms, qualified(st, nm));
        free(nm);
        return scan_body(s, tl, a + 1, b, err);
    }
    if (tok_is(s, t0, "example")) {
        char nm[64];
        (*nexamples)++;
        snprintf(nm, sizeof nm, "mathilda_example_%d", *nexamples);
        names_push(thms, copy_n(nm, strlen(nm)));
        examples[*nexamples - 1] = t0->start;
        return scan_body(s, tl, a + 1, b, err);
    }
    if (certify && tok_is(s, t0, "def")) {
        size_t e = name_at(s, tl, a + 1, b);
        char* nm = e > a + 1 ? tok_text(s, tl, a + 1, e) : copy_n("", 0);
        bool ok = cert_aux_name(nm);
        free(nm);
        if (!ok) { set_err_s(err, "cmd", line, "def"); return false; }
        return scan_body(s, tl, a + 1, b, err);
    }
    if (tok_is(s, t0, "namespace")) {
        size_t e = name_at(s, tl, a + 1, b);
        if (e == a + 1 || e != b) { set_err_s(err, "cmd", line, "namespace"); return false; }
        ss_push(st, tok_text(s, tl, a + 1, e));
        return scan_body(s, tl, a + 1, b, err);
    }
    if (tok_is(s, t0, "section") ||
        (tok_is(s, t0, "noncomputable") && a + 1 < b && tok_is(s, &tl->v[a + 1], "section"))) {
        size_t k = tok_is(s, t0, "section") ? a + 1 : a + 2;
        size_t e = name_at(s, tl, k, b);
        if (e != b) { set_err_s(err, "cmd", line, "section"); return false; }
        ss_push(st, copy_n("", 0));
        return scan_body(s, tl, k, b, err);
    }
    if (tok_is(s, t0, "end")) {
        if (name_at(s, tl, a + 1, b) != b) { set_err_s(err, "cmd", line, "end"); return false; }
        if (st->n > 0) free(st->v[--st->n]);
        return scan_body(s, tl, a + 1, b, err);
    }
    if (tok_is(s, t0, "variable") || tok_is(s, t0, "universe"))
        return scan_body(s, tl, a + 1, b, err);
    if (tok_is(s, t0, "open")) {
        size_t k = a + 1;
        if (k < b && tok_is(s, &tl->v[k], "scoped")) k++;
        size_t first = k;
        while (k < b) {
            if (tok_is(s, &tl->v[k], "in") && k + 1 == b && k > first) { k++; break; }
            size_t e = name_at(s, tl, k, b);
            if (e == k) break;
            k = e;
        }
        if (k != b || first == b) {
            set_err_s(err, "cmd", line, "open (only namespace names may follow)");
            return false;
        }
        for (size_t j = first; j < b; j++)
            if (tl->v[j].kind == TK_WORD && !tok_is(s, &tl->v[j], "in") && refused_word(s, &tl->v[j])) {
                set_err(err, "token", line, s + tl->v[j].start, tl->v[j].end - tl->v[j].start);
                return false;
            }
        return true;
    }
    if (tok_is(s, t0, "set_option")) {
        if (option_ok(s, tl, a, b)) return true;
        char* what = tok_text(s, tl, a, b);
        set_err_s(err, "cmd", line, what);
        free(what);
        return false;
    }
    char* what = tok_text(s, tl, a, a + 1);
    set_err_s(err, "cmd", line, what);
    free(what);
    return false;
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

static Expr* builtin_lean_checksource(Expr* res) {
    if (res->data.function.arg_count != 4) return NULL;
    Expr* const* a = res->data.function.args;
    if (a[0]->type != EXPR_STRING) return NULL;
    if (!is_symbol(a[1], SYM_True) && !is_symbol(a[1], SYM_False)) return NULL;
    if (!is_symbol(a[2], SYM_True) && !is_symbol(a[2], SYM_False)) return NULL;
    if (!(a[3]->type == EXPR_STRING || is_symbol(a[3], SYM_None))) return NULL;
    const char* s = a[0]->data.string;
    size_t n = strlen(s);
    bool mathlib = is_symbol(a[1], SYM_True), certify = is_symbol(a[2], SYM_True);
    const char* prelude = a[3]->type == EXPR_STRING ? a[3]->data.string : NULL;

    CheckErr err = { NULL, NULL, 0 };
    TokList tl = { NULL, 0, 0 };
    StrStack st = { NULL, 0, 0 };
    Names thms = { NULL, 0, 0 };
    size_t pre_at, pre_len;
    size_t* examples = NULL;
    int nexamples = 0;
    /* Import lines first, so that the reason reported for a source that
     * imports a meta-level library is the import (not whatever the lexer
     * trips over later).  A line that starts with "import " at column 0 must
     * be a single Mathlib module (and only with `mathlib`); the same rule is
     * enforced again on the tokens below. */
    bool ok = true;
    for (size_t i = 0; ok && i < n; ) {
        size_t e = i;
        while (e < n && s[e] != '\n') e++;
        if (e - i >= 7 && memcmp(s + i, "import", 6) == 0 && (s[i + 6] == ' ')) {
            size_t k = i + 6;
            while (k < e && s[k] == ' ') k++;
            size_t m = k;
            while (m < e && (word_char((unsigned char)s[m]) || s[m] == '.')) m++;
            size_t t = m;
            while (t < e && s[t] == ' ') t++;
            bool good = mathlib && t == e && m > k && s[m - 1] != '.' &&
                        ((m - k == 7 && memcmp(s + k, "Mathlib", 7) == 0) ||
                         (m - k > 8 && memcmp(s + k, "Mathlib.", 8) == 0));
            if (!good) {
                int line = 1;
                for (size_t q = 0; q < i; q++) if (s[q] == '\n') line++;
                size_t len = e - i;
                while (len > 0 && s[i + len - 1] == ' ') len--;
                set_err(&err, "import", line, s + i, len);
                ok = false;
            }
        }
        i = e + 1;
    }
    if (ok) ok = lex(s, n, prelude, &tl, &err, &pre_at, &pre_len);
    if (ok) {
        examples = malloc((tl.n + 1) * sizeof *examples);
        if (!examples) abort();
        if (tl.n > 0 && !tl.v[0].col0) {
            set_err_s(&err, "cmd", tl.v[0].line, "an indented first command (commands start "
                                                 "at column 0)");
            ok = false;
        }
        bool past_imports = false;
        for (size_t k = 0; ok && k < tl.n; ) {
            size_t e = k + 1;
            while (e < tl.n && !tl.v[e].col0) e++;
            ok = check_command(s, &tl, k, e, mathlib, certify, &past_imports, &st, &thms,
                               examples, &nexamples, &err);
            k = e;
        }
    }
    Expr* out;
    if (ok) {
        /* Rename the examples (same line; the source grows). */
        const char* rep = "noncomputable def _root_.mathilda_example_";
        size_t extra = (size_t)nexamples * (strlen(rep) + 24);
        char* buf = malloc(n + extra + 1);
        if (!buf) abort();
        size_t w = 0, r = 0;
        for (int k = 0; k < nexamples; k++) {
            memcpy(buf + w, s + r, examples[k] - r);
            w += examples[k] - r;
            w += (size_t)sprintf(buf + w, "%s%d", rep, k + 1);
            r = examples[k] + strlen("example");
        }
        memcpy(buf + w, s + r, n - r);
        w += n - r;
        buf[w] = '\0';
        Expr** names = malloc((thms.n + 1) * sizeof *names);
        if (!names) abort();
        for (size_t k = 0; k < thms.n; k++) names[k] = expr_new_string(thms.v[k]);
        Expr* rules[3] = {
            rule("OK", expr_new_symbol(SYM_True)),
            rule("Source", expr_new_string(buf)),
            rule("Theorems", expr_new_function(expr_new_symbol(SYM_List), names, thms.n)),
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
    free(tl.v);
    for (size_t k = 0; k < st.n; k++) free(st.v[k]);
    free(st.v);
    for (size_t k = 0; k < thms.n; k++) free(thms.v[k]);
    free(thms.v);
    free(examples);
    return out;
}

void lean_check_init(void) {
    symtab_add_builtin("Lean`Internal`CheckSource", builtin_lean_checksource);
    symtab_set_docstring("Lean`Internal`CheckSource",
        "Lean`Internal`CheckSource[src, mathlib, certify, prelude]\n"
        "\tchecks the Lean 4 source src against LeanVerify's allowlist grammar (no\n"
        "\tuser-defined code: only theorem / example / open / namespace / section /\n"
        "\tvariable / universe / numeric set_option commands, imports of Mathlib only\n"
        "\twhen mathlib is True) and gives <|\"OK\" -> True, \"Source\", \"Theorems\"|>\n"
        "\tor <|\"OK\" -> False, \"Kind\", \"What\", \"Line\"|>. certify -> True also\n"
        "\taccepts Lean.m's generator auxiliaries (def mathilda_cert_k_x); prelude (or\n"
        "\tNone) is a trusted text skipped once when it occurs byte-identical.");
    symtab_get_def("Lean`Internal`CheckSource")->attributes |= ATTR_PROTECTED;
}
