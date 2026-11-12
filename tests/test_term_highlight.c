/*
 * test_term_highlight.c -- the editor's display widths (src/term/term_width.c)
 * and its incremental line tokenizer (src/term/term_highlight.c): token
 * classes, state carried across lines, symbol classes from the symbol table
 * (without creating symbols), and the coloured fputs used for results.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "term/term_width.h"
#include "term/term_highlight.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void test_width(void) {
    ASSERT(term_cp_width('a') == 1);
    ASSERT(term_cp_width(0x01) == -1 && term_cp_width(0x7F) == -1);
    ASSERT(term_cp_width(0x4E2D) == 2);          /* CJK 中 */
    ASSERT(term_cp_width(0xAC00) == 2);          /* Hangul */
    ASSERT(term_cp_width(0xFF21) == 2);          /* fullwidth A */
    ASSERT(term_cp_width(0x1F600) == 2);         /* emoji */
    ASSERT(term_cp_width(0x0301) == 0);          /* combining acute */
    ASSERT(term_cp_width(0x200D) == 0);          /* ZWJ */
    ASSERT(term_cp_width(0x03B1) == 1);          /* alpha */
    ASSERT(term_cp_width(0xF7C0) == 1);          /* private use (linear boxes) */

    const char* s = "a\xE4\xB8\xAD" "e\xCC\x81" "\x01";   /* a 中 é(combining) ^A */
    size_t n = strlen(s);
    ASSERT(term_text_width(s, n) == 1 + 2 + 1 + 0 + 2);

    uint32_t cp;
    size_t i = term_utf8_next(s, n, 1, &cp);
    ASSERT(cp == 0x4E2D && i == 4);
    ASSERT(term_utf8_prev(s, 4) == 1);
    ASSERT(term_utf8_prev(s, 1) == 0);
    /* A stray continuation byte and a truncated sequence are one char each. */
    const char* bad = "x\x80\xE4\xB8";
    i = term_utf8_next(bad, 4, 1, &cp);
    ASSERT(i == 2 && cp == 0x80);
    i = term_utf8_next(bad, 4, 2, &cp);
    ASSERT(i == 3 && cp == 0xE4);
    ASSERT(term_utf8_prev(bad, 4) == 3);
}

/* Styles of `src` (one line) as a string of class letters. */
static const char LETTER[] = ".BUNSCPK";
static void styles(const char* src, TermHlState* st, char* out) {
    size_t n = strlen(src);
    uint8_t sty[256];
    *st = term_hl_line(src, n, *st, sty, term_hl_classify_symtab, NULL);
    for (size_t i = 0; i < n; i++) out[i] = LETTER[sty[i]];
    out[n] = '\0';
}

static void expect(const char* src, const char* want) {
    TermHlState st = term_hl_state_initial();
    char got[256];
    styles(src, &st, got);
    ASSERT_MSG(strcmp(got, want) == 0, "%s\n  got  %s\n  want %s", src, got, want);
}

static void test_tokens(void) {
    expect("Sin[x] + 12", "BBBK.K...NN");
    expect("1.5`20 2*^-3 16^^ff .5", "NNNNNN.NNNNN.NNNNNN.NN");
    expect("f[x_, _Integer, y__., #2 &]", ".KPP..PPPPPPPP..PPPP..PP..K");
    expect("\"a(*b\" (* \"c *)", "SSSSSS.CCCCCCCC");
    expect("(* a (* b *) c *) 1", "CCCCCCCCCCCCCCCCC.N");
    expect("{a, b}", "K....K");
    expect("\"x\\\"y\"", "SSSSSS");
    expect("1..3", "N..N");
}

static void test_carried_state(void) {
    TermHlState st = term_hl_state_initial();
    char got[256];
    styles("x = \"open", &st, got);
    ASSERT(st.in_string == 1);
    ASSERT(strcmp(got, "....SSSSS") == 0);
    styles("still\" + 1", &st, got);
    ASSERT(st.in_string == 0);
    ASSERT_MSG(strcmp(got, "SSSSSS...N") == 0, "%s", got);

    styles("(* one (* two", &st, got);
    ASSERT(st.comment == 2);
    styles("*) *) Sin", &st, got);
    ASSERT(st.comment == 0);
    ASSERT_MSG(strcmp(got, "CCCCC.BBB") == 0, "%s", got);

    TermHlState a = term_hl_state_initial(), b = term_hl_state_initial();
    ASSERT(term_hl_state_eq(a, b));
    b.invalid = 1;
    ASSERT(!term_hl_state_eq(a, b));
}

static size_t names_length(void) {
    Expr* p = parse_expression("Length[Names[]]");
    Expr* r = evaluate(p);
    expr_free(p);
    ASSERT(r && r->type == EXPR_INTEGER);
    size_t n = (size_t)r->data.integer;
    expr_free(r);
    return n;
}

static void test_symbol_classes(void) {
    Expr* p = parse_expression("hlUserF[x_] := x^2");
    Expr* r = evaluate(p);
    expr_free(p); expr_free(r);
    expect("hlUserF[1] + hlUnknownG", "UUUUUUUKNK.............");
    /* Lexing never creates a symbol. */
    size_t before = names_length();
    for (int i = 0; i < 50; i++) expect("neverSeenHlSymbol123", "....................");
    ASSERT(names_length() == before);
    /* A NULL classifier leaves every symbol plain. */
    uint8_t sty[8];
    (void)term_hl_line("Sin", 3, term_hl_state_initial(), sty, NULL, NULL);
    ASSERT(sty[0] == TERM_HL_PLAIN);
}

static void test_fputs(void) {
    char buf[512];
    memset(buf, 0, sizeof(buf));
    FILE* f = tmpfile();
    ASSERT(f != NULL);
    const char* text = "1/2 x^2 + Sin[y]";
    ASSERT(term_hl_fputs(text, strlen(text), f, term_hl_classify_symtab, NULL));
    rewind(f);
    size_t got = fread(buf, 1, sizeof(buf) - 1, f);
    fclose(f);
    buf[got] = '\0';
    ASSERT(strstr(buf, "\033[35m1") != NULL);         /* number colour */
    ASSERT(strstr(buf, "\033[36mSin") != NULL);       /* builtin colour */
    /* The text itself is unchanged once the SGR is removed. */
    char plain[512];
    size_t k = 0;
    for (size_t i = 0; buf[i]; i++) {
        if (buf[i] == '\033') { while (buf[i] && buf[i] != 'm') i++; continue; }
        plain[k++] = buf[i];
    }
    plain[k] = '\0';
    ASSERT_MSG(strcmp(plain, text) == 0, "%s", plain);
    /* The last SGR is a reset: no colour leaks into what follows. */
    char* last = strrchr(buf, '\033');
    ASSERT(last && strncmp(last, "\033[0m", 4) == 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_width);
    TEST(test_tokens);
    TEST(test_carried_state);
    TEST(test_symbol_classes);
    TEST(test_fputs);
    printf("All term_highlight tests passed.\n");
    return 0;
}
