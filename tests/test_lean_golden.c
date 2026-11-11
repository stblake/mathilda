/*
 * test_lean_golden.c -- byte identity of everything Certify generates.
 *
 * Loads tests/lean_golden/dump.m and corpus.m, evaluates goldenDump[] (one
 * record per corpus item: the messages issued, the result's head, and every
 * ProofObject property except "Timing", the inlined prelude elided) and
 * compares it with tests/lean_golden/golden.txt.  No Lean process is run.
 *
 * On a mismatch it writes the actual record to lean_golden_actual.txt in the
 * current directory, prints the first differing line (number and both
 * texts) and exits 1.  With LEAN_GOLDEN_UPDATE=1 in the environment it
 * rewrites golden.txt instead and exits 0; that is done only when a change
 * to the generated certificates is intended.
 */
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef LEAN_GOLDEN_DIR
#define LEAN_GOLDEN_DIR "lean_golden/"
#endif

/* The whole file at path as a NUL-terminated buffer (caller frees), or NULL. */
static char* read_file(const char* path, size_t* len_out) {
    FILE* f = fopen(path, "rb");
    if (!f) return NULL;
    size_t cap = 1 << 16, len = 0;
    char* buf = malloc(cap);
    if (!buf) { fclose(f); return NULL; }
    for (;;) {
        if (len + 4096 + 1 > cap) {
            char* nb = realloc(buf, cap * 2);
            if (!nb) { free(buf); fclose(f); return NULL; }
            buf = nb; cap *= 2;
        }
        size_t n = fread(buf + len, 1, 4096, f);
        len += n;
        if (n < 4096) break;
    }
    fclose(f);
    buf[len] = '\0';
    if (len_out) *len_out = len;
    return buf;
}

static int write_file(const char* path, const char* data, size_t len) {
    FILE* f = fopen(path, "wb");
    if (!f) return 0;
    int ok = fwrite(data, 1, len, f) == len;
    if (fclose(f) != 0) ok = 0;
    return ok;
}

/* Evaluate input and discard the value. */
static void run(const char* input) {
    Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "FAIL: does not parse: %s\n", input); exit(1); }
    Expr* v = evaluate(p);
    expr_free(p);
    if (v) expr_free(v);
}

/* The 1-based line number and the start of the first line where a and b
 * differ; *la / *lb receive the lengths of those lines. */
static size_t first_diff(const char* a, const char* b, const char** pa, const char** pb,
                         size_t* la, size_t* lb) {
    size_t line = 1;
    const char* sa = a;
    const char* sb = b;
    while (*a && *a == *b) {
        if (*a == '\n') { line++; sa = a + 1; sb = b + 1; }
        a++; b++;
    }
    *pa = sa; *pb = sb;
    *la = strcspn(sa, "\n");
    *lb = strcspn(sb, "\n");
    return line;
}

int main(void) {
    symtab_init();
    core_init();

    char in[1024];
    snprintf(in, sizeof in, "Get[\"%sdump.m\"]; Get[\"%scorpus.m\"]; Null",
             LEAN_GOLDEN_DIR, LEAN_GOLDEN_DIR);
    run(in);

    /* The corpus issues many messages on purpose; $MessageList records
     * their names, so their printed text is not needed here. */
    int saved = test_stderr_redirect("/dev/null");
    Expr* p = parse_expression("goldenDump[]");
    Expr* v = p ? evaluate(p) : NULL;
    if (p) expr_free(p);
    test_stderr_restore(saved);
    if (!v || v->type != EXPR_STRING) {
        char* s = v ? expr_to_string(v) : NULL;
        fprintf(stderr, "FAIL: goldenDump[] did not return a String: %.300s\n",
                s ? s : "(null)");
        free(s);
        if (v) expr_free(v);
        return 1;
    }
    const char* actual = v->data.string;
    size_t alen = strlen(actual);

    char golden_path[1024];
    snprintf(golden_path, sizeof golden_path, "%sgolden.txt", LEAN_GOLDEN_DIR);

    const char* upd = getenv("LEAN_GOLDEN_UPDATE");
    if (upd && strcmp(upd, "1") == 0) {
        int ok = write_file(golden_path, actual, alen);
        printf("%s %s (%zu bytes)\n", ok ? "wrote" : "FAIL: could not write", golden_path,
               alen);
        expr_free(v);
        return ok ? 0 : 1;
    }

    size_t glen = 0;
    char* golden = read_file(golden_path, &glen);
    if (!golden) {
        fprintf(stderr, "FAIL: cannot read %s\n", golden_path);
        expr_free(v);
        return 1;
    }
    int same = glen == alen && memcmp(golden, actual, alen) == 0;
    if (!same) {
        const char *pa, *pb;
        size_t la, lb;
        size_t line = first_diff(golden, actual, &pa, &pb, &la, &lb);
        write_file("lean_golden_actual.txt", actual, alen);
        fprintf(stderr, "FAIL: the generated certificates differ from %s at line %zu\n"
                "  golden: %.*s\n  actual: %.*s\n"
                "  (full record in lean_golden_actual.txt; LEAN_GOLDEN_UPDATE=1 rewrites "
                "golden.txt when the change is intended)\n",
                golden_path, line, (int)(la > 2000 ? 2000 : la), pa,
                (int)(lb > 2000 ? 2000 : lb), pb);
    } else {
        size_t records = 0;
        for (const char* s = actual; (s = strstr(s, "=== ")) != NULL; s += 4)
            if (s == actual || s[-1] == '\n') records++;
        printf("lean_golden: %zu records identical to %s\n", records, golden_path);
    }
    free(golden);
    expr_free(v);
    return same ? 0 : 1;
}
