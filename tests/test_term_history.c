/*
 * test_term_history.c -- the persistent REPL history (src/term/term_history.c):
 * the escaped one-line-per-entry format, the environment rules, loading,
 * appending, filtering, truncation at load, foreign files and the 0600 mode.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "test_utils.h"
#include "term/term_history.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

static char g_dir[128];

static void path_in_dir(char* out, size_t cap, const char* leaf) {
    snprintf(out, cap, "%s/%s", g_dir, leaf);
}

static char* slurp(const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* b = malloc((size_t)n + 1);
    size_t got = fread(b, 1, (size_t)n, f);
    b[got] = '\0';
    fclose(f);
    return b;
}

static void roundtrip(const char* s) {
    char* e = term_history_encode(s);
    ASSERT(e && !strchr(e, '\n') && !strchr(e, '\r'));
    char* d = term_history_decode(e);
    ASSERT_MSG(strcmp(d, s) == 0, "round trip of %s gave %s", s, d);
    free(e);
    free(d);
}

static void test_encoding(void) {
    /* 1. Round trips. */
    roundtrip("a");
    roundtrip("f[\n1]");
    roundtrip("\"a\\nb\"");            /* a literal backslash-n */
    roundtrip("x\r\ny");
    roundtrip("\\");
    roundtrip("");
    roundtrip("#^2 &");
    char* e = term_history_encode("f[\n1]");
    ASSERT_STR_EQ(e, "f[\\n1]");
    free(e);
    /* 2. An unknown escape is kept literally; so is a trailing backslash. */
    char* d = term_history_decode("\\q");
    ASSERT_STR_EQ(d, "\\q");
    free(d);
    d = term_history_decode("ab\\");
    ASSERT_STR_EQ(d, "ab\\");
    free(d);
}

static void test_env_rules(void) {
    /* 3. Paths. */
    char* p = term_history_resolve_path("", "/home/u");
    ASSERT(p == NULL);
    p = term_history_resolve_path(NULL, NULL);
    ASSERT(p == NULL);
    p = term_history_resolve_path(NULL, "");
    ASSERT(p == NULL);
    p = term_history_resolve_path("/x/h", "/home/u");
    ASSERT_STR_EQ(p, "/x/h");
    free(p);
    p = term_history_resolve_path(NULL, "/home/u");
    ASSERT_STR_EQ(p, "/home/u/.mathilda_history");
    free(p);
    /* 4. Sizes. */
    ASSERT(term_history_size_from_env("0") == 0);
    ASSERT(term_history_size_from_env("abc") == 1000);
    ASSERT(term_history_size_from_env("5") == 5);
    ASSERT(term_history_size_from_env("99999999") == 1000000);
    ASSERT(term_history_size_from_env("-3") == 1000);
    ASSERT(term_history_size_from_env("") == 1000);
    ASSERT(term_history_size_from_env(NULL) == 1000);
    ASSERT(term_history_open("/tmp/never-created-by-this-test", 0) == NULL);
}

static void test_store(void) {
    char path[256];
    path_in_dir(path, sizeof(path), "h1");

    /* 5. A fresh path: header written, adds appended, reopening loads them. */
    TermHistory* h = term_history_open(path, 100);
    ASSERT(h && term_history_persistent(h));
    char* text = slurp(path);
    ASSERT_STR_EQ(text, "#mathilda-history v1\n");
    free(text);
    ASSERT(term_history_add(h, "a1 = 11"));
    ASSERT(term_history_add(h, "f[\n1]"));
    /* 6. Repeats, blanks and comments are skipped. */
    ASSERT(!term_history_add(h, "f[\n1]"));
    ASSERT(!term_history_add(h, "   "));
    ASSERT(!term_history_add(h, "(* c *)"));
    ASSERT(term_history_add(h, "a1 = 11"));     /* not consecutive: kept */
    /* 10. A '#' entry. */
    ASSERT(term_history_add(h, "#^2 &"));
    term_history_close(h);
    text = slurp(path);
    ASSERT_STR_EQ(text, "#mathilda-history v1\na1 = 11\nf[\\n1]\na1 = 11\n#^2 &\n");
    free(text);

    h = term_history_open(path, 100);
    ASSERT(term_history_count(h) == 4);
    ASSERT_STR_EQ(term_history_get(h, 1), "f[\n1]");
    ASSERT_STR_EQ(term_history_get(h, 3), "#^2 &");
    ASSERT(term_history_get(h, 4) == NULL);
    /* The repeat check sees the loaded tail. */
    ASSERT(!term_history_add(h, "#^2 &"));
    term_history_close(h);

#ifndef _WIN32
    /* 9. Created 0600. */
    struct stat st;
    ASSERT(stat(path, &st) == 0);
    ASSERT_MSG((st.st_mode & 0777) == 0600, "mode %o", (unsigned)(st.st_mode & 0777));
#endif
    remove(path);
}

static void test_truncation(void) {
    /* 7. max 3 with 5 entries: the file keeps the header and the last 3. */
    char path[256], tmp[300];
    path_in_dir(path, sizeof(path), "h2");
    snprintf(tmp, sizeof(tmp), "%s.tmp", path);
    TermHistory* h = term_history_open(path, 10);
    const char* in[] = { "e1", "e2", "e3", "e4", "e5" };
    for (int i = 0; i < 5; i++) ASSERT(term_history_add(h, in[i]));
    term_history_close(h);
    h = term_history_open(path, 3);
    ASSERT(term_history_count(h) == 3);
    ASSERT_STR_EQ(term_history_get(h, 0), "e3");
    term_history_close(h);
    char* text = slurp(path);
    ASSERT_STR_EQ(text, "#mathilda-history v1\ne3\ne4\ne5\n");
    free(text);
    ASSERT(access(tmp, F_OK) != 0);
    /* In memory too, the store holds at most max. */
    h = term_history_open(path, 3);
    ASSERT(term_history_add(h, "e6"));
    ASSERT(term_history_count(h) == 3);
    ASSERT_STR_EQ(term_history_get(h, 0), "e4");
    term_history_close(h);
    remove(path);
}

static void test_foreign_and_memory(void) {
    /* 8. A foreign file is never written to. */
    char path[256];
    path_in_dir(path, sizeof(path), "h3");
    FILE* f = fopen(path, "w");
    fputs("hi\n", f);
    fclose(f);
    int saved = test_stderr_redirect("/dev/null");   /* the one note */
    TermHistory* h = term_history_open(path, 100);
    test_stderr_restore(saved);
    ASSERT(h && !term_history_persistent(h));
    ASSERT(term_history_count(h) == 0);
    ASSERT(term_history_add(h, "x"));               /* kept for the session */
    term_history_close(h);
    char* text = slurp(path);
    ASSERT_STR_EQ(text, "hi\n");
    free(text);
    remove(path);

    /* An empty existing file is adopted. */
    f = fopen(path, "w");
    fclose(f);
    h = term_history_open(path, 100);
    ASSERT(term_history_persistent(h));
    term_history_close(h);
    text = slurp(path);
    ASSERT_STR_EQ(text, "#mathilda-history v1\n");
    free(text);
    remove(path);

    /* A pathless store filters in memory only. */
    h = term_history_open(NULL, 2);
    ASSERT(h && !term_history_persistent(h));
    ASSERT(term_history_add(h, "a"));
    ASSERT(!term_history_add(h, "a"));
    ASSERT(term_history_add(h, "b"));
    ASSERT(term_history_add(h, "c"));
    ASSERT(term_history_count(h) == 2);
    term_history_close(h);
    term_history_close(NULL);
}

int main(void) {
    snprintf(g_dir, sizeof(g_dir), "/tmp/mth_hist_test_%ld", (long)getpid());
    mkdir(g_dir, 0700);
    TEST(test_encoding);
    TEST(test_env_rules);
    TEST(test_store);
    TEST(test_truncation);
    TEST(test_foreign_and_memory);
    rmdir(g_dir);
    printf("All term_history tests passed.\n");
    return 0;
}
