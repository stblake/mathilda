/*
 * term_config.c -- see term_config.h.
 */
/* mkdir is POSIX, hidden by glibc under -std=c99. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "term_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

static char* dup_str(const char* s) {
    size_t n = strlen(s);
    char* p = malloc(n + 1);
    if (p) memcpy(p, s, n + 1);
    return p;
}

char* term_config_dir(const char* env_dir, const char* env_home) {
    if (env_dir) return *env_dir ? dup_str(env_dir) : NULL;
    if (!env_home || !*env_home) return NULL;
    return term_config_join(env_home, ".mathilda");
}

char* term_config_join(const char* dir, const char* name) {
    if (!dir) return NULL;
    size_t a = strlen(dir), b = strlen(name);
    while (a > 1 && dir[a - 1] == '/') a--;
    char* p = malloc(a + 1 + b + 1);
    if (!p) return NULL;
    memcpy(p, dir, a);
    p[a] = '/';
    memcpy(p + a + 1, name, b + 1);
    return p;
}

static bool is_space(char c) { return c == ' ' || c == '\t' || c == '\r'; }

/* Does line [s, e) set `key`? On a match *val and *vlen get the trimmed value. */
static bool line_key(const char* s, const char* e, const char* key,
                     const char** val, size_t* vlen) {
    while (s < e && is_space(*s)) s++;
    size_t k = strlen(key);
    if ((size_t)(e - s) < k || memcmp(s, key, k) != 0) return false;
    s += k;
    while (s < e && is_space(*s)) s++;
    if (s >= e || *s != '=') return false;
    s++;
    while (s < e && is_space(*s)) s++;
    const char* v = s;
    for (const char* h = v; h < e; h++)
        if (*h == '#') { e = h; break; }
    while (e > v && is_space(e[-1])) e--;
    *val = v;
    *vlen = (size_t)(e - v);
    return true;
}

static char* read_file(const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) return NULL;
    size_t cap = 1024, n = 0;
    char* b = malloc(cap);
    while (b) {
        size_t r = fread(b + n, 1, cap - n - 1, f);
        n += r;
        if (n + 1 < cap) break;
        if (cap >= (1u << 20)) break;              /* a config file is small */
        char* nb = realloc(b, cap * 2);
        if (!nb) { free(b); b = NULL; break; }
        b = nb;
        cap *= 2;
    }
    fclose(f);
    if (b) b[n] = '\0';
    return b;
}

bool term_config_get(const char* path, const char* key, char* buf, size_t cap) {
    if (!path || !cap) return false;
    char* t = read_file(path);
    if (!t) return false;
    bool found = false;
    const char* s = t;
    while (*s && !found) {
        const char* e = strchr(s, '\n');
        if (!e) e = s + strlen(s);
        const char* v;
        size_t vl;
        if (line_key(s, e, key, &v, &vl)) {
            if (vl >= cap) vl = cap - 1;
            memcpy(buf, v, vl);
            buf[vl] = '\0';
            found = true;
        }
        s = *e ? e + 1 : e;
    }
    free(t);
    return found;
}

char* term_config_edit(const char* text, const char* key, const char* value) {
    if (!text) text = "";
    size_t cap = strlen(text) + strlen(key) + strlen(value) + 8;
    char* out = malloc(cap);
    if (!out) return NULL;
    size_t n = 0;
    bool done = false;
    const char* s = text;
    while (*s) {
        const char* e = strchr(s, '\n');
        const char* next = e ? e + 1 : s + strlen(s);
        if (!e) e = next;
        const char* v;
        size_t vl;
        if (line_key(s, e, key, &v, &vl)) {
            if (!done) {
                n += (size_t)sprintf(out + n, "%s = %s\n", key, value);
                done = true;
            }
        } else {
            memcpy(out + n, s, (size_t)(e - s));
            n += (size_t)(e - s);
            out[n++] = '\n';
        }
        s = next;
    }
    if (!done) n += (size_t)sprintf(out + n, "%s = %s\n", key, value);
    out[n] = '\0';
    return out;
}

static void make_parent_dir(const char* path) {
    const char* slash = strrchr(path, '/');
    if (!slash || slash == path) return;
    size_t n = (size_t)(slash - path);
    char* d = malloc(n + 1);
    if (!d) return;
    memcpy(d, path, n);
    d[n] = '\0';
#ifdef _WIN32
    (void)_mkdir(d);
#else
    (void)mkdir(d, 0700);
#endif
    free(d);
}

bool term_config_set(const char* path, const char* key, const char* value) {
    if (!path) return false;
    char* old = read_file(path);
    char* text = term_config_edit(old, key, value);
    free(old);
    if (!text) return false;
    make_parent_dir(path);
    size_t pl = strlen(path);
    char* tmp = malloc(pl + 5);
    if (!tmp) { free(text); return false; }
    memcpy(tmp, path, pl);
    memcpy(tmp + pl, ".tmp", 5);
    bool ok = false;
    FILE* f = fopen(tmp, "wb");
    if (f) {
        size_t n = strlen(text);
        ok = fwrite(text, 1, n, f) == n;
        ok = (fclose(f) == 0) && ok;
        if (ok) {
#ifdef _WIN32
            remove(path);                          /* rename does not replace */
#endif
            ok = rename(tmp, path) == 0;
        }
        if (!ok) remove(tmp);
    }
    free(tmp);
    free(text);
    return ok;
}
