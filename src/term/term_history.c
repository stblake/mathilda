/*
 * term_history.c -- see term_history.h.
 */
/* open()/fdopen() for the 0600 create are POSIX, hidden by glibc under -std=c99. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "term_history.h"

#include "dialog.h"   /* mth_is_blank_or_comment_only */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#define TERM_HISTORY_HEADER "#mathilda-history v1"

struct TermHistory {
    char*  path;         /* NULL: in-memory only                         */
    size_t max;
    char** items;        /* ring of the newest entries, oldest first     */
    size_t count, cap;
    bool   writable;     /* entries are appended to `path`               */
};

static char* dup_str(const char* s) {
    size_t n = strlen(s);
    char* r = malloc(n + 1);
    if (r) memcpy(r, s, n + 1);
    return r;
}

/* ---- environment ----------------------------------------------------------- */

char* term_history_resolve_path(const char* env_histfile, const char* env_home) {
    if (env_histfile) return *env_histfile ? dup_str(env_histfile) : NULL;
    if (!env_home || !*env_home) return NULL;
    static const char leaf[] = "/.mathilda_history";
    size_t n = strlen(env_home);
    char* r = malloc(n + sizeof(leaf));
    if (!r) return NULL;
    memcpy(r, env_home, n);
    memcpy(r + n, leaf, sizeof(leaf));
    return r;
}

size_t term_history_size_from_env(const char* v) {
    if (!v || !*v) return TERM_HISTORY_DEFAULT_SIZE;
    size_t n = 0;
    for (const char* p = v; *p; p++) {
        if (*p < '0' || *p > '9') return TERM_HISTORY_DEFAULT_SIZE;
        if (n <= TERM_HISTORY_MAX_SIZE) n = n * 10 + (size_t)(*p - '0');
    }
    return n > TERM_HISTORY_MAX_SIZE ? TERM_HISTORY_MAX_SIZE : n;
}

/* ---- encoding ---------------------------------------------------------------- */

char* term_history_encode(const char* entry) {
    size_t n = 0;
    for (const char* p = entry; *p; p++)
        n += (*p == '\\' || *p == '\n' || *p == '\r') ? 2 : 1;
    char* r = malloc(n + 1);
    if (!r) return NULL;
    char* o = r;
    for (const char* p = entry; *p; p++) {
        if (*p == '\\')      { *o++ = '\\'; *o++ = '\\'; }
        else if (*p == '\n') { *o++ = '\\'; *o++ = 'n'; }
        else if (*p == '\r') { *o++ = '\\'; *o++ = 'r'; }
        else *o++ = *p;
    }
    *o = '\0';
    return r;
}

char* term_history_decode(const char* line) {
    char* r = malloc(strlen(line) + 1);
    if (!r) return NULL;
    char* o = r;
    for (const char* p = line; *p; p++) {
        if (*p == '\\' && p[1]) {
            char c = p[1];
            if (c == '\\')     { *o++ = '\\'; p++; }
            else if (c == 'n') { *o++ = '\n'; p++; }
            else if (c == 'r') { *o++ = '\r'; p++; }
            else *o++ = '\\';          /* unknown escape: kept literally */
        } else {
            *o++ = *p;
        }
    }
    *o = '\0';
    return r;
}

/* ---- the in-memory list ------------------------------------------------------ */

static void items_push(TermHistory* h, char* owned) {
    if (h->count == h->max) {            /* drop the oldest */
        free(h->items[0]);
        memmove(h->items, h->items + 1, (h->count - 1) * sizeof(char*));
        h->count--;
    }
    if (h->count == h->cap) {
        size_t nc = h->cap ? h->cap * 2 : 16;
        if (nc > h->max) nc = h->max;
        char** ni = realloc(h->items, nc * sizeof(char*));
        if (!ni) { free(owned); return; }
        h->items = ni;
        h->cap = nc;
    }
    h->items[h->count++] = owned;
}

/* ---- the file ------------------------------------------------------------------ */

/* Read one line of any length (no getline in C99). Returns malloc'd text
 * without the trailing LF / CRLF, or NULL at end of file. */
static char* read_line(FILE* f) {
    size_t cap = 256, n = 0;
    char* buf = malloc(cap);
    if (!buf) return NULL;
    for (;;) {
        if (!fgets(buf + n, (int)(cap - n), f)) {
            if (n == 0) { free(buf); return NULL; }
            break;
        }
        n += strlen(buf + n);
        if (n > 0 && buf[n - 1] == '\n') { buf[--n] = '\0'; break; }
        if (cap - n < 2) {
            char* nb = realloc(buf, cap * 2);
            if (!nb) { free(buf); return NULL; }
            buf = nb;
            cap *= 2;
        }
    }
    if (n > 0 && buf[n - 1] == '\r') buf[--n] = '\0';
    return buf;
}

/* Open `path` for appending, creating it 0600 when absent. */
static FILE* open_append(const char* path) {
#ifndef _WIN32
    int fd = open(path, O_WRONLY | O_APPEND | O_CREAT, 0600);
    if (fd < 0) return NULL;
    FILE* f = fdopen(fd, "a");
    if (!f) close(fd);
    return f;
#else
    return fopen(path, "a");
#endif
}

static void note(const char* what, const char* path) {
    fprintf(stderr, "Mathilda: %s %s; history is kept for this session only.\n", what, path);
}

static bool write_header_file(const char* path) {
    FILE* f = open_append(path);
    if (!f) return false;
    bool ok = fputs(TERM_HISTORY_HEADER "\n", f) >= 0;
    if (fclose(f) != 0) ok = false;
    return ok;
}

/* Rewrite the file with the entries now in memory (the newest `max`). */
static bool rewrite_file(TermHistory* h) {
    size_t n = strlen(h->path);
    char* tmp = malloc(n + 5);
    if (!tmp) return false;
    memcpy(tmp, h->path, n);
    memcpy(tmp + n, ".tmp", 5);
    remove(tmp);
    FILE* f = open_append(tmp);
    bool ok = f != NULL;
    if (ok) ok = fputs(TERM_HISTORY_HEADER "\n", f) >= 0;
    for (size_t i = 0; ok && i < h->count; i++) {
        char* enc = term_history_encode(h->items[i]);
        ok = enc && fputs(enc, f) >= 0 && fputc('\n', f) != EOF;
        free(enc);
    }
    if (f && fclose(f) != 0) ok = false;
    if (ok) ok = rename(tmp, h->path) == 0;
    if (!ok) remove(tmp);
    free(tmp);
    return ok;
}

TermHistory* term_history_open(const char* path, size_t max) {
    if (max == 0) return NULL;
    if (max > TERM_HISTORY_MAX_SIZE) max = TERM_HISTORY_MAX_SIZE;
    TermHistory* h = calloc(1, sizeof(*h));
    if (!h) return NULL;
    h->max = max;
    if (!path) return h;                      /* in-memory store */
    h->path = dup_str(path);
    if (!h->path) return h;

    FILE* f = fopen(path, "r");
    if (!f) {
        if (errno == ENOENT && write_header_file(path)) h->writable = true;
        else note("cannot create history file", path);
        return h;
    }
    char* first = read_line(f);
    if (!first) {                             /* an empty file: adopt it */
        fclose(f);
        if (write_header_file(path)) h->writable = true;
        else note("cannot write history file", path);
        return h;
    }
    bool ours = strcmp(first, TERM_HISTORY_HEADER) == 0;
    free(first);
    if (!ours) {
        fclose(f);
        note("not a Mathilda history file:", path);
        return h;
    }
    size_t total = 0;
    for (char* line; (line = read_line(f)) != NULL; ) {
        char* e = term_history_decode(line);
        free(line);
        if (e) { items_push(h, e); total++; }
    }
    fclose(f);
    h->writable = true;
    if (total > max && !rewrite_file(h)) {
        /* Still appendable; the file is only longer than it needs to be. */
    }
    return h;
}

bool term_history_add(TermHistory* h, const char* entry) {
    if (!h || !entry) return false;
    if (mth_is_blank_or_comment_only(entry)) return false;
    if (h->count > 0 && strcmp(h->items[h->count - 1], entry) == 0) return false;
    char* copy = dup_str(entry);
    if (!copy) return false;
    items_push(h, copy);
    if (h->writable) {
        char* enc = term_history_encode(entry);
        FILE* f = enc ? open_append(h->path) : NULL;
        bool ok = f != NULL;
        if (ok) ok = fputs(enc, f) >= 0 && fputc('\n', f) != EOF;
        if (f && fclose(f) != 0) ok = false;
        free(enc);
        if (!ok) {
            h->writable = false;
            note("cannot write history file", h->path);
        }
    }
    return true;
}

size_t term_history_count(const TermHistory* h) { return h ? h->count : 0; }

const char* term_history_get(const TermHistory* h, size_t i) {
    return (h && i < h->count) ? h->items[i] : NULL;
}

size_t term_history_max(const TermHistory* h) { return h ? h->max : 0; }

bool term_history_persistent(const TermHistory* h) { return h && h->writable; }

void term_history_close(TermHistory* h) {
    if (!h) return;
    for (size_t i = 0; i < h->count; i++) free(h->items[i]);
    free(h->items);
    free(h->path);
    free(h);
}
