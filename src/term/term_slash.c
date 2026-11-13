/*
 * term_slash.c -- see term_slash.h.
 */
#include "term_slash.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#define NAME_MAX_LEN 32

static TermSlashCmd* g_cmds = NULL;
static size_t        g_n = 0, g_cap = 0;

static char* dup_n(const char* s, size_t n) {
    char* p = malloc(n + 1);
    if (!p) return NULL;
    memcpy(p, s, n);
    p[n] = '\0';
    return p;
}

bool term_slash_valid_name(const char* s, size_t n) {
    if (!s || n == 0 || n > NAME_MAX_LEN || !isalpha((unsigned char)s[0])) return false;
    for (size_t i = 1; i < n; i++)
        if (!isalnum((unsigned char)s[i])) return false;
    return true;
}

static bool name_eq(const char* a, const char* b, size_t n) {
    return strlen(a) == n && memcmp(a, b, n) == 0;
}

const TermSlashCmd* term_slash_find(const char* name, size_t n) {
    if (n && name[0] == '/') { name++; n--; }
    for (size_t i = 0; i < g_n; i++)
        if (name_eq(g_cmds[i].name, name, n)) return &g_cmds[i];
    return NULL;
}

static void cmd_free(TermSlashCmd* c) {
    if (c->payload && c->payload_free) c->payload_free(c->payload);
    free(c->name);
    free(c->summary);
    free(c->owner);
    memset(c, 0, sizeof(*c));
}

/* Insert at position `at`, shifting the rest. */
static TermSlashCmd* insert_at(size_t at) {
    if (g_n == g_cap) {
        size_t nc = g_cap ? g_cap * 2 : 16;
        TermSlashCmd* nv = realloc(g_cmds, nc * sizeof(*nv));
        if (!nv) return NULL;
        g_cmds = nv;
        g_cap = nc;
    }
    memmove(&g_cmds[at + 1], &g_cmds[at], (g_n - at) * sizeof(*g_cmds));
    g_n++;
    memset(&g_cmds[at], 0, sizeof(*g_cmds));
    return &g_cmds[at];
}

static bool fill(TermSlashCmd* c, const char* name, const char* summary, const char* owner) {
    c->name = dup_n(name, strlen(name));
    c->summary = dup_n(summary ? summary : "", summary ? strlen(summary) : 0);
    c->owner = owner ? dup_n(owner, strlen(owner)) : NULL;
    return c->name && c->summary && (!owner || c->owner);
}

bool term_slash_add_builtin(const char* name, const char* summary, int id) {
    if (!term_slash_valid_name(name, strlen(name)) || term_slash_find(name, strlen(name))) return false;
    size_t at = 0;
    while (at < g_n && g_cmds[at].kind == TERM_SLASH_BUILTIN) at++;
    TermSlashCmd* c = insert_at(at);
    if (!c) return false;
    if (!fill(c, name, summary, NULL)) {
        cmd_free(c);
        memmove(&g_cmds[at], &g_cmds[at + 1], (g_n - at - 1) * sizeof(*g_cmds));
        g_n--;
        return false;
    }
    c->kind = TERM_SLASH_BUILTIN;
    c->id = id;
    return true;
}

bool term_slash_add_mod(const char* name, const char* summary, const char* owner,
                        void* payload, void (*payload_free)(void*)) {
    size_t n = name ? strlen(name) : 0;
    if (!term_slash_valid_name(name, n)) return false;
    const TermSlashCmd* old = term_slash_find(name, n);
    if (old && old->kind == TERM_SLASH_BUILTIN) return false;
    TermSlashCmd tmp;
    memset(&tmp, 0, sizeof tmp);
    if (!fill(&tmp, name, summary, owner)) { tmp.payload = NULL; cmd_free(&tmp); return false; }
    tmp.kind = TERM_SLASH_MOD;
    tmp.payload = payload;
    tmp.payload_free = payload_free;
    if (old) {
        TermSlashCmd* c = &g_cmds[old - g_cmds];
        cmd_free(c);
        *c = tmp;
        return true;
    }
    size_t at = 0;
    while (at < g_n && g_cmds[at].kind == TERM_SLASH_BUILTIN) at++;
    while (at < g_n && strcmp(g_cmds[at].name, name) < 0) at++;
    TermSlashCmd* c = insert_at(at);
    if (!c) { tmp.payload = NULL; cmd_free(&tmp); return false; }
    *c = tmp;
    return true;
}

bool term_slash_remove(const char* name) {
    const TermSlashCmd* c = term_slash_find(name, strlen(name));
    if (!c || c->kind != TERM_SLASH_MOD) return false;
    size_t at = (size_t)(c - g_cmds);
    cmd_free(&g_cmds[at]);
    memmove(&g_cmds[at], &g_cmds[at + 1], (g_n - at - 1) * sizeof(*g_cmds));
    g_n--;
    return true;
}

void term_slash_clear(void) {
    for (size_t i = 0; i < g_n; i++) cmd_free(&g_cmds[i]);
    free(g_cmds);
    g_cmds = NULL;
    g_n = g_cap = 0;
}

size_t term_slash_count(void) { return g_n; }
const TermSlashCmd* term_slash_at(size_t i) { return i < g_n ? &g_cmds[i] : NULL; }

static bool prefix_ci(const char* name, const char* p, size_t n) {
    for (size_t i = 0; i < n; i++) {
        if (!name[i]) return false;
        if (tolower((unsigned char)name[i]) != tolower((unsigned char)p[i])) return false;
    }
    return true;
}

void term_slash_complete(const char* buf, size_t cursor, TermCompletions* out) {
    out->items = NULL;
    out->count = 0;
    out->word_start = 0;
    if (!buf || cursor == 0 || buf[0] != '/') return;
    for (size_t i = 1; i < cursor; i++)
        if (!isalnum((unsigned char)buf[i])) return;
    if (g_n == 0) return;
    out->items = malloc(g_n * sizeof(char*));
    if (!out->items) return;
    for (size_t i = 0; i < g_n; i++) {
        if (!prefix_ci(g_cmds[i].name, buf + 1, cursor - 1)) continue;
        size_t ln = strlen(g_cmds[i].name);
        char* it = malloc(ln + 2);
        if (!it) break;
        it[0] = '/';
        memcpy(it + 1, g_cmds[i].name, ln + 1);
        out->items[out->count++] = it;
    }
    if (out->count == 0) { free(out->items); out->items = NULL; }
}

const char* term_slash_summary(const char* item, size_t n) {
    const TermSlashCmd* c = term_slash_find(item, n);
    return c && c->summary[0] ? c->summary : NULL;
}

static bool is_ws(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }

bool term_slash_parse(const char* line, const char** name, size_t* nlen,
                      const char** args, size_t* alen) {
    if (!line || line[0] != '/') return false;
    const char* p = line + 1;
    const char* q = p;
    while (*q && !is_ws(*q)) q++;
    *name = p;
    *nlen = (size_t)(q - p);
    while (*q && is_ws(*q)) q++;
    const char* e = q + strlen(q);
    while (e > q && is_ws(e[-1])) e--;
    *args = q;
    *alen = (size_t)(e - q);
    return true;
}
