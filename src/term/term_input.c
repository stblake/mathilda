/*
 * term_input.c -- see term_input.h. The reader seam and the stdio backend.
 */
#include "term_input.h"
#include "term_interrupt.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TERM_STDIO_MAX_LINE 10240   /* repl.c's MAX_INPUT_LEN, unchanged */

TermReadStatus term_read(TermReader* rd, const char* prompt, char** out) {
    *out = NULL;
    if (!rd || !rd->ops || !rd->ops->read) return TERM_READ_EOF;
    (void)term_interrupt_take_reading();            /* a stale prompt Ctrl-C */
    TermState prev = term_interrupt_set_state(TERM_STATE_READING);
    TermReadStatus st = rd->ops->read(rd, prompt, out);
    term_interrupt_set_state(prev);
    if (st != TERM_READ_OK && *out) { free(*out); *out = NULL; }
    return st;
}

void term_reader_history_add(TermReader* rd, const char* entry) {
    if (rd && rd->ops && rd->ops->history_add && entry) rd->ops->history_add(rd, entry);
}

TermBackendKind term_reader_kind(const TermReader* rd) {
    return rd ? rd->kind : TERM_BACKEND_STDIO;
}

void term_reader_close(TermReader* rd) {
    if (!rd) return;
    if (rd->ops && rd->ops->close) rd->ops->close(rd);
    free(rd);
}

/* ---- stdio backend ------------------------------------------------------- */

static TermReadStatus stdio_read(TermReader* rd, const char* prompt, char** out) {
    (void)rd;
    char line[TERM_STDIO_MAX_LINE];
    printf("%s", prompt);
    fflush(stdout);
    if (!fgets(line, sizeof(line), stdin)) return TERM_READ_EOF;
    size_t len = strlen(line);
    while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = '\0';
    char* copy = malloc(len + 1);
    if (!copy) return TERM_READ_EOF;
    memcpy(copy, line, len + 1);
    *out = copy;
    return TERM_READ_OK;
}

static const TermReaderOps STDIO_OPS = { stdio_read, NULL, NULL };

TermReader* term_reader_open_stdio(void) {
    TermReader* rd = malloc(sizeof(*rd));
    if (!rd) return NULL;
    rd->ops = &STDIO_OPS;
    rd->kind = TERM_BACKEND_STDIO;
    rd->impl = NULL;
    return rd;
}
