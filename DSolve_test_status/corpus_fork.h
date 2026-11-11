/* corpus_fork.h -- the fork-per-case core shared by the DSolve corpus runners
 * (test_dsolve_corpus.c for ODEs, pde/test_dsolve_pde_corpus.c for PDEs).
 *
 * corpus_run_case() evaluates one verifier call in a forked child under a
 * wall-clock alarm, so a crash, a hang, or leaked state in one record cannot
 * touch the next.  The child maps the verifier's integer verdict through the
 * caller's table to ONE pipe byte; the parent reaps it with its own deadline
 * (alarm + 2 s) and SIGKILLs a child that ignores SIGALRM.
 *
 * Header-only (static functions) so each runner stays a single translation
 * unit, as CMake builds it.  POSIX: fork/pipe/waitpid -- the corpus runners are
 * nightly Linux/macOS tools, never part of the Windows build.
 */
#ifndef DSOLVE_CORPUS_FORK_H
#define DSOLVE_CORPUS_FORK_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "expr.h"

typedef enum {
    CORPUS_RAN,       /* the child answered; *code holds its byte        */
    CORPUS_TIMEOUT,   /* the child overran its alarm (either deadline)   */
    CORPUS_CRASH      /* the child died, or fork/pipe failed             */
} CorpusOutcome;

/* Map an integer verdict to its pipe byte; `bad` for anything off the table
 * (including a non-integer verdict). */
static char corpus_verdict_byte(const Expr* v, const char* table, int ntable, char bad) {
    if (!v || v->type != EXPR_INTEGER) return bad;
    if (v->data.integer < 0 || v->data.integer >= ntable) return bad;
    return table[v->data.integer];
}

/* Evaluate `call` (borrowed) in a forked child.  table[k] is the byte for
 * verdict k; `bad` the byte for anything else. */
static CorpusOutcome corpus_run_case(Expr* call, int timeout_sec,
                                     const char* table, int ntable, char bad,
                                     char* code) {
    int pipefd[2];
    if (pipe(pipefd) != 0) return CORPUS_CRASH;
    pid_t pid = fork();
    if (pid < 0) { close(pipefd[0]); close(pipefd[1]); return CORPUS_CRASH; }

    if (pid == 0) {                                   /* --- child --- */
        close(pipefd[0]);
        alarm((unsigned)timeout_sec);
        Expr* verdict = evaluate(call);
        char c = corpus_verdict_byte(verdict, table, ntable, bad);
        alarm(0);
        if (verdict) expr_free(verdict);
        ssize_t w = write(pipefd[1], &c, 1);
        (void)w;
        close(pipefd[1]);
        _exit(0);
    }

    close(pipefd[1]);                                 /* --- parent --- */
    time_t deadline = time(NULL) + timeout_sec + 2;
    int status = 0;
    bool reaped = false;
    while (time(NULL) < deadline) {
        pid_t r = waitpid(pid, &status, WNOHANG);
        if (r == pid) { reaped = true; break; }
        if (r < 0) break;
        usleep(50000);
    }
    if (!reaped) {
        kill(pid, SIGKILL);
        waitpid(pid, &status, 0);
        close(pipefd[0]);
        return CORPUS_TIMEOUT;
    }
    ssize_t r = read(pipefd[0], code, 1);
    close(pipefd[0]);
    if (r != 1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        if (WIFSIGNALED(status) && WTERMSIG(status) == SIGALRM) return CORPUS_TIMEOUT;
        fprintf(stderr, "      -> CRASH (exit %d, sig %d)\n",
                WIFEXITED(status) ? WEXITSTATUS(status) : -1,
                WIFSIGNALED(status) ? WTERMSIG(status) : 0);
        return CORPUS_CRASH;
    }
    return CORPUS_RAN;
}

/* Read a whole file into a NUL-terminated buffer (caller frees), or NULL. */
static char* corpus_slurp(const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (sz < 0) { fclose(f); return NULL; }
    char* buf = (char*)malloc((size_t)sz + 1);
    if (!buf) { fclose(f); return NULL; }
    size_t rd = fread(buf, 1, (size_t)sz, f);
    buf[rd] = '\0';
    fclose(f);
    return buf;
}

#endif /* DSOLVE_CORPUS_FORK_H */
