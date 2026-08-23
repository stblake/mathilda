/* test_reduce_corpus.c
 *
 * Semantic stress-corpus runner for Reduce[].
 *
 * Mirrors test_solve_corpus.c: load the verifier prelude
 * (tests/reduce_check_prelude.m) once in the parent via Get[], read the corpus
 * (tests/reduce_corpus.m) as a HELD List literal so each case's statement stays
 * unevaluated until its own forked child evaluates it, then fork-per-case with a
 * wall-clock alarm.  Each child evaluates
 *
 *     reduceCheckCode[label, expr, vars, dom, expected]
 *
 * which -- unlike Solve, whose output is a list of rules -- verifies a LOGICAL
 * FORMULA by SAMPLING: at a bounded rational/integer grid over the variables and
 * free parameters the input and the output must agree wherever both decide, and
 * every fully-determined solution branch of the output must satisfy the input
 * (witness check).  "decline" cases require Reduce to bubble back unevaluated
 * (the soundness invariant).  It returns an integer verdict: 0 = PASS, 1 = FAIL
 * (a wrong constant, an equivalence violation, or a formula where a decline was
 * required), 2 = UNEVALUATED (bubbled back when a result was expected).  The
 * child maps the verdict to a single pipe byte; the parent tallies.
 *
 * The corpus ships fully green, so the baseline is 0: any non-PASS is a real
 * regression.  A genuine equivalence failure (not a mere non-minimal form) is a
 * soundness bug and must be fixed, never baselined away.
 */

#include "test_utils.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>

void symtab_init(void);
void core_init(void);

/* Full Expr definition + constructors (test_utils.h only forward-declares). */
#include "expr.h"

extern char* expr_to_string(struct Expr*);

/* Per-case wall-clock cap.  Reduce + the sampling verifier (grid + witnesses,
 * each substituting into the input and the output formula) is comparable to the
 * Solve corpus's back-substitution, so allow the same headroom. */
#define REDUCE_PER_CASE_TIMEOUT_SEC 15

/* Single-byte verdict codes returned from each child over the pipe. */
enum {
    REDUCE_CODE_PASS   = 'P',  /* verdict 0: verified                        */
    REDUCE_CODE_FAIL   = 'F',  /* verdict 1: wrong constant / equivalence /  */
                               /*            formula where decline required  */
    REDUCE_CODE_UNEVAL = 'U',  /* verdict 2: Reduce bubbled back             */
    REDUCE_CODE_BADVAL = 'X',  /* reduceCheckCode returned a non-integer     */
};

/* Read an entire file into a NUL-terminated heap buffer, or NULL. */
static char* slurp(const char* path) {
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

int main(int argc, char** argv) {
    /* The default 60 s test_utils alarm is per-binary; disable it and let each
     * forked case carry its own alarm instead. */
    alarm(0);
    setvbuf(stderr, NULL, _IONBF, 0);

    const char* corpus_file  = (argc > 1) ? argv[1] : "../reduce_corpus.m";
    const char* prelude_file = (argc > 2) ? argv[2] : "../reduce_check_prelude.m";

    symtab_init();
    core_init();

    /* Load the verifier prelude in the PARENT so its definitions
     * (reduceCheckCode, reduceVerdict, ...) are in the symbol table before we
     * fork; children inherit them copy-on-write. */
    fprintf(stderr, "==> prelude: %s\n", prelude_file);
    Expr* get_call = expr_new_function(expr_new_symbol("Get"),
        (Expr*[]){ expr_new_string(prelude_file) }, 1);
    Expr* get_res = evaluate(get_call);
    expr_free(get_call);
    if (get_res) expr_free(get_res);

    /* Read the corpus WITHOUT evaluating it, so every case (constant relations,
     * degenerate equations, and everything else) reaches Reduce exactly as
     * written inside the forked child rather than collapsing at load. */
    fprintf(stderr, "==> corpus:  %s\n", corpus_file);
    char* fbuf = slurp(corpus_file);
    if (!fbuf) {
        fprintf(stderr, "FAIL: cannot read %s\n", corpus_file);
        ASSERT(false);
        return 1;
    }
    Expr* cases = parse_expression(fbuf);
    free(fbuf);

    if (!cases || cases->type != EXPR_FUNCTION
        || !cases->data.function.head
        || cases->data.function.head->type != EXPR_SYMBOL
        || strcmp(cases->data.function.head->data.symbol.name, "List") != 0) {
        fprintf(stderr, "FAIL: %s did not parse as a List literal.\n", corpus_file);
        if (cases) expr_free(cases);
        ASSERT(false);
        return 1;
    }

    size_t n = cases->data.function.arg_count;
    fprintf(stderr, "==> Reduce corpus: %zu cases\n", n);

    int passed = 0, failed = 0, unevaluated = 0;
    int malformed = 0, timed_out = 0, crashed = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* rec = cases->data.function.args[i];
        if (rec->type != EXPR_FUNCTION
            || !rec->data.function.head
            || rec->data.function.head->type != EXPR_SYMBOL
            || strcmp(rec->data.function.head->data.symbol.name, "List") != 0
            || rec->data.function.arg_count < 5) {
            malformed++;
            continue;
        }
        Expr* label = rec->data.function.args[0];
        char* label_str = expr_to_string(label);
        fprintf(stderr, "  [%3zu/%zu] %s\n", i + 1, n,
                label_str ? label_str : "?");
        free(label_str);

        int pipefd[2];
        if (pipe(pipefd) != 0) { crashed++; continue; }
        pid_t pid = fork();
        if (pid < 0) { close(pipefd[0]); close(pipefd[1]); crashed++; continue; }

        if (pid == 0) {
            /* --- Child --- */
            close(pipefd[0]);
            char code = REDUCE_CODE_UNEVAL;
            alarm(REDUCE_PER_CASE_TIMEOUT_SEC);

            Expr* call = expr_new_function(expr_new_symbol("reduceCheckCode"),
                (Expr*[]){
                    expr_copy(rec->data.function.args[0]),
                    expr_copy(rec->data.function.args[1]),
                    expr_copy(rec->data.function.args[2]),
                    expr_copy(rec->data.function.args[3]),
                    expr_copy(rec->data.function.args[4])
                }, 5);
            Expr* verdict = evaluate(call);
            expr_free(call);

            if (verdict && verdict->type == EXPR_INTEGER) {
                switch (verdict->data.integer) {
                    case 0:  code = REDUCE_CODE_PASS;   break;
                    case 1:  code = REDUCE_CODE_FAIL;   break;
                    case 2:  code = REDUCE_CODE_UNEVAL; break;
                    default: code = REDUCE_CODE_BADVAL; break;
                }
            } else {
                code = REDUCE_CODE_BADVAL;
            }
            alarm(0);
            if (verdict) expr_free(verdict);

            /* reduceCheckCode Prints a "REDUCE-FAIL/UNEVAL <label> <why>"
             * diagnostic on any non-PASS; _exit() would drop the buffered
             * stdio, so flush it before signalling the verdict. */
            fflush(stdout);
            fflush(stderr);

            ssize_t w = write(pipefd[1], &code, 1);
            (void)w;
            close(pipefd[1]);
            _exit(0);
        }

        /* --- Parent --- */
        close(pipefd[1]);
        time_t deadline = time(NULL) + REDUCE_PER_CASE_TIMEOUT_SEC;
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
            timed_out++;
            fprintf(stderr, "      -> TIMEOUT (>%ds)\n", REDUCE_PER_CASE_TIMEOUT_SEC);
            continue;
        }

        char code = 0;
        ssize_t r = read(pipefd[0], &code, 1);
        close(pipefd[0]);

        if (r != 1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
            if (WIFSIGNALED(status) && WTERMSIG(status) == SIGALRM) {
                timed_out++;
                fprintf(stderr, "      -> TIMEOUT (child alarm)\n");
                continue;
            }
            crashed++;
            fprintf(stderr, "      -> CRASH (exit %d, sig %d)\n",
                    WIFEXITED(status) ? WEXITSTATUS(status) : -1,
                    WIFSIGNALED(status) ? WTERMSIG(status)  : 0);
            continue;
        }

        switch (code) {
            case REDUCE_CODE_PASS:   passed++;      break;
            case REDUCE_CODE_FAIL:   failed++;      break;
            case REDUCE_CODE_UNEVAL: unevaluated++; break;
            default:                 crashed++;
                fprintf(stderr, "      -> bad verdict 0x%02x\n",
                        (unsigned char)code);
                break;
        }
    }
    expr_free(cases);

    int regressions = failed + unevaluated + timed_out + crashed;

    fprintf(stderr, "\n=== Reduce corpus result ===\n");
    fprintf(stderr, "  Total cases:   %zu\n", n);
    fprintf(stderr, "  Passed:        %d\n", passed);
    fprintf(stderr, "  Failed:        %d\n", failed);
    fprintf(stderr, "  Unevaluated:   %d\n", unevaluated);
    fprintf(stderr, "  Timed out:     %d\n", timed_out);
    fprintf(stderr, "  Crashed:       %d\n", crashed);
    fprintf(stderr, "  Malformed:     %d\n", malformed);
    fprintf(stderr, "  Non-PASS:      %d\n", regressions);
    fprintf(stderr, "============================\n");

    /* Regression baseline: the corpus ships fully green (every case verified by
     * back-sampling), so this is 0.  A case that newly fails, bubbles back,
     * crashes, or times out pushes it above baseline and fails the test.  A
     * genuine equivalence failure is a Reduce soundness bug -- fix it, never
     * raise this number to hide it. */
    const int REDUCE_FAIL_BASELINE = 0;
    if (regressions > REDUCE_FAIL_BASELINE) {
        fprintf(stderr,
            "FAIL: %d non-PASS case(s) exceed the baseline of %d "
            "(new Reduce regression).\n", regressions, REDUCE_FAIL_BASELINE);
        return 1;
    }
    if (malformed > 0) {
        fprintf(stderr, "FAIL: %d malformed corpus record(s).\n", malformed);
        return 1;
    }

    fprintf(stderr, "Reduce corpus within baseline (%d/%d non-PASS <= %d).\n",
            regressions, (int)n, REDUCE_FAIL_BASELINE);
    return 0;
}
