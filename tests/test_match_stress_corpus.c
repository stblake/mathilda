/* test_match_stress_corpus.c
 *
 * Conformance + robustness runner for the pattern matcher.
 *
 * Loads a corpus `.m` file: a List of {inputString, expected} pairs, where
 * `inputString` is the exact text a user would type at the REPL (a String,
 * so it is NOT evaluated when Get[] builds the list) and `expected` is the
 * documented Wolfram-Language result (a bare expression, so Get[] evaluates
 * it to its canonical form). For every pair the runner:
 *
 *   1. parses inputString to an Expr
 *   2. evaluates it
 *   3. checks the result is structurally equal (expr_eq) to `expected`,
 *      after normalising any packed NDArray back to a plain List so the
 *      comparison is on value, not on a storage decision.
 *
 * Each case runs in a forked child with a per-case wall-clock timeout, the
 * same isolation test_crc_corpus.c uses: a pattern that infinite-loops
 * (adversarial //. / backtracking) or overflows the C stack (deep nesting,
 * the matcher has no depth guard) is reported as TIMEOUT / CRASH for that
 * one case instead of taking down the whole run.
 *
 * Two modes, selected by argv[2]:
 *   (default)   assert mode  -- nonzero exit if any FAIL / TIMEOUT / CRASH.
 *   --observe   report mode  -- prints MATCH / DIVERGE per case and always
 *                              exits 0. Used for the ordering-sensitive /
 *                              genuinely-subtle cases whose "expected" value
 *                              is documented-from-memory and wants a human
 *                              eye rather than a hard gate.
 *
 * argv[1] overrides the corpus path (default ../match_stress_corpus.m, the
 * layout used when running from tests/build).
 */

#include "test_utils.h"   /* parse_expression, evaluate, expr_to_string,
                             expr_free, test_delist, ASSERT */
#include "expr.h"

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

#define CORPUS_PER_CASE_TIMEOUT_SEC 10

enum {
    CASE_PASS    = 'P',   /* result matched expected                     */
    CASE_FAIL    = 'F',   /* result differed from expected               */
    CASE_MALFORM = 'M',   /* pair not {String, expr} shaped              */
};

static int is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, "List") == 0;
}

/* Structural equality that ignores packed-vs-nested storage. Operates on
 * copies so the shared corpus tree (inherited across fork) is never mutated. */
static bool value_eq(Expr* result, Expr* expected) {
    Expr* r = test_delist(expr_copy(result));
    Expr* e = test_delist(expr_copy(expected));
    bool ok = expr_eq(r, e);
    expr_free(r);
    expr_free(e);
    return ok;
}

int main(int argc, char** argv) {
    alarm(0);                       /* cancel test_utils.h's 60s constructor alarm */
    setvbuf(stderr, NULL, _IONBF, 0);

    const char* corpus_file =
        (argc > 1) ? argv[1] : "../match_stress_corpus.m";
    bool observe = (argc > 2) && strcmp(argv[2], "--observe") == 0;

    symtab_init();
    core_init();

    fprintf(stderr, "==> match stress corpus: %s  (%s mode)\n",
            corpus_file, observe ? "observe" : "assert");

    char get_buf[1024];
    snprintf(get_buf, sizeof(get_buf), "Get[\"%s\"]", corpus_file);
    Expr* get_call = parse_expression(get_buf);
    ASSERT(get_call != NULL);
    Expr* tests = evaluate(get_call);
    expr_free(get_call);

    if (!is_list(tests)) {
        fprintf(stderr, "FAIL: could not load %s as a List.\n", corpus_file);
        if (tests) expr_free(tests);
        ASSERT(false);
        return 1;
    }

    size_t n = tests->data.function.arg_count;
    fprintf(stderr, "==> %zu cases\n", n);

    int passed = 0, failed = 0, malformed = 0, timed_out = 0, crashed = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* pair = tests->data.function.args[i];
        Expr* input    = is_list(pair) && pair->data.function.arg_count >= 2
                       ? pair->data.function.args[0] : NULL;
        Expr* expected = input ? pair->data.function.args[1] : NULL;
        if (!input || input->type != EXPR_STRING) {
            malformed++;
            fprintf(stderr, "  [%3zu/%zu] MALFORMED (input must be a String)\n",
                    i + 1, n);
            continue;
        }

        int pipefd[2];
        if (pipe(pipefd) != 0) { crashed++; continue; }
        pid_t pid = fork();
        if (pid < 0) { close(pipefd[0]); close(pipefd[1]); crashed++; continue; }

        if (pid == 0) {
            /* --- Child --- */
            close(pipefd[0]);
            alarm(CORPUS_PER_CASE_TIMEOUT_SEC);

            char code = CASE_FAIL;
            Expr* parsed = parse_expression(input->data.string);
            if (!parsed) {
                fprintf(stderr, "  [%3zu/%zu] PARSE ERROR: %s\n",
                        i + 1, n, input->data.string);
            } else {
                Expr* result = evaluate(parsed);
                expr_free(parsed);
                bool ok = value_eq(result, expected);
                code = ok ? CASE_PASS : CASE_FAIL;
                if (observe || !ok) {
                    char* ex = expr_to_string(expected);
                    char* rs = expr_to_string(result);
                    fprintf(stderr, "  [%3zu/%zu] %s  %s\n"
                                    "        expected: %s\n"
                                    "        actual:   %s\n",
                            i + 1, n, ok ? "MATCH  " : "DIVERGE",
                            input->data.string, ex, rs);
                    free(ex); free(rs);
                }
                expr_free(result);
            }
            alarm(0);
            ssize_t w = write(pipefd[1], &code, 1); (void)w;
            close(pipefd[1]);
            _exit(0);
        }

        /* --- Parent --- */
        close(pipefd[1]);
        time_t deadline = time(NULL) + CORPUS_PER_CASE_TIMEOUT_SEC + 1;
        int status = 0;
        bool reaped = false;
        while (time(NULL) < deadline) {
            pid_t r = waitpid(pid, &status, WNOHANG);
            if (r == pid) { reaped = true; break; }
            if (r < 0) break;
            usleep(20000);
        }
        if (!reaped) {
            kill(pid, SIGKILL);
            waitpid(pid, &status, 0);
            close(pipefd[0]);
            timed_out++;
            fprintf(stderr, "  [%3zu/%zu] TIMEOUT (>%ds): %s\n",
                    i + 1, n, CORPUS_PER_CASE_TIMEOUT_SEC, input->data.string);
            continue;
        }

        char code = 0;
        ssize_t r = read(pipefd[0], &code, 1);
        close(pipefd[0]);

        if (r != 1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
            if (WIFSIGNALED(status) && WTERMSIG(status) == SIGALRM) {
                timed_out++;
                fprintf(stderr, "  [%3zu/%zu] TIMEOUT (child alarm): %s\n",
                        i + 1, n, input->data.string);
            } else {
                crashed++;
                fprintf(stderr, "  [%3zu/%zu] CRASH (exit %d, sig %d): %s\n",
                        i + 1, n,
                        WIFEXITED(status) ? WEXITSTATUS(status) : -1,
                        WIFSIGNALED(status) ? WTERMSIG(status) : 0,
                        input->data.string);
            }
            continue;
        }

        switch (code) {
            case CASE_PASS:    passed++;    break;
            case CASE_MALFORM: malformed++; break;
            default:           failed++;    break;
        }
    }
    expr_free(tests);

    fprintf(stderr, "\n=== match stress corpus result ===\n");
    fprintf(stderr, "  Total cases:  %zu\n", n);
    fprintf(stderr, "  Passed:       %d\n", passed);
    fprintf(stderr, "  Failed:       %d\n", failed);
    fprintf(stderr, "  Timed out:    %d\n", timed_out);
    fprintf(stderr, "  Crashed:      %d\n", crashed);
    fprintf(stderr, "  Malformed:    %d\n", malformed);
    fprintf(stderr, "==================================\n");

    if (observe) return 0;   /* report-only tier never gates */
    return (failed + timed_out + crashed + malformed) == 0 ? 0 : 1;
}
