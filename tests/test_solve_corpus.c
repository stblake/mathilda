/* test_solve_corpus.c
 *
 * Form-invariant stress-corpus runner for Solve[].
 *
 * Approach mirrors test_intrat_corpus.c: load the verifier prelude
 * (tests/solve_check_prelude.m) once in the parent via Get[], read the
 * corpus (tests/solve_corpus.m) as a HELD List literal so each case's
 * equation stays unevaluated until its own forked child evaluates it,
 * then fork-per-case with a wall-clock alarm.  Each child evaluates
 *
 *     solveCheckCode[label, eqn, vars, dom, expected]
 *
 * which back-substitutes every returned solution into the equation(s)
 * and returns an integer verdict: 0 = PASS, 1 = FAIL (wrong count or a
 * residual that does not vanish or an out-of-domain value), 2 =
 * UNEVALUATED (Solve bubbled back).  The child maps the verdict to a
 * single pipe byte; the parent tallies.
 *
 * The corpus is a progress dashboard, not an all-green assertion: the
 * runner fails only when the number of non-PASS cases exceeds
 * SOLVE_FAIL_BASELINE, the checked-in high-water mark of known-missing
 * capabilities (the Modulus / Rationals / single-eq-multivar /
 * poly-in-kernel gaps).  Each landed fix must lower the baseline; any
 * NEW regression (a case that newly fails, crashes, or times out) trips
 * the test.  See tasks/todo.md and the plan for the phase list.
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

/* Per-case wall-clock cap.  Solve + the back-substitution verifier
 * (which may call Simplify / N over Root objects) is heavier than a
 * single Integrate, so allow more headroom than the rational-integral
 * corpus while still bounding a runaway case. */
#define SOLVE_PER_CASE_TIMEOUT_SEC 15

/* Single-byte verdict codes returned from each child over the pipe. */
enum {
    SOLVE_CODE_PASS   = 'P',  /* verdict 0: verified                     */
    SOLVE_CODE_FAIL   = 'F',  /* verdict 1: wrong count / residual / dom */
    SOLVE_CODE_UNEVAL = 'U',  /* verdict 2: Solve bubbled back           */
    SOLVE_CODE_BADVAL = 'X',  /* solveCheckCode returned a non-integer   */
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
    /* The default 60 s test_utils alarm is per-binary; disable it and
     * let each forked case carry its own alarm instead. */
    alarm(0);
    setvbuf(stderr, NULL, _IONBF, 0);

    const char* corpus_file  = (argc > 1) ? argv[1] : "../solve_corpus.m";
    const char* prelude_file = (argc > 2) ? argv[2] : "../solve_check_prelude.m";

    symtab_init();
    core_init();

    /* Load the verifier prelude in the PARENT so its definitions
     * (solveCheckCode, solveVerdict, ...) are in the symbol table before
     * we fork; children inherit them copy-on-write. */
    fprintf(stderr, "==> prelude: %s\n", prelude_file);
    Expr* get_call = expr_new_function(expr_new_symbol("Get"),
        (Expr*[]){ expr_new_string(prelude_file) }, 1);
    Expr* get_res = evaluate(get_call);
    expr_free(get_call);
    if (get_res) expr_free(get_res);

    /* Read the corpus WITHOUT evaluating it, so degenerate equations
     * (0==0, 1==0, x==x) and every other case reach Solve exactly as
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
    fprintf(stderr, "==> Solve corpus: %zu cases\n", n);

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
            char code = SOLVE_CODE_UNEVAL;
            alarm(SOLVE_PER_CASE_TIMEOUT_SEC);

            Expr* call = expr_new_function(expr_new_symbol("solveCheckCode"),
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
                    case 0:  code = SOLVE_CODE_PASS;   break;
                    case 1:  code = SOLVE_CODE_FAIL;   break;
                    case 2:  code = SOLVE_CODE_UNEVAL; break;
                    default: code = SOLVE_CODE_BADVAL; break;
                }
            } else {
                code = SOLVE_CODE_BADVAL;
            }
            alarm(0);
            if (verdict) expr_free(verdict);

            ssize_t w = write(pipefd[1], &code, 1);
            (void)w;
            close(pipefd[1]);
            _exit(0);
        }

        /* --- Parent --- */
        close(pipefd[1]);
        time_t deadline = time(NULL) + SOLVE_PER_CASE_TIMEOUT_SEC;
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
            fprintf(stderr, "      -> TIMEOUT (>%ds)\n", SOLVE_PER_CASE_TIMEOUT_SEC);
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
            case SOLVE_CODE_PASS:   passed++;      break;
            case SOLVE_CODE_FAIL:   failed++;      break;
            case SOLVE_CODE_UNEVAL: unevaluated++; break;
            default:                crashed++;
                fprintf(stderr, "      -> bad verdict 0x%02x\n",
                        (unsigned char)code);
                break;
        }
    }
    expr_free(cases);

    int regressions = failed + unevaluated + timed_out + crashed;

    fprintf(stderr, "\n=== Solve corpus result ===\n");
    fprintf(stderr, "  Total cases:   %zu\n", n);
    fprintf(stderr, "  Passed:        %d\n", passed);
    fprintf(stderr, "  Failed:        %d\n", failed);
    fprintf(stderr, "  Unevaluated:   %d\n", unevaluated);
    fprintf(stderr, "  Timed out:     %d\n", timed_out);
    fprintf(stderr, "  Crashed:       %d\n", crashed);
    fprintf(stderr, "  Malformed:     %d\n", malformed);
    fprintf(stderr, "  Non-PASS:      %d\n", regressions);
    fprintf(stderr, "===========================\n");

    /* Regression baseline: the checked-in high-water mark of non-PASS
     * cases.  The corpus landed at 10 empirically confirmed gaps --
     * Modulus (4), Rationals (2), single-equation multi-variable (2),
     * polynomial-in-transcendental-kernel (2).  Each phase that lands a
     * fix LOWERS this number; a case that newly fails/crashes/times out
     * pushes it above baseline and fails the test.  Drive it to 0.
     *   10 -> 8 : Phase 2 wired the Rationals domain (D3).
     *    8 -> 4 : Phase 4 wired the Modulus option (D1, 4 cases).
     *    4 -> 2 : Phase 6 wired single-equation multi-variable (D4, 2).
     *    2 -> 0 : Phase 7 wired polynomial-in-transcendental-kernel (D2).
     * The corpus is fully green; any new non-PASS is a real regression. */
    const int SOLVE_FAIL_BASELINE = 0;
    if (regressions > SOLVE_FAIL_BASELINE) {
        fprintf(stderr,
            "FAIL: %d non-PASS case(s) exceed the baseline of %d "
            "(new Solve regression).\n", regressions, SOLVE_FAIL_BASELINE);
        return 1;
    }
    if (malformed > 0) {
        fprintf(stderr, "FAIL: %d malformed corpus record(s).\n", malformed);
        return 1;
    }

    fprintf(stderr, "Solve corpus within baseline (%d/%d non-PASS <= %d).\n",
            regressions, (int)n, SOLVE_FAIL_BASELINE);
    return 0;
}
