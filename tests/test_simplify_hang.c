/* test_simplify_hang.c
 *
 * Regression + feature tests for two things that landed together:
 *
 *  1. The Simplify[(x+Sqrt[u]) u^(-1/2)] HANG. The heuristic search finished
 *     instantly (best = 1 + x/Sqrt[u]); the hang was downstream, in Simplify's
 *     radical-fraction polish calling Factor, which spun forever inside the
 *     multivariate GCD's pseudo_rem -- an unguarded while(true) that never
 *     terminates when Factor treats the algebraically-dependent generators
 *     Sqrt[u] and u as independent (the pseudo-remainder degree in the main
 *     variable never drops). Fixes: a degree-monotonicity guard in pseudo_rem
 *     (src/poly/poly.c) and a narrowed polish gate (src/simp/simp_builtins.c,
 *     has_compound_radicand). The bug was latent for EVERY Factor/GCD caller.
 *
 *  2. The new Simplify TimeConstraint option: a synchronous, per-sub-expression
 *     wall-clock budget that fails gracefully (returns the best-so-far form)
 *     with no memory leak.
 *
 * The termination battery runs each previously-dangerous input in a FORKED
 * child under a hard per-case alarm, so a genuine regression reports the exact
 * input that hung instead of freezing the whole suite. Correctness and option
 * checks run in-process (fast; backstopped by the test harness's global alarm).
 * Failures are counted and exit(1) is used so they survive -DNDEBUG builds.
 */

/* fork/waitpid/alarm/_exit/SIGALRM are POSIX; expose them under -std=c99. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "test_utils.h"     /* parse_expression/evaluate/expr_to_string/expr_free + unistd/signal */
#include <sys/wait.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void symtab_init(void);
extern void core_init(void);

static int failures = 0;

/* ---- in-process correctness helpers (fast; global alarm is the backstop) --- */

/* Evaluate `input`, returning the printed result (caller frees) or NULL. */
static char* eval_to_str(const char* input) {
    struct Expr* parsed = parse_expression(input);
    if (!parsed) return NULL;
    struct Expr* res = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(res);
    expr_free(res);
    return s;
}

static void check_eq(const char* input, const char* expected) {
    char* s = eval_to_str(input);
    if (!s) { fprintf(stderr, "FAIL (parse): %s\n", input); failures++; return; }
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
        failures++;
    }
    free(s);
}

/* Assert the result is a real value, not the TimeConstrained/$Aborted sentinel
 * (used for the tiny-budget bail, whose exact form is timing-dependent). */
static void check_not_aborted(const char* input) {
    char* s = eval_to_str(input);
    if (!s) { fprintf(stderr, "FAIL (parse): %s\n", input); failures++; return; }
    if (strstr(s, "$Aborted") != NULL) {
        fprintf(stderr, "FAIL (aborted): %s -> %s\n", input, s);
        failures++;
    }
    free(s);
}

/* ---- termination harness: evaluate in a forked child under a hard alarm ----
 * The child arms alarm(timeout); if the evaluation runs away, SIGALRM (default
 * disposition) terminates the child even mid-CPU-loop. The parent blocks in
 * waitpid and classifies the exit. Returns 0 = terminated cleanly, 3 = parse
 * failure, 100+signal = killed (100+SIGALRM => timed out), -1 = fork failure. */
static int run_terminates(const char* input, unsigned timeout_sec) {
    fflush(stdout);
    fflush(stderr);
    pid_t pid = fork();
    if (pid < 0) return -1;
    if (pid == 0) {
        alarm(timeout_sec);
        struct Expr* p = parse_expression(input);
        if (p) {
            struct Expr* r = evaluate(p);
            expr_free(p);
            if (r) expr_free(r);
        }
        alarm(0);
        _exit(p ? 0 : 3);
    }
    int status = 0;
    waitpid(pid, &status, 0);
    if (WIFEXITED(status))   return (int)WEXITSTATUS(status);
    if (WIFSIGNALED(status)) return 100 + WTERMSIG(status);
    return -2;
}

static void expect_terminates(const char* input, unsigned timeout_sec) {
    int rc = run_terminates(input, timeout_sec);
    if (rc == 0) return;                       /* terminated cleanly */
    if (rc == 100 + SIGALRM)
        fprintf(stderr, "FAIL (HANG >%us): %s\n", timeout_sec, input);
    else if (rc == 3)
        fprintf(stderr, "FAIL (parse): %s\n", input);
    else if (rc == -1)
        fprintf(stderr, "FAIL (fork): %s\n", input);
    else
        fprintf(stderr, "FAIL (rc=%d): %s\n", rc, input);
    failures++;
}

int main(void) {
    symtab_init();
    core_init();
    /* Cancel the inherited global alarm from test_utils.h -- each forked case
     * arms its own bound, and the in-process checks are fast. */
    alarm(0);

    /* ------------------------------------------------------------------ */
    /* Termination battery: 24 radical/rationalisation inputs that route   */
    /* through the Factor/pseudo_rem path. Bases a,b,c,u,v,w sort BEFORE x  */
    /* (the actual trigger); y,z sort after. Each must finish under 5 s.    */
    /* ------------------------------------------------------------------ */
    const char* battery[] = {
        "Simplify[(x+Sqrt[u]) u^(-1/2)]",              /* the reported case */
        "Simplify[(x+Sqrt[a])/Sqrt[a]]",
        "Simplify[(x+Sqrt[b]) b^(-1/2)]",
        "Simplify[(x+Sqrt[c]) c^(-1/2)]",
        "Simplify[(x+Sqrt[w])^2 w^(-1/2)]",
        "Simplify[(2 x + 3 Sqrt[a]) a^(-1/2)]",
        "Simplify[(x + y Sqrt[a]) a^(-1/2)]",
        "Simplify[(x + Sqrt[a] + Sqrt[b]) a^(-1/2)]",
        "Simplify[(x + Sqrt[a b]) (a b)^(-1/2)]",
        "Simplify[1 + x/Sqrt[a]]",
        "Simplify[1 + x/Sqrt[u] + y/Sqrt[u]]",
        "Simplify[(x + a^(1/3)) a^(-1/3)]",
        "Simplify[(x + a^(2/3)) a^(-1/3)]",
        "Simplify[(p + q Sqrt[a]) a^(-1/2)]",
        "Simplify[(x - Sqrt[a]) a^(-1/2)]",
        "Simplify[(x + Sqrt[a])^3 a^(-1/2)]",
        "Simplify[Sqrt[a]/(x + Sqrt[a])]",
        "Simplify[(x + Sqrt[a])/(y + Sqrt[a])]",
        "Simplify[-x + Sqrt[x^2 + a]]",                /* Lagrange radical branch */
        "Simplify[(-x + Sqrt[x^2 + a])/a]",
        "Factor[1 + x/Sqrt[a]]",                       /* direct Factor path */
        "Factor[1 + x/Sqrt[u]]",
        "Factor[(x + Sqrt[a])/Sqrt[a]]",
        "Simplify[(v + Sqrt[u]) u^(-1/2)]",
    };
    for (size_t i = 0; i < sizeof(battery) / sizeof(battery[0]); i++)
        expect_terminates(battery[i], 5);

    /* ------------------------------------------------------------------ */
    /* Exact-output regressions (deterministic forms, verified post-fix).  */
    /* ------------------------------------------------------------------ */
    check_eq("Simplify[(x+Sqrt[u]) u^(-1/2)]", "1 + x/Sqrt[u]");
    check_eq("Simplify[(x+Sqrt[a])/Sqrt[a]]", "1 + x/Sqrt[a]");
    check_eq("Simplify[(x+Sqrt[y]) y^(-1/2)]", "1 + x/Sqrt[y]");
    check_eq("Factor[1 + x/Sqrt[a]]", "(Sqrt[a] + x)/Sqrt[a]");
    /* The narrowed polish must still fire for a COMPOUND radicand. */
    check_eq("Simplify[(Sqrt[6] Sqrt[6 + x^2])/(6 x + x^3)]",
             "Sqrt[6]/(x Sqrt[6 + x^2])");
    /* Unrelated simplifications must be unchanged. */
    check_eq("Simplify[Sin[x]^2+Cos[x]^2]", "1");
    check_eq("Simplify[(x^2-1)/(x-1)]", "1 + x");

    /* ------------------------------------------------------------------ */
    /* TimeConstraint option.                                              */
    /* ------------------------------------------------------------------ */
    /* A generous / Infinite budget is transparent on easy inputs. */
    check_eq("Simplify[(x+Sqrt[u]) u^(-1/2), TimeConstraint -> 5]",
             "1 + x/Sqrt[u]");
    check_eq("Simplify[(x+Sqrt[u]) u^(-1/2), TimeConstraint -> Infinity]",
             "1 + x/Sqrt[u]");
    /* The option registers a default and shows up in Options[]. */
    check_eq("Options[Simplify]",
             "{Assumptions -> Automatic, ComplexityFunction -> Automatic, "
             "TransformationFunctions -> Automatic, TimeConstraint -> Infinity}");
    /* TimeConstraint must NOT be swallowed as a positional assumption: a real
     * positional assumption still applies alongside it. */
    check_eq("Simplify[Sqrt[x^2], x > 0]", "x");
    check_eq("Simplify[Sqrt[x^2], x > 0, TimeConstraint -> 5]", "x");
    check_eq("Simplify[Sqrt[x^2], Assumptions -> x > 0, TimeConstraint -> 2]", "x");
    /* A tiny budget bails gracefully: terminates fast, returns a real form
     * (never $Aborted, never a crash). The exact form is timing-dependent, so
     * we only assert termination + non-abort. */
    expect_terminates("Simplify[(x+Sqrt[u]) u^(-1/2), TimeConstraint -> 0.0001]", 5);
    check_not_aborted("Simplify[(x+Sqrt[u]) u^(-1/2), TimeConstraint -> 0.0001]");

    if (failures) {
        fprintf(stderr, "\n%d Simplify-hang test(s) FAILED\n", failures);
        return 1;
    }
    printf("All Simplify-hang tests passed!\n");
    return 0;
}
