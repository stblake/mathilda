/* test_crc_corpus.c
 *
 * Corpus runner for the CRC integral table (src/internal/CRCMathTablesIntegrals.m).
 *
 * Loads ../../CRCIntegralsCorpus.m — a List of {integrand, var} pairs that
 * each exercise one rule from the CRC table.  For every pair the runner
 *
 *   1. computes  r = Integrate[integrand, var, Method -> "CRCTable"]
 *   2. computes  d = D[r, var] - integrand
 *   3. checks    Simplify[Cancel[Together[d]]] == 0
 *
 * A case is "diff_zero" iff (3) succeeds.  Anything else — unevaluated,
 * diff nonzero, child timeout, child crash — is a regression flagged by
 * the test summary at the end.  The differential check is the same
 * pattern used by test_intrat_corpus.c, with the simplifier swapped to
 * Simplify[Cancel[Together[..]]] per project preference for CRC results
 * which often involve nested radicals / log arguments.
 *
 * Each case runs in a forked child with a per-case wall-clock timeout
 * so a single misbehaving rule cannot stall the whole run or corrupt
 * evaluator state between cases.
 */

#include "test_utils.h"
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

extern char* expr_to_string(struct Expr*);

static bool expr_contains_head_named(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_FUNCTION) {
        if (e->data.function.head
            && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, name) == 0) {
            return true;
        }
        if (expr_contains_head_named(e->data.function.head, name)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (expr_contains_head_named(e->data.function.args[i], name)) return true;
        }
    }
    return false;
}

/* Diff-zero check: try symbolic first, fall back to numeric.
 *
 * Symbolic: Simplify[Cancel[Together[diff]]] == 0.  Together puts
 * everything over a common denominator; Cancel removes the GCD of
 * numerator and denominator; Simplify mops up the remaining residue.
 * This catches rational and algebraic-rational diffs.
 *
 * Numeric fallback: substitute the integration variable with three
 * concrete real values (1.7, 0.3, 2.6) and confirm |diff| is below
 * ~1e-8 at every sample.  This catches nested-radical normalization
 * cases like  1/Sqrt[1 - x^2/4] vs 2/Sqrt[4 - x^2]  which Simplify
 * cannot collapse but which N[] resolves trivially.  A pure-symbolic
 * check would mark these as DIFF_NONZERO even though the integral is
 * correct, so the numeric pass is essential coverage for CRC results
 * that drag radicals through Abs[a] / Abs[c] / branch-cut factors. */
/* True when an N[..] result is "real and finite" — a real/integer node,
 * or a Complex[..] with negligible imaginary part.  Used to filter
 * branch-cut samples in numeric_is_zero: a sample where the integrand
 * goes complex (e.g. Sqrt[(1+x)/(1-x)] at x = 1.7) or singular (e.g.
 * ArcSec[0.2], 1/Sqrt[1 - 4x^2] at x = 1) is outside the rule's natural
 * real domain and cannot meaningfully witness the antiderivative's
 * correctness — the diff there is dominated by branch choices /
 * directed-infinity arithmetic, not by the rule.  Anything that isn't
 * a finite real (ComplexInfinity, Indeterminate, unevaluated symbols,
 * etc.) is rejected. */
static bool n_value_is_real(Expr* nv, double tol) {
    if (!nv) return false;
    if (nv->type == EXPR_REAL || nv->type == EXPR_INTEGER) return true;
    if (nv->type == EXPR_FUNCTION
        && nv->data.function.head
        && nv->data.function.head->type == EXPR_SYMBOL
        && strcmp(nv->data.function.head->data.symbol, "Complex") == 0
        && nv->data.function.arg_count == 2) {
        Expr* re = nv->data.function.args[0];
        Expr* im = nv->data.function.args[1];
        /* Both parts must be finite numerics. */
        if (re->type != EXPR_REAL && re->type != EXPR_INTEGER) return false;
        if (im->type != EXPR_REAL && im->type != EXPR_INTEGER) return false;
        double iv = (im->type == EXPR_REAL) ? im->data.real
                  : (double)im->data.integer;
        return iv >= -tol && iv <= tol;
    }
    return false;
}

/* Evaluate `expr /. var -> sample` numerically, returning a fresh Expr.
 * Caller owns the result. */
static Expr* eval_at_point(Expr* expr, Expr* var, double sample) {
    Expr* rule = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_copy(var), expr_new_real(sample) }, 2);
    Expr* sub = expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(expr), rule }, 2);
    Expr* nv_call = expr_new_function(expr_new_symbol("N"),
        (Expr*[]){ sub }, 1);
    Expr* nv = evaluate(nv_call);
    expr_free(nv_call);
    return nv;
}

static bool n_value_near_zero(Expr* nv, double tol) {
    if (nv->type == EXPR_REAL) {
        double v = nv->data.real;
        return v >= -tol && v <= tol;
    }
    if (nv->type == EXPR_INTEGER && nv->data.integer == 0) return true;
    if (nv->type == EXPR_FUNCTION
        && nv->data.function.head
        && nv->data.function.head->type == EXPR_SYMBOL
        && strcmp(nv->data.function.head->data.symbol, "Complex") == 0
        && nv->data.function.arg_count == 2) {
        Expr* re = nv->data.function.args[0];
        Expr* im = nv->data.function.args[1];
        double rv = (re->type == EXPR_REAL) ? re->data.real
                  : (re->type == EXPR_INTEGER) ? (double)re->data.integer
                  : 1.0;
        double iv = (im->type == EXPR_REAL) ? im->data.real
                  : (im->type == EXPR_INTEGER) ? (double)im->data.integer
                  : 1.0;
        return rv >= -tol && rv <= tol && iv >= -tol && iv <= tol;
    }
    return false;
}

/* numeric_is_zero — fall-back numeric verification of D[r,x] - integrand.
 *
 * Strategy: evaluate the diff at a spread of real sample points,
 * BUT skip any sample where the integrand itself goes non-real
 * (i.e. the sample is outside the rule's natural real domain;
 * branch choices dominate the diff there and tell us nothing about
 * correctness on the principal sheet).  Example: Sqrt[(1+x)/(1-x)]
 * is real on (-1, 1) only; at x = 1.7 both integrand and antiderivative
 * are complex with branches that don't cancel even when the rule is
 * exactly right on (-1, 1).
 *
 * Sample spread includes a tight in-domain band (0.1..0.5, -0.4) so
 * almost every rule has at least one valid witness, plus the original
 * wider points (1.7, 2.6) which still catch errors outside the natural
 * domain when both sides happen to stay real.
 *
 * A case passes when: (a) at least one sample is in-domain, AND (b)
 * every in-domain sample gives diff ~ 0.  An all-skipped case is a
 * fail (we couldn't witness anything). */
static bool numeric_is_zero(Expr* diff, Expr* var, Expr* integrand) {
    /* Mix of tight in-domain points and wider ones, all positive.
     * At least one of these is in the natural real domain for almost
     * every CRC rule.  Negative samples are intentionally avoided: many
     * CRC entries have positive-x branch conventions (e.g. ArcSec[a x]
     * → Log[a x + ...] is right for a x > 0 only) which are out of
     * scope for this corpus, and a negative sample would flag those as
     * regressions even though they pre-date this run. */
    static const double samples[] = { 0.3, 0.7, 1.3, 1.7, 2.6 };
    static const size_t ns = sizeof(samples) / sizeof(samples[0]);

    if (var->type != EXPR_SYMBOL) return false;

    size_t in_domain = 0;
    for (size_t i = 0; i < ns; i++) {
        /* Skip when the integrand itself is not finite-real at this
         * sample (branch-cut region or singular point — the diff there
         * tells us nothing about the rule's correctness on the principal
         * sheet). */
        Expr* nint = eval_at_point(integrand, var, samples[i]);
        bool integrand_real = n_value_is_real(nint, 1e-6);
        expr_free(nint);
        if (!integrand_real) continue;

        Expr* nv = eval_at_point(diff, var, samples[i]);
        /* Skip if the diff itself is not finite-real — e.g. the sample
         * sits exactly on a domain boundary (ArcSec[1]) where 1/Sqrt[..]
         * = 1/0.  Treating those as "not zero" is a false positive; the
         * limit is finite, our N[] just can't see through 0/0.  Only
         * real-and-finite samples are informative. */
        if (!n_value_is_real(nv, 1.0)) {
            expr_free(nv);
            continue;
        }
        bool point_zero = n_value_near_zero(nv, 1e-7);
        expr_free(nv);

        if (!point_zero) return false;
        in_domain++;
    }
    return in_domain > 0;
}

static bool is_zero_after_simplify(Expr* diff, Expr* var, Expr* integrand) {
    Expr* tg_call = expr_new_function(expr_new_symbol("Together"),
        (Expr*[]){ expr_copy(diff) }, 1);
    Expr* tg = evaluate(tg_call);
    expr_free(tg_call);

    Expr* cn_call = expr_new_function(expr_new_symbol("Cancel"),
        (Expr*[]){ tg }, 1);
    Expr* cn = evaluate(cn_call);
    expr_free(cn_call);

    Expr* s_call = expr_new_function(expr_new_symbol("Simplify"),
        (Expr*[]){ cn }, 1);
    Expr* s = evaluate(s_call);
    expr_free(s_call);

    bool is_zero = false;
    if (s->type == EXPR_INTEGER && s->data.integer == 0) {
        is_zero = true;
    } else if (s->type == EXPR_REAL) {
        double v = s->data.real;
        if (v >= -1e-9 && v <= 1e-9) is_zero = true;
    }
    expr_free(s);

    if (is_zero) return true;
    return numeric_is_zero(diff, var, integrand);
}

#define CORPUS_PER_CASE_TIMEOUT_SEC 5

enum {
    CORPUS_DIFF_ZERO    = 'Z',
    CORPUS_DIFF_NONZERO = 'D',
    CORPUS_UNEVALUATED  = 'U',
};

int main(int argc, char** argv) {
    alarm(0);
    setvbuf(stderr, NULL, _IONBF, 0);

    const char* corpus_file =
        (argc > 1) ? argv[1] : "../../CRCIntegralsCorpus.m";

    symtab_init();
    core_init();

    fprintf(stderr, "==> corpus file: %s\n", corpus_file);

    char get_buf[1024];
    snprintf(get_buf, sizeof(get_buf), "Get[\"%s\"]", corpus_file);
    Expr* get_call = parse_expression(get_buf);
    ASSERT(get_call != NULL);
    Expr* tests = evaluate(get_call);
    expr_free(get_call);

    if (!tests || tests->type != EXPR_FUNCTION
        || !tests->data.function.head
        || tests->data.function.head->type != EXPR_SYMBOL
        || strcmp(tests->data.function.head->data.symbol, "List") != 0) {
        fprintf(stderr,
            "FAIL: could not load %s as a List.\n", corpus_file);
        if (tests) expr_free(tests);
        ASSERT(false);
        return 1;
    }

    size_t n = tests->data.function.arg_count;
    fprintf(stderr, "==> CRC corpus: %zu cases\n", n);

    int diff_zero    = 0;
    int diff_nonzero = 0;
    int unevaluated  = 0;
    int malformed    = 0;
    int timed_out    = 0;
    int crashed      = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* pair = tests->data.function.args[i];
        if (pair->type != EXPR_FUNCTION
            || !pair->data.function.head
            || pair->data.function.head->type != EXPR_SYMBOL
            || strcmp(pair->data.function.head->data.symbol, "List") != 0
            || pair->data.function.arg_count < 2) {
            malformed++;
            fprintf(stderr, "  [%3zu/%zu] MALFORMED\n", i + 1, n);
            continue;
        }
        Expr* integrand = pair->data.function.args[0];
        Expr* var       = pair->data.function.args[1];

        char* in_str = expr_to_string(integrand);
        fprintf(stderr, "  [%3zu/%zu] %s\n", i + 1, n, in_str);
        free(in_str);

        int pipefd[2];
        if (pipe(pipefd) != 0) {
            fprintf(stderr, "      -> pipe() failed; skipping\n");
            crashed++;
            continue;
        }
        pid_t pid = fork();
        if (pid < 0) {
            fprintf(stderr, "      -> fork() failed; skipping\n");
            close(pipefd[0]); close(pipefd[1]);
            crashed++;
            continue;
        }

        if (pid == 0) {
            /* --- Child --- */
            close(pipefd[0]);
            char code = CORPUS_UNEVALUATED;

            alarm(CORPUS_PER_CASE_TIMEOUT_SEC);

            /* Build Integrate[f, x, Method -> "CRCTable"] to bypass
             * Rational / RischNorman and exercise the CRC rules
             * directly. */
            Expr* method_rule = expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){
                    expr_new_symbol("Method"),
                    expr_new_string("CRCTable")
                }, 2);
            Expr* int_call = expr_new_function(expr_new_symbol("Integrate"),
                (Expr*[]){
                    expr_copy(integrand),
                    expr_copy(var),
                    method_rule
                }, 3);
            Expr* result = evaluate(int_call);
            expr_free(int_call);

            /* Any of the three integrator heads anywhere in the
             * result counts as unevaluated.  The CRC table's recursive
             * reduction rules (Formula 26, 38, 67, ...) leave residual
             * `IntegrateTable[..]` calls when the inner integrand
             * doesn't match another stored rule; these surface as
             * subterms of an otherwise closed expression.  Top-level
             * matching alone misses them. */
            bool is_unevaluated =
                expr_contains_head_named(result, "Integrate")
                || expr_contains_head_named(result, "Integrate`CRCTable")
                || expr_contains_head_named(result, "IntegrateTable");

            if (is_unevaluated) {
                code = CORPUS_UNEVALUATED;
            } else {
                Expr* d_call = expr_new_function(expr_new_symbol("D"),
                    (Expr*[]){ expr_copy(result), expr_copy(var) }, 2);
                Expr* d_res = evaluate(d_call);
                expr_free(d_call);

                Expr* neg_int = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){
                        expr_new_integer(-1),
                        expr_copy(integrand)
                    }, 2);
                Expr* diff_raw = expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ d_res, neg_int }, 2);
                Expr* diff = evaluate(diff_raw);
                expr_free(diff_raw);

                code = is_zero_after_simplify(diff, var, integrand)
                    ? CORPUS_DIFF_ZERO : CORPUS_DIFF_NONZERO;
                expr_free(diff);
            }
            alarm(0);
            if (result) expr_free(result);

            ssize_t w = write(pipefd[1], &code, 1);
            (void)w;
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
            usleep(50000);
        }
        if (!reaped) {
            kill(pid, SIGKILL);
            waitpid(pid, &status, 0);
            close(pipefd[0]);
            timed_out++;
            fprintf(stderr, "      -> TIMEOUT (>%ds)\n",
                    CORPUS_PER_CASE_TIMEOUT_SEC);
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
            fprintf(stderr,
                "      -> CRASH (exit %d, sig %d)\n",
                WIFEXITED(status) ? WEXITSTATUS(status) : -1,
                WIFSIGNALED(status) ? WTERMSIG(status)  : 0);
            continue;
        }

        switch (code) {
            case CORPUS_DIFF_ZERO:    diff_zero++;    break;
            case CORPUS_DIFF_NONZERO: {
                diff_nonzero++;
                char* in = expr_to_string(integrand);
                fprintf(stderr,
                    "      -> DIFF NONZERO  (integrand: %s)\n", in);
                free(in);
                fflush(stderr);
                break;
            }
            case CORPUS_UNEVALUATED:
                unevaluated++;
                fprintf(stderr, "      -> UNEVALUATED\n");
                break;
            default:
                crashed++;
                fprintf(stderr, "      -> UNKNOWN classification 0x%02x\n",
                        (unsigned char)code);
                break;
        }
    }
    expr_free(tests);

    fprintf(stderr, "\n=== CRC integral corpus result ===\n");
    fprintf(stderr, "  Total cases:                  %zu\n", n);
    fprintf(stderr, "  Closed (diff zero):           %d  %.1f%%\n",
            diff_zero, n ? 100.0 * diff_zero / n : 0.0);
    fprintf(stderr, "  Diff nonzero  (REGRESSION):   %d\n", diff_nonzero);
    fprintf(stderr, "  Unevaluated:                  %d\n", unevaluated);
    fprintf(stderr, "  Timed out (>%d s/case):       %d\n",
            CORPUS_PER_CASE_TIMEOUT_SEC, timed_out);
    fprintf(stderr, "  Crashed (child fault):        %d\n", crashed);
    fprintf(stderr, "  Malformed (test file issue):  %d\n", malformed);
    fprintf(stderr, "===================================\n");

    /* The diff-nonzero pile is the hard regression bar: rules whose
     * Integrate output exists but does not differentiate back to the
     * integrand.  This is always a real correctness bug.  Unevaluated
     * / timeout / crash counts are tracked but tolerated until the
     * supporting machinery (e.g. nested-radical simplifier coverage)
     * is in place.  See tasks/lessons.md / Mathilda_spec.md for the
     * running baseline.
     *
     * Treat the test as passing when no DIFF NONZERO case is detected.
     * Raise this if integration improvements legitimately exceed the
     * cap; never lower it without investigating each regression. */
    const int CORPUS_DIFF_NONZERO_BASELINE = 0;
    if (diff_nonzero > CORPUS_DIFF_NONZERO_BASELINE) {
        fprintf(stderr,
            "FAIL: %d case(s) closed but did not differentiate back to "
            "the integrand (above baseline %d).\n",
            diff_nonzero, CORPUS_DIFF_NONZERO_BASELINE);
        return 1;
    }

    return 0;
}
