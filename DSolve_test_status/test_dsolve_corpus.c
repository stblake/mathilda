/* test_dsolve_corpus.c
 *
 * Self-verifying stress-corpus runner for DSolve[], modelled on
 * test_solve_corpus.c.  Loads the verifier prelude
 * (tests/dsolve_corpus_prelude.m) once in the parent via Get[], reads the
 * corpus (DE_examples_2.m -- the 12000.org section 2.1.2 ODEs) as a List
 * literal WITHOUT evaluating it, then forks per case with a wall-clock
 * alarm.  Each child evaluates
 *
 *     dsolveCheckCode[label, eqn, fn, iv, classif]
 *
 * which runs DSolve under TimeConstrained and numerically back-substitutes
 * every explicit branch, returning an integer verdict:
 *   0 = PASS, 1 = FAIL (wrong closed form), 2 = UNEVAL (declined/timeout/{}),
 *   3 = SKIP (system_of_ODEs -- this is the scalar-first campaign).
 *
 * The child maps the verdict to a single pipe byte; the parent tallies and
 * prints one TSV line per case to STDOUT ("<P|F|U|S>\t<label>\t<classif>")
 * so tools/dsolve_corpus_report.py can bucket the gaps by Maple class.
 *
 * The corpus is a PROGRESS DASHBOARD, not an all-green assertion: the runner
 * fails only when the number of non-PASS scalar cases exceeds the checked-in
 * high-water mark DSOLVE_CORPUS_FAIL_BASELINE.  Each landed method-wave must
 * LOWER the baseline; a case that newly fails/crashes/times out trips the
 * test.  Overridable at run time via env DSOLVE_CORPUS_BASELINE (for the
 * initial measurement) and DSOLVE_CORPUS_LIMIT (to run a prefix of the
 * corpus during iteration).
 */

#include "test_utils.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void symtab_init(void);
void core_init(void);

#include "expr.h"
#include "corpus_fork.h"   /* fork-per-case core shared with the PDE runner */

extern char* expr_to_string(struct Expr*);

/* Per-case wall-clock cap.  The prelude runs DSolve under an 8 s
 * TimeConstrained plus a numeric back-substitution sweep; allow headroom
 * for a method that ignores TimeConstrained while still bounding a hang. */
#define DSOLVE_PER_CASE_TIMEOUT_SEC 20

/* Single-byte verdict codes returned from each child over the pipe. */
enum {
    DS_CODE_PASS   = 'P',  /* verdict 0: verified closed form              */
    DS_CODE_FAIL   = 'F',  /* verdict 1: a branch is demonstrably nonzero  */
    DS_CODE_UNEVAL = 'U',  /* verdict 2: declined / timed out / {}         */
    DS_CODE_SKIP   = 'S',  /* verdict 3: system_of_ODEs (not counted)      */
    DS_CODE_BADVAL = 'X',  /* dsolveCheckCode returned a non-integer       */
};

/* extract a plain string field (label / classif) for the TSV, best-effort. */
static char* field_str(Expr* rec, size_t idx) {
    if (rec->type == EXPR_FUNCTION && rec->data.function.arg_count > idx) {
        Expr* a = rec->data.function.args[idx];
        char* s = expr_to_string(a);
        return s ? s : NULL;
    }
    return NULL;
}

int main(int argc, char** argv) {
    /* The default test_utils alarm is per-binary; disable it -- each forked
     * case carries its own alarm instead. */
    alarm(0);
    setvbuf(stderr, NULL, _IONBF, 0);

    const char* corpus_file  = (argc > 1) ? argv[1] : "../../DSolve_test_status/DE_examples_2.m";
    const char* prelude_file = (argc > 2) ? argv[2] : "../../DSolve_test_status/dsolve_corpus_prelude.m";

    int baseline = -1;
    if (argc > 3) baseline = atoi(argv[3]);            /* per-section baseline (CMake) */
    const char* be = getenv("DSOLVE_CORPUS_BASELINE"); /* manual measurement override */
    if (be) baseline = atoi(be);
    long limit = -1;
    const char* le = getenv("DSOLVE_CORPUS_LIMIT");
    if (le) limit = atol(le);

    symtab_init();
    core_init();

    fprintf(stderr, "==> prelude: %s\n", prelude_file);
    Expr* get_call = expr_new_function(expr_new_symbol("Get"),
        (Expr*[]){ expr_new_string(prelude_file) }, 1);
    Expr* get_res = evaluate(get_call);
    expr_free(get_call);
    if (get_res) expr_free(get_res);

    fprintf(stderr, "==> corpus:  %s\n", corpus_file);
    char* fbuf = corpus_slurp(corpus_file);
    if (!fbuf) {
        fprintf(stderr, "FAIL: cannot read %s\n", corpus_file);
        ASSERT(false);
        return 1;
    }
    Expr* cases = parse_expression(fbuf);
    free(fbuf);

    if (!cases || cases->type != EXPR_FUNCTION || !cases->data.function.head
        || cases->data.function.head->type != EXPR_SYMBOL
        || strcmp(cases->data.function.head->data.symbol.name, "List") != 0) {
        fprintf(stderr, "FAIL: %s did not parse as a List literal.\n", corpus_file);
        if (cases) expr_free(cases);
        ASSERT(false);
        return 1;
    }

    size_t n = cases->data.function.arg_count;
    if (limit >= 0 && (size_t)limit < n) n = (size_t)limit;
    fprintf(stderr, "==> DSolve corpus: %zu cases\n", n);

    int passed = 0, failed = 0, unevaluated = 0, skipped = 0;
    int malformed = 0, timed_out = 0, crashed = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* rec = cases->data.function.args[i];
        if (rec->type != EXPR_FUNCTION || !rec->data.function.head
            || rec->data.function.head->type != EXPR_SYMBOL
            || strcmp(rec->data.function.head->data.symbol.name, "List") != 0
            || rec->data.function.arg_count < 5) {
            malformed++;
            continue;
        }
        char* label_str = field_str(rec, 0);
        char* classif_str = field_str(rec, 4);
        fprintf(stderr, "  [%4zu/%zu] %s\n", i + 1, n, label_str ? label_str : "?");

        static const char verdict_bytes[] = {
            DS_CODE_PASS, DS_CODE_FAIL, DS_CODE_UNEVAL, DS_CODE_SKIP };
        Expr* call = expr_new_function(expr_new_symbol("dsolveCheckCode"),
            (Expr*[]){
                expr_copy(rec->data.function.args[0]),
                expr_copy(rec->data.function.args[1]),
                expr_copy(rec->data.function.args[2]),
                expr_copy(rec->data.function.args[3]),
                expr_copy(rec->data.function.args[4])
            }, 5);
        char code = 0;
        CorpusOutcome oc = corpus_run_case(call, DSOLVE_PER_CASE_TIMEOUT_SEC,
                                           verdict_bytes, 4, DS_CODE_BADVAL, &code);
        expr_free(call);
        if (oc != CORPUS_RAN) {
            if (oc == CORPUS_TIMEOUT) {
                timed_out++;
                fprintf(stderr, "      -> TIMEOUT (>%ds)\n", DSOLVE_PER_CASE_TIMEOUT_SEC);
            } else {
                crashed++;
            }
            printf("U\t%s\t%s\n", label_str ? label_str : "?", classif_str ? classif_str : "");
            fflush(stdout);
            free(label_str); free(classif_str);
            continue;
        }

        char tag = 'U';
        switch (code) {
            case DS_CODE_PASS:   passed++;      tag = 'P'; break;
            case DS_CODE_FAIL:   failed++;      tag = 'F'; break;
            case DS_CODE_UNEVAL: unevaluated++; tag = 'U'; break;
            case DS_CODE_SKIP:   skipped++;     tag = 'S'; break;
            default: crashed++;  tag = 'U';
                fprintf(stderr, "      -> bad verdict 0x%02x\n", (unsigned char)code);
                break;
        }
        printf("%c\t%s\t%s\n", tag, label_str ? label_str : "?", classif_str ? classif_str : "");
        fflush(stdout);
        free(label_str); free(classif_str);
    }
    expr_free(cases);

    int scalar_total = passed + failed + unevaluated + timed_out + crashed;
    int nonpass = failed + unevaluated + timed_out + crashed;

    fprintf(stderr, "\n=== DSolve corpus result ===\n");
    fprintf(stderr, "  Total cases:   %zu\n", n);
    fprintf(stderr, "  Skipped (sys): %d\n", skipped);
    fprintf(stderr, "  Scalar cases:  %d\n", scalar_total);
    fprintf(stderr, "  Passed:        %d\n", passed);
    fprintf(stderr, "  Failed:        %d\n", failed);
    fprintf(stderr, "  Unevaluated:   %d\n", unevaluated);
    fprintf(stderr, "  Timed out:     %d\n", timed_out);
    fprintf(stderr, "  Crashed:       %d\n", crashed);
    fprintf(stderr, "  Malformed:     %d\n", malformed);
    fprintf(stderr, "  Non-PASS:      %d / %d scalar\n", nonpass, scalar_total);
    fprintf(stderr, "============================\n");

    /* Regression baseline: checked-in high-water mark of non-PASS scalar
     * cases.  Set from the first full measurement; each method-wave LOWERS
     * it.  Override at run time with DSOLVE_CORPUS_BASELINE. */
    const int DSOLVE_CORPUS_FAIL_BASELINE = 605;   /* 2.1.2 fallback; per-section via argv[3] */
    int eff_baseline = (baseline >= 0) ? baseline : DSOLVE_CORPUS_FAIL_BASELINE;

    if (malformed > 0) {
        fprintf(stderr, "FAIL: %d malformed corpus record(s).\n", malformed);
        return 1;
    }
    if (nonpass > eff_baseline) {
        fprintf(stderr, "FAIL: %d non-PASS scalar case(s) exceed baseline %d "
                        "(new DSolve regression).\n", nonpass, eff_baseline);
        return 1;
    }
    fprintf(stderr, "DSolve corpus within baseline (%d non-PASS <= %d).\n",
            nonpass, eff_baseline);
    return 0;
}
