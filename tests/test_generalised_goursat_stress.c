/* test_generalised_goursat_stress.c -- nightly: every GeneralisedGoursat corpus.
 *
 * Runs all of tests/data/generalised_goursat/*.m -- 1,183 integrals exported
 * from the research directory by tools/goursat/export_corpora.py: the
 * test_goursat.py checks and Chebyshev sweep, run_examples.py, the literature
 * list, the 30 welz_cases.py benchmarks, stress_goursat.py's nine constructed
 * families at three seeds, and the two pseudo-elliptic example lists -- through
 * tests/data/generalised_goursat/check_corpus.m, with a 120 s per-case budget
 * (as stress_goursat.py).  Every answer is checked by the package's exact
 * verifier; a verdict is judged against the constructed ground truth.
 *
 * Fails on any WRONG, DECLINE, UNEVAL or TIMEOUT.  MISSED ("undecided" where
 * the truth is known) is incompleteness, not a bug; the cases allowed to be
 * MISSED are listed in ALLOWED_MISSED with the reason.
 */

#define _POSIX_C_SOURCE 200809L   /* access() */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

static const char* CORPORA[] = {
    "unit", "examples", "literature", "welz",
    "stress_2026", "stress_7", "stress_11",
    "pseudoelliptic_sqrt", "pseudoelliptic_cbrt",
};

/* "undecided" with a known truth, accepted: the residual is of the third kind
 * and its elementarity is a torsion question the symmetry reductions do not
 * decide (paper, "What remains undecided"). */
static const char* ALLOWED_MISSED[] = {
    "order 4: 1/((t^4+3) sqrt(t^4-1))",          /* truth "undecided" in test_goursat.py */
    NULL,
};

static bool allowed_missed(const char* label) {
    for (int i = 0; ALLOWED_MISSED[i]; i++)
        if (strstr(label, ALLOWED_MISSED[i])) return true;
    return false;
}

static Expr* eval_str(const char* s) {
    Expr* p = parse_expression(s);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    expr_free(p);
    return r;
}

/* Locate tests/data/generalised_goursat from the build directory. */
static const char* data_dir(void) {
    static const char* cands[] = { "../data/generalised_goursat", "../../tests/data/generalised_goursat",
                                   "tests/data/generalised_goursat", "data/generalised_goursat" };
    char probe[512];
    for (size_t i = 0; i < sizeof(cands) / sizeof(cands[0]); i++) {
        snprintf(probe, sizeof(probe), "%s/check_corpus.m", cands[i]);
        if (access(probe, R_OK) == 0) return cands[i];
    }
    return NULL;
}

int main(void) {
    alarm(0);                              /* nightly: no 120 s unit watchdog */
    symtab_init();
    core_init();
    expr_free(eval_str("Integrate`GeneralisedGoursat[1/Sqrt[t^3 - 1], t]"));  /* lazy load */

    const char* dir = data_dir();
    ASSERT_MSG(dir != NULL, "cannot find tests/data/generalised_goursat");
    char buf[1024];
    snprintf(buf, sizeof(buf), "Get[\"%s/check_corpus.m\"]", dir);
    expr_free(eval_str(buf));

    int bad = 0, total_ok = 0;
    for (size_t c = 0; c < sizeof(CORPORA) / sizeof(CORPORA[0]); c++) {
        snprintf(buf, sizeof(buf), "ggCheckCorpus[\"%s/%s.m\", 120]", dir, CORPORA[c]);
        Expr* r = eval_str(buf);
        ASSERT(r && r->type == EXPR_FUNCTION && r->data.function.arg_count == 2);
        Expr* okc = r->data.function.args[0];
        Expr* rows = r->data.function.args[1];
        int nok = okc->type == EXPR_INTEGER ? (int)okc->data.integer : -1;
        total_ok += nok;
        printf("%-22s OK=%d  not OK=%zu\n", CORPORA[c], nok, rows->data.function.arg_count);
        for (size_t i = 0; i < rows->data.function.arg_count; i++) {
            Expr* row = rows->data.function.args[i];
            const char* label = row->data.function.args[0]->data.string;
            const char* verdict = row->data.function.args[1]->data.string;
            bool fine = strcmp(verdict, "MISSED") == 0 && allowed_missed(label);
            char* st = expr_to_string(row->data.function.args[2]);
            char* sec = expr_to_string(row->data.function.args[3]);
            printf("  %-8s %s  [%s, %s s]%s\n", verdict, label, st, sec, fine ? "  (allowed)" : "");
            free(st); free(sec);
            if (!fine) bad++;
        }
        expr_free(r);
    }
    printf("total OK=%d, failures=%d\n", total_ok, bad);
    ASSERT_MSG(bad == 0, "%d GeneralisedGoursat corpus cases are not OK", bad);
    printf("All GeneralisedGoursat stress corpora passed!\n");
    return 0;
}
