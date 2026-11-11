/* test_dsolve_pde_corpus.c -- self-verifying PDE corpus runner for DSolve[].
 *
 * The PDE twin of ../test_dsolve_corpus.c.  Loads the verifier prelude
 * (dsolve_pde_corpus_prelude.m) once via Get[], reads one DE_PDE_*.m corpus file
 * (N. Abbasi's "PDE in CAS", imported by tools/pde_corpus/import_abbasi_pde.py)
 * as a List literal WITHOUT evaluating it, and forks per record (corpus_fork.h)
 * to evaluate
 *
 *     pdeCheckCode[label, section, eqns, fn, vars, opts]
 *
 * whose integer verdict maps to one TSV line on STDOUT,
 *
 *     <tag>\t<label>\t<section>\t<colA>\t<colB>
 *
 *     P  PASS            explicit answer, every residual proved zero
 *     I  PASS-IMPLICIT   implicit relation verified by the implicit-function rule
 *     S  PASS-SERIES     eigenfunction series verified mode by mode + projection
 *     K  UNK             neither proved nor refuted (counted non-PASS)
 *     F  FAIL            a residual proved NONZERO -- a wrong answer
 *     N  UNFIT           an IVP/BVP answer still carrying C[k] (unsolved)
 *     U  UNEVAL          declined / timed out / crashed
 *
 * for tools/dsolve_pde_corpus_report.py.  Two gates, both from argv (CMake):
 *
 *     argv[3]  nonpass baseline   K + F + N + U + timeouts + crashes may not exceed
 *                                 it (a high-water mark each method wave lowers)
 *     argv[4]  FAIL baseline      F may not exceed it, SEPARATELY -- unlike the
 *                                 ODE gate, a new wrong answer cannot hide behind
 *                                 a new solve.  Driven to 0, then asserted empty.
 *
 * Overrides: DSOLVE_CORPUS_BASELINE / DSOLVE_CORPUS_FAIL_BASELINE (measurement),
 * DSOLVE_CORPUS_LIMIT (run a prefix while iterating).
 */

#include "test_utils.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void symtab_init(void);
void core_init(void);

#include "expr.h"
#include "../corpus_fork.h"

extern char* expr_to_string(struct Expr*);

/* The prelude runs DSolve under an 8 s TimeConstrained, then a handful of 4 s
 * PossibleZeroQ probes; the alarm bounds a method that ignores the budget. */
#define PDE_PER_CASE_TIMEOUT_SEC 20

/* verdict k -> pipe byte, in pdeCheckCode's numbering */
static const char PDE_VERDICT_BYTES[] = {
    'P',   /* 0 PASS          */
    'F',   /* 1 FAIL          */
    'U',   /* 2 UNEVAL        */
    'X',   /* 3 SKIP (reserved; counted as malformed) */
    'I',   /* 4 PASS-IMPLICIT */
    'S',   /* 5 PASS-SERIES   */
    'K',   /* 6 UNK           */
    'N',   /* 7 UNFIT         */
};
#define PDE_BAD_BYTE '?'
#define PDE_NFIELDS 8   /* {label, section, eqns, fn, vars, opts, colA, colB} */

/* A string field, unquoted for the TSV (the record holds "pde-1" etc.). */
static char* field_plain(Expr* rec, size_t idx) {
    Expr* a = rec->data.function.args[idx];
    if (a->type == EXPR_STRING) {
        const char* s = a->data.string;
        char* out = (char*)malloc(strlen(s) + 1);
        if (out) strcpy(out, s);
        return out;
    }
    return expr_to_string(a);
}

static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, "List") == 0;
}

static int env_int(const char* name, int dflt) {
    const char* v = getenv(name);
    return v ? atoi(v) : dflt;
}

int main(int argc, char** argv) {
    alarm(0);                         /* each forked case carries its own alarm */
    setvbuf(stderr, NULL, _IONBF, 0);

    if (argc < 3) {
        fprintf(stderr, "usage: %s <DE_PDE_*.m> <prelude.m> [nonpass-baseline] [fail-baseline]\n",
                argv[0]);
        return 2;
    }
    const char* corpus_file  = argv[1];
    const char* prelude_file = argv[2];
    int baseline      = env_int("DSOLVE_CORPUS_BASELINE", argc > 3 ? atoi(argv[3]) : -1);
    int fail_baseline = env_int("DSOLVE_CORPUS_FAIL_BASELINE", argc > 4 ? atoi(argv[4]) : 0);
    long limit = env_int("DSOLVE_CORPUS_LIMIT", -1);

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
    if (!fbuf) { fprintf(stderr, "FAIL: cannot read %s\n", corpus_file); return 1; }
    Expr* cases = parse_expression(fbuf);
    free(fbuf);
    if (!is_list(cases)) {
        fprintf(stderr, "FAIL: %s did not parse as a List literal.\n", corpus_file);
        if (cases) expr_free(cases);
        return 1;
    }

    size_t n = cases->data.function.arg_count;
    if (limit >= 0 && (size_t)limit < n) n = (size_t)limit;
    fprintf(stderr, "==> DSolve PDE corpus: %zu cases\n", n);

    int pass = 0, impl = 0, series = 0, unk = 0, fail = 0, unfit = 0, uneval = 0;
    int timed_out = 0, crashed = 0, malformed = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* rec = cases->data.function.args[i];
        if (!is_list(rec) || rec->data.function.arg_count != PDE_NFIELDS) {
            malformed++;
            continue;
        }
        char* label   = field_plain(rec, 0);
        char* section = field_plain(rec, 1);
        char* colA    = field_plain(rec, 6);   /* Abbasi's two per-problem */
        char* colB    = field_plain(rec, 7);   /* result columns (Y/N/?)   */
        fprintf(stderr, "  [%4zu/%zu] %s\n", i + 1, n, label ? label : "?");

        Expr* call = expr_new_function(expr_new_symbol("pdeCheckCode"),
            (Expr*[]){
                expr_copy(rec->data.function.args[0]), expr_copy(rec->data.function.args[1]),
                expr_copy(rec->data.function.args[2]), expr_copy(rec->data.function.args[3]),
                expr_copy(rec->data.function.args[4]), expr_copy(rec->data.function.args[5])
            }, 6);
        char code = 0;
        CorpusOutcome oc = corpus_run_case(call, PDE_PER_CASE_TIMEOUT_SEC,
                                           PDE_VERDICT_BYTES,
                                           (int)sizeof PDE_VERDICT_BYTES, PDE_BAD_BYTE, &code);
        expr_free(call);

        char tag = 'U';
        if (oc == CORPUS_TIMEOUT) {
            timed_out++;
            fprintf(stderr, "      -> TIMEOUT (>%ds)\n", PDE_PER_CASE_TIMEOUT_SEC);
        } else if (oc == CORPUS_CRASH) {
            crashed++;
        } else {
            switch (code) {
                case 'P': pass++;   tag = 'P'; break;
                case 'I': impl++;   tag = 'I'; break;
                case 'S': series++; tag = 'S'; break;
                case 'K': unk++;    tag = 'K'; break;
                case 'F': fail++;   tag = 'F'; break;
                case 'N': unfit++;  tag = 'N'; break;
                case 'U': uneval++; tag = 'U'; break;
                default:
                    malformed++;
                    fprintf(stderr, "      -> bad verdict 0x%02x\n", (unsigned char)code);
                    break;
            }
        }
        printf("%c\t%s\t%s\t%s\t%s\n", tag, label ? label : "?", section ? section : "?",
               colA ? colA : "?", colB ? colB : "?");
        fflush(stdout);
        free(label); free(section); free(colA); free(colB);
    }
    expr_free(cases);

    int passed  = pass + impl + series;
    int nonpass = unk + fail + unfit + uneval + timed_out + crashed;

    fprintf(stderr, "\n=== DSolve PDE corpus result ===\n");
    fprintf(stderr, "  Cases:          %zu\n", n);
    fprintf(stderr, "  PASS:           %d\n", pass);
    fprintf(stderr, "  PASS-IMPLICIT:  %d\n", impl);
    fprintf(stderr, "  PASS-SERIES:    %d\n", series);
    fprintf(stderr, "  UNK:            %d\n", unk);
    fprintf(stderr, "  FAIL:           %d\n", fail);
    fprintf(stderr, "  UNFIT:          %d\n", unfit);
    fprintf(stderr, "  UNEVAL:         %d\n", uneval);
    fprintf(stderr, "  Timed out:      %d\n", timed_out);
    fprintf(stderr, "  Crashed:        %d\n", crashed);
    fprintf(stderr, "  Malformed:      %d\n", malformed);
    fprintf(stderr, "  Passed:         %d / %zu   Non-PASS: %d\n", passed, n, nonpass);
    fprintf(stderr, "================================\n");

    int rc = 0;
    if (malformed > 0) {
        fprintf(stderr, "FAIL: %d malformed record(s) / verdict(s).\n", malformed);
        rc = 1;
    }
    if (fail > fail_baseline) {
        fprintf(stderr, "FAIL: %d wrong answer(s) exceed the FAIL baseline %d.\n",
                fail, fail_baseline);
        rc = 1;
    }
    if (baseline >= 0 && nonpass > baseline) {
        fprintf(stderr, "FAIL: %d non-PASS case(s) exceed baseline %d "
                        "(new DSolve regression).\n", nonpass, baseline);
        rc = 1;
    }
    if (rc == 0)
        fprintf(stderr, "DSolve PDE corpus within baselines (non-PASS %d <= %d, FAIL %d <= %d).\n",
                nonpass, baseline, fail, fail_baseline);
    return rc;
}
