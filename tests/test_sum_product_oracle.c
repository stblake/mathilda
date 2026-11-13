/* test_sum_product_oracle.c
 *
 * Ratchet gate for the Sum / Product oracle corpus
 * (tests/sum_product_oracle/corpus.m, verifier prelude.m).
 *
 *   sum_product_oracle_tests <prelude.m> <corpus.m> <baseline.txt>
 *
 * Every corpus record is {label, Hold[input], kind, data, subs}; the child
 * process evaluates spoCheck @@ record, which returns "PASS" / "WRONG" /
 * "UNEVAL" / "NOCHECK".  Each record runs in its own forked child under a
 * wall-clock alarm, so a hang or a crash costs one record, not the run.
 *
 * baseline.txt lists, one per line, the labels known NOT to pass (with the
 * verdict and triage ID after a '#').  The gate:
 *   - FAILS when a label outside the baseline does not PASS (a regression,
 *     or a new record that was never baselined);
 *   - FAILS when a baseline line names a label that is not in the corpus
 *     (a stale line hides nothing and must be deleted);
 *   - REPORTS a baseline label that now passes, so its line can be deleted.
 * Each fix campaign therefore shrinks the baseline; it never grows except
 * when a newly added record is deliberately baselined.  SPO_VERBOSE=1 also
 * prints the current verdict of every baselined record.
 */
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "expr.h"

void symtab_init(void);
void core_init(void);

#define SPO_PER_CASE_TIMEOUT_SEC 20
#define SPO_MAX_BASELINE 1024
#define SPO_LABEL_MAX 128

/* Evaluate Get["path"] and discard the result; returns the result Expr
 * (caller frees) when keep is true. */
static Expr* spo_get(const char* path, int keep) {
    Expr* call = expr_new_function(expr_new_symbol("Get"),
        (Expr*[]){ expr_new_string(path) }, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    if (!keep && r) { expr_free(r); r = NULL; }
    return r;
}

static int spo_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, "List") == 0;
}

/* Baseline: labels known not to pass.  "label   # note" per line. */
typedef struct {
    char label[SPO_LABEL_MAX];
    int seen;
} SpoBase;

static int spo_read_baseline(const char* path, SpoBase* out, int max) {
    FILE* f = fopen(path, "r");
    if (!f) return -1;
    char line[512];
    int n = 0;
    while (fgets(line, sizeof line, f) && n < max) {
        char* p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\0') continue;
        size_t len = strcspn(p, " \t#\r\n");
        if (len == 0 || len >= SPO_LABEL_MAX) continue;
        memcpy(out[n].label, p, len);
        out[n].label[len] = '\0';
        out[n].seen = 0;
        n++;
    }
    fclose(f);
    return n;
}

static SpoBase* spo_find(SpoBase* b, int nb, const char* label) {
    for (int i = 0; i < nb; i++)
        if (strcmp(b[i].label, label) == 0) return &b[i];
    return NULL;
}

/* Run one record in a forked child; returns the verdict as a short static
 * string ("PASS", "WRONG", "UNEVAL", "NOCHECK", "TIMEOUT", "CRASH"). */
static const char* spo_run_record(Expr* rec) {
    int pipefd[2];
    if (pipe(pipefd) != 0) return "CRASH";
    pid_t pid = fork();
    if (pid < 0) { close(pipefd[0]); close(pipefd[1]); return "CRASH"; }

    if (pid == 0) {
        close(pipefd[0]);
        alarm(SPO_PER_CASE_TIMEOUT_SEC);
        Expr* call = expr_new_function(expr_new_symbol("Apply"),
            (Expr*[]){ expr_new_symbol("spoCheck"), expr_copy(rec) }, 2);
        Expr* v = evaluate(call);
        expr_free(call);
        char code = 'B';
        if (v && v->type == EXPR_STRING) {
            const char* s = v->data.string;
            if      (strcmp(s, "PASS") == 0)    code = 'P';
            else if (strcmp(s, "WRONG") == 0)   code = 'W';
            else if (strcmp(s, "UNEVAL") == 0)  code = 'U';
            else if (strcmp(s, "NOCHECK") == 0) code = 'N';
        }
        alarm(0);
        if (write(pipefd[1], &code, 1) != 1) _exit(3);
        close(pipefd[1]);
        _exit(0);
    }

    close(pipefd[1]);
    int status = 0;
    waitpid(pid, &status, 0);
    char code = 0;
    ssize_t r = read(pipefd[0], &code, 1);
    close(pipefd[0]);
    if (WIFSIGNALED(status) && WTERMSIG(status) == SIGALRM) return "TIMEOUT";
    if (r != 1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) return "CRASH";
    switch (code) {
        case 'P': return "PASS";
        case 'W': return "WRONG";
        case 'U': return "UNEVAL";
        case 'N': return "NOCHECK";
        default:  return "CRASH";
    }
}

int main(int argc, char** argv) {
    alarm(0);                       /* per-record alarms instead */
    setvbuf(stderr, NULL, _IONBF, 0);
    if (argc < 4) {
        fprintf(stderr, "usage: %s prelude.m corpus.m baseline.txt\n", argv[0]);
        return 2;
    }

    static SpoBase base[SPO_MAX_BASELINE];
    int nb = spo_read_baseline(argv[3], base, SPO_MAX_BASELINE);
    if (nb < 0) { fprintf(stderr, "FAIL: cannot read %s\n", argv[3]); return 1; }

    symtab_init();
    core_init();
    spo_get(argv[1], 0);
    Expr* cases = spo_get(argv[2], 1);
    if (!spo_is_list(cases)) {
        fprintf(stderr, "FAIL: %s did not load as a List\n", argv[2]);
        return 1;
    }

    /* SPO_VERBOSE=1 also lists every baselined record's current verdict,
     * to keep the baseline's annotations in step. */
    const char* venv = getenv("SPO_VERBOSE");
    int verbose = venv && *venv && strcmp(venv, "0") != 0;

    size_t n = cases->data.function.arg_count;
    int counts[6] = {0};            /* PASS WRONG UNEVAL NOCHECK TIMEOUT CRASH */
    const char* names[6] = {"PASS", "WRONG", "UNEVAL", "NOCHECK", "TIMEOUT", "CRASH"};
    int regressions = 0, fixed = 0, malformed = 0;

    fprintf(stderr, "==> Sum/Product oracle: %zu records, %d baselined\n", n, nb);
    for (size_t i = 0; i < n; i++) {
        Expr* rec = cases->data.function.args[i];
        if (!spo_is_list(rec) || rec->data.function.arg_count != 5
            || rec->data.function.args[0]->type != EXPR_STRING) {
            malformed++;
            continue;
        }
        const char* label = rec->data.function.args[0]->data.string;
        const char* v = spo_run_record(rec);
        for (int k = 0; k < 6; k++) if (strcmp(v, names[k]) == 0) counts[k]++;

        SpoBase* b = spo_find(base, nb, label);
        int pass = strcmp(v, "PASS") == 0;
        if (b) {
            b->seen = 1;
            if (!pass && verbose)
                fprintf(stderr, "  KNOWN   %-24s %s\n", label, v);
            if (pass) {
                fixed++;
                fprintf(stderr, "  FIXED   %-24s (delete it from the baseline)\n", label);
            }
        } else if (!pass) {
            regressions++;
            fprintf(stderr, "  REGRESS %-24s %s\n", label, v);
        }
    }
    expr_free(cases);

    int stale = 0;
    for (int i = 0; i < nb; i++)
        if (!base[i].seen) {
            stale++;
            fprintf(stderr, "  STALE   %-24s (baseline names no corpus record)\n",
                    base[i].label);
        }

    fprintf(stderr, "\n=== Sum/Product oracle ===\n");
    for (int k = 0; k < 6; k++) fprintf(stderr, "  %-8s %d\n", names[k], counts[k]);
    fprintf(stderr, "  baseline %d (fixed %d, stale %d), regressions %d, malformed %d\n",
            nb, fixed, stale, regressions, malformed);

    if (regressions || stale || malformed) {
        fprintf(stderr, "FAIL\n");
        return 1;
    }
    fprintf(stderr, "OK\n");
    return 0;
}
