/* test_roach_corpus.c
 *
 * Corpus runner for tests/RoachCorpus.m: the definite integrals and pFq
 * special values / representations of Kelly Roach's 2008 Quantum paper and
 * his ISSAC '96 and '97 galleries (see the corpus header for the list).
 *
 * Each record is  {id, Hold[expr], {param -> val, ...}, reference, {heads}}.
 * For every record, in a forked child with a per-case timeout:
 *
 *   1. r = ReleaseHold[expr];
 *   2. r still contains one of the record's uneval heads (Integrate; the
 *      pFq heads for a FunctionExpand record)            -> UNEVAL
 *   3. ConditionalExpression[v, _] -> v; v = N[r /. params] (machine, else
 *      30 digits); |v - reference| <= 1e-9 (1 + |reference|) -> CORRECT
 *   4. a finite value that misses the reference, or an Indeterminate /
 *      infinite one                                        -> WRONG
 *      a value with symbolic leftovers N cannot evaluate   -> UNEVAL
 *   5. the child outlives the budget                        -> TIMEOUT
 *      the child dies                                      -> CRASH
 *
 * The references come from tools/roach_corpus_refs.py (mpmath, 40 digits,
 * from the definition of each quantity), never from Mathilda.
 *
 * GATES. WRONG is assert-empty. CORRECT ratchets against a baseline file of
 * record ids (tests/RoachCorpus.baseline, the second argument): an id listed
 * there that is no longer CORRECT fails the run; an id CORRECT now and not
 * listed is reported so that the line can be added. With the environment
 * variable ROACH_CORPUS_WRITE_BASELINE=1 the runner rewrites the baseline
 * with the ids CORRECT in this run instead of gating on it.
 *
 * Modeled on test_extreme_definite_corpus.c.
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

#define ROACH_PER_CASE_TIMEOUT_SEC 30
#define ROACH_MAX_CASES 512
#define ROACH_ID_LEN 32

/* One-letter verdicts sent from the child over a pipe. */
enum {
    V_CORRECT = 'C',
    V_WRONG   = 'W',
    V_UNEVAL  = 'U',
};

/* The scorer, in the language: the comparison is numeric and Mathilda's N
 * does it better than a C walk over Real/Complex/MPFR nodes. A symbolic
 * leftover N cannot evaluate is UNEVAL ("N"); a non-finite value is WRONG. */
static const char* SCORER =
    "RoachCorpus`num[e_] := Module[{v = Quiet[N[e]]},"
    "  If[!NumericQ[v], v = Quiet[N[N[e, 30]]]];"
    "  If[NumericQ[v], {Re[v], Im[v]}, $Failed]];"
    "RoachCorpus`score[h_, inst_, ref_, heads_] := Module[{r, v, w},"
    "  r = ReleaseHold[h];"
    "  Print[\"      result: \", StringTake[ToString[r, InputForm], UpTo[240]]];"
    "  If[!FreeQ[r, Alternatives @@ (ToExpression /@ heads)], Return[\"U\"]];"
    "  r = r //. ConditionalExpression[c_, _] :> c;"
    "  If[!FreeQ[r, Indeterminate | ComplexInfinity | DirectedInfinity], Return[\"W\"]];"
    "  v = RoachCorpus`num[r /. inst]; w = RoachCorpus`num[ref];"
    "  If[v === $Failed, Return[\"N\"]];"
    "  If[Abs[(v[[1]] - w[[1]]) + I (v[[2]] - w[[2]])] <= 10.^-9 (1 + Abs[w[[1]] + I w[[2]]]),"
    "     \"C\", Print[\"      value: \", v, \"  reference: \", w]; \"W\"]]";

static bool is_list_expr(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, "List") == 0;
}

/* ---- baseline: one record id per line, '#' starts a comment ---- */

static int load_baseline(const char* path, char ids[][ROACH_ID_LEN], int max) {
    FILE* f = path ? fopen(path, "r") : NULL;
    if (!f) return -1;
    int n = 0;
    char line[256];
    while (fgets(line, sizeof line, f) && n < max) {
        char* p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\0') continue;
        size_t k = strcspn(p, " \t\r\n#");
        if (k == 0) continue;
        if (k >= ROACH_ID_LEN) k = ROACH_ID_LEN - 1;
        memcpy(ids[n], p, k);
        ids[n][k] = '\0';
        n++;
    }
    fclose(f);
    return n;
}

static bool id_in(const char* id, char ids[][ROACH_ID_LEN], int n) {
    for (int i = 0; i < n; i++)
        if (strcmp(ids[i], id) == 0) return true;
    return false;
}

int main(int argc, char** argv) {
    alarm(0);
    setvbuf(stderr, NULL, _IONBF, 0);
    setvbuf(stdout, NULL, _IONBF, 0);

    const char* corpus_file = (argc > 1) ? argv[1] : "../RoachCorpus.m";
    const char* baseline_file = (argc > 2) ? argv[2] : NULL;
    const char* wb = getenv("ROACH_CORPUS_WRITE_BASELINE");
    bool write_baseline = wb && strcmp(wb, "1") == 0;

    symtab_init();
    core_init();
    /* test_load_init_m chdir's to the repo root to find init.m; come back,
     * so that the corpus and baseline paths are read relative to where the
     * runner was started. */
    char start_dir[4096];
    bool have_dir = getcwd(start_dir, sizeof start_dir) != NULL;
    test_load_init_m();          /* the production environment, as -file sees it */
    if (have_dir) (void)!chdir(start_dir);
    test_mute_messages();

    Expr* sc = parse_expression(SCORER);
    ASSERT(sc != NULL);
    Expr* scr = evaluate(sc);
    expr_free(sc);
    if (scr) expr_free(scr);

    fprintf(stderr, "==> corpus file: %s\n", corpus_file);
    char get_buf[1024];
    snprintf(get_buf, sizeof(get_buf), "Get[\"%s\"]", corpus_file);
    Expr* get_call = parse_expression(get_buf);
    ASSERT(get_call != NULL);
    Expr* tests = evaluate(get_call);
    expr_free(get_call);
    if (!is_list_expr(tests)) {
        fprintf(stderr, "FAIL: could not load %s as a List.\n", corpus_file);
        if (tests) expr_free(tests);
        ASSERT(false);
        return 1;
    }

    size_t n = tests->data.function.arg_count;
    if (n > ROACH_MAX_CASES) n = ROACH_MAX_CASES;
    fprintf(stderr, "==> Roach corpus: %zu cases\n", n);

    static char correct_ids[ROACH_MAX_CASES][ROACH_ID_LEN];
    int correct = 0, wrong = 0, uneval = 0, timed_out = 0, crashed = 0, malformed = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* rec = tests->data.function.args[i];
        if (!is_list_expr(rec) || rec->data.function.arg_count != 5 ||
            rec->data.function.args[0]->type != EXPR_STRING) {
            malformed++;
            fprintf(stderr, "  [%3zu/%zu] MALFORMED\n", i + 1, n);
            continue;
        }
        const char* id = rec->data.function.args[0]->data.string;
        char* in_str = expr_to_string(rec->data.function.args[1]);
        fprintf(stderr, "  [%3zu/%zu] %s  %s\n", i + 1, n, id, in_str);
        free(in_str);

        int pipefd[2];
        if (pipe(pipefd) != 0) { crashed++; continue; }
        fflush(stdout);
        pid_t pid = fork();
        if (pid < 0) { close(pipefd[0]); close(pipefd[1]); crashed++; continue; }

        if (pid == 0) {
            /* --- child --- */
            close(pipefd[0]);
            alarm(ROACH_PER_CASE_TIMEOUT_SEC);
            Expr* call = expr_new_function(expr_new_symbol("RoachCorpus`score"),
                (Expr*[]){ expr_copy(rec->data.function.args[1]),
                           expr_copy(rec->data.function.args[2]),
                           expr_copy(rec->data.function.args[3]),
                           expr_copy(rec->data.function.args[4]) }, 4);
            Expr* r = evaluate(call);
            expr_free(call);
            alarm(0);
            char code = V_UNEVAL;
            if (r && r->type == EXPR_STRING) {
                switch (r->data.string[0]) {
                    case 'C': code = V_CORRECT; break;
                    case 'W': code = V_WRONG; break;
                    case 'N':
                        fprintf(stderr, "      (no numeric value)\n");
                        code = V_UNEVAL; break;
                    default: code = V_UNEVAL; break;
                }
            }
            if (r) expr_free(r);
            ssize_t w = write(pipefd[1], &code, 1); (void)w;
            close(pipefd[1]);
            _exit(0);
        }

        /* --- parent --- */
        close(pipefd[1]);
        time_t deadline = time(NULL) + ROACH_PER_CASE_TIMEOUT_SEC + 3;
        int status = 0;
        bool reaped = false;
        while (time(NULL) < deadline) {
            pid_t rr = waitpid(pid, &status, WNOHANG);
            if (rr == pid) { reaped = true; break; }
            if (rr < 0) break;
            usleep(20000);
        }
        if (!reaped) {
            kill(pid, SIGKILL);
            waitpid(pid, &status, 0);
            close(pipefd[0]);
            timed_out++;
            fprintf(stderr, "      -> TIMEOUT (>%ds)\n", ROACH_PER_CASE_TIMEOUT_SEC);
            continue;
        }
        char code = 0;
        ssize_t rd = read(pipefd[0], &code, 1);
        close(pipefd[0]);
        if (rd != 1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
            if (WIFSIGNALED(status) && WTERMSIG(status) == SIGALRM) {
                timed_out++;
                fprintf(stderr, "      -> TIMEOUT (child alarm)\n");
                continue;
            }
            crashed++;
            fprintf(stderr, "      -> CRASH (exit %d, sig %d)\n",
                    WIFEXITED(status) ? WEXITSTATUS(status) : -1,
                    WIFSIGNALED(status) ? WTERMSIG(status) : 0);
            continue;
        }
        switch (code) {
            case V_CORRECT:
                snprintf(correct_ids[correct], ROACH_ID_LEN, "%s", id);
                correct++;
                fprintf(stderr, "      -> CORRECT\n");
                break;
            case V_WRONG:
                wrong++;
                fprintf(stderr, "      -> WRONG (%s: value != reference)\n", id);
                break;
            default:
                uneval++;
                fprintf(stderr, "      -> UNEVAL\n");
                break;
        }
    }
    expr_free(tests);

    fprintf(stderr, "\n=== Roach corpus result ===\n");
    fprintf(stderr, "  Total cases:   %zu\n", n);
    fprintf(stderr, "  CORRECT:       %d\n", correct);
    fprintf(stderr, "  WRONG:         %d\n", wrong);
    fprintf(stderr, "  UNEVAL:        %d\n", uneval);
    fprintf(stderr, "  TIMEOUT (>%ds): %d\n", ROACH_PER_CASE_TIMEOUT_SEC, timed_out);
    fprintf(stderr, "  CRASH:         %d\n", crashed);
    fprintf(stderr, "  MALFORMED:     %d\n", malformed);
    fprintf(stderr, "===========================\n");

    int fail = 0;
    if (wrong > 0) {
        fprintf(stderr, "FAIL: %d case(s) returned a WRONG value (WRONG is assert-empty).\n", wrong);
        fail = 1;
    }
    if (malformed > 0) {
        fprintf(stderr, "FAIL: %d malformed record(s).\n", malformed);
        fail = 1;
    }

    if (write_baseline && baseline_file) {
        FILE* f = fopen(baseline_file, "w");
        if (!f) { fprintf(stderr, "FAIL: cannot write %s\n", baseline_file); return 1; }
        fprintf(f, "# tests/RoachCorpus.baseline -- the record ids of tests/RoachCorpus.m\n"
                   "# that are CORRECT. The runner (test_roach_corpus.c) fails when one of\n"
                   "# them is no longer CORRECT; add an id when a case starts passing.\n"
                   "# Regenerate with ROACH_CORPUS_WRITE_BASELINE=1.\n");
        for (int k = 0; k < correct; k++) fprintf(f, "%s\n", correct_ids[k]);
        fclose(f);
        fprintf(stderr, "==> wrote %d ids to %s\n", correct, baseline_file);
        return fail;
    }

    static char base[ROACH_MAX_CASES][ROACH_ID_LEN];
    int nb = load_baseline(baseline_file, base, ROACH_MAX_CASES);
    if (nb < 0) {
        fprintf(stderr, "note: no baseline file; CORRECT is not gated.\n");
        return fail;
    }
    int lost = 0, gained = 0;
    for (int k = 0; k < nb; k++) {
        if (!id_in(base[k], correct_ids, correct)) {
            fprintf(stderr, "FAIL: %s is in the baseline but no longer CORRECT.\n", base[k]);
            lost++;
        }
    }
    for (int k = 0; k < correct; k++) {
        if (!id_in(correct_ids[k], base, nb)) {
            fprintf(stderr, "note: %s is now CORRECT -- add it to %s.\n",
                    correct_ids[k], baseline_file);
            gained++;
        }
    }
    fprintf(stderr, "==> baseline: %d ids, %d lost, %d gained\n", nb, lost, gained);
    if (lost > 0) fail = 1;
    return fail;
}
