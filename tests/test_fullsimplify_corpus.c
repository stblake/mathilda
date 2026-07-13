/*
 * test_fullsimplify_corpus.c -- per-rule corpus runner for FullSimplify.
 *
 * Loads tests/fullsimplify_corpus.m (a List of {input, expected} pairs) and
 * checks that FullSimplify[input] is structurally equal to expected for every
 * pair. Modeled on test_crc_corpus.c but in-process (the cases are small and
 * fast, and the transform rules are non-recursive), so no per-case fork is
 * needed. The corpus path may be overridden on the command line; it defaults
 * to ../fullsimplify_corpus.m (the layout used when running from tests/build).
 *
 * The kernel's FullSimplify driver is loaded by core_init (it resolves
 * src/internal/simp/FullSimplify.m independently of the working directory), so
 * no extra setup is required here.
 */

#include "expr.h"
#include "symtab.h"
#include "parse.h"
#include "eval.h"
#include "print.h"
#include "context.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void core_init(void);

static int is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, "List") == 0;
}

int main(int argc, char** argv) {
    const char* corpus_file =
        (argc > 1) ? argv[1] : "../fullsimplify_corpus.m";

    symtab_init();
    core_init();

    char get_buf[1024];
    snprintf(get_buf, sizeof(get_buf), "Get[\"%s\"]", corpus_file);
    Expr* get_call = parse_expression(get_buf);
    if (!get_call) {
        fprintf(stderr, "FAIL: could not parse Get for %s\n", corpus_file);
        return 1;
    }
    Expr* tests = evaluate(get_call);
    expr_free(get_call);

    if (!is_list(tests)) {
        fprintf(stderr, "FAIL: could not load %s as a List.\n", corpus_file);
        if (tests) expr_free(tests);
        return 1;
    }

    size_t n = tests->data.function.arg_count;
    fprintf(stderr, "==> FullSimplify corpus: %zu cases\n", n);

    int passed = 0, failed = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* pair = tests->data.function.args[i];
        if (!is_list(pair) || pair->data.function.arg_count < 2) {
            fprintf(stderr, "  [%3zu/%zu] MALFORMED\n", i + 1, n);
            failed++;
            continue;
        }
        Expr* input    = pair->data.function.args[0];
        Expr* expected = pair->data.function.args[1];

        /* Build FullSimplify[input] from a copy and evaluate it. The head name
         * is resolved through the context system exactly as the parser would,
         * so the bare "FullSimplify" maps to the package's FullSimplify`
         * FullSimplify symbol (where the DownValue lives) via $ContextPath. */
        Expr** args = (Expr**)malloc(sizeof(Expr*));
        args[0] = expr_copy(input);
        char* head_name = context_resolve_name("FullSimplify");
        Expr* call = expr_new_function(expr_new_symbol(head_name), args, 1);
        free(head_name);
        free(args);
        Expr* result = evaluate(call);
        expr_free(call); /* evaluate does not free its argument */

        if (expr_eq(result, expected)) {
            passed++;
        } else {
            failed++;
            char* in_s  = expr_to_string(input);
            char* ex_s  = expr_to_string(expected);
            char* rs_s  = expr_to_string(result);
            fprintf(stderr,
                "  [%3zu/%zu] FAIL: FullSimplify[%s]\n"
                "        expected: %s\n"
                "        actual:   %s\n",
                i + 1, n, in_s, ex_s, rs_s);
            free(in_s); free(ex_s); free(rs_s);
        }
        expr_free(result);
    }

    expr_free(tests);
    fprintf(stderr, "==> FullSimplify corpus: %d passed, %d failed\n",
            passed, failed);
    return failed == 0 ? 0 : 1;
}
