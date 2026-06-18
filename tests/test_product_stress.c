/*
 * test_product_stress.c -- generated, oracle-verified stress families for the
 * symbolic Product (Stages 0-3), plus a curated regression corpus.
 *
 * Every generated case is checked by finite expansion: the symbolic closed form
 * Product[body, {k, lo, n}] is substituted at several integer n and compared to
 * the direct product Product[body, {k, lo, N}] (the multiplicative ground
 * truth).  No hand-written expected strings -- the oracle is the only judge.
 *
 * Run FOREGROUND only (redirect to /tmp); never background or poll.
 */

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;
static int checks = 0;
static int cases = 0;

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "PARSE FAIL: %s\n", input); exit(1); }
    Expr* v = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(v);
    expr_free(v);
    return s;
}

/* Oracle: compare Product[body,{k,lo,n}] /. n->N against Product[body,{k,lo,N}]
 * for N in {Ns...}.  `extra_cond`, if non-NULL, is an additional predicate
 * string (e.g. a FreeQ check) asserted == "True" on the symbolic closed form. */
static void oracle_body(const char* body, int lo, const int* Ns, int nN,
                        const char* extra_cond) {
    cases++;
    char closed[512];
    snprintf(closed, sizeof closed, "Product[%s, {k, %d, n}]", body, lo);

    if (extra_cond) {
        char cbuf[700];
        snprintf(cbuf, sizeof cbuf, "%s", extra_cond);
        char* s = eval_str(cbuf);
        checks++;
        if (strcmp(s, "True") != 0) {
            fprintf(stderr, "FAIL (cond): %s => %s\n", extra_cond, s);
            failures++;
        }
        free(s);
    }

    for (int i = 0; i < nN; i++) {
        int N = Ns[i];
        char buf[1400];
        snprintf(buf, sizeof buf,
                 "Simplify[((%s) /. n -> %d) - Product[%s, {k, %d, %d}]]",
                 closed, N, body, lo, N);
        char* s = eval_str(buf);
        checks++;
        if (strcmp(s, "0") != 0) {
            fprintf(stderr, "FAIL (oracle): %s @ n=%d  => %s\n", closed, N, s);
            failures++;
        }
        free(s);
    }
}

int main(void) {
    symtab_init();
    core_init();

    const int Ns[] = { 0, 1, 2, 4, 6, 8 };
    const int nN = (int)(sizeof Ns / sizeof Ns[0]);
    /* Avoid poles: lower bound 1, so k+b != 0 requires b not in {-1,-2,...}.
     * The shift values below keep every factor nonzero for k >= 1. */
    const char* shifts[] = { "0", "1", "2", "3", "1/2", "3/2", "5/2", "a" };
    const int nsh = (int)(sizeof shifts / sizeof shifts[0]);

    /* --- Family 1: rational  prod (k+a)/(k+b)  --> Pochhammer / telescoping --- */
    for (int ia = 0; ia < nsh; ia++) {
        for (int ib = 0; ib < nsh; ib++) {
            if (ia == ib) continue;
            char body[128];
            snprintf(body, sizeof body, "(k + %s)/(k + %s)", shifts[ia], shifts[ib]);
            oracle_body(body, 1, Ns, nN, NULL);
        }
    }

    /* --- Family 2: polynomial  prod (k+c1)(k+c2)  --> Pochhammer products --- */
    for (int ia = 0; ia < nsh; ia++) {
        for (int ib = ia; ib < nsh; ib++) {
            char body[128];
            snprintf(body, sizeof body, "(k + %s) (k + %s)", shifts[ia], shifts[ib]);
            oracle_body(body, 1, Ns, nN, NULL);
        }
    }

    /* --- Family 3: geometric  prod p(k) r^k  --> base^Sum * Pochhammer --- */
    {
        const char* rs[] = { "2", "3", "1/2", "a" };
        const char* ps[] = { "1", "k", "k^2", "(k + 1)", "k (k + 1)" };
        for (int ir = 0; ir < 4; ir++)
            for (int ip = 0; ip < 5; ip++) {
                char body[160];
                snprintf(body, sizeof body, "%s %s^k", ps[ip], rs[ir]);
                oracle_body(body, 1, Ns, nN, NULL);
            }
    }

    /* --- Family 4: telescoping ladders  prod h(k+1)/h(k)  (Gamma-free) --- */
    {
        const char* hs[] = { "k", "(k + 1)", "(k + 1) (k + 2)", "k (k + 2)",
                             "(k + 1/2)", "(2 k + 1)" };
        for (int ih = 0; ih < 6; ih++) {
            char body[200];
            /* (h /. k -> k+1) / h  -- a clean telescoping multiplicand */
            snprintf(body, sizeof body, "(%s /. k -> k + 1)/(%s)", hs[ih], hs[ih]);
            /* Telescoping output must be free of Gamma AND Pochhammer. */
            char cond[700];
            snprintf(cond, sizeof cond,
                "FreeQ[Product[%s, {k, 1, n}, Method -> \"Telescoping\"], Gamma] && "
                "FreeQ[Product[%s, {k, 1, n}, Method -> \"Telescoping\"], Pochhammer]",
                body, body);
            oracle_body(body, 1, Ns, nN, cond);
        }
    }

    /* --- Regression corpus: known shapes + edge cases --- */
    /* empty / reversed ranges -> 1 */
    { char* s = eval_str("Product[k, {k, 1, 0}]"); checks++; cases++;
      if (strcmp(s, "1")) { fprintf(stderr, "FAIL empty\n"); failures++; } free(s); }
    { char* s = eval_str("Product[k^3, {k, 7, 2}]"); checks++; cases++;
      if (strcmp(s, "1")) { fprintf(stderr, "FAIL reversed\n"); failures++; } free(s); }
    /* half-integer shifts */
    oracle_body("(2 k - 1)/(2 k)", 1, Ns, nN, NULL);
    oracle_body("(k - 1/2)", 1, Ns, nN, NULL);
    /* geometric x rational mixes */
    oracle_body("k^3 2^k", 1, Ns, nN, NULL);
    oracle_body("(k + 1)/(k + 3) 3^k", 1, Ns, nN, NULL);

    /* symbolic-exponent hang guard: must return promptly and not hang. */
    { char* s = eval_str("FreeQ[Product[a^(k^2), {k, 1, n}], Product]");
      checks++; cases++;
      if (strcmp(s, "True")) { fprintf(stderr, "FAIL hang-guard a^(k^2)\n"); failures++; }
      free(s); }

    printf("\n%d cases, %d checks, %d failures\n", cases, checks, failures);
    if (failures) { fprintf(stderr, "STRESS FAILED\n"); return 1; }
    printf("All product stress checks passed!\n");
    return 0;
}
