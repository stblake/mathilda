/* test_ratcanon_build.c — Phase 2 of RATCANON_REWRITE_PLAN.md.
 *
 * Validates rat_canon_build: it turns a rational expression into one ordered
 * generator tower (transcendental free, algebraic with an explicit relation),
 * with num/den a plain rational function over the generators.  Tests:
 *   - structure: generator counts and kinds;
 *   - relations: every algebraic gen's relation vanishes at its kernel, every
 *     transcendental gen has no relation;
 *   - independence: commensurate exponentials collapse to one fundamental,
 *     Log[x^2] expands to Log[x];
 *   - round-trip: substituting the generators back reproduces the input value.
 */
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "sym_names.h"
#include "ratcanon.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static RatCanonForm* build_from(const char* input) {
    Expr* e = parse_expression(input);
    Expr* v = evaluate(e);
    RatCanonForm* f = rat_canon_build(v);
    expr_free(e); expr_free(v);
    return f;
}

static int count_kind(const RatCanonForm* f, RcGenKind k) {
    int n = 0;
    for (size_t i = 0; i < f->n; i++) if (f->gens[i].kind == k) n++;
    return n;
}

/* Generator counts. */
static void gens(const char* input, int n_alg, int n_trans) {
    RatCanonForm* f = build_from(input);
    int a = count_kind(f, RCG_ALGEBRAIC), t = count_kind(f, RCG_TRANSCENDENTAL);
    if (a != n_alg || t != n_trans) {
        printf("FAIL gens: %s\n  expected alg=%d trans=%d, got alg=%d trans=%d\n",
               input, n_alg, n_trans, a, t);
    }
    ASSERT(a == n_alg);
    ASSERT(t == n_trans);
    /* algebraic generators must be ordered first (LEX-leading). */
    int seen_trans = 0;
    for (size_t i = 0; i < f->n; i++) {
        if (f->gens[i].kind == RCG_TRANSCENDENTAL) seen_trans = 1;
        else if (seen_trans) { printf("FAIL order: %s (alg after trans)\n", input); ASSERT(0); }
    }
    rat_canon_free(f);
}

/* Invariants: every algebraic relation vanishes at its kernel; every
 * transcendental gen has relation == NULL. */
static void relations_ok(const char* input) {
    RatCanonForm* f = build_from(input);
    for (size_t i = 0; i < f->n; i++) {
        if (f->gens[i].kind == RCG_TRANSCENDENTAL) {
            if (f->gens[i].relation) { printf("FAIL: %s trans gen has relation\n", input); ASSERT(0); }
            continue;
        }
        ASSERT(f->gens[i].relation != NULL);
        /* relation[sym -> kernel] must evaluate to 0 */
        Expr* subst = rat_canon_subst_back(f, f->gens[i].relation);
        Expr* val = evaluate(subst);
        int zero = (val->type == EXPR_INTEGER && val->data.integer == 0);
        if (!zero) {
            char* s = expr_to_string(val);
            printf("FAIL relation: %s gen %zu relation != 0 at kernel: %s\n", input, i, s);
            free(s);
        }
        ASSERT(zero);
        expr_free(subst); expr_free(val);
    }
    rat_canon_free(f);
}

/* Round-trip: Together[roundtrip(build(e)) - e] == 0. */
static void roundtrip(const char* input) {
    Expr* e = parse_expression(input);
    Expr* v = evaluate(e);
    RatCanonForm* f = rat_canon_build(v);
    Expr* rt = rat_canon_roundtrip(f);
    /* diff = Together[rt - v] */
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ expr_new_integer(-1), expr_copy(v) }, 2);
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus),
                  (Expr*[]){ rt, neg }, 2);
    Expr* tg = expr_new_function(expr_new_symbol(SYM_Together), (Expr*[]){ sum }, 1);
    Expr* diff = evaluate(tg);
    int zero = (diff->type == EXPR_INTEGER && diff->data.integer == 0);
    if (!zero) {
        char* s = expr_to_string(diff);
        printf("FAIL roundtrip: %s -> nonzero diff %s\n", input, s);
        free(s);
    }
    ASSERT(zero);
    expr_free(e); expr_free(v); expr_free(tg); expr_free(diff);
    rat_canon_free(f);
}

/* -------------------------------------------------------------------------- */

static void test_structure(void) {
    gens("(1+x)/(1-x^2)", 0, 0);                       /* plain Q            */
    gens("1/(1+Log[x])+1/Log[x]", 0, 1);               /* one Log            */
    gens("Log[x]+Log[y]", 0, 2);                       /* two independent    */
    gens("1/(x-I)+1/(x+I)", 1, 0);                     /* I                  */
    gens("1/(b-a Sqrt[k])+1/(b+a Sqrt[k])", 1, 0);     /* Sqrt[k]            */
    gens("1/(x-Sqrt[2])+1/(x+Sqrt[2])", 1, 0);         /* Sqrt[2]            */
    gens("1/(x-Sqrt[3])+1/(x-Sqrt[5])", 2, 0);         /* Sqrt3, Sqrt5       */
    gens("1/(x-(-1)^(1/3))", 1, 0);                    /* root of unity      */
    gens("1/(1+Tan[x])+1/(1-Tan[x])", 0, 1);           /* Tan                */
    gens("1/(1+ArcTan[x])", 0, 1);                     /* inverse-trig       */
    gens("Log[x]/(x - Sqrt[2])", 1, 1);                /* mixed alg + trans  */
}

static void test_independence(void) {
    /* commensurate exponentials collapse to one fundamental */
    gens("(E^(2x)-1)/(E^x-1)", 0, 1);
    gens("1/(E^x-1)+1/(E^(3x)-1)", 0, 1);
    /* distinct exp bases stay independent */
    gens("E^x + E^y", 0, 2);
    /* Log[x^2] expands to Log[x] (one gen) */
    gens("1/(1+Log[x^2])+1/Log[x]", 0, 1);
}

static void test_relations(void) {
    relations_ok("1/(x-I)+1/(x+I)");
    relations_ok("1/(b-a Sqrt[k])+1/(b+a Sqrt[k])");
    relations_ok("1/(x-Sqrt[2])+1/(x+Sqrt[2])");
    relations_ok("1/(x-Sqrt[3])+1/(x-Sqrt[5])");
    relations_ok("1/(x-(-1)^(1/3))");
    relations_ok("(x^2-2)/(x-Sqrt[2])");
    relations_ok("1/(y^(1/3)-1)");
    relations_ok("(E^(2x)-1)/(E^x-1)");
    relations_ok("1/(1+Log[x])+1/Log[x]");
}

static void test_roundtrip(void) {
    roundtrip("(1+x)/(1-x^2)");
    roundtrip("1/(1+Log[x])+1/Log[x]");
    roundtrip("1/(x-I)+1/(x+I)");
    roundtrip("1/(b-a Sqrt[k])+1/(b+a Sqrt[k])");
    roundtrip("(E^(2x)-1)/(E^x-1)");
    roundtrip("1/(x-Sqrt[2])+1/(x+Sqrt[2])");
    roundtrip("1/(x-Sqrt[3])+1/(x-Sqrt[5])");
    roundtrip("1/(x-(-1)^(1/3))+1/(x+(-1)^(1/3))");
    roundtrip("1/(1+Tan[x])+1/(1-Tan[x])");
    roundtrip("Log[x]/(x-Sqrt[2])+1/(x+Sqrt[2])");
    roundtrip("(x^2-2)/(x-Sqrt[2])");
    roundtrip("a/b+c/d");
    roundtrip("1/(E^x-1)+1/(E^(3x)-1)");
}

int main(void) {
    symtab_init();
    core_init();
    printf("Running ratcanon build (Phase 2) tests...\n");
    TEST(test_structure);
    TEST(test_independence);
    TEST(test_relations);
    TEST(test_roundtrip);
    printf("All ratcanon build tests passed!\n");
    return 0;
}
