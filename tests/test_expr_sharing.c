/* test_expr_sharing.c -- M3 phase-2 (function-node sharing) stress tests.
 *
 * The Phase-2 invariant: expr_copy() inc-refs ANY node (atom or
 * FUNCTION), so callers receive a shared pointer. Mutators must call
 * expr_unshare() before writing to fields, or work on freshly
 * constructed (refcount==1) nodes only.
 *
 * These tests directly probe the refcount mechanics, the unshare
 * semantics, and the safety of the audited mutation hot-paths in
 * eval.c, print.c, plus.c, times.c, core.c, match.c, etc.
 */

#include "core.h"
#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>

/* ------------------------------------------------------------------ */
/* 1. Refcount mechanics                                               */
/* ------------------------------------------------------------------ */

static void test_atom_refcount_inc_dec(void) {
    Expr* a = expr_new_integer(42);
    ASSERT(a->refcount == 1);

    Expr* b = expr_copy(a);
    ASSERT(a == b);                /* atom shares the same pointer */
    ASSERT(a->refcount == 2);

    Expr* c = expr_copy(a);
    ASSERT(c == a);
    ASSERT(a->refcount == 3);

    expr_free(b);
    ASSERT(a->refcount == 2);
    expr_free(c);
    ASSERT(a->refcount == 1);
    expr_free(a);                  /* finally physically frees */
}

static void test_function_refcount_inc_dec(void) {
    Expr* args[2] = { expr_new_integer(2), expr_new_integer(3) };
    Expr* f = expr_new_function(expr_new_symbol("Plus"), args, 2);
    ASSERT(f->refcount == 1);

    Expr* g = expr_copy(f);
    ASSERT(g == f);                /* Phase-2: FUNCTION shares pointer */
    ASSERT(f->refcount == 2);

    /* Children stay singly-owned by the (still single) function node;
     * sharing the parent does NOT bump children's refcount. */
    ASSERT(f->data.function.args[0]->refcount == 1);
    ASSERT(f->data.function.args[1]->refcount == 1);
    ASSERT(f->data.function.head->refcount == 1);

    expr_free(g);
    ASSERT(f->refcount == 1);
    expr_free(f);
}

static void test_atom_share_pointer_identity(void) {
    /* Across many copies and frees, the pointer stays stable. */
    Expr* a = expr_new_real(3.14);
    Expr* refs[100];
    for (int i = 0; i < 100; i++) {
        refs[i] = expr_copy(a);
        ASSERT(refs[i] == a);
    }
    ASSERT(a->refcount == 101);
    for (int i = 99; i >= 0; i--) expr_free(refs[i]);
    ASSERT(a->refcount == 1);
    expr_free(a);
}

/* ------------------------------------------------------------------ */
/* 2. expr_unshare semantics                                           */
/* ------------------------------------------------------------------ */

static void test_unshare_unique_is_noop(void) {
    Expr* a = expr_new_integer(7);
    ASSERT(a->refcount == 1);
    Expr* b = expr_unshare(a);
    ASSERT(b == a);                /* fast path: same pointer */
    ASSERT(b->refcount == 1);
    expr_free(b);
}

static void test_unshare_atom_creates_private(void) {
    Expr* a = expr_new_integer(7);
    Expr* shared = expr_copy(a);   /* a->refcount == 2 */
    ASSERT(a->refcount == 2);

    /* unshare consumes one ref of `shared`. Since refcount > 1, it
     * allocates a fresh private copy and dec-refs the original. */
    Expr* priv = expr_unshare(shared);
    ASSERT(priv != a);             /* different node */
    ASSERT(priv->refcount == 1);
    ASSERT(priv->data.integer == 7);
    ASSERT(a->refcount == 1);      /* original now back to refcount 1 */

    expr_free(priv);
    expr_free(a);
}

static void test_unshare_function_creates_private_args_array(void) {
    Expr* args[2] = { expr_new_integer(10), expr_new_integer(20) };
    Expr* f = expr_new_function(expr_new_symbol("Plus"), args, 2);

    Expr* shared = expr_copy(f);   /* f->refcount == 2 */
    ASSERT(f->refcount == 2);

    Expr* priv = expr_unshare(shared);
    ASSERT(priv != f);             /* new top-level node */
    ASSERT(priv->refcount == 1);
    ASSERT(f->refcount == 1);      /* dec'd by unshare */

    /* Args are shared (inc-ref'd from `f`), not deep-copied. */
    ASSERT(priv->data.function.args[0] == f->data.function.args[0]);
    ASSERT(f->data.function.args[0]->refcount == 2);

    /* But the args ARRAY is private — different malloc. */
    ASSERT(priv->data.function.args != f->data.function.args);

    /* Mutating priv's args[0] slot must NOT affect f. */
    expr_free(priv->data.function.args[0]);
    priv->data.function.args[0] = expr_new_integer(99);

    /* f sees its original 10 unchanged. */
    ASSERT(f->data.function.args[0]->type == EXPR_INTEGER);
    ASSERT(f->data.function.args[0]->data.integer == 10);

    expr_free(priv);
    expr_free(f);
}

static void test_unshare_bigint_payload_independent(void) {
    Expr* a = expr_new_bigint_from_str("99999999999999999999");
    Expr* shared = expr_copy(a);
    ASSERT(a->refcount == 2);

    Expr* priv = expr_unshare(shared);
    ASSERT(priv != a);
    ASSERT(priv->refcount == 1);
    ASSERT(a->refcount == 1);

    /* Each owns its own mpz_t. Mutate `a`'s mpz to verify priv unaffected
     * (we never do this in production, but the test validates the
     * payload is genuinely independent). */
    mpz_set_si(a->data.bigint, 0);
    ASSERT(mpz_cmp_si(priv->data.bigint, 0) != 0);

    expr_free(priv);
    expr_free(a);
}

static void test_unshare_string_payload_independent(void) {
    Expr* a = expr_new_string("hello");
    Expr* shared = expr_copy(a);
    ASSERT(a->refcount == 2);

    Expr* priv = expr_unshare(shared);
    ASSERT(priv != a);
    ASSERT(priv->data.string != a->data.string);
    ASSERT(strcmp(priv->data.string, "hello") == 0);

    expr_free(priv);
    expr_free(a);
}

/* ------------------------------------------------------------------ */
/* 3. End-to-end correctness on the audited hot paths                  */
/* ------------------------------------------------------------------ */

/* Plus does in-place numeric contagion (plus.c:233 args[i] = numed[i])
 * on the freshly-built `res`. Verify the result is correct AND repeats
 * are stable (no aliasing accumulating across runs). */
static void test_plus_numeric_contagion_repeated(void) {
    for (int i = 0; i < 200; i++) {
        assert_eval_eq("1.0 + Pi", "4.14159", 0);
    }
}

/* Times does the same args[i] mutation. */
static void test_times_numeric_contagion_repeated(void) {
    for (int i = 0; i < 200; i++) {
        assert_eval_eq("2.0 * Pi", "6.28319", 0);
    }
}

/* eval_flatten_args mutates res->data.function.args (eval.c:191-192).
 * Verify Flat works under heavy repetition. */
static void test_flat_flatten_repeated(void) {
    for (int i = 0; i < 200; i++) {
        assert_eval_eq("Plus[1, Plus[2, Plus[3, 4]], 5]", "15", 0);
    }
}

/* Orderless triggers qsort on res->data.function.args (eval.c:588). */
static void test_orderless_sort_repeated(void) {
    for (int i = 0; i < 200; i++) {
        assert_eval_eq("Plus[c, b, a, d]", "a + b + c + d", 0);
    }
}

/* flatten_sequences mutates res->data.function.args (eval.c:445-446). */
static void test_flatten_sequences_repeated(void) {
    for (int i = 0; i < 100; i++) {
        assert_eval_eq("Plus[1, Sequence[2, 3], 4]", "10", 0);
        assert_eval_eq("List[a, Sequence[b, c], d]", "{a, b, c, d}", 0);
    }
}

/* QuotientRemainder zeroes arg_count on intermediate FUNCTION nodes
 * (core.c:1498-1499). Verify correctness under repetition. */
static void test_quotient_remainder_repeated(void) {
    for (int i = 0; i < 100; i++) {
        assert_eval_eq("QuotientRemainder[17, 5]", "{3, 2}", 0);
        assert_eval_eq("QuotientRemainder[100, 7]", "{14, 2}", 0);
    }
}

/* print.c mutates `t_copy = expr_unshare(expr_copy(arg))` to negate
 * literal-prefix coefficients before printing. Verify the original
 * expression is not corrupted. */
static void test_print_does_not_mutate_input(void) {
    /* Build an expression "Times[-2, x]" by parsing then re-print
     * many times. The print path used to mutate args[0] in place;
     * we now unshare. After many print calls, the parsed tree must
     * still reproduce the same string. */
    Expr* parsed = parse_expression("-2 x");
    ASSERT(parsed != NULL);
    char* expected = expr_to_string(parsed);
    for (int i = 0; i < 200; i++) {
        char* s = expr_to_string(parsed);
        ASSERT_STR_EQ(s, expected);
        free(s);
    }
    free(expected);
    expr_free(parsed);
}

static void test_print_negative_rational_does_not_mutate(void) {
    /* Times[Rational[-1, 2], x] prints as "-x/2". The print path
     * unshares two levels (t_copy AND the inner Rational). */
    Expr* parsed = parse_expression("-x/2");
    ASSERT(parsed != NULL);
    Expr* eval = evaluate(parsed);
    expr_free(parsed);

    char* expected = expr_to_string(eval);
    for (int i = 0; i < 200; i++) {
        char* s = expr_to_string(eval);
        ASSERT_STR_EQ(s, expected);
        free(s);
    }
    free(expected);
    expr_free(eval);
}

/* ------------------------------------------------------------------ */
/* 4. Pattern matching with shared bindings                            */
/* ------------------------------------------------------------------ */

/* When a binding variable appears multiple times in the RHS, the
 * substitute path expr_copy()'s the binding once per slot. Phase-2
 * makes those copies share. The result tree contains multiple aliased
 * pointers; later evaluation (e.g., Plus combining identical terms)
 * must not corrupt anything. */
static void test_repeated_binding_substitution(void) {
    /* x_ binds to (a + b) which is itself a FUNCTION node.
     * The RHS produces (a + b) * (a + b). */
    for (int i = 0; i < 50; i++) {
        assert_eval_eq("(a + b) /. x_ -> x * x", "(a + b)^2", 0);
    }
}

static void test_pattern_binding_evaluates_shared_subtree(void) {
    /* The binding for `expr` is a shared subtree. After substitution,
     * the result is `f[expr, expr]` -- two pointers into the same
     * tree. Evaluating this still has to work. */
    for (int i = 0; i < 50; i++) {
        assert_eval_eq("(1 + 2) /. x_ -> {x, x, x}", "{3, 3, 3}", 0);
    }
}

/* ------------------------------------------------------------------ */
/* 5. Heavy use-the-system stress tests                                */
/* ------------------------------------------------------------------ */

static void test_repeated_complex_evaluation(void) {
    /* Expand triggers many subexpression copies and rewrites.
     * Run many times; any aliasing bug typically surfaces as a
     * crash or wrong answer within a few iterations. */
    for (int i = 0; i < 100; i++) {
        assert_eval_eq("Expand[(x + y + z)^4]",
                       "x^4 + 4 x^3 y + 6 x^2 y^2 + 4 x y^3 + y^4 + "
                       "4 x^3 z + 12 x^2 y z + 12 x y^2 z + 4 y^3 z + "
                       "6 x^2 z^2 + 12 x y z^2 + 6 y^2 z^2 + 4 x z^3 + "
                       "4 y z^3 + z^4", 0);
    }
}

static void test_repeated_polynomial_factoring(void) {
    for (int i = 0; i < 50; i++) {
        assert_eval_eq("Factor[x^4 - 1]", "(-1 + x) (1 + x) (1 + x^2)", 0);
        assert_eval_eq("PolynomialGCD[x^3 - 1, x^2 - 1]", "-1 + x", 0);
    }
}

static void test_repeated_simplify(void) {
    for (int i = 0; i < 50; i++) {
        assert_eval_eq("Simplify[(x^2 - 1)/(x - 1)]", "1 + x", 0);
    }
}

static void test_repeated_derivative(void) {
    for (int i = 0; i < 50; i++) {
        assert_eval_eq("D[Sin[x^2], x]", "2 x Cos[x^2]", 0);
        assert_eval_eq("D[x^3 + 2 x^2 + x + 1, x]", "1 + 4 x + 3 x^2", 0);
    }
}

/* ------------------------------------------------------------------ */
/* 6. Random refcount-stress walk                                      */
/* ------------------------------------------------------------------ */

/* Build a population of random FUNCTION trees, then on each step:
 * randomly pick one and either:
 *   - clone it via expr_copy() (share), or
 *   - unshare a clone, mutate args[0], and verify original unchanged, or
 *   - drop a ref via expr_free.
 * Repeat for many iterations. Goal: catch any path where sharing
 * corrupts data or refcounts go astray.
 *
 * We use deterministic pseudo-randomness so failures are reproducible.
 */
static unsigned int rstate = 0xC0FFEE;
static unsigned int rrand(void) {
    /* xorshift32 */
    unsigned int x = rstate;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    rstate = x;
    return x;
}

static Expr* build_random_tree(int depth) {
    if (depth == 0 || (rrand() & 0x7) == 0) {
        return expr_new_integer((int64_t)(rrand() % 1000));
    }
    int n = 2 + (int)(rrand() % 3);  /* arity 2..4 */
    Expr** kids = malloc(sizeof(Expr*) * n);
    for (int i = 0; i < n; i++) kids[i] = build_random_tree(depth - 1);
    static const char* heads[] = { "Plus", "Times", "List", "f" };
    Expr* fn = expr_new_function(expr_new_symbol(heads[rrand() % 4]), kids, n);
    free(kids);
    return fn;
}

static void test_random_share_unshare_walk(void) {
    enum { POP = 32, ITERS = 4000 };
    Expr* pop[POP];
    for (int i = 0; i < POP; i++) pop[i] = build_random_tree(3);

    for (int it = 0; it < ITERS; it++) {
        unsigned int op = rrand() & 0x3;
        unsigned int idx = rrand() % POP;

        switch (op) {
            case 0: {
                /* Share: clone + immediately drop. Net no-op for refcount. */
                Expr* c = expr_copy(pop[idx]);
                ASSERT(c == pop[idx]);
                ASSERT(pop[idx]->refcount >= 2);
                expr_free(c);
                break;
            }
            case 1: {
                /* Share twice, then unshare one clone, mutate it, verify
                 * original tree string is unchanged. */
                Expr* shared1 = expr_copy(pop[idx]);
                Expr* shared2 = expr_copy(pop[idx]);
                ASSERT(pop[idx]->refcount >= 3);

                char* before = expr_to_string(pop[idx]);

                /* Unshare shared2 (consumes its ref). */
                Expr* priv = expr_unshare(shared2);
                if (priv == pop[idx]) {
                    /* refcount-1 fast path can't happen here because we
                     * just bumped it. Only triggered if shared2 was
                     * already private. With pop[idx] having refcount >= 3
                     * before this branch, expr_unshare must allocate. */
                    ASSERT(0 && "unexpected fast path");
                }
                ASSERT(priv->refcount == 1);

                /* If it's a FUNCTION, replace args[0] with something
                 * recognizable. Original must NOT see the change. */
                if (priv->type == EXPR_FUNCTION && priv->data.function.arg_count > 0) {
                    expr_free(priv->data.function.args[0]);
                    priv->data.function.args[0] = expr_new_integer(-77777);
                }

                char* after = expr_to_string(pop[idx]);
                ASSERT_STR_EQ(before, after);
                free(before);
                free(after);

                expr_free(priv);
                expr_free(shared1);
                break;
            }
            case 2: {
                /* Replace pop[idx] with a fresh random tree.
                 * Drops the old ref. */
                expr_free(pop[idx]);
                pop[idx] = build_random_tree(2 + (int)(rrand() % 3));
                break;
            }
            case 3: {
                /* Round-trip through evaluate() (which will call
                 * expr_copy under the hood). The returned tree may
                 * differ structurally; we only check it parses. */
                Expr* c = expr_copy(pop[idx]);
                Expr* ev = evaluate(c);
                ASSERT(ev != NULL);
                expr_free(ev);
                break;
            }
        }
    }

    for (int i = 0; i < POP; i++) {
        ASSERT(pop[i]->refcount == 1);
        expr_free(pop[i]);
    }
}

/* ------------------------------------------------------------------ */
/* 7. Refcount accounting on deep evaluation                           */
/* ------------------------------------------------------------------ */

/* After Phase-2, evaluating a complex expression and freeing the
 * result must leave the symbol table's still-reachable atoms with
 * sane refcounts. We can't probe internal refcounts directly here,
 * but we can verify that evaluating, freeing, and re-evaluating
 * gives identical results indefinitely (no slow drift, no crash). */
static void test_evaluate_free_roundtrip_stable(void) {
    const char* prog = "Apart[1/((x-1)(x-2)(x-3))]";
    Expr* parsed = parse_expression(prog);
    char* canon = expr_to_string(evaluate(expr_copy(parsed)));

    for (int i = 0; i < 100; i++) {
        Expr* p = expr_copy(parsed);
        Expr* r = evaluate(p);
        char* s = expr_to_string(r);
        ASSERT_STR_EQ(s, canon);
        free(s);
        expr_free(r);
    }
    expr_free(parsed);
    free(canon);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    /* 1. Refcount mechanics */
    TEST(test_atom_refcount_inc_dec);
    TEST(test_function_refcount_inc_dec);
    TEST(test_atom_share_pointer_identity);

    /* 2. expr_unshare */
    TEST(test_unshare_unique_is_noop);
    TEST(test_unshare_atom_creates_private);
    TEST(test_unshare_function_creates_private_args_array);
    TEST(test_unshare_bigint_payload_independent);
    TEST(test_unshare_string_payload_independent);

    /* 3. Hot-path mutators */
    TEST(test_plus_numeric_contagion_repeated);
    TEST(test_times_numeric_contagion_repeated);
    TEST(test_flat_flatten_repeated);
    TEST(test_orderless_sort_repeated);
    TEST(test_flatten_sequences_repeated);
    TEST(test_quotient_remainder_repeated);
    TEST(test_print_does_not_mutate_input);
    TEST(test_print_negative_rational_does_not_mutate);

    /* 4. Pattern matching */
    TEST(test_repeated_binding_substitution);
    TEST(test_pattern_binding_evaluates_shared_subtree);

    /* 5. End-to-end stress */
    TEST(test_repeated_complex_evaluation);
    TEST(test_repeated_polynomial_factoring);
    TEST(test_repeated_simplify);
    TEST(test_repeated_derivative);

    /* 6. Random share/unshare walk */
    TEST(test_random_share_unshare_walk);

    /* 7. Roundtrip stability */
    TEST(test_evaluate_free_roundtrip_stable);

    printf("\nAll expression-sharing stress tests passed!\n");
    return 0;
}
