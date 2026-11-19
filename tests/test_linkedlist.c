/* Unit tests for LinkedList (src/linkedlist.c; design in LINKED_LIST_DESIGN.md).
 *
 * Two layers:
 *   - the storage core, driven directly through the ll_* C API: the two-run
 *     deque, chunk growth, rebase under queue traffic, typed -> boxed
 *     conversion, copies and bulk ingest;
 *   - the language surface: the LinkedList[] constructor, the atom's identity
 *     semantics (Head, AtomQ, SameQ, ordering, patterns), printing, Normal,
 *     Length and the constructor's messages.
 */
#include "print.h"
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "parse.h"
#include "eval.h"
#include "linkedlist.h"
#include "sym_intern.h"
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ helpers */

static void expect_list(const LinkedList* ll, const char* expected) {
    Expr* l = ll_to_list(ll);
    char* s = expr_to_string(l);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: list contents\n  Expected: %s\n  Actual:   %s\n", expected, s);
        abort();
    }
    free(s);
    expr_free(l);
}

static void expect_true(bool c, const char* what) {
    if (!c) { fprintf(stderr, "FAIL: %s\n", what); abort(); }
}

static Expr* eval_str(const char* in) {
    Expr* p = parse_expression(in);
    Expr* r = evaluate(p);
    expr_free(p);
    return r;
}

/* --------------------------------------------------------- storage: deque */

void test_push_both_ends() {
    LinkedList* ll = ll_new();
    expect_true(ll_length(ll) == 0 && ll_kind(ll) == LL_EMPTY, "new list is empty");
    ll_push_back(ll, expr_new_integer(3));
    ll_push_back(ll, expr_new_integer(4));
    ll_push_front(ll, expr_new_integer(2));
    ll_push_front(ll, expr_new_integer(1));
    expect_true(ll_kind(ll) == LL_I64, "machine integers stay typed");
    expect_list(ll, "{1, 2, 3, 4}");
    Expr* e = ll_get(ll, 2);
    expect_true(e->type == EXPR_INTEGER && e->data.integer == 3, "ll_get index 2");
    expr_free(e);
    ll_release(ll);
}

void test_pop_both_ends() {
    LinkedList* ll = ll_new();
    for (int i = 1; i <= 5; i++) ll_push_back_i64(ll, i);
    ll_push_front_i64(ll, 0);
    Expr* a = ll_pop_back(ll);   /* 5 */
    Expr* b = ll_pop_front(ll);  /* 0, from the front run */
    Expr* c = ll_pop_front(ll);  /* 1, from the back run's inner end */
    expect_true(a->data.integer == 5 && b->data.integer == 0 && c->data.integer == 1,
                "pops return the right ends");
    expr_free(a); expr_free(b); expr_free(c);
    expect_list(ll, "{2, 3, 4}");
    /* Exhaust through the front run's inner end too. */
    LinkedList* f = ll_new();
    for (int i = 1; i <= 3; i++) ll_push_front_i64(f, i);       /* {3, 2, 1} */
    Expr* x = ll_pop_back(f);                                   /* 1 */
    expect_true(x->data.integer == 1, "pop_back from front-only list");
    expr_free(x);
    expect_list(f, "{3, 2}");
    while (ll_length(f)) expr_free(ll_pop_back(f));
    expect_true(ll_pop_back(f) == NULL && ll_pop_front(f) == NULL, "pop of empty is NULL");
    ll_release(f);
    ll_release(ll);
}

void test_large_growth_and_index() {
    LinkedList* ll = ll_new();
    const int64_t n = 100000;
    for (int64_t i = 0; i < n; i++) ll_push_back_i64(ll, i);
    for (int64_t i = 1; i <= n; i++) ll_push_front_i64(ll, -i);
    expect_true(ll_length(ll) == 2 * n, "length after 2n pushes");
    /* Every index maps to the right element: logical i holds i - n. */
    for (int64_t i = 0; i < 2 * n; i += 997) {
        int64_t v;
        expect_true(ll_get_i64(ll, i, &v) && v == i - n, "random access across both runs");
    }
    ll_release(ll);
}

void test_queue_rebase_bounds_memory() {
    /* Steady-state queue of 10 elements through 10^6 operations: the dead
     * prefix must be reclaimed (rebase), so memory stays O(live). */
    LinkedList* ll = ll_new();
    for (int i = 0; i < 10; i++) ll_push_back_i64(ll, i);
    int64_t expect_front = 0;
    for (int64_t t = 0; t < 1000000; t++) {
        ll_push_back_i64(ll, 10 + t);
        Expr* e = ll_pop_front(ll);
        expect_true(e->data.integer == expect_front, "FIFO order through rebases");
        expect_front++;
        expr_free(e);
    }
    expect_true(ll_length(ll) == 10, "queue length steady");
    expect_true(ll_payload_bytes(ll) < 4096, "queue memory stays O(live)");
    ll_release(ll);
}

void test_typed_to_boxed_conversion() {
    LinkedList* ll = ll_new();
    ll_push_back_real(ll, 1.5);
    ll_push_front_real(ll, 0.5);
    expect_true(ll_kind(ll) == LL_REAL, "reals stay typed");
    ll_push_back(ll, expr_new_integer(7));        /* Integer into Real: boxed */
    expect_true(ll_kind(ll) == LL_BOXED, "mixed Integer/Real boxes");
    ll_push_back(ll, expr_new_symbol(intern_symbol("x")));
    expect_list(ll, "{0.5, 1.5, 7, x}");
    ll_release(ll);
    /* An emptied list retypes for free. */
    LinkedList* r = ll_new();
    ll_push_back(r, expr_new_symbol(intern_symbol("y")));
    expr_free(ll_pop_back(r));
    ll_push_back_real(r, 2.0);
    expect_true(ll_kind(r) == LL_REAL, "empty boxed list retypes");
    ll_release(r);
}

void test_complex_typed() {
    LinkedList* ll = ll_new();
    ll_push_back_cplx(ll, 1.0, 2.0);
    ll_push_back(ll, eval_str("Complex[3., -1.]"));
    expect_true(ll_kind(ll) == LL_CPLX, "machine complex stays typed");
    expect_list(ll, "{1.0 + 2.0*I, 3.0 - 1.0*I}");
    double re, im;
    expect_true(ll_get_cplx(ll, 1, &re, &im) && re == 3.0 && im == -1.0, "typed complex read");
    ll_release(ll);
}

void test_set_and_copy() {
    LinkedList* ll = ll_new();
    for (int i = 1; i <= 4; i++) ll_push_back_i64(ll, i);
    ll_push_front_i64(ll, 0);
    LinkedList* c = ll_copy(ll);
    ll_set(ll, 0, expr_new_integer(100));
    ll_set(ll, 4, expr_new_symbol(intern_symbol("z")));   /* converts ll to boxed */
    expect_list(ll, "{100, 1, 2, 3, z}");
    expect_list(c, "{0, 1, 2, 3, 4}");                     /* copy unaffected */
    expect_true(ll_kind(c) == LL_I64, "copy keeps typed storage");
    LinkedList* c2 = ll_copy(ll);                          /* boxed copy */
    ll_release(ll);                                        /* elements shared by refcount */
    expect_list(c2, "{100, 1, 2, 3, z}");
    ll_release(c2);
    ll_release(c);
}

void test_bulk_ingest() {
    LinkedList* ll = ll_new();
    Expr* packed = eval_str("Developer`ToPackedArray[Range[1000] / 2.]");
    expect_true(ll_append_elements(ll, packed), "packed real ingest");
    expect_true(ll_kind(ll) == LL_REAL && ll_length(ll) == 1000, "packed ingest stays typed");
    double v;
    expect_true(ll_get_real(ll, 999, &v) && v == 500.0, "last packed element");
    expr_free(packed);
    Expr* plain = eval_str("{a, {b, c}}");
    expect_true(ll_append_elements(ll, plain), "boxed ingest");
    expect_true(ll_length(ll) == 1002 && ll_kind(ll) == LL_BOXED, "mixed ingest boxes");
    expr_free(plain);
    Expr* notlist = eval_str("f[1]");
    expect_true(!ll_append_elements(ll, notlist), "non-list source declines");
    expr_free(notlist);
    ll_release(ll);
}

/* ---------------------------------------------------------- the surface */

void test_constructor_and_printing() {
    assert_eval_eq("LinkedList[]", "LinkedList[<0>]", 0);
    assert_eval_eq("LinkedList[{1, 2, 3}]", "LinkedList[<3>]", 0);
    assert_eval_eq("LinkedList[{1, 2, 3}]", "LinkedList[<3>]", 1);           /* FullForm */
    assert_eval_eq("InputForm[LinkedList[{a}]]", "LinkedList[<1>]", 0);
    assert_eval_eq("LinkedList[Range[10^5]]", "LinkedList[<100000>]", 0);
    assert_eval_eq("LinkedList[NDArray[{1., 2.}]]", "LinkedList[<2>]", 0);
    assert_eval_eq("Normal[LinkedList[{1, x, 2.5}]]", "{1, x, 2.5}", 0);
    assert_eval_eq("Normal[LinkedList[]]", "{}", 0);
    assert_eval_eq("Length[LinkedList[{a, b, c}]]", "3", 0);
    assert_eval_eq("Normal[LinkedList[Range[5]]] === Range[5]", "True", 0);
}

void test_identity_semantics() {
    assert_eval_eq("Head[LinkedList[]]", "LinkedList", 0);
    assert_eval_eq("AtomQ[LinkedList[{1}]]", "True", 0);
    assert_eval_eq("MatchQ[LinkedList[], _LinkedList]", "True", 0);
    assert_eval_eq("ll1 = LinkedList[{1}]; ll2 = ll1; ll1 === ll2", "True", 0);
    /* Two lists with equal contents are different objects. */
    assert_eval_eq("LinkedList[{1}] === LinkedList[{1}]", "False", 0);
    /* Each LinkedList[] call creates a new object (no memoised constructor). */
    assert_eval_eq("a1 = LinkedList[]; a2 = LinkedList[]; a1 === a2", "False", 0);
    /* The copy constructor gives an independent object with the same contents. */
    assert_eval_eq("c0 = LinkedList[{1, 2}]; c1 = LinkedList[c0]; {c0 === c1, Normal[c1]}",
                   "{False, {1, 2}}", 0);
    /* Ordering is stable (creation order) and Union keeps distinct objects. */
    assert_eval_eq("u = LinkedList[]; Length[Union[{u, u, LinkedList[]}]]", "2", 0);
    assert_eval_eq("Attributes[LinkedList]", "{Protected}", 0);
    assert_eval_eq("NumericQ[LinkedList[{1}]]", "False", 0);
    assert_eval_eq("N[LinkedList[{1}]]", "LinkedList[<1>]", 0);
}

void test_constructor_messages() {
    char path[] = "/tmp/mth_ll_msgXXXXXX";
    int fd = mkstemp(path);
    if (fd >= 0) close(fd);
    int saved = test_stderr_redirect(path);
    assert_eval_eq("LinkedList[5]", "LinkedList[5]", 0);
    assert_eval_eq("LinkedList[{1}, {2}]", "LinkedList[{1}, {2}]", 0);
    test_stderr_restore(saved);
    FILE* f = fopen(path, "r");
    char buf[1024] = {0};
    size_t got = f ? fread(buf, 1, sizeof buf - 1, f) : 0;
    if (f) fclose(f);
    unlink(path);
    buf[got] = 0;
    expect_true(strstr(buf, "LinkedList::list") != NULL, "LinkedList::list emitted");
    expect_true(strstr(buf, "LinkedList::argt") != NULL, "LinkedList::argt emitted");
    /* Routed through the funnel: Quiet silences, Check sees it. */
    assert_eval_eq("Check[LinkedList[5], failed]", "failed", 0);
    assert_eval_eq("Quiet[Check[LinkedList[5], failed]]", "failed", 0);
}

void test_bytecount() {
    assert_eval_eq("ByteCount[LinkedList[Range[1000]]] > 8000", "True", 0);
    assert_eval_eq("ByteCount[LinkedList[]] < 2000", "True", 0);
}

/* ------------------------------------------------- Push / Pop / PopFirst */

void test_stack_and_queue() {
    assert_eval_eq("st = LinkedList[]; Push[st, a]; Push[st, b]; Push[st, c]; "
                   "{Pop[st], Pop[st], Normal[st]}", "{c, b, {a}}", 0);
    assert_eval_eq("qu = LinkedList[{1, 2}]; Push[qu, 3]; "
                   "{PopFirst[qu], PopFirst[qu], PopFirst[qu], Length[qu]}", "{1, 2, 3, 0}", 0);
    /* Push returns the list itself (the same object). */
    assert_eval_eq("pr = LinkedList[]; Push[pr, 1] === pr", "True", 0);
    /* Aliases see the mutation: reference semantics. */
    assert_eval_eq("al = LinkedList[]; al2 = al; Push[al2, x]; Normal[al]", "{x}", 0);
    /* A helper mutates its caller's list. */
    assert_eval_eq("enq[q_, v_] := Push[q, v]; hq = LinkedList[]; enq[hq, 7]; enq[hq, 8]; "
                   "Normal[hq]", "{7, 8}", 0);
    /* The copy constructor is independent of later pushes. */
    assert_eval_eq("o1 = LinkedList[{1}]; o2 = LinkedList[o1]; Push[o1, 2]; "
                   "{Normal[o1], Normal[o2]}", "{{1, 2}, {1}}", 0);
    /* Plain Lists: rebind like AppendTo. */
    assert_eval_eq("pl = {1, 2, 3}; {Push[pl, 4], Pop[pl], PopFirst[pl], pl}",
                   "{{1, 2, 3, 4}, 4, 1, {2, 3}}", 0);
    /* Typed storage survives pushes of machine reals. */
    assert_eval_eq("tr = LinkedList[]; Do[Push[tr, N[i]], {i, 1000}]; "
                   "Developer`PackedArrayQ[Normal[tr]]", "True", 0);
    assert_eval_eq("Attributes[{Push, Pop, PopFirst}]",
                   "{{HoldFirst, Protected}, {HoldFirst, Protected}, {HoldFirst, Protected}}", 0);
}

void test_mutation_invalidates_caches() {
    /* Length must be recomputed after every push (no memoised fixed point). */
    assert_eval_eq("mc = LinkedList[]; Table[Push[mc, i]; Length[mc], {3}]", "{1, 2, 3}", 0);
    assert_eval_eq("mc2 = LinkedList[]; l1 = Length[mc2]; Push[mc2, 1]; {l1, Length[mc2]}",
                   "{0, 1}", 0);
    /* A rule whose condition reads the list must re-fire after a mutation,
     * including through a stored unevaluated call (r1 holds big[<list>]). */
    assert_eval_eq("big[q_] /; Length[q] > 1 := yes; mq = LinkedList[{1}]; r1 = big[mq]; "
                   "Push[mq, 2]; {r1, big[mq]}", "{yes, yes}", 0);
}

void test_stack_messages() {
    char path[] = "/tmp/mth_ll_msgXXXXXX";
    int fd = mkstemp(path);
    if (fd >= 0) close(fd);
    int saved = test_stderr_redirect(path);
    assert_eval_eq("em = LinkedList[]; Pop[em]", "Pop[em]", 0);
    assert_eval_eq("PopFirst[em]", "PopFirst[em]", 0);
    assert_eval_eq("Pop[neverAssignedSym]", "Pop[neverAssignedSym]", 0);
    assert_eval_eq("ix = 5; Push[ix, 1]", "Push[ix, 1]", 0);
    assert_eval_eq("cy = LinkedList[]; Push[cy, cy]", "Push[cy, LinkedList[<0>]]", 0);
    assert_eval_eq("ca = LinkedList[]; cb = LinkedList[{ca}]; Push[ca, {f[cb]}]",
                   "Push[ca, {f[LinkedList[<1>]]}]", 0);
    assert_eval_eq("Push[em]", "Push[em]", 0);
    test_stderr_restore(saved);
    FILE* f = fopen(path, "r");
    char buf[4096] = {0};
    size_t got = f ? fread(buf, 1, sizeof buf - 1, f) : 0;
    if (f) fclose(f);
    unlink(path);
    buf[got] = 0;
    expect_true(strstr(buf, "Pop::empty") != NULL, "Pop::empty");
    expect_true(strstr(buf, "PopFirst::empty") != NULL, "PopFirst::empty");
    expect_true(strstr(buf, "Pop::rvalue") != NULL, "Pop::rvalue");
    expect_true(strstr(buf, "Push::normal") != NULL, "Push::normal");
    expect_true(strstr(buf, "Push::llcyc") != NULL, "Push::llcyc");
    expect_true(strstr(buf, "Push::argrx") != NULL || strstr(buf, "Push::argr") != NULL,
                "Push arity");
    /* The cycle-rejected list is unchanged, and the messages reach Check. */
    assert_eval_eq("Length[cy]", "0", 0);
    assert_eval_eq("Quiet[Check[Pop[em], failed]]", "failed", 0);
    assert_eval_eq("Quiet[Check[Push[cy, cy], failed]]", "failed", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_push_both_ends);
    TEST(test_pop_both_ends);
    TEST(test_large_growth_and_index);
    TEST(test_queue_rebase_bounds_memory);
    TEST(test_typed_to_boxed_conversion);
    TEST(test_complex_typed);
    TEST(test_set_and_copy);
    TEST(test_bulk_ingest);
    TEST(test_constructor_and_printing);
    TEST(test_identity_semantics);
    TEST(test_constructor_messages);
    TEST(test_bytecount);
    TEST(test_stack_and_queue);
    TEST(test_mutation_invalidates_caches);
    TEST(test_stack_messages);

    printf("All LinkedList tests passed!\n");
    symtab_clear();
    return 0;
}
