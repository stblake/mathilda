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

/* ------------------------------------------------- stress regressions */

void test_deep_nesting_no_recursion() {
    /* A 10^5-deep chain of lists: freeing it, the cycle check and ByteCount
     * all overflowed the C stack when they recursed per level. */
    assert_eval_eq("dc = LinkedList[]; dcur = dc; "
                   "Do[dn = LinkedList[]; Push[dcur, dn]; dcur = dn, {100000}]; "
                   "{Quiet[Check[Push[dcur, dc], rejected]], ByteCount[dc] > 10^6}",
                   "{rejected, True}", 0);
    assert_eval_eq("Clear[dc, dcur, dn]; done", "done", 0);
}

void test_abort_and_memory_limits_leave_list_consistent() {
    assert_eval_eq("ta = LinkedList[]; "
                   "{TimeConstrained[Do[Push[ta, i], {i, 10^9}], 0.1, aborted], "
                   "Normal[ta] === Range[Length[ta]], Length[ta] > 0}",
                   "{aborted, True, True}", 0);
    assert_eval_eq("tm = LinkedList[]; "
                   "{MemoryConstrained[Do[Push[tm, N[i]], {i, 10^8}], 10^7, overmem], "
                   "Normal[tm] === N[Range[Length[tm]]], Length[tm] > 0}",
                   "{overmem, True, True}", 0);
}

void test_identity_keys_and_sharing() {
    assert_eval_eq("ik1 = LinkedList[]; ik2 = LinkedList[]; ias = <|ik1 -> 1, ik2 -> 2|>; "
                   "Do[Push[ik1, i], {i, 100}]; {ias[ik1], ias[ik2]}", "{1, 2}", 0);
    assert_eval_eq("sin = LinkedList[{1}]; sout = LinkedList[{sin, {sin}}]; Push[sin, 2]; "
                   "{Normal[Normal[sout][[1]]], Normal[Normal[sout][[2, 1]]]}",
                   "{{1, 2}, {1, 2}}", 0);
}


/* ------------------------------------------- phase 2: reads, Part, cursor */

/* The span cursor visits every element exactly once, in logical order, across
 * both runs and every chunk boundary. */
void test_span_cursor_covers_both_runs() {
    LinkedList* ll = ll_new();
    for (int64_t i = 0; i < 1000; i++) ll_push_back_i64(ll, i);
    for (int64_t i = 1; i <= 777; i++) ll_push_front_i64(ll, -i);
    int64_t n = ll_length(ll), i = 0, calls = 0;
    while (i < n) {
        const void* p;
        int dir;
        int64_t k = ll_span_at(ll, i, &p, &dir);
        expect_true(k >= 1 && i + k <= n, "span length in range");
        const int64_t* d = (const int64_t*)p;
        for (int64_t t = 0; t < k; t++) {
            int64_t want;
            ll_get_i64(ll, i + t, &want);
            expect_true(d[dir > 0 ? t : -t] == want, "span slot = logical element");
        }
        i += k;
        calls++;
    }
    expect_true(calls < 40, "O(log n) cursor calls");
    ll_release(ll);
}

/* A machine float Total over a list grown at BOTH ends -- leaves that
 * straddle chunks and the downward-growing front run -- bit-matches the
 * packed Total of Normal[ll], serial and threaded sizes alike. */
static void check_sum_parity(LinkedList* ll, const char* what) {
    Expr* fast = ll_sum_machine(ll);
    Expr* lst = ll_to_list(ll);
    Expr* a[1] = { lst };
    Expr* call = expr_new_function(expr_new_symbol(intern_symbol("Total")), a, 1);
    Expr* ref = evaluate(call);
    expr_free(call);
    if (!fast || !expr_eq(fast, ref)) {
        char* x = fast ? expr_to_string(fast) : NULL;
        char* y = expr_to_string(ref);
        fprintf(stderr, "FAIL: %s: ll %s, packed %s\n", what, x ? x : "NULL", y);
        abort();
    }
    expr_free(fast);
    expr_free(ref);
}

void test_total_bit_parity() {
    int64_t sizes[] = { 1, 7, 128, 129, 1000, 70001, 300000 };
    uint64_t g = 12345;
    for (size_t s = 0; s < sizeof sizes / sizeof sizes[0]; s++) {
        LinkedList* r = ll_new();
        LinkedList* c = ll_new();
        for (int64_t i = 0; i < sizes[s]; i++) {
            g = g * UINT64_C(6364136223846793005) + 1;
            double v = (double)(g >> 11) / 9007199254740992.0 - 0.3;
            if (g & 1) ll_push_front_real(r, v); else ll_push_back_real(r, v);
            if (g & 2) ll_push_front_cplx(c, v, 1.0 - v); else ll_push_back_cplx(c, v * 3, v);
        }
        check_sum_parity(r, "real Total parity");
        check_sum_parity(c, "complex Total parity");
        ll_release(r);
        ll_release(c);
    }
    /* int64 overflow promotes to an exact bignum. */
    LinkedList* b = ll_new();
    for (int i = 0; i < 9; i++) ll_push_front_i64(b, INT64_MAX - i);
    ll_push_back_i64(b, -3);
    check_sum_parity(b, "int64 overflow -> bignum");
    ll_release(b);
}

void test_part_reads() {
    assert_eval_eq("ll = LinkedList[Range[10^5]]; {ll[[1234]], ll[[-1]], ll[[0]]}",
                   "{1234, 100000, LinkedList}", 0);
    assert_eval_eq("Head[ll[[{1, -1, 7}]]]", "LinkedList", 0);
    assert_eval_eq("Normal[ll[[{1, -1, 7}]]]", "{1, 100000, 7}", 0);
    assert_eval_eq("Normal[ll[[3 ;; 12 ;; 3]]]", "{3, 6, 9, 12}", 0);
    assert_eval_eq("Normal[ll[[-1 ;; -5 ;; -2]]]", "{100000, 99998, 99996}", 0);
    assert_eval_eq("Normal[ll[[5 ;; 4]]]", "{}", 0);
    assert_eval_eq("Length[ll[[All]]]", "100000", 0);
    assert_eval_eq("ll[[All]] === ll", "False", 0);                 /* a new list */
    /* Deeper subscripts apply to the selected elements. */
    assert_eval_eq("m = LinkedList[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]; "
                   "{m[[2, 3]], Normal[m[[All, 2]]], Normal[m[[{1, 3}, 2 ;; 3]]], m[[2]]}",
                   "{6, {2, 5, 8}, {{2, 3}, {8, 9}}, {4, 5, 6}}", 0);
    /* Packed rows and nested LinkedLists as elements. */
    assert_eval_eq("p = LinkedList[]; Push[p, Range[400]]; Push[p, Range[401, 800]]; "
                   "{p[[2, -1]], Normal[p[[All, 1]]]}", "{800, {1, 401}}", 0);
    assert_eval_eq("n = LinkedList[{LinkedList[{a, b}], LinkedList[{c, d}]}]; "
                   "{n[[2, 1]], {n}[[1, 1, 2]], Normal[n[[All, -1]]]}", "{c, b, {b, d}}", 0);
    /* Typed selections stay typed (and exact heads survive). */
    assert_eval_eq("r = LinkedList[N[Range[10]]]; {Normal[r[[2 ;; 4]]], Head[r[[2]]]}",
                   "{{2.0, 3.0, 4.0}, Real}", 0);
    assert_eval_eq("Extract[m, {2, 1}]", "4", 0);
    /* Reads see mutations (no stale memo). */
    assert_eval_eq("q = LinkedList[{1}]; Table[Push[q, i]; {q[[-1]], Last[q]}, {i, 3}]",
                   "{{1, 1}, {2, 2}, {3, 3}}", 0);
}

void test_part_messages() {
    assert_eval_eq("Check[LinkedList[{1, 2}][[5]], err]", "err", 0);
    assert_eval_eq("Check[LinkedList[{1, 2}][[{1, -3}]], err]", "err", 0);
    assert_eval_eq("Check[LinkedList[{1, 2}][[1 ;; 3]], err]", "err", 0);
    assert_eval_eq("Check[LinkedList[{1, 2, 3}][[3 ;; 1]], err]", "err", 0);
    assert_eval_eq("Quiet[LinkedList[{1, 2}][[5]]]", "Part[LinkedList[<2>], 5]", 0);
    assert_eval_eq("Check[First[LinkedList[]], err]", "err", 0);
    assert_eval_eq("Check[Last[LinkedList[]], err]", "err", 0);
    assert_eval_eq("{First[LinkedList[], d], Last[LinkedList[], d]}", "{d, d}", 0);
}

void test_reading_heads() {
    assert_eval_eq("ll = LinkedList[Range[10^6]]; {First[ll], Last[ll], Total[ll]}",
                   "{1, 1000000, 500000500000}", 0);
    assert_eval_eq("{MemberQ[ll, 500], MemberQ[ll, -1], MemberQ[ll, 5.], Count[ll, 7], "
                   "Count[ll, _Integer], Count[ll, _Real], Count[ll, _?EvenQ], "
                   "MemberQ[ll, _LinkedList, {0}], MemberQ[ll, _LinkedList]}",
                   "{True, False, False, 1, 1000000, 0, 500000, True, False}", 0);
    assert_eval_eq("Total[LinkedList[{2^62, 2^62, 2^62, -5}]]", "13835058055282163707", 0);
    assert_eval_eq("m = LinkedList[{{1, 2, 3}, {4, 5, 6}}]; "
                   "{Total[m], Total[m, 2], Normal[Total[m, {2}]], Total[m, Infinity], Total[m, -1]}",
                   "{{5, 7, 9}, 21, {6, 15}, 21, 21}", 0);
    assert_eval_eq("s = LinkedList[{x, y, x, f[x]}]; "
                   "{Total[s], Count[s, x], Count[s, x, Infinity], MemberQ[s, f[_]], Count[s, _, {2}]}",
                   "{2 x + y + f[x], 2, 3, True, 1}", 0);
    assert_eval_eq("{Total[LinkedList[]], Total[LinkedList[{1.*^308, 1.*^308}]]}",
                   "{0, 2.0*^308}", 0);
    /* Nested LinkedLists are walked like nested Lists, at every level spec. */
    assert_eval_eq("n = LinkedList[{LinkedList[{1, 2}], LinkedList[{3, x}]}]; "
                   "{Count[n, _Integer, 2], Count[n, _Integer, {-1}], MemberQ[n, x, Infinity], "
                   "MemberQ[{n}, 3, 3], Count[n, _LinkedList], Count[{1, LinkedList[{1, 1}]}, 1, Infinity], "
                   "Normal[Total[n, {2}]]}",
                   "{3, 3, True, True, 2, 3, {3, 3 + x}}", 0);
    /* The same answers as the List, element for element. */
    assert_eval_eq("v = {3, x, 2.5, {1, 2}, f[3], 3}; w = LinkedList[v]; "
                   "Table[{Count[w, pt, lv], MemberQ[w, pt, lv]} === {Count[v, pt, lv], MemberQ[v, pt, lv]}, "
                   "{pt, {3, _Integer, _Real, _List, _, f[_]}}, {lv, {1, 2, {2}, Infinity, {-1}}}] // Flatten // Union",
                   "{True}", 0);
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
    TEST(test_deep_nesting_no_recursion);
    TEST(test_abort_and_memory_limits_leave_list_consistent);
    TEST(test_identity_keys_and_sharing);
    TEST(test_span_cursor_covers_both_runs);
    TEST(test_total_bit_parity);
    TEST(test_part_reads);
    TEST(test_part_messages);
    TEST(test_reading_heads);

    printf("All LinkedList tests passed!\n");
    symtab_clear();
    return 0;
}
