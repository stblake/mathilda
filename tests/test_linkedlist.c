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
#include "compile/autocompile.h"
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


/* ------------------------------------------- phase 3: in-place mutators */

void test_append_prepend_to() {
    assert_eval_eq("a = LinkedList[]; AppendTo[a, 1]; PrependTo[a, 0]; AppendTo[a, 2]; "
                   "{Normal[a], AppendTo[a, 3] === a, PrependTo[a, -1] === a}",
                   "{{0, 1, 2}, True, True}", 0);
    /* Through an alias and through a function parameter: the same object. */
    assert_eval_eq("a2 = a; AppendTo[a2, 4]; add[l_, x_] := AppendTo[l, x]; add[a, 5]; Normal[a]",
                   "{-1, 0, 1, 2, 3, 4, 5}", 0);
    assert_eval_eq("Quiet[Check[AppendTo[a, a], failed]]", "failed", 0);
    assert_eval_eq("Quiet[Check[PrependTo[a, {f[a]}], failed]]", "failed", 0);
    assert_eval_eq("Length[a]", "7", 0);
    /* Typed storage survives in-place growth at both ends. */
    assert_eval_eq("b = LinkedList[]; Do[AppendTo[b, i]; PrependTo[b, -i], {i, 1000}]; "
                   "{Length[b], b[[1]], b[[-1]], Total[b], Normal[b] === Join[Range[-1000, -1], Range[1000]]}",
                   "{2000, -1000, 1000, 0, True}", 0);
    /* Reads after each mutation see it (no stale memo). */
    assert_eval_eq("t = LinkedList[{1, 2, 3}]; Table[AppendTo[t, i]; {Length[t], Last[t]}, {i, 3}]",
                   "{{4, 1}, {5, 2}, {6, 3}}", 0);
}

void test_part_assignment() {
    assert_eval_eq("ll = LinkedList[Range[10]]; ll2 = ll; ll[[3]] = 30; ll[[-1]] = 100; Normal[ll2]",
                   "{1, 2, 30, 4, 5, 6, 7, 8, 9, 100}", 0);
    assert_eval_eq("ll[[1 ;; 4]] = {pa, pb, pc, pd}; ll[[{5, 7}]] = 0; ll[[2 ;; 10 ;; 4]] = pz; Normal[ll]",
                   "{pa, pz, pc, pd, 0, pz, 0, 8, 9, pz}", 0);
    assert_eval_eq("(ll[[2]] = pw) === pw", "True", 0);
    /* Errors leave the list unchanged and the Set unevaluated. */
    assert_eval_eq("{Quiet[Check[ll[[11]] = 1, failed]], Quiet[Check[ll[[9 ;; 12]] = 1, failed]], "
                   "Quiet[Check[ll[[{1, 20}]] = 1, failed]], Normal[ll][[1]]}",
                   "{failed, failed, failed, pa}", 0);
    /* Typed lists: scalar broadcast and a packed buffer written slot to slot. */
    assert_eval_eq("r = LinkedList[N[Range[5]]]; r[[All]] = 0.5; r[[2 ;; 4]] = {1., 2., 3.}; Normal[r]",
                   "{0.5, 1.0, 2.0, 3.0, 0.5}", 0);
    assert_eval_eq("q = LinkedList[Range[400]]; q[[1 ;; 300]] = Range[300, 1, -1]; Normal[q[[{1, 300, 301}]]]",
                   "{300, 1, 301}", 0);
    assert_eval_eq("q[[1]] = 2.5; {q[[1]], q[[2]], Head[q[[2]]]}", "{2.5, 299, Integer}", 0);
    /* A LinkedList right-hand side is distributed like a List. */
    assert_eval_eq("q[[1 ;; 3]] = LinkedList[{px, py, pz}]; Normal[q[[1 ;; 4]]]", "{px, py, pz, 297}", 0);
    /* Deeper subscripts: copy-on-write for an immutable element ... */
    assert_eval_eq("pm = LinkedList[{{1, 2}, {3, 4}}]; row = pm[[2]]; pm[[2, 1]] = px; pm[[All, 2]] = {pp, ps}; "
                   "{Normal[pm], row}", "{{{1, pp}, {px, ps}}, {3, 4}}", 0);
    /* ... and in place for a LinkedList element, at any depth. */
    assert_eval_eq("pn = LinkedList[{LinkedList[{1, 2}]}]; inner = pn[[1]]; pn[[1, 2]] = 99; "
                   "outer = {pn}; outer[[1, 1, 1]] = 7; Normal[inner]", "{7, 99}", 0);
    /* Cycles are refused, directly or through a value that contains the list. */
    assert_eval_eq("pcy = LinkedList[{1}]; {Quiet[Check[pcy[[1]] = pcy, failed]], "
                   "Quiet[Check[pcy[[1]] = {pcy}, failed]], Quiet[Check[pn[[1, 1]] = pn, failed]], Normal[pcy]}",
                   "{failed, failed, failed, {1}}", 0);
    /* The object itself as target: a function parameter writes the caller's list. */
    assert_eval_eq("setFirst[l_, v_] := (l[[1]] = v); setFirst[pcy, \"s\"]; Normal[pcy]", "{\"s\"}", 0);
}

void test_mutators_linear_time() {
    /* 10^6 appends and prepends and 10^5 Part writes in well under a second
     * each: a regression to O(n) per call would take minutes. */
    assert_eval_eq("b = LinkedList[]; Do[AppendTo[b, i], {i, 10^6}]; "
                   "p = LinkedList[]; Do[PrependTo[p, N[i]], {i, 10^6}]; "
                   "Do[b[[i]] = -i, {i, 10^5}]; {Length[b], b[[5]], b[[-1]], p[[1]], Length[p]}",
                   "{1000000, -5, 1000000, 1.0*^6, 1000000}", 0);
}


/* ------------------------------------------------ phase 4: copying heads */

/* Every copying head on a LinkedList agrees with the same head on the List,
 * for typed (integer, real, complex) and boxed storage, and for a list built
 * from both ends (front-run spans run downward). */
void test_copying_heads_match_list() {
    assert_eval_eq("c4mk[v_] := Module[{l = LinkedList[], h = Floor[Length[v]/2]}, "
                   "Scan[PrependTo[l, #] &, Reverse[Take[v, h]]]; Scan[AppendTo[l, #] &, Drop[v, h]]; l]; "
                   "c4ops = {Take[#, 3] &, Take[#, -2] &, Take[#, {2, -2, 2}] &, Take[#, {-1, 1, -3}] &, "
                   "Take[#, UpTo[100]] &, Take[#, All] &, Drop[#, 2] &, Drop[#, -3] &, Drop[#, {2, 8, 3}] &, "
                   "Drop[#, {-2, 1, -2}] &, Drop[#, None] &, Most, Rest, Reverse, Append[#, c4q] &, "
                   "Prepend[#, 1] &, Insert[#, c4q, 2] &, Insert[#, c4q, -1] &, "
                   "Insert[#, c4q, {{1}, {4}, {4}, {13}}] &, Delete[#, 3] &, Delete[#, -1] &, "
                   "Delete[#, {{1}, {3}, {3}, {12}}] &, Delete[#, {}] &, ReplacePart[#, {2 -> c4q, -1 -> 5}] &, "
                   "ReplacePart[#, i_?OddQ :> i] &, ReplacePart[#, {{1}, {4}} -> 0] &, Select[#, NumericQ] &, "
                   "Select[#, NumericQ, 2] &, Map[c4f, #] &, Map[#^2 &, #] &, Map[c4f, #, {1}] &, "
                   "Join[#, #] &, Join[{c4q}, #, {2}] &, Sin, Exp, # + 1 &, # + 1.5 &, 2 # &, #^2 &, # + # &, # # &}; "
                   "c4vs = {Range[12], N[Range[12]]/3, Range[12] + 0.5 I, "
                   "{c4x, 1, 2.5, c4f[c4y], \"s\", 3, 4, 5, 6, 7, 8, 9}}; "
                   "Union[Flatten[Table[With[{r = op[mk[v]]}, Head[r] === LinkedList && Normal[r] === op[v]], "
                   "{op, c4ops}, {v, c4vs}, {mk, {LinkedList, c4mk}}]]]",
                   "{True}", 0);
    /* Deeper specifications act inside the elements, as for a List. */
    assert_eval_eq("c4m = LinkedList[{{1, 2}, {3, 4}, {5, 6}}]; c4n = Normal[c4m]; "
                   "Union[Table[Normal[op[c4m]] === op[c4n], {op, {Delete[#, {2, 1}] &, "
                   "Delete[#, {{2, 1}, {3}, {1, 2}}] &, Insert[#, c4x, {3, 2}] &, Reverse[#, {1, 2}] &, "
                   "Reverse[#, 2] &, Take[#, 2, 1] &, Take[#, All, -1] &, Drop[#, 1, -1] &, "
                   "ReplacePart[#, {2, 2} -> c4z] &, ReplacePart[#, {_, 1} -> 0] &, Map[c4f, #, {2}] &, "
                   "Map[c4f, #, 2] &, Map[c4f, #, {-1}] &, Join[#, #, 2] &, # + {10, 20, 30} &}}]]",
                   "{True}", 0);
    assert_eval_eq("{Map[c4f, c4m, {0, 1}] // Head, Map[c4f, c4m, {0}][[1]] === c4m, c4h[Delete[c4m, 0]], "
                   "c4h[Delete[c4m, {{}}]]}",
                   "{c4f, True, c4h[{1, 2}, {3, 4}, {5, 6}], c4h[]}", 0);
    /* Nested LinkedLists take the deeper specifications their own way. */
    assert_eval_eq("c4nn = LinkedList[{LinkedList[{1, 2, 3}], LinkedList[{4, 5}]}]; "
                   "{Head[Take[c4nn, All, 1][[1]]], Normal[Normal /@ Normal[Take[c4nn, All, -1]]], "
                   "Normal[Delete[c4nn, {1, 2}][[1]]], Normal[Normal[Reverse[c4nn, 2]][[2]]]}",
                   "{LinkedList, {{3}, {5}}, {1, 3}, {5, 4}}", 0);
}

void test_copying_heads_leave_source() {
    assert_eval_eq("c4s = LinkedList[Range[5]]; c4a = c4s; "
                   "{Append[c4s, 6], Prepend[c4s, 0], Take[c4s, 2], Drop[c4s, 2], Rest[c4s], Most[c4s], "
                   "Reverse[c4s], Join[c4s, c4s], Insert[c4s, 9, 2], Delete[c4s, 2], ReplacePart[c4s, 1 -> 0], "
                   "Map[Minus, c4s], Select[c4s, EvenQ], c4s + 1, Sin[c4s]}; "
                   "{Normal[c4a], Append[c4s, 6] === c4s, Join[c4s] === c4s}",
                   "{{1, 2, 3, 4, 5}, False, False}", 0);
    /* The result is a fresh, independent object. */
    assert_eval_eq("c4t = Take[c4s, 3]; AppendTo[c4t, 99]; {Normal[c4t], Length[c4s]}",
                   "{{1, 2, 3, 99}, 5}", 0);
}

void test_copying_messages() {
    assert_eval_eq("c4e = LinkedList[{1, 2, 3}]; "
                   "{Check[Take[c4e, 5], err], Check[Take[c4e, -4], err], Check[Take[c4e, {2, 7}], err], "
                   "Check[Drop[c4e, 4], err], Check[Take[c4e, c4x], err], Check[Most[LinkedList[]], err], "
                   "Check[Rest[LinkedList[]], err], Check[Insert[c4e, 0, 5], err], Check[Insert[c4e, 0, {5, 1}], err], "
                   "Check[Delete[c4e, 4], err], Check[Delete[c4e, {{1}, {-4}}], err], "
                   "Check[Join[c4e, c4f[1]], err], Check[c4e + {1, 2}, err]}",
                   "{err, err, err, err, err, err, err, err, err, err, err, err, err}", 0);
    assert_eval_eq("Quiet[{Take[c4e, 5], Delete[c4e, 4], Insert[c4e, 0, 0], Join[c4e, c4f[1]], c4e + {1, 2}}]",
                   "{Take[LinkedList[<3>], 5], Delete[LinkedList[<3>], 4], Insert[LinkedList[<3>], 0, 0], "
                   "Join[LinkedList[<3>], c4f[1]], LinkedList[<3>] + {1, 2}}", 0);
    /* In range at the boundaries: no message. */
    assert_eval_eq("Check[{Normal[Take[c4e, {4, 3}]], Normal[Drop[c4e, 3]], Normal[Insert[c4e, 0, 4]], "
                   "Normal[Insert[c4e, 0, -4]], Normal[Take[c4e, 0]]}, err]",
                   "{{}, {}, {1, 2, 3, 0}, {0, 1, 2, 3}, {}}", 0);
}

/* Listable threading over machine reals runs the kernel from chunk to chunk
 * and must agree with the packed List to the bit -- including where the
 * packed path declines to the List (overflow to MPFR, exact zeros). */
void test_listable_threading_bit_parity() {
    assert_eval_eq("SeedRandom[4]; c4r = Join[RandomReal[{-5, 5}, 70000], {0., -0., 800., -800.}]; "
                   "c4lr = LinkedList[c4r]; c4fr = LinkedList[]; Scan[PrependTo[c4fr, #] &, Reverse[c4r]]; "
                   "Union[Flatten[Table[Normal[op[l]] === op[c4r], "
                   "{op, {Sin, Cos, Exp, Tanh, Abs, Sqrt, Log, # + 1.5 &, 3 # &, # # &, # + # &, # 1.*^307 &, "
                   "#/3 &, # - 1 &, # + 0 &, 0 # &, #^2 &, Floor, Round}}, {l, {c4lr, c4fr}}]]]",
                   "{True}", 0);
    assert_eval_eq("c4i = LinkedList[{1, 2, 2^62, -3}]; "
                   "{Normal[c4i + 1], Normal[2 c4i], Normal[c4i^2] === {1, 4, 2^124, 9}, Normal[Mod[c4i, 3]]}",
                   "{{2, 3, 4611686018427387905, -2}, {2, 4, 9223372036854775808, -6}, True, {1, 2, 1, 0}}", 0);
    /* A user function with Listable threads too, and threading reads every
     * kind of operand: List, packed List, LinkedList. */
    assert_eval_eq("SetAttributes[c4g, Listable]; {Normal[c4g[LinkedList[{1, c4x}]]], "
                   "Normal[LinkedList[{1, 2}] + {c4x, 3} + Range[2] + LinkedList[{0.5, 0.5}]]}",
                   "{{c4g[1], c4g[c4x]}, {2.5 + c4x, 7.5}}", 0);
}

void test_copying_heads_scale() {
    /* 10^6 elements: Insert, Delete and ReplacePart copy spans (a List copies
     * Expr pointers one by one), so twenty of each stay far under a second. */
    assert_eval_eq("c4b = LinkedList[N[Range[10^6]]]; "
                   "c4t = First[AbsoluteTiming[Do[Insert[c4b, 0., 500000]; Delete[c4b, 500000]; "
                   "ReplacePart[c4b, 7 -> 0.]; Take[c4b, {10, -10}]; Reverse[c4b], {20}]]]; "
                   "{c4t < 1.5, Length[Insert[c4b, 0., 500000]], Total[Delete[c4b, -1]] == 499999500000}",
                   "{True, 1000001, True}", 0);
}


/* ------------------------------- everywhere: List heads on the chunks */

/* llchk[h, l]: h on a LinkedList of l agrees with h[Normal[ll]] -- held, so an
 * unevaluated h[ll] is never rescued by evaluating h[{...}] afterwards --
 * without materialising the list (only Normal may) and without changing it. */
static const char* LLCHK =
    "llchk[h_, l0_] := Module[{ll = LinkedList[l0], l, c0, r, c1}, l = Normal[ll]; "
    "c0 = Internal`LinkedListConversionCount[]; r = h[ll]; "
    "c1 = Internal`LinkedListConversionCount[]; "
    "{c1 == c0, (Hold[Evaluate[r]] //. x_LinkedList :> With[{v = Normal[x]}, v /; True]) === "
    "Hold[Evaluate[h[l]]], Normal[ll] === l}]; ";

static void expect_llchk(const char* heads, const char* fixtures) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "%sUnion[Flatten[Quiet[Table[llchk[h, l], {h, %s}, {l, %s}]]]]", LLCHK, heads, fixtures);
    assert_eval_eq(buf, "{True}", 0);
}

/* Map, Select and Listable threading read the slots: compiled functions and
 * predicates per slot, every ND kernel and the buffer arithmetic through the
 * engines' container-free cores, else f[e] per element -- never a List. */
void test_map_select_listable_direct() {
    const char* fx = "{N[Range[-6, 13]]/7, Range[-5, 14], N[Range[12]] - 0.5 I, {m1a, 2, m1b, 3.5}, {}}";
    expect_llchk("{Sin, Sqrt, Log, Floor, Abs, Round, # + 1 &, # 2 &, #^2 &, #^0.5 &, 2^# &, #^# &, "
                 "Mod[#, 3] &, ArcTan[#, #] &, BesselJ[1, #] &, LerchPhi[0.5, 2, #] &, # - 1 &, #/3 &, "
                 "-# &, # + # &, # # &, # + 2.5 &, IntegerLength, Map[#^2 + 1 &, #] &, Map[# + 0.5 &, #] &, "
                 "Map[m1f, #] &, Select[#, # > 0 &] &, Select[#, EvenQ] &, Select[#, EvenQ, 2] &}", fx);
    /* bignum promotion, exact Rationals and messages, as on the List */
    expect_llchk("{# + 2^62 &, #^-1 &, # 2^62 &, #^3 &}", "{{2^62, 1, 0, -2}}");
    /* the callback may change the list: the elements come from a snapshot */
    assert_eval_eq("m1x = LinkedList[{1, 2, 3}]; {Normal[Map[(AppendTo[m1x, 0]; #) &, m1x]], Length[m1x], "
                   "Normal[Select[m1x, (PrependTo[m1x, 9]; True) &]] === {1, 2, 3, 0, 0, 0}, Length[m1x]}",
                   "{{1, 2, 3}, 6, True, 12}", 0);
    /* a Throw or an Abort leaves the head */
    assert_eval_eq("{Catch[Map[If[# == 2, Throw[m1t], #] &, LinkedList[{1, 2, 3}]]], "
                   "Catch[Select[LinkedList[{1, 2, 3}], If[# == 2, Throw[m1u], True] &]], "
                   "CheckAbort[Map[If[# == 2, Abort[], #] &, LinkedList[{1, 2, 3}]], m1ab]}",
                   "{m1t, m1u, m1ab}", 0);
}

/* The level walkers and functional heads see a LinkedList's elements one level
 * down, as a List's: patterns, replacement, Apply, Fold, Scan, the threading
 * heads, MapIndexed, MapAt, Pick, Part/Extract with a LinkedList of indices. */
void test_walkers_and_functional() {
    const char* fx = "{Range[-3, 4], N[Range[5]]/3, {m2a, m2f[m2b], 2, {m2c, 2}}, {}}";
    expect_llchk("{Cases[#, 2] &, Cases[#, _Integer, Infinity] &, Position[#, 2] &, DeleteCases[#, 2] &, "
                 "DeleteCases[#, 2, Infinity] &, FreeQ[#, 2] &, FreeQ[#, m2c] &, Level[#, {1}] &, Level[#, {-1}] &, "
                 "Depth, LeafCount, FirstCase[#, _m2f] &, FirstPosition[#, 2] &, # /. 2 -> 5 &, "
                 "Replace[#, 2 -> 5, {1}] &, # //. 1 -> 2 &, m2g @@ # &, Apply[m2g, #, {1}] &, "
                 "Fold[m2h, 0, #] &, FoldList[m2h, 0, #] &, FoldList[Plus, #] &, Map[m2g, #, {2}] &, MapAll[m2g, #] &, "
                 "MapIndexed[m2g, #] &, MapAt[m2g, #, -1] &, MapThread[m2h, {#, #}] &, Outer[m2h, #, {1, 2}] &, "
                 "Inner[m2h, #, #, m2g] &, Thread[m2g[#, 1]] &, AllTrue[#, IntegerQ] &, AnyTrue[#, EvenQ] &, "
                 "NoneTrue[#, EvenQ] &, SelectFirst[#, EvenQ] &, LengthWhile[#, IntegerQ] &, TakeWhile[#, IntegerQ] &, "
                 "Pick[#, EvenQ /@ #] &, Part[{9, 8, 7, 6, 5, 4, 3, 2}, LinkedList[{1, 3}]] &, "
                 "Module[{m2s = {}}, Scan[(m2s = Append[m2s, #]) &, #]; m2s] &}", fx);
    /* ReplaceAll that matches nothing keeps the object; a change makes a new one */
    assert_eval_eq("m2x = LinkedList[{1, 2}]; {(m2x /. 7 -> 8) === m2x, (m2x /. 2 -> 8) === m2x}",
                   "{True, False}", 0);
    /* nested LinkedLists are walked at every depth */
    assert_eval_eq("m2n = LinkedList[{LinkedList[{1, 2}], LinkedList[{3, m2y}]}]; "
                   "{Position[m2n, m2y], Normal /@ Normal[m2n /. m2y -> 4], Depth[m2n], Level[m2n, {2}]}",
                   "{{{2, 2}}, {{1, 2}, {3, 4}}, 3, {1, 2, 3, m2y}}", 0);
}

/* Structure queries and the discrete iterators. */
void test_structure_and_iterators() {
    const char* fx = "{Range[-3, 4], N[Range[5]]/3, {m3a, 2}, {{1, 2}, {3, 4}}, {}}";
    expect_llchk("{VectorQ, VectorQ[#, NumberQ] &, MatrixQ, ArrayQ, Dimensions, ArrayDepth, Precision, Accuracy, "
                 "Table[m3i^2, {m3i, #}] &, Sum[m3i, {m3i, #}] &, Product[m3i, {m3i, #}] &, "
                 "Module[{m3s = 0}, Do[m3s += m3i, {m3i, #}]; m3s] &}", fx);
    /* a LinkedList of rows counts as a matrix; a LinkedList level of LinkedLists too */
    assert_eval_eq("{Dimensions[LinkedList[{LinkedList[{1, 2}], LinkedList[{3, 4}]}]], "
                   "MatrixQ[LinkedList[{LinkedList[{1, 2}], {3, 4}}]], VectorQ[LinkedList[{1, {2}}]]}",
                   "{{2, 2}, True, False}", 0);
    /* the iterator runs over the values the list had: the body may grow it */
    assert_eval_eq("m3y = LinkedList[{1, 2, 3}]; Do[AppendTo[m3y, m3i], {m3i, m3y}]; Normal[m3y]",
                   "{1, 2, 3, 1, 2, 3}", 0);
}

/* ---------------------------------------- phase 5: Compile and autocompile */

/* Module-local lists compile, at every element type, and answer exactly as the
 * interpreter does (the same bodies run through CompileDiagnostics and the
 * interpreter). */
void test_compile_module_local_lists() {
    static const char* bodies[] = {
        "Module[{b = LinkedList[]}, Do[AppendTo[b, i^2.], {i, n}]; Normal[b]]",
        "Module[{b = LinkedList[]}, Do[AppendTo[b, i^2], {i, n}]; Total[b]]",
        "Module[{b = LinkedList[]}, Do[PrependTo[b, i + 2. I], {i, n}]; b[[1]] + Last[b] + Length[b]]",
        "Module[{st = LinkedList[], acc = 0}, Do[Push[st, i]; If[Mod[i, 3] == 0, acc = acc + Pop[st]], {i, n}]; acc + 100 Length[st]]",
        "Module[{q = LinkedList[], s = 0.}, Do[Push[q, N[i]]; If[i > 3, s = s + PopFirst[q]], {i, n}]; s]",
        "Module[{b = LinkedList[Range[n]]}, Do[b[[i]] = b[[i]] b[[-i]], {i, n}]; Total[b]]",
        "Module[{b = LinkedList[N[Range[n]]]}, Normal[Join[Reverse[Rest[b]], Join[Take[b, 2], Drop[b, -3]]]]]",
        "Module[{b = LinkedList[N[Range[n]]]}, Normal[Join[Append[Most[b], 0.5], Prepend[b, -1.]]]]",
        "Module[{b = LinkedList[Range[n]], a}, a = 0; Module[{c = b}, Push[c, 100]]; Length[b] + Last[b]]",
    };
    char buf[1024];
    for (size_t k = 0; k < sizeof bodies / sizeof bodies[0]; k++) {
        snprintf(buf, sizeof buf, "CompileDiagnostics[{{n, _Integer}}, %s][[1, 2]]", bodies[k]);
        assert_eval_eq(buf, "True", 0);
        snprintf(buf, sizeof buf, "Compile[{{n, _Integer}}, %s][9] === Block[{n = 9}, %s]",
                 bodies[k], bodies[k]);
        assert_eval_eq(buf, "True", 0);
    }
    /* An Integer push into a guessed-Real local retypes it (above); pushes of
     * two machine types make a mixed list, which is not a machine list, so the
     * body runs through the interpreter -- with the interpreter's answer. */
    assert_eval_eq("c5mb = Module[{b = LinkedList[]}, Do[AppendTo[b, i]; AppendTo[b, N[i]], {i, n}]; Normal[b]]; "
                   "{CompileDiagnostics[{{n, _Integer}}, Module[{b = LinkedList[]}, "
                   "Do[AppendTo[b, i]; AppendTo[b, N[i]], {i, n}]; Normal[b]]][[1, 2]], "
                   "Compile[{{n, _Integer}}, Module[{b = LinkedList[]}, Do[AppendTo[b, i]; AppendTo[b, N[i]], {i, n}]; Normal[b]]][2]}",
                   "{False, {1, 1.0, 2, 2.0}}", 0);
}

/* A {b, _LinkedList} argument is the caller's own list: mutated in place. */
void test_compile_list_arguments() {
    assert_eval_eq("c5h = Compile[{{b, _LinkedList}, {x, _Real}}, AppendTo[b, x]; PrependTo[b, -x]; Length[b]]; "
                   "c5l = LinkedList[{1., 2.}]; {c5h[c5l, 7.], Normal[c5l], "
                   "CompileDiagnostics[{{b, _LinkedList}, {x, _Real}}, AppendTo[b, x]; Length[b]][[1, 2]]}",
                   "{4, {-7.0, 1.0, 2.0, 7.0}, True}", 0);
    /* The element type is found (Integer here) or given; a list of another kind
     * at call time is answered by the interpreter, correctly. */
    assert_eval_eq("c5i = Compile[{{b, _LinkedList}}, b[[1]] = 99; b[[-1]] + First[b]]; "
                   "c5li = LinkedList[{1, 2, 3}]; c5lr = LinkedList[{1., 2.}]; "
                   "{c5i[c5li], Normal[c5li], c5i[c5lr], Normal[c5lr]}",
                   "{102, {99, 2, 3}, 101.0, {99, 2.0}}", 0);
    assert_eval_eq("c5k = Compile[{{b, _LinkedList, _Integer}}, Pop[b] + Pop[b]]; "
                   "{c5k[LinkedList[{4, 5, 6}]], Head[Compile[{{b, _LinkedList, _Integer}}, b][LinkedList[{1}]]]}",
                   "{11, LinkedList}", 0);
    /* Returning the argument returns the same object. */
    assert_eval_eq("c5r = Compile[{{b, _LinkedList, _Real}}, Push[b, 1.]]; c5o = LinkedList[]; "
                   "{c5r[c5o] === c5o, Length[c5o]}", "{True, 1}", 0);
}

/* A compiled call that fails half way leaves no trace on a list it did not
 * create: the interpreter's re-run then makes every change exactly once and
 * emits its message exactly once. */
void test_compile_failure_is_undone() {
    assert_eval_eq("c5f = Compile[{{b, _LinkedList, _Integer}}, (AppendTo[b, 1]; PrependTo[b, 0]; b[[1]] = 7; Pop[b]; b[[10]])]; "
                   "c5q = LinkedList[{1, 2, 3}]; {Quiet[c5f[c5q]], Normal[c5q]}",
                   "{Part[LinkedList[<4>], 10], {7, 1, 2, 3}}", 0);
    assert_eval_eq("c5p = Compile[{{b, _LinkedList, _Integer}}, Pop[b]]; c5e = LinkedList[]; "
                   "{Check[c5p[c5e], failed], Length[c5e]}", "{failed, 0}", 0);
    /* An Integer Total past int64 is a bignum: the interpreter answers. */
    assert_eval_eq("Compile[{{b, _LinkedList, _Integer}}, Total[b]][LinkedList[{2^62, 2^62, 2^62}]]",
                   "13835058055282163712", 0);
}

/* Do / While / For over a GLOBAL list compile as one program. Same answer as
 * the interpreter; a failure part way runs the loop through the interpreter,
 * with the list as the interpreter would leave it. */
void test_autocompile_loops_over_global_lists() {
    long runs0 = autocompile_closed_stmt_runs();
    assert_eval_eq("c5g = LinkedList[]; Do[AppendTo[c5g, N[i]], {i, 10^5}]; "
                   "c5gi = LinkedList[]; Do[AppendTo[c5gi, i^2], {i, 1000}]; Do[PrependTo[c5gi, -i], {i, 3}]; "
                   "c5gz = LinkedList[{1. + I}]; Do[AppendTo[c5gz, c5gz[[-1]] I], {4}]; "
                   "{Length[c5g], Total[c5g], c5gi[[1]], c5gi[[-1]], Normal[c5gz]}",
                   "{100000, 5.00005*^9, -3, 1000000, {1.0 + 1.0*I, -1.0 + 1.0*I, -1.0 - 1.0*I, 1.0 - 1.0*I, 1.0 + 1.0*I}}", 0);
    expect_true(autocompile_closed_stmt_runs() >= runs0 + 4, "the loops ran compiled");
    /* The same loops interpreted give the same lists. */
    autocompile_set_enabled(false);
    assert_eval_eq("c5h1 = LinkedList[]; Do[AppendTo[c5h1, N[i]], {i, 10^5}]; Normal[c5h1] === Normal[c5g]", "True", 0);
    autocompile_set_enabled(true);
    /* Overflow part way (2^63): undone, then interpreted from the start. */
    assert_eval_eq("c5ov = LinkedList[{2^61}]; Do[AppendTo[c5ov, 2 c5ov[[-1]]], {i, 4}]; Normal[c5ov]",
                   "{2305843009213693952, 4611686018427387904, 9223372036854775808, "
                   "18446744073709551616, 36893488147419103232}", 0);
    /* Branches of two types keep the interpreter's types (no widening in a
     * loop that stands in for the interpreter: COMPILE_STRICT_JOIN). */
    assert_eval_eq("c5mx = LinkedList[]; Do[AppendTo[c5mx, If[EvenQ[i], i, N[i]]], {i, 4}]; Normal[c5mx]",
                   "{1.0, 2, 3.0, 4}", 0);
    /* A body that sets a global (acc) is not compiled; still correct. */
    assert_eval_eq("c5st = LinkedList[Range[5]]; c5acc = 0; Do[c5acc = c5acc + Pop[c5st], {3}]; {c5acc, Normal[c5st]}",
                   "{12, {1, 2}}", 0);
    /* A per-element autocompiled body may READ a global list ... */
    assert_eval_eq("c5rd = LinkedList[N[Range[10]]]; Table[c5rd[[i]]^2, {i, 3}]", "{1.0, 4.0, 9.0}", 0);
    /* ... but never mutate one (Table re-runs from the start on a decline). */
    assert_eval_eq("c5tm = LinkedList[]; Table[AppendTo[c5tm, N[i]]; i, {i, 5}]; Normal[c5tm]",
                   "{1.0, 2.0, 3.0, 4.0, 5.0}", 0);
}

void test_compile_lists_are_fast() {
    /* 10^6 appends compiled: tens of milliseconds (interpreted ~0.5 s). */
    assert_eval_eq("c5big = LinkedList[]; c5t = First[AbsoluteTiming[Do[AppendTo[c5big, N[i]], {i, 10^6}]]]; "
                   "{c5t < 0.25, Length[c5big]}", "{True, 1000000}", 0);
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
    TEST(test_append_prepend_to);
    TEST(test_part_assignment);
    TEST(test_mutators_linear_time);
    TEST(test_copying_heads_match_list);
    TEST(test_copying_heads_leave_source);
    TEST(test_copying_messages);
    TEST(test_listable_threading_bit_parity);
    TEST(test_map_select_listable_direct);
    TEST(test_walkers_and_functional);
    TEST(test_structure_and_iterators);
    TEST(test_copying_heads_scale);
    TEST(test_compile_module_local_lists);
    TEST(test_compile_list_arguments);
    TEST(test_compile_failure_is_undone);
    TEST(test_autocompile_loops_over_global_lists);
    TEST(test_compile_lists_are_fast);

    printf("All LinkedList tests passed!\n");
    symtab_clear();
    return 0;
}
