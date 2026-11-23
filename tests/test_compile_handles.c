/* Compile[] handle kinds through the CompileHandleOps registry
 * (src/compile/handle_ops.{h,c}).
 *
 * Two parts.
 *
 *  1. The registry rows: window, element accessor, display name, flags and
 *     retain opcode of Array / Association / LinkedList / PVector, and that
 *     scalars and managed types are not handles.
 *
 *  2. Coverage probes: one body per lowered shape of each handle kind.  Every
 *     probe must (a) really compile (CompileDiagnostics -> Compiled True) and
 *     (b) equal the interpreter, compared as the printed InputForm of
 *     {Head, Normal} of the result (a list argument is rebuilt for each side, so
 *     in-place mutation by one side never reaches the other).  Only shapes that
 *     already lowered when this file was written are listed; a shape that does
 *     not compile belongs in the lowering work, not here.
 */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include "compile/compile.h"
#include "compile/handle_ops.h"
#include "compile/compile_emit.h"   /* num_common */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

/* ---- 1. registry rows ---------------------------------------------------- */

static void test_registry_rows(void) {
    CompileType arr2 = CT_ARRAY(CT_REAL, 2);
    CompileType as = CT_ASSOC_TYPE(CT_INT);
    CompileType ll = CT_LLIST_TYPE(CT_COMPLEX);
    CompileType pv = CT_PVEC_TYPE(CT_REAL);

    ASSERT(CT_IS_HANDLE(CT_ARR) && CT_IS_HANDLE(arr2) && CT_IS_HANDLE(as)
           && CT_IS_HANDLE(ll) && CT_IS_HANDLE(pv));
    ASSERT(!CT_IS_HANDLE(CT_BOOL) && !CT_IS_HANDLE(CT_INT) && !CT_IS_HANDLE(CT_REAL)
           && !CT_IS_HANDLE(CT_COMPLEX) && !CT_IS_HANDLE(CT_BIGINT)
           && !CT_IS_HANDLE(CT_BIGREAL) && !CT_IS_HANDLE(CT_BIGCOMPLEX)
           && !CT_IS_HANDLE(CT_ERR));
    /* every encoding in every window belongs to exactly the macro's kind */
    for (int t = 0; t < 64; t++) {
        const CompileHandleOps* h = ct_handle_ops((CompileType)t);
        bool any = CT_IS_ARRAY(t) || CT_IS_ASSOC(t) || CT_IS_LLIST(t) || CT_IS_PVEC(t);
        ASSERT_MSG((h != NULL) == any, "window mismatch at %d", t);
        if (!h) continue;
        ASSERT(ct_handle_elem((CompileType)t) ==
               (CT_IS_ARRAY(t) ? CT_ELEM(t) : CT_IS_ASSOC(t) ? CT_ASSOC_VALTYPE(t)
                : CT_IS_LLIST(t) ? CT_LLIST_ELEM(t) : CT_PVEC_ELEM(t)));
    }
    ASSERT(strcmp(ct_handle_name(arr2), "Array") == 0);
    ASSERT(strcmp(ct_handle_name(as), "Association") == 0);
    ASSERT(strcmp(ct_handle_name(ll), "LinkedList") == 0);
    ASSERT(strcmp(ct_handle_name(pv), "PVector") == 0);
    ASSERT(ct_handle_name(CT_REAL) == NULL);

    ASSERT(!ct_is_opaque_handle(arr2) && ct_is_opaque_handle(as)
           && ct_is_opaque_handle(ll) && ct_is_opaque_handle(pv) && !ct_is_opaque_handle(CT_INT));
    ASSERT(ct_has_value_semantics(pv) && !ct_has_value_semantics(ll)
           && !ct_has_value_semantics(as) && !ct_has_value_semantics(arr2));
    ASSERT(ct_is_scope_local(arr2) && ct_is_scope_local(ll) && ct_is_scope_local(pv)
           && !ct_is_scope_local(as));
    ASSERT(ct_retain_op(ll) != 0 && ct_retain_op(pv) != 0 && ct_retain_op(ll) != ct_retain_op(pv)
           && ct_retain_op(as) == 0 && ct_retain_op(arr2) == 0 && ct_retain_op(CT_INT) == 0);
    /* an open element (CT_BOOL) is only legal for the list kinds, with a reason */
    ASSERT(ct_handle_ops(CT_LLIST)->open_choices == 3 && ct_handle_ops(CT_PVEC)->open_choices == 2
           && ct_handle_ops(CT_ASSOC)->open_choices == 0 && ct_handle_ops(CT_ARR)->open_choices == 0);
    ASSERT(ct_handle_ops(CT_LLIST)->open_bail && ct_handle_ops(CT_PVEC)->open_bail);
    /* a handle is never a numeric operand */
    ASSERT((int)num_common(as, CT_INT) < 0 && (int)num_common(CT_REAL, ll) < 0
           && (int)num_common(pv, pv) < 0 && (int)num_common(arr2, CT_INT) >= 0);
}

/* ---- 2. coverage probes --------------------------------------------------- */

typedef struct {
    const char* kind;   /* for the log */
    const char* spec;   /* Compile argument spec */
    const char* names;  /* the same parameters as a Function slot list */
    const char* body;
    const char* args;   /* call arguments, rebuilt per side */
} Probe;

#define AS  "{{p, _Association, _Real}}", "{p}"
#define ASI "{{p, _Association, _Integer}}", "{p}"
#define LLI "{{b, _LinkedList, _Integer}}", "{b}"
#define LLR "{{b, _LinkedList, _Real}}", "{b}"
#define PVI "{{v, _PVector, _Integer}}", "{v}"
#define PVR "{{v, _PVector, _Real}}", "{v}"
#define NI  "{{n, _Integer}}", "{n}"

static const Probe PROBES[] = {
    /* Association (read-only bag, owned KeyDrop/KeyTake results) */
    {"Association", AS,  "Lookup[p, \"a\"]",                      "<|\"a\" -> 1.5, \"b\" -> 2.|>"},
    {"Association", AS,  "Lookup[p, \"zz\", -1.]",                "<|\"a\" -> 1.5|>"},
    {"Association", AS,  "KeyExistsQ[p, \"a\"]",                  "<|\"a\" -> 1.5|>"},
    {"Association", AS,  "Length[p]",                             "<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>"},
    {"Association", AS,  "Total[Values[p]]",                      "<|\"a\" -> 1., \"b\" -> 2.5|>"},
    {"Association", AS,  "Values[p]",                             "<|\"a\" -> 1., \"b\" -> 2.5|>"},
    {"Association", AS,  "KeyDrop[p, \"b\"]",                     "<|\"a\" -> 1., \"b\" -> 2.5|>"},
    {"Association", AS,  "KeyTake[p, {\"a\", \"c\"}]",            "<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>"},
    {"Association", ASI, "Lookup[p, 1] - Lookup[p, 3]",           "<|1 -> 10, 3 -> 4|>"},
    {"Association", ASI, "Module[{s = 0}, Do[s += Lookup[p, k, 0], {k, 1, 5}]; s]", "<|1 -> 10, 3 -> 4, 5 -> 7|>"},
    /* LinkedList (reference semantics, in place) */
    {"LinkedList",  NI,  "Module[{b = LinkedList[]}, Do[AppendTo[b, i*i], {i, n}]; b]",          "6"},
    {"LinkedList",  NI,  "Module[{b = LinkedList[]}, Do[PrependTo[b, i], {i, n}]; b]",           "6"},
    {"LinkedList",  NI,  "Module[{b = LinkedList[]}, Do[AppendTo[b, i], {i, n}]; Length[b] + Total[b]]", "7"},
    {"LinkedList",  NI,  "Module[{b = LinkedList[]}, Do[AppendTo[b, N[i]], {i, n}]; Normal[b]]", "5"},
    {"LinkedList",  NI,  "Module[{b = LinkedList[]}, Do[AppendTo[b, i], {i, n}]; b[[1]] = 99; b[[-1]] + First[b]]", "4"},
    {"LinkedList",  NI,  "Module[{b = LinkedList[]}, Do[AppendTo[b, i], {i, n}]; Pop[b] + Pop[b]]", "5"},
    {"LinkedList",  NI,  "Module[{a = LinkedList[]}, Do[AppendTo[a, i], {i, n}]; Module[{b = a}, AppendTo[b, 100]; Length[a]]]", "4"},
    {"LinkedList",  LLI, "AppendTo[b, 7]; PrependTo[b, 0]; Total[b]",                              "LinkedList[{1, 2, 3}]"},
    {"LinkedList",  LLI, "b",                                                                       "LinkedList[{4, 5}]"},
    {"LinkedList",  LLR, "AppendTo[b, 2.5]; Length[b]",                                             "LinkedList[{1., 2.}]"},
    /* PVector (persistent value) */
    {"PVector",     NI,  "Module[{v = PVector[]}, Do[v = Append[v, 3*i - 7], {i, n}]; v]",         "40"},
    {"PVector",     NI,  "Module[{v = PVector[]}, Do[AppendTo[v, i/7.], {i, n}]; v]",              "40"},
    {"PVector",     NI,  "Module[{v = PVector[]}, Do[AppendTo[v, i], {i, n}]; Total[v] + Length[v]]", "33"},
    {"PVector",     NI,  "Module[{v = PVector[Table[i*2, {i, n}]]}, v[[1]] = -1; v]",              "33"},
    {"PVector",     NI,  "Module[{v = PVector[]}, Do[AppendTo[v, i], {i, n}]; Module[{old = v}, v = Append[v, 99]; v[[1]] = -1; old]]", "32"},
    {"PVector",     NI,  "Module[{v = PVector[]}, Do[v = Append[v, i], {i, n}]; Do[v = Most[v], {k, 1, n - 3}]; v]", "35"},
    {"PVector",     PVI, "Total[v] + Length[v] + v[[1]] + Last[v]",                                 "PVector[Range[40]]"},
    {"PVector",     PVI, "Append[v, 5]",                                                            "PVector[Range[5]]"},
    {"PVector",     PVI, "v",                                                                       "PVector[Range[5]]"},
    {"PVector",     PVR, "Module[{w = v}, w[[1]] = 0.; Total[v]*1000. + Total[w]]",                 "PVector[{1.5, 2.5, 3.}]"},
};

static void run_probe(const Probe* pr) {
    char buf[4096];
    /* the body really compiled */
    snprintf(buf, sizeof buf, "(\"Compiled\" /. CompileDiagnostics[%s, %s]) === True", pr->spec, pr->body);
    Expr* p = parse_expression(buf);
    Expr* r = evaluate(p);
    bool lowers = r->type == EXPR_SYMBOL && strcmp(r->data.symbol.name, "True") == 0;
    expr_free(p); expr_free(r);
    ASSERT_MSG(lowers, "[%s] expected to lower: %s", pr->kind, pr->body);
    /* equal to the interpreter, as printed InputForm */
    snprintf(buf, sizeof buf, "hp[Compile[%s, %s][%s]] === hp[Function[%s, %s][%s]]",
             pr->spec, pr->body, pr->args, pr->names, pr->body, pr->args);
    Expr* q = parse_expression(buf);
    Expr* s = evaluate(q);
    bool same = s->type == EXPR_SYMBOL && strcmp(s->data.symbol.name, "True") == 0;
    expr_free(q); expr_free(s);
    ASSERT_MSG(same, "[%s] compiled != interpreter: %s", pr->kind, pr->body);
}

static void test_probes(void) {
    /* hp: {Head, Normal} printed -- exact where SameQ would tolerate reals */
    assert_eval_eq("hp[x_] := ToString[{Head[x], Normal[x]}, InputForm]; 1", "1", 0);
    for (size_t i = 0; i < sizeof PROBES / sizeof *PROBES; i++) {
        printf("  lowers + parity: [%s] %s\n", PROBES[i].kind, PROBES[i].body);
        run_probe(&PROBES[i]);
    }
}

/* A body that mixes a handle into scalar arithmetic or a function application
 * declines (the interpreter runs it); the registry must keep refusing. */
static void test_handle_gates(void) {
    assert_eval_eq("(\"Compiled\" /. CompileDiagnostics[{{p, _Association, _Real}}, p + 1.])", "False", 0);
    assert_eval_eq("(\"Compiled\" /. CompileDiagnostics[{{b, _LinkedList, _Integer}}, b + 1])", "False", 0);
    assert_eval_eq("(\"Compiled\" /. CompileDiagnostics[{{v, _PVector, _Integer}}, v + 1])", "False", 0);
    assert_eval_eq("(\"Compiled\" /. CompileDiagnostics[{{v, _PVector, _Integer}}, Module[{k = 1}, Map[Function[x, x], v]]])", "False", 0);
    assert_eval_eq("(\"Compiled\" /. CompileDiagnostics[{{v, _PVector, _Integer}}, With[{w = PVector[]}, w]])", "False", 0);
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    TEST(test_registry_rows);
    TEST(test_probes);
    TEST(test_handle_gates);
    printf("All compile-handles tests passed!\n");
    return 0;
}
