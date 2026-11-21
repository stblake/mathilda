/* test_scoping_semantics.c -- the scoping invariants every substitution fast
 * path must keep (docs/design/scoping_substitution.md, EVAL_CHANGE_PROPOSAL.md
 * section 8.1).
 *
 * Lexical renaming (Module / With / named Function, into held code too),
 * dynamic Block and iterators, shadowing, simultaneous binding, capture
 * avoidance in rule application (A11), the $ModuleNumber / Unique counter,
 * and repeated application after per-node summaries have been filled.  Each
 * expected value is the output of the binary before the summary fast paths
 * (v0.6652), so a fast path that changes any of them -- even to an
 * alpha-equivalent form -- fails here.
 *
 * Run it three ways: as is, with MATHILDA_NO_SYMSUM=1, and with
 * MATHILDA_SYMSUM_DIFF=1 (every fast walk re-checked against the slow one). */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fail = 0;

/* Evaluate `input` (definitions may precede the probe, separated by ';') and
 * compare ToString[InputForm[...]] of the value with `expected`. */
static void probe(const char* input, const char* expected) {
    size_t n = strlen(input) + 64;
    char* src = malloc(n);
    snprintf(src, n, "ToString[InputForm[%s]]", input);
    Expr* e = parse_expression(src);
    if (!e) { fprintf(stderr, "FAIL: does not parse: %s\n", input); abort(); }
    Expr* r = evaluate(e);
    const char* got = (r && r->type == EXPR_STRING) ? r->data.string : "<not a string>";
    if (strcmp(got, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, got);
        g_fail++;
    }
    expr_free(r);
    expr_free(e);
    free(src);
}

/* Run definitions whose value is not checked. */
static void run(const char* input) {
    Expr* e = parse_expression(input);
    if (!e) { fprintf(stderr, "FAIL: does not parse: %s\n", input); abort(); }
    expr_free(evaluate(e));
    expr_free(e);
}

static void test_capture_and_shadowing(void) {
    run("g[v_] := Module[{e = 1}, v + e]");
    probe("g[e + 1]", "2 + e");                     /* A11: the local is renamed, not captured */
    run("sh[e_] := Module[{e = 1}, e + 1]");
    probe("sh[99]", "2");
    run("gb[v_] := Block[{e = 1}, v + e]");
    probe("gb[e + 1]", "3");                        /* Block is dynamic: no renaming */
    run("cap[v_] := Module[{t}, t = v; {t, Hold[v]}]");
    probe("cap[t]", "{t, Hold[t]}");
    /* A scope with no pattern variable inside is still renamed when its local
     * collides with an injected value: the summary skip must not change that. */
    run("fc[x_] := {x, Module[{e}, e]}");
    probe("StringMatchQ[SymbolName[Last[fc[e + 1]]], \"e$*$*\"]", "True");
    run("fc2[x_] := Module[{e}, {x, e}]");
    probe("fc2[e + 1][[1]]", "1 + e");
    probe("fc2[e + 1][[2]] =!= e", "True");
}

static void test_binding_forms(void) {
    probe("With[{q = 12}, With[{k = q}, k]]", "12");
    probe("With[{k = 2}, Hold[Table[k, {k, 3}]]]", "Hold[Table[2, {2, 3}]]");
    run("gs[f_] := Sum[f, {k, 1, 3}]");
    probe("gs[1/k]", "11/6");                       /* iterators are not scoping constructs */
    probe("Function[x, (x + #) &][a]", "a + #1 &");
    probe("Module[{c = 0}, Scan[(c = c + 1) &, {1, 2, 3}]; c]", "3");
    probe("Module[{x}, Hold[x]] /. s_Symbol :> StringMatchQ[SymbolName[s], \"x$*\"]", "False[True]");
    probe("With[{x = 5}, Function[x, x + 1][2]]", "3");
    probe("With[{x = 5}, Function[{y}, x + y][2]]", "7");
    probe("With[{x = 5}, Module[{y = x}, y + x]]", "10");
    probe("Module[{a = 1}, With[{b = a}, Module[{a = 2}, a + b]]]", "3");
    probe("Function[{x}, Function[{x}, x][2] + x][10]", "12");
    probe("With[{a = 1}, Hold[With[{a = 2}, a], a]]", "Hold[With[{a = 2}, a], 1]");
    probe("With[{a = 1}, Hold[Module[{a}, a], a]]", "Hold[Module[{a}, a], 1]");
    probe("With[{Plus = Times}, Plus[2, 3]]", "6");
    probe("Module[{List = 7}, List]", "7");
    probe("Catch[Module[{a = Throw[7]}, a]]", "7");
    probe("Module[{a = 1}, Return[a + 1]; 5]", "2");
}

static void test_large_bodies_and_repetition(void) {
    run("big = Hold @@ {Table[h[k, u, v], {k, 200}]}");
    run("f1[x_] := Module[{t = x}, If[t > 0, t, big; ReleaseHold[big]; u + t]]");
    probe("f1[1]", "1");
    probe("Length[f1[-1]]", "2");
    run("f2[x_] := With[{y = x}, Hold[{y, z, big}]]");
    probe("f2[3][[1, 1 ;; 2]]", "{3, z}");
    run("k1[x_] := Module[{r = x}, Do[r = r + 1, {3}]; If[r > 100, big, r]]");
    probe("Table[k1[i], {i, 5}]", "{4, 5, 6, 7, 8}");
    probe("Table[k1[i], {i, 5}]", "{4, 5, 6, 7, 8}");   /* again, summaries now filled */
    /* Redefinition: summaries of the old right-hand side must not leak. */
    run("k2[x_] := Module[{r = x}, {r, w}]");
    probe("k2[1]", "{1, w}");
    run("k2[x_] := Module[{w = x}, {w, r}]");
    probe("k2[1] /. s_Symbol :> SymbolName[s]", "\"List\"[1, \"r\"]");
}

/* A name seen as a free symbol inside a summarised value BEFORE it is ever
 * used as a local: the epoch must invalidate the old summary. */
static void test_epoch(void) {
    run("vv = Hold[{uu, ww}]; vv; vv");
    run("useV[x_] := {x, vv}; useV[1]; useV[2]");      /* summaries of vv's value get filled */
    run("fe[] := Module[{uu = 1}, {uu, vv}]");
    probe("fe[] /. s_Symbol :> SymbolName[s]", "\"List\"[1, \"Hold\"[\"List\"[\"uu\", \"ww\"]]]");
    run("fe2[] := With[{ww = 9}, {ww, vv}]");
    probe("fe2[]", "{9, Hold[{uu, ww}]}");
    probe("ToExpression[\"Module[{zq9 = 3}, zq9^2]\"]", "9");
    probe("ToExpression[\"Module[{zq10}, Hold[zq10]]\"] /. s_Symbol :> StringMatchQ[SymbolName[s], \"zq10$*\"]",
          "False[True]");
}

/* $ModuleNumber and Unique share one monotone counter, and a user Set of
 * $ModuleNumber is honoured by both. */
static void test_module_number(void) {
    run("$ModuleNumber = 5000");
    probe("Module[{x}, x] /. s_Symbol :> SymbolName[s]", "\"x$5000\"");
    probe("$ModuleNumber", "5001");
    probe("SymbolName[Unique[x]]", "\"x$5001\"");
    probe("Module[{x}, SymbolName[x]]", "\"x$5002\"");
    probe("$ModuleNumber", "5003");
    probe("{SymbolName[Unique[y]], Module[{y}, SymbolName[y]], SymbolName[Unique[y]]}",
          "{\"y$5003\", \"y$5004\", \"y$5005\"}");
    probe("$ModuleNumber", "5006");
    probe("IntegerQ[$ModuleNumber]", "True");
    /* $ModuleNumber read inside a value and a condition sees the live count. */
    run("mnq[] := $ModuleNumber");
    probe("Module[{a}, 1]; m1 = mnq[]; Module[{a}, 1]; mnq[] - m1", "1");
    run("Clear[$ModuleNumber]");
    probe("IntegerQ[$ModuleNumber]", "False");
    probe("Module[{x}, StringMatchQ[SymbolName[x], \"x$*\"]]", "True");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_capture_and_shadowing);
    TEST(test_binding_forms);
    TEST(test_large_bodies_and_repetition);
    TEST(test_epoch);
    TEST(test_module_number);

    if (g_fail) { fprintf(stderr, "%d scoping probe(s) failed\n", g_fail); return 1; }
    printf("All scoping semantics probes passed.\n");
    return 0;
}
