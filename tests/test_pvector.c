/* PVector at the interpreter level.  Interpreter parity: every operation applied
 * to a PVector must equal the same operation on the List holding the same
 * elements -- compared as Normal[op[PVector[l]]] === op[l] for the heads that
 * preserve structure, and op[PVector[l]] === op[l] for reductions and queries.
 * Run over unboxed integer and real storage, boxed mixed storage, an empty
 * vector, and a vector long enough to need several trie levels.  Then the
 * persistence contract (an update never changes another holder's value), value
 * equality/hash (usable as an Association key), pattern matching, and printing. */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

typedef struct { const char* tpl; int wrap; } Case;   /* '@' stands for the argument */
static const Case CASES[] = {
    /* structure-preserving: the answer is a PVector, compared via Normal */
    {"Reverse[@]", 1}, {"Rest[@]", 1}, {"Most[@]", 1}, {"Sort[@]", 1},
    {"Take[@, 2]", 1}, {"Take[@, -2]", 1}, {"Take[@, {2, 3}]", 1}, {"Drop[@, 1]", 1}, {"Drop[@, -2]", 1},
    {"Join[@, @]", 1}, {"Append[@, 7]", 1}, {"Prepend[@, 7]", 1}, {"Append[@, {1, 2}]", 1},
    {"Insert[@, q, 2]", 1}, {"Delete[@, 2]", 1}, {"ReplacePart[@, 1 -> z]", 1},
    {"Map[f, @]", 1}, {"f /@ @", 1}, {"Select[@, Positive]", 1}, {"DeleteDuplicates[@]", 1},
    {"RotateLeft[@, 1]", 1}, {"Riffle[@, 0]", 1}, {"DeleteCases[@, 1]", 1},
    {"Accumulate[@]", 1}, {"Differences[@]", 1}, {"Union[@]", 1},
    {"@ + 1", 1}, {"2 @", 1}, {"@^2", 1}, {"@ + @", 1}, {"Sin[@]", 1},
    /* reductions and queries: the answer is whatever the List gives */
    {"Total[@]", 0}, {"Max[@]", 0}, {"Min[@]", 0}, {"Mean[@]", 0}, {"Count[@, 1]", 0},
    {"Position[@, 1]", 0}, {"MemberQ[@, 1]", 0}, {"FreeQ[@, 1]", 0}, {"Flatten[@]", 0},
    {"Apply[Plus, @]", 0}, {"Fold[Plus, 0, @]", 0}, {"Tally[@]", 0}, {"Length[@]", 0},
    {"First[@]", 0}, {"Last[@]", 0}, {"Dimensions[@]", 0}, {"Normal[@]", 0},
};

static void build(char* out, size_t cap, const char* tpl, const char* arg) {
    size_t k = 0;
    for (const char* p = tpl; *p && k + 1 < cap; p++) {
        if (*p == '@') { size_t l = strlen(arg); memcpy(out + k, arg, l); k += l; }
        else out[k++] = *p;
    }
    out[k] = 0;
}

/* One data set: a plain List literal.  Cases whose List answer is unevaluated or
 * a message (e.g. Rest[{}]) still agree, since both sides evaluate to the same
 * held expression -- Normal of an unevaluated PVector call is the call itself, so
 * only compare when the List side evaluates cleanly (checked in-language). */
static void run_dataset(const char* name, const char* lst) {
    char vbuf[600], lbuf[600], test[1600];
    snprintf(vbuf, sizeof vbuf, "PVector[%s]", lst);
    for (size_t i = 0; i < sizeof CASES / sizeof *CASES; i++) {
        char vexpr[900], lexpr[900];
        build(vexpr, sizeof vexpr, CASES[i].tpl, vbuf);
        build(lexpr, sizeof lexpr, CASES[i].tpl, lst);
        (void)lbuf;
        if (CASES[i].wrap)
            snprintf(test, sizeof test, "Module[{r = %s, l = %s}, If[ListQ[l], PVectorQ[r] && Normal[r] === l, True]]", vexpr, lexpr);
        else
            snprintf(test, sizeof test, "%s === %s", vexpr, lexpr);
        (void)name;
        assert_true(test);
    }
}

static void test_parity_ints(void)   { run_dataset("ints",  "{5, 3, 8, 1, 9, 2}"); }
static void test_parity_reals(void)  { run_dataset("reals", "{1.5, 2.5, -3.5, 1., 4.25}"); }
static void test_parity_mixed(void)  { run_dataset("mixed", "{1, a, 2.5, {1, 2}, \"s\", x + y}"); }
static void test_parity_empty(void)  { run_dataset("empty", "{}"); }
static void test_parity_long(void)   { run_dataset("long",  "Range[1, 200]"); }
static void test_parity_leaf(void)   { run_dataset("leaf",  "Table[Mod[7 i, 11], {i, 33}]"); }   /* 33: one past a leaf */

static void test_core(void) {
    assert_true("PVectorQ[PVector[{1, 2}]] && !PVectorQ[{1, 2}] && !PVectorQ[5]");
    assert_true("Head[PVector[{1}]] === PVector && Length[PVector[{1, 2, 3}]] == 3");
    assert_true("Normal[PVector[]] === {} && Length[PVector[]] == 0");
    assert_true("Normal[PVector[PVector[{1, 2}]]] === {1, 2}");
    assert_true("PVector[Range[5]][[3]] == 3 && PVector[{1., 2.}][[2]] === 2.");
    assert_true("PVector[Range[100]][[-1]] == 100 && Normal[PVector[Range[100]][[10 ;; 12]]] === {10, 11, 12}");
    assert_true("Normal[PVector[Range[10]][[{1, 3, 5}]]] === {1, 3, 5}");
    assert_true("PVector[{{1, 2}, {3, 4}}][[2, 1]] == 3");
    assert_true("Dimensions[PVector[{1, 2, 3}]] === {3}");
    assert_true("MatchQ[PVector[{1}], _PVector] && !MatchQ[PVector[{1}], _List] && !ListQ[PVector[{1}]]");
    assert_true("MatchQ[Range[3][[7]], _Part]");
    assert_true("MatchQ[PVector[Range[3]][[7]], _Part]");                 /* out of range: unevaluated */
    assert_true("MatchQ[PVector[{1, 2}][[0]], PVector]");                  /* head extraction */
}

static void test_persistence(void) {
    assert_true("Module[{v = PVector[Range[100]], w}, w = v; w[[50]] = -1; "
                "v[[50]] == 50 && w[[50]] == -1 && Normal[v] === Range[100]]");
    assert_true("Module[{v = PVector[Range[100]], w}, w = Append[v, 0]; Length[v] == 100 && Length[w] == 101]");
    assert_true("Module[{v = PVector[Range[100]], w}, w = v; AppendTo[w, 7]; "
                "Length[v] == 100 && Length[w] == 101 && w[[101]] == 7]");
    assert_true("Module[{v = PVector[Range[100]], w}, w = v; PrependTo[w, 7]; "
                "Length[v] == 100 && w[[1]] == 7 && w[[2]] == 1]");
    assert_true("Module[{v = PVector[Range[40]], snaps = {}}, "
                "Do[v[[i]] = 0; AppendTo[snaps, v], {i, 40}]; "
                "Total[Normal[snaps[[10]]]] == Total[Range[11, 40]] && Total[Normal[snaps[[40]]]] == 0]");
    assert_true("Module[{v = PVector[{{1, 2}, {3, 4}}], w}, w = v; w[[2, 1]] = 9; Normal[v] === {{1, 2}, {3, 4}} && Normal[w] === {{1, 2}, {9, 4}}]");
    assert_true("Module[{v = PVector[{1, 2, 3}]}, v[[2]] = x; Normal[v] === {1, x, 3}]");        /* typed -> boxed */
    assert_true("Module[{v = PVector[{1, 2, 3}]}, v[[2 ;; 3]] = 0; Normal[v] === {1, 0, 0}]");     /* non-integer spec */
    assert_true("Module[{v = PVector[Range[5]]}, v[[9]] = 1; Normal[v] === Range[5]]");           /* out of range */
    /* a PVector stored inside another structure is unaffected by later updates */
    assert_true("Module[{v = PVector[{1, 2}], h}, h = <|\"k\" -> v|>; v[[1]] = 9; Normal[h[\"k\"]] === {1, 2}]");
}

static void test_equality_hash(void) {
    assert_true("PVector[{1, 2, 3}] === PVector[{1, 2, 3}] && PVector[{1, 2, 3}] =!= PVector[{1, 2, 4}]");
    assert_true("PVector[{1, 2, 3}] =!= {1, 2, 3} && PVector[{1, 2}] =!= PVector[{1, 2, 3}]");
    assert_true("Append[PVector[{1, 2}], 3] === PVector[{1, 2, 3}]");                       /* built two ways */
    assert_true("Module[{v = PVector[{1, 2}]}, v[[1]] = x; v[[1]] = 1; v === PVector[{1, 2}]]");   /* boxed vs typed */
    assert_true("<|PVector[{1, 2}] -> \"a\"|>[PVector[{1, 2}]] === \"a\"");                  /* usable as a key */
    assert_true("Length[DeleteDuplicates[{PVector[{1}], PVector[{1}], PVector[{2}]}]] == 2");
    assert_true("Length[Union[{PVector[{1}], PVector[{1}], PVector[{2}]}]] == 2");
    assert_true("(PVector[{1, 2}] /. PVector[{1, 2}] -> ok) === ok");
    assert_true("Sort[{PVector[{2}], PVector[{1}]}] === {PVector[{1}], PVector[{2}]}");
}

static void test_printing(void) {
    assert_true("ToString[PVector[{1, 2, 3}]] === \"PVector[<3>]\"");
    assert_true("ToString[PVector[]] === \"PVector[<0>]\"");
    assert_true("ToString[PVector[{1, 2, 3}], InputForm] === \"PVector[<3>]\"");
    assert_true("StringContainsQ[ToString[PVector[{1, 2, 3}], TeXForm], \"PVector\"] && StringContainsQ[ToString[PVector[{1, 2, 3}], TeXForm], \"3\"]");
    assert_true("StringContainsQ[ToString[LinkedList[{1, 2}], TeXForm], \"LinkedList\"]");
}

static void test_nesting(void) {
    assert_true("Normal[PVector[{PVector[{1, 2}], 3}][[1]]] === {1, 2}");
    assert_true("Length[{PVector[{1}], PVector[{2}]}] == 2 && Head[{PVector[{1}]}[[1]]] === PVector");
    assert_true("Head[f[PVector[{1}]]] === f && FreeQ[f[PVector[{1}]], Plus]");
    assert_true("Module[{g}, g[x_] := Length[x]; g[PVector[Range[9]]] == 9]");
    assert_true("(Function[v, v[[2]]])[PVector[{4, 5, 6}]] == 5");
    assert_true("Normal[Map[PVector, {{1}, {2, 3}}][[2]]] === {2, 3}");
}
/* Compile[] parity: a PVector argument, result or in-body use must give the answer
 * the interpreter gives (the compiler declines what it cannot lower and the
 * interpreter runs it; it must never answer differently or crash). */
static void test_compile_parity(void) {
    assert_true("Compile[{{a, _Integer, 1}}, Total[a]][PVector[Range[10]]] === 55");
    assert_true("Compile[{x}, Length[x]][PVector[Range[7]]] === 7");
    assert_true("Compile[{x}, x][PVector[{1, 2}]] === PVector[{1, 2}]");
    assert_true("Compile[{{n, _Integer}}, Module[{p = PVector[]}, Do[p = Append[p, i^2], {i, n}]; Total[Normal[p]]]][20] === Total[Range[20]^2]");
    assert_true("Compile[{{n, _Integer}}, Module[{p = PVector[Range[n]], s = 0}, Do[s += p[[i]], {i, n}]; s]][50] === 1275");
    assert_true("Compile[{{n, _Integer}}, Module[{p = PVector[Range[n]]}, Do[p[[i]] = 0, {i, n}]; Total[Normal[p]]]][30] === 0");
    assert_true("Normal[Compile[{{n, _Integer}}, PVector[Range[n]]][5]] === Range[5]");
    assert_true("Compile[{{n, _Integer}}, Module[{p = PVector[Range[n]]}, p[[n]]]][12] === 12");
    /* the same compiled function, many times: no state leaks between calls */
    assert_true("Module[{cf = Compile[{{n, _Integer}}, Module[{p = PVector[]}, Do[p = Append[p, i], {i, n}]; Length[p]]]}, "
                "And @@ Table[cf[k] === k, {k, 0, 60}]]");
    assert_true("Module[{v = PVector[Range[40]], cf = Compile[{{a, _Integer, 1}}, Total[a]]}, "
                "v[[3]] = 100; cf[v] === Total[Normal[v]]]");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_core);
    TEST(test_parity_ints);
    TEST(test_parity_reals);
    TEST(test_parity_mixed);
    TEST(test_parity_empty);
    TEST(test_parity_long);
    TEST(test_parity_leaf);
    TEST(test_persistence);
    TEST(test_equality_hash);
    TEST(test_printing);
    TEST(test_nesting);
    TEST(test_compile_parity);
    printf("All pvector tests passed!\n");
    return 0;
}
