/* Compile[] lowering of PVector (src/compile/compile_pvec.c).
 *
 * The contract under test: a compiled body gives EXACTLY the interpreter's
 * answer, or the VM declines and the interpreter runs the body from the start.
 * Every case therefore runs the same body twice -- through Compile[] and through
 * Function[] (the interpreter) -- and compares the two as printed InputForm
 * (SameQ is tolerant for reals; the printed form is not), after Normal for a
 * PVector result.
 *
 *   - every lowered op, over Integer and Real elements, at the trie boundaries
 *     0, 1, 31, 32, 33, 1056, 1057 elements
 *   - the body really compiled (CompileDiagnostics -> Compiled True)
 *   - persistence inside a compiled body (an old value survives an update)
 *   - arguments: borrowed, never changed; returned as the same value
 *   - declines give the interpreter's answer (mixed type, boxed argument,
 *     out of range, string element, Integer overflow, With constant)
 *   - 200 calls of one compiled function: no state leaks between calls
 *   - the optimiser A/B switches (NO_OPT / NO_CSE / NO_FUSE / NO_PAR) do not
 *     change a result
 */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include "compile/compile.h"
#include "pvec.h"
#include "sym_intern.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

/* The exact comparison used throughout: a PVector is read through Normal, and
 * two answers are the same when their InputForm strings are. */
static void define_helpers(void) {
    assert_eval_eq("norm[x_] := If[PVectorQ[x], Normal[x], x]; 1", "1", 0);
    assert_eval_eq("exact[a_, b_] := ToString[norm[a], InputForm] === ToString[norm[b], InputForm]; 1", "1", 0);
}

/* Expressions are spliced as TEXT: a pattern-substituted Function[{n}, body] is
 * renamed by the evaluator's scoping rules and would no longer bind the n in body. */
static char* fmt_unused(const char* f, const char* a, const char* b) {
    size_t cap = strlen(f) + strlen(a) + (b ? strlen(b) : 0) + 64;
    char* buf = (char*)malloc(cap);
    snprintf(buf, cap, f, a, a, b ? b : "");
    return buf;
}
/* compiled body (one Integer argument n) === the interpreter's, at n = k */
static void cmp_ok(const char* body, const char* k) {
    size_t cap = 2 * strlen(body) + 200;
    char* buf = (char*)malloc(cap);
    snprintf(buf, cap, "exact[Compile[{{n, _Integer}}, %s][%s], Function[{n}, %s][%s]]", body, k, body, k);
    assert_true(buf);
    free(buf);
}
static bool lowers_q(const char* body) {
    size_t cap = strlen(body) + 200;
    char* buf = (char*)malloc(cap);
    snprintf(buf, cap, "(\"Compiled\" /. CompileDiagnostics[{{n, _Integer}}, %s]) === True", body);
    Expr* p = parse_expression(buf);
    Expr* r = evaluate(p);
    bool t = r->type == EXPR_SYMBOL && strcmp(r->data.symbol.name, "True") == 0;
    expr_free(p); expr_free(r); free(buf);
    return t;
}
static void lowers_ok(const char* body)  { ASSERT_MSG(lowers_q(body), "expected to lower: %s", body); }
static void lowers_not(const char* body) { ASSERT_MSG(!lowers_q(body), "expected NOT to lower: %s", body); }
static const char* SIZES[] = {"0", "1", "31", "32", "33", "1056", "1057"};
static void sweep_ok(const char* body) {
    for (size_t i = 0; i < sizeof SIZES / sizeof *SIZES; i++) cmp_ok(body, SIZES[i]);
}
/* Bodies over n.  The Integer-element builder and the Real-element builder are
 * spliced into each shape. */
#define BUILD_I "Module[{v = PVector[]}, Do[v = Append[v, 3*i - 7], {i, n}]; "
#define BUILD_R "Module[{v = PVector[]}, Do[v = Append[v, i/7.], {i, n}]; "

static const struct { const char* name; const char* body; } BODIES[] = {
    {"build ints",            BUILD_I "v]"},
    {"build reals",           BUILD_R "v]"},
    {"AppendTo ints",         "Module[{v = PVector[]}, Do[AppendTo[v, 3*i - 7], {i, n}]; v]"},
    {"AppendTo reals",        "Module[{v = PVector[]}, Do[AppendTo[v, i/7.], {i, n}]; v]"},
    {"Length",                BUILD_I "Length[v]]"},
    {"Part/First/Last ints",  BUILD_I "If[n > 0, v[[n]] + v[[-1]] + First[v] + Last[v] + v[[Quotient[n + 1, 2]]], 0]]"},
    {"Part/First/Last reals", BUILD_R "If[n > 0, v[[n]] + v[[-1]] + First[v] + Last[v], 0.]]"},
    {"Total ints",            BUILD_I "Total[v]]"},
    {"Total reals",           BUILD_R "Total[v]]"},
    {"Normal ints",           BUILD_I "Normal[v]]"},
    {"Normal reals",          BUILD_R "Normal[v]]"},
    {"Part set ints",         BUILD_I "Do[v[[i]] = v[[i]] + 1, {i, n}]; v]"},
    {"Part set reals",        BUILD_R "Do[v[[i]] = v[[i]]*2., {i, n}]; v]"},
    {"Part set negative",     BUILD_I "If[n > 0, v[[-1]] = 5; v[[1]] = -5, 0]; v]"},
    {"Most loop",             BUILD_I "Do[v = Most[v], {k, 1, n - 3}]; v]"},
    {"Most/Append values",    BUILD_I "If[n > 0, Length[Most[Append[v, 5]]] + Total[Append[Most[v], 7]], 0]]"},
    {"PVector[Table] ints",   "Module[{v = PVector[Table[i*2, {i, n}]]}, v = Append[v, 1]; v]"},
    {"PVector[Table] reals",  "Module[{v = PVector[Table[i/3., {i, n}]]}, v = Append[v, 1.5]; v]"},
    {"two vectors",           "Module[{a = PVector[], b = PVector[]}, Do[AppendTo[a, i]; AppendTo[b, 2.5*i], {i, n}]; Total[a] + Total[b] + Length[a]]"},
    {"scalar around vector",  "Module[{v = PVector[], s = 0}, Do[v = Append[v, i]; s += v[[i]], {i, n}]; s + Length[v]]"},
    {"While/For loops",       "Module[{v = PVector[], k = 1}, While[k <= n, AppendTo[v, k*k]; k++]; Do[v[[j]] = v[[j]] - 1, {j, Length[v]}]; v]"},
};

static void test_parity_every_op(void) {
    define_helpers();
    for (size_t i = 0; i < sizeof BODIES / sizeof *BODIES; i++) {
        printf("  lowers + parity: %s\n", BODIES[i].name);
        lowers_ok(BODIES[i].body);
        sweep_ok(BODIES[i].body);
    }
}

/* Compiled results are PVectors, as the interpreter's own representation, and
 * the element storage is the interpreter's: typed, not boxed. */
static void test_result_representation(void) {
    assert_true("PVectorQ[Compile[{{n, _Integer}}, " BUILD_I "v]][40]]");
    assert_true("PVectorQ[Compile[{{n, _Integer}}, " BUILD_I "v]][0]]");
    assert_true("Head[Compile[{{n, _Integer}}, " BUILD_R "v]][40]] === PVector");
    /* Normal presents as a plain List, as pv_to_list does, packed or not */
    assert_true("PackedArrayQ[Compile[{{n, _Integer}}, " BUILD_I "Normal[v]]][2000]] === PackedArrayQ[Function[{n}, " BUILD_I "Normal[v]]][2000]]");
    assert_true("PackedArrayQ[Compile[{{n, _Integer}}, " BUILD_I "Normal[v]]][5]] === PackedArrayQ[Function[{n}, " BUILD_I "Normal[v]]][5]]");
    /* the result is a value of its own: updating it later does not disturb a copy */
    assert_true("Module[{r = Compile[{{n, _Integer}}, " BUILD_I "v]][100], s}, s = r; r[[1]] = 0; Normal[s] === Normal[Function[{n}, " BUILD_I "v]][100]] && r[[1]] === 0]");
}

/* Persistence inside a compiled body: keep an old vector, update, old unchanged. */
static void test_persistence(void) {
    const char* keep =
        "Module[{v = PVector[]}, Do[AppendTo[v, i], {i, n}]; "
        "Module[{old = v}, v = Append[v, 99]; v[[1]] = -1; Do[v = Append[v, i], {i, 5}]; ";
    char buf[1500];
    snprintf(buf, sizeof buf, "%s old]]", keep); lowers_ok(buf); sweep_ok(buf);
    snprintf(buf, sizeof buf, "%s v]]", keep); sweep_ok(buf);
    snprintf(buf, sizeof buf, "%s Total[old]*1000000 + Length[old]*1000 + Length[v] + Total[v]]]", keep);
    lowers_ok(buf); sweep_ok(buf);
    /* the old vector equals what it was, built independently */
    assert_true("Compile[{{n, _Integer}}, Module[{v = PVector[]}, Do[AppendTo[v, i], {i, n}]; "
                "Module[{old = v}, v = Append[v, 99]; v[[1]] = -1; old]]][100] === PVector[Range[100]]");
    /* Most of a copy leaves the original */
    assert_true("Compile[{{n, _Integer}}, Module[{v = PVector[]}, Do[AppendTo[v, i], {i, n}]; "
                "Module[{w = v}, w = Most[w]; w = Most[w]; v]]][40] === PVector[Range[40]]");
    /* snapshots in a loop: each retained copy keeps its own length */
    assert_true("Compile[{{n, _Integer}}, Module[{v = PVector[], s = 0}, Do[v = Append[v, i]; "
                "Module[{snap = v}, v = Append[v, 0]; s += Length[snap]*1000 + Length[v]], {i, n}]; s]][50] === "
                "Function[{n}, Module[{v = PVector[], s = 0}, Do[v = Append[v, i]; "
                "Module[{snap = v}, v = Append[v, 0]; s += Length[snap]*1000 + Length[v]], {i, n}]; s]][50]");
}

/* A PVector argument: borrowed, never changed, returned as the same value. */
static void test_arguments(void) {
    const char* sizes[] = {"1", "31", "32", "33", "1056", "1057"};
    for (size_t i = 0; i < sizeof sizes / sizeof *sizes; i++) {
        char buf[1500];
        snprintf(buf, sizeof buf,
            "Module[{pv = PVector[Table[3 i - 7, {i, %s}]]}, "
            "cf = Compile[{{v, _PVector, _Integer}}, Total[v] + Length[v] + v[[1]] + Last[v] + First[v] + v[[-1]]]; "
            "exact[cf[pv], Total[Normal[pv]] + Length[pv] + pv[[1]] + Last[pv] + First[pv] + pv[[-1]]]]", sizes[i]);
        assert_true(buf);
        snprintf(buf, sizeof buf,
            "Module[{pv = PVector[Table[i/3., {i, %s}]]}, "
            "cf = Compile[{{v, _PVector, _Real}}, Total[v] + Length[v] + v[[1]] + Last[v]]; "
            "exact[cf[pv], Total[Normal[pv]] + Length[pv] + pv[[1]] + Last[pv]]]", sizes[i]);
        assert_true(buf);
        /* the element type is found for an untyped declaration */
        snprintf(buf, sizeof buf,
            "Module[{pv = PVector[Table[i/3., {i, %s}]]}, "
            "exact[Compile[{{v, _PVector}}, Total[v] + v[[1]]][pv], Total[Normal[pv]] + pv[[1]]]]", sizes[i]);
        assert_true(buf);
        /* an update through a copy leaves the caller's vector alone */
        snprintf(buf, sizeof buf,
            "Module[{pv = PVector[Table[i, {i, %s}]], r, c0, cf2}, c0 = Normal[pv]; "
            "cf2 = Compile[{{v, _PVector, _Integer}}, Module[{w = v}, AppendTo[w, 7]; w[[1]] = 0; Total[v]*1000 + Total[w]]]; "
            "r = cf2[pv]; r === Total[c0]*1000 + Total[c0] - c0[[1]] + 7 && Normal[pv] === c0]", sizes[i]);
        assert_true(buf);
        /* returned: the argument itself, and an update of it */
        snprintf(buf, sizeof buf,
            "Module[{pv = PVector[Range[%s]]}, Compile[{{v, _PVector, _Integer}}, v][pv] === pv "
            "&& Compile[{{v, _PVector, _Integer}}, Append[v, 5]][pv] === Append[pv, 5] "
            "&& Compile[{{v, _PVector, _Integer}}, Most[v]][pv] === Most[pv] && Normal[pv] === Range[%s]]", sizes[i], sizes[i]);
        assert_true(buf);
    }
    assert_true("(\"Compiled\" /. CompileDiagnostics[{{v, _PVector, _Integer}}, Total[v] + Length[v]]) === True");
    assert_true("(\"Compiled\" /. CompileDiagnostics[{{v, _PVector}}, Total[v] + Length[v]]) === True");
    assert_true("(\"Compiled\" /. CompileDiagnostics[{{v, _PVector, _Real}}, Append[v, 2.5]]) === True");
    /* an array argument next to a vector */
    assert_true("Module[{pv = PVector[Range[20]], a = Table[1.5 i, {i, 20}]}, "
                "exact[Compile[{{v, _PVector, _Integer}, {a, _Real, 1}}, Total[v]*Total[a]][pv, a], Total[Range[20]]*Total[a]]]");
}

/* Every decline must give the interpreter's answer. */
static void test_declines(void) {
    /* mixed element types: not lowered, answered by the interpreter (boxed) */
    lowers_not("Module[{v = PVector[]}, v = Append[v, 1]; v = Append[v, 2.5]; Length[v]]");
    cmp_ok("Module[{v = PVector[]}, v = Append[v, 1]; v = Append[v, 2.5]; Total[v]]", "3");
    cmp_ok("Module[{v = PVector[{n, 2.5}]}, Total[v]]", "3");
    lowers_not("Module[{v = PVector[{n, 2.5}]}, Total[v]]");
    /* a boxed vector or a wrong-kind vector for the declared element type */
    /* the compiler has no List-literal lowering: PVector[{...}] in a body is the interpreter's */
    lowers_not("Module[{v = PVector[{n, 2, 3}]}, v = Append[v, 4]; Total[v] + v[[1]] + Length[v]]");
    cmp_ok("Module[{v = PVector[{n, 2, 3}]}, v = Append[v, 4]; Total[v] + v[[1]] + Length[v]]", "5");
    assert_true("exact[Compile[{{v, _PVector, _Real}}, Total[v]][PVector[{1, 2.5}]], 3.5]");
    assert_true("exact[Compile[{{v, _PVector, _Real}}, v[[1]]][PVector[{1, 2}]], 1]");
    assert_true("exact[Compile[{{v, _PVector, _Integer}}, Total[v]][PVector[{1., 2.}]], 3.]");
    assert_true("exact[Compile[{{v, _PVector, _Integer}}, Total[v]][PVector[{1, 2.5}]], 3.5]");
    /* a List where a PVector is declared, and a scalar */
    assert_true("exact[Compile[{{v, _PVector, _Integer}}, Length[v]][{1, 2, 3}], Length[{1, 2, 3}]]");
    /* string elements */
    assert_true("exact[Compile[{{v, _PVector, _Integer}}, Length[v]][PVector[{\"a\", \"b\"}]], 2]");
    assert_true("exact[Compile[{{v, _PVector}}, Length[v]][PVector[{\"a\", \"b\"}]], 2]");
    /* out of range: the interpreter's Part::partw answer */
    assert_true("exact[Quiet[Compile[{{v, _PVector, _Integer}, {i, _Integer}}, v[[i]]][PVector[{1, 2, 3}], 7]], "
                "Quiet[Function[{v, i}, v[[i]]][PVector[{1, 2, 3}], 7]]]");
    assert_true("exact[Quiet[Compile[{{v, _PVector, _Integer}, {i, _Integer}}, v[[i]]][PVector[{1, 2, 3}], 0]], "
                "Quiet[Function[{v, i}, v[[i]]][PVector[{1, 2, 3}], 0]]]");
    assert_true("exact[Quiet[Compile[{{v, _PVector, _Integer}, {i, _Integer}}, v[[i]]][PVector[{1, 2, 3}], -4]], "
                "Quiet[Function[{v, i}, v[[i]]][PVector[{1, 2, 3}], -4]]]");
    assert_true("exact[Quiet[Compile[{{n, _Integer}}, Module[{v = PVector[{1, 2}]}, v[[n]] = 9; v]][5]], "
                "Quiet[Function[{n}, Module[{v = PVector[{1, 2}]}, v[[n]] = 9; v]][5]]]");
    assert_true("exact[Quiet[Compile[{{n, _Integer}}, Module[{v = PVector[]}, First[v]]][1]], "
                "Quiet[Function[{n}, Module[{v = PVector[]}, First[v]]][1]]]");
    assert_true("exact[Quiet[Compile[{{n, _Integer}}, Module[{v = PVector[]}, Last[v]]][1]], "
                "Quiet[Function[{n}, Module[{v = PVector[]}, Last[v]]][1]]]");
    assert_true("exact[Quiet[Compile[{{n, _Integer}}, Module[{v = PVector[]}, v = Most[v]; Length[v]]][1]], "
                "Quiet[Function[{n}, Module[{v = PVector[]}, v = Most[v]; Length[v]]][1]]]");
    /* the sum of Integers past int64 is a bignum: the interpreter's answer */
    assert_true("exact[Compile[{{v, _PVector, _Integer}}, Total[v]][PVector[{9223372036854775807, 5}]], 9223372036854775812]");
    /* Total of an empty Real vector is the Integer 0, as in the interpreter */
    assert_true("exact[Compile[{{n, _Integer}}, Total[PVector[]]][1], Function[{n}, Total[PVector[]]][1]]");
    cmp_ok("Module[{v = PVector[]}, Do[v = Append[v, 1.5 i], {i, n}]; Total[v]]", "0");
    cmp_ok("Module[{v = PVector[]}, Normal[v]]", "0");
    /* a With constant cannot be assigned in the interpreter: not lowered */
    lowers_not("With[{v = PVector[{1, 2}]}, v = Append[v, 3]; Length[v]]");
    /* a PVector through a function application, a conditional join, a loop
     * value or a higher-order head: declined, never miscompiled */
    cmp_ok("Module[{v = PVector[{1, 2}]}, Length[If[n > 1, v, Append[v, 1]]]]", "3");
    cmp_ok("Module[{v = PVector[{1, 2}]}, Length[Nest[Append[#, 1] &, v, n]]]", "3");
    cmp_ok("Module[{v = PVector[{1, 2}]}, Total[Map[#^2 &, v]]]", "3");
    cmp_ok("Module[{v = PVector[{1, 2}]}, Total[Table[v, {2}][[1]]]]", "3");
    /* an assignment inside another operand changes the evaluation order: declined */
    cmp_ok("Module[{v = PVector[{1, 2, 3}]}, Append[v, (v = Most[v]; 1)]; Length[v]]", "0");
    cmp_ok("Module[{v = PVector[{1, 2, 3}]}, v[[(v = Most[v]; 1)]] + Length[v]]", "0");
    /* a global vector is not folded: answered by the interpreter */
    assert_true("Module[{}, gv = PVector[Range[10]]; exact[Compile[{{n, _Integer}}, Length[gv] + n][3], 13]]");
    /* arithmetic on a vector is not a machine scalar */
    cmp_ok("Module[{v = PVector[{1, 2}]}, Total[v + 1]]", "0");
}

/* The same compiled function over many calls: no state leaks between them. */
static void test_many_calls(void) {
    assert_true("Module[{cf = Compile[{{n, _Integer}}, Module[{v = PVector[], old = PVector[{7}]}, "
                "Do[v = Append[v, i*i], {i, n}]; old = Append[old, n]; Total[v] + Length[old]*1000 + old[[2]]]]}, "
                "And @@ Table[cf[k] === Function[{n}, Module[{v = PVector[], old = PVector[{7}]}, "
                "Do[v = Append[v, i*i], {i, n}]; old = Append[old, n]; Total[v] + Length[old]*1000 + old[[2]]]][k], {k, 1, 200}]]");
    /* a vector argument handed to 200 calls is the same afterwards */
    assert_true("Module[{pv = PVector[Range[64]], cf = Compile[{{v, _PVector, _Integer}, {k, _Integer}}, "
                "Module[{w = v}, w[[Mod[k, 64] + 1]] = k; AppendTo[w, k]; Total[w] + Length[w]]], r}, "
                "r = Table[cf[pv, k], {k, 1, 200}]; "
                "Normal[pv] === Range[64] && r === Table[Total[ReplacePart[Range[64], Mod[k, 64] + 1 -> k]] + k + 65, {k, 1, 200}]]");
    /* and the results of earlier calls are untouched by later ones */
    assert_true("Module[{cf = Compile[{{n, _Integer}}, Module[{v = PVector[]}, Do[v = Append[v, i], {i, n}]; v]], rs}, "
                "rs = Table[cf[k], {k, 0, 60}]; And @@ Table[Normal[rs[[k + 1]]] === Range[k], {k, 0, 60}]]");
}

/* ---- the optimiser A/B ---------------------------------------------------- */

static bool run_flags(const char* src, unsigned flags, long long arg, Expr** result) {
    Expr* body = parse_expression(src);
    const char* names[1] = { intern_symbol("n") };
    CompileType types[1] = { CT_INT };
    CompiledProgram* p = compile_expr_ex(body, names, types, 1, flags);
    if (!p) { expr_free(body); return false; }
    CompileValue av, out;
    memset(&av, 0, sizeof av); memset(&out, 0, sizeof out);
    av.type = CT_INT; av.v.i = arg;
    bool ok = compiled_eval(p, &av, &out);
    *result = NULL;
    if (ok) {
        if (out.type == CT_INT) *result = expr_new_integer(out.v.i);
        else if (out.type == CT_REAL) *result = expr_new_real(out.v.r);
        else if (CT_IS_PVEC(out.type) || CT_IS_ARRAY(out.type)) *result = out.v.a;
    }
    compiled_free(p);
    expr_free(body);
    return ok;
}

static void test_optimiser_ab(void) {
    static const char* bodies[] = {
        BUILD_I "Total[v] + Length[v]]",
        BUILD_R "Total[v]]",
        BUILD_I "v]",
        BUILD_R "Do[v[[i]] = v[[i]]*2. + v[[i]], {i, n}]; v]",
        "Module[{v = PVector[], s = 0}, Do[v = Append[v, i]; s += v[[i]] + Length[v], {i, n}]; s]",
        "Module[{v = PVector[]}, Do[AppendTo[v, i], {i, n}]; Module[{old = v}, v = Append[v, 9]; Length[old]*100 + Length[v] + Length[v] + Length[v]]]",
        "Module[{v = PVector[Table[i, {i, n}]]}, Do[v = Most[v], {k, 1, n - 3}]; v]",
    };
    static const unsigned flags[] = { COMPILE_NO_OPT, COMPILE_NO_CSE, COMPILE_NO_FUSE, COMPILE_NO_PAR,
                                      COMPILE_NO_OPT | COMPILE_NO_CSE | COMPILE_NO_FUSE | COMPILE_NO_PAR };
    static const long long args[] = { 0, 1, 31, 32, 33, 100, 1057 };
    int compared = 0;
    for (size_t b = 0; b < sizeof bodies / sizeof *bodies; b++) {
        for (size_t a = 0; a < sizeof args / sizeof *args; a++) {
            Expr* ref = NULL;
            bool rok = run_flags(bodies[b], 0, args[a], &ref);
            ASSERT_MSG(rok || args[a] < 4, "default flags did not compile/run body %zu at n=%lld", b, args[a]);
            for (size_t f = 0; f < sizeof flags / sizeof *flags; f++) {
                Expr* got = NULL;
                bool gok = run_flags(bodies[b], flags[f], args[a], &got);
                ASSERT_MSG(gok == rok, "body %zu n=%lld flags %u: compiled-or-not differs", b, args[a], flags[f]);
                if (rok) {
                    ASSERT_MSG(got && expr_eq(ref, got), "body %zu n=%lld flags %u: result differs", b, args[a], flags[f]);
                    compared++;
                }
                if (got) expr_free(got);
            }
            if (ref) expr_free(ref);
        }
    }
    ASSERT_MSG(compared > 100, "A/B compared only %d results", compared);
    printf("  optimiser A/B: %d results identical\n", compared);
}

/* The body not only compiles: at run time the VM does NOT decline (compiled_eval
 * returns true), except for the documented empty-vector cases the interpreter
 * answers itself. */
static void test_runs_compiled(void) {
    static const long long ns[] = {0, 1, 31, 32, 33, 1056, 1057};
    int ran = 0;
    for (size_t b = 0; b < sizeof BODIES / sizeof *BODIES; b++) {
        for (size_t i = 0; i < sizeof ns / sizeof *ns; i++) {
            Expr* r = NULL;
            bool ok = run_flags(BODIES[b].body, 0, ns[i], &r);
            bool empty_decline = ns[i] == 0 &&
                (strcmp(BODIES[b].name, "Total reals") == 0 || strcmp(BODIES[b].name, "two vectors") == 0 || strncmp(BODIES[b].name, "Normal", 6) == 0);
            ASSERT_MSG(ok != empty_decline, "%s at n=%lld: compiled run %s", BODIES[b].name, ns[i],
                       ok ? "succeeded but was expected to decline" : "declined");
            if (r) expr_free(r);
            ran += ok;
        }
    }
    printf("  %d compiled runs completed without declining\n", ran);
}

/* Vector ARGUMENTS reach the VM (no decline), borrowed and unchanged. */
static void test_args_run_compiled(void) {
    static const char* srcs[] = {
        "Total[v] + Length[v] + v[[1]] + Last[v]",
        "Module[{w = v}, AppendTo[w, @3]; w[[1]] = @0; Total[w] + Total[v]]",
        "Append[v, @9]",
        "Most[v]",
    };
    static const long long sizes[] = {1, 2, 31, 32, 33, 1056, 1057};
    for (int kind = 0; kind < 2; kind++) {
        for (size_t b = 0; b < sizeof srcs / sizeof *srcs; b++) {
            char src[300];
            {   /* '@' stands for the element kind: 3 -> 3 or 3. */
                size_t o = 0;
                for (const char* q = srcs[b]; *q; q++) {
                    if (*q == '@') continue;
                    src[o++] = *q;
                    if (kind && q > srcs[b] && q[-1] == '@') src[o++] = '.';
                }
                src[o] = 0;
            }
            Expr* body = parse_expression(src);
            const char* names[1] = { intern_symbol("v") };
            CompileType types[1] = { CT_PVEC_TYPE(kind ? CT_REAL : CT_INT) };
            CompiledProgram* p = compile_expr_ex(body, names, types, 1, 0);
            ASSERT_MSG(p != NULL, "argument body %zu (%s) did not compile", b, kind ? "Real" : "Integer");
            for (size_t i = 0; i < sizeof sizes / sizeof *sizes; i++) {
                char mk[200];
                snprintf(mk, sizeof mk, kind ? "PVector[Table[i/3., {i, %lld}]]" : "PVector[Table[3 i - 7, {i, %lld}]]", sizes[i]);
                Expr* mp = parse_expression(mk);
                Expr* pv = evaluate(mp);
                expr_free(mp);
                Expr* before = expr_copy(pv);
                Expr* snap = pv_to_list(pv->data.pvec);
                CompileValue av, out;
                memset(&av, 0, sizeof av); memset(&out, 0, sizeof out);
                av.type = types[0]; av.v.a = pv;
                bool ok = compiled_eval(p, &av, &out);
                ASSERT_MSG(ok, "argument body %zu n=%lld: the VM declined", b, sizes[i]);
                Expr* after = pv_to_list(pv->data.pvec);
                ASSERT_MSG(expr_eq(snap, after), "argument body %zu n=%lld: the argument changed", b, sizes[i]);
                if (CT_IS_PVEC(out.type) || CT_IS_ARRAY(out.type)) expr_free(out.v.a);
                expr_free(snap); expr_free(after); expr_free(before); expr_free(pv);
            }
            compiled_free(p);
            expr_free(body);
        }
    }
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    TEST(test_parity_every_op);
    TEST(test_result_representation);
    TEST(test_persistence);
    TEST(test_arguments);
    TEST(test_declines);
    TEST(test_many_calls);
    TEST(test_optimiser_ab);
    TEST(test_runs_compiled);
    TEST(test_args_run_compiled);
    printf("All compile-pvec tests passed!\n");
    return 0;
}
