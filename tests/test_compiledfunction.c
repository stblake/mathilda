/* Tests for the user-facing Compile[] / CompiledFunction object (M1b).
 *
 * These drive the whole pipeline (parse → evaluate → apply), unlike
 * test_compile.c which unit-tests the compile_expr engine directly.  Float
 * results are checked via a Mathilda equality/tolerance expression that
 * evaluates to True, so no fragile decimal-string comparisons are needed;
 * integer and symbolic results are compared by exact printed form. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/* Integer-typed compiled functions return exact Integers. */
void test_cf_integer(void) {
    assert_eval_eq("Compile[{{n, _Integer}}, n^2 + 1][5]", "26", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, a b - b][6, 7]", "35", 0);
}

/* Real path, default (bare-symbol) type, and a machine-precision identity. */
void test_cf_real(void) {
    assert_eval_eq("Compile[{{x, _Real}}, x^2 + 1][3.0] == 10", "True", 0);
    assert_eval_eq("Compile[{x}, x + x][2.5] == 5", "True", 0);          /* default _Real */
    assert_eval_eq("Compile[{{x, _Real}}, x^2][3] == 9", "True", 0);      /* Integer arg → Real */
    assert_eval_eq("Abs[Compile[{{x, _Real}}, Sin[x]^2 + Cos[x]^2][1.234] - 1] < 10^-12", "True", 0);
}

/* Complex arguments and results. */
void test_cf_complex(void) {
    assert_eval_eq("Chop[Compile[{{z, _Complex}}, z^2][1.0 + 2.0 I] - (-3 + 4 I)] == 0", "True", 0);
}

/* Symbolic argument → interpreter fallback (still produces the right value). */
void test_cf_symbolic_fallback(void) {
    assert_eval_eq("Compile[{{x, _Real}}, x^2 + 1][a]", "1 + a^2", 0);
}

/* Body outside the compilable subset → the object is still built and
 * application falls back to the interpreter.
 *
 * The marker is a user DownValue, not a head that happens to lack a kernel:
 * this test used `Zeta` until Zeta got one, at which point it was quietly
 * exercising the COMPILED path and no longer testing the fallback at all.
 * CompileDiagnostics below asserts the marker really does bail. */
void test_cf_uncompilable_fallback(void) {
    assert_eval_eq("uncid[t_] := t; Compile[{{x, _Real}}, uncid[x^2]][3] == 9", "True", 0);
    assert_eval_eq("uncid2[t_] := t; "
                   "CompileDiagnostics[{{x, _Real}}, uncid2[x^2]][[1]] == (\"Compiled\" -> False)",
                   "True", 0);
}

/* Indexed Part through the user-facing object (M3c).
 *
 * The engine-level battery in test_compile.c drives compile_expr directly on
 * NDArrays; what is only reachable from here is the BOUNDARY — a List argument
 * packed on the way in and unpacked on the way out, an array built inside the
 * body with no array argument to take its kind from, and the interpreter
 * fallback when a subscript is out of range. */
void test_cf_part(void) {
    /* Scalar subscripts, including from the end and computed. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, v[[2]]][{10., 20., 30.}] == 20", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, v[[-1]]][{10., 20., 30.}] == 30", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}, {k, _Integer}}, v[[k]]][{10., 20., 30.}, 3] == 30",
                   "True", 0);
    assert_eval_eq("Compile[{{m, _Real, 2}}, m[[2, 1]]][{{1., 2.}, {3., 4.}}] == 3", "True", 0);

    /* Slices, All, position lists, partial indexing.  A List goes in, so a
     * List must come back. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, v[[2 ;; 3]]][{10., 20., 30., 40.}]",
                   "{20.0, 30.0}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, v[[{1, 3}]]][{10., 20., 30., 40.}]",
                   "{10.0, 30.0}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, v[[1 ;; 4 ;; 2]]][{10., 20., 30., 40.}]",
                   "{10.0, 30.0}", 0);
    assert_eval_eq("Compile[{{m, _Real, 2}}, m[[2]]][{{1., 2.}, {3., 4.}}]", "{3.0, 4.0}", 0);
    assert_eval_eq("Compile[{{m, _Real, 2}}, m[[All, 1]]][{{1., 2.}, {3., 4.}}]", "{1.0, 3.0}", 0);

    /* Out of range: the compiled call fails and the INTERPRETER answers, which
     * for Part means leaving it unevaluated rather than inventing an element. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, v[[9]]][{1., 2.}]", "Part[{1.0, 2.0}, 9]", 0);
    assert_eval_eq("Compile[{{m, _Real, 2}}, m[[1, 5]]][{{1., 2.}, {3., 4.}}]",
                   "Part[{{1.0, 2.0}, {3.0, 4.0}}, 1, 5]", 0);

    /* Assignment into a local, returned as the result.  With no array ARGUMENT
     * the result has no kind to inherit, and the interpreter running the same
     * body returns a List, so this must too. */
    assert_eval_eq("Compile[{{n, _Integer}}, Module[{u = ConstantArray[0., n]}, "
                   "Do[u[[i]] = 1. i, {i, 1, n}]; u]][4]", "{1.0, 2.0, 3.0, 4.0}", 0);
    assert_eval_eq("Head[Compile[{{n, _Integer}}, Module[{u = ConstantArray[0., n]}, u]][3]]",
                   "List", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Module[{u = v}, u[[2]] = 99.; u]][{1., 2., 3.}]",
                   "{1.0, 99.0, 3.0}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Module[{u = v}, u[[2 ;; 3]] = 0.; u]]"
                   "[{1., 2., 3.}]", "{1.0, 0.0, 0.0}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Module[{u = v}, u[[1]] += 5.; u]][{1., 2.}]",
                   "{6.0, 2.0}", 0);

    /* The local is a COPY: an NDArray passed in must come back untouched, or a
     * compiled body would be mutating a value its caller still owns. */
    assert_eval_eq("nd = NDArray[{1., 2., 3.}]; "
                   "Compile[{{v, _Real, 1}}, Module[{u = v}, u[[1]] = 99.; u]][nd]; Normal[nd]",
                   "{1.0, 2.0, 3.0}", 0);

    /* Writing through an argument is not in the subset, so the object still
     * builds and the interpreter answers — which is where value semantics live. */
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Real, 1}}, v[[1]] = 0.]",
                   "False", 0);

    /* A 5-point stencil: the shape this feature exists for.  Compiled and
     * interpreted must agree element for element. */
    assert_eval_eq(
        "stencil[a_] := Module[{n = Length[a], b = a}, "
        "  Do[b[[i, j]] = (a[[i - 1, j]] + a[[i + 1, j]] + a[[i, j - 1]] + a[[i, j + 1]])/4, "
        "     {i, 2, n - 1}, {j, 2, n - 1}]; b]; "
        "cs = Compile[{{a, _Real, 2}}, Module[{n = Length[a], b = a}, "
        "  Do[b[[i, j]] = (a[[i - 1, j]] + a[[i + 1, j]] + a[[i, j - 1]] + a[[i, j + 1]])/4, "
        "     {i, 2, n - 1}, {j, 2, n - 1}]; b]]; "
        "g = Table[1.0 (10 i + j), {i, 1, 5}, {j, 1, 5}]; "
        "cs[g] === stencil[g]", "True", 0);
}

/* CompileDiagnostics: the answer to "why is this slow".  A bail is otherwise
 * invisible — same result, 10-40x slower — and because the compilable subset is
 * a cliff, the useful report is the single INNERMOST subexpression that stopped
 * the whole body, not the body itself. */
void test_compile_diagnostics(void) {
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, Sin[x] + x^2]", "True", 0);
    assert_eval_eq("\"ResultType\" /. CompileDiagnostics[{{x, _Real}}, Sin[x]]", "\"Real\"", 0);
    assert_eval_eq("\"ResultType\" /. CompileDiagnostics[{{x, _Real}}, x > 0]", "\"Boolean\"", 0);
    assert_eval_eq("\"ResultType\" /. CompileDiagnostics[{{z, _Complex}}, z^2]", "\"Complex\"", 0);
    assert_eval_eq("(\"Instructions\" /. CompileDiagnostics[{{x, _Real}}, Sin[x] + x^2]) > 0",
                   "True", 0);

    /* The innermost cause, not the enclosing construct: the report must name
     * BarnesG[x], not the Plus that contains it. */
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, Sin[x] + BarnesG[x]]",
                   "False", 0);
    assert_eval_eq("\"Subexpression\" /. CompileDiagnostics[{{x, _Real}}, Sin[x] + BarnesG[x]]",
                   "\"BarnesG[x]\"", 0);
    /* A free symbol is diagnosed as a symbol, not as an unsupported head. */
    assert_eval_eq("\"Subexpression\" /. CompileDiagnostics[{{x, _Real}}, Sin[x] + freevar]",
                   "\"freevar\"", 0);

    /* The optimiser's effect is reportable: a body with a repeated subtree
     * compiles to fewer instructions than the same body with the passes off. */
    assert_eval_eq("With[{d = CompileDiagnostics[{{x, _Real}}, (x + 1)^2 (x + 1)^3 + Sin[x + 1]]}, "
                   "(\"Instructions\" /. d) <= (\"InstructionsUnoptimized\" /. d)]", "True", 0);
    assert_eval_eq("(\"CommonSubexpressions\" /. "
                   "CompileDiagnostics[{{x, _Real}}, (x + 1)^2 (x + 1)^3 + Sin[x + 1]]) > 0",
                   "True", 0);

    /* A malformed argspec leaves the call unevaluated, as Compile does. */
    assert_eval_startswith("CompileDiagnostics[x, x^2]", "CompileDiagnostics[");
}

/* Procedural / control-flow bodies (M2c) reachable through Compile. */
void test_cf_procedural(void) {
    assert_eval_eq("Abs[Compile[{{n, _Integer}}, Module[{s = 0.}, Do[s = s + i, {i, 1, n}]; s]][10] - 55] < 10^-9", "True", 0);
    assert_eval_eq("Compile[{{x, _Real}}, Nest[Function[u, u/2], x, 3]][8.0] == 1", "True", 0);
}

/* The object prints as CompiledFunction[...]; a wrong-arity application stays
 * unevaluated (the object is not consumed). */
void test_cf_object_and_arity(void) {
    assert_eval_startswith("Compile[{{x, _Real}}, x^2]", "CompiledFunction[{x}");
    assert_eval_startswith("Compile[{{x, _Real}}, x][1.0, 2.0]", "CompiledFunction[{x}");
}

/* Array argspec: `{v, _Real, r}` declares a rank-r machine array, the same
 * spelling the Wolfram Language uses.  Until this landed the whole array engine
 * — fusion, any-rank, the reductions — was reachable only from C. */
void test_cf_array_argspec(void) {
    /* A List argument is packed at the boundary and comes back as a List: the
     * result KIND must follow the argument kind, or the compiled path would
     * answer with a different head from the interpreter fallback on the same
     * input.  The two lines below are exactly that pair — the second has a
     * symbolic element, so it falls back — and they must agree in shape. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, v^2 + 2 v + 1][{1., 2., 3.}]",
                   "{4.0, 9.0, 16.0}", 0);
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, 2 v][{1., 2.}]]", "List", 0);
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, 2 v][{1., 2., x}]]", "List", 0);

    /* An NDArray argument stays an NDArray — borrowed in, new one out. */
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, 2 v][NDArray[{1., 2.}]]]", "NDArray", 0);

    /* Any rank. */
    assert_eval_eq("Compile[{{m, _Real, 2}}, m + 1][{{1., 2.}, {3., 4.}}]",
                   "{{2.0, 3.0}, {4.0, 5.0}}", 0);
    assert_eval_eq("Compile[{{t, _Real, 3}}, 2 t][{{{1., 2.}}, {{3., 4.}}}]",
                   "{{{2.0, 4.0}}, {{6.0, 8.0}}}", 0);

    /* Reduction to a scalar, and several array parameters. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Total[v^2]][{1., 2., 3.}] == 14", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}, {u, _Real, 1}}, v u + v][{1., 2.}, {3., 4.}]",
                   "{4.0, 10.0}", 0);

    /* Complex elements. */
    assert_eval_eq("Chop[Total[Compile[{{v, _Complex, 1}}, v v][{1. + 2. I, 3.}]"
                   " - {-3 + 4 I, 9}]] == 0", "True", 0);

    /* Parity with the interpreter over a libm chain — same body, same data. */
    assert_eval_eq("Max[Abs[Compile[{{v, _Real, 1}}, Sin[v] Exp[-v] + Sqrt[v]][{0.5, 1.5, 2.5}]"
                   " - (Sin[#] Exp[-#] + Sqrt[#] & /@ {0.5, 1.5, 2.5})]] < 10^-12", "True", 0);

    /* Shapes that must NOT take the fast path: a ragged list cannot be packed,
     * a rank mismatch is not this function's signature, and a symbolic element
     * is not a machine number.  All three fall back and stay correct. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, 2 v][{{1., 2.}, {3.}}]",
                   "{{2.0, 4.0}, {6.0}}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, 2 v][{{1., 2.}, {3., 4.}}]",
                   "{{2.0, 4.0}, {6.0, 8.0}}", 0);

    /* Malformed argspecs must leave Compile[] unevaluated rather than guess. */
    assert_eval_eq("Head[Compile[{{v, _Real, 0.5}}, v]]", "Compile", 0);
    assert_eval_eq("Head[Compile[{{v, _Integer, 2}}, v]]", "Compile", 0);  /* no integer dtype */
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_cf_integer);
    TEST(test_cf_real);
    TEST(test_cf_complex);
    TEST(test_cf_symbolic_fallback);
    TEST(test_cf_uncompilable_fallback);
    TEST(test_cf_procedural);
    TEST(test_cf_object_and_arity);
    TEST(test_cf_array_argspec);
    TEST(test_cf_part);
    TEST(test_compile_diagnostics);

    printf("All CompiledFunction tests passed!\n");
    return 0;
}
