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
#include "expr.h"
#include "sym_intern.h"
#include "compile/compile.h"
#include "compile/compiled_function.h"

/* Integer-typed compiled functions return exact Integers. */
void test_cf_integer(void) {
    assert_eval_eq("Compile[{{n, _Integer}}, n^2 + 1][5]", "26", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, a b - b][6, 7]", "35", 0);
}

/* Order[a,b] -> the canonical comparison {1,0,-1}, lowered as Sign[b - a] from
 * the existing SUB + SIGN opcodes. Always an Integer result, matching the
 * interpreter's Order (result-HEAD parity), for real, integer and mixed args. */
void test_cf_order(void) {
    /* Real args: 1 when a<b, -1 when a>b, 0 when equal. */
    assert_eval_eq("Compile[{{a, _Real}, {b, _Real}}, Order[a, b]][1., 2.]", "1", 0);
    assert_eval_eq("Compile[{{a, _Real}, {b, _Real}}, Order[a, b]][2., 1.]", "-1", 0);
    assert_eval_eq("Compile[{{a, _Real}, {b, _Real}}, Order[a, b]][1., 1.]", "0", 0);
    /* Integer args. */
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Order[a, b]][3, 7]", "1", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Order[a, b]][7, 3]", "-1", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Order[a, b]][5, 5]", "0", 0);
    /* Mixed Integer/Real args unify to Real. */
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Real}}, Order[a, b]][2, 3.5]", "1", 0);
    /* The body really lowers -- guards against the silent interpreter-bail cliff. */
    assert_eval_eq("Lookup[CompileDiagnostics[{{a, _Real}, {b, _Real}}, Order[a, b]], "
                   "\"Compiled\"]", "True", 0);
    /* Result-HEAD parity: the compiled result is an Integer, not a Real. */
    assert_eval_eq("Head[Compile[{{a, _Real}, {b, _Real}}, Order[a, b]][1., 2.]]",
                   "Integer", 0);
    /* Agreement with the interpreter on numeric args. */
    assert_eval_eq("Compile[{{a, _Real}, {b, _Real}}, Order[a, b]][3., 2.] === Order[3., 2.]",
                   "True", 0);
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

    /* ...and it must STILL be a List when an unrelated NDArray argument is
     * present.  Deciding the result kind from the arguments alone got this
     * wrong: `from_nd` was true because a genuine NDArray came in, so a freshly
     * BUILT array came back as an NDArray where the interpreter's own
     * ConstantArray gives a List.  Every array-constructing head would have
     * inherited that, so the program now reports whether its result was built.
     *
     * The two cases below are the fault line: the first builds its result and
     * ignores the argument's kind, the second DERIVES its result and must keep
     * following it. */
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, ConstantArray[1., 3]][NDArray[{1., 2.}]]]",
                   "List", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, ConstantArray[1., 3]][NDArray[{1., 2.}]]",
                   "{1.0, 1.0, 1.0}", 0);
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, v + 1.][NDArray[{1., 2.}]]]", "NDArray", 0);
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, Module[{u = v}, u + 1.]][NDArray[{1., 2.}]]]",
                   "NDArray", 0);
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, v[[1 ;; 2]]][NDArray[{1., 2., 3.}]]]",
                   "NDArray", 0);
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
    assert_eval_eq("Head[Compile[{{v, _Real, 99}}, v]]", "Compile", 0);   /* rank out of range */

    /* An _Integer array argspec is ACCEPTED (NDT_INT64); it used to be rejected
     * for want of an integer dtype.  A List of Integers goes in and a List of
     * Integers comes back — never an NDArray, because no user syntax can build
     * an int64 one, so handing one back would be a value the rest of the system
     * has no way to have produced. */
    assert_eval_eq("Head[Compile[{{v, _Integer, 1}}, Total[v]][{1, 2, 3}]]", "Integer", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, Total[v]][{1, 2, 3}]", "6", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, v[[2]]][{10, 20, 30}]", "20", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, Length[v]][{1, 2, 3}]", "3", 0);
    assert_eval_eq("Compile[{{m, _Integer, 2}}, m[[2, 1]]][{{1, 2}, {3, 4}}]", "3", 0);
    /* Exact past 2^53, which is where a float64 buffer would start lying and
     * where the double-based ndt_get/ndt_set pair would too. */
    assert_eval_eq("Compile[{{v, _Integer, 1}}, v[[1]]][{9007199254740993}]",
                   "9007199254740993", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, Total[v]][{9007199254740993, 1}]"
                   " === 9007199254740994", "True", 0);
    /* A float64 NDArray is not an integer array: different element types to the
     * interpreter too, so the call goes back rather than rounding. */
    assert_eval_eq("Compile[{{v, _Integer, 1}}, Total[v]][NDArray[{1., 2., 3.}]] == 6",
                   "True", 0);
}

/* The delegated array heads (src/compile/compile.c, ND_FNS and ND_REDS).
 *
 * Each of these has an NDArray entry point in the interpreter and now a
 * lowering that calls that same entry point from the VM, so the compiled
 * answer is the interpreted answer by construction — which is exactly what
 * these assertions pin, element heads included.  The gap they close is not
 * the speed of the operation (it was already a buffer walk) but that a body
 * CONTAINING one no longer bails wholesale: the compilable subset is a cliff,
 * so one unlowerable head used to cost every other head in the body too. */
static void test_cf_delegated_array_heads(void) {
    /* Array -> array, one argument. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Differences[v]][{3., 1., 7.}]",
                   "{-2.0, 6.0}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Ratios[v]][{2., 4., 8.}]",
                   "{2.0, 2.0}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Most[v]][{3., 1., 7.}]", "{3.0, 1.0}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Rest[v]][{3., 1., 7.}]", "{1.0, 7.0}", 0);

    /* Array -> array with a trailing integer. */
    assert_eval_eq("Compile[{{v, _Real, 1}, {k, _Integer}}, RotateLeft[v, k]][{1., 2., 3.}, 1]",
                   "{2.0, 3.0, 1.0}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}, {k, _Integer}}, RotateRight[v, k]][{1., 2., 3.}, 1]",
                   "{3.0, 1.0, 2.0}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}, {k, _Integer}}, MovingAverage[v, k]][{1., 2., 3.}, 2]",
                   "{1.5, 2.5}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}, {k, _Integer}}, TakeLargest[v, k]][{3., 1., 7.}, 2]",
                   "{7.0, 3.0}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}, {k, _Integer}}, TakeSmallest[v, k]][{3., 1., 7.}, 2]",
                   "{1.0, 3.0}", 0);

    /* The delegated reductions.  Bit-identical to the interpreter because both
     * run the same summation, which is the reason for delegating rather than
     * re-implementing the loop in the VM. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Mean[v]][{3., 1., 7., 2., 5.}]"
                   " === Mean[{3., 1., 7., 2., 5.}]", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Median[v]][{3., 1., 7., 2., 5.}]"
                   " === Median[{3., 1., 7., 2., 5.}]", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Variance[v]][{3., 1., 7., 2., 5.}]"
                   " === Variance[{3., 1., 7., 2., 5.}]", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, StandardDeviation[v]][{3., 1., 7., 2.}]"
                   " === StandardDeviation[{3., 1., 7., 2.}]", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, RootMeanSquare[v]][{3., 1., 7., 2.}]"
                   " === RootMeanSquare[{3., 1., 7., 2.}]", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Max[v]][{3., 1., 7., 2.}]", "7.0", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Min[v]][{3., 1., 7., 2.}]", "1.0", 0);

    /* THE REGRESSION.  Max/Min used to reach the scalar pairwise fold with an
     * array operand, and an empty fold returns its accumulator — so Max[v]
     * lowered to the IDENTITY and Max[v] + 1. answered {4., 2., 8., 3.} where
     * the interpreter answers 8.  Alone, Max[v] happened to be rejected further
     * down, which is why this needed the head to sit inside a larger expression
     * to show.  Found by tools/compile_coverage.py, 2026-08-02. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Max[v] + 1.][{3., 1., 7., 2.}]", "8.0", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Min[v] * 2.][{3., 1., 7., 2.}]", "2.0", 0);
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, Max[v] + 1.][{3., 1., 7., 2.}]]",
                   "Real", 0);
    /* The scalar fold itself is untouched. */
    assert_eval_eq("Compile[{{x, _Real}, {y, _Real}}, Max[x, y, 3.]][1., 9.]", "9.0", 0);
    assert_eval_eq("Compile[{{x, _Integer}, {y, _Integer}}, Min[x, y, 3]][1, 9]", "1", 0);

    /* Max/Min SELECT an element, so an integer vector answers with an Integer;
     * Mean and friends AVERAGE, so an integer vector answers with a Rational
     * that no machine slot holds and the body must decline rather than round. */
    assert_eval_eq("Compile[{{v, _Integer, 1}}, Max[v]][{3, 1, 7}]", "7", 0);
    assert_eval_eq("Head[Compile[{{v, _Integer, 1}}, Min[v]][{3, 1, 7}]]", "Integer", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, Mean[v]][{1, 2}]", "3/2", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, Ratios[v]][{3, 1, 7}]",
                   "{1/3, 7}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, Differences[v]][{3, 1, 7}]",
                   "{-2, 6}", 0);

    /* Rank 2 is a DIFFERENT operation for the reductions — Mean of a matrix is
     * the vector of column means — so the rank-1-only rule must hold. */
    assert_eval_eq("Compile[{{m, _Real, 2}}, Mean[m]][{{1., 2.}, {3., 4.}}]",
                   "{2.0, 3.0}", 0);

    /* Degenerate shapes the ND path declines: the interpreter answers and the
     * two still agree. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Most[v]][{4.}]", "{}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Differences[v]][{4.}]", "{}", 0);

    /* And the point of the exercise: a body that MIXES a delegated head with
     * ordinary arithmetic compiles as one program rather than falling off the
     * cliff at the first unlowered head. */
    assert_eval_eq("Lookup[CompileDiagnostics[{{v, _Real, 1}}, "
                   "(Max[v] - Min[v]) / Mean[v]], \"Compiled\"]", "True", 0);
    assert_eval_eq("Lookup[CompileDiagnostics[{{v, _Real, 1}}, "
                   "Total[Differences[Sort[v]]]], \"Compiled\"]", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, (Max[v] - Min[v]) / Mean[v]][{3., 1., 7., 2., 5.}]"
                   " === (Max[#] - Min[#]) / Mean[#] &[{3., 1., 7., 2., 5.}]", "True", 0);
}

/* The NARROWING kernels over an array — real in, exact Integer out.
 *
 * Two separate defects, both found by tools/compile_coverage.py on 2026-08-02
 * and both invisible from the head on its own.
 *
 * 1. The declared element type LIED.  emit_arr_unary computed the result
 *    element as "to_real ? Real : the input's", which for Floor over a Real
 *    array says Real — while ndarray_map_unary, preferring the to_int path,
 *    writes an NDT_INT64 buffer.  Standalone the caller re-reads the real
 *    dtype and all is well; a CONSUMER inside the same program reads the slot
 *    as the declared type and gets the integer's bits back as a double.
 *    Total[Floor[v]] answered 2.96439*10^-323, which is the int64 6.
 * 2. IntegerPart and UnitStep had no array lowering at all: both have a
 *    dedicated scalar branch (so try_kernel is never reached) and a kernel
 *    with neither a real nor a complex arm (so the array interception list
 *    skipped them).  Each took its whole enclosing body down to the
 *    interpreter. */
/* Ordering[vector] lowers through the same ND_FNS delegation as Sort, with the
 * one wrinkle that its result element type is int64 whatever the operand's dtype
 * (a permutation is integer) -- the NdFnSpec.int_result path. */
static void test_cf_ordering(void) {
    assert_eval_eq("Compile[{{v, _Real, 1}}, Ordering[v]][{3., 1., 2.}]", "{2, 3, 1}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, Ordering[v]][{3, 1, 2}]", "{2, 3, 1}", 0);
    /* The result element is Integer even for a _Real input. */
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, Ordering[v]][{3., 1., 2.}][[1]]]",
                   "Integer", 0);
    /* Bit-identical to the interpreter, ties (stability) included. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Ordering[v]][{2., 6., 1., 9., 1., 2., 3.}]"
                   " === Ordering[{2., 6., 1., 9., 1., 2., 3.}]", "True", 0);
    /* CompileDiagnostics confirms it lowers at a rank-1 array shape. */
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Real, 1}}, Ordering[v]]",
                   "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, Ordering[v]]",
                   "True", 0);
    /* As a subexpression of a scalar body: Ordering[v][[1]] is the position of the
     * minimum (correct whether or not this exact shape compiles). */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Ordering[v][[1]]][{3., 1., 2.}]", "2", 0);
}

static void test_cf_narrowing_kernels_over_arrays(void) {
    /* The reinterpretation bug: a narrowing head inside a reduction. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Total[Floor[v]]][{1.5, 2.5, 3.5}]", "6", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Total[Sign[v]]][{1.5, -2.5, 3.5}]", "1", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Total[Ceiling[v]]][{1.5, 2.5}]", "5", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Total[Round[v]]][{1.4, 2.6}]", "4", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Total[IntegerPart[v]]][{1.9, -2.9}]", "-1", 0);

    /* ... and the element HEADS it made wrong: exact Integers, not Reals. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Floor[v] + 1][{1.5, 2.5, 3.5}]",
                   "{2, 3, 4}", 0);
    assert_eval_eq("Head[First[Compile[{{v, _Real, 1}}, Floor[v] + 1][{1.5, 2.5}]]]",
                   "Integer", 0);

    /* The two heads with no array lowering at all. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, UnitStep[v]][{-1.5, 0., 2.5}]",
                   "{0, 1, 1}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, IntegerPart[v]][{-1.5, 0., 2.5}]",
                   "{-1, 0, 2}", 0);
    assert_eval_eq("Lookup[CompileDiagnostics[{{v, _Real, 1}}, UnitStep[v]], \"Compiled\"]",
                   "True", 0);
    assert_eval_eq("Lookup[CompileDiagnostics[{{v, _Real, 1}}, IntegerPart[v]], \"Compiled\"]",
                   "True", 0);
    /* The scalar and multi-argument UnitStep forms are untouched. */
    assert_eval_eq("Compile[{{x, _Real}}, UnitStep[x]][-0.5]", "0", 0);
    assert_eval_eq("Compile[{{x, _Real}, {y, _Real}}, UnitStep[x, y]][1., -1.]", "0", 0);

    /* Abs is `to_int` too, but with only the int64 arm — a real array is a
     * PROJECTION and must stay Real, so the same rule has to give a different
     * answer for it than for Floor. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Total[Abs[v]]][{1.5, -2.5}]", "4.0", 0);
    assert_eval_eq("Head[First[Compile[{{v, _Real, 1}}, Abs[v]][{-1.5, 2.5}]]]", "Real", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, Abs[v]][{-3, 4}]", "{3, 4}", 0);

    /* Every one of them still agrees with the interpreter element for element. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Floor[v]][{-1.5, 0., 2.5}]"
                   " === Floor[{-1.5, 0., 2.5}]", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Sign[v]][{-1.5, 0., 2.5}]"
                   " === Sign[{-1.5, 0., 2.5}]", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, UnitStep[v]][{-1.5, 0., 2.5}]"
                   " === UnitStep[{-1.5, 0., 2.5}]", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, IntegerPart[v]][{-1.5, 0., 2.5}]"
                   " === IntegerPart[{-1.5, 0., 2.5}]", "True", 0);
}

/* The exact-integer BINARY kernels over an array — Mod, Quotient, GCD, LCM and
 * the two-argument ArcTan.
 *
 * These have no complex arm at all (Mod and Quotient are narrowing, GCD and LCM
 * are defined only on Z), and the compiler read a missing complex arm as the
 * degrade sentinel — so the array spelling had no lowering while the SCALAR one
 * had a dedicated opcode.  `Mod[x, 3]` compiled and `Mod[v, 3]` did not, and
 * one unlowered head costs the whole body. */
static void test_cf_integer_binary_kernels_over_arrays(void) {
    assert_eval_eq("Compile[{{v, _Real, 1}, {k, _Integer}}, Mod[v, k]][{7.5, -3.5, 10.}, 3]"
                   " === Mod[{7.5, -3.5, 10.}, 3]", "True", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}, {k, _Integer}}, Mod[v, k]][{7, -3, 10}, 3]",
                   "{1, 0, 1}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}, {k, _Integer}}, Quotient[v, k]][{7.5, -3.5, 10.}, 3]",
                   "{2, -2, 3}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}, {k, _Integer}}, Quotient[v, k]][{7, -3, 10}, 3]",
                   "{2, -1, 3}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}, {k, _Integer}}, GCD[v, k]][{12, 18, 25}, 6]",
                   "{6, 6, 1}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}, {k, _Integer}}, LCM[v, k]][{4, 6, 9}, 6]",
                   "{12, 6, 18}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}, {y, _Real}}, ArcTan[v, y]][{1., -1.}, 2.]"
                   " === ArcTan[{1., -1.}, 2.]", "True", 0);

    /* Element heads: a real Mod is Real, a real Quotient is an exact Integer. */
    assert_eval_eq("Head[First[Compile[{{v, _Real, 1}, {k, _Integer}}, Mod[v, k]]"
                   "[{7.5}, 3]]]", "Real", 0);
    assert_eval_eq("Head[First[Compile[{{v, _Real, 1}, {k, _Integer}}, Quotient[v, k]]"
                   "[{7.5}, 3]]]", "Integer", 0);

    /* Nested, which is the point: the head no longer takes its body down. */
    assert_eval_eq("Lookup[CompileDiagnostics[{{v, _Real, 1}, {k, _Integer}}, "
                   "Total[Mod[v, k]]], \"Compiled\"]", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}, {k, _Integer}}, Total[Mod[v, k]]]"
                   "[{7.5, -3.5, 10.}, 3] === Total[Mod[{7.5, -3.5, 10.}, 3]]", "True", 0);

    /* Abandon, never wrap: a zero divisor has no exact int64 answer, so the
     * kernel declines the WHOLE array and the interpreter answers — with the
     * unevaluated Mod[k, 0] the interpreter itself gives. */
    assert_eval_eq("Compile[{{v, _Integer, 1}, {k, _Integer}}, Mod[v, k]][{7, -3}, 0]"
                   " === Mod[{7, -3}, 0]", "True", 0);
    /* And exactness past 2^53, where a double round-trip would lie. */
    assert_eval_eq("Compile[{{v, _Integer, 1}, {k, _Integer}}, Quotient[v, k]]"
                   "[{9007199254740993}, 2]", "{4503599627370496}", 0);

    /* The scalar forms keep their dedicated opcodes. */
    assert_eval_eq("Compile[{{x, _Real}, {k, _Integer}}, Mod[x, k]][7.5, 3]", "1.5", 0);
    assert_eval_eq("Compile[{{x, _Integer}, {k, _Integer}}, Mod[x, k]][7, 3]", "1", 0);

    /* The INTEGER-ONLY unary kernels (MoebiusMu, EulerPhi, IntegerLength) have
     * `to_int_i` and nothing else -- MoebiusMu of a real is not a machine
     * question -- so a guard keyed on the REAL narrowing arm read exactly the
     * heads it was written for as sentinels.  Integer arrays only; a real one
     * still declines. */
    assert_eval_eq("Compile[{{v, _Integer, 1}}, MoebiusMu[v]][{6, 12, 30}]",
                   "{1, 0, -1}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, EulerPhi[v]][{6, 12, 30}]",
                   "{2, 4, 8}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, IntegerLength[v]][{100, 255, 7}]",
                   "{3, 3, 1}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, Total[EulerPhi[v]]][{6, 12, 30}]",
                   "14", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, EulerPhi[v]][{6., 12.}]"
                   " === EulerPhi[{6., 12.}]", "True", 0);

    /* MovingMedian joins MovingAverage in the delegated table. */
    assert_eval_eq("Compile[{{v, _Real, 1}, {k, _Integer}}, MovingMedian[v, k]]"
                   "[{3., 1., 7., 2., 5.}, 3] === MovingMedian[{3., 1., 7., 2., 5.}, 3]",
                   "True", 0);
}

/* Integer ARRAYS, built rather than passed in.  Each construct here used to bail
 * for want of an integer dtype, so each is checked against the interpreter
 * evaluating the same expression — element heads included, since that is the
 * property the restriction existed to protect. */
void test_cf_integer_arrays(void) {
    assert_eval_eq("Compile[{{n, _Integer}}, Range[n]][5] === Range[5]", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Range[n]][0]", "{}", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Range[n]][-3]", "{}", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Range[a, b]][3, 7] === Range[3, 7]",
                   "True", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Range[a, b]][7, 3]", "{}", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Range[a, b, 2]][1, 9]"
                   " === Range[1, 9, 2]", "True", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Range[a, b, -2]][9, 1]"
                   " === Range[9, 1, -2]", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Total[Range[n]]][100]", "5050", 0);
    /* A REAL Range must still bail: the interpreter walks a real iterator by
     * repeated addition against a slack, which a closed-form step does not
     * reproduce — in the last bits, and at the endpoint in the COUNT. */
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, Range[x]]", "False", 0);

    assert_eval_eq("Compile[{{n, _Integer}}, Table[i, {i, 1, n}]][5] === Table[i, {i, 1, 5}]",
                   "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Table[i^2, {i, 1, n}]][4]"
                   " === Table[i^2, {i, 1, 4}]", "True", 0);
    assert_eval_eq("Head[First[Compile[{{n, _Integer}}, Table[i, {i, 1, n}]][3]]]", "Integer", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, ConstantArray[7, n]][4] === ConstantArray[7, 4]",
                   "True", 0);
    assert_eval_eq("Head[First[Compile[{{n, _Integer}}, ConstantArray[0, n]][3]]]", "Integer", 0);

    /* Build-and-fill through an integer local: A_STORE_I on one side and
     * A_LOAD_I on the other. */
    assert_eval_eq("Compile[{{n, _Integer}}, Module[{u = ConstantArray[0, n]}, "
                   "Do[u[[i]] = i*i, {i, 1, n}]; u]][5] === Table[i^2, {i, 1, 5}]", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Module[{u = Range[n], s = 0}, "
                   "Do[s = s + u[[i]], {i, 1, n}]; s]][10]", "55", 0);

    /* The overflow rule reaches into array ELEMENTS too: a few elements here
     * leave the int64 range, so the whole call defers and the interpreter
     * produces the exact bigints.  Deliberately a SHORT range with large
     * elements — a long range would make the fallback build a huge list. */
    assert_eval_eq("Compile[{{n, _Integer}}, Table[(i*4000000000)*4000000000, {i, 1, n}]][3]"
                   " === Table[(i*4000000000)*4000000000, {i, 1, 3}]", "True", 0);
}

/* ------------------------------------------------------------------ *
 *  CompilePrint / disassembler                                         *
 * ------------------------------------------------------------------ *
 * Asserted through the C entry point rather than through the builtin,
 * because CompilePrint[] writes to stdout and returns Null — the return
 * value proves nothing about the listing.  The listing is deterministic
 * by construction (kernels are resolved to names, callee programs and
 * parallel loops are numbered, no address is ever printed), which is
 * exactly what makes these substring assertions stable. */

static char* disasm_of(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* v = evaluate(parsed);
    expr_free(parsed);
    ASSERT(v != NULL && v->type == EXPR_COMPILED);
    char* s = compiled_function_disassemble(v->data.compiled);
    expr_free(v);
    ASSERT(s != NULL);
    return s;
}

static void has(const char* text, const char* needle, const char* what) {
    if (!strstr(text, needle)) {
        fprintf(stderr, "FAIL: %s — expected \"%s\" in disassembly:\n%s\n", what, needle, text);
        exit(1);
    }
}

static void lacks(const char* text, const char* needle, const char* what) {
    if (strstr(text, needle)) {
        fprintf(stderr, "FAIL: %s — unexpected \"%s\" in disassembly:\n%s\n", what, needle, text);
        exit(1);
    }
}

/* Header, scalar opcodes, and the constants folded INTO the instructions:
 * `_RK` is the whole point of looking at a polynomial body. */
void test_disasm_scalar(void) {
    char* d = disasm_of("Compile[{{x, _Real}}, x^2 + 2.5 x + 1]");
    has(d, "R0   : Real", "argument register and type");
    has(d, "Result", "result line");
    has(d, "all-Real fast path", "all-Real signature reported");
    has(d, "POWI_R", "integer power opcode");
    has(d, "MUL_RK", "constant folded into the multiply");
    has(d, "ADD_RK", "constant folded into the add");
    has(d, "R0 * 2.5", "rendered constant operand");
    has(d, "return R", "RET rendering");
    /* Every CONST that survives must have had its type recovered; raw bits
     * mean the inference gave up. */
    lacks(d, "0x", "no un-typed constants");
    free(d);
}

/* A special function lowers to a machine kernel whose immediate is a bare
 * function pointer; the listing must name it, not print an address. */
void test_disasm_kernel_names(void) {
    char* d = disasm_of("Compile[{{x, _Real}}, Gamma[x] + Erf[x]]");
    has(d, "<Gamma>", "unary kernel resolved to its symbol");
    has(d, "Gamma[R", "kernel rendered as a call");
    free(d);

    char* n = disasm_of("Compile[{{a, _Real}, {b, _Real}, {c, _Real}, {z, _Real}}, "
                        "Hypergeometric2F1[a, b, c, z]]");
    has(n, "KERNN", "n-ary kernel opcode");
    has(n, "<Hypergeometric2F1>", "n-ary kernel resolved");
    has(n, "[R", "operand range printed");
    free(n);
}

/* Control flow: branch targets are marked and rendered as jumps. */
void test_disasm_control_flow(void) {
    char* d = disasm_of("Compile[{{x, _Real}}, If[x > 0, Sqrt[x], -Sqrt[-x]]]");
    has(d, "JZ", "conditional branch");
    has(d, "goto", "rendered jump");
    has(d, "\n>", "branch target marked in the gutter");
    free(d);

    /* A unit-step Do closes with OP_LOOP, which increments, tests and branches
     * in ONE instruction — the whole point of it, since the general shape spends
     * four instructions per iteration on control alone.  Asserting the opcode
     * rather than "some increment happens" is what makes a regression to the
     * four-instruction form visible here rather than only in a benchmark. */
    char* l = disasm_of("Compile[{{n, _Integer}}, Module[{s = 0.}, "
                        "Do[s = s + 1./k, {k, 1, n}]; s]]");
    has(l, "LOOP", "counted loop closes with the fused increment/test/branch");
    has(l, "if ++", "OP_LOOP rendered as its increment-and-test");
    has(l, "R2 = 0.", "Real constant typed through the MOVE that initialises it");
    free(l);

    /* A non-unit step cannot use it (OP_LOOP steps by one), so that form keeps
     * the explicit increment — and both must still be correct, which the
     * engine's own iterator tests cover. */
    char* l2 = disasm_of("Compile[{{n, _Integer}}, Module[{s = 0.}, "
                         "Do[s = s + 1./k, {k, 1, n, 2}]; s]]");
    has(l2, "INC_I", "explicit increment on a non-unit step");
    free(l2);
}

/* Arrays: the three register banks, the strip-mined loop and the fan-out
 * marker are the facts a fused body is read for. */
void test_disasm_arrays(void) {
    char* d = disasm_of("Compile[{{v, _Real, 1}}, v^2 + 2 v + 1]");
    has(d, "V0   : Real[1]", "array argument keeps element type and rank");
    has(d, "APAR", "parallel fan-out marker emitted");
    has(d, "parallel loop", "fan-out reported in the header");
    has(d, "VLOAD_R", "strip-mined load");
    has(d, "VSTORE_R", "strip-mined store");
    has(d, "T", "tile registers named");
    /* argdep must count a fused leaf as a read of the argument. */
    lacks(d, "(unused)", "fused array argument is not reported unused");
    free(d);

    char* p = disasm_of("Compile[{{v, _Real, 1}}, v[[2 ;; 3]]]");
    has(p, "A_PART", "general Part opcode");
    has(p, "Span[2, 3]", "literal subscript rendered from the PartSpec");
    free(p);

    /* A delegated head must say WHICH head.  Before the name accessors an
     * A_NDFN rendered as bare "A_NDFN" and a V_NDRED as "V_NDRED(V0, V0)" —
     * telling Mean from Median is the point of having a dump. */
    char* n = disasm_of("Compile[{{v, _Real, 1}}, Total[Differences[v]]]");
    has(n, "A_NDFN", "delegated structural opcode");
    has(n, "Differences[", "A_NDFN names the head it delegates to");
    free(n);

    char* r = disasm_of("Compile[{{v, _Real, 1}}, (Max[v] - Min[v])/Mean[v]]");
    has(r, "V_NDRED", "delegated reduction opcode");
    has(r, "Max[", "V_NDRED names its head");
    has(r, "Mean[", "and distinguishes the three reductions in one body");
    has(r, "Min[", "and distinguishes the three reductions in one body");
    free(r);

    char* t = disasm_of("Compile[{{v, _Real, 1}, {k, _Integer}}, "
                        "Total[RotateLeft[v, k]]]");
    has(t, "RotateLeft[", "A_NDFN renders its trailing integer operands");
    free(t);
}

/* A body outside the compilable subset has no bytecode, so the useful
 * answer is why — the same question CompileDiagnostics answers about a
 * body, asked about an object. */
void test_disasm_uncompiled(void) {
    char* d = disasm_of("Compile[{x}, Integrate[x, x]]");
    has(d, "not compiled", "reports that there is no program");
    has(d, "Reason", "bail reason reported");
    has(d, "Integrate[x, x]", "offending subexpression reported");
    free(d);
}

/* OP_CALL's immediate is a borrowed callee program, listed once in its own
 * numbered section.  A user Compile[] never emits one (that needs
 * COMPILE_FOLD_GLOBALS, which only the autocompile paths set), so drive the
 * engine directly rather than leave the worklist untested. */
void test_disasm_callee_program(void) {
    /* Big enough that the inliner declines and a real CALL is emitted. */
    Expr* def = parse_expression(
        "gBig = Compile[{{x, _Real}}, Sin[x] + Cos[x] + Sin[2 x] + Cos[2 x] + "
        "Sin[3 x] + Cos[3 x] + Sin[4 x] + Cos[4 x] + Sin[5 x] + Cos[5 x] + "
        "Sin[6 x] + Cos[6 x] + Sin[7 x] + Cos[7 x]]");
    ASSERT(def != NULL);
    expr_free(evaluate(def));
    expr_free(def);

    Expr* body = parse_expression("gBig[y] + 1.");
    ASSERT(body != NULL);
    const char* names[1]; names[0] = intern_symbol("y");
    CompileType types[1]; types[0] = CT_REAL;
    CompiledProgram* p = compile_expr_ex(body, names, types, 1, COMPILE_FOLD_GLOBALS);
    expr_free(body);
    ASSERT(p != NULL);

    char* d = compiled_disassemble(p, names);
    ASSERT(d != NULL);
    has(d, "CALL", "call opcode");
    has(d, "<call #1>", "callee numbered, not printed as an address");
    has(d, "--- called program #1 ---", "callee listed in its own section");
    free(d);
    compiled_free(p);
}

/* ------------------------------------------------------------------ *
 *  RuntimeAttributes -> Listable                                       *
 * ------------------------------------------------------------------ */

/* Threading over scalar parameters.
 *
 * The discriminating body is `If[x > 0, ...]`, NOT something like `x^2`: with a
 * List substituted for x, Power/Plus are themselves Listable, so an arithmetic
 * body gives the threaded answer through the interpreter fallback whether the
 * object is Listable or not — a test built on one would pass with the feature
 * removed.  `Greater` does not thread, so the non-Listable object provably
 * cannot produce a List here. */
void test_cf_runtime_attributes(void) {
    assert_eval_eq("Compile[{{x, _Real}}, If[x > 0, 1., -1.], RuntimeAttributes -> Listable]"
                   "[{1., -2., 3.}]", "{1.0, -1.0, 1.0}", 0);
    /* The control: same body, default RuntimeAttributes -> {}, no threading. */
    assert_eval_eq("Head[Compile[{{x, _Real}}, If[x > 0, 1., -1.]][{1., -2., 3.}]]", "If", 0);
    assert_eval_eq("Head[Compile[{{x, _Real}}, If[x > 0, 1., -1.], "
                   "RuntimeAttributes -> {}][{1., -2., 3.}]]", "If", 0);

    /* Both spellings of the setting; scalar arguments are unaffected. */
    assert_eval_eq("Compile[{{x, _Real}}, x^2, RuntimeAttributes -> {Listable}][{1., 2., 3.}]",
                   "{1.0, 4.0, 9.0}", 0);
    assert_eval_eq("Compile[{{x, _Real}}, x^2, RuntimeAttributes -> Listable][3.] == 9",
                   "True", 0);

    /* Every result type survives threading: Integer, Complex, Boolean. */
    assert_eval_eq("Compile[{{n, _Integer}}, n^2 + 1, RuntimeAttributes -> Listable][{1, 2, 3}]",
                   "{2, 5, 10}", 0);
    assert_eval_eq("Chop[Compile[{{z, _Complex}}, z^2, RuntimeAttributes -> Listable]"
                   "[{1. + 2. I, 3.}] - {-3 + 4 I, 9}] == {0, 0}", "True", 0);
    assert_eval_eq("Compile[{{x, _Real}}, x > 0, RuntimeAttributes -> Listable][{1., -1.}]",
                   "{True, False}", 0);

    /* Nested lists thread level by level; an empty list threads to an empty one. */
    assert_eval_eq("Compile[{{x, _Real}}, x^2, RuntimeAttributes -> Listable]"
                   "[{{1., 2.}, {3., 4.}}]", "{{1.0, 4.0}, {9.0, 16.0}}", 0);
    assert_eval_eq("Compile[{{x, _Real}}, x^2, RuntimeAttributes -> Listable][{}]", "{}", 0);

    /* Several parameters: equal-length lists thread in step, a scalar is reused
     * for every element, and unequal lengths leave the application unevaluated
     * (after a Thread::tdlen message), exactly as for a Listable symbol. */
    assert_eval_eq("Compile[{{x, _Real}, {y, _Real}}, x - y, RuntimeAttributes -> Listable]"
                   "[{1., 2.}, {10., 20.}]", "{-9.0, -18.0}", 0);
    assert_eval_eq("Compile[{{x, _Real}, {y, _Real}}, x - y, RuntimeAttributes -> Listable]"
                   "[{1., 2.}, 10.]", "{-9.0, -8.0}", 0);
    assert_eval_startswith("Compile[{{x, _Real}, {y, _Real}}, x - y, "
                           "RuntimeAttributes -> Listable][{1., 2.}, {10., 20., 30.}]",
                           "CompiledFunction[{x, y}");

    /* Listable belongs to the OBJECT, not to its bytecode: a body outside the
     * compilable subset threads too, and each element still falls back on its
     * own (the second element here is symbolic). */
    assert_eval_eq("uncra[t_] := t; Compile[{{x, _Real}}, uncra[x^2], "
                   "RuntimeAttributes -> Listable][{2., 3.}]", "{4.0, 9.0}", 0);
    assert_eval_eq("Compile[{{x, _Real}}, x^2, RuntimeAttributes -> Listable][{a, 2.}]",
                   "{a^2, 4.0}", 0);

    /* An unusable option leaves Compile[] unevaluated rather than silently
     * building an object that ignores what was asked for. */
    assert_eval_eq("Head[Compile[{{x, _Real}}, x, RuntimeAttributes -> Orderless]]",
                   "Compile", 0);
    assert_eval_eq("Head[Compile[{{x, _Real}}, x, RuntimeAttributes -> {Listable, Flat}]]",
                   "Compile", 0);
    assert_eval_eq("Head[Compile[{{x, _Real}}, x, ThisIsNotAnOption -> True]]", "Compile", 0);

    /* The registered defaults, and SetOptions changing them. */
    assert_eval_eq("Options[Compile]",
                   "{RuntimeAttributes -> {}, "
                   "RuntimeOptions -> {\"CatchMachineIntegerOverflow\" -> True}}", 0);
    assert_eval_eq("SetOptions[Compile, RuntimeAttributes -> Listable]; "
                   "Compile[{{x, _Real}}, If[x > 0, 1., -1.]][{1., -2.}]", "{1.0, -1.0}", 0);
    /* Restore, and prove the restore took: leaving this on would silently make
     * every later test in this file run against a Listable Compile. */
    assert_eval_eq("SetOptions[Compile, RuntimeAttributes -> {}]; "
                   "Head[Compile[{{x, _Real}}, If[x > 0, 1., -1.]][{1., -2.}]]", "If", 0);

    /* CompilePrint reports it — a call that threads before running any bytecode
     * should not be invisible in the listing. */
    char* d = disasm_of("Compile[{{x, _Real}}, x^2, RuntimeAttributes -> Listable]");
    has(d, "Attributes  Listable", "RuntimeAttributes reported in the listing");
    free(d);
    char* p = disasm_of("Compile[{{x, _Real}}, x^2]");
    lacks(p, "Attributes", "no attribute line for a plain object");
    free(p);
}

/* Machine-integer overflow.
 *
 * The engine's contract is that a compiled body answers IDENTICALLY to the
 * interpreter or does not answer at all.  int64 arithmetic breaks that on its
 * own — the interpreter promotes to a bigint where the machine wraps — so every
 * integer opcode that can overflow abandons the call and lets the interpreter
 * redo it.  Each case below is written as a comparison against the interpreter
 * evaluating the SAME expression, so the test states the contract rather than a
 * particular number. */
void test_cf_integer_overflow(void) {
    /* One case per checked opcode, at or just past the int64 boundary. */
    assert_eval_eq("Compile[{{n, _Integer}}, n^3][3000000] === 3000000^3", "True", 0);        /* POWI_I */
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, a*b][4000000000, 4000000000]"
                   " === 4000000000*4000000000", "True", 0);                                  /* MUL_I  */
    assert_eval_eq("Compile[{{a, _Integer}}, a + 1][9223372036854775807]"
                   " === 9223372036854775807 + 1", "True", 0);                                /* ADD_IK */
    assert_eval_eq("Compile[{{a, _Integer}}, a - 1][-9223372036854775808]"
                   " === -9223372036854775808 - 1", "True", 0);                               /* SUB_IK */
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, a - b][-9223372036854775808, 1]"
                   " === -9223372036854775808 - 1", "True", 0);                               /* SUB_I  */
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, a + b][9223372036854775807, 1]"
                   " === 9223372036854775807 + 1", "True", 0);                                /* ADD_I  */
    assert_eval_eq("Compile[{{a, _Integer}}, -a][-9223372036854775808]"
                   " === -(-9223372036854775808)", "True", 0);                                /* NEG_I  */
    assert_eval_eq("Compile[{{a, _Integer}}, Abs[a]][-9223372036854775808]"
                   " === Abs[-9223372036854775808]", "True", 0);                              /* ABS_I  */
    assert_eval_eq("Compile[{{a, _Integer}}, a*7][2000000000000000000]"
                   " === 2000000000000000000*7", "True", 0);                                  /* MUL_IK */

    /* An accumulator inside a loop is the same opcodes, so it inherits the
     * property — worth pinning because the overflow happens mid-run, after the
     * loop has already written to the accumulator many times. */
    assert_eval_eq("Compile[{{n, _Integer}}, Sum[i^2, {i, 1, n}]][3000000]"
                   " === Sum[i^2, {i, 1, 3000000}]", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Module[{s = 1}, Do[s = s*3, {i, 1, n}]; s]][50]"
                   " === 3^50", "True", 0);

    /* Values that FIT must not be disturbed by any of this. */
    assert_eval_eq("Compile[{{n, _Integer}}, n^2 + 3 n + 1][10]", "131", 0);
    assert_eval_eq("Compile[{{a, _Integer}}, Abs[a]][-5]", "5", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, n^3][2000000] === 2000000^3", "True", 0);

    /* Integer division traps the hardware rather than merely answering wrongly:
     * a zero divisor and INT64_MIN/-1 both used to take the process down with
     * SIGFPE.  The interpreter leaves Mod[5, 0] unevaluated, so the compiled
     * path must hand it back and let it do that. */
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Mod[a, b]][5, 0]", "Mod[5, 0]", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Quotient[a, b]][5, 0]",
                   "Quotient[5, 0]", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Quotient[a, b]]"
                   "[-9223372036854775808, -1] === Quotient[-9223372036854775808, -1]",
                   "True", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Mod[a, b]][17, 5]", "2", 0);
}

/* Integer-CLOSED heads: given integer arguments the interpreter returns an
 * Integer, so the compiled object must too.
 *
 * These all have a registered REAL kernel, and reaching it first is how they
 * came back as `35.` and `120.` — values that compare equal to the interpreter's
 * but carry a different head, which IntegerQ, exact arithmetic and the printer
 * all notice.  So each case asserts Head as well as value, and the
 * out-of-domain cases assert that the interpreter's own answer comes back
 * (ComplexInfinity, a Rational, a bigint) rather than a machine approximation
 * of it. */
void test_cf_integer_closed_heads(void) {
    /* Power: Integer for a non-negative exponent, Rational below, Indeterminate
     * at 0^0, exact bigint past the machine range. */
    assert_eval_eq("Head[Compile[{{a, _Integer}, {b, _Integer}}, a^b][2, 10]]", "Integer", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, a^b][2, 10] === 2^10", "True", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, a^b][7, -3]", "1/343", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, a^b][0, 0]", "Indeterminate", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, a^b][2, 64] === 2^64", "True", 0);

    /* Factorial / Gamma: ComplexInfinity outside the domain, bigint past 20!. */
    assert_eval_eq("Head[Compile[{{n, _Integer}}, Factorial[n]][5]]", "Integer", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Factorial[n]][20] === 20!", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Factorial[n]][21] === 21!", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Factorial[n]][-1]", "ComplexInfinity", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Factorial[n]][0]", "1", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Gamma[n]][3]", "2", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Gamma[n]][0]", "ComplexInfinity", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Gamma[n]][-2]", "ComplexInfinity", 0);

    /* Binomial: the multiplicative recurrence must stay exact well past the
     * point where n!/(k!(n-k)!) would overflow — C(60,30) fits an int64 while
     * 60! does not, so a naive lowering fails this one and not the small ones. */
    assert_eval_eq("Head[Compile[{{n, _Integer}, {k, _Integer}}, Binomial[n, k]][7, 3]]",
                   "Integer", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, Binomial[n, k]][60, 30]"
                   " === Binomial[60, 30]", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, Binomial[n, k]][70, 35]"
                   " === Binomial[70, 35]", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, Binomial[n, k]][3, 7]", "0", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {k, _Integer}}, Binomial[n, k]][7, -1]", "0", 0);

    /* Pochhammer: Rational for a negative second argument. */
    assert_eval_eq("Compile[{{a, _Integer}, {n, _Integer}}, Pochhammer[a, n]][7, 3]", "504", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {n, _Integer}}, Pochhammer[a, n]][7, 0]", "1", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {n, _Integer}}, Pochhammer[a, n]][7, -3]",
                   "1/120", 0);

    /* Fibonacci / LucasL, including the negative index, where the two have
     * OPPOSITE sign rules — F[-10] = -55 but L[-10] = +123. */
    assert_eval_eq("Compile[{{n, _Integer}}, Fibonacci[n]][10]", "55", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Fibonacci[n]][-10]", "-55", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Fibonacci[n]][0]", "0", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Fibonacci[n]][92] === Fibonacci[92]", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Fibonacci[n]][93] === Fibonacci[93]", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, LucasL[n]][10]", "123", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, LucasL[n]][-10]", "123", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, LucasL[n]][-1]", "-1", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, LucasL[n]][0]", "2", 0);

    /* Re / Im / Arg / FractionalPart over the integers. Arg is 0 at or above
     * zero and the SYMBOL Pi below, which no machine type holds. */
    assert_eval_eq("Head[Compile[{{n, _Integer}}, Im[n]][3]]", "Integer", 0);
    assert_eval_eq("Head[Compile[{{n, _Integer}}, Re[n]][3]]", "Integer", 0);
    assert_eval_eq("Head[Compile[{{n, _Integer}}, FractionalPart[n]][3]]", "Integer", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Arg[n]][3]", "0", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Arg[n]][0]", "0", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, Arg[n]][-3]", "Pi", 0);

    /* The REAL lowerings of the same heads must be untouched: these branches sit
     * in front of the kernel dispatch, so a mistake there silently removes the
     * real fast path rather than producing a wrong answer. */
    assert_eval_eq("Compile[{{x, _Real}}, Gamma[x]][3.5] == Gamma[3.5]", "True", 0);
    assert_eval_eq("Compile[{{x, _Real}, {y, _Real}}, Binomial[x, y]][5.5, 2.] == 12.375",
                   "True", 0);
    assert_eval_eq("Head[Compile[{{x, _Real}}, FractionalPart[x]][2.5]]", "Real", 0);
    assert_eval_eq("Head[Compile[{{x, _Real}}, Im[x]][2.5]]", "Real", 0);
    assert_eval_eq("Compile[{{x, _Real}, {y, _Real}}, x^y][2., 0.5] == Sqrt[2.]", "True", 0);
}

/* Integer-ONLY heads: no real counterpart, so they compile over CT_INT or not at
 * all.  Every case is checked against the interpreter's own answer, including
 * the ones where the interpreter declines (Infinity, an unevaluated call) — that
 * is where a compiled path is most tempted to invent something. */
void test_cf_integer_only_heads(void) {
    /* GCD / LCM: non-negative, 0 handled the interpreter's way, and n-ary
     * because both are Flat. */
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, GCD[a, b]][12, 18]", "6", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, GCD[a, b]][-12, 18]", "6", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, GCD[a, b]][0, 0]", "0", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, GCD[a, b]][0, 5]", "5", 0);
    assert_eval_eq("Head[Compile[{{a, _Integer}, {b, _Integer}}, GCD[a, b]][12, 18]]",
                   "Integer", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}, {c, _Integer}}, GCD[a, b, c]]"
                   "[12, 18, 30]", "6", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, LCM[a, b]][4, 6]", "12", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, LCM[a, b]][-4, 6]", "12", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, LCM[a, b]][0, 5]", "0", 0);
    /* Divide-before-multiply: this product fits an int64 only because the gcd is
     * taken out first, so it fails a naive |a b|/g. */
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, LCM[a, b]]"
                   "[3037000499, 3037000493] === LCM[3037000499, 3037000493]", "True", 0);

    /* Predicates, which are Mod-and-compare rather than opcodes of their own —
     * so they inherit MOD_I's zero-divisor guard. */
    assert_eval_eq("Compile[{{n, _Integer}}, EvenQ[n]][4]", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, EvenQ[n]][5]", "False", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, EvenQ[n]][-4]", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, EvenQ[n]][0]", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, OddQ[n]][-3]", "True", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Divisible[a, b]][12, 4]", "True", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Divisible[a, b]][12, 5]", "False", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, Divisible[a, b]][12, 0]"
                   " === Divisible[12, 0]", "True", 0);

    /* IntegerLength / IntegerExponent, with and without a base. */
    assert_eval_eq("Compile[{{n, _Integer}}, IntegerLength[n]][12345]", "5", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, IntegerLength[n]][0]", "0", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, IntegerLength[n]][-12345]", "5", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {b, _Integer}}, IntegerLength[n, b]][255, 16]",
                   "2", 0);
    /* Base 1 is not a base; the interpreter says so and the compiled path must
     * let it. */
    assert_eval_eq("Compile[{{n, _Integer}, {b, _Integer}}, IntegerLength[n, b]][255, 1]"
                   " === IntegerLength[255, 1]", "True", 0);
    assert_eval_eq("Compile[{{n, _Integer}, {b, _Integer}}, IntegerExponent[n, b]][24, 2]",
                   "3", 0);
    assert_eval_eq("Compile[{{n, _Integer}}, IntegerExponent[n]][100]", "2", 0);
    /* IntegerExponent[0, b] is Infinity — not a machine integer. */
    assert_eval_eq("Compile[{{n, _Integer}, {b, _Integer}}, IntegerExponent[n, b]][0, 2]",
                   "Infinity", 0);

    /* PowerMod, the one ternary opcode.  Includes the modular INVERSE (negative
     * exponent), a modulus too large for an int64 intermediate — which defers to
     * the interpreter and must still agree — and a non-invertible case. */
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}, {m, _Integer}}, PowerMod[a, b, m]]"
                   "[3, 100, 7]", "4", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}, {m, _Integer}}, PowerMod[a, b, m]]"
                   "[2, 10, 1000]", "24", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}, {m, _Integer}}, PowerMod[a, b, m]]"
                   "[7, 100, 1000000007] === PowerMod[7, 100, 1000000007]", "True", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}, {m, _Integer}}, PowerMod[a, b, m]]"
                   "[3, -1, 7]", "5", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}, {m, _Integer}}, PowerMod[a, b, m]]"
                   "[3, 100, 4000000000] === PowerMod[3, 100, 4000000000]", "True", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}, {m, _Integer}}, PowerMod[a, b, m]]"
                   "[2, -1, 4] === PowerMod[2, -1, 4]", "True", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}, {m, _Integer}}, PowerMod[a, b, m]]"
                   "[2, 3, 0] === PowerMod[2, 3, 0]", "True", 0);

    /* REAL arguments must not reach any of this: the interpreter has no real
     * GCD either, so the body simply bails and it answers. */
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{x, _Real}, {y, _Real}}, GCD[x, y]]",
                   "False", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, EvenQ[x]]", "False", 0);
}

/* RuntimeOptions -> {"CatchMachineIntegerOverflow" -> False} keeps the wrapped
 * int64 instead of deferring to the interpreter.  It is opt-in precisely because
 * it makes the object answer differently from the interpreter, so the test
 * asserts that difference rather than a value: that is the whole feature. */
void test_cf_runtime_options(void) {
    assert_eval_eq("Options[Compile]",
                   "{RuntimeAttributes -> {}, "
                   "RuntimeOptions -> {\"CatchMachineIntegerOverflow\" -> True}}", 0);

    /* Default: promote via the interpreter. */
    assert_eval_eq("Compile[{{a, _Integer}}, a*a][4000000000] === 4000000000^2", "True", 0);
    /* Off: wrap.  Not equal to the exact answer, and equal to it mod 2^64. */
    assert_eval_eq("Compile[{{a, _Integer}}, a*a, "
                   "RuntimeOptions -> {\"CatchMachineIntegerOverflow\" -> False}]"
                   "[4000000000] === 4000000000^2", "False", 0);
    assert_eval_eq("Mod[Compile[{{a, _Integer}}, a*a, "
                   "RuntimeOptions -> {\"CatchMachineIntegerOverflow\" -> False}]"
                   "[4000000000], 2^64] == Mod[4000000000^2, 2^64]", "True", 0);
    /* The Wolfram shorthands. */
    assert_eval_eq("Compile[{{a, _Integer}}, a*a, RuntimeOptions -> \"Speed\"]"
                   "[4000000000] === 4000000000^2", "False", 0);
    assert_eval_eq("Compile[{{a, _Integer}}, a*a, RuntimeOptions -> \"Quality\"]"
                   "[4000000000] === 4000000000^2", "True", 0);
    /* Turning it off must not change any result that FITS. */
    assert_eval_eq("Compile[{{a, _Integer}}, a*a, RuntimeOptions -> \"Speed\"][7]", "49", 0);

    /* Composes with the other option, in either order. */
    assert_eval_eq("Compile[{{a, _Integer}}, a*a, RuntimeOptions -> \"Speed\", "
                   "RuntimeAttributes -> Listable][{3, 4}]", "{9, 16}", 0);
    assert_eval_eq("Compile[{{a, _Integer}}, a*a, RuntimeAttributes -> Listable, "
                   "RuntimeOptions -> \"Speed\"][{3, 4}]", "{9, 16}", 0);

    /* A setting we cannot honour leaves Compile[...] unevaluated rather than
     * being quietly dropped — the same rule RuntimeAttributes follows. */
    assert_eval_eq("Head[Compile[{{a, _Integer}}, a, RuntimeOptions -> {\"Nope\" -> False}]]",
                   "Compile", 0);
    assert_eval_eq("Head[Compile[{{a, _Integer}}, a, "
                   "RuntimeOptions -> {\"CatchMachineIntegerOverflow\" -> 7}]]", "Compile", 0);

    /* SetOptions moves the default, and the restore takes. */
    assert_eval_eq("SetOptions[Compile, RuntimeOptions -> "
                   "{\"CatchMachineIntegerOverflow\" -> False}]; "
                   "Compile[{{a, _Integer}}, a*a][4000000000] === 4000000000^2", "False", 0);
    assert_eval_eq("SetOptions[Compile, RuntimeOptions -> "
                   "{\"CatchMachineIntegerOverflow\" -> True}]; "
                   "Compile[{{a, _Integer}}, a*a][4000000000] === 4000000000^2", "True", 0);
}

/* Threading and the ARRAY signature meet in one place: a rank-r parameter
 * consumes r levels, so threading is over the levels above it — and a packed
 * NDArray must answer the same as the List it packs, or the object would
 * contradict itself depending on how its input happened to be stored. */
void test_cf_runtime_attributes_arrays(void) {
    /* Depth == rank: one whole argument, no threading. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Total[v], RuntimeAttributes -> Listable]"
                   "[{1., 2., 3.}] == 6", "True", 0);
    /* Depth > rank: thread over the rows. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Total[v], RuntimeAttributes -> Listable]"
                   "[{{1., 2.}, {3., 4.}}]", "{3.0, 7.0}", 0);
    assert_eval_eq("Compile[{{m, _Real, 2}}, Total[Flatten[m]], RuntimeAttributes -> Listable]"
                   "[{{{1., 2.}, {3., 4.}}, {{5., 6.}, {7., 8.}}}]", "{10.0, 26.0}", 0);

    /* The same three answers from packed arrays, which come back packed. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Total[v], RuntimeAttributes -> Listable]"
                   "[NDArray[{1., 2., 3.}]] == 6", "True", 0);
    assert_eval_eq("Normal[Compile[{{v, _Real, 1}}, Total[v], RuntimeAttributes -> Listable]"
                   "[NDArray[{{1., 2.}, {3., 4.}}]]]", "{3.0, 7.0}", 0);
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, Total[v], RuntimeAttributes -> Listable]"
                   "[NDArray[{{1., 2.}, {3., 4.}}]]]", "NDArray", 0);

    /* A rank-2 array over a SCALAR parameter threads twice and is re-packed once
     * at the top — one rank-2 array back, not a list of rows. */
    assert_eval_eq("Normal[Compile[{{x, _Real}}, x^2, RuntimeAttributes -> Listable]"
                   "[NDArray[{{1., 2.}, {3., 4.}}]]]", "{{1.0, 4.0}, {9.0, 16.0}}", 0);

    /* A result that cannot be packed stays a List. */
    assert_eval_eq("Compile[{{x, _Real}}, x > 0, RuntimeAttributes -> Listable]"
                   "[NDArray[{1., -1.}]]", "{True, False}", 0);

    /* The argument array is borrowed, not consumed: threading must not disturb
     * the caller's value. */
    assert_eval_eq("ndra = NDArray[{{1., 2.}, {3., 4.}}]; "
                   "Compile[{{x, _Real}}, x^2, RuntimeAttributes -> Listable][ndra]; "
                   "Normal[ndra]", "{{1.0, 2.0}, {3.0, 4.0}}", 0);
}

/* The language-level contract: prints and returns Null, and leaves anything
 * that is not a CompiledFunction unevaluated. */
void test_compile_print_builtin(void) {
    assert_eval_eq("CompilePrint[Compile[{x}, x^2]]", "Null", 0);
    assert_eval_eq("CompilePrint[3]", "CompilePrint[3]", 0);
    assert_eval_eq("CompilePrint[]", "CompilePrint[]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_cf_integer);
    TEST(test_cf_order);
    TEST(test_cf_real);
    TEST(test_cf_complex);
    TEST(test_cf_symbolic_fallback);
    TEST(test_cf_uncompilable_fallback);
    TEST(test_cf_procedural);
    TEST(test_cf_object_and_arity);
    TEST(test_cf_array_argspec);
    TEST(test_cf_part);
    TEST(test_cf_runtime_attributes);
    TEST(test_cf_runtime_attributes_arrays);
    TEST(test_cf_integer_overflow);
    TEST(test_cf_integer_closed_heads);
    TEST(test_cf_integer_only_heads);
    TEST(test_cf_integer_arrays);
    TEST(test_cf_delegated_array_heads);
    TEST(test_cf_ordering);
    TEST(test_cf_narrowing_kernels_over_arrays);
    TEST(test_cf_integer_binary_kernels_over_arrays);
    TEST(test_cf_runtime_options);
    TEST(test_compile_diagnostics);
    TEST(test_disasm_scalar);
    TEST(test_disasm_kernel_names);
    TEST(test_disasm_control_flow);
    TEST(test_disasm_arrays);
    TEST(test_disasm_uncompiled);
    TEST(test_disasm_callee_program);
    TEST(test_compile_print_builtin);

    printf("All CompiledFunction tests passed!\n");
    return 0;
}
