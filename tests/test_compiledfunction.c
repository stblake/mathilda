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
    assert_eval_eq("Head[Compile[{{v, _Integer, 2}}, v]]", "Compile", 0);  /* no integer dtype */
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

    /* The registered default, and SetOptions changing it. */
    assert_eval_eq("Options[Compile]", "{RuntimeAttributes -> {}}", 0);
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
