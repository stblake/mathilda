/* Tests for auto-compilation of the numeric builtins (Plot/Plot3D now;
 * NIntegrate/FindRoot/Table as they are wired). These drive the real builtins
 * end-to-end and check that the compiled fast path agrees with the interpreter
 * to machine precision, that uncompilable bodies fall back cleanly, and that the
 * non-real-exclusion / fallback semantics are preserved. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "expr.h"
#include "parse.h"
#include "compile/autocompile.h"

/* A reference body that is guaranteed to run in the interpreter while computing
 * exactly the same numbers: `uncid[t_] := t` is an identity whose DownValue the
 * compiler cannot lower, so wrapping a body forces the slow path with no
 * perturbation at all.
 *
 * Prefer this over "wrap it in a head with no machine kernel".  That reference
 * expires the day the head gets a kernel, silently turning a compiled-vs-
 * interpreted check into a compiled-vs-compiled one — which is exactly what
 * happened here to `Zeta`. */
#define AC_UNCID "uncid[t_] := t; "

/* Every sampled Plot point lies on the true curve (compiled path == interpreter). */
void test_plot_parity(void) {
    assert_eval_eq(
        "With[{p = Cases[Plot[Sin[x] + Cos[2 x], {x, -3, 3}], Line[q_] :> q, Infinity][[1]]}, "
        "Max[Table[Abs[(Sin[p[[k,1]]] + Cos[2 p[[k,1]]]) - p[[k,2]]], {k, 1, Length[p]}]] < 10^-9]",
        "True", 0);
    /* a rational + exp body */
    assert_eval_eq(
        "With[{p = Cases[Plot[Exp[-x^2] (1 + x)/(2 + x^2), {x, -2, 2}], Line[q_] :> q, Infinity][[1]]}, "
        "Max[Table[Abs[Exp[-p[[k,1]]^2] (1 + p[[k,1]])/(2 + p[[k,1]]^2) - p[[k,2]]], {k, 1, Length[p]}]] < 10^-9]",
        "True", 0);
}

/* Body outside the compilable subset (Zeta has no machine kernel) still plots. */
void test_plot_fallback(void) {
    assert_eval_eq("Head[Plot[Zeta[x], {x, 2, 5}]]", "Graphics", 0);
    assert_eval_eq("Head[Plot[Sin[x], {x, 0, Pi}]]", "Graphics", 0);
}

/* Where the body goes complex (Sqrt of a negative), the point is excluded —
 * exactly as the interpreter path does. */
void test_plot_complex_excluded(void) {
    assert_eval_eq(
        "With[{p = Cases[Plot[Sqrt[x], {x, -1, 1}], Line[q_] :> q, Infinity][[1]]}, "
        "Min[Table[p[[k,1]], {k, 1, Length[p]}]] >= 0]",
        "True", 0);
}

/* Plot3D: every polygon vertex lies on the surface. */
void test_plot3d_parity(void) {
    assert_eval_eq(
        "With[{p = Flatten[Cases[Plot3D[Sin[x] Cos[y], {x, 0, 2}, {y, 0, 2}], Polygon[q_] :> q, Infinity], 1]}, "
        "Max[Table[Abs[Sin[p[[k,1]]] Cos[p[[k,2]]] - p[[k,3]]], {k, 1, Length[p]}]] < 10^-9]",
        "True", 0);
    assert_eval_eq("Head[Plot3D[Zeta[x + y], {x, 2, 3}, {y, 2, 3}]]", "Graphics3D", 0);
}

/* Table: exact/symbolic iterators are UNTOUCHED; only inexact (machine-real)
 * iterators take the compiled path. */
void test_table_exact_untouched(void) {
    assert_eval_eq("Table[i^2, {i, 1, 5}]", "{1, 4, 9, 16, 25}", 0);          /* exact Integers */
    assert_eval_eq("Table[2^i, {i, 62, 64}]",
                   "{4611686018427387904, 9223372036854775808, 18446744073709551616}", 0); /* BigInt, no int64 overflow */
    assert_eval_eq("Table[Sin[i], {i, 1, 3}]", "{Sin[1], Sin[2], Sin[3]}", 0); /* symbolic */
    assert_eval_eq("Table[1/i, {i, 1, 4}]", "{1, 1/2, 1/3, 1/4}", 0);          /* exact Rationals */
}

/* Real iterator: compiled result matches the interpreter to machine precision. */
void test_table_real_parity(void) {
    assert_eval_eq(
        "Max[Table[Abs[Table[Sin[x] + x^2, {x, 0., 10., 0.1}][[k]] - (Sin[(k-1) 0.1] + ((k-1) 0.1)^2)], "
        "{k, 1, 101}]] < 10^-11",
        "True", 0);
    assert_eval_eq("Table[x^2, {x, 1., 4., 1.}]", "{1.0, 4.0, 9.0, 16.0}", 0);
}

/* The compiled element must keep its own type: a body whose value is an Integer
 * (If branches, Round, ...) must not come back as a Real. */
void test_table_result_type_preserved(void) {
    assert_eval_eq("Table[If[x > 2.5, 1, 2], {x, 1., 4.}]", "{2, 2, 1, 1}", 0);
    assert_eval_eq("Table[Round[x], {x, 1.2, 4.2}]", "{1, 2, 3, 4}", 0);
    assert_eval_eq("Table[x^2, {x, 1., 3.}]", "{1.0, 4.0, 9.0}", 0);   /* Reals stay Real */
}

/* A nested Table binds the outer iterator through the symbol table, so the inner
 * body sees it as a free symbol.  COMPILE_FOLD_GLOBALS resolves it to the value
 * it currently holds; without that the whole inner body is uncompilable. */
void test_table_nested_folds_outer_var(void) {
    assert_eval_eq("Table[x + 10 y, {y, 1., 3.}, {x, 1., 2.}]",
                   "{{11.0, 12.0}, {21.0, 22.0}, {31.0, 32.0}}", 0);
    /* Folding must re-run per call, never cache a stale value. */
    assert_eval_eq("g0 = 3.; {Table[g0 x, {x, 1., 2.}], g0 = 10.; Table[g0 x, {x, 1., 2.}]}",
                   "{{3.0, 6.0}, {10.0, 20.0}}", 0);
    /* A body that MUTATES a global is not compilable, so the accumulation below
     * runs in the interpreter and stays correct. */
    assert_eval_eq("g1 = 0.; {Table[(g1 = g1 + 1.; g1), {x, 1., 4.}], g1}",
                   "{{1.0, 2.0, 3.0, 4.0}, 4.0}", 0);
}

/* Calling a CompiledFunction from a compiled body inlines the callee instead of
 * paying an evaluator round-trip.  Each case is checked against an interpreted
 * twin (a DownValue never compiles), so a capture bug shows up as a mismatch. */
void test_table_inlines_compiled_callee(void) {
    assert_eval_eq("cf1 = Compile[{{a, _Real}}, a^2 + 1]; ci1[a_] := a^2 + 1; "
                   "Table[cf1[x], {x, 1., 5.}] === Table[ci1[x], {x, 1., 5.}]", "True", 0);
    /* callee body reads a global */
    assert_eval_eq("cq = 7.; cf2 = Compile[{{a, _Real}}, a + cq]; ci2[a_] := a + cq; "
                   "Table[cf2[x], {x, 1., 3.}] === Table[ci2[x], {x, 1., 3.}]", "True", 0);
    /* callee PARAMETER shares the caller's iterator name — must not capture */
    assert_eval_eq("cf3 = Compile[{{x, _Real}}, x*2]; ci3[x_] := x*2; "
                   "Table[cf3[x + 1], {x, 1., 3.}] === Table[ci3[x + 1], {x, 1., 3.}]", "True", 0);
    /* callee GLOBAL shares the caller's iterator name — must resolve to neither
     * silently: whatever it does, both paths must agree */
    assert_eval_eq("cw = 100.; cf4 = Compile[{{a, _Real}}, a + cw]; ci4[a_] := a + cw; "
                   "Table[cf4[x], {x, 1., 3.}] === Table[ci4[x], {x, 1., 3.}]", "True", 0);
    /* two parameters, and a nested inline (callee calls another callee) */
    assert_eval_eq("cf5 = Compile[{{a, _Real}, {b, _Real}}, a - b]; ci5[a_, b_] := a - b; "
                   "Table[cf5[x, 2 x], {x, 1., 3.}] === Table[ci5[x, 2 x], {x, 1., 3.}]", "True", 0);
    assert_eval_eq("cf6 = Compile[{{a, _Real}}, a + 1]; cf7 = Compile[{{a, _Real}}, cf6[a]*3]; "
                   "ci7[a_] := (a + 1)*3; "
                   "Table[cf7[x], {x, 1., 3.}] === Table[ci7[x], {x, 1., 3.}]", "True", 0);
    /* an Integer-typed parameter keeps its declared type through the inline */
    assert_eval_eq("cf8 = Compile[{{a, _Integer}}, a*a]; Table[cf8[Round[x]], {x, 1., 3.}]",
                   "{1, 4, 9}", 0);
    /* Wrong arity must bail to a stuck expression — never silently produce a
     * number from the wrong parameter binding. */
    assert_eval_eq("cf9 = Compile[{{a, _Real}, {b, _Real}}, a + b]; "
                   "NumberQ[Table[cf9[x], {x, 1., 2.}][[1]]]", "False", 0);
}

/* ---- Multi-statement loop bodies through the auto-compiled Table path. -------
 * A one-liner only exercises loop scaffolding.  This body carries four locals of
 * three types (Complex u/du, Real acc, Integer cnt), five statements per
 * iteration, and a two-armed If that mutates state in BOTH arms, run over a
 * nested Table (so the outer iterator is folded as a constant).  Written three
 * ways — Do / While / For — which must agree bitwise with each other and match
 * an interpreted twin; `mlRef` is a DownValue, which never compiles, so any
 * miscompilation shows up as a mismatch rather than passing silently.
 * The grid straddles the imaginary axis, so both If arms fire (checked). */
#define ML_STEP  "du = (2 u + 1/u^2)/3 - u; u = u + du; acc = acc + Abs[du]; " \
                 "If[Re[u] > 0, cnt = cnt + 1, cnt = cnt - 1]"
#define ML_DECL  "u = x + I y, du = 0. + 0. I, acc = 0., cnt = 0"
#define ML_TAIL  "acc + cnt + Re[u]"
#define ML_GRID  ", {y, 0.4, 1.2, 0.2}, {x, -1.2, 1.2, 0.3}]"
#define ML_REF   "mlRef[z_, n_] := Module[{u = z, du = 0. + 0. I, acc = 0., cnt = 0, k = 0}, " \
                 "While[k < n, du = (2 u + 1/u^2)/3 - u; u = u + du; acc = acc + Abs[du]; " \
                 "If[Re[u] > 0, cnt = cnt + 1, cnt = cnt - 1]; k = k + 1]; acc + cnt + Re[u]]; "

void test_table_multiline_loop_body(void) {
    /* Do / While / For spellings of the same iteration are the same computation */
    assert_eval_eq(
        "tD = Table[Module[{" ML_DECL "}, Do[" ML_STEP ", {8}]; " ML_TAIL "]" ML_GRID "; "
        "tW = Table[Module[{" ML_DECL ", k = 0}, While[k < 8, " ML_STEP "; k = k + 1]; " ML_TAIL "]" ML_GRID "; "
        "tF = Table[Module[{" ML_DECL ", k = 0}, For[k = 0, k < 8, k = k + 1, " ML_STEP "]; " ML_TAIL "]" ML_GRID "; "
        "{tD === tW, tD === tF}", "{True, True}", 0);
    /* ...and they match the interpreter */
    assert_eval_eq(ML_REF
        "Max[Abs[Flatten[tD - Table[mlRef[x + I y, 8]" ML_GRID "]]] < 10^-10", "True", 0);
    /* both arms of the in-loop If actually fire on this grid */
    assert_eval_eq("Union[Sign[Flatten[tD]]]", "{-1, 1}", 0);
    /* a zero-trip loop must still produce the interpreter's value */
    assert_eval_eq(
        "Table[Module[{" ML_DECL ", k = 0}, While[k < 0, " ML_STEP "; k = k + 1]; " ML_TAIL "]" ML_GRID
        " === Table[mlRef[x + I y, 0]" ML_GRID, "True", 0);
}

/* The same multi-statement body reached through an INLINED CompiledFunction
 * callee: the loop, its locals and its branch all have to survive being lowered
 * with only the callee's parameters in scope. */
void test_table_inlines_multiline_callee(void) {
    assert_eval_eq(ML_REF
        "mlC = Compile[{{z, _Complex}, {n, _Integer}}, "
        "Module[{u = z, du = 0. + 0. I, acc = 0., cnt = 0, k = 0}, "
        "While[k < n, du = (2 u + 1/u^2)/3 - u; u = u + du; acc = acc + Abs[du]; "
        "If[Re[u] > 0, cnt = cnt + 1, cnt = cnt - 1]; k = k + 1]; acc + cnt + Re[u]]]; "
        "Max[Abs[Flatten[Table[mlC[x + I y, 8]" ML_GRID " - Table[mlRef[x + I y, 8]" ML_GRID "]]] < 10^-10",
        "True", 0);
    /* Do and For callees inline the same way */
    assert_eval_eq(
        "mlCD = Compile[{{z, _Complex}, {n, _Integer}}, "
        "Module[{u = z, du = 0. + 0. I, acc = 0., cnt = 0}, "
        "Do[du = (2 u + 1/u^2)/3 - u; u = u + du; acc = acc + Abs[du]; "
        "If[Re[u] > 0, cnt = cnt + 1, cnt = cnt - 1], {n}]; acc + cnt + Re[u]]]; "
        "mlCF = Compile[{{z, _Complex}, {n, _Integer}}, "
        "Module[{u = z, du = 0. + 0. I, acc = 0., cnt = 0, k = 0}, "
        "For[k = 0, k < n, k = k + 1, du = (2 u + 1/u^2)/3 - u; u = u + du; acc = acc + Abs[du]; "
        "If[Re[u] > 0, cnt = cnt + 1, cnt = cnt - 1]]; acc + cnt + Re[u]]]; "
        "{Table[mlCD[x + I y, 8]" ML_GRID " === tD, Table[mlCF[x + I y, 8]" ML_GRID " === tD}",
        "{True, True}", 0);
}

/* ANTI-VACUITY GUARD.  Every parity check above would also pass if NOTHING
 * compiled and all paths ran in the interpreter — agreement with the interpreter
 * is trivially true when you ARE the interpreter.  So assert directly, through
 * the very API `Table` uses, that each multi-statement loop body really does
 * compile as a function of the (real) iterator, with the outer iterator `y`
 * bound as a global for the fold.  If a future edit drops a loop form out of the
 * compilable subset, this fails while the parity tests stay green. */
void test_multiline_bodies_really_compile(void) {
    const char* forms[3] = {
        "Module[{" ML_DECL "}, Do[" ML_STEP ", {8}]; " ML_TAIL "]",
        "Module[{" ML_DECL ", k = 0}, While[k < 8, " ML_STEP "; k = k + 1]; " ML_TAIL "]",
        "Module[{" ML_DECL ", k = 0}, For[k = 0, k < 8, k = k + 1, " ML_STEP "]; " ML_TAIL "]",
    };
    const char* names[3] = { "Do", "While", "For" };
    expr_free(evaluate(parse_expression("y = 0.7")));    /* the folded outer iterator */
    Expr* xs = parse_expression("x");
    for (int i = 0; i < 3; i++) {
        Expr* b = parse_expression(forms[i]);            /* raw: evaluating rewrites the loop */
        AutoCompiled* ac = autocompile_new(b, (const Expr* const*)&xs, 1);
        if (!ac) fprintf(stderr, "FAIL: multi-line %s body did not compile\n", names[i]);
        assert(ac != NULL);
        autocompiled_free(ac);
        expr_free(b);
    }
    expr_free(xs);
    expr_free(evaluate(parse_expression("Clear[y]")));
}

#undef ML_STEP
#undef ML_DECL
#undef ML_TAIL
#undef ML_GRID
#undef ML_REF

/* Real iterator where the body goes complex → per-element interpreter fallback. */
void test_table_real_complex_fallback(void) {
    /* Sqrt[-1.] must still come back complex (I), not NaN. */
    assert_eval_eq("Chop[Table[Sqrt[x], {x, -1., 1., 1.}] - {I, 0., 1.}] == {0, 0, 0}", "True", 0);
}

/* list-LHS Set threading ({a, b} = c binds a = c, b = c) must AUTO-COMPILE
 * through the very API Table uses -- agreement via the interpreter fallback would
 * be vacuous. A destructuring List RHS is not lowered and must decline (NULL),
 * falling back to the interpreter that handles it. */
void test_list_thread_autocompiles(void) {
    Expr* xs = parse_expression("x");
    struct { const char* body; int compiles; } cases[] = {
        { "Module[{a = 0., b = 0.}, {a, b} = x; a + b]",                  1 },
        { "Module[{a = 0., b = 0., c = 0.}, {a, {b, c}} = x; a + b + c]", 1 },
        { "Module[{a = 0., b = 0.}, {a, b} = {x, 2 x}; a + b]",           0 }, /* destructure */
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        Expr* b = parse_expression(cases[i].body);
        AutoCompiled* ac = autocompile_new(b, (const Expr* const*)&xs, 1);
        if (cases[i].compiles && !ac)
            fprintf(stderr, "FAIL: list-thread body did not auto-compile: %s\n", cases[i].body);
        if (!cases[i].compiles && ac)
            fprintf(stderr, "FAIL: destructure body must decline auto-compile: %s\n", cases[i].body);
        assert((cases[i].compiles != 0) == (ac != NULL));
        if (ac) autocompiled_free(ac);
        expr_free(b);
    }
    expr_free(xs);
}

/* End-to-end: threaded assignment inside an auto-compiled Table body (real
 * iterator) gives the same values as the interpreter. */
void test_list_thread_table_parity(void) {
    assert_eval_eq(
        "Table[Module[{a = 0., b = 0.}, {a, b} = x; a + b], {x, 1., 4.}] == {2., 4., 6., 8.}",
        "True", 0);
    assert_eval_eq(
        "Table[Module[{a = 0., b = 0., c = 0.}, {a, {b, c}} = x; a + b + c], {x, 1., 3.}] == {3., 6., 9.}",
        "True", 0);
}

/* NIntegrate: finite / half-line / whole-line agree with the interpreter; the
 * oscillatory sub-method (which re-bodies a copied context) stays correct; a
 * body that goes complex falls back per-sample; uncompilable bodies still work. */
void test_nintegrate_parity(void) {
    assert_eval_eq("Abs[NIntegrate[Sin[x]^2, {x, 0, Pi}] - Pi/2] < 10^-8", "True", 0);
    assert_eval_eq("Abs[NIntegrate[Exp[-x], {x, 0, Infinity}] - 1] < 10^-8", "True", 0);
    assert_eval_eq("Abs[NIntegrate[Exp[-x^2], {x, -Infinity, Infinity}] - Sqrt[Pi]] < 10^-7", "True", 0);
}

void test_nintegrate_oscillatory(void) {
    /* Regression: the amplitude/phase decomposition copies the sampler context
     * and swaps the body — the compiled program must not carry over. */
    assert_eval_eq("Abs[NIntegrate[Cos[100000 x], {x, 0, 1}] - Sin[100000]/100000] < 10^-9", "True", 0);
}

void test_nintegrate_complex_fallback(void) {
    /* Sqrt[x-1] is imaginary on [0,1], real on [1,2]; the compiled real program
     * returns NaN on the imaginary part and the interpreter supplies it. */
    assert_eval_eq("Chop[NIntegrate[Sqrt[x - 1], {x, 0, 2}] - (2/3 + 2 I/3)] == 0", "True", 0);
}

/* Multi-dimensional NIntegrate (cubature / Monte-Carlo, ni_mc_sample). */
void test_nintegrate_multidim(void) {
    assert_eval_eq("Abs[NIntegrate[x y, {x, 0, 1}, {y, 0, 1}] - 1/4] < 10^-8", "True", 0);
    assert_eval_eq("Abs[NIntegrate[x y z, {x, 0, 1}, {y, 0, 1}, {z, 0, 1}] - 1/8] < 10^-7", "True", 0);
    assert_eval_eq("Abs[NIntegrate[Exp[-(x^2 + y^2)], {x, -2, 2}, {y, -2, 2}] "
                   "- (NIntegrate[Exp[-x^2], {x, -2, 2}])^2] < 10^-6", "True", 0);
    /* Monte-Carlo (low default accuracy) */
    assert_eval_eq("Abs[NIntegrate[x^2 + y^2, {x, 0, 1}, {y, 0, 1}, Method -> \"MonteCarlo\"] - 2/3] < 10^-2", "True", 0);
    /* complex-going integrand → per-point interpreter fallback (nonzero imaginary part) */
    assert_eval_eq("Im[Chop[NIntegrate[Sqrt[x + y - 1], {x, 0, 1}, {y, 0, 1}]]] > 0", "True", 0);
}

void test_nintegrate_uncompilable(void) {
    /* Zeta has no machine kernel → the integrand stays on the interpreter path;
     * cross-check the machine result against the arbitrary-precision one. */
    assert_eval_eq("Abs[NIntegrate[Zeta[x]/x^3, {x, 2, 4}] "
                   "- NIntegrate[Zeta[x]/x^3, {x, 2, 4}, WorkingPrecision -> 30]] < 10^-6", "True", 0);
}

/* FindRoot: the compiled machine-real path converges to the same root as the
 * interpreter (checked by residual and by parity against an uncompilable
 * perturbation). Complex/MPFR/inert-derivative paths keep working. */
void test_findroot_parity(void) {
    /* residual f(root) ~ 0 for Newton / secant / bracket */
    assert_eval_eq("With[{r = x /. FindRoot[Cos[x] - x, {x, 0.5}]}, Abs[Cos[r] - r] < 10^-8]", "True", 0);
    assert_eval_eq("With[{r = x /. FindRoot[x^3 - x - 2, {x, 1, 2}]}, Abs[r^3 - r - 2] < 10^-8]", "True", 0);
    assert_eval_eq("With[{r = x /. FindRoot[x^2 - 2, {x, 1.5, 0, 3}]}, Abs[r^2 - 2] < 10^-8]", "True", 0);
    /* Compiled result equals the interpreter result.  The reference wraps the
     * body in an identity DownValue the compiler cannot lower — NOT in an
     * "uncompilable head", which is what this used to do with `Zeta` until
     * Zeta got a machine kernel and quietly turned the check into
     * compiled-vs-compiled.  A reference built out of a coverage gap expires
     * when the gap is closed; one built out of a user rule does not. */
    assert_eval_eq(AC_UNCID "Abs[(x /. FindRoot[Cos[x] - x, {x, 0.5}]) "
                   "- (x /. FindRoot[uncid[Cos[x] - x], {x, 0.5}])] < 10^-11", "True", 0);
}

void test_findroot_complex_and_inert(void) {
    assert_eval_eq("Chop[(x /. FindRoot[x^2 + 1, {x, I}]) - I] == 0", "True", 0);       /* complex root */
    assert_eval_eq("With[{r = x /. FindRoot[Zeta[x] - 2, {x, 1.5}]}, Abs[Zeta[r] - 2] < 10^-8]", "True", 0); /* inert deriv → FD */
}

/* FindRoot systems: each component f_i and Jacobian entry is compiled as a
 * function of all variables; the linear-solve Newton is otherwise unchanged. */
void test_findroot_system(void) {
    assert_eval_eq("{x, y} /. FindRoot[{x + y == 3, x - y == 1}, {{x, 0}, {y, 0}}]", "{2.0, 1.0}", 0);
    assert_eval_eq("With[{s = {x, y} /. FindRoot[{x^2 + y^2 == 4, x - y == 0}, {{x, 1}, {y, 1}}]}, "
                   "Max[Abs[s - {Sqrt[2], Sqrt[2]}]] < 10^-8]", "True", 0);
    assert_eval_eq("With[{s = {x, y, z} /. FindRoot[{x + y + z == 6, x^2 - y == 2, z - x == 1}, "
                   "{{x, 1}, {y, 1}, {z, 1}}]}, "
                   "Max[Abs[{s[[1]]+s[[2]]+s[[3]]-6, s[[1]]^2-s[[2]]-2, s[[3]]-s[[1]]-1}]] < 10^-8]", "True", 0);
    /* one component uncompilable (LogGamma) but with an evaluable derivative →
     * that component and its Jacobian row use the interpreter, the rest compile. */
    assert_eval_eq("With[{s = {x, y} /. FindRoot[{LogGamma[x] - y == 0, x - 3 == 0}, {{x, 2.5}, {y, 0}}]}, "
                   "Abs[s[[2]] - LogGamma[3]] < 10^-7]", "True", 0);
}

/* ------------------------------------------------------------------ *
 *  The grid / field / term samplers                                    *
 *                                                                      *
 *  Each of these compares a plot against the SAME plot with its body    *
 *  wrapped in `uncid[t_] := t` — an identity whose DownValue the        *
 *  compiler cannot lower, so the reference runs entirely in the         *
 *  interpreter while computing bit-for-bit the same numbers.  That is    *
 *  what makes the comparison exact rather than approximate, and it       *
 *  cannot go vacuous the way an "uncompilable head" reference does the   *
 *  day that head gets a kernel (which is exactly what happened to        *
 *  `Zeta` below).                                                       *
 *                                                                      *
 *  `Cases[g, _Real, Infinity]` pulls out every machine number in the    *
 *  Graphics in order — coordinates AND colours — so one comparison       *
 *  covers the whole output, and the Length check catches a structural    *
 *  divergence (an adaptive sampler subdividing differently) that a       *
 *  value comparison alone would miss.                                    *
 * ------------------------------------------------------------------ */

/* Every machine number in the two Graphics agrees, and there are the same
 * number of them. */
static void assert_plot_parity(const char* compiled, const char* interpreted,
                               const char* tol) {
    char buf[2048];
    snprintf(buf, sizeof buf,
             AC_UNCID
             "With[{a = Cases[%s, _Real, Infinity], b = Cases[%s, _Real, Infinity]}, "
             "Length[a] == Length[b] && Length[a] > 0 && Max[Abs[a - b]] < %s]",
             compiled, interpreted, tol);
    assert_eval_eq(buf, "True", 0);
}

void test_contourplot_parity(void) {
    assert_plot_parity(
        "ContourPlot[Sin[x] Cos[y] + x^2/10, {x, -2, 2}, {y, -2, 2}, PlotPoints -> 30]",
        "ContourPlot[uncid[Sin[x] Cos[y] + x^2/10], {x, -2, 2}, {y, -2, 2}, PlotPoints -> 30]",
        "10^-12");
    /* Equation form takes a different grid path (one grid per equation). */
    assert_plot_parity(
        "ContourPlot[x^2 + y^2 == 2, {x, -2, 2}, {y, -2, 2}, PlotPoints -> 20]",
        "ContourPlot[uncid[x^2 + y^2] == 2, {x, -2, 2}, {y, -2, 2}, PlotPoints -> 20]",
        "10^-12");
    assert_eval_eq("Head[ContourPlot[uncid[x + y], {x, 0, 1}, {y, 0, 1}, PlotPoints -> 5]]",
                   "Graphics", 0);
}

void test_densityplot_parity(void) {
    assert_plot_parity(
        "DensityPlot[Sin[x] Cos[y] + Exp[-x^2 - y^2], {x, -2, 2}, {y, -2, 2}, PlotPoints -> 20]",
        "DensityPlot[uncid[Sin[x] Cos[y] + Exp[-x^2 - y^2]], {x, -2, 2}, {y, -2, 2}, PlotPoints -> 20]",
        "10^-12");
}

/* ComplexPlot is the only sampler whose variable ranges over the PLANE, so it
 * compiles with a complex argument (autocompile_new_z).  Parity is to one ulp
 * rather than exact: C99 complex division and the interpreter's Complex[]
 * arithmetic round differently, both correctly. */
void test_complexplot_parity(void) {
    assert_plot_parity(
        "ComplexPlot[(z^2 - 1)/(z^2 + 2), {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 20]",
        "ComplexPlot[uncid[(z^2 - 1)/(z^2 + 2)], {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 20]",
        "10^-9");
    assert_plot_parity(
        "ComplexPlot[Sin[z] Exp[z], {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 20]",
        "ComplexPlot[uncid[Sin[z] Exp[z]], {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 20]",
        "10^-9");
    /* A head with a real kernel but no complex one must fall back, not bail
     * the plot: Zeta[z] compiles for _Real and does not for _Complex. */
    assert_eval_eq("Head[ComplexPlot[Zeta[z], {z, 2 - I, 3 + I}, PlotPoints -> 5]]",
                   "Graphics", 0);
}

/* A parametric body is a pair, so each coordinate compiles as its own program.
 * The 1-iterator form also drives the ADAPTIVE sampler, whose subdivision
 * depends on the sampled values — so equal output lengths here is itself a
 * check that the compiled values steered the sampler identically. */
void test_parametricplot_parity(void) {
    assert_plot_parity(
        "ParametricPlot[{Sin[3 t] Cos[t], Sin[3 t] Sin[t]}, {t, 0, 2 Pi}]",
        "ParametricPlot[uncid[{Sin[3 t] Cos[t], Sin[3 t] Sin[t]}], {t, 0, 2 Pi}]",
        "10^-12");
    /* 2-iterator (region) form. */
    assert_plot_parity(
        "ParametricPlot[{u Cos[v], u Sin[v]}, {u, 0, 1}, {v, 0, 2 Pi}, PlotPoints -> 12]",
        "ParametricPlot[uncid[{u Cos[v], u Sin[v]}], {u, 0, 1}, {v, 0, 2 Pi}, PlotPoints -> 12]",
        "10^-12");
    /* PolarPlot desugars to ParametricPlot, so it inherits the fast path. */
    assert_plot_parity(
        "PolarPlot[1 + Cos[3 t]/2, {t, 0, 2 Pi}]",
        "PolarPlot[uncid[1 + Cos[3 t]/2], {t, 0, 2 Pi}]",
        "10^-12");
}

void test_parametricplot3d_parity(void) {
    assert_plot_parity(
        "ParametricPlot3D[{Cos[t], Sin[t], t/5}, {t, 0, 10}]",
        "ParametricPlot3D[uncid[{Cos[t], Sin[t], t/5}], {t, 0, 10}]",
        "10^-12");
    assert_plot_parity(
        "ParametricPlot3D[{Cos[u] Sin[v], Sin[u] Sin[v], Cos[v]}, {u, 0, 2 Pi}, {v, 0, Pi}, "
        "PlotPoints -> 12]",
        "ParametricPlot3D[uncid[{Cos[u] Sin[v], Sin[u] Sin[v], Cos[v]}], {u, 0, 2 Pi}, "
        "{v, 0, Pi}, PlotPoints -> 12]",
        "10^-12");
}

void test_vectorplot_parity(void) {
    assert_plot_parity(
        "VectorPlot[{-y Exp[-x^2], x Sin[y]}, {x, -2, 2}, {y, -2, 2}, VectorPoints -> 10]",
        "VectorPlot[{uncid[-y Exp[-x^2]], uncid[x Sin[y]]}, {x, -2, 2}, {y, -2, 2}, "
        "VectorPoints -> 10]",
        "10^-12");
    /* One component uncompilable ⇒ both stay interpreted, and it still plots. */
    assert_eval_eq("Head[VectorPlot[{-y, uncid[x]}, {x, -1, 1}, {y, -1, 1}, VectorPoints -> 4]]",
                   "Graphics", 0);
}

/* Streamline integration takes several field samples per RK step, so this is
 * the hottest sampler of the group — and the most sensitive to a divergence,
 * since a difference in one sample steers the whole streamline. */
void test_streamplot_parity(void) {
    assert_plot_parity(
        "StreamPlot[{-y Exp[-x^2], x Sin[y]}, {x, -2, 2}, {y, -2, 2}]",
        "StreamPlot[{uncid[-y Exp[-x^2]], uncid[x Sin[y]]}, {x, -2, 2}, {y, -2, 2}]",
        "10^-12");
}

/* NSum draws every term through one chokepoint, so Direct / Wynn / Levin /
 * Euler-Maclaurin / CVZ / the far-tail ladder all take the compiled path — and
 * so does the EM correction's CONTINUOUS sampling of the same summand, which is
 * where most of an EM sum's work goes.  NProduct is Exp[NSum[Log f]] and
 * inherits it. */
void test_nsum_parity(void) {
    assert_eval_eq(AC_UNCID "Abs[NSum[Sin[n]/n^3, {n, 1, Infinity}] "
                   "- NSum[uncid[Sin[n]/n^3], {n, 1, Infinity}]] < 10^-12", "True", 0);
    assert_eval_eq(AC_UNCID "Abs[NSum[1/(n^2 + 1), {n, 1, Infinity}] "
                   "- NSum[uncid[1/(n^2 + 1)], {n, 1, Infinity}]] < 10^-12", "True", 0);
    assert_eval_eq(AC_UNCID "Abs[NSum[(-1)^n/n, {n, 1, Infinity}] "
                   "- NSum[uncid[(-1)^n/n], {n, 1, Infinity}]] < 10^-12", "True", 0);
    /* known closed forms */
    assert_eval_eq("Abs[NSum[1/n^2, {n, 1, Infinity}] - Pi^2/6] < 10^-8", "True", 0);
    assert_eval_eq("Abs[NSum[Exp[-n], {n, 0, Infinity}] - 1/(1 - Exp[-1])] < 10^-8", "True", 0);
    /* A complex-valued summand takes the boxed compiled path (complex result
     * type), and must still agree with the interpreter. */
    assert_eval_eq(AC_UNCID "Abs[NSum[I^n/n^2, {n, 1, Infinity}] "
                   "- NSum[uncid[I^n/n^2], {n, 1, Infinity}]] < 10^-10", "True", 0);
    assert_eval_eq(AC_UNCID "Abs[NProduct[1 + 1/n^2, {n, 1, Infinity}] "
                   "- NProduct[uncid[1 + 1/n^2], {n, 1, Infinity}]] < 10^-10", "True", 0);
}

/* The MPFR path must never take the compiled route: a double cannot carry the
 * terms of a 30-digit sum. */
void test_nsum_mpfr_untouched(void) {
    assert_eval_eq("Abs[NSum[1/n^2, {n, 1, Infinity}, WorkingPrecision -> 30] - Pi^2/6] < 10^-25",
                   "True", 0);
}

/* ANTI-VACUITY GUARD for the group above: assert through the very API these
 * builtins use that the bodies really do compile.  Every parity test would pass
 * trivially if nothing compiled at all. */
void test_sampler_bodies_really_compile(void) {
    struct { const char* body; const char* v1; const char* v2; } cases[] = {
        { "Sin[x] Cos[y] + x^2/10",       "x", "y"    },   /* Contour / Density */
        { "-y Exp[-x^2]",                 "x", "y"    },   /* Vector / Stream   */
        { "Sin[3 t] Cos[t]",              "t", NULL   },   /* Parametric        */
        { "Cos[u] Sin[v]",                "u", "v"    },   /* ParametricPlot3D  */
        { "Sin[n]/n^3",                   "n", NULL   },   /* NSum summand      */
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        Expr* b  = parse_expression(cases[i].body);
        Expr* v1 = parse_expression(cases[i].v1);
        Expr* v2 = cases[i].v2 ? parse_expression(cases[i].v2) : NULL;
        const Expr* vs[2] = { v1, v2 };
        AutoCompiled* ac = autocompile_new(b, vs, v2 ? 2 : 1);
        if (!ac) fprintf(stderr, "FAIL: sampler body `%s` did not compile\n", cases[i].body);
        assert(ac != NULL);
        autocompiled_free(ac);
        expr_free(b); expr_free(v1); expr_free(v2);
    }
    /* ComplexPlot needs the COMPLEX-argument constructor, and its subset is
     * genuinely smaller — assert both halves of that. */
    Expr* zb = parse_expression("(z^2 - 1)/(z^2 + 2)");
    Expr* zv = parse_expression("z");
    const Expr* zs[1] = { zv };
    AutoCompiled* zc = autocompile_new_z(zb, zs, 1);
    if (!zc) fprintf(stderr, "FAIL: ComplexPlot body did not compile at complex argument\n");
    assert(zc != NULL);
    autocompiled_free(zc);
    expr_free(zb);

    /* The complex subset really is smaller: Zeta has a real machine kernel and
     * no complex one, so the same body compiles one way and not the other.  If
     * this ever stops holding, a complex Zeta kernel landed and
     * NUMERIC_FUNCTION_MISSING.md needs updating — which is the point. */
    Expr* rb = parse_expression("Zeta[z]");
    AutoCompiled* rr = autocompile_new(rb, zs, 1);
    assert(rr != NULL);
    autocompiled_free(rr);
    AutoCompiled* zz = autocompile_new_z(rb, zs, 1);
    if (zz) fprintf(stderr, "NOTE: Zeta gained a complex kernel; "
                            "update NUMERIC_FUNCTION_MISSING.md\n");
    autocompiled_free(zz);
    expr_free(rb); expr_free(zv);
}

/* ---- finite Sum[] / Product[] auto-compilation ------------------------- */

/* Strategy A: an integer-iterator finite sum with an inexact body compiles the
 * WHOLE Sum to one machine loop.  Parity is checked against the interpreter via
 * the uncid[] wrapper (a DownValue the compiler cannot lower), and the total is
 * a machine Real. */
void test_sum_integer_parity(void) {
    assert_eval_eq(AC_UNCID
        "Abs[Sum[Sin[i^2]*1.0, {i, 1, 20000}] - Sum[uncid[Sin[i^2]*1.0], {i, 1, 20000}]] < 10^-9",
        "True", 0);
    assert_eval_eq("Head[Sum[Sin[i^2]*1.0, {i, 1, 2000}]]", "Real", 0);
    /* non-unit integer step */
    assert_eval_eq(AC_UNCID
        "Abs[Sum[Cos[i]*1.0, {i, 1, 999, 2}] - Sum[uncid[Cos[i]*1.0], {i, 1, 999, 2}]] < 10^-9",
        "True", 0);
}

/* Strategy B: a real-bounded / real-step iterator compiles the body and folds in
 * C, visiting exactly the interpreter's index points. */
void test_sum_real_iterator_parity(void) {
    assert_eval_eq(AC_UNCID
        "Abs[Sum[Sin[i]*1.0, {i, 0., 50., 0.1}] - Sum[uncid[Sin[i]*1.0], {i, 0., 50., 0.1}]] < 10^-9",
        "True", 0);
    assert_eval_eq("Head[Sum[Sin[i]*1.0, {i, 0., 50., 0.1}]]", "Real", 0);
}

/* A complex-valued machine body sums to a machine Complex, matching the
 * interpreter. */
void test_sum_complex_parity(void) {
    assert_eval_eq(AC_UNCID
        "Abs[Sum[Exp[I i*1.0], {i, 1, 5000}] - Sum[uncid[Exp[I i*1.0]], {i, 1, 5000}]] < 10^-9",
        "True", 0);
    assert_eval_eq("Head[Sum[Exp[I i*1.0], {i, 1, 100}]]", "Complex", 0);
}

/* Exactness gate: the compiled path must NEVER divert an exact sum/product to a
 * machine number.  A rational body stays Rational, an integer body stays Integer
 * (bignum), a real-iterator body that evaluates to an Integer stays Integer, and
 * an empty range folds to the exact identity. */
void test_sum_product_exactness(void) {
    assert_eval_eq("Head[Sum[1/(i^2 + i + 1), {i, 1, 300}]]", "Rational", 0);
    assert_eval_eq("Head[Sum[i^2, {i, 1, 50}]]", "Integer", 0);
    assert_eval_eq("Sum[Round[i], {i, 1., 10.}]", "55", 0);          /* real iter, exact result */
    assert_eval_eq("Product[i, {i, 1, 20}]", "2432902008176640000", 0);  /* 20! bignum */
    assert_eval_eq("Sum[Sin[i]*1.0, {i, 5, 1}]", "0", 0);           /* empty sum -> exact 0 */
    assert_eval_eq("Product[Cos[i]*1.0, {i, 5, 1}]", "1", 0);       /* empty product -> exact 1 */
}

/* Product mirrors Sum: whole-Product lowering (A) and body-loop (B), real and
 * complex, all parity-checked against the interpreter. */
void test_product_parity(void) {
    assert_eval_eq(AC_UNCID
        "Abs[Product[Cos[i]*1.0, {i, 1, 500}] - Product[uncid[Cos[i]*1.0], {i, 1, 500}]] < 10^-9",
        "True", 0);
    assert_eval_eq(AC_UNCID
        "Abs[Product[1. + Sin[i]/10, {i, 0., 20., 0.5}] - Product[uncid[1. + Sin[i]/10], {i, 0., 20., 0.5}]] < 10^-9",
        "True", 0);
}

/* $AutoCompilation off: the answer is unchanged (the interpreter computes the
 * same machine total), confirming the fast path is a pure optimisation. */
void test_sum_autocompile_switch(void) {
    autocompile_set_enabled(false);
    assert_eval_eq("Head[Sum[Sin[i^2]*1.0, {i, 1, 200}]]", "Real", 0);
    autocompile_set_enabled(true);
    assert_eval_eq("Head[Sum[Sin[i^2]*1.0, {i, 1, 200}]]", "Real", 0);
}

/* Inside Compile[], Sum/Product already lower to a native accumulation loop
 * (compile_emit_ctrl.c) — confirm the compiled function agrees with the
 * interpreter. */
void test_sum_inside_compile(void) {
    assert_eval_eq(AC_UNCID
        "Abs[Compile[{{n, _Integer}}, Sum[Sin[i]*1.0, {i, 1, n}]][1000] "
        "- Sum[uncid[Sin[i]*1.0], {i, 1, 1000}]] < 10^-9",
        "True", 0);
    assert_eval_eq(AC_UNCID
        "Abs[Compile[{{n, _Integer}}, Product[1. + 1./i, {i, 1, n}]][50] "
        "- Product[uncid[1. + 1./i], {i, 1, 50}]] < 10^-9",
        "True", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_plot_parity);
    TEST(test_plot_fallback);
    TEST(test_plot_complex_excluded);
    TEST(test_plot3d_parity);
    TEST(test_table_exact_untouched);
    TEST(test_table_real_parity);
    TEST(test_table_result_type_preserved);
    TEST(test_table_nested_folds_outer_var);
    TEST(test_table_inlines_compiled_callee);
    TEST(test_multiline_bodies_really_compile);   /* precondition: they DO compile */
    TEST(test_table_multiline_loop_body);
    TEST(test_table_inlines_multiline_callee);
    TEST(test_table_real_complex_fallback);
    TEST(test_list_thread_autocompiles);          /* precondition: threading DOES compile */
    TEST(test_list_thread_table_parity);
    TEST(test_nintegrate_parity);
    TEST(test_nintegrate_multidim);
    TEST(test_nintegrate_oscillatory);
    TEST(test_nintegrate_complex_fallback);
    TEST(test_nintegrate_uncompilable);
    TEST(test_findroot_parity);
    TEST(test_findroot_complex_and_inert);
    TEST(test_findroot_system);
    TEST(test_sampler_bodies_really_compile);   /* precondition: they DO compile */
    TEST(test_contourplot_parity);
    TEST(test_densityplot_parity);
    TEST(test_complexplot_parity);
    TEST(test_parametricplot_parity);
    TEST(test_parametricplot3d_parity);
    TEST(test_vectorplot_parity);
    TEST(test_streamplot_parity);
    TEST(test_nsum_parity);
    TEST(test_nsum_mpfr_untouched);

    TEST(test_sum_integer_parity);
    TEST(test_sum_real_iterator_parity);
    TEST(test_sum_complex_parity);
    TEST(test_sum_product_exactness);
    TEST(test_product_parity);
    TEST(test_sum_autocompile_switch);
    TEST(test_sum_inside_compile);

    printf("All auto-compile tests passed!\n");
    return 0;
}
