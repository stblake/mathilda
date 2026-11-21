/* Unit tests for LinkedList arguments to the vector-numeric, structural and
 * symbolic heads (worker C of the "LinkedList everywhere" campaign).
 *
 * The contract (LINKED_LIST_DESIGN.md section 3.2): H[ll, ...] gives what
 * H[Normal[ll], ...] gives -- held, with every LinkedList in the answer
 * replaced by its elements -- read from the chunks (the materialisation
 * counter does not move) and leaving the argument unchanged. Machine answers
 * agree with the packed List (Normal of a typed list) to the bit.
 */
#include "print.h"
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "parse.h"
#include "eval.h"
#include "linkedlist.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* nchk[h, l]: h on a LinkedList of l agrees with h[Normal[ll]], held, without
 * a materialisation and without changing the list (as test_linkedlist.c's
 * llchk). nchk2 makes every element of l a LinkedList first: a List of
 * LinkedList rows. */
static const char* NCHK =
    "nchk[h_, l0_] := Module[{ll = LinkedList[l0], l, c0, r, c1}, l = Normal[ll]; "
    "c0 = Internal`LinkedListConversionCount[]; r = h[ll]; "
    "c1 = Internal`LinkedListConversionCount[]; "
    "{c1 == c0, (Hold[Evaluate[r]] //. x_LinkedList :> With[{v = Normal[x]}, v /; True]) === "
    "Hold[Evaluate[h[l]]], Normal[ll] === l}]; "
    "nchk2[h_, l0_] := Module[{ll = LinkedList /@ l0, l, c0, r, c1}, l = Normal /@ ll; "
    "c0 = Internal`LinkedListConversionCount[]; r = h[ll]; "
    "c1 = Internal`LinkedListConversionCount[]; "
    "{c1 == c0, (Hold[Evaluate[r]] //. x_LinkedList :> With[{v = Normal[x]}, v /; True]) === "
    "Hold[Evaluate[h[l]]], (Normal /@ ll) === l}]; ";

static void expect_chk(const char* fn, const char* heads, const char* fixtures) {
    char buf[8192];
    snprintf(buf, sizeof buf, "%sUnion[Flatten[Quiet[Table[%s[h, l], {h, %s}, {l, %s}]]]]",
             NCHK, fn, heads, fixtures);
    assert_eval_eq(buf, "{True}", 0);
}

/* Vectors: machine reals (packed Normal), machine integers, complex, mixed and
 * symbolic, and short ones (a plain Normal). */
#define VEC_FIX "{N[Range[-6, 13]]/7, Range[-5, 14], N[Range[12]] - 0.5 I, {1.5, -2., 3., 0.5}, " \
                "{3, 1, 2, 2}, {n4a, n4b, 2, 1/3}, {n4a, 2.5, 1}}"

void test_vector_numerics() {
    expect_chk("nchk", "{Norm, Norm[#, 1] &, Norm[#, 3] &, Norm[#, Infinity] &, Norm[#, 2.5] &, "
               "Normalize, Total, # . # &, Dot[#, Range[Length[#]]] &, Tr, Tr[#, Times] &, "
               "Signature, Clip, Clip[#, {0, 2}] &, Clip[#, {-1.5, 0.75}] &, Chop, Chop[#, 0.6] &, "
               "Rationalize, Rationalize[#, 1/10] &, N, N[#, 20] &, Transpose, ConjugateTranspose, "
               "TensorProduct[#, {1, 2}] &, ListConvolve[{1, -1}, #] &, "
               "ListCorrelate[{1, 2, 1}, #, {2, 2}] &, ListConvolve[#, #, 1] &, Variables, "
               "PossibleZeroQ, Simplify, ComplexExpand, D[#, n4a] &, Dt, Integrate[#, n4a] &, "
               "NextPrime, PolynomialMod[#, 2] &}", VEC_FIX);
}

/* Matrices given as a LinkedList of rows and as a List of LinkedList rows. */
#define MAT_FIX_SMALL "{{{1, 2, 3}, {4, 5, 6}}, {{1.5, 2.}, {-3., 0.25}, {7., 8.}}, {{n4a, 1}, {2, n4b}}}"
#define MAT_FIX "{{{1, 2, 3}, {4, 5, 6}}, {{1.5, 2.}, {-3., 0.25}, {7., 8.}}, {{n4a, 1}, {2, n4b}}, " \
                "N[Table[Sin[i + 2 j], {i, 12}, {j, 9}]]}"

void test_matrix_shapes() {
    const char* heads = "{Transpose, ConjugateTranspose, Diagonal, Diagonal[#, 1] &, Tr, "
                        "UpperTriangularize, LowerTriangularize[#, -1] &, Norm, Norm[#, 1] &, "
                        "Norm[#, Infinity] &, Norm[#, \"Frobenius\"] &, "
                        "# . Table[i - j, {i, Dimensions[#][[2]]}, {j, 2}] &, Transpose[#, {2, 1}] &}";
    expect_chk("nchk", heads, MAT_FIX);
    expect_chk("nchk2", heads, MAT_FIX_SMALL);
}

/* The spectral transforms, to the value (their List forms box a zero
 * imaginary part differently for packed and plain Lists). */
void test_transforms() {
    assert_eval_eq("n4t[h_, l_] := Module[{ll = LinkedList[l], c0, r}, c0 = Internal`LinkedListConversionCount[]; "
                   "r = h[ll]; {Internal`LinkedListConversionCount[] == c0, "
                   "Max[Abs[Flatten[Normal /@ Normal[r] - h[Normal[ll]]]]] == 0}]; "
                   "Union[Flatten[Table[n4t[h, l], {h, {Fourier, InverseFourier, Fourier[#, FourierParameters -> {1, -1}] &, "
                   "FourierDCT, FourierDST, FourierDCT[#, 3] &, FourierDST[#, 1] &}}, "
                   "{l, {N[Range[-6, 13]]/7, Range[-5, 14], N[Range[12]] - 0.5 I, {1.5, -2., 3., 0.5}, {3, 1, 2, 2}}}]]]",
                   "{True}", 0);
    assert_eval_eq("Normal[Fourier[LinkedList[{n4a, n4b}]]]", "{(n4a + n4b)/Sqrt[2], (n4a - n4b)/Sqrt[2]}", 0);
}

/* Machine answers are the packed List's to the bit: BLAS dot / gemv / gemm,
 * FFTW, the scaled 2-norm, N's int64 -> double. */
void test_machine_bit_parity() {
    assert_eval_eq("n4v = N[Sin[Range[2000]]]; n4w = N[Cos[Range[2000]]]; "
                   "n4m = N[Table[Sin[i j], {i, 60}, {j, 60}]]; "
                   "{LinkedList[n4v] . LinkedList[n4w] === n4v . n4w, "
                   "Normal[LinkedList[n4m] . n4v[[;; 60]]] === n4m . n4v[[;; 60]], "
                   "Normal[LinkedList[n4v[[;; 60]]] . n4m] === n4v[[;; 60]] . n4m, "
                   "Normal /@ Normal[LinkedList[n4m] . n4m] === n4m . n4m, "
                   "Norm[LinkedList[n4v]] === Norm[n4v], Norm[LinkedList[n4v], 3] === Norm[n4v, 3], "
                   "Normal[Normalize[LinkedList[n4v]]] === Normalize[n4v], "
                   "Normal[Fourier[LinkedList[n4v]]] === Fourier[n4v], "
                   "Normal[FourierDST[LinkedList[n4v]]] === FourierDST[n4v], "
                   "Normal[N[LinkedList[Range[3000]]]] === N[Range[3000]], "
                   "Normal[ListConvolve[{1., -2., 1.}, LinkedList[n4v]]] === ListConvolve[{1., -2., 1.}, n4v], "
                   "Normal[Clip[LinkedList[n4v], {-0.5, 0.5}]] === Clip[n4v, {-0.5, 0.5}], "
                   "Normal[Clip[LinkedList[n4v], {-1/2, 1/2}]] === Clip[n4v, {-1/2, 1/2}]}",
                   "{True, True, True, True, True, True, True, True, True, True, True, True, True}", 0);
    /* nothing was materialised for those (only the explicit Normal calls) */
    assert_eval_eq("n4c0 = Internal`LinkedListConversionCount[]; n4l = LinkedList[n4v]; "
                   "{n4l . n4l, Norm[n4l], Fourier[n4l], N[n4l], Clip[n4l], Chop[n4l], Tr[n4l]}; "
                   "Internal`LinkedListConversionCount[] - n4c0", "0", 0);
    /* a machine-integer list is numericalised slot for slot into a real list */
    assert_eval_eq("n4n = N[LinkedList[{1, 2, 3}]]; {Head[n4n], Normal[n4n]}", "{LinkedList, {1.0, 2.0, 3.0}}", 0);
}

/* Heads that take a LinkedList as a sequence of data: digits, terms,
 * residues, points. */
void test_sequence_data_heads() {
    expect_chk("nchk", "{FromDigits, FromDigits[#, 2] &, FromDigits[#, n4b] &, FromContinuedFraction, "
               "ChineseRemainder[#, {3, 5, 7, 11}] &, ChineseRemainder[{1, 2, 3, 4}, #] &, "
               "DigitCount, DigitCount[#, 2] &, DigitCount[#, 10, 1] &, FindIntegerNullVector, "
               "Interpolation, Interpolation[#][2.5] &, ListInterpolation, ListInterpolation[#][1.5] &, "
               "InterpolatingPolynomial[#, n4x] &, ToNumberField}",
               "{{3, 1, 4, 1}, {1, 0, 1, 1}, {2.5, 1.5, -1., 4.}, {1, Sqrt[2], 3, 2 + Sqrt[2]}, {n4a, n4b, 1, 0}}");
    /* pairs, a mantissa-exponent form, a LinkedList of exponents */
    assert_eval_eq("{Normal[Interpolation[LinkedList[{{1, 2.}, {2, 3.5}, {4, 1.}}]][3]] === "
                   "Interpolation[{{1, 2.}, {2, 3.5}, {4, 1.}}][3], "
                   "FromDigits[LinkedList[{{1, 2, 3}, 2}]], Normal[FromDigits[{{1, 2}, LinkedList[{3, 4}]}]]}",
                   "{True, 123/10, {120, 1200}}", 0);
}

/* Symbolic heads: LinkedLists as expressions mapped over, and as lists of
 * variables, assumptions, moduli, equations. */
void test_symbolic_heads() {
    expect_chk("nchk", "{Simplify, FullSimplify, Simplify[#, n4a > 0] &, Simplify[Sqrt[n4a^2], #] &, "
                       "FullSimplify[Sqrt[n4a^2], #] &, ComplexExpand[Exp[I n4a] + n4b, #] &, "
                       "PossibleZeroQ, Variables, D[#, n4a] &, D[#, {n4a, 2}] &, Dt, Integrate[#, n4a] &, "
                       "Series[#, {n4a, 0, 2}] &, Limit[#, n4a -> 1] &, Grad[n4a^2 n4b, #] &, "
                       "Grad[#, {n4a, n4b}] &, Laplacian[n4a^2 n4b^3, #] &, Div[#, {n4a, n4b}] &, "
                       "Curl[#, {n4a, n4b}] &, CoefficientList[(1 + n4a)^2 n4b, #] &, "
                       "CoefficientRules[(1 + n4a)^2 n4b, #] &, CoefficientArrays[#, {n4a, n4b}] &, "
                       "CoefficientArrays[n4a + 2 n4b, #] &, PolynomialQ[n4a^2 + n4b, #] &, "
                       "PolynomialQ[#, n4a] &, PolynomialMod[#, 2] &, PolynomialMod[n4a^3 + 3 n4b, #] &, "
                       "PolynomialReduce[n4a^2 n4b + n4b^2, #, {n4a, n4b}] &, "
                       "GroebnerBasis[#, {n4a, n4b}] &, GroebnerBasis[{n4a^2 - n4b, n4b^2 - 1}, #] &, "
                       "Minimize[(n4a - 1)^2 + (n4b + 2)^2, #] &, (SeedRandom[7]; NMinimize[(n4a - 1)^2 + n4b^2 + 1, #]) &, "
                       "(SeedRandom[7]; NMaximize[-(n4a - 1)^2 - n4b^2, #]) &, NSolve[#, {n4a, n4b}] &, NSolve}",
               "{{n4a, n4b}, {n4a^2 + 2 n4b, 3 n4a n4b - 1}, {n4a^2 - 1, n4b - 2}, {n4a > 0, n4b > 1}}");
}

/* A LinkedList derivative spec is read as the List one ({x, n}, {{x, y}}):
 * never as an opaque "variable" the expression is free of (D used to answer 0). */
void test_derivative_specs() {
    assert_eval_eq("{D[n4x^2 n4y, LinkedList[{n4x, 2}]], Normal[D[n4x^2 n4y, {LinkedList[{n4x, n4y}]}]], "
                   "Normal /@ Normal[Grad[LinkedList[{n4x^2 + n4y, n4x n4y}], {n4x, n4y}]], "
                   "HornerForm[n4x^2 + 2 n4x n4y + 1, LinkedList[{n4x, n4y}]] === "
                   "HornerForm[n4x^2 + 2 n4x n4y + 1, {n4x, n4y}]}",
                   "{2 n4y, {2 n4x n4y, n4x^2}, {{2 n4x, 1}, {n4y, n4x}}, True}", 0);
}

/* The answer is a new LinkedList; the argument is never changed, and a
 * callback that changes it sees a snapshot. */
void test_result_containers() {
    assert_eval_eq("n4r = LinkedList[{n4a^2, n4b}]; {Head[D[n4r, n4a]], Head[N[LinkedList[{1, 2}]]], "
                   "Head[Clip[LinkedList[{1.5, -3.}]]], Head[Transpose[LinkedList[{{1, 2}, {3, 4}}]]], "
                   "Head[LinkedList[{1., 2.}] . {{1, 2}, {3, 4}}], Head[Fourier[LinkedList[{1, 2}]]], "
                   "Head[FindIntegerNullVector[LinkedList[{1, 2, 3}]]], Normal[n4r]}",
                   "{LinkedList, LinkedList, LinkedList, LinkedList, LinkedList, LinkedList, LinkedList, "
                   "{n4a^2, n4b}}", 0);
    /* nested LinkedLists are walked at every depth by the tree walkers */
    assert_eval_eq("n4q = LinkedList[{LinkedList[{0.5, 1.}], 1.5}]; "
                   "{Normal /@ Normal[Rationalize[n4q]], Normal /@ Normal[N[LinkedList[{LinkedList[{1}], 2}]]]}",
                   "{{{1/2, 1}, 3/2}, {{1.0}, 2.0}}", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_vector_numerics);
    TEST(test_matrix_shapes);
    TEST(test_transforms);
    TEST(test_machine_bit_parity);
    TEST(test_sequence_data_heads);
    TEST(test_symbolic_heads);
    TEST(test_derivative_specs);
    TEST(test_result_containers);

    printf("All LinkedList numeric tests passed!\n");
    symtab_clear();
    return 0;
}
