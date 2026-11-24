/* Builtin-audit regressions: image processing (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* RandomImage#1 (critical): w*h*ch overflowed (2^40 x 2^40 is 0 mod 2^64), so
 * the buffer was under-allocated (SIGBUS), or the request ran for minutes.
 * Checked, an impossible size aborts at once with General::nomem. */
static void test_random_image_huge_size(void) {
    audit_chk("RandomImage#1", "Quiet[CheckAbort[RandomImage[1, {2^40, 2^40}], aborted]]", "aborted");
    audit_chk("RandomImage#1",
              "Quiet[Check[CheckAbort[RandomImage[1, {4611686018427387904, 4}], 0], fired, General::nomem]]",
              "fired");
    audit_chk("RandomImage#1",
              "Quiet[CheckAbort[RandomImage[1, {2^31, 2^31}, ColorSpace -> \"RGB\"], aborted]]", "aborted");
    audit_chk("RandomImage#1", "ImageDimensions[RandomImage[1, {4, 3}]]", "{4, 3}");
    audit_chk("RandomImage#1", "ImageDimensions[RandomImage[2, 5]]", "{5, 5}");
}

/* Opening#1 / Closing#1 (high): the second pass of an opening or closing reused the
 * UNREFLECTED element, so for any element that is not centrally symmetric the pair
 * was not adjoint: the result was shifted, Opening exceeded the image and Closing fell
 * below it, and neither was idempotent.  The second pass now uses the reflection.
 * Checked by the defining laws on random data with asymmetric and even-sized
 * elements: Opening <= f <= Closing pointwise, and both idempotent. */
static void test_opening_closing_asymmetric(void) {
    /* a plateau that exactly fits a 1x2 element is preserved by the opening; the
     * closing of it is itself too (nothing to fill) */
    audit_chk("Opening#1", "ImageData[Opening[Image[{{0., 0., 1., 1., 0., 0.}}, \"Real\"], {{1, 1}}]]",
              "{{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}}");
    audit_true("Opening#1",
               "Module[{g = Image[{{0., 0., 0., 0., 0.}, {0., 0., 1., 0., 0.}, {0., 0., 0., 1., 0.},"
               " {0., 0., 0., 0., 0.}}, \"Real\"]}, ImageData[Opening[g, {{1, 0}, {0, 1}}]] == ImageData[g]]");
    audit_chk("Closing#1", "ImageData[Closing[Image[{{0., 0., 1., 1., 0., 0.}}, \"Real\"], {{1, 1}}]]",
              "{{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}}");
    /* the audit's closing: the bright pixel at {2, 2} survives and the gap between the
     * two pixels is filled.  At the right edge the element restricted to the image reads
     * only the dilation there, max(f[4], f[5]) = 1, so column 5 is 1 as well. */
    audit_chk("Closing#1",
              "ImageData[Closing[Image[{{0., 0., 0., 0., 0.}, {0., 1., 0., 1., 0.}, {0., 0., 0., 0., 0.}}],"
              " {{1, 1, 0}}]]",
              "{{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 1.0, 1.0, 1.0}, {0.0, 0.0, 0.0, 0.0, 0.0}}");
    audit_true("Closing#1",
               "Module[{c = Closing[Image[{{0., 0., 0., 0., 0.}, {0., 1., 0., 1., 0.}, {0., 0., 0., 0., 0.}}],"
               " {{1, 1, 0}}]}, ImageData[Closing[c, {{1, 1, 0}}]] == ImageData[c]]");
    const char* laws =
        "SeedRandom[%d]; Module[{m = RandomReal[1, {12, 13}], k = %s, op, cl},"
        " op = ImageData[Opening[Image[m], k]]; cl = ImageData[Closing[Image[m], k]];"
        " And @@ Flatten[MapThread[LessEqual, {op, m}, 2]]"
        " && And @@ Flatten[MapThread[GreaterEqual, {cl, m}, 2]]"
        " && ImageData[Opening[Image[op], k]] == op && ImageData[Closing[Image[cl], k]] == cl]";
    const char* elems[] = { "{{1, 1}}", "{{1, 1, 0}}", "{{1, 0, 1}, {1, 1, 0}}", "{{0, 0, 1}, {0, 0, 0}, {0, 0, 0}}",
                            "{{1, 1}, {1, 1}}", "{{1, 1, 1, 1}, {1, 1, 1, 1}}" };
    char buf[1024];
    for (int i = 0; i < 6; i++) {
        snprintf(buf, sizeof buf, laws, 3 + i, elems[i]);
        audit_true("Opening#1", buf);
        audit_true("Closing#1", buf);
    }
    /* sanity: a symmetric element is unchanged by the reflection -- the radius form
     * agrees with its explicit square, and a Bit image stays Bit */
    audit_true("Opening#1",
               "SeedRandom[5]; Module[{m = Image[RandomReal[1, {9, 9}]]},"
               " ImageData[Opening[m, 1]] == ImageData[Opening[m, BoxMatrix[1]]]"
               " && ImageData[Closing[m, 1]] == ImageData[Closing[m, BoxMatrix[1]]]]");
    audit_chk("Closing#1", "ImageType[Closing[Image[{{0, 1, 0}, {1, 0, 1}}], {{1, 1}}]]", "\"Bit\"");
}

/* ColorQuantize#1 (high): the median cut split a box at len/2 even inside a run of
 * one colour, so identical pixels were given different colours and box means invented
 * colours (red + blue -> purple).  The cut is now at the value boundary nearest the
 * median, and a grey level is the mean of its identical channels (no luminance
 * round trip), so an image with at most n colours comes back exactly. */
static void test_color_quantize_value_boundary(void) {
    audit_chk("ColorQuantize#1",
              "ImageData[ColorQuantize[Image[{{{1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}, {0., 0., 1.}}}, \"Real\"], 2]]",
              "{{{1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}}}");
    audit_chk("ColorQuantize#1", "ImageData[ColorQuantize[Image[{{0., 0.}, {0., 1.}}, \"Real\"], 2]]",
              "{{0.0, 0.0}, {0.0, 1.0}}");
    audit_chk("ColorQuantize#1", "ImageData[ColorQuantize[Image[{{0., 0.5}}, \"Real\"], 2]]", "{{0.0, 0.5}}");
    /* identical input pixels always map to one output colour, for any n */
    audit_true("ColorQuantize#1",
               "SeedRandom[4]; Module[{px = RandomChoice[{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}, {1., 1., 0.},"
               " {0.2, 0.4, 0.6}}, {6, 7}], q},"
               " q = ImageData[ColorQuantize[Image[px, \"Real\"], 3]];"
               " AllTrue[GatherBy[Transpose[{Flatten[px, 1], Flatten[q, 1]}], First],"
               " Length[Union[#[[All, 2]]]] == 1 &] && Length[Union[Flatten[q, 1]]] <= 3"
               " && Sort[Union[Flatten[ImageData[ColorQuantize[Image[px, \"Real\"], 5]], 1]]] =="
               " Sort[Union[Flatten[px, 1]]]]");
    /* sanity: a ramp still quantises to at most n levels */
    audit_true("ColorQuantize#1",
               "Length[Union[Flatten[ImageData[ColorQuantize[Image[Table[N[(i + j)/32], {i, 16}, {j, 16}], \"Real\"], 4]]]]] == 4");
}

/* ImageCorrelate#1 (high): for an even-sized template the NCC window statistics were
 * anchored at k/2 while the cross term (a convolution with the reversed template)
 * anchors at (k-1)/2, so scores left [-1, 1].  Both now use one anchor: every score is
 * a correlation coefficient, affine-invariant in the template, and a crop scores 1 at
 * its own position. */
static void test_image_correlate_ncc_even(void) {
    audit_true("ImageCorrelate#1",
               "Module[{a = ImageData[ImageCorrelate[Image[{{0., 1., 0., 1.}}], {{2., 7.}}, \"NormalizedCrossCorrelation\"]],"
               " b = ImageData[ImageCorrelate[Image[{{0., 1., 0., 1.}}], {{0., 1.}}, \"NormalizedCrossCorrelation\"]]},"
               " Max[Abs[a - b]] < 1*^-12 && Max[Abs[a]] <= 1 + 1*^-12]");
    audit_true("ImageCorrelate#1",
               "SeedRandom[2]; Module[{m = RandomReal[1, {20, 20}], t, r},"
               " t = m[[5 ;; 8, 7 ;; 10]];"
               " r = ImageData[ImageCorrelate[Image[m], t, \"NormalizedCrossCorrelation\"]];"
               " Max[Abs[r]] <= 1 + 1*^-12 && Abs[Max[r] - 1] < 1*^-12 && Position[r, Max[r]] == {{6, 8}}]");
    /* the score at an interior position is the textbook coefficient of that window */
    audit_true("ImageCorrelate#1",
               "SeedRandom[6]; Module[{m = RandomReal[1, {10, 11}], t = RandomReal[1, {2, 4}], r, w, c},"
               " r = ImageData[ImageCorrelate[Image[m], t, \"NormalizedCrossCorrelation\"]];"
               " w = Flatten[m[[4 ;; 5, 5 ;; 8]]]; c = Flatten[t];"
               " Abs[r[[4, 6]] - (w - Mean[w]) . (c - Mean[c])/(Norm[w - Mean[w]] Norm[c - Mean[c]])] < 1*^-12]");
    /* sanity: odd templates unchanged (a 5x5 crop peaks at its centre) */
    audit_true("ImageCorrelate#1",
               "SeedRandom[1]; Module[{m = RandomReal[1, {20, 20}], r},"
               " r = ImageData[ImageCorrelate[Image[m], m[[5 ;; 9, 7 ;; 11]], \"NormalizedCrossCorrelation\"]];"
               " Position[r, Max[r]] == {{7, 9}} && Abs[Max[r] - 1] < 1*^-12]");
}

/* GradientFilter#1 (high): GradientFilter[image, r] -- the derivative-of-Gaussian
 * form -- was unevaluated.  Its kernels are normalised so a linear ramp of slope s has
 * magnitude s away from the border, and a constant image has gradient 0. */
static void test_gradient_filter_radius(void) {
    audit_true("GradientFilter#1",
               "Module[{g = GradientFilter[Image[{{0., 0., 1., 1.}, {0., 0., 1., 1.}}], 1]},"
               " ImageQ[g] && ImageDimensions[g] == {4, 2}]");
    audit_true("GradientFilter#1",
               "Module[{ramp = Image[Table[0.01 j, {i, 11}, {j, 11}]]},"
               " Max[Abs[ImageData[GradientFilter[ramp, 2]][[4 ;; 8, 4 ;; 8]] - 0.01]] < 1*^-14"
               " && Max[Abs[ImageData[GradientFilter[ramp, {2, 1.5}]][[4 ;; 8, 4 ;; 8]] - 0.01]] < 1*^-14]");
    /* rotation invariance: a diagonal ramp of slope 1/8 per axis has magnitude Sqrt[2]/8 */
    audit_true("GradientFilter#1",
               "Module[{d = Image[Table[(i + j)/80., {i, 13}, {j, 13}]]},"
               " Max[Abs[ImageData[GradientFilter[d, 3]][[5 ;; 9, 5 ;; 9]] - Sqrt[2]/80]] < 1*^-14]");
    audit_true("GradientFilter#1",
               "Max[ImageData[GradientFilter[Image[ConstantArray[0.3, {6, 7}]], 2]]] < 1*^-15");
    /* a step edge peaks on the edge */
    audit_true("GradientFilter#1",
               "Module[{r = ImageData[GradientFilter[Image[Table[If[j > 6, 1., 0.], {i, 9}, {j, 12}]], 2]][[5]]},"
               " MemberQ[{{6}, {7}}, First[Position[r, Max[r]]]]]");
    audit_chk("GradientFilter#1", "GradientFilter[Image[{{0., 1.}}], -1] // Head", "GradientFilter");
    /* sanity: the one-argument Sobel form is unchanged */
    audit_true("GradientFilter#1",
               "Module[{ramp = Image[Table[0.01 j, {i, 9}, {j, 9}]]},"
               " Max[Abs[ImageData[GradientFilter[ramp]][[3 ;; 7, 3 ;; 7]] - 0.01]] < 1*^-14]");
}

int main(void) {
    audit_init();
    TEST(test_random_image_huge_size);
    TEST(test_opening_closing_asymmetric);
    TEST(test_color_quantize_value_boundary);
    TEST(test_image_correlate_ncc_even);
    TEST(test_gradient_filter_radius);
    return audit_finish("audit_image_processing");
}
