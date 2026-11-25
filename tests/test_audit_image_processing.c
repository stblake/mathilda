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
     * two pixels is filled.  Dilation is the standard max f(x - b) (it stamps the element,
     * audit Dilation#1), so the element {{1, 1, 0}} reaches one column to the RIGHT of each
     * pixel: the dilation is 1 on columns 1..4 and the erosion min g(x + b) keeps them; the
     * last column, read only at itself, stays 0. */
    audit_chk("Closing#1",
              "ImageData[Closing[Image[{{0., 0., 0., 0., 0.}, {0., 1., 0., 1., 0.}, {0., 0., 0., 0., 0.}}],"
              " {{1, 1, 0}}]]",
              "{{0.0, 0.0, 0.0, 0.0, 0.0}, {1.0, 1.0, 1.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}}");
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

/* RemoveAlphaChannel#1 (medium): an exact or arbitrary-precision background. */
static void test_s14_removealpha_exact(void) {
    audit_run("auditT = SetAlphaChannel[Image[{{1.0, 1.0}, {1.0, 1.0}}, \"Real\"], 0.5]");
    audit_true("RemoveAlphaChannel#1",
               "Max[Abs[Flatten[ImageData[RemoveAlphaChannel[auditT, 1/2]]] - 0.75]] < 10^-12");
    audit_true("RemoveAlphaChannel#1",
               "Max[Abs[Flatten[ImageData[RemoveAlphaChannel[auditT, N[1/2, 30]]]] - 0.75]] < 10^-12");
    audit_run("ClearAll[auditT]");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Round#2", "Quiet[Check[Round[1, 2, 3], fired, Round::argt]]", "fired");
    audit_chk("Round#2", "Quiet[Round[1, 2, 3]]", "Round[1, 2, 3]");
    audit_chk("AlphaChannel#1", "Quiet[Check[AlphaChannel[], fired, AlphaChannel::argx]]", "fired");
    audit_chk("AlphaChannel#1", "Quiet[AlphaChannel[]]", "AlphaChannel[]");
    audit_chk("Erosion#3", "Quiet[Check[Erosion[Image[{{0.5, 0.2}}]], fired, Erosion::argr]]", "fired");
    audit_chk("FindThreshold#4", "Quiet[Check[FindThreshold[], fired, FindThreshold::argx]]", "fired");
    audit_chk("FindThreshold#4", "Quiet[FindThreshold[]]", "FindThreshold[]");
    audit_chk("RandomImage#6", "Quiet[Check[RandomImage[1, {2, 3}, Foo -> 1], fired, RandomImage::optx]]", "fired");
    audit_chk("RandomImage#6", "Quiet[RandomImage[1, {2, 3}, Foo -> 1]]", "RandomImage[1, {2, 3}, Foo -> 1]");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* Round#3: Round[Indeterminate] */
    audit_true("Round#3", "Round[Indeterminate] === Indeterminate");
    audit_true("Round#3", "Round[Indeterminate, 2] === Indeterminate");
}

/* ================================================================== */
/* Wave 4 G3: image processing                                        */
/* ================================================================== */

/* Round of an exact complex numeric rounds both parts (src/piecewise.c). */
static void test_wave4_round_complex(void) {
    audit_true("Round#1", "Round[Sqrt[2] + I Pi] === 1 + 3 I");
    audit_true("Round#1", "Round[I Sqrt[2]] === I");
    audit_true("Round#1", "Round[E + I] === 3 + I");
    audit_true("Round#1", "Floor[E + I Pi] === 2 + 3 I && Ceiling[E - I Pi] === 3 - 3 I");
    audit_true("Round#1", "Round[N[Sqrt[2] + I Pi]] === 1 + 3 I");
}

/* Exact numeric parameters (1/2, Pi/4, {{Pi}}) read as their values
 * (image_real_value, src/image.c). */
static void test_wave4_exact_parameters(void) {
    audit_run("img = Image[{{{1.,0.,0.},{0.,1.,0.}}}, \"Real\"]");
    audit_true("ColorReplace#2", "ImageData[ColorReplace[img, Red -> Blue, 1/10]] == {{{0., 0., 1.}, {0., 1., 0.}}}");
    audit_true("ColorReplace#2", "ImageData[ColorReplace[img, RGBColor[1/2,0,0] -> Blue, 0.6]] == {{{0., 0., 1.}, {0., 1., 0.}}}");
    audit_true("ImageCompose#1", "ImageData[ImageCompose[Image[{{0.5}}], {Image[{{1.}}], 1/2}]] == {{0.75}}");
    audit_true("ImageConvolve#1", "ImageData[ImageConvolve[Image[{{0.,1.,0.}}], {{Pi}}]] == N[{{0, Pi, 0}}]");
    audit_run("b = Image[{{0.1, 0.2}, {0.3, 0.4}}, \"Real\"]");
    audit_true("SetAlphaChannel#1", "ImageData[SetAlphaChannel[b, 1/2]] == {{{0.1, 0.5}, {0.2, 0.5}}, {{0.3, 0.5}, {0.4, 0.5}}}");
    audit_true("SetAlphaChannel#1", "ImageData[SetAlphaChannel[b, Pi/4]][[1, 1, 2]] == N[Pi/4]");
}

/* Index permutations keep the pixel type. */
static void test_wave4_type_preserved(void) {
    audit_chk("ImageCrop#1", "ImageType[ImageCrop[Image[{{10, 20}, {30, 40}}, \"Byte\"], {1, 1}]]", "\"Byte\"");
    audit_chk("ImageCrop#1", "ImageType[ImageCrop[Image[{{0, 1}, {1, 0}}, \"Bit\"], {1, 1}]]", "\"Bit\"");
    audit_chk("ImageCrop#1", "ImageType[ImageCrop[Image[{{0, 1, 0}, {1, 1, 1}, {0, 1, 0}}, \"Bit\"]]]", "\"Bit\"");
    audit_run("bi = Image[{{0, 128}, {255, 7}}, \"Byte\"]");
    audit_chk("ImageReflect#1", "ImageType[ImageReflect[bi]]", "\"Byte\"");
    audit_true("ImageReflect#1", "ImageData[ImageReflect[bi], \"Byte\"] == {{255, 7}, {0, 128}}");
    audit_chk("ColorConvert#1", "ImageType[ColorConvert[Image[{{{255, 0, 0}}}, \"Byte\"], \"Grayscale\"]]", "\"Byte\"");
    audit_chk("ColorConvert#1", "ImageType[ColorConvert[Image[{{3, 200}}, \"Byte\"], \"Grayscale\"]]", "\"Byte\"");
}

static void test_wave4_reflect_pad_rotate(void) {
    audit_run("img = Image[{{1., 2., 3.}, {4., 5., 6.}}]");
    audit_true("ImageReflect#2", "ImageData[ImageReflect[img, Top -> Bottom]] == ImageData[ImageReflect[img]] == {{4., 5., 6.}, {1., 2., 3.}}");
    audit_true("ImageReflect#2", "ImageData[ImageReflect[img, Left -> Top]] == {{1., 4.}, {2., 5.}, {3., 6.}}");
    audit_true("ImageReflect#2", "ImageData[ImageReflect[img, Left -> Right]] == {{3., 2., 1.}, {6., 5., 4.}}");
    audit_true("ImageReflect#2", "ImageData[ImageReflect[img, Left -> Left]] == ImageData[img]");
    audit_run("i = Image[{{0.1, 0.2, 0.3}, {0.4, 0.5, 0.6}}]");
    audit_true("ImagePad#1", "ImageData[ImagePad[i, 1, \"Periodic\"]] == "
                             "{{0.6, 0.4, 0.5, 0.6, 0.4}, {0.3, 0.1, 0.2, 0.3, 0.1}, {0.6, 0.4, 0.5, 0.6, 0.4}, {0.3, 0.1, 0.2, 0.3, 0.1}}");
    audit_true("ImagePad#1", "ImageData[ImagePad[i, 1, \"Reversed\"]] == "
                             "{{0.1, 0.1, 0.2, 0.3, 0.3}, {0.1, 0.1, 0.2, 0.3, 0.3}, {0.4, 0.4, 0.5, 0.6, 0.6}, {0.4, 0.4, 0.5, 0.6, 0.6}}");
    audit_run("c = Image[{{{1., 0., 0.}}}]");
    audit_true("ImagePad#2", "ImageData[ImagePad[c, 1, {0., 1., 0.}]] == "
                             "{{{0., 1., 0.}, {0., 1., 0.}, {0., 1., 0.}}, {{0., 1., 0.}, {1., 0., 0.}, {0., 1., 0.}}, "
                             "{{0., 1., 0.}, {0., 1., 0.}, {0., 1., 0.}}}");
    audit_true("ImagePad#2", "ImageData[ImagePad[c, 1, Blue]][[1, 1]] == {0., 0., 1.} && ImageData[ImagePad[c, 1, Blue]][[2, 2]] == {1., 0., 0.}");
    audit_run("i3 = Image[{{1., 2., 3.}, {4., 5., 6.}}]");
    audit_true("ImageRotate#1", "ImageDimensions[ImageRotate[i3, 1.*^20]] == {3, 2}");
    audit_true("ImageRotate#1", "ImageData[ImageRotate[i3, 1.*^20]] != {{4., 1.}, {5., 2.}, {6., 3.}}");
}

/* The filters take plain numeric data (builtin_filter_data, src/imagefilter.c). */
static void test_wave4_filter_data(void) {
    audit_true("Erosion#1", "Erosion[{{1,2,3},{4,5,6},{7,8,9}}, 1] === {{1, 1, 2}, {1, 1, 2}, {4, 4, 5}}");
    audit_true("Closing#2", "Closing[{{0.1,0.9},{0.5,0.2}}, 1] == {{0.9, 0.9}, {0.9, 0.9}}");
    audit_true("GaussianFilter#1", "Dimensions[GaussianFilter[{1., 2., 3., 4., 5.}, 1]] == {5} && "
                                   "VectorQ[GaussianFilter[{1., 2., 3., 4., 5.}, 1], NumericQ]");
    audit_true("GaussianFilter#1", "GaussianFilter[{1., 2., 3., 4., 5.}, 1][[3]] == 3.");
    audit_true("GaussianFilter#1", "Dimensions[GaussianFilter[{{1., 2.}, {3., 4.}}, 1]] == {2, 2}");
    audit_true("MeanFilter#1", "MeanFilter[{1,2,3,4,5}, 1] === {4/3, 2, 3, 4, 14/3}");
    audit_true("MedianFilter#1", "MedianFilter[{1, 5, 2, 8, 3}, 1][[2 ;; 4]] === {2, 5, 3}");
    audit_true("MedianFilter#1", "MedianFilter[{{1, 2, 9}, {3, 4, 5}, {6, 7, 8}}, 1][[2, 2]] === 5");
    audit_true("DerivativeFilter#1", "DerivativeFilter[{{0.,1.,2.},{0.,1.,2.},{0.,1.,2.}}, {0,1}][[All, 2]] == {1., 1., 1.}");
    audit_true("GradientFilter#2", "Dimensions[GradientFilter[{{0,0,1,1},{0,0,1,1},{0,0,1,1}}, 1]] == {3, 4}");
    audit_true("GradientFilter#2", "Dimensions[GradientFilter[{{0,0,1,1},{0,0,1,1}}]] == {2, 4}");
    audit_true("MorphologicalComponents#1", "MorphologicalComponents[{{1,0},{0,1}}] === {{1, 0}, {0, 1}}");
    audit_true("MorphologicalComponents#1",
               "MorphologicalComponents[{{1,0,1},{1,0,1},{1,1,1}}, CornerNeighbors->False] === {{1, 0, 1}, {1, 0, 1}, {1, 1, 1}}");
    /* an Image argument still takes the image operator */
    audit_true("Erosion#1", "ImageQ[Erosion[Image[{{0.1, 0.9}, {0.5, 0.2}}], 1]]");
}

/* Out[] is Out[-1]. */
static void test_wave4_out_noarg(void) {
    audit_true("Out#1", "Unprotect[Out]; Out[2] = 7; Out[3] = 8; Protect[Out]; $Line = 4; Out[] === Out[-1] === 8");
}

static void test_wave4_threshold_forms(void) {
    audit_true("BoxMatrix#1", "BoxMatrix[{1, 2}] === ConstantArray[1, {3, 5}]");
    audit_true("BoxMatrix#1", "BoxMatrix[1] === ConstantArray[1, {3, 3}]");
    audit_true("Binarize#1", "b1 = Binarize[Image[{{0.5,0.5}}]]; ImageQ[b1] && ImageType[b1] === \"Bit\" && ImageDimensions[b1] == {2, 1}");
    audit_true("Binarize#2", "ImageData[Binarize[Image[{{0.1,0.9}}], {0.2,0.95}]] == {{0, 1}}");
    audit_true("Binarize#2", "ImageData[Binarize[Image[{{0.1,0.5,0.9}}], {0.2,0.6}]] == {{0, 1, 0}}");
    audit_true("FindThreshold#1", "FindThreshold[Image[{{0.2,0.2},{0.8,0.8}}], Method -> \"Mean\"] == 0.5");
    audit_true("FindThreshold#1", "FindThreshold[Image[{{0.2,0.2},{0.8,0.8}}], Method -> \"Cluster\"] == "
                                  "FindThreshold[Image[{{0.2,0.2},{0.8,0.8}}]]");
    audit_true("EdgeDetect#1", "ImageData[EdgeDetect[Image[ConstantArray[0.5,{3,3}]],0,0]] == ConstantArray[0, {3, 3}]");
    audit_true("ImageAdjust#1", "ImageData[ImageAdjust[Image[{{0.2,0.6}}], 0.5]] == ImageData[ImageAdjust[Image[{{0.2,0.6}}], {0.5, 0}]]");
    audit_true("ImageAdjust#1", "ImageQ[ImageAdjust[Image[{{0.2,0.6}}], 0.5]]");
}

static void test_wave4_resize_random_compose(void) {
    audit_run("img = Image[Table[N[(i+j)/10], {i, 3}, {j, 6}]]");
    audit_true("ImageResize#1", "ImageDimensions[ImageResize[img, {Automatic, 6}]] == {12, 6}");
    audit_true("ImageResize#1", "ImageDimensions[ImageResize[img, {12, Automatic}]] == {12, 6}");
    audit_true("ImageResize#2", "ImageDimensions[ImageResize[img, Scaled[2]]] == {12, 6}");
    audit_true("ImageResize#2", "ImageDimensions[ImageResize[img, Scaled[1/3]]] == {2, 1}");
    audit_true("RandomImage#2", "r2 = RandomImage[{2,3},{4,4}]; Head[r2] === Image && ImageDimensions[r2] == {4, 4} && "
                                "2 <= Min[ImageData[r2]] && Max[ImageData[r2]] <= 3");
    audit_run("t = SetAlphaChannel[Image[{{1.0, 1.0}, {1.0, 1.0}}, \"Real\"], 0.5]");
    audit_true("RemoveAlphaChannel#2", "ImageData[RemoveAlphaChannel[t, GrayLevel[0.2]]] == {{0.6, 0.6}, {0.6, 0.6}}");
    audit_true("RemoveAlphaChannel#2", "ImageData[RemoveAlphaChannel[t, Red]][[1, 1]] == {1., 0.5, 0.5}");
    audit_true("ImageCompose#2",
               "bb=SetAlphaChannel[Image[{{{0.,0.,0.}}}], 0.]; ImageData[ImageCompose[bb, SetAlphaChannel[Image[{{{1.,0.,0.}}}], 0.5]]] == {{{1., 0., 0., 0.5}}}");
    /* an opaque base is unchanged by the over operator */
    audit_true("ImageCompose#2", "ImageData[ImageCompose[Image[{{0.5}}], {Image[{{1.}}], 0.5}]] == {{0.75}}");
}

static void test_wave4_levels_colorreplace(void) {
    audit_true("ImageLevels#1", "ImageLevels[Image[{{{1,0,0},{0,1,0}}}]] == {{0., {1, 1, 2}}, {1., {1, 1, 0}}}");
    audit_true("ImageLevels#1", "ImageLevels[Image[{{0, 1}, {1, 0}}]] == {{0., 2}, {1., 2}}");
    audit_run("ga = Image[{{{0.5, 1.}, {0.2, 0.5}}}, \"Real\"]");
    audit_true("ColorReplace#1", "ImageChannels[ColorReplace[ga, GrayLevel[0.5] -> RGBColor[1,0,0]]] == 4");
    audit_true("ColorReplace#1", "ImageData[ColorReplace[ga, GrayLevel[0.5] -> RGBColor[1,0,0]]] == {{{1., 0., 0., 1.}, {0.2, 0.2, 0.2, 0.5}}}");
    audit_true("ColorReplace#1", "ImageData[ColorReplace[ga, GrayLevel[0.5] -> GrayLevel[0.1]]] == {{{0.1, 1.}, {0.2, 0.5}}}");
}

/* Dilation stamps the element; erosion reads x + b; opening <= f <= closing. */
static void test_wave4_dilation_element(void) {
    audit_true("Dilation#1", "ImageData[Dilation[Image[{{0.,0.,0.},{0.,1.,0.},{0.,0.,0.}}], {{1,0,0},{0,1,0},{0,0,0}}]] == "
                             "{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}");
    audit_run("f = Image[{{0.,0.,0.,0.,0.},{0.,1.,1.,1.,0.},{0.,1.,1.,1.,0.},{0.,1.,1.,1.,0.},{0.,0.,0.,0.,0.}}]; k = {{1,1,0},{1,1,0},{0,0,0}}");
    audit_true("Dilation#1", "ImageData[Dilation[Erosion[f, k], k]] == ImageData[f]");
    audit_true("Dilation#1", "And @@ Flatten[MapThread[LessEqual, {ImageData[Opening[f, k]], ImageData[f]}, 2]]");
    audit_true("Dilation#1", "And @@ Flatten[MapThread[GreaterEqual, {ImageData[Closing[f, k]], ImageData[f]}, 2]]");
    audit_true("Dilation#1", "SeedRandom[3]; Module[{g = RandomImage[1, {7, 6}], e = {{1, 1, 0}, {0, 1, 1}, {0, 0, 0}}},"
                             " Max[ImageData[Opening[g, e]] - ImageData[g]] <= 0 && Min[ImageData[Closing[g, e]] - ImageData[g]] >= 0"
                             " && ImageData[Opening[Opening[g, e], e]] == ImageData[Opening[g, e]]]");
}

/* Thinning keeps a pixel of a 2x2 block; Pruning removes a spur up to the loop it meets. */
static void test_wave4_thinning_pruning(void) {
    audit_run("sq = Image[{{1,1,0,0,0},{1,1,0,0,0},{0,0,0,0,0},{0,0,0,1,0},{0,0,0,0,0}}]");
    audit_true("Thinning#1", "Max[MorphologicalComponents[Thinning[sq]]] == 2");
    audit_true("Thinning#1", "Total[ImageData[Thinning[Image[{{1,1},{1,1}}]]], 2] >= 1");
    audit_true("Thinning#1", "ImageData[Thinning[Image[ConstantArray[0, {3, 3}]]]] == ConstantArray[0, {3, 3}]");
    audit_run("sp = Image[Table[If[(i==2||i==6)&&2<=j<=6 || (j==2||j==6)&&2<=i<=6 || (i==4 && 7<=j<=9),1,0],{i,1,8},{j,1,11}]]");
    audit_true("Pruning#1", "Round[ImageData[Pruning[sp, 10]]][[4]] == {0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0}");
    /* the loop itself is untouched */
    audit_true("Pruning#1", "Total[Round[ImageData[Pruning[sp, 10]]], 2] == 16");
    /* a straight segment loses one pixel from each end per pass */
    audit_true("Pruning#1", "Total[Round[ImageData[Pruning[Image[{{0,0,0,0,0,0,0},{0,1,1,1,1,1,0},{0,0,0,0,0,0,0}}], 1]]], 2] == 3");
}

/* Image presentation options are accepted; Image3D stacks 2-D images. */
static void test_wave4_image_options_stack(void) {
    audit_true("Image#1", "ImageQ[Image[{{0, 1}}, ColorSpace -> \"RGB\"]]");
    audit_true("Image#1", "ImageQ[Image[{{0, 1}}, \"Byte\", ColorSpace -> \"Grayscale\"]] && "
                          "ImageType[Image[{{0, 1}}, \"Byte\", ColorSpace -> \"Grayscale\"]] === \"Byte\"");
    audit_true("Image#1", "ImageQ[Image[{{0.5}}, Interleaving -> True]]");
    audit_true("Image#1", "Head[Image[{{0.5}}, NoSuchImageOption -> 1]] === Image && !ImageQ[Image[{{0.5}}, NoSuchImageOption -> 1]]");
    audit_true("ImageDimensions#1", "ImageDimensions[Image[{{0.,1.}}, ColorSpace->\"Grayscale\"]] == {2, 1}");
    audit_true("ImageDimensions#1", "ImageDimensions[Image[Table[0.,{4},{6}], \"Real\", ColorSpace->\"Grayscale\"]] == {6, 4}");
    audit_true("Image3D#1", "Image3D[{Image[{{0.,1.}}],Image[{{1.,0.}}]}] === Image3D[{{{0., 1.}}, {{1., 0.}}}, \"Real\"]");
    audit_true("Image3D#1", "ImageDimensions[Image3D[{Image[{{0.,1.}}],Image[{{1.,0.}}]}]] == {2, 1, 2}");
}

static void test_wave4_gaussian_derivative(void) {
    audit_true("GaussianMatrix#1", "Dimensions[GaussianMatrix[1, {0, 1}]] == {3, 3}");
    audit_true("GaussianMatrix#1", "GaussianMatrix[1, {0, 1}] == -Reverse /@ GaussianMatrix[1, {0, 1}]");
    audit_true("GaussianMatrix#1", "GaussianMatrix[1, {0, 1}][[2, 2]] == 0 && GaussianMatrix[1, {0, 1}][[2, 1]] > 0");
    audit_true("GaussianMatrix#1", "GaussianMatrix[1, {0, 0}] == GaussianMatrix[1]");
    audit_true("GaussianMatrix#1", "GaussianMatrix[2, {1, 0}] == Transpose[GaussianMatrix[2, {0, 1}]]");
}

static void test_wave4_palette_padding_metrics(void) {
    audit_true("ColorQuantize#2", "ImageData[ColorQuantize[Image[{{0., 1.}}, \"Real\"], {Black, White}]] == {{0., 1.}}");
    audit_true("ColorQuantize#2", "ImageData[ColorQuantize[Image[{{0.2, 0.7}}, \"Real\"], {Black, White}]] == {{0., 1.}}");
    audit_true("ColorQuantize#2", "ImageData[ColorQuantize[Image[{{{0.9, 0.1, 0.}}}], {Red, Blue}]] == {{{1., 0., 0.}}}");
    audit_true("ImageConvolve#2", "ImageData[ImageConvolve[Image[{{0.,1.,0.}}], {{1,1,1}}, Padding -> 0]] == {{1., 1., 1.}}");
    audit_true("ImageConvolve#2", "ImageData[ImageConvolve[Image[{{0.,1.,0.}}], {{1, 2, 3}}, Padding -> \"Fixed\"]] == {{1., 2., 3.}}");
    audit_true("ImageConvolve#2", "ImageData[ImageConvolve[Image[{{1.,1.,1.}}], {{1,1,1}}, Padding -> 0]] == {{2., 3., 2.}}");
    audit_true("DistanceTransform#2",
               "ImageData[DistanceTransform[Image[{{0,1,1}}, \"Bit\"], DistanceFunction -> ManhattanDistance]] == {{0., 1., 2.}}");
    audit_true("DistanceTransform#2",
               "ImageData[DistanceTransform[Image[{{0,1,1},{1,1,1},{1,1,1}}, \"Bit\"], DistanceFunction -> ChessboardDistance]] == "
               "{{0., 1., 2.}, {1., 1., 2.}, {2., 2., 2.}}");
    audit_true("DistanceTransform#2",
               "ImageData[DistanceTransform[Image[{{0,1,1},{1,1,1},{1,1,1}}, \"Bit\"], DistanceFunction -> ManhattanDistance]] == "
               "{{0., 1., 2.}, {1., 2., 3.}, {2., 3., 4.}}");
    audit_true("DistanceTransform#2",
               "ImageData[DistanceTransform[Image[{{0,1,1}}, \"Bit\"], DistanceFunction -> EuclideanDistance]] == "
               "ImageData[DistanceTransform[Image[{{0,1,1}}, \"Bit\"]]]");
}

static void test_wave4_histogram_match_correlate(void) {
    audit_true("HistogramTransform#1", "ImageData[HistogramTransform[Image[{{0.2, 0.8}}], Image[{{0.5, 0.6}}]]] == {{0.5, 0.6}}");
    audit_true("HistogramTransform#1", "ImageQ[HistogramTransform[Image[{{0.2, 0.8}}], NormalDistribution[0.5, 0.1]]]");
    audit_true("HistogramTransform#1", "Module[{d = ImageData[HistogramTransform[Image[{{0.2, 0.8}}], NormalDistribution[0.5, 0.1]]]},"
                                       " d[[1, 1]] < 0.5 < d[[1, 2]] && Abs[d[[1, 1]] + d[[1, 2]] - 1] < 10^-9]");
    audit_true("ImageCorrelate#2", "ImageData[ImageCorrelate[Image[{{0.,1.}}],{{1}},EuclideanDistance]] == {{1., 0.}}");
    audit_true("ImageCorrelate#2", "ImageData[ImageCorrelate[Image[{{0.,1.,0.5}}],{{1, 0}},ManhattanDistance]] == {{1., 2., 0.5}}");
}

int main(void) {
    audit_init();
    TEST(test_random_image_huge_size);
    TEST(test_opening_closing_asymmetric);
    TEST(test_color_quantize_value_boundary);
    TEST(test_image_correlate_ncc_even);
    TEST(test_gradient_filter_radius);
    TEST(test_s14_removealpha_exact);
    TEST(test_wave3_arity_options);
    TEST(test_s11_special_values);
    TEST(test_wave4_round_complex);
    TEST(test_wave4_exact_parameters);
    TEST(test_wave4_type_preserved);
    TEST(test_wave4_reflect_pad_rotate);
    TEST(test_wave4_filter_data);
    TEST(test_wave4_out_noarg);
    TEST(test_wave4_threshold_forms);
    TEST(test_wave4_resize_random_compose);
    TEST(test_wave4_levels_colorreplace);
    TEST(test_wave4_dilation_element);
    TEST(test_wave4_thinning_pruning);
    TEST(test_wave4_image_options_stack);
    TEST(test_wave4_gaussian_derivative);
    TEST(test_wave4_palette_padding_metrics);
    TEST(test_wave4_histogram_match_correlate);
    return audit_finish("audit_image_processing");
}
