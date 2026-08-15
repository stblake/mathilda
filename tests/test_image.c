/* test_image.c -- the Image representation and its accessors.
 *
 * Three of these assertions are the ones that would actually catch a bug, and they are the
 * reason the others are here at all:
 *
 *   - THE TRANSPOSITION. ImageDimensions is {width, height} while ImageData is height x width.
 *     Every row that touches it uses a NON-SQUARE image, because a square one cannot tell the
 *     two apart and would pass with the axes swapped;
 *   - EXACT BYTE SCALING. 255 must come back as exactly 1.0 and 0 as exactly 0.0, so these are
 *     equalities rather than tolerances -- the values are chosen to be exactly representable;
 *   - THE FIXED POINT. Image[data, type] is canonical and must evaluate to itself, or the
 *     evaluator would rewrite it forever.
 */
#include <stdio.h>
#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

/* Two rows of three: non-square, so width and height cannot be confused. */
#define I23 "Image[{{0., 0.1, 0.2}, {0.3, 0.4, 0.5}}]"

static void test_dimensions_are_transposed_from_data(void) {
    /* THE row. Two rows of three pixels is width 3, height 2 -- and ImageData is the other way
     * round. Getting this backwards is the classic image bug, and a square test image would
     * hide it completely. */
    assert_eval_eq("ImageDimensions[" I23 "]", "{3, 2}", 0);
    assert_eval_eq("Dimensions[ImageData[" I23 "]]", "{2, 3}", 0);
    /* And the pixel at data[[y, x]] is the one the dimensions imply: row 2, column 3 is 0.5. */
    assert_eval_eq("Part[ImageData[" I23 "], 2, 3]", "0.5", 0);
    /* A colour image keeps the same convention with channels last. */
    assert_eval_eq("ImageDimensions[Image[{{{1.,0.,0.}, {0.,1.,0.}}}]]", "{2, 1}", 0);
    assert_eval_eq("Dimensions[ImageData[Image[{{{1.,0.,0.}, {0.,1.,0.}}}]]]", "{1, 2, 3}", 0);
}

static void test_byte_scaling_is_exact(void) {
    /* Equalities, not tolerances: 0, 255 and the ratio 128/255 are all exactly what double
     * arithmetic gives, so a wrong divisor (254, or 256) fails outright rather than by a
     * hair. */
    assert_eval_eq("ImageData[Image[{{0, 255}}]]", "{{0.0, 1.0}}", 0);
    assert_eval_eq("Last[First[ImageData[Image[{{0, 255}}]]]] == 1.0", "True", 0);
    assert_eval_eq("First[First[ImageData[Image[{{0, 255}}]]]] == 0.0", "True", 0);
    assert_eval_eq("Part[ImageData[Image[{{0, 128, 255}}]], 1, 2] == 128./255.", "True", 0);
    /* A Bit image is already unit-scaled, so scaling must be the identity rather than /255. */
    assert_eval_eq("ImageData[Image[{{0, 1}, {1, 0}}]]",
                   "{{0.0, 1.0}, {1.0, 0.0}}", 0);
    /* Real data passes through untouched, including values outside [0,1] -- storing faithfully
     * beats clamping silently, because a clamp would destroy data the caller may want back. */
    assert_eval_eq("ImageData[Image[{{-0.5, 1.5}}]]", "{{-0.5, 1.5}}", 0);
}

static void test_stored_values_round_trip(void) {
    /* Asking for the image's own type gives back exactly what went in. This is what makes the
     * scaling a property of ImageData alone rather than something baked into construction. */
    assert_eval_eq("ImageData[Image[{{0, 128, 255}}], \"Byte\"] === {{0, 128, 255}}",
                   "True", 0);
    assert_eval_eq("ImageData[Image[{{0, 1}}], \"Bit\"] === {{0, 1}}", "True", 0);
    /* Asking for a DIFFERENT type declines: conversion has its own rounding decisions and is
     * not something an accessor should do silently. */
    assert_eval_eq("Head[ImageData[Image[{{0, 255}}], \"Real\"]]", "ImageData", 0);
    assert_eval_eq("Head[ImageData[Image[{{0., 0.5}}], \"Byte\"]]", "ImageData", 0);
}

static void test_type_is_inferred_from_the_values(void) {
    /* Inference looks at the values and nothing else, so it is predictable from the data. The
     * integer/real distinction is the interesting one: the same numbers written differently
     * give different types, and a caller can rely on that. */
    assert_eval_eq("ImageType[Image[{{0, 1}, {1, 0}}]]", "\"Bit\"", 0);
    assert_eval_eq("ImageType[Image[{{0., 1.}, {1., 0.}}]]", "\"Real\"", 0);
    assert_eval_eq("ImageType[Image[{{0, 255}, {128, 64}}]]", "\"Byte\"", 0);
    assert_eval_eq("ImageType[Image[{{0, 2}}]]", "\"Byte\"", 0);
    assert_eval_eq("ImageType[Image[{{0, 256}}]]", "\"Real\"", 0);
    assert_eval_eq("ImageType[Image[{{-1, 1}}]]", "\"Real\"", 0);
    /* A stated type is honoured where the data fits it, and refused where it does not --
     * reinterpreting 300 as a byte would corrupt every later scaling. */
    assert_eval_eq("ImageType[Image[{{0, 1}}, \"Byte\"]]", "\"Byte\"", 0);
    assert_eval_eq("Head[Image[{{0, 300}}, \"Byte\"]]", "Image", 0);
    assert_eval_eq("Head[Image[{{0, 5}}, \"Bit\"]]", "Image", 0);
    assert_eval_eq("Head[Image[{{0, 1}}, \"Nonsense\"]]", "Image", 0);
}

static void test_canonical_form_is_a_fixed_point(void) {
    /* Image[data, type] must evaluate to ITSELF. Normalising unconditionally would rewrite it
     * on every pass and the evaluator would never settle -- this row is what fails if the
     * already-canonical early return is ever removed. */
    assert_eval_eq("Image[{{0, 1}}, \"Bit\"]", "Image[{{0, 1}}, \"Bit\"]", 0);
    assert_eval_eq("Image[Image[{{0, 1}}][[1]], Image[{{0, 1}}][[2]]] === Image[{{0, 1}}]",
                   "True", 0);
    /* Normalisation added the type, so the one-argument form is gone after evaluation. */
    assert_eval_eq("Length[Image[{{0, 1}}]]", "2", 0);
}

static void test_channels(void) {
    assert_eval_eq("ImageChannels[" I23 "]", "1", 0);
    assert_eval_eq("ImageChannels[Image[{{{1., 0., 0.}}}]]", "3", 0);
    assert_eval_eq("ImageChannels[Image[{{{1., 0., 0., 0.5}}}]]", "4", 0);
    /* A pixel list of a different length in one place is not a channel count, so it declines
     * rather than reporting the first pixel's length as if it were the image's. */
    assert_eval_eq("Head[Image[{{{1., 0., 0.}, {0., 1.}}}]]", "Image", 0);
}

static void test_malformed_input_declines(void) {
    /* Ragged data would be indexed as rectangular by every filter downstream, so a clear
     * refusal here beats an out-of-bounds read somewhere far away. */
    assert_eval_eq("Head[Image[{{1, 2}, {3}}]]", "Image", 0);
    assert_eval_eq("Head[Image[{}]]", "Image", 0);
    assert_eval_eq("Head[Image[{{}}]]", "Image", 0);
    assert_eval_eq("Head[Image[\"hello\"]]", "Image", 0);
    assert_eval_eq("Head[Image[{{1, x}}]]", "Image", 0);
    assert_eval_eq("Head[Image[{1, 2, 3}]]", "Image", 0);            /* rank 1 is not an image */
    assert_eval_eq("Head[Image[{{1, 2}}, \"Byte\", 3]]", "Image", 0);
    /* Complex pixels are not brightnesses. */
    assert_eval_eq("Head[Image[{{1, I}}]]", "Image", 0);
    /* The accessors decline on a non-image rather than inventing an answer. */
    assert_eval_eq("Head[ImageDimensions[{{1, 2}}]]", "ImageDimensions", 0);
    assert_eval_eq("Head[ImageData[{{1, 2}}]]", "ImageData", 0);
    assert_eval_eq("Head[ImageType[42]]", "ImageType", 0);
}

static void test_imageq_and_attributes(void) {
    /* Malformed input to Image stays unevaluated, so ImageQ is the way validity is tested --
     * Head alone cannot distinguish a valid image from a refused one, since both are Image. */
    assert_eval_eq("ImageQ[" I23 "]", "True", 0);
    assert_eval_eq("ImageQ[Image[{{0, 1}}]]", "True", 0);
    assert_eval_eq("ImageQ[{{1, 2}}]", "False", 0);
    assert_eval_eq("ImageQ[Image[{{1, 2}, {3}}]]", "False", 0);
    assert_eval_eq("ImageQ[Image[{{1, 2}}, \"Nonsense\"]]", "False", 0);
    /* A HAND-TYPED canonical form is not trusted blindly. image_info reads the shape without
     * re-checking every pixel -- a benchmark showed the full walk made ImageDimensions cost
     * 0.59 ms on a 512x512 image -- but it still verifies rectangularity in O(height), because
     * every filter downstream indexes the array as rectangular. This row is what that cheap
     * check exists for, and what fails if it is ever dropped as redundant. */
    assert_eval_eq("ImageQ[Image[{{1, 2}, {3}}, \"Bit\"]]", "False", 0);
    assert_eval_eq("Head[ImageDimensions[Image[{{1, 2}, {3}}, \"Bit\"]]]",
                   "ImageDimensions", 0);
    assert_eval_eq("And @@ Map[MemberQ[Attributes[#], Protected] &,"
                   " {Image, ImageQ, ImageDimensions, ImageChannels, ImageType, ImageData}]",
                   "True", 0);
}

static void test_convolution_reflects_the_kernel(void) {
    /* THE discriminating row, and the only one that can tell convolution from correlation.
     * They agree exactly on every symmetric kernel -- a Gaussian, a box -- so a smoothing test
     * passes either way and the mistake goes unnoticed. With a delta at the centre of a 1x3
     * image and the asymmetric kernel {{1,2,3}}, convolution gives {1,2,3} and correlation gives
     * {3,2,1}: exact integers, no tolerance, and no way to pass by accident. */
    assert_eval_eq("ImageData[ImageConvolve[Image[{{0., 1., 0.}}], {{1, 2, 3}}]]",
                   "{{1.0, 2.0, 3.0}}", 0);
    /* Vertically too, so a transposed index cannot hide in the row case. */
    assert_eval_eq("ImageData[ImageConvolve[Image[{{0.}, {1.}, {0.}}], {{1}, {2}, {3}}]]",
                   "{{1.0}, {2.0}, {3.0}}", 0);
}

static void test_identity_kernel_is_the_identity(void) {
    /* A 1x1 kernel of 1 must return the pixels EXACTLY -- not close, identical. Any stray
     * scaling, any off-by-one in the centre calculation, and this fails. */
    assert_eval_eq("ImageData[ImageConvolve[" I23 ", {{1}}]] === ImageData[" I23 "]",
                   "True", 0);
    /* And dimensions survive, which is what fails if the output is built transposed. */
    assert_eval_eq("ImageDimensions[ImageConvolve[" I23 ", GaussianMatrix[2]]]", "{3, 2}", 0);
}

static void test_constant_image_survives_a_normalised_kernel(void) {
    /* An exact property that specifically tests the BORDER. A kernel summing to 1 over a
     * constant image must give that same constant everywhere -- and at the edges that is only
     * true because out-of-range reads clamp to the edge pixel. Zero padding would darken the
     * border, which looks exactly like a real vignetting bug, so this row is what stands between
     * that and going unnoticed. Chop, because the sum of 25 scaled doubles need not be bit-exact. */
    assert_eval_eq("Chop[Max[Abs[Flatten[ImageData[ImageConvolve["
                   "Image[Table[0.25, {6}, {6}]], GaussianMatrix[2]]] - 0.25]]]]", "0", 0);
    /* A box kernel is NOT normalised, so the same image comes back nine times brighter. That is
     * Mathematica's definition and this row pins that it was not "helpfully" rescaled. */
    assert_eval_eq("Chop[Part[ImageData[ImageConvolve[Image[Table[0.1, {4}, {4}]],"
                   " BoxMatrix[1]]], 2, 2] - 0.9]", "0", 0);
}

static void test_kernel_constructors(void) {
    /* A Gaussian matrix must sum to exactly 1, because it is normalised by the REALISED sum
     * rather than the analytic 2 pi sigma^2 -- the analytic constant is right only for an
     * infinite kernel, and on a truncated one it leaves the sum under 1, darkening an image a
     * little on every pass. */
    assert_eval_eq("Chop[Total[Flatten[GaussianMatrix[2]]] - 1.]", "0", 0);
    assert_eval_eq("Chop[Total[Flatten[GaussianMatrix[4]]] - 1.]", "0", 0);
    assert_eval_eq("Dimensions[GaussianMatrix[3]]", "{7, 7}", 0);
    /* Symmetric about the centre, and the centre is the largest entry. */
    assert_eval_eq("GaussianMatrix[1] === Transpose[GaussianMatrix[1]]", "True", 0);
    assert_eval_eq("Part[GaussianMatrix[2], 3, 3] == Max[Flatten[GaussianMatrix[2]]]",
                   "True", 0);
    /* r = 0 is a single tap of 1, which makes it the identity kernel. */
    assert_eval_eq("GaussianMatrix[0]", "{{1.0}}", 0);
    /* An explicit sigma is honoured: a larger sigma is flatter, so its centre is smaller. */
    assert_eval_eq("Part[GaussianMatrix[{2, 5.}], 3, 3] < Part[GaussianMatrix[{2, 0.5}], 3, 3]",
                   "True", 0);
    assert_eval_eq("BoxMatrix[0]", "{{1}}", 0);
    assert_eval_eq("BoxMatrix[2] === Table[1, {5}, {5}]", "True", 0);
    /* A fractional or negative radius has no matrix size, so it declines. */
    assert_eval_eq("Head[GaussianMatrix[1.5]]", "GaussianMatrix", 0);
    assert_eval_eq("Head[GaussianMatrix[-1]]", "GaussianMatrix", 0);
    assert_eval_eq("Head[BoxMatrix[1.5]]", "BoxMatrix", 0);
    assert_eval_eq("Head[GaussianMatrix[{2, -1.}]]", "GaussianMatrix", 0);
}

static void test_gaussianfilter_equals_imageconvolve(void) {
    /* Mathematica DOCUMENTS these as equal, so the identity is asserted rather than assumed.
     * GaussianFilter is implemented by building the same matrix and calling the same convolution
     * for exactly this reason -- two independent implementations of one identity is how the
     * identity quietly stops holding. */
    assert_eval_eq("ImageData[GaussianFilter[" I23 ", 1]] === "
                   "ImageData[ImageConvolve[" I23 ", GaussianMatrix[1]]]", "True", 0);
    assert_eval_eq("ImageData[GaussianFilter[" I23 ", 2]] === "
                   "ImageData[ImageConvolve[" I23 ", GaussianMatrix[2]]]", "True", 0);
}

static void test_filters_handle_colour_and_decline_junk(void) {
    /* Each channel is convolved independently, so the channel count survives. */
    assert_eval_eq("Module[{c = Image[{{{1.,0.,0.}, {0.,1.,0.}}, {{0.,0.,1.}, {1.,1.,0.}}}]},"
                   " {ImageChannels[GaussianFilter[c, 1]], ImageDimensions[GaussianFilter[c, 1]],"
                   "  ImageType[GaussianFilter[c, 1]]}]", "{3, {2, 2}, \"Real\"}", 0);
    /* The result is always Real, because a Gaussian of bytes is not a byte and rounding back
     * would discard precision nobody asked to lose. */
    assert_eval_eq("ImageType[ImageConvolve[Image[{{0, 255}, {255, 0}}], GaussianMatrix[1]]]",
                   "\"Real\"", 0);
    /* Declines: a non-image, a ragged kernel, a rank-1 kernel, a symbolic entry. */
    assert_eval_eq("Head[ImageConvolve[{{1, 2}}, {{1}}]]", "ImageConvolve", 0);
    assert_eval_eq("Head[ImageConvolve[" I23 ", {{1, 2}, {3}}]]", "ImageConvolve", 0);
    assert_eval_eq("Head[ImageConvolve[" I23 ", {1, 2, 3}]]", "ImageConvolve", 0);
    assert_eval_eq("Head[ImageConvolve[" I23 ", {{1, x}}]]", "ImageConvolve", 0);
    assert_eval_eq("Head[GaussianFilter[{{1, 2}}, 1]]", "GaussianFilter", 0);
    assert_eval_eq("And @@ Map[MemberQ[Attributes[#], Protected] &,"
                   " {ImageConvolve, GaussianMatrix, BoxMatrix, GaussianFilter}]", "True", 0);
}

#define BIMODAL "Image[{{0.2,0.2,0.8,0.8},{0.2,0.2,0.8,0.8}," \
                "{0.2,0.2,0.8,0.8},{0.2,0.2,0.8,0.8}}]"

static void test_otsu_splits_a_bimodal_image_perfectly(void) {
    /* The absolute property, and the analogue of "reproduces every training label": on data with
     * two clean clusters there is a ground truth, so Otsu must recover it EXACTLY -- every pixel
     * on the right side, not merely most. */
    assert_eval_eq("ImageData[Binarize[" BIMODAL "], \"Bit\"] === "
                   "{{0,0,1,1},{0,0,1,1},{0,0,1,1},{0,0,1,1}}", "True", 0);
    /* And the threshold lands strictly between the two clusters. */
    assert_eval_eq("0.2 < FindThreshold[" BIMODAL "] < 0.8", "True", 0);
    /* Result is a Bit image at the original dimensions. */
    assert_eval_eq("{ImageType[Binarize[" BIMODAL "]], ImageDimensions[Binarize[" BIMODAL "]]}",
                   "{\"Bit\", {4, 4}}", 0);
}

static void test_otsu_is_the_argmax_of_between_class_variance(void) {
    /* Verified against an INDEPENDENT computation of the objective, written from the pixel values
     * and sharing no code with the implementation: between-class variance is
     * w0 w1 (mu0 - mu1)^2, and the returned threshold must maximise it. This is much stronger
     * than comparing against another library's reported number, which would only pin a bin
     * convention -- skimage returns the winning bin's CENTRE and bins over [min, max], where this
     * returns the upper EDGE and bins over [0, 1], so the two legitimately differ by a bin or two
     * while agreeing on which split is best. */
    assert_eval_eq("Module[{r, px, sb, t0, best},"
                   " r = Image[Table[N[Mod[x*7 + y*13, 251]/251], {y, 48}, {x, 48}]];"
                   " px = Flatten[ImageData[r]];"
                   " sb[t_] := Module[{lo, hi, w0, w1},"
                   "   lo = Select[px, # <= t &]; hi = Select[px, # > t &];"
                   "   w0 = Length[lo]/Length[px]; w1 = Length[hi]/Length[px];"
                   "   If[w0 == 0 || w1 == 0, 0., w0 w1 (Mean[lo] - Mean[hi])^2]];"
                   " t0 = FindThreshold[r];"
                   " best = Max[Table[sb[t], {t, 0.02, 0.98, 0.02}]];"
                   " sb[t0] >= 0.99 best]", "True", 0);
}

static void test_binarize_boundary_is_strictly_above(void) {
    /* "Above" versus "at or above" differ on exactly the pixels a threshold was chosen to sit
     * between, so the boundary is pinned rather than left to chance. At the threshold: 0. */
    assert_eval_eq("ImageData[Binarize[Image[{{0.4, 0.5, 0.6}}], 0.5], \"Bit\"]",
                   "{{0, 0, 1}}", 0);
    /* A threshold below everything gives all 1s; above everything, all 0s. */
    assert_eval_eq("ImageData[Binarize[Image[{{0.4, 0.6}}], -1.], \"Bit\"]", "{{1, 1}}", 0);
    assert_eval_eq("ImageData[Binarize[Image[{{0.4, 0.6}}], 2.], \"Bit\"]", "{{0, 0}}", 0);
}

static void test_greyscale_uses_rec601_luminance(void) {
    /* EXACT weights, asserted as the values themselves: pure red, green and blue must come back
     * as 0.299, 0.587 and 0.114. An unweighted mean would give 1/3 for all three, which is the
     * mistake this row exists to catch -- and it matters for thresholding, since an average puts
     * saturated red and saturated blue on the same side of any threshold when perceptually they
     * are far apart. */
    assert_eval_eq("ImageData[ColorConvert[Image[{{{1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}}}],"
                   " \"Grayscale\"]]", "{{0.299, 0.587, 0.114}}", 0);
    /* One channel out, whatever went in. */
    assert_eval_eq("ImageChannels[ColorConvert[Image[{{{1.,0.,0.}}}], \"Grayscale\"]]", "1", 0);
    /* Already grey is a no-op on the values. */
    assert_eval_eq("ImageData[ColorConvert[Image[{{0.25, 0.75}}], \"Grayscale\"]]",
                   "{{0.25, 0.75}}", 0);
    /* Binarize on colour goes through luminance, so red (0.299) and blue (0.114) separate. */
    assert_eval_eq("ImageData[Binarize[Image[{{{1.,0.,0.}, {0.,0.,1.}}}], 0.2], \"Bit\"]",
                   "{{1, 0}}", 0);
    /* Other colour spaces are declined rather than approximated. */
    assert_eval_eq("Head[ColorConvert[Image[{{{1.,0.,0.}}}], \"LAB\"]]", "ColorConvert", 0);
}

static void test_degenerate_and_malformed_decline(void) {
    /* A constant image is ONE cluster: no threshold splits it in two, and inventing one would be
     * a fiction. Declining is the honest answer. */
    assert_eval_eq("Head[FindThreshold[Image[Table[0.5, {3}, {3}]]]]", "FindThreshold", 0);
    assert_eval_eq("Head[Binarize[Image[Table[0.5, {3}, {3}]]]]", "Binarize", 0);
    /* But an EXPLICIT threshold on a constant image is perfectly well defined. */
    assert_eval_eq("ImageData[Binarize[Image[Table[0.5, {2}, {2}]], 0.4], \"Bit\"]",
                   "{{1, 1}, {1, 1}}", 0);
    assert_eval_eq("Head[FindThreshold[{{1, 2}}]]", "FindThreshold", 0);
    assert_eval_eq("Head[Binarize[{{1, 2}}]]", "Binarize", 0);
    assert_eval_eq("Head[Binarize[" BIMODAL ", x]]", "Binarize", 0);
    assert_eval_eq("And @@ Map[MemberQ[Attributes[#], Protected] &,"
                   " {Binarize, FindThreshold, ColorConvert}]", "True", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_dimensions_are_transposed_from_data);
    TEST(test_byte_scaling_is_exact);
    TEST(test_stored_values_round_trip);
    TEST(test_type_is_inferred_from_the_values);
    TEST(test_canonical_form_is_a_fixed_point);
    TEST(test_channels);
    TEST(test_malformed_input_declines);
    TEST(test_imageq_and_attributes);
    TEST(test_convolution_reflects_the_kernel);
    TEST(test_identity_kernel_is_the_identity);
    TEST(test_constant_image_survives_a_normalised_kernel);
    TEST(test_kernel_constructors);
    TEST(test_gaussianfilter_equals_imageconvolve);
    TEST(test_filters_handle_colour_and_decline_junk);
    TEST(test_otsu_splits_a_bimodal_image_perfectly);
    TEST(test_otsu_is_the_argmax_of_between_class_variance);
    TEST(test_binarize_boundary_is_strictly_above);
    TEST(test_greyscale_uses_rec601_luminance);
    TEST(test_degenerate_and_malformed_decline);

    printf("All image tests passed.\n");
    return 0;
}
