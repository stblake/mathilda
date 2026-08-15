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

#define CHK4 "Image[{{0.,1.,0.,1.},{1.,0.,1.,0.},{0.,1.,0.,1.},{1.,0.,1.,0.}}]"

static void test_downsampling_does_not_alias(void) {
    /* THE row, and the analogue of the convolution-reflection test: it separates a correct
     * resampler from one that merely looks correct on smooth data.
     *
     * A 4x4 checkerboard reduced to 2x2 is sampled at exactly the frequency that annihilates the
     * pattern. Area averaging returns its mean, 0.5 everywhere. Nearest-neighbour returns a FLAT
     * FIELD -- the pattern is simply gone -- which is what aliasing looks like and what no amount
     * of interpolation afterwards could undo. Exact values either way, so nothing can land
     * between them by accident. */
    assert_eval_eq("ImageData[ImageResize[" CHK4 ", {2, 2}]]",
                   "{{0.5, 0.5}, {0.5, 0.5}}", 0);
    assert_eval_eq("ImageData[ImageResize[" CHK4 ", {2, 2}, Resampling -> \"Nearest\"]]",
                   "{{0.0, 0.0}, {0.0, 0.0}}", 0);
    /* Automatic must CHOOSE area averaging when shrinking; if the default ever flipped to
     * bilinear or nearest this fails, which is the point. */
    assert_eval_eq("ImageData[ImageResize[" CHK4 ", {2, 2}]] === "
                   "ImageData[ImageResize[" CHK4 ", {2, 2}, Resampling -> \"Average\"]]",
                   "True", 0);
}

static void test_area_averaging_is_an_exact_block_mean(void) {
    /* At an integer reduction factor the fractional coverage terms vanish and the result is an
     * exact mean of each source block -- so these are equalities. Block (1,1) covers
     * {0, 0.25, 0, 0.25} with mean 0.125; block (1,2) covers {0.5, 0.75, 0.5, 0.75} with mean
     * 0.625. Chosen so both are exactly representable. */
    assert_eval_eq("ImageData[ImageResize[Image[{{0.,0.25,0.5,0.75},{0.,0.25,0.5,0.75},"
                   "{1.,1.,1.,1.},{1.,1.,1.,1.}}], {2, 2}]]",
                   "{{0.125, 0.625}, {1.0, 1.0}}", 0);
    /* FRACTIONAL coverage, on a 3 -> 2 reduction where no block boundary lines up. Destination
     * pixel 0 spans source [0, 1.5): all of pixel 0 (value 0) plus half of pixel 1 (value 0.5),
     * so (0*1 + 0.5*0.5)/1.5 = 1/6. Restricting the area path to integer factors would have
     * silently fallen back to something worse on exactly the sizes people ask for. */
    assert_eval_eq("Module[{r = ImageData[ImageResize[Image[{{0., 0.5, 1.}}], {2, 1}]]},"
                   " {Chop[Part[r,1,1] - 1./6.], Chop[Part[r,1,2] - 5./6.]}]", "{0, 0}", 0);
}

static void test_resizing_to_the_same_size_is_the_identity(void) {
    /* Every method must be the identity at 1:1 -- exactly, not approximately. This is what
     * catches a half-pixel shift: the naive sx = i * scale is right at 1:1 and wrong everywhere
     * else, so a test only at 1:1 would MISS it, which is why the fractional and checkerboard
     * rows above exist too. */
    assert_eval_eq("ImageData[ImageResize[" I23 ", {3, 2}, Resampling -> \"Bilinear\"]] === "
                   "ImageData[" I23 "]", "True", 0);
    assert_eval_eq("ImageData[ImageResize[" I23 ", {3, 2}, Resampling -> \"Average\"]] === "
                   "ImageData[" I23 "]", "True", 0);
    assert_eval_eq("ImageData[ImageResize[" I23 ", {3, 2}, Resampling -> \"Nearest\"]] === "
                   "ImageData[" I23 "]", "True", 0);
    /* A constant image survives ENLARGEMENT unchanged, which bilinear must give exactly. */
    assert_eval_eq("Union[Flatten[ImageData[ImageResize[Image[{{0.25,0.25},{0.25,0.25}}],"
                   " {4, 4}]]]]", "{0.25}", 0);
}

static void test_size_specification(void) {
    /* A single number is a WIDTH and the height follows the aspect ratio: 4x2 to width 2 is
     * height 1. */
    assert_eval_eq("ImageDimensions[ImageResize[Image[Table[0.5, {2}, {4}]], 2]]", "{2, 1}", 0);
    assert_eval_eq("ImageDimensions[ImageResize[" I23 ", {7, 5}]]", "{7, 5}", 0);
    /* Enlarging and shrinking both land on exactly the requested size. */
    assert_eval_eq("ImageDimensions[ImageResize[" I23 ", {1, 1}]]", "{1, 1}", 0);
    /* Colour channels and the Real result type survive a resize. */
    assert_eval_eq("{ImageChannels[ImageResize[Image[{{{1.,0.,0.},{0.,1.,0.}}}], {1, 1}]],"
                   " ImageType[ImageResize[" I23 ", {2, 2}]]}", "{3, \"Real\"}", 0);
}

static void test_resize_declines_bad_input(void) {
    /* A fractional pixel count has no meaning, and rounding one silently would make
     * ImageResize[img, 10.5] quietly mean something the caller did not say. */
    assert_eval_eq("Head[ImageResize[" I23 ", 2.5]]", "ImageResize", 0);
    assert_eval_eq("Head[ImageResize[" I23 ", 0]]", "ImageResize", 0);
    assert_eval_eq("Head[ImageResize[" I23 ", -3]]", "ImageResize", 0);
    assert_eval_eq("Head[ImageResize[" I23 ", {2, 0}]]", "ImageResize", 0);
    assert_eval_eq("Head[ImageResize[" I23 ", {2, 3, 4}]]", "ImageResize", 0);
    assert_eval_eq("Head[ImageResize[{{1, 2}}, 2]]", "ImageResize", 0);
    /* An unknown resampling name declines rather than falling back to a default -- a silent
     * fallback would make a typo look like it worked. */
    assert_eval_eq("Head[ImageResize[" I23 ", 2, Resampling -> \"Bogus\"]]", "ImageResize", 0);
    assert_eval_eq("Head[ImageResize[" I23 ", 2, Bogus -> \"Nearest\"]]", "ImageResize", 0);
    assert_eval_eq("MemberQ[Attributes[ImageResize], Protected]", "True", 0);
}

static void test_buffer_backed_storage_agrees_with_the_nested_path(void) {
    /* Above the packing threshold a filter's output is stored as a machine BUFFER rather than
     * nested expressions. The two representations must be indistinguishable through the API, and
     * the identity kernel makes that an EXACT equality rather than a tolerance: convolving with
     * {{1}} goes in through the nested path and comes out through the buffer path, so if the two
     * disagreed anywhere this would fail.
     *
     * The Head assertion is what stops this test passing vacuously. Every other row in this file
     * uses a small image that stays nested, so none of them touch the buffer path at all -- which
     * is exactly how the first version of this change shipped an ImageData that could not read its
     * own storage. */
    assert_eval_eq("Module[{d, img, out},"
                   " d = Table[N[Mod[x*7 + y*13, 251]/251], {y, 64}, {x, 64}];"
                   " img = Image[d];"
                   " out = ImageConvolve[img, {{1}}];"
                   " {Head[Part[out, 1]], ImageData[out] === ImageData[img]}]",
                   "{NDArray, True}", 0);
    /* Storage is buffer-backed UNIFORMLY, at every size -- ndbuild_open packs even six elements,
     * so there is no size threshold at which the representation changes. That is the better
     * outcome: one representation for every computed image, and no class of image that only the
     * large-input tests would ever reach. What stays nested is the data a CALLER handed to
     * Image[...], which is copied as given rather than converted. */
    assert_eval_eq("Head[Part[ImageConvolve[" I23 ", {{1}}], 1]]", "NDArray", 0);
    assert_eval_eq("Head[Part[" I23 ", 1]]", "List", 0);
    /* And the two agree exactly, small as well as large: an identity convolution of the nested
     * form gives back the same pixels through the buffer form. */
    assert_eval_eq("ImageData[ImageConvolve[" I23 ", {{1}}]] === ImageData[" I23 "]", "True", 0);
    /* ImageData always answers with a List whatever the storage -- a caller writing
     * Part[ImageData[img], y, x] must not have to care which side of the threshold it landed. */
    assert_eval_eq("Module[{out = ImageConvolve[Image[Table[0.5, {64}, {64}]], {{1}}]},"
                   " {Head[Part[out, 1]], Head[ImageData[out]],"
                   "  Dimensions[ImageData[out]], Union[Flatten[ImageData[out]]]}]",
                   "{NDArray, List, {64, 64}, {0.5}}", 0);
}

static void test_every_accessor_reads_buffer_storage(void) {
    /* Each accessor and each filter must work on buffer-backed storage, not just on nested. These
     * all go through image_info's O(1) dims read or image_load's buffer read. */
    assert_eval_eq("Module[{g = ImageConvolve[Image[Table[N[x/70.], {y, 64}, {x, 64}]],"
                   " GaussianMatrix[1]]},"
                   " {ImageQ[g], ImageDimensions[g], ImageChannels[g], ImageType[g],"
                   "  Dimensions[ImageData[g]]}]",
                   "{True, {64, 64}, 1, \"Real\", {64, 64}}", 0);
    /* Colour keeps its rank-3 buffer and its per-channel values. */
    assert_eval_eq("Module[{c = ImageConvolve[Image[Table[{0.5, 0.25, 0.75}, {40}, {40}]],"
                   " GaussianMatrix[1]]},"
                   " {Head[Part[c, 1]], ImageChannels[c], Dimensions[ImageData[c]],"
                   "  Part[ImageData[c], 20, 20]}]",
                   "{NDArray, 3, {40, 40, 3}, {0.5, 0.25, 0.75}}", 0);
    /* Thresholding and resizing consume buffer storage too. */
    assert_eval_eq("Module[{g = ImageConvolve[Image[Table[N[x/70.], {y, 64}, {x, 64}]],"
                   " GaussianMatrix[1]]},"
                   " {ImageType[Binarize[g]], ImageDimensions[ImageResize[g, {16, 16}]],"
                   "  NumberQ[FindThreshold[g]]}]",
                   "{\"Bit\", {16, 16}, True}", 0);
}

static void test_separable_kernels_are_detected_and_decomposed(void) {
    /* An outer product is rank 1, so ImageConvolve factors it and runs two 1-D passes. The
     * INDEPENDENT check that the factorisation is right: {{1,2},{2,4}} is u (x) v for u = {1,2} and
     * v = {1,2}, so convolving with it must equal convolving with {{1,2}} and then {{1},{2}} --
     * two genuine one-dimensional convolutions, arrived at by a different route through the code.
     * If the decomposition picked the wrong pivot or scaled a factor wrongly, these would differ. */
    assert_eval_eq("Module[{img = Image[Table[N[Mod[x*3 + y*5, 17]/17.], {y, 9}, {x, 7}]]},"
                   " Chop[Max[Abs[Flatten["
                   "   ImageData[ImageConvolve[img, {{1, 2}, {2, 4}}]]"
                   " - ImageData[ImageConvolve[ImageConvolve[img, {{1, 2}}], {{1}, {2}}]]]]]]]",
                   "0", 0);
    /* A kernel with a ZERO in the corner is still separable, and pivoting on the LARGEST entry
     * rather than on K[[1,1]] is what makes that work -- pivoting on a zero would divide by it.
     * {{0,0,0},{1,2,1},{0,0,0}} is {0,1,0} (x) {1,2,1}, so it must equal the same two 1-D
     * convolutions. This row is what fails if the pivot search is ever simplified away. */
    assert_eval_eq("Module[{img = Image[Table[N[Mod[x*3 + y*5, 17]/17.], {y, 6}, {x, 6}]]},"
                   " Chop[Max[Abs[Flatten["
                   "   ImageData[ImageConvolve[img, {{0,0,0}, {1,2,1}, {0,0,0}}]]"
                   " - ImageData[ImageConvolve[ImageConvolve[img, {{1,2,1}}],"
                   "                           {{0}, {1}, {0}}]]]]]]]", "0", 0);
}

static void test_non_separable_kernels_take_the_direct_path(void) {
    /* THE correctness row for separability. Treating a rank-2 kernel as rank 1 does not make the
     * answer slightly wrong -- it computes a completely different filter -- so the detector must
     * reject one, and the tolerance is relative and tight (1e-12) for exactly that reason.
     *
     * {{1,0},{0,1}} is rank 2. Its centre is (1,1) by floor division, so
     * out[y][x] = s[y+1][x+1] + s[y][x]; a delta at the centre of a 3x3 therefore gives
     * {{1,0,0},{0,1,0},{0,0,0}}. Computed by hand, and no rank-1 approximation of that kernel
     * produces it. */
    assert_eval_eq("ImageData[ImageConvolve[Image[{{0.,0.,0.},{0.,1.,0.},{0.,0.,0.}}],"
                   " {{1, 0}, {0, 1}}]]",
                   "{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}}", 0);
    /* An all-zero kernel has no factorisation at all and must not divide by its pivot. */
    assert_eval_eq("Union[Flatten[ImageData[ImageConvolve[" I23 ", {{0, 0}, {0, 0}}]]]]",
                   "{0.0}", 0);
}

static void test_separability_preserves_the_documented_identity(void) {
    /* GaussianFilter and ImageConvolve[.., GaussianMatrix[..]] must remain BIT-identical, not
     * merely close. Both route through the same separable path when the kernel is rank 1, which is
     * precisely why separability was put in ImageConvolve rather than only in GaussianFilter --
     * doing it in one of them would have made the documented identity approximate. */
    assert_eval_eq("ImageData[GaussianFilter[" I23 ", 2]] === "
                   "ImageData[ImageConvolve[" I23 ", GaussianMatrix[2]]]", "True", 0);
    /* And the earlier exact answers survive the new path: a 1xN kernel is trivially rank 1. */
    assert_eval_eq("ImageData[ImageConvolve[Image[{{0., 1., 0.}}], {{1, 2, 3}}]]",
                   "{{1.0, 2.0, 3.0}}", 0);
    /* A constant image through a normalised separable kernel is still unchanged at the border,
     * which is what fails if the two passes clamp differently from the direct form. */
    assert_eval_eq("Chop[Max[Abs[Flatten[ImageData[ImageConvolve["
                   "Image[Table[0.25, {6}, {6}]], GaussianMatrix[2]]] - 0.25]]]]", "0", 0);
}

#define RAMP8 "Image[Table[N[x/8.], {y, 6}, {x, 6}]]"

static void test_derivative_of_a_ramp_is_exactly_its_slope(void) {
    /* THE row, and it is the one that caught a real sign bug. A ramp rising 1/8 per pixel has
     * derivative exactly +1/8, and 0.125 is exactly representable, so this is an EQUALITY.
     *
     * The sign is the whole point. ImageConvolve reflects its kernel, so a central-difference
     * stencil written in natural reading order computes the NEGATED derivative. The gradient
     * MAGNITUDE squares that away and looked perfectly correct; so would any edge-detection
     * result, since only the ranking of edges matters there. Only an exact signed value could see
     * it -- an absolute-property test on the magnitude would have shipped the bug. The convention
     * now matches scipy's correlate to the digit. */
    assert_eval_eq("Part[ImageData[DerivativeFilter[" RAMP8 ", {0, 1}]], 3, 3] == 0.125",
                   "True", 0);
    /* An x-ramp has NO y-derivative, exactly zero everywhere including the border. */
    assert_eval_eq("Union[Flatten[ImageData[DerivativeFilter[" RAMP8 ", {1, 0}]]]]", "{0.0}", 0);
    /* A y-ramp is the mirror case, so a transposed stencil index cannot hide. */
    assert_eval_eq("Part[ImageData[DerivativeFilter["
                   "Image[Table[N[y/8.], {y, 6}, {x, 6}]], {1, 0}]], 3, 3] == 0.125", "True", 0);
    /* Every derivative of a constant is exactly zero -- including the second, and including the
     * gradient magnitude. */
    assert_eval_eq("Module[{c = Image[Table[0.4, {5}, {5}]]},"
                   " {Union[Flatten[ImageData[DerivativeFilter[c, {0, 1}]]]],"
                   "  Union[Flatten[ImageData[DerivativeFilter[c, {1, 0}]]]],"
                   "  Union[Flatten[ImageData[DerivativeFilter[c, {0, 2}]]]],"
                   "  Union[Flatten[ImageData[GradientFilter[c]]]]}]",
                   "{{0.0}, {0.0}, {0.0}, {0.0}}", 0);
}

static void test_second_derivative_is_exact_on_a_quadratic(void) {
    /* The second difference {1, -2, 1} on f(x) = x^2 gives exactly 2, for the same reason the
     * first is exact on a ramp: the stencil is the exact finite difference, not an approximation
     * of one. A wrong normalisation here would show as 1 or 4 rather than 2. */
    assert_eval_eq("Part[ImageData[DerivativeFilter["
                   "Image[Table[N[x*x], {y, 6}, {x, 6}]], {0, 2}]], 3, 3] == 2.", "True", 0);
    /* Order 0 is a smoothing, and being symmetric and normalised it preserves a LINEAR ramp
     * exactly -- which is also what makes the mixed orders exact. */
    assert_eval_eq("Part[ImageData[DerivativeFilter[" RAMP8 ", {0, 0}]], 3, 3] == 0.375",
                   "True", 0);
}

static void test_gradient_magnitude_is_rotation_invariant(void) {
    /* Rotation invariance is why the magnitude is Sqrt[dx^2 + dy^2] and not |dx| + |dy|. A
     * DIAGONAL ramp x + y has gradient of magnitude slope*Sqrt[2]; the absolute sum would report
     * 2*slope, biasing every downstream threshold by orientation. With slope 1/8 the answer is
     * 0.125 Sqrt[2], asserted to rounding. */
    assert_eval_eq("Chop[Part[ImageData[GradientFilter["
                   "Image[Table[N[(x + y)/8.], {y, 6}, {x, 6}]]]], 3, 3] - 0.125 Sqrt[2.]]",
                   "0", 0);
    /* On an axis-aligned ramp the magnitude is just the slope. */
    assert_eval_eq("Part[ImageData[GradientFilter[" RAMP8 "]], 3, 3] == 0.125", "True", 0);
    /* Magnitude is non-negative everywhere, whatever the sign of the derivatives. */
    assert_eval_eq("Module[{g = Flatten[ImageData[GradientFilter["
                   "Image[Table[N[Mod[x*3 + y*5, 7]/7.], {y, 8}, {x, 8}]]]]]},"
                   " And @@ Map[# >= 0. &, g]]", "True", 0);
}

static void test_derivative_filter_matches_an_explicit_kernel(void) {
    /* DerivativeFilter must be exactly ImageConvolve with the kernel it claims to use -- the
     * pre-flipped outer product of {1,2,1}/4 and {1/2, 0, -1/2}. Writing the kernel out here is
     * what pins the SIGN CONVENTION in a place a reader can check, rather than leaving it implicit
     * in a stencil table. */
    assert_eval_eq("ImageData[DerivativeFilter[" RAMP8 ", {0, 1}]] === "
                   "ImageData[ImageConvolve[" RAMP8 ", "
                   "{{0.125, 0., -0.125}, {0.25, 0., -0.25}, {0.125, 0., -0.125}}]]", "True", 0);
    /* Orders outside 0..2 and a malformed spec decline rather than guessing a stencil. */
    assert_eval_eq("Head[DerivativeFilter[" RAMP8 ", {0, 3}]]", "DerivativeFilter", 0);
    assert_eval_eq("Head[DerivativeFilter[" RAMP8 ", {-1, 0}]]", "DerivativeFilter", 0);
    assert_eval_eq("Head[DerivativeFilter[" RAMP8 ", 1]]", "DerivativeFilter", 0);
    assert_eval_eq("Head[DerivativeFilter[" RAMP8 ", {0, 1.5}]]", "DerivativeFilter", 0);
    assert_eval_eq("Head[DerivativeFilter[{{1, 2}}, {0, 1}]]", "DerivativeFilter", 0);
    assert_eval_eq("Head[GradientFilter[{{1, 2}}]]", "GradientFilter", 0);
    /* Colour goes through luminance, so the result is single-channel. */
    assert_eval_eq("ImageChannels[GradientFilter[Image[{{{1.,0.,0.},{0.,1.,0.}}}]]]", "1", 0);
    assert_eval_eq("And @@ Map[MemberQ[Attributes[#], Protected] &,"
                   " {DerivativeFilter, GradientFilter}]", "True", 0);
}

#define STEPV "Image[Table[If[x >= 5, 1., 0.], {y, 8}, {x, 8}]]"
#define STEPH "Image[Table[If[y >= 5, 1., 0.], {y, 8}, {x, 8}]]"

static void test_canny_thins_an_edge_to_one_pixel(void) {
    /* THE property that separates an edge detector from a thick mask, and it is exact.
     *
     * A clean step does NOT give a single-pixel gradient peak: the central difference responds 0.5
     * at both k-1 and k, an even plateau. Testing `mag >= both neighbours` keeps both and yields a
     * two-pixel edge -- which is what scikit-image's canny does on this input, giving width 2. The
     * asymmetric test `mag > backward && mag >= forward` resolves the tie to the lower-index side
     * and thins to one. Radius 0 keeps the geometry exact so this is a pixel pattern, not a
     * tendency. */
    assert_eval_eq("Map[Total, ImageData[EdgeDetect[" STEPV ", 0], \"Bit\"]]",
                   "{1, 1, 1, 1, 1, 1, 1, 1}", 0);
    /* Same for a horizontal step, transposed -- so a swapped index cannot hide in one axis. */
    assert_eval_eq("Map[Total, Transpose[ImageData[EdgeDetect[" STEPH ", 0], \"Bit\"]]]",
                   "{1, 1, 1, 1, 1, 1, 1, 1}", 0);
    /* And it lands in the same column on every row: a straight edge stays straight. */
    assert_eval_eq("Union[Map[First[Flatten[Position[#, 1]]] &, "
                   "ImageData[EdgeDetect[" STEPV ", 0], \"Bit\"]]]", "{4}", 0);
    /* Border rows are KEPT, not discarded: off-the-edge counts as zero magnitude, so a border
     * pixel is kept if it beats the one neighbour it has. Clamping instead would compare a pixel
     * against itself and keep every border pixel unconditionally; discarding the border, as
     * scikit-image does, loses real edges that run along it. */
    assert_eval_eq("{First[Map[Total, ImageData[EdgeDetect[" STEPV ", 0], \"Bit\"]]],"
                   " Last[Map[Total, ImageData[EdgeDetect[" STEPV ", 0], \"Bit\"]]]}",
                   "{1, 1}", 0);
}

static void test_canny_finds_no_edges_where_there_are_none(void) {
    /* A constant image has no gradient anywhere, so Otsu finds no two classes and the honest answer
     * is no edges -- not an arbitrary threshold applied to numerical noise. */
    assert_eval_eq("Union[Flatten[ImageData[EdgeDetect[Image[Table[0.5, {6}, {6}]], 0],"
                   " \"Bit\"]]]", "{0}", 0);
    /* A threshold above any possible gradient suppresses everything, which is what shows the
     * explicit-threshold argument is really used rather than ignored. */
    assert_eval_eq("Union[Flatten[ImageData[EdgeDetect[" STEPV ", 0, 10.], \"Bit\"]]]",
                   "{0}", 0);
    /* A threshold below the step's gradient keeps the edge, still one pixel wide. */
    assert_eval_eq("Map[Total, ImageData[EdgeDetect[" STEPV ", 0, 0.2], \"Bit\"]]",
                   "{1, 1, 1, 1, 1, 1, 1, 1}", 0);
}

static void test_hysteresis_drops_isolated_weak_edges(void) {
    /* Hysteresis is the reason for two thresholds. A step of height 0.3 has gradient 0.15, which
     * sits between 0.4*t and t for t = 0.25 -- so it is a WEAK edge. Isolated, with no strong pixel
     * anywhere to be 8-connected to, it must be dropped entirely. If the implementation collapsed
     * to a single threshold at the low value, this would come back as a full edge. */
    assert_eval_eq("Union[Flatten[ImageData[EdgeDetect["
                   "Image[Table[If[x >= 5, 0.3, 0.], {y, 8}, {x, 8}]], 0, 0.25], \"Bit\"]]]",
                   "{0}", 0);
    /* The same weak step at a threshold it exceeds outright is kept -- so the row above is about
     * hysteresis and not about the edge being invisible. */
    assert_eval_eq("Map[Total, ImageData[EdgeDetect["
                   "Image[Table[If[x >= 5, 0.3, 0.], {y, 8}, {x, 8}]], 0, 0.1], \"Bit\"]]",
                   "{1, 1, 1, 1, 1, 1, 1, 1}", 0);
}

static void test_edgedetect_shape_and_declines(void) {
    assert_eval_eq("{ImageType[EdgeDetect[" STEPV ", 0]], ImageDimensions[EdgeDetect[" STEPV ", 0]]}",
                   "{\"Bit\", {8, 8}}", 0);
    /* Colour goes through luminance, so the result is single channel. */
    assert_eval_eq("ImageChannels[EdgeDetect[Image[Table[{0.2, 0.6, 0.9}, {8}, {8}]], 0]]",
                   "1", 0);
    assert_eval_eq("Head[EdgeDetect[{{1, 2}}]]", "EdgeDetect", 0);
    assert_eval_eq("Head[EdgeDetect[" STEPV ", 1.5]]", "EdgeDetect", 0);
    assert_eval_eq("Head[EdgeDetect[" STEPV ", -1]]", "EdgeDetect", 0);
    assert_eval_eq("Head[EdgeDetect[" STEPV ", 0, -0.5]]", "EdgeDetect", 0);
    assert_eval_eq("MemberQ[Attributes[EdgeDetect], Protected]", "True", 0);
}

#define VOL234 "Image3D[Table[N[(x + 10 y + 100 z)/1000.], {z, 2}, {y, 3}, {x, 4}]]"

static void test_image3d_dimensions_are_fully_reversed(void) {
    /* THE row for volumes, and the 3-D version of the 2-D transposition trap -- only worse: with
     * three axes there are SIX possible orderings, and a cubic test volume validates none of them.
     * So this uses 2 slices of 3 rows of 4 columns, every extent distinct.
     *
     * Storage is depth x height x width, indexed data[[z, y, x]]; ImageDimensions reports
     * {width, height, depth}, fully reversed. Both are asserted, because either one alone would
     * pass with two axes swapped. */
    assert_eval_eq("ImageDimensions[" VOL234 "]", "{4, 3, 2}", 0);
    assert_eval_eq("Dimensions[ImageData[" VOL234 "]]", "{2, 3, 4}", 0);
    /* And the voxel the indices imply really is the one stored: (x=4) + 10(y=3) + 100(z=2) over
     * 1000 is 0.234, exactly representable. */
    assert_eval_eq("Part[ImageData[" VOL234 "], 2, 3, 4] == 0.234", "True", 0);
    assert_eval_eq("Part[ImageData[" VOL234 "], 1, 1, 1] == 0.111", "True", 0);
}

static void test_image3d_is_distinct_from_image(void) {
    /* A volume is not a 2-D image and must not answer to ImageQ, or every filter written for a
     * plane would silently accept one and index it wrongly. */
    assert_eval_eq("{Image3DQ[" VOL234 "], ImageQ[" VOL234 "]}", "{True, False}", 0);
    /* And conversely a plane is not a volume. */
    assert_eval_eq("{Image3DQ[" I23 "], ImageQ[" I23 "]}", "{False, True}", 0);
    /* Malformed input stays unevaluated, so Image3DQ is how validity is tested. */
    assert_eval_eq("Image3DQ[Image3D[{{{1, 2}}, {{3}}}]]", "False", 0);
    assert_eval_eq("Image3DQ[{{{1}}}]", "False", 0);
    assert_eval_eq("Head[Image3D[{{1, 2}, {3, 4}}]]", "Image3D", 0);
    assert_eval_eq("Head[Image3D[{}]]", "Image3D", 0);
    assert_eval_eq("Head[Image3D[{{{1, x}}}]]", "Image3D", 0);
}

static void test_image3d_shares_the_type_rules(void) {
    /* Type inference and ImageData's scaling are the same code for both ranks, and this pins that:
     * a byte volume scales 255 to exactly 1.0 and 128 to 128/255, as a byte plane does. */
    assert_eval_eq("ImageType[Image3D[{{{0, 255}}, {{128, 64}}}]]", "\"Byte\"", 0);
    assert_eval_eq("Part[ImageData[Image3D[{{{0, 255}}, {{128, 64}}}]], 1, 1, 2] == 1.0",
                   "True", 0);
    assert_eval_eq("ImageType[Image3D[{{{0, 1}}, {{1, 0}}}]]", "\"Bit\"", 0);
    assert_eval_eq("ImageType[Image3D[{{{0., 0.5}}}]]", "\"Real\"", 0);
    /* A stated type inconsistent with the data declines, as in 2-D. */
    assert_eval_eq("Head[Image3D[{{{0, 300}}}, \"Byte\"]]", "Image3D", 0);
    /* Canonical form is a fixed point. */
    assert_eval_eq("Length[Image3D[{{{0, 1}}}]]", "2", 0);
}

static void test_image3d_colour_volumes(void) {
    /* A colour volume is RANK 4 -- depth x height x width x channels -- which is the case that
     * forced ImageData's nested rebuild to become a general recursion over dims rather than the
     * unrolled two-and-a-bit levels the plane case uses. */
    assert_eval_eq("Module[{cv = Image3D[Table[{0.5, 0.25, 0.75}, {z, 2}, {y, 2}, {x, 3}]]},"
                   " {ImageChannels[cv], ImageDimensions[cv], Dimensions[ImageData[cv]],"
                   "  Part[ImageData[cv], 1, 1, 1]}]",
                   "{3, {3, 2, 2}, {2, 2, 3, 3}, {0.5, 0.25, 0.75}}", 0);
    /* Grey volumes report one channel. */
    assert_eval_eq("ImageChannels[" VOL234 "]", "1", 0);
    assert_eval_eq("And @@ Map[MemberQ[Attributes[#], Protected] &, {Image3D, Image3DQ}]",
                   "True", 0);
}

static void test_volume_convolution_reflects_on_every_axis(void) {
    /* The z-axis version of the 2-D reflection test. A delta in a 3-slice column with the kernel
     * varying only in z gives {1, 2, 3} under convolution and {3, 2, 1} under correlation -- exact
     * integers, and it covers the axis the 2-D tests cannot reach. */
    assert_eval_eq("Flatten[ImageData[ImageConvolve[Image3D[{{{0.}}, {{1.}}, {{0.}}}],"
                   " {{{1}}, {{2}}, {{3}}}]]]", "{1.0, 2.0, 3.0}", 0);
    /* Identity kernel returns the voxels exactly, and the dimensions survive. */
    assert_eval_eq("ImageData[ImageConvolve[" VOL234 ", {{{1}}}]] === ImageData[" VOL234 "]",
                   "True", 0);
    assert_eval_eq("ImageDimensions[ImageConvolve[" VOL234 ", {{{1}}}]]", "{4, 3, 2}", 0);
    /* A constant volume through a normalised kernel is unchanged everywhere, borders included --
     * which is only true because the per-axis clamping composes. */
    assert_eval_eq("Chop[Max[Abs[Flatten[ImageData[GaussianFilter["
                   "Image3D[Table[0.25, {4}, {4}, {4}]], 1]]] - 0.25]]]", "0", 0);
}

static void test_volume_separability_equals_three_1d_passes(void) {
    /* The independent check on the rank-3 factorisation, and the analogue of the 2-D outer-product
     * row: a separable 3x3x3 kernel must give exactly what three successive 1-D convolutions give.
     * Reached by a different route through the code, so a mis-scaled factor or a wrong pivot shows
     * up here rather than being invisible. */
    assert_eval_eq("Module[{w = Image3D[Table[N[Mod[x*3 + y*5 + z*7, 11]/11.],"
                   " {z, 4}, {y, 5}, {x, 6}]]},"
                   " Chop[Max[Abs[Flatten["
                   "   ImageData[ImageConvolve[w, Table[a*b*c, {a, {1,2,1}}, {b, {1,2,1}},"
                   "                                          {c, {1,2,1}}]]]"
                   " - ImageData[ImageConvolve[ImageConvolve[ImageConvolve[w,"
                   "     {{{1,2,1}}}], {{{1},{2},{1}}}], {{{1}},{{2}},{{1}}}]]]]]]]",
                   "0", 0);
    /* A non-separable rank-3 kernel must still run, on the direct path, and return a volume. */
    assert_eval_eq("Module[{r = ImageConvolve[Image3D[{{{0.}}, {{1.}}, {{0.}}}],"
                   " {{{1, 0}, {0, 0}}, {{0, 0}, {0, 1}}}]},"
                   " {Image3DQ[r], ImageDimensions[r]}]", "{True, {1, 1, 3}}", 0);
}

static void test_volume_gaussian_and_small_volumes(void) {
    /* THE row for the nested fallback, which was a live bug: ndbuild_open declines an array under
     * the packing threshold, so a small volume returned NULL and the whole convolution DECLINED
     * instead of answering. The 2-D builder always had a nested fallback; the 3-D one was written
     * without it, and only large volumes had been tried, so it looked correct. Three voxels and
     * eight voxels are both below the threshold. */
    assert_eval_eq("Image3DQ[ImageConvolve[Image3D[{{{0.}}, {{1.}}, {{0.}}}], {{{1}}}]]",
                   "True", 0);
    assert_eval_eq("Image3DQ[GaussianFilter[Image3D[Table[0.5, {2}, {2}, {2}]], 1]]", "True", 0);
    /* Shape, type and channels through a 3-D Gaussian, including a colour volume. */
    assert_eval_eq("Module[{w = Image3D[Table[N[x/8.], {z, 4}, {y, 5}, {x, 6}]]},"
                   " {ImageDimensions[GaussianFilter[w, 1]], ImageType[GaussianFilter[w, 1]],"
                   "  ImageChannels[GaussianFilter[w, 1]]}]", "{{6, 5, 4}, \"Real\", 1}", 0);
    assert_eval_eq("ImageChannels[GaussianFilter["
                   "Image3D[Table[{0.5, 0.25, 0.75}, {2}, {3}, {4}]], 1]]", "3", 0);
    /* A rank-2 kernel on a volume declines rather than being reinterpreted as a stack of planes. */
    assert_eval_eq("Head[ImageConvolve[" VOL234 ", {{1, 2}, {3, 4}}]]", "ImageConvolve", 0);
    assert_eval_eq("Head[GaussianFilter[" VOL234 ", 1.5]]", "GaussianFilter", 0);
    assert_eval_eq("Head[GaussianFilter[" VOL234 ", -1]]", "GaussianFilter", 0);
}

#define MIMG "Image[Table[N[Mod[x*3 + y*5, 7]/7.], {y, 7}, {x, 8}]]"

static void test_morphology_ordering_chain(void) {
    /* Erosion <= Opening <= f <= Closing <= Dilation, POINTWISE and everywhere including the
     * border. One assertion covering four operators and their relationship to the original, which
     * no per-operator check would catch: a wrong element centre shifts one of them off the chain
     * while leaving each individually plausible. It holds at the border only because the padding
     * replicates rather than zero-fills. */
    assert_eval_eq("Module[{f = ImageData[" MIMG "], d, e, o, c},"
                   " d = ImageData[Dilation[" MIMG ", 1]]; e = ImageData[Erosion[" MIMG ", 1]];"
                   " o = ImageData[Opening[" MIMG ", 1]]; c = ImageData[Closing[" MIMG ", 1]];"
                   " And @@ MapThread[#1 <= #2 <= #3 <= #4 <= #5 &,"
                   "   Map[Flatten, {e, o, f, c, d}]]]", "True", 0);
}

static void test_morphology_duality_is_exact(void) {
    /* Erosion[f, k] == 1 - Dilation[1 - f, k]. An EXACT identity, and it fails for a swapped
     * min/max, a mis-centred element, or a padding rule that is not self-dual -- zero padding would
     * break it at the border specifically, which is why replicate is used. */
    assert_eval_eq("Module[{f = ImageData[" MIMG "], neg},"
                   " neg = Image[Map[(1. - #) &, f, {2}]];"
                   " Chop[Max[Abs[Flatten[ImageData[Erosion[" MIMG ", 1]]"
                   "   - Map[(1. - #) &, ImageData[Dilation[neg, 1]], {2}]]]]]]", "0", 0);
}

static void test_opening_and_closing_are_idempotent(void) {
    /* THE defining property of an opening, and the reason opening twice is not a sharpening loop.
     * Exact, not approximate. A mis-composed pair -- dilate-then-erode where erode-then-dilate was
     * meant, or a different element on the second pass -- still smooths and still looks reasonable,
     * and fails here. */
    assert_eval_eq("Module[{o = ImageData[Opening[" MIMG ", 1]]},"
                   " Chop[Max[Abs[Flatten[o - ImageData[Opening[Image[o], 1]]]]]]]", "0", 0);
    assert_eval_eq("Module[{c = ImageData[Closing[" MIMG ", 1]]},"
                   " Chop[Max[Abs[Flatten[c - ImageData[Closing[Image[c], 1]]]]]]]", "0", 0);
}

static void test_dilating_a_point_gives_the_element(void) {
    /* The cleanest statement of what dilation IS: a single bright voxel spreads to exactly the
     * structuring element's footprint, so the output is the element itself. Exact pixel pattern. */
    assert_eval_eq("ImageData[Dilation[Image[Table[If[x == 3 && y == 3, 1., 0.],"
                   " {y, 5}, {x, 5}]], 1]]",
                   "{{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 1.0, 1.0, 0.0},"
                   " {0.0, 1.0, 1.0, 1.0, 0.0}, {0.0, 1.0, 1.0, 1.0, 0.0},"
                   " {0.0, 0.0, 0.0, 0.0, 0.0}}", 0);
    /* A radius of 0 is a single-pixel element, hence the identity for both operators. */
    assert_eval_eq("ImageData[Dilation[" MIMG ", 0]] === ImageData[" MIMG "]", "True", 0);
    assert_eval_eq("ImageData[Erosion[" MIMG ", 0]] === ImageData[" MIMG "]", "True", 0);
}

static void test_element_forms_and_declines(void) {
    /* An integer radius and the equivalent BoxMatrix must be the SAME operation -- which is what
     * flat morphology means: only the support enters, never the values. Both take the separable
     * path, so this also pins that the two paths agree. */
    assert_eval_eq("ImageData[Dilation[" MIMG ", BoxMatrix[1]]] === "
                   "ImageData[Dilation[" MIMG ", 1]]", "True", 0);
    /* Values are ignored, only the support: an element of 5s is the same as an element of 1s. */
    assert_eval_eq("ImageData[Dilation[" MIMG ", {{5, 5, 5}, {5, 5, 5}, {5, 5, 5}}]] === "
                   "ImageData[Dilation[" MIMG ", 1]]", "True", 0);
    /* A non-rectangular element takes the direct path and must still bracket the image. */
    assert_eval_eq("Module[{f = ImageData[" MIMG "],"
                   " cr = ImageData[Dilation[" MIMG ", {{0,1,0},{1,1,1},{0,1,0}}]]},"
                   " And @@ MapThread[#1 <= #2 &, {Flatten[f], Flatten[cr]}]]", "True", 0);
    /* An all-zero element marks no neighbourhood at all, so it declines rather than returning the
     * infinities an empty max would produce. */
    assert_eval_eq("Head[Dilation[" MIMG ", {{0, 0}, {0, 0}}]]", "Dilation", 0);
    assert_eval_eq("Head[Dilation[" MIMG ", -1]]", "Dilation", 0);
    assert_eval_eq("Head[Dilation[" MIMG ", 1.5]]", "Dilation", 0);
    assert_eval_eq("Head[Erosion[{{1, 2}}, 1]]", "Erosion", 0);
    assert_eval_eq("And @@ Map[MemberQ[Attributes[#], Protected] &,"
                   " {Dilation, Erosion, Opening, Closing}]", "True", 0);
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
    TEST(test_downsampling_does_not_alias);
    TEST(test_area_averaging_is_an_exact_block_mean);
    TEST(test_resizing_to_the_same_size_is_the_identity);
    TEST(test_size_specification);
    TEST(test_resize_declines_bad_input);
    TEST(test_buffer_backed_storage_agrees_with_the_nested_path);
    TEST(test_every_accessor_reads_buffer_storage);
    TEST(test_separable_kernels_are_detected_and_decomposed);
    TEST(test_non_separable_kernels_take_the_direct_path);
    TEST(test_separability_preserves_the_documented_identity);
    TEST(test_derivative_of_a_ramp_is_exactly_its_slope);
    TEST(test_second_derivative_is_exact_on_a_quadratic);
    TEST(test_gradient_magnitude_is_rotation_invariant);
    TEST(test_derivative_filter_matches_an_explicit_kernel);
    TEST(test_canny_thins_an_edge_to_one_pixel);
    TEST(test_canny_finds_no_edges_where_there_are_none);
    TEST(test_hysteresis_drops_isolated_weak_edges);
    TEST(test_edgedetect_shape_and_declines);
    TEST(test_image3d_dimensions_are_fully_reversed);
    TEST(test_image3d_is_distinct_from_image);
    TEST(test_image3d_shares_the_type_rules);
    TEST(test_image3d_colour_volumes);
    TEST(test_volume_convolution_reflects_on_every_axis);
    TEST(test_volume_separability_equals_three_1d_passes);
    TEST(test_volume_gaussian_and_small_volumes);
    TEST(test_morphology_ordering_chain);
    TEST(test_morphology_duality_is_exact);
    TEST(test_opening_and_closing_are_idempotent);
    TEST(test_dilating_a_point_gives_the_element);
    TEST(test_element_forms_and_declines);

    printf("All image tests passed.\n");
    return 0;
}
