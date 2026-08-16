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
    /* Storage is buffer-backed UNIFORMLY, and -- since the constructor was made to canonicalise --
     * REGARDLESS OF PROVENANCE. It did not always: data a caller handed to Image[...] used to be
     * copied as given, so a constructed image and a computed one with identical pixels were never
     * SameQ, and `===` was quietly useless on images. Both agreeing on representation is what makes
     * every exactness claim in this file assertable at all. */
    assert_eval_eq("Head[Part[ImageConvolve[" I23 ", {{1}}], 1]]", "NDArray", 0);
    assert_eval_eq("Head[Part[" I23 ", 1]]", "NDArray", 0);
    assert_eval_eq("ImageConvolve[" I23 ", {{1}}] === " I23, "True", 0);
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

static void test_connectivity_is_the_discriminating_property(void) {
    /* THE row. Two pixels touching only at a corner are ONE component under 8-connectivity and TWO
     * under 4. Every other property here -- background staying 0, labels contiguous, a single blob
     * labelling 1 -- holds under either rule, so a wrong default would pass all of them and only
     * this can tell them apart. Exact label patterns, not counts. */
    assert_eval_eq("MorphologicalComponents[Image[{{1., 0.}, {0., 1.}}]]",
                   "{{1, 0}, {0, 1}}", 0);
    assert_eval_eq("MorphologicalComponents[Image[{{1., 0.}, {0., 1.}}],"
                   " CornerNeighbors -> False]", "{{1, 0}, {0, 2}}", 0);
    /* And the counts follow, since labels are contiguous from 1. */
    assert_eval_eq("{Max[Flatten[MorphologicalComponents[Image[{{1., 0.}, {0., 1.}}]]]],"
                   " Max[Flatten[MorphologicalComponents[Image[{{1., 0.}, {0., 1.}}],"
                   "   CornerNeighbors -> False]]]}", "{1, 2}", 0);
}

static void test_a_u_shape_needs_the_union_find(void) {
    /* The case a single pass cannot handle, and the reason for union-find at all. Scanning in raster
     * order, the two arms of a U get DIFFERENT provisional labels -- nothing has connected them yet
     * -- and only the base at the bottom reveals they are one component. A one-pass implementation
     * returns two components here and looks perfectly reasonable on every convex shape.
     *
     * Asserted under 4-connectivity, so the arms cannot be joined by a diagonal instead. */
    assert_eval_eq("MorphologicalComponents[Image[{{1.,0.,1.},{1.,0.,1.},{1.,1.,1.}}],"
                   " CornerNeighbors -> False]",
                   "{{1, 0, 1}, {1, 0, 1}, {1, 1, 1}}", 0);
}

static void test_labels_are_contiguous_in_raster_order(void) {
    /* Relabelling to 1..k in raster order of first appearance is what makes Max the component count
     * and makes a label pattern assertable at all. Without it the labels would be whatever the first
     * pass allocated, with GAPS where two provisional labels were later merged. */
    assert_eval_eq("MorphologicalComponents[Image[{{1.,0.,1.},{0.,0.,0.},{1.,0.,0.}}],"
                   " CornerNeighbors -> False]", "{{1, 0, 2}, {0, 0, 0}, {3, 0, 0}}", 0);
    assert_eval_eq("Module[{lb = Flatten[MorphologicalComponents["
                   "Image[{{1.,0.,1.},{0.,0.,0.},{1.,0.,0.}}], CornerNeighbors -> False]]},"
                   " Union[Select[lb, # > 0 &]] === Range[Max[lb]]]", "True", 0);
    /* An all-background image has no components and every label is 0. */
    assert_eval_eq("MorphologicalComponents[Image[Table[0., {3}, {3}]]]",
                   "{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}", 0);
    /* A single blob is all label 1. */
    assert_eval_eq("Union[Flatten[MorphologicalComponents[Image[Table[1., {3}, {3}]]]]]",
                   "{1}", 0);
}

static void test_components_interact_with_morphology(void) {
    /* Dilation can only MERGE components, never split them, so the count cannot increase -- an
     * absolute relationship between the two features rather than a property of either alone. Four
     * isolated corners become one blob at radius 1. */
    assert_eval_eq("Module[{sep = Image[{{1.,0.,0.,1.},{0.,0.,0.,0.},{1.,0.,0.,1.}}]},"
                   " {Max[Flatten[MorphologicalComponents[sep]]],"
                   "  Max[Flatten[MorphologicalComponents[Dilation[sep, 1]]]]}]", "{4, 1}", 0);
    /* The threshold argument selects the foreground: above t, not at or above. */
    assert_eval_eq("MorphologicalComponents[Image[{{0.3, 0.7}, {0.7, 0.3}}], 0.5]",
                   "{{0, 1}, {1, 0}}", 0);
    assert_eval_eq("Head[MorphologicalComponents[{{1, 2}}]]", "MorphologicalComponents", 0);
    assert_eval_eq("Head[MorphologicalComponents[Image[{{1., 0.}}],"
                   " CornerNeighbors -> Maybe]]", "MorphologicalComponents", 0);
    assert_eval_eq("MemberQ[Attributes[MorphologicalComponents], Protected]", "True", 0);
}

static void test_median_removes_an_outlier_exactly(void) {
    /* THE discriminating property, and why a median exists alongside a Gaussian: an isolated outlier
     * is removed EXACTLY, not attenuated. One bright pixel in a constant field leaves NOTHING behind
     * -- the whole result is the background value -- where the mean leaves a visible bump at the same
     * place. Both halves are asserted, because "the median smooths" is true of the mean too and
     * proves nothing. */
    assert_eval_eq("Union[Flatten[ImageData[MedianFilter["
                   "Image[Table[If[x == 3 && y == 3, 1., 0.25], {y, 5}, {x, 5}]], 1]]]]",
                   "{0.25}", 0);
    assert_eval_eq("Part[ImageData[MeanFilter["
                   "Image[Table[If[x == 3 && y == 3, 1., 0.25], {y, 5}, {x, 5}]], 1]], 3, 3]"
                   " > 0.25", "True", 0);
    /* A median of a constant is that constant, exactly -- and 0.375 is exactly representable. */
    assert_eval_eq("Union[Flatten[ImageData[MedianFilter["
                   "Image[Table[0.375, {4}, {4}]], 1]]]]", "{0.375}", 0);
    /* Radius 0 is a one-element window, hence the identity. */
    assert_eval_eq("ImageData[MedianFilter[" MIMG ", 0]] === ImageData[" MIMG "]", "True", 0);
}

static void test_median_is_not_separable(void) {
    /* The median is the ONE operator in this file that does not decompose, and this pins that the
     * implementation does the real thing rather than the fast wrong thing.
     *
     * A sum, a maximum and a minimum all separate because they ignore grouping. A median depends on a
     * value's RANK within the whole window, and grouping destroys rank. On {{1,2,9},{3,4,5},{6,7,8}}
     * the true median of all nine values is 5, while the median of the row medians {2,4,7} is 4. A
     * separable implementation would return 4 here -- fast, plausible, and wrong. */
    assert_eval_eq("Part[ImageData[MedianFilter[Image[{{1.,2.,9.},{3.,4.,5.},{6.,7.,8.}}], 1]],"
                   " 2, 2] == 5.", "True", 0);
    /* The output of a rank filter is always one of its INPUTS -- an even window takes the lower
     * middle rather than averaging the two, which would invent a value not present in the window. */
    assert_eval_eq("Module[{v = Flatten[ImageData[MedianFilter["
                   "Image[{{0.25, 0.5}, {0.75, 1.}}], 1]]], src = {0.25, 0.5, 0.75, 1.}},"
                   " And @@ Map[MemberQ[src, #] &, v]]", "True", 0);
}

static void test_meanfilter_is_a_box_convolution(void) {
    /* MeanFilter IS a convolution with a normalised box, and is implemented as one -- so this asserts
     * the identity rather than a reimplementation agreeing with itself. Two implementations of one
     * identity is how the identity quietly stops holding. */
    assert_eval_eq("Module[{img = Image[Table[N[Mod[x*3 + y*5, 7]/7.], {y, 6}, {x, 7}]]},"
                   " Chop[Max[Abs[Flatten[ImageData[MeanFilter[img, 1]]"
                   " - ImageData[ImageConvolve[img, Table[1./9., {3}, {3}]]]]]]]]", "0", 0);
    /* A constant survives a mean exactly, borders included, because the box sums to 1 and the padding
     * replicates. */
    assert_eval_eq("Chop[Max[Abs[Flatten[ImageData[MeanFilter["
                   "Image[Table[0.25, {5}, {5}]], 2]]] - 0.25]]]", "0", 0);
    /* Both rank filters bracket nothing in general, but a mean must lie between the min and the max
     * of the image -- a weak property, and the one that fails if the kernel is not normalised. */
    assert_eval_eq("Module[{f = Flatten[ImageData[" MIMG "]],"
                   " m = Flatten[ImageData[MeanFilter[" MIMG ", 1]]]},"
                   " And @@ Map[(Min[f] <= # <= Max[f]) &, m]]", "True", 0);
    assert_eval_eq("Head[MedianFilter[" MIMG ", -1]]", "MedianFilter", 0);
    assert_eval_eq("Head[MedianFilter[" MIMG ", 1.5]]", "MedianFilter", 0);
    assert_eval_eq("Head[MeanFilter[{{1, 2}}, 1]]", "MeanFilter", 0);
    assert_eval_eq("And @@ Map[MemberQ[Attributes[#], Protected] &,"
                   " {MedianFilter, MeanFilter}]", "True", 0);
}

#define DTONE "Image[Table[If[x == 1 && y == 1, 0., 1.], {y, 6}, {x, 6}]]"

static void test_distance_transform_is_exactly_euclidean(void) {
    /* THE row separating an exact transform from the classic chamfer approximation, and it is an
     * EQUALITY because a 3-4-5 triangle has an integer hypotenuse. One background pixel at the
     * corner; the pixel three across and four down must read exactly 5. A chamfer transform gives
     * about 5.03 -- invisible on any picture, and the reason a tolerance here would have hidden the
     * wrong algorithm entirely. */
    assert_eval_eq("Part[ImageData[DistanceTransform[" DTONE "]], 5, 4] == 5.", "True", 0);
    /* A background pixel is exactly 0, and its immediate neighbour exactly 1. */
    assert_eval_eq("Part[ImageData[DistanceTransform[" DTONE "]], 1, 1] == 0.", "True", 0);
    assert_eval_eq("Part[ImageData[DistanceTransform[" DTONE "]], 1, 2] == 1.", "True", 0);
    /* The diagonal neighbour is Sqrt[2], which is what a chamfer transform cannot represent with
     * integer steps -- the root of its inexactness. */
    assert_eval_eq("Chop[Part[ImageData[DistanceTransform[" DTONE "]], 2, 2] - Sqrt[2.]]",
                   "0", 0);
    /* And a 5-12-13 triangle, so the 3-4-5 case cannot be passing by coincidence. */
    assert_eval_eq("Module[{one = Image[Table[If[x == 1 && y == 1, 0., 1.], {y, 14}, {x, 14}]]},"
                   " Part[ImageData[DistanceTransform[one]], 13, 6] == 13.]", "True", 0);
}

static void test_distance_transform_degenerate_cases(void) {
    /* All background: every distance is 0, and nothing divides by an empty parabola stack. */
    assert_eval_eq("Union[Flatten[ImageData[DistanceTransform["
                   "Image[Table[0., {4}, {4}]]]]]]", "{0.0}", 0);
    /* A checkerboard has every foreground pixel adjacent to a background one, so the only distances
     * are 0 and 1 -- which also confirms the transform is over the FOUR-neighbour distance to the
     * nearest zero rather than something looser. */
    assert_eval_eq("Union[Flatten[ImageData[DistanceTransform["
                   "Image[Table[If[Mod[x + y, 2] == 0, 1., 0.], {y, 4}, {x, 4}]]]]]]",
                   "{0.0, 1.0}", 0);
    /* The threshold argument chooses the foreground: at t = 0.5 a 0.3 pixel is background, so its
     * own distance is 0. */
    assert_eval_eq("Part[ImageData[DistanceTransform["
                   "Image[{{0.3, 0.7}, {0.7, 0.7}}], 0.5]], 1, 1] == 0.", "True", 0);
    assert_eval_eq("Head[DistanceTransform[{{1, 2}}]]", "DistanceTransform", 0);
    assert_eval_eq("Head[DistanceTransform[" DTONE ", x]]", "DistanceTransform", 0);
    assert_eval_eq("MemberQ[Attributes[DistanceTransform], Protected]", "True", 0);
}

static void test_distance_transform_relates_to_morphology(void) {
    /* An absolute relationship between two features rather than a property of either alone, and the
     * OFF-BY-ONE in it is the interesting part. A 3x3 erosion removes every pixel that has a
     * background pixel in its 8-neighbourhood -- that is, every pixel at Chebyshev distance 1 -- so
     * surviving ONE erosion means distance >= 2, and surviving k means distance >= k + 1, not >= k.
     * The first version of this row asserted >= 2 against TWO erosions and failed, correctly.
     *
     * The Euclidean-to-Chebyshev step is exact rather than approximate: Chebyshev 1 bounds Euclidean
     * by Sqrt[2] < 2, and Euclidean >= Chebyshev always, so "Euclidean >= 2" and "Chebyshev >= 2"
     * pick out the same pixels. Checked as an equality over the whole image, not a spot check.
     *
     * SameQ rather than Equal, and for a reason worth recording because it is NOT about booleans.
     * `False == False` reduces to True perfectly well. But `(2.0 >= 2.) == (1.0 > 0.)` evaluates to
     * `2.0 == True` -- the parenthesised comparison is being swallowed into a chained comparison
     * rather than evaluated first. That is a comparison bug independent of anything here, filed to
     * look at separately; SameQ sidesteps it.
     *
     * The route to it is the lesson: the assertion failed while the underlying arrays matched
     * exactly, and printing the two arrays is what settled it. Two attempts at re-deriving the claim
     * had already gone wrong before that -- once genuinely (the off-by-one) and once not (the claim
     * was right and the test was broken). Print the data before arguing with it a third time. */
    assert_eval_eq("Module[{img = Image[Table[If[2 <= x <= 7 && 2 <= y <= 7, 1., 0.],"
                   " {y, 9}, {x, 9}]], dt, er},"
                   " dt = ImageData[DistanceTransform[img]];"
                   " er = ImageData[Erosion[img, 1]];"
                   " And @@ MapThread[((#1 >= 2.) === (#2 > 0.)) &, {Flatten[dt], Flatten[er]}]]",
                   "True", 0);
    /* And two erosions correspond to distance >= 3, confirming the k + 1 relation rather than a
     * coincidence at k = 1. */
    assert_eval_eq("Module[{img = Image[Table[If[2 <= x <= 7 && 2 <= y <= 7, 1., 0.],"
                   " {y, 9}, {x, 9}]], dt, er},"
                   " dt = ImageData[DistanceTransform[img]];"
                   " er = ImageData[Erosion[Erosion[img, 1], 1]];"
                   " And @@ MapThread[((#1 >= 3.) === (#2 > 0.)) &, {Flatten[dt], Flatten[er]}]]",
                   "True", 0);
}

static void test_volume_downsampling_does_not_alias(void) {
    /* The 3-D version of the checkerboard row. A 2x2x2-periodic pattern halved on every axis is
     * sampled at exactly the frequency that annihilates it: area averaging returns its mean, 0.5
     * everywhere, and nearest returns a FLAT FIELD. Exact values on both sides. */
    assert_eval_eq("Union[Flatten[ImageData[ImageResize["
                   "Image3D[Table[N[Mod[x + y + z, 2]], {z, 4}, {y, 4}, {x, 4}]], {2, 2, 2}]]]]",
                   "{0.5}", 0);
    assert_eval_eq("Union[Flatten[ImageData[ImageResize["
                   "Image3D[Table[N[Mod[x + y + z, 2]], {z, 4}, {y, 4}, {x, 4}]], {2, 2, 2},"
                   " Resampling -> \"Nearest\"]]]]", "{0.0}", 0);
    /* An exact block mean over EIGHT voxels: a 2x2x2 volume of slice values 1 and 2 averages to
     * exactly 1.5. */
    assert_eval_eq("ImageData[ImageResize[Image3D[Table[N[z], {z, 2}, {y, 2}, {x, 2}]],"
                   " {1, 1, 1}]]", "{{{1.5}}}", 0);
    /* Area averaging preserves the mean at an integer factor, since every source voxel carries
     * equal total weight -- a conservation law in three dimensions as in two. */
    assert_eval_eq("Module[{big = Image3D[Table[N[Mod[x*7 + y*13 + z*3, 251]/251],"
                   " {z, 8}, {y, 8}, {x, 8}]]},"
                   " Chop[Mean[Flatten[ImageData[ImageResize[big, {4, 4, 4}]]]]"
                   " - Mean[Flatten[ImageData[big]]]]]", "0", 0);
}

static void test_volume_resize_respects_the_axis_order(void) {
    /* THE reversal, and it needs a non-cubic volume AND a non-cubic TARGET. The spec is
     * {width, height, depth} while the storage is depth x height x width, so a resize has to reverse
     * the spec before indexing. A cubic source or a cubic target validates none of that -- and with
     * three axes there are six ways to get it wrong, each of which a cube would hide.
     *
     * Source is 4 wide, 3 high, 2 deep; asked for 2 wide, 6 high, 4 deep. Both the reported
     * dimensions and the data shape are asserted, since either alone would pass with axes swapped. */
    assert_eval_eq("ImageDimensions[ImageResize["
                   "Image3D[Table[N[(x + 10 y + 100 z)/1000.], {z, 2}, {y, 3}, {x, 4}]],"
                   " {2, 6, 4}]]", "{2, 6, 4}", 0);
    assert_eval_eq("Dimensions[ImageData[ImageResize["
                   "Image3D[Table[N[(x + 10 y + 100 z)/1000.], {z, 2}, {y, 3}, {x, 4}]],"
                   " {2, 6, 4}]]]", "{4, 6, 2}", 0);
    /* Every method is the identity at 1:1, which is what fails on a half-voxel shift in any axis --
     * though only the aliasing and block-mean rows above can catch a shift at other scales. */
    assert_eval_eq("Module[{v = " VOL234 "},"
                   " {ImageData[ImageResize[v, {4, 3, 2}, Resampling -> \"Average\"]]"
                   "    === ImageData[v],"
                   "  ImageData[ImageResize[v, {4, 3, 2}, Resampling -> \"Bilinear\"]]"
                   "    === ImageData[v],"
                   "  ImageData[ImageResize[v, {4, 3, 2}, Resampling -> \"Nearest\"]]"
                   "    === ImageData[v]}]", "{True, True, True}", 0);
    /* A constant volume survives enlargement exactly. */
    assert_eval_eq("Union[Flatten[ImageData[ImageResize["
                   "Image3D[Table[0.25, {2}, {2}, {2}]], {4, 4, 4}]]]]", "{0.25}", 0);
    /* A two-element spec on a volume declines rather than being guessed at as a plane. */
    assert_eval_eq("Head[ImageResize[" VOL234 ", {2, 3}]]", "ImageResize", 0);
    assert_eval_eq("Head[ImageResize[" VOL234 ", {2, 3, 0}]]", "ImageResize", 0);
    assert_eval_eq("Head[ImageResize[" VOL234 ", {2, 3, 1.5}]]", "ImageResize", 0);
    /* Colour volumes keep their channels through a resize. */
    assert_eval_eq("Module[{cv = Image3D[Table[{0.5, 0.25, 0.75}, {2}, {3}, {4}]]},"
                   " {ImageChannels[ImageResize[cv, {2, 2, 1}]],"
                   "  Dimensions[ImageData[ImageResize[cv, {2, 2, 1}]]]}]",
                   "{3, {1, 2, 2, 3}}", 0);
}

#define LIMG "Image[Table[N[Mod[x*3 + y*5, 7]/7.], {y, 6}, {x, 8}]]"

static void test_levels_counts_sum_to_the_pixel_count(void) {
    /* THE property for a histogram, and it is exact: every pixel lands in exactly one bin, so a total
     * that disagrees means a bin boundary is wrong or a pixel was dropped. 6 x 8 is 48. */
    assert_eval_eq("Total[Map[Last, ImageLevels[" LIMG "]]]", "48", 0);
    assert_eval_eq("Total[Map[Last, ImageLevels[" LIMG ", 8]]]", "48", 0);
    assert_eval_eq("{Length[ImageLevels[" LIMG "]], Length[ImageLevels[" LIMG ", 8]]}",
                   "{256, 8}", 0);
    /* A Bit image has TWO natural levels, not 256 -- those are its distinct values, and binning it
     * into anything else would invent structure. Both counts are pinned. */
    assert_eval_eq("ImageLevels[Image[{{0, 1}, {1, 0}}]]", "{{0.0, 2}, {1.0, 2}}", 0);
    /* Levels are on ImageData's unit scale, so the last level is exactly 1. */
    assert_eval_eq("First[Last[ImageLevels[" LIMG "]]] == 1.", "True", 0);
    /* Volumes count too -- 2 x 3 x 4 is 24 voxels. */
    assert_eval_eq("Total[Map[Last, ImageLevels["
                   "Image3D[Table[N[z/4.], {z, 2}, {y, 3}, {x, 4}]]]]]", "24", 0);
    assert_eval_eq("Head[ImageLevels[" LIMG ", 1]]", "ImageLevels", 0);
}

static void test_imageadjust_stretches_exactly_and_is_idempotent(void) {
    /* The stretch puts the darkest pixel on exactly 0 and the brightest on exactly 1 -- equalities,
     * not tolerances. */
    assert_eval_eq("Module[{a = ImageData[ImageAdjust[Image[{{0.25, 0.5}, {0.5, 0.75}}]]]},"
                   " {Min[Flatten[a]] == 0., Max[Flatten[a]] == 1.}]", "{True, True}", 0);
    /* IDEMPOTENT, and this is the row that would catch an off-by-one in the range: a slightly wrong
     * divisor still produces a plausible-looking contrast curve, but a second stretch would then move
     * the pixels again. After a correct stretch the second pass is the identity. */
    assert_eval_eq("Module[{a = ImageData[ImageAdjust[Image[{{0.25, 0.5}, {0.5, 0.75}}]]]},"
                   " ImageData[ImageAdjust[Image[a]]] === a]", "True", 0);
    /* A constant image has no range: unchanged, rather than dividing by zero or mapping the single
     * value to an arbitrary end. */
    assert_eval_eq("Union[Flatten[ImageData[ImageAdjust[Image[Table[0.4, {3}, {3}]]]]]]",
                   "{0.4}", 0);
    /* Monotone -- a stretch may not reorder pixels. */
    assert_eval_eq("OrderedQ[Flatten[ImageData[ImageAdjust[Image[{{0.1, 0.2, 0.3, 0.9}}]]]]]",
                   "True", 0);
}

static void test_imageadjust_parametric_curve(void) {
    /* c = 0, b = 0, g = 1 must be the EXACT identity: contrast pivots about mid-grey so zero contrast
     * leaves the pivot arithmetic a no-op, and gamma 1 is skipped. If any step were misordered this
     * would drift. */
    assert_eval_eq("Chop[Max[Abs[Flatten[ImageData[ImageAdjust[" LIMG ", {0., 0., 1.}]]"
                   " - ImageData[" LIMG "]]]]]", "0", 0);
    /* Brightness is a plain offset: 0.5 + 0.25 is exactly 0.75. */
    assert_eval_eq("Part[ImageData[ImageAdjust[Image[{{0.5}}], {0., 0.25}]], 1, 1] == 0.75",
                   "True", 0);
    /* Contrast pivots about mid-grey, so mid-grey itself is FIXED by any contrast change -- that is
     * what "pivots" means, and it is the property that fails if the pivot is 0 instead of 1/2. */
    assert_eval_eq("Part[ImageData[ImageAdjust[Image[{{0.5}}], {3., 0.}]], 1, 1] == 0.5",
                   "True", 0);
    /* Clipping happens before gamma, so a value driven past 1 by brightness stays 1 rather than
     * becoming a power of something out of range. */
    assert_eval_eq("Part[ImageData[ImageAdjust[Image[{{0.9}}], {0., 0.5, 2.}]], 1, 1] == 1.",
                   "True", 0);
    assert_eval_eq("Head[ImageAdjust[" LIMG ", {1.}]]", "ImageAdjust", 0);
    assert_eval_eq("Head[ImageAdjust[" LIMG ", {0., 0., -1.}]]", "ImageAdjust", 0);
    assert_eval_eq("Head[ImageAdjust[{{1, 2}}]]", "ImageAdjust", 0);
    /* Volumes go through both forms. */
    assert_eval_eq("ImageDimensions[ImageAdjust["
                   "Image3D[Table[N[z/4.], {z, 2}, {y, 3}, {x, 4}]]]]", "{4, 3, 2}", 0);
    assert_eval_eq("And @@ Map[MemberQ[Attributes[#], Protected] &, {ImageLevels, ImageAdjust}]",
                   "True", 0);
}

static void test_vanherk_agrees_with_an_independent_reference(void) {
    /* The 1-D max/min is now van Herk--Gil-Werman: three comparisons per pixel whatever the radius.
     * It is EXACT rather than approximate, because max and min are associative and idempotent, so
     * splitting a window at a block boundary and recombining loses nothing -- which means the fast
     * path must agree BIT-EXACTLY with a direct computation.
     *
     * The reference is written in Mathilda and shares no code with the C implementation: replicate
     * padding means an out-of-range read clamps to the edge, so the neighbourhood max equals the max
     * over the CLAMPED index range, which Span with Max/Min expresses directly.
     *
     * Comparing against Dilation with an all-ones matrix would NOT have been a check: an all-ones
     * element is a full rectangle, so both sides would take van Herk and agree with themselves. */
    assert_eval_eq("Module[{img = Image[Table[N[Mod[x*11 + y*7, 23]/23.], {y, 9}, {x, 11}]], d, ref},"
                   " d = ImageData[img];"
                   " ref[dd_, r_] := Module[{h = Length[dd], w = Length[First[dd]]},"
                   "   Table[Max[Flatten[dd[[Max[1, y - r] ;; Min[h, y + r],"
                   "                        Max[1, x - r] ;; Min[w, x + r]]]]], {y, h}, {x, w}]];"
                   " And @@ Table[ImageData[Dilation[img, r]] === ref[d, r], {r, {1, 2, 3, 4}}]]",
                   "True", 0);
    assert_eval_eq("Module[{img = Image[Table[N[Mod[x*11 + y*7, 23]/23.], {y, 9}, {x, 11}]], d, ref},"
                   " d = ImageData[img];"
                   " ref[dd_, r_] := Module[{h = Length[dd], w = Length[First[dd]]},"
                   "   Table[Min[Flatten[dd[[Max[1, y - r] ;; Min[h, y + r],"
                   "                        Max[1, x - r] ;; Min[w, x + r]]]]], {y, h}, {x, w}]];"
                   " And @@ Table[ImageData[Erosion[img, r]] === ref[d, r], {r, {1, 2, 3, 4}}]]",
                   "True", 0);
    /* A radius LARGER than the image is the case the block decomposition can get wrong, since the
     * last block is short and a window then spans the whole padded line. */
    assert_eval_eq("Module[{img = Image[{{0.25, 0.75}, {0.5, 1.}}]},"
                   " {Union[Flatten[ImageData[Dilation[img, 4]]]],"
                   "  Union[Flatten[ImageData[Erosion[img, 4]]]]}]", "{{1.0}, {0.25}}", 0);
}

static void test_correlation_is_convolution_reflected(void) {
    /* The reflection asserted from the OTHER side: on a delta with {{1,2,3}}, correlation gives
     * {3,2,1} where convolution gives {1,2,3}. Both are pinned, so neither can drift toward the
     * other. */
    assert_eval_eq("ImageData[ImageCorrelate[Image[{{0., 1., 0.}}], {{1, 2, 3}}]]",
                   "{{3.0, 2.0, 1.0}}", 0);
    /* THE identity: correlation equals convolution with the kernel reversed on both axes -- and it is
     * BIT-EXACT, not merely close, because correlation is now implemented that way rather than as a
     * second loop. The identity holds by construction instead of by two implementations agreeing.
     *
     * It was a tolerance first, when correlation had its own loop: the difference was 3.6e-15, the
     * same products summed in a different order. Deriving it from convolution removed the difference
     * along with the duplicate code, so the assertion tightened from Chop to ===. */
    assert_eval_eq("Module[{img = Image[Table[N[Mod[x*3 + y*5, 7]/7.], {y, 6}, {x, 7}]],"
                   " k = {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}}},"
                   " ImageData[ImageCorrelate[img, k]]"
                   " === ImageData[ImageConvolve[img, Reverse[Reverse[k], 2]]]]", "True", 0);
    /* On a SYMMETRIC kernel the two must agree -- again to rounding, since the Gaussian is rank 1 and
     * so takes the separable path on the convolution side. */
    assert_eval_eq("Module[{img = Image[Table[N[Mod[x*3 + y*5, 7]/7.], {y, 6}, {x, 7}]]},"
                   " Chop[Max[Abs[Flatten[ImageData[ImageCorrelate[img, GaussianMatrix[1]]]"
                   " - ImageData[ImageConvolve[img, GaussianMatrix[1]]]]]]]]", "0", 0);
}

static void test_ncc_scores_exactly_one_at_an_exact_match(void) {
    /* NCC's defining properties, and the bound is the strong one: |NCC| <= 1 by Cauchy-Schwarz, with
     * equality exactly when the window is an affine image of the template. So where the template is a
     * CROP of the image the score is exactly 1, and nothing anywhere exceeds 1.
     *
     * NOTE WHAT IS *NOT* ASSERTED: that the peak is unique. The first version of this test expected a
     * single maximum at the crop's location and failed -- because the test image was periodic
     * (Mod[x*11 + y*7, 23]) and had five exact matches. Uniqueness is a property of the IMAGE, not of
     * NCC, so asserting it would have been testing the wrong thing. The claim is that an exact match
     * scores exactly 1 and that 1 is the ceiling. */
    assert_eval_eq("Module[{big = Image[Table[N[Mod[x*11 + y*7, 23]/23.], {y, 12}, {x, 12}]],"
                   " bd, tmpl, nc},"
                   " bd = ImageData[big];"
                   " tmpl = bd[[5 ;; 7, 6 ;; 8]];"
                   " nc = ImageData[ImageCorrelate[big, tmpl, \"NormalizedCrossCorrelation\"]];"
                   " {Chop[Part[nc, 6, 7] - 1.], Chop[Max[Flatten[nc]] - 1.],"
                   "  Min[Flatten[nc]] >= -1.000001}]", "{0, 0, True}", 0);
    /* Brightness and contrast invariance is the whole reason for normalising: the same template scaled
     * and offset still scores 1, where plain correlation would rank by brightness instead. */
    assert_eval_eq("Module[{big = Image[Table[N[Mod[x*11 + y*7, 23]/23.], {y, 12}, {x, 12}]],"
                   " bd, tmpl, nc},"
                   " bd = ImageData[big];"
                   " tmpl = Map[(0.4 # + 0.1) &, bd[[5 ;; 7, 6 ;; 8]], {2}];"
                   " nc = ImageData[ImageCorrelate[big, tmpl, \"NormalizedCrossCorrelation\"]];"
                   " Chop[Part[nc, 6, 7] - 1.]]", "0", 0);
    /* A flat window has no shape to compare, so it scores 0 rather than dividing by zero. Scoring 1
     * would be worse than wrong: every flat region would then match every template perfectly. */
    assert_eval_eq("Union[Flatten[ImageData[ImageCorrelate[Image[Table[0.5, {5}, {5}]],"
                   " {{1., 2.}, {3., 4.}}, \"NormalizedCrossCorrelation\"]]]]", "{0.0}", 0);
    assert_eval_eq("Head[ImageCorrelate[" LIMG ", {{1., 2.}}, \"Bogus\"]]", "ImageCorrelate", 0);
    assert_eval_eq("Head[ImageCorrelate[{{1, 2}}, {{1.}}]]", "ImageCorrelate", 0);
    assert_eval_eq("MemberQ[Attributes[ImageCorrelate], Protected]", "True", 0);
}

#define RIMG "Image[{{0.1, 0.2, 0.3, 0.4}, {0.5, 0.6, 0.7, 0.8}}]"

static void test_right_angle_rotation_is_exact(void) {
    /* A quarter turn is a pure index permutation -- every pixel lands on another pixel's position, so
     * nothing is interpolated and FOUR of them are exactly the identity. Asserted with ===, which is
     * only available because right angles are kept off the resampler; routed through bilinear the
     * identity would become approximate for no reason. */
    assert_eval_eq("ImageData[ImageRotate[ImageRotate[ImageRotate[ImageRotate[" RIMG "]]]]]"
                   " === ImageData[" RIMG "]", "True", 0);
    /* A quarter turn SWAPS the dimensions, which a square image could not show. 4 wide by 2 high
     * becomes 2 wide by 4 high. */
    assert_eval_eq("{ImageDimensions[" RIMG "], ImageDimensions[ImageRotate[" RIMG "]]}",
                   "{{4, 2}, {2, 4}}", 0);
    /* Two quarter turns are exactly a half turn, and a half turn of a single row reverses it. */
    assert_eval_eq("ImageData[ImageRotate[ImageRotate[" RIMG "]]]"
                   " === ImageData[ImageRotate[" RIMG ", Pi]]", "True", 0);
    assert_eval_eq("ImageData[ImageRotate[Image[{{0.1, 0.2, 0.3}}], Pi]]",
                   "{{0.3, 0.2, 0.1}}", 0);
    /* A full turn is the image itself. Pi and 2 Pi are EXACT SYMBOLIC values, not machine reals, so
     * these rows also pin that the angle is numericalised -- ImageRotate[img, Pi] declined outright
     * before that was added, which is a poor answer to the most natural way of writing a half turn. */
    assert_eval_eq("ImageData[ImageRotate[" RIMG ", 2 Pi]] === ImageData[" RIMG "]", "True", 0);
    assert_eval_eq("ImageData[ImageRotate[" RIMG ", 90 Degree]] === ImageData[ImageRotate[" RIMG "]]",
                   "True", 0);
}

static void test_reflection_is_exact_and_self_inverse(void) {
    /* Also a pure permutation, so reflecting twice is exactly the identity -- on both axes. */
    assert_eval_eq("ImageData[ImageReflect[ImageReflect[" RIMG "]]] === ImageData[" RIMG "]",
                   "True", 0);
    assert_eval_eq("ImageData[ImageReflect[ImageReflect[" RIMG ", Left], Left]]"
                   " === ImageData[" RIMG "]", "True", 0);
    /* And the two directions do different things, which a self-inverse test alone cannot show. */
    assert_eval_eq("ImageData[ImageReflect[" RIMG "]]",
                   "{{0.5, 0.6, 0.7, 0.8}, {0.1, 0.2, 0.3, 0.4}}", 0);
    assert_eval_eq("ImageData[ImageReflect[" RIMG ", Left]]",
                   "{{0.4, 0.3, 0.2, 0.1}, {0.8, 0.7, 0.6, 0.5}}", 0);
    /* Dimensions are unchanged by a reflection, unlike a quarter turn. */
    assert_eval_eq("ImageDimensions[ImageReflect[" RIMG "]]", "{4, 2}", 0);
    assert_eval_eq("Head[ImageReflect[" RIMG ", Sideways]]", "ImageReflect", 0);
    assert_eval_eq("Head[ImageRotate[{{1, 2}}]]", "ImageRotate", 0);
    assert_eval_eq("And @@ Map[MemberQ[Attributes[#], Protected] &, {ImageRotate, ImageReflect}]",
                   "True", 0);
}

static void test_arbitrary_angle_rotation_round_trips_smooth_content(void) {
    /* An arbitrary angle CANNOT be exact -- a rotated pixel grid does not land on a pixel grid -- so
     * the property is that theta then -theta recovers the interior, and it holds to 1.1e-16 on SMOOTH
     * content. The interior specifically, because the corners rotate out of frame and back as
     * background.
     *
     * SMOOTH content specifically, and that distinction is the point. The first version of this row
     * used the noise-like Mod[x*7 + y*13, 251] pattern and recovered only to 0.377 -- not a bug but
     * two bilinear interpolations destroying energy at the Nyquist limit, which is what interpolation
     * does. A tolerance loose enough to admit that would have accepted an outright wrong rotation
     * too; using band-limited content instead keeps the assertion tight enough to mean something. */
    assert_eval_eq("Module[{sm = Image[Table[N[(x + 2 y)/128.], {y, 32}, {x, 32}]], rt, d0},"
                   " rt = ImageData[ImageRotate[ImageRotate[sm, 0.3], -0.3]];"
                   " d0 = ImageData[sm];"
                   " Chop[Max[Abs[Flatten[rt[[9 ;; 24, 9 ;; 24]] - d0[[9 ;; 24, 9 ;; 24]]]]]]]",
                   "0", 0);
    /* Out-of-frame area reads as 0, not the replicated edge: a rotation exposes area that was never
     * photographed, and smearing the border across it would invent content. A half turn of a
     * one-pixel-tall image by a small angle pulls background into the corners. */
    assert_eval_eq("Module[{r = ImageData[ImageRotate[Image[Table[1., {8}, {8}]], 0.4]]},"
                   " Min[Flatten[r]] < 0.5]", "True", 0);
    /* Dimensions are preserved by an arbitrary-angle rotation, unlike a quarter turn. */
    assert_eval_eq("ImageDimensions[ImageRotate[" RIMG ", 0.3]]", "{4, 2}", 0);
}


/* ImagePad and ImageCrop -- the geometry pair whose composition is exactly the identity. */
static void test_pad_then_crop_is_the_identity(void) {
    /* THE property. Both are index arithmetic with no interpolation, so cropping back to the
     * original size after a symmetric pad returns the original pixels bit for bit. It is also the
     * test that catches an off-by-one on either side INDEPENDENTLY: a pad that adds one row too many
     * at the top and one too few at the bottom still has the right total size, and only a round trip
     * notices. */
    assert_eval_eq("Module[{img = Image[Table[N[Sin[i*0.7] + Cos[j*0.4]], {i, 5}, {j, 7}]]},"
                   " ImageCrop[ImagePad[img, 3], ImageDimensions[img]] === img]", "True", 0);
    /* Three channels too, since the pad indexes pixels and the channel stride is the easy thing to
     * drop. */
    assert_eval_eq("Module[{img = Image[Table[N[i + j + k]/20, {i, 3}, {j, 4}, {k, 3}]]},"
                   " ImageCrop[ImagePad[img, 2], ImageDimensions[img]] === img]", "True", 0);
    /* A CENTRED crop can only invert a SYMMETRIC pad -- with pt=4, pb=3 the content starts one row
     * below where a centred crop begins looking. The exact inverse of an asymmetric pad is a
     * negative pad, and that is the property worth asserting. */
    assert_eval_eq("Module[{img = Image[Table[N[i*j]/30, {i, 5}, {j, 7}]]},"
                   " ImagePad[ImagePad[img, {{1, 2}, {3, 4}}], {{-1, -2}, {-3, -4}}] === img]",
                   "True", 0);
    assert_eval_eq("ImageDimensions[ImagePad[Image[Table[N[i*j]/30, {i, 5}, {j, 7}]], 3]]",
                   "{13, 11}", 0);
    /* Negative padding IS cropping, and agrees with the centred crop that removes the same border. */
    assert_eval_eq("Module[{img = Image[Table[N[i*j]/30, {i, 5}, {j, 7}]]},"
                   " ImagePad[img, -1] === ImageCrop[img, {5, 3}]]", "True", 0);
    /* Padding may not erase the image, and a crop may not enlarge it: both decline. */
    assert_eval_eq("Head[ImagePad[Image[Table[N[i*j]/30, {i, 5}, {j, 7}]], -9]]", "ImagePad", 0);
    assert_eval_eq("Head[ImageCrop[Image[Table[N[i*j]/30, {i, 5}, {j, 7}]], {99, 99}]]",
                   "ImageCrop", 0);
}

static void test_pad_side_convention_and_modes(void) {
    /* THE CONVENTION TRAP. {{left, right}, {bottom, top}} names the pair in VISUAL order, while the
     * data's first row is the TOP -- so `top` padding adds rows at the START of the array, the
     * reverse of how the spec reads. A symmetric pad cannot tell; this pads one side only. */
    assert_eval_eq("Module[{img = Image[Table[N[i*j]/30, {i, 5}, {j, 7}]], q},"
                   " q = ImagePad[img, {{0, 0}, {0, 2}}];"
                   " {ImageDimensions[q], First[ImageData[q]] === Table[0., {7}],"
                   "  Chop[ImageData[q][[3]] - First[ImageData[img]]] === Table[0, {7}]}]",
                   "{{7, 7}, True, True}", 0);
    /* The three fill modes, on a row where each gives a different answer. Reflection does NOT repeat
     * the edge pixel: doubling the edge sample would bias any later average toward the border. */
    assert_eval_eq("First[ImageData[ImagePad[Image[{{0.1, 0.2, 0.3}}], {{1, 1}, {0, 0}}, "
                   "\"Reflected\"]]]", "{0.2, 0.1, 0.2, 0.3, 0.2}", 0);
    assert_eval_eq("First[ImageData[ImagePad[Image[{{0.1, 0.2, 0.3}}], {{1, 1}, {0, 0}}, "
                   "\"Fixed\"]]]", "{0.1, 0.1, 0.2, 0.3, 0.3}", 0);
    assert_eval_eq("First[ImageData[ImagePad[Image[{{0.1, 0.2, 0.3}}], {{1, 1}, {0, 0}}, 0.9]]]",
                   "{0.9, 0.1, 0.2, 0.3, 0.9}", 0);
    /* Reflection deeper than the image works, because the period is 2n-2 rather than one width. */
    assert_eval_eq("First[ImageData[ImagePad[Image[{{0.1, 0.2, 0.3}}], {{4, 0}, {0, 0}}, "
                   "\"Reflected\"]]]", "{0.1, 0.2, 0.3, 0.2, 0.1, 0.2, 0.3}", 0);
    /* Every mode is invertible by a negative pad, since none of them touches the interior. */
    assert_eval_eq("Module[{img = Image[Table[N[i*j]/30, {i, 5}, {j, 7}]]},"
                   " {ImagePad[ImagePad[img, 1, \"Fixed\"], -1] === img,"
                   "  ImagePad[ImagePad[img, 2, \"Reflected\"], -2] === img}]",
                   "{True, True}", 0);
}

static void test_crop_trims_a_uniform_border(void) {
    /* ImageCrop with no size asks how much of the frame carries no information. */
    assert_eval_eq("Module[{b = Image[{{0., 0., 0., 0., 0.}, {0., 0.4, 0.6, 0., 0.},"
                   " {0., 0., 0., 0., 0.}}]},"
                   " {ImageDimensions[ImageCrop[b]], ImageData[ImageCrop[b]]}]",
                   "{{2, 1}, {{0.4, 0.6}}}", 0);
    /* The border colour is read from a corner rather than assumed black -- a scanned page's margin
     * is white, and assuming black would trim nothing at all. */
    assert_eval_eq("Module[{w = Image[{{1., 1., 1.}, {1., 0.2, 1.}, {1., 1., 1.}}]},"
                   " {ImageDimensions[ImageCrop[w]], ImageData[ImageCrop[w]]}]",
                   "{{1, 1}, {{0.2}}}", 0);
    /* An entirely uniform image has no content to keep, and a zero-sized image is not an image, so
     * it comes back unchanged. */
    assert_eval_eq("Module[{u = Image[Table[0.5, {3}, {4}]]}, ImageCrop[u] === u]", "True", 0);
}


/* Volumetric pad and crop. Every property the planar pair has is asserted again at rank 3, because
 * "same as the other rank" is a claim to verify: the volumetric paths here have twice dropped
 * something the planar ones had. */
#define VOL3 "Image3D[Table[N[100 z + 10 y + x], {z, 1, 3}, {y, 1, 4}, {x, 1, 5}]]"

static void test_volume_pad_then_crop_is_the_identity(void) {
    /* NON-CUBIC extents on purpose -- 3 slices of 4 rows of 5 columns -- so any axis swap in the
     * pad or the crop shows up as a wrong size or wrong voxels rather than passing by symmetry. */
    assert_eval_eq("ImageDimensions[" VOL3 "]", "{5, 4, 3}", 0);
    assert_eval_eq("ImageDimensions[ImagePad[" VOL3 ", 2]]", "{9, 8, 7}", 0);
    assert_eval_eq("Module[{v = " VOL3 "},"
                   " ImageCrop[ImagePad[v, 2], ImageDimensions[v]] === v]", "True", 0);
    /* Negative padding crops, and agrees with the centred crop that removes the same shell. */
    assert_eval_eq("Module[{v = " VOL3 "},"
                   " {ImageDimensions[ImagePad[v, -1]], ImagePad[v, -1] === ImageCrop[v, {3, 2, 1}]}]",
                   "{{3, 2, 1}, True}", 0);
    /* An asymmetric pad is inverted by the negated pad, not by a centred crop. */
    assert_eval_eq("Module[{v = " VOL3 "},"
                   " ImagePad[ImagePad[v, {{1, 2}, {3, 0}, {0, 1}}],"
                   "          {{-1, -2}, {-3, 0}, {0, -1}}] === v]", "True", 0);
    /* Padding may not erase the volume; a crop may not enlarge it; and trimming a uniform border in
     * three dimensions is a question this declines rather than answering silently. */
    assert_eval_eq("Head[ImagePad[" VOL3 ", -9]]", "ImagePad", 0);
    assert_eval_eq("Head[ImageCrop[" VOL3 ", {99, 99, 99}]]", "ImageCrop", 0);
    assert_eval_eq("Head[ImageCrop[" VOL3 "]]", "ImageCrop", 0);
}

static void test_volume_pad_axis_conventions(void) {
    /* THE THREE CONVENTIONS, pinned one axis at a time, since a symmetric pad cannot tell any of
     * them apart. The spec is {{left, right}, {bottom, top}, {first slice, last slice}} in the order
     * ImageDimensions reports, and the HEIGHT pair is the reversed one -- Mathematica names it
     * {bottom, top} while row 1 is the top of the image. */
    assert_eval_eq("Module[{v = " VOL3 ", q},"
                   " q = ImagePad[v, {{0, 0}, {0, 0}, {0, 2}}];"
                   " {ImageDimensions[q], First[ImageData[q]] === First[ImageData[v]],"
                   "  ImageData[q][[5]] === Table[0., {4}, {5}]}]",
                   "{{5, 4, 5}, True, True}", 0);
    assert_eval_eq("Module[{v = " VOL3 ", r},"
                   " r = ImagePad[v, {{0, 0}, {0, 2}, {0, 0}}];"
                   " {ImageDimensions[r], ImageData[r][[1, 1]] === Table[0., {5}],"
                   "  ImageData[r][[1, 3]] === ImageData[v][[1, 1]]}]",
                   "{{5, 6, 3}, True, True}", 0);
    assert_eval_eq("Module[{v = " VOL3 ", u},"
                   " u = ImagePad[v, {{2, 0}, {0, 0}, {0, 0}}];"
                   " {ImageDimensions[u], ImageData[u][[1, 1, 1 ;; 2]] === {0., 0.}}]",
                   "{{7, 4, 3}, True}", 0);
}

static void test_volume_pad_fill_modes(void) {
    /* Exact voxel values, which is the only way to tell the modes apart: the volume holds
     * 100z + 10y + x, so slice 1 row 1 is {111, 112, 113, 114, 115}. */
    assert_eval_eq("ImageData[ImagePad[" VOL3 ", {{1, 0}, {0, 0}, {0, 0}}, \"Fixed\"]][[1, 1, 1 ;; 3]]",
                   "{111.0, 111.0, 112.0}", 0);
    /* Reflection does NOT repeat the edge voxel: padding by 2 gives 113, 112 before 111. */
    assert_eval_eq("ImageData[ImagePad[" VOL3 ", {{2, 0}, {0, 0}, {0, 0}}, "
                   "\"Reflected\"]][[1, 1, 1 ;; 4]]", "{113.0, 112.0, 111.0, 112.0}", 0);
    /* The modes apply per axis, not only across rows: Fixed in DEPTH repeats the first slice. */
    assert_eval_eq("Module[{v = " VOL3 "},"
                   " First[ImageData[ImagePad[v, {{0, 0}, {0, 0}, {1, 0}}, \"Fixed\"]]]"
                   " === First[ImageData[v]]]", "True", 0);
    assert_eval_eq("ImageData[ImagePad[" VOL3 ", {{1, 0}, {0, 0}, {0, 0}}, 0.5]][[1, 1, 1]]",
                   "0.5", 0);
    /* Every mode is invertible by a negative pad, since none touches the interior. */
    assert_eval_eq("Module[{v = " VOL3 "},"
                   " {ImagePad[ImagePad[v, 1, \"Fixed\"], -1] === v,"
                   "  ImagePad[ImagePad[v, 2, \"Reflected\"], -2] === v}]", "{True, True}", 0);
}


/* NCC on summed-area tables. The properties that matter are the argmax (exact, an integer) and
 * agreement with the definition; the peak is 1.0 only to a few ulp now, which is the documented
 * trade for taking the window statistics out of the inner loop. */
#define NIMG "Image[Table[N[Mod[i*13 + j*7, 97]]/97, {i, 1, 12}, {j, 1, 14}]]"

static void test_ncc_finds_the_patch_it_was_given(void) {
    /* A patch OF THE IMAGE, so the correct answer is known: the score peaks where the patch came
     * from. The argmax is asserted exactly -- it is the property template matching rests on. */
    assert_eval_eq("Module[{img = " NIMG ", t, r},"
                   " t = ImageData[img][[4 ;; 6, 5 ;; 7]];"
                   " r = ImageData[ImageCorrelate[img, t, \"NormalizedCrossCorrelation\"]];"
                   " Position[r, Max[Flatten[r]]]]", "{{5, 6}}", 0);
    /* And the peak is 1 to within 1e-12. Not exactly 1: the numerator comes from the correlation and
     * the variance from the tables, so they no longer share a summation order. */
    assert_eval_eq("Module[{img = " NIMG ", t, r},"
                   " t = ImageData[img][[4 ;; 6, 5 ;; 7]];"
                   " r = ImageData[ImageCorrelate[img, t, \"NormalizedCrossCorrelation\"]];"
                   " Abs[1 - Max[Flatten[r]]] < 1.*^-12]", "True", 0);
    /* Every score is a cosine of two deviation vectors, so it cannot leave [-1, 1]. */
    assert_eval_eq("Module[{img = " NIMG ", t, r},"
                   " t = ImageData[img][[4 ;; 6, 5 ;; 7]];"
                   " r = Flatten[ImageData[ImageCorrelate[img, t, "
                   "     \"NormalizedCrossCorrelation\"]]];"
                   " Min[r] >= -1.000000001 && Max[r] <= 1.000000001]", "True", 0);
}

static void test_ncc_agrees_with_its_definition(void) {
    /* The reference is the definition itself, written out in Mathilda with the same clamped border:
     * the tables are an optimisation, and an optimisation is only correct if it agrees with the thing
     * it replaced. A summed-area table that is off by a row still produces plausible scores. */
    assert_eval_eq("Module[{img = " NIMG ", d, t, r, ref, cl},"
                   " d = ImageData[img]; t = d[[4 ;; 6, 5 ;; 7]];"
                   " r = ImageData[ImageCorrelate[img, t, \"NormalizedCrossCorrelation\"]];"
                   " cl[v_, n_] := Max[1, Min[n, v]];"
                   " ref = Table[Module[{win, wb, tb, nu, sn, tn},"
                   "   win = Table[d[[cl[y + i - 2, 12], cl[x + j - 2, 14]]], {i, 3}, {j, 3}];"
                   "   wb = Mean[Flatten[win]]; tb = Mean[Flatten[t]];"
                   "   nu = Total[Flatten[(win - wb)*(t - tb)]];"
                   "   sn = Sqrt[Total[Flatten[(win - wb)^2]]];"
                   "   tn = Sqrt[Total[Flatten[(t - tb)^2]]];"
                   "   If[sn > 0 && tn > 0, nu/(sn*tn), 0.]], {y, 12}, {x, 14}];"
                   " Max[Abs[Flatten[r - ref]]] < 1.*^-11]", "True", 0);
}

static void test_ncc_degenerate_and_invariant(void) {
    /* A uniform image has no variance at any position. Every score is 0 -- and specifically not NaN,
     * which is what a sqrt of a slightly-negative variance would give and what would then spread
     * through Max and Position without ever looking like an error. */
    assert_eval_eq("Union[Flatten[ImageData[ImageCorrelate[Image[Table[0.5, {8}, {8}]],"
                   " {{1., 2.}, {3., 4.}}, \"NormalizedCrossCorrelation\"]]]]", "{0.0}", 0);
    /* NCC compares shape, so an affine change of the template must not move any score. */
    assert_eval_eq("Module[{img = " NIMG ", t, r},"
                   " t = ImageData[img][[4 ;; 6, 5 ;; 7]];"
                   " r = ImageData[ImageCorrelate[img, t, \"NormalizedCrossCorrelation\"]];"
                   " Max[Abs[Flatten[ImageData[ImageCorrelate[img, 2 t + 3,"
                   "   \"NormalizedCrossCorrelation\"]] - r]]] < 1.*^-11]", "True", 0);
    /* An unknown method name declines rather than silently correlating. */
    assert_eval_eq("Head[ImageCorrelate[" NIMG ", {{1.}}, \"NoSuchMethod\"]]", "ImageCorrelate", 0);
}


/* Convolution through the transform. None of the groups above reaches it: every kernel they use is
 * either separable or smaller than the crossover, so before these the FFT path shipped untested. */
#define FIMG "Image[Table[N[Mod[i*13 + j*7, 97]]/97, {i, 1, 40}, {j, 1, 48}]]"
/* The definition, longhand, with the clamped border and ci = Floor[kh/2] -- the same convention the
 * direct loop uses. Written once as a Mathilda function the tests below call.
 *
 * The first version of this reference was WRONG, off by one in both axes, and it accused the working
 * separable path of a 376-unit error before it accused the new one. Which is the argument for having
 * it: a reference that disagrees with two independent implementations is the thing that is broken, and
 * that is only visible when there ARE two. */
#define FREF \
  "cl[v_, n_] := Max[1, Min[n, v]];" \
  "refOf[kk_, kh_, kw_] := Module[{ci = Floor[kh/2], cj = Floor[kw/2]}," \
  " Table[Sum[d[[cl[y - i + ci, 40], cl[x - j + cj, 48]]] * kk[[i + 1, j + 1]]," \
  "       {i, 0, kh - 1}, {j, 0, kw - 1}], {y, 1, 40}, {x, 1, 48}]];"
#define FKER(n) "Table[N[Mod[i*5 + j*3, 11]] - 5, {i, 1, " n "}, {j, 1, " n "}]"

static void test_fft_convolution_agrees_with_the_definition(void) {
    /* Straddling the crossover: 6x6 is still the dense loop, 7x7 and up take the transform, and both
     * must answer the same question. The kernel is deliberately full-rank so the separable path
     * cannot claim it. */
    assert_eval_eq("Module[{d = ImageData[" FIMG "], k}," FREF
                   " Max[Table[k = " FKER("n") ";"
                   "   Max[Abs[Flatten[ImageData[ImageConvolve[" FIMG ", k]] - refOf[k, n, n]]]],"
                   "   {n, {6, 7, 9, 12}}]] < 1.*^-12]", "True", 0);
    /* A NON-SQUARE kernel, both orientations. The transform has a different length on each axis, and
     * a transposed pad or a transposed output crop is exactly the bug that would survive a square
     * test. */
    assert_eval_eq("Module[{d = ImageData[" FIMG "], k}," FREF
                   " k = Table[N[Mod[i*5 + j*3, 11]] - 5, {i, 1, 7}, {j, 1, 11}];"
                   " Max[Abs[Flatten[ImageData[ImageConvolve[" FIMG ", k]] - refOf[k, 7, 11]]]]"
                   " < 1.*^-12]", "True", 0);
    assert_eval_eq("Module[{d = ImageData[" FIMG "], k}," FREF
                   " k = Table[N[Mod[i*5 + j*3, 11]] - 5, {i, 1, 11}, {j, 1, 7}];"
                   " Max[Abs[Flatten[ImageData[ImageConvolve[" FIMG ", k]] - refOf[k, 11, 7]]]]"
                   " < 1.*^-12]", "True", 0);
    /* An EVEN extent, where Floor[kh/2] puts the centre off-centre and the pad offset has to follow
     * it rather than assuming symmetry. */
    assert_eval_eq("Module[{d = ImageData[" FIMG "], k}," FREF
                   " k = " FKER("8") ";"
                   " Max[Abs[Flatten[ImageData[ImageConvolve[" FIMG ", k]] - refOf[k, 8, 8]]]]"
                   " < 1.*^-12]", "True", 0);
}

static void test_fft_convolution_across_channels(void) {
    /* The kernel's transform is computed once and reused for every channel, so a colour image is
     * where a stale reused buffer would appear -- as identical channels, which is why the assertion
     * is that they DIFFER rather than that the call returned something. */
    assert_eval_eq("Module[{img, o},"
                   " img = Image[Table[N[Mod[i*13 + j*7 + k*5, 97]]/97,"
                   "   {i, 1, 20}, {j, 1, 24}, {k, 1, 3}]];"
                   " o = ImageData[ImageConvolve[img, " FKER("9") "]];"
                   " {Dimensions[o], o[[1, 1, 1]] != o[[1, 1, 2]],"
                   "  o[[1, 1, 2]] != o[[1, 1, 3]]}]", "{{20, 24, 3}, True, True}", 0);
    /* Below the crossover nothing changed, and the identity kernel is still EXACT -- the transform
     * would have cost that, which is one reason the switch is where it is. */
    assert_eval_eq("Module[{img = " FIMG "},"
                   " ImageData[ImageConvolve[img, {{1}}]] === ImageData[img]]", "True", 0);
}


/* The rank-3 transform. It matters more than the planar one: a kd*kh*kw kernel is CUBIC in the
 * radius, so 9x9x9 is 729 taps per voxel where 9x9 is 81 per pixel. */
#define F3VOL "Image3D[Table[N[Mod[z*13 + y*7 + x*3, 97]]/97, {z, 1, 10}, {y, 1, 12}, {x, 1, 14}]]"
#define F3REF \
  "cl[a_, n_] := Max[1, Min[n, a]];" \
  "refOf[kk_, kd_, kh_, kw_] := Module[{cz = Floor[kd/2], cy = Floor[kh/2], cx = Floor[kw/2]}," \
  " Table[Sum[dd[[cl[z - m + cz, 10], cl[y - i + cy, 12], cl[x - j + cx, 14]]]" \
  "           * kk[[m + 1, i + 1, j + 1]]," \
  "       {m, 0, kd - 1}, {i, 0, kh - 1}, {j, 0, kw - 1}], {z, 1, 10}, {y, 1, 12}, {x, 1, 14}]];" \
  "mk[a_, b_, cc_] := Table[N[Mod[m*5 + i*3 + j*2, 11]] - 5, {m, 1, a}, {i, 1, b}, {j, 1, cc}];"

static void test_volume_fft_convolution_agrees_with_the_definition(void) {
    /* Straddling the crossover, which at rank 3 sits at 3x3x3: that size stays dense (the two paths
     * measured equal there, so staying dense avoids the transform's scratch memory) and 5x5x5 up take
     * the transform. Both must answer the same question. */
    assert_eval_eq("Module[{dd = ImageData[" F3VOL "], k}," F3REF
                   " Max[Table[k = mk[n, n, n];"
                   "   Max[Abs[Flatten[ImageData[ImageConvolve[" F3VOL ", k]]"
                   "                   - refOf[k, n, n, n]]]], {n, {3, 5, 7}}]] < 1.*^-12]",
                   "True", 0);
    /* NON-CUBIC extents with all three axes distinct, in both orderings. The transform has a
     * different length on each axis, so a transposed pad or a transposed output window is the bug
     * that a cubic kernel could never reveal -- and the volumetric paths in this file have twice
     * shipped exactly that class of mistake. */
    assert_eval_eq("Module[{dd = ImageData[" F3VOL "], k}," F3REF
                   " k = mk[3, 4, 5];"
                   " Max[Abs[Flatten[ImageData[ImageConvolve[" F3VOL ", k]]"
                   "                 - refOf[k, 3, 4, 5]]]] < 1.*^-12]", "True", 0);
    assert_eval_eq("Module[{dd = ImageData[" F3VOL "], k}," F3REF
                   " k = mk[5, 4, 3];"
                   " Max[Abs[Flatten[ImageData[ImageConvolve[" F3VOL ", k]]"
                   "                 - refOf[k, 5, 4, 3]]]] < 1.*^-12]", "True", 0);
    /* EVEN extents on every axis, where Floor[k/2] is off-centre and the pad offset must follow. */
    assert_eval_eq("Module[{dd = ImageData[" F3VOL "], k}," F3REF
                   " k = mk[4, 4, 4];"
                   " Max[Abs[Flatten[ImageData[ImageConvolve[" F3VOL ", k]]"
                   "                 - refOf[k, 4, 4, 4]]]] < 1.*^-12]", "True", 0);
}

static void test_volume_fft_convolution_across_channels(void) {
    /* A colour VOLUME. The kernel transform is hoisted out of the channel loop here too, so a stale
     * reused buffer would show as identical channels. */
    assert_eval_eq("Module[{cv, k, o},"
                   " cv = Image3D[Table[N[Mod[z*13 + y*7 + x*3 + ch*29, 97]]/97,"
                   "   {z, 1, 6}, {y, 1, 8}, {x, 1, 10}, {ch, 1, 3}]];"
                   " k = Table[N[Mod[m*5 + i*3 + j*2, 11]] - 5, {m, 1, 5}, {i, 1, 5}, {j, 1, 5}];"
                   " o = ImageData[ImageConvolve[cv, k]];"
                   " {Dimensions[o], o[[1, 1, 1, 1]] != o[[1, 1, 1, 2]],"
                   "  o[[1, 1, 1, 2]] != o[[1, 1, 1, 3]]}]", "{{6, 8, 10, 3}, True, True}", 0);
    /* A separable volumetric kernel still takes the factorisation, which beats any transform:
     * kd + kh + kw taps against a transform of the whole padded extent. */
    assert_eval_eq("Module[{v = " F3VOL ", u, k, o},"
                   " u = Table[N[i], {i, 1, 5}];"
                   " k = Table[u[[m]] u[[i]] u[[j]], {m, 1, 5}, {i, 1, 5}, {j, 1, 5}];"
                   " o = ImageData[ImageConvolve[v, k]];"
                   " Dimensions[o]]", "{10, 12, 14}", 0);
}


/* Corner detection. The properties here are absolute, not scores: an edge must be EXACTLY zero, a
 * uniform region exactly zero, and the response must rotate with the image. */
#define CK24 "Image[Table[If[Mod[Quotient[i - 1, 6] + Quotient[j - 1, 6], 2] == 0, 0., 1.]," \
             " {i, 1, 24}, {j, 1, 24}]]"

static void test_corner_response_is_zero_where_there_is_no_corner(void) {
    /* No gradient at all: exactly zero, and not merely small. */
    assert_eval_eq("Union[Flatten[ImageData[CornerFilter[Image[Table[0.5, {16}, {16}]]]]]]",
                   "{0.0}", 0);
    assert_eval_eq("ImageCorners[Image[Table[0.5, {16}, {16}]]]", "{}", 0);
    /* THE DISCRIMINATING PROPERTY. Along a straight edge every gradient in the window is parallel, so
     * the structure tensor has rank 1 and both det and lambda_min vanish. A corner detector that
     * fires on an edge is an edge detector, and no visual check of a response map reveals it. */
    assert_eval_eq("Module[{e = Image[Table[If[j <= 8, 0., 1.], {i, 1, 24}, {j, 1, 24}]]},"
                   " {Max[Flatten[ImageData[CornerFilter[e]]]],"
                   "  Max[Flatten[ImageData[CornerFilter[e, 2, \"Harris\"]]]],"
                   "  Length[ImageCorners[e]]}]", "{0.0, 0.0, 0}", 0);
    /* The tensor is positive semi-definite, so the smaller eigenvalue is never negative -- which is
     * also why it is computed from (Sxx-Syy)^2 + 4 Sxy^2 rather than trace^2 - 4 det, where
     * cancellation could produce a negative discriminant and a NaN. */
    assert_eval_eq("Min[Flatten[ImageData[CornerFilter[" CK24 "]]]] >= -1.*^-15", "True", 0);
}

static void test_corner_response_finds_corners_and_rotates_with_them(void) {
    /* One square corner, at a known place. */
    assert_eval_eq("Module[{sq, rs},"
                   " sq = Image[Table[If[i >= 8 && j >= 8, 1., 0.], {i, 1, 24}, {j, 1, 24}]];"
                   " rs = ImageData[CornerFilter[sq]];"
                   " {Position[rs, Max[Flatten[rs]]], ImageCorners[sq]}]",
                   "{{{8, 8}}, {{8, 8}}}", 0);
    /* A 24x24 checkerboard of 6-pixel blocks has interior corners at every multiple of 6, so 3x3 of
     * them -- an exact count, not an approximate one. */
    assert_eval_eq("Length[ImageCorners[" CK24 "]]", "9", 0);
    /* ROTATIONAL COVARIANCE. A quarter turn is an exact index permutation, so the response of the
     * rotated image must be the rotation of the response. This is what catches a swapped Ix/Iy, a
     * transposed smoothing pass, or a sign error in one derivative -- all of which leave a response
     * map that looks entirely plausible. */
    assert_eval_eq("Module[{ck = " CK24 ", r0},"
                   " r0 = CornerFilter[ck];"
                   " Max[Abs[Flatten[ImageData[CornerFilter[ImageRotate[ck]]]"
                   "                 - ImageData[ImageRotate[r0]]]]] < 1.*^-15]", "True", 0);
    assert_eval_eq("Module[{ck = " CK24 "},"
                   " CornerFilter[Nest[ImageRotate, ck, 4]] === CornerFilter[ck]]", "True", 0);
}

static void test_corner_options_are_wired(void) {
    /* Both parameters have to change the answer, or they are decoration. The radius is the window the
     * tensor averages over and the method is what is done with the eigenvalues. */
    assert_eval_eq("Module[{ck = " CK24 "},"
                   " {Max[Abs[Flatten[ImageData[CornerFilter[ck, 1]]"
                   "                  - ImageData[CornerFilter[ck, 4]]]]] > 1.*^-6,"
                   "  Max[Abs[Flatten[ImageData[CornerFilter[ck, 2, \"Harris\"]]"
                   "                  - ImageData[CornerFilter[ck, 2, \"MinimumEigenvalue\"]]]]]"
                   "  > 1.*^-6}]", "{True, True}", 0);
    /* A position indexes ImageData directly -- {row, column}, which is the documented convention. */
    assert_eval_eq("Module[{sq},"
                   " sq = Image[Table[If[i >= 8 && j >= 8, 1., 0.], {i, 1, 24}, {j, 1, 24}]];"
                   " Head[ImageData[sq][[Sequence @@ First[ImageCorners[sq]]]]]]", "Real", 0);
    /* A raised threshold keeps no more corners than a lower one. */
    assert_eval_eq("Module[{ck = " CK24 "},"
                   " Length[ImageCorners[ck, 2, 0.9]] <= Length[ImageCorners[ck, 2, 0.01]]]",
                   "True", 0);
    assert_eval_eq("Head[CornerFilter[" CK24 ", 2, \"Nope\"]]", "CornerFilter", 0);
    assert_eval_eq("Head[CornerFilter[" CK24 ", 0]]", "CornerFilter", 0);
}


/* ImageCorners' separation and feature limit. The properties here are absolute: a minimum separation
 * d means no two returned positions are within d, which is a statement about the output alone. */
#define CIMG "Image[Table[N[Mod[i*7 + j*13, 251]]/251, {i, 1, 128}, {j, 1, 128}]]"
/* Smallest pairwise distance in a list of {row, column} positions; Infinity for fewer than two. */
#define MIND \
  "minD[l_] := If[Length[l] < 2, Infinity," \
  " Min[Table[Sqrt[N[(l[[a, 1]] - l[[b, 1]])^2 + (l[[a, 2]] - l[[b, 2]])^2]]," \
  "           {a, 1, Length[l]}, {b, a + 1, Length[l]}]]];"

static void test_corner_separation_is_respected(void) {
    /* THE property, at three separations. Not "roughly spread out" -- no pair closer than d. */
    assert_eval_eq("Module[{img = " CIMG "}," MIND
                   " And @@ Table[minD[ImageCorners[img, 2, 0.05, d]] >= d, {d, {2., 5., 10.}}]]",
                   "True", 0);
    /* Asking for more separation cannot return more corners. */
    assert_eval_eq("Module[{img = " CIMG "},"
                   " Length[ImageCorners[img, 2, 0.05, 2.]]"
                   " >= Length[ImageCorners[img, 2, 0.05, 10.]]]", "True", 0);
    /* Separation removes a great deal on a busy image: the first two filters leave clusters a pixel
     * or two apart, which is what makes the raw list unusable as a feature set. */
    assert_eval_eq("Module[{img = " CIMG "},"
                   " Length[ImageCorners[img, 2, 0.05, 10.]] < Length[ImageCorners[img]]]",
                   "True", 0);
    /* The same image must always give the same list -- ties are broken by position for exactly
     * this reason. */
    assert_eval_eq("Module[{img = " CIMG "},"
                   " ImageCorners[img, 2, 0.05, 5.] === ImageCorners[img, 2, 0.05, 5.]]", "True", 0);
}

static void test_corner_feature_limit_keeps_the_strongest(void) {
    /* Exactly n, and they are the strongest n. The second half is asserted as "the weakest kept is at
     * least as strong as the strongest dropped", which is what "the top n" means -- and NOT by
     * comparing against Sort[..., Greater], which cannot serve as an oracle here: Mathilda's real
     * comparison is tolerant in Mathematica's way, so two responses a single ulp apart are Equal and
     * neither is Greater, and Sort then permutes them freely. */
    assert_eval_eq("Module[{img = " CIMG ", rs, all, top, kept, rest},"
                   " rs = ImageData[CornerFilter[img]];"
                   " all = ImageCorners[img]; top = ImageCorners[img, 2, 0.05, 0, 5];"
                   " kept = Table[rs[[p[[1]], p[[2]]]], {p, top}];"
                   " rest = Table[rs[[p[[1]], p[[2]]]], {p, Drop[all, 5]}];"
                   " {Length[top], Min[kept] >= Max[rest]}]", "{5, True}", 0);
    /* Strongest first, stated as the list being non-increasing. */
    assert_eval_eq("Module[{img = " CIMG ", rs, v},"
                   " rs = ImageData[CornerFilter[img]];"
                   " v = Table[rs[[p[[1]], p[[2]]]], {p, Take[ImageCorners[img], 12]}];"
                   " And @@ Table[v[[i]] >= v[[i + 1]], {i, 1, Length[v] - 1}]]", "True", 0);
    /* The limit is applied AFTER separation, so a limited request still separates: applied first it
     * would return n positions out of a single cluster. */
    assert_eval_eq("Module[{img = " CIMG ", l}," MIND
                   " l = ImageCorners[img, 2, 0.05, 10., 3];"
                   " {Length[l] <= 3, minD[l] >= 10.}]", "{True, True}", 0);
    /* The default path is unchanged by any of this. */
    assert_eval_eq("Length[ImageCorners[" CK24 "]]", "9", 0);
    assert_eval_eq("Head[ImageCorners[" CK24 ", 2, 0.05, -1.]]", "ImageCorners", 0);
    assert_eval_eq("Head[ImageCorners[" CK24 ", 2, 0.05, 0, 0]]", "ImageCorners", 0);
}


/* Volumetric corner detection. The eigenvalue HIERARCHY is the test, and it is what distinguishes a
 * real 3-D detector from a 2-D one applied slicewise: in three dimensions the rank of the structure
 * tensor says what the neighbourhood holds -- 1 is a plane, 2 is an edge, 3 is a corner -- so
 * lambda_min is exactly zero for BOTH a plane and an edge, and positive only at a corner. */
#define V16(cond) "Image3D[Table[If[" cond ", 0., 1.], {z, 1, 16}, {y, 1, 16}, {x, 1, 16}]]"
#define VPLANE  V16("x <= 6")
#define VEDGE   V16("x <= 8 && y <= 8")
#define VCORNER V16("x <= 8 && y <= 8 && z <= 8")
#define VCK "Image3D[Table[If[Mod[Quotient[z - 1, 6] + Quotient[y - 1, 6] + Quotient[x - 1, 6], 2]" \
            " == 0, 0., 1.], {z, 1, 24}, {y, 1, 24}, {x, 1, 24}]]"

static void test_volume_corner_eigenvalue_hierarchy(void) {
    /* Flat: exactly zero, not merely small. */
    assert_eval_eq("Union[Flatten[ImageData[CornerFilter[Image3D[Table[0.5, {12}, {12}, {12}]]]]]]",
                   "{0.0}", 0);
    /* A PLANE is rank 1 and an EDGE is rank 2, so both have lambda_min = 0. A detector that fires on
     * either has not been written for three dimensions, and a response map cannot show the
     * difference. */
    assert_eval_eq("Max[Flatten[ImageData[CornerFilter[" VPLANE "]]]] < 1.*^-15", "True", 0);
    assert_eval_eq("Max[Flatten[ImageData[CornerFilter[" VEDGE "]]]] < 1.*^-15", "True", 0);
    /* Only rank 3 responds, and it responds at the octant's corner. */
    assert_eval_eq("Module[{rs = ImageData[CornerFilter[" VCORNER "]]},"
                   " {Max[Flatten[rs]] > 1.*^-3, Position[rs, Max[Flatten[rs]]]}]",
                   "{True, {{7, 7, 7}}}", 0);
    /* The tensor is positive semi-definite, so lambda_min is never negative. */
    assert_eval_eq("Min[Flatten[ImageData[CornerFilter[" VCK "]]]] >= 0.", "True", 0);
    /* Harris in three dimensions subtracts k*trace^3, not trace^2: det scales as lambda^3 here, and
     * the two terms must share a dimension or the constant means nothing. On a plane det is 0, so the
     * response is negative. */
    assert_eval_eq("Max[Flatten[ImageData[CornerFilter[" VPLANE ", 2, \"Harris\"]]]] <= 1.*^-18",
                   "True", 0);
    assert_eval_eq("ImageDimensions[CornerFilter[" VCORNER "]]", "{16, 16, 16}", 0);
}

static void test_volume_corner_axis_covariance(void) {
    /* Swapping two axes must swap them in the response. There is no ImageRotate for volumes, so the
     * permutation is done on the data -- the property is the same one the planar tests get from a
     * quarter turn, and it is what catches a transposed gradient, a transposed smoothing pass, or a
     * mis-indexed tensor entry, none of which spoil the look of the output. */
    assert_eval_eq("Module[{v, a},"
                   " v = Image3D[Table[N[Mod[z*13 + y*7 + x*3, 97]]/97,"
                   "   {z, 1, 10}, {y, 1, 12}, {x, 1, 14}]];"
                   " a = ImageData[CornerFilter[v]];"
                   " Max[Abs[Flatten[ImageData[CornerFilter[Image3D[Transpose[ImageData[v],"
                   "   {1, 3, 2}]]]] - Transpose[a, {1, 3, 2}]]]] < 1.*^-15]", "True", 0);
}

static void test_volume_imagecorners(void) {
    /* One corner, found once, as a triple that indexes ImageData. */
    assert_eval_eq("ImageCorners[" VCORNER "]", "{{7, 7, 7}}", 0);
    /* Neither a plane nor an edge is a corner -- the same hierarchy, now in the position list. */
    assert_eval_eq("{ImageCorners[" VPLANE "], Length[ImageCorners[" VEDGE "]]}", "{{}, 0}", 0);
    /* A 24^3 checkerboard of 6-voxel blocks has interior corners at every multiple of 6: 3^3 of
     * them, an exact count. */
    assert_eval_eq("Length[ImageCorners[" VCK "]]", "27", 0);
    assert_eval_eq("Length[First[ImageCorners[" VCK "]]]", "3", 0);
    assert_eval_eq("Head[ImageData[" VCK "][[Sequence @@ First[ImageCorners[" VCK "]]]]]", "Real", 0);
    /* Separation is Euclidean in three dimensions, and the same absolute property holds: no two
     * returned positions within d. */
    assert_eval_eq("Module[{ck = " VCK ", minD},"
                   " minD[l_] := If[Length[l] < 2, Infinity, Min[Table["
                   "   Sqrt[N[Total[(l[[a]] - l[[b]])^2]]], {a, 1, Length[l]}, {b, a + 1, Length[l]}]]];"
                   " And @@ Table[minD[ImageCorners[ck, 2, 0.05, dd]] >= dd, {dd, {4., 8.}}]]",
                   "True", 0);
    assert_eval_eq("Length[ImageCorners[" VCK ", 2, 0.05, 0, 3]]", "3", 0);
    assert_eval_eq("ImageCorners[" VCK ", 2, 0.05, 4.] === ImageCorners[" VCK ", 2, 0.05, 4.]",
                   "True", 0);
    /* The planar path is untouched by sharing the peak-finding code. */
    assert_eval_eq("Length[ImageCorners[" CK24 "]]", "9", 0);
}


/* Named options, through the generic reader in options.c. */
static void test_corner_named_options(void) {
    /* The registered defaults are what Options[head] reports, and they are the single source the
     * reader consults -- so this is not decoration, it is the mechanism. */
    /* expr_to_string quotes strings, unlike Print -- so the expected form carries the quotes. */
    assert_eval_eq("Options[CornerFilter]", "{Method -> \"MinimumEigenvalue\"}", 0);
    assert_eval_eq("Options[ImageCorners]", "{MaxFeatures -> Infinity}", 0);

    /* THE EQUIVALENCE THAT MATTERS, and the one that caught the bug: the option form and the
     * positional form must agree. The first version of the reader could not distinguish an option the
     * caller passed from a default it had filled in, so the default overrode the positional argument
     * and CornerFilter[img, 2, "Harris"] silently computed the other response. Both spellings of both
     * settings are asserted equal here. */
    assert_eval_eq("Module[{img = " CIMG "},"
                   " ImageCorners[img, MaxFeatures -> 5] === ImageCorners[img, 2, 0.05, 0, 5]]",
                   "True", 0);
    assert_eval_eq("Module[{ck = " CK24 "},"
                   " ImageData[CornerFilter[ck, 2, Method -> \"Harris\"]]"
                   " === ImageData[CornerFilter[ck, 2, \"Harris\"]]]", "True", 0);
    /* And the default still applies when nothing positional was given. */
    assert_eval_eq("Module[{ck = " CK24 "},"
                   " ImageData[CornerFilter[ck]]"
                   " === ImageData[CornerFilter[ck, 2, Method -> \"MinimumEigenvalue\"]]]",
                   "True", 0);

    /* Options are trailing, may accompany positional arguments, and the last duplicate wins. */
    assert_eval_eq("Module[{img = " CIMG "},"
                   " {Length[ImageCorners[img, 2, 0.05, 10., MaxFeatures -> 3]],"
                   "  Length[ImageCorners[img, MaxFeatures -> 9, MaxFeatures -> 2]]}]",
                   "{3, 2}", 0);
    /* An unknown option DECLINES rather than being ignored. Mathematica warns and continues; refusing
     * is the more conservative reading, and a typo'd option name that silently does nothing is the
     * failure this avoids. */
    assert_eval_eq("Module[{img = " CIMG "}, Head[ImageCorners[img, Nonsense -> 3]]]",
                   "ImageCorners", 0);
    assert_eval_eq("Module[{img = " CIMG "}, Head[ImageCorners[img, MaxFeatures -> -2]]]",
                   "ImageCorners", 0);
    assert_eval_eq("Head[CornerFilter[" CK24 ", Method -> \"Nope\"]]", "CornerFilter", 0);

    /* SetOptions works with no code in the builtin beyond reading the registry -- which is the
     * argument for the defaults living there rather than as constants in C. */
    assert_eval_eq("Module[{img = " CIMG ", a, b, c},"
                   " SetOptions[ImageCorners, MaxFeatures -> 4];"
                   " a = Length[ImageCorners[img]];"
                   " b = Length[ImageCorners[img, MaxFeatures -> 7]];"
                   " SetOptions[ImageCorners, MaxFeatures -> Infinity];"
                   " c = Length[ImageCorners[img]];"
                   " {a, b, c > 100}]", "{4, 7, True}", 0);

    /* Options reach the volumetric path too, since they are stripped before the rank dispatch. */
    assert_eval_eq("Length[ImageCorners[" VCK ", MaxFeatures -> 5]]", "5", 0);
}


/* LocalAdaptiveBinarize. A global threshold cannot binarize unevenly lit content -- not as a tuning
 * problem but in principle -- and the discriminating test below is what shows the feature earning its
 * place rather than merely running. */
#define LIMG "Image[Table[N[Mod[i*13 + j*7, 97]]/97, {i, 1, 14}, {j, 1, 16}]]"
/* The definition, longhand, with the same clamped border. */
#define LREF \
  "cl[v_, n_] := Max[1, Min[n, v]];" \
  "winOf[y_, x_, r_] := Flatten[Table[d[[cl[y + a, 14], cl[x + b, 16]]], {a, -r, r}, {b, -r, r}]];" \
  "muOf[y_, x_, r_] := Mean[winOf[y, x, r]];" \
  "sdOf[y_, x_, r_] := Sqrt[Max[0., Mean[winOf[y, x, r]^2] - muOf[y, x, r]^2]];"

static void test_local_adaptive_binarize_matches_the_definition(void) {
    /* With c2 or c1 moving the threshold off the pixel's own mean, agreement is EXACT. */
    assert_eval_eq("Module[{d = ImageData[" LIMG "]}," LREF
                   " ImageData[LocalAdaptiveBinarize[" LIMG ", 2, {1, -0.3, 0.02}]]"
                   " === N[Table[If[d[[y, x]] > muOf[y, x, 2] - 0.3 sdOf[y, x, 2] + 0.02, 1, 0],"
                   "        {y, 1, 14}, {x, 1, 16}]]]", "True", 0);
    assert_eval_eq("Module[{d = ImageData[" LIMG "]}," LREF
                   " ImageData[LocalAdaptiveBinarize[" LIMG ", 2, {0.9, 0, 0.}]]"
                   " === N[Table[If[d[[y, x]] > 0.9 muOf[y, x, 2], 1, 0],"
                   "        {y, 1, 14}, {x, 1, 16}]]]", "True", 0);

    /* MEAN-ONLY IS THE BOUNDARY CASE, and it is worth stating precisely rather than excluding. With
     * the threshold equal to the window mean, a pixel can tie it -- and on this periodic image many
     * do exactly. The summed-area mean (four table lookups, then sum/area) and a direct sum of the
     * window differ in the last bit, and a BINARY decision amplifies that into a visible flip: 15 of
     * 224 pixels here, every one of them a tie.
     *
     * So the assertion is not "the outputs are equal" but the stronger and true statement: every
     * disagreement is a pixel within one ulp of its own threshold. Where the definition is
     * numerically determined, the implementation matches it. */
    assert_eval_eq("Module[{d = ImageData[" LIMG "], mine, ref, bad}," LREF
                   " mine = ImageData[LocalAdaptiveBinarize[" LIMG ", 2]];"
                   " ref = N[Table[If[d[[y, x]] > muOf[y, x, 2], 1, 0], {y, 1, 14}, {x, 1, 16}]];"
                   " bad = Position[mine - ref, _?(# != 0 &)];"
                   " And @@ Table[Abs[d[[p[[1]], p[[2]]]] - muOf[p[[1]], p[[2]], 2]] < 1.*^-15,"
                   "              {p, bad}]]", "True", 0);
}

static void test_local_adaptive_binarize_properties(void) {
    /* A uniform image: the window mean IS the value, so nothing exceeds it -- all zero, exactly. */
    assert_eval_eq("Union[Flatten[ImageData[LocalAdaptiveBinarize[Image[Table[0.5, {16}, {16}]], 3]]]]",
                   "{0.0}", 0);
    /* Typed Bit, because the result is binary by construction. */
    assert_eval_eq("Part[LocalAdaptiveBinarize[" LIMG ", 3], 2]", "\"Bit\"", 0);
    assert_eval_eq("Union[Flatten[ImageData[LocalAdaptiveBinarize[" LIMG ", 3]]]]", "{0.0, 1.0}", 0);
    /* Raising the offset can only turn pixels off. */
    assert_eval_eq("Module[{img = " LIMG ", n},"
                   " n = Table[Total[Flatten[ImageData[LocalAdaptiveBinarize[img, 3, {1, 0, c}]]]],"
                   "      {c, {0., 0.05, 0.2}}];"
                   " And @@ Table[n[[i]] >= n[[i + 1]], {i, 1, 2}]]", "True", 0);
    assert_eval_eq("Head[LocalAdaptiveBinarize[" LIMG ", 0]]", "LocalAdaptiveBinarize", 0);
    assert_eval_eq("Head[LocalAdaptiveBinarize[" LIMG ", 3, {1, 0}]]", "LocalAdaptiveBinarize", 0);
}

static void test_local_adaptive_beats_global_under_a_ramp(void) {
    /* THE REASON THE FUNCTION EXISTS. A checkerboard under a strong lighting ramp: no single number
     * separates the two tones in both halves, so a global threshold must collapse one half to a
     * single value. The local form keeps the pattern in both. Asserted as set membership -- "this half
     * contains both values" -- which is absolute, rather than as a percentage recovered. */
    assert_eval_eq("Module[{ramp, gl, lo},"
                   " ramp = Image[Table[N[(0.1 + 0.8 (j - 1)/63)"
                   "   * If[Mod[Quotient[i - 1, 4] + Quotient[j - 1, 4], 2] == 0, 0.55, 1.]],"
                   "   {i, 1, 64}, {j, 1, 64}]];"
                   " gl = ImageData[Binarize[ramp]]; lo = ImageData[LocalAdaptiveBinarize[ramp, 5]];"
                   " {Length[Union[Flatten[Take[gl, All, 20]]]],"
                   "  Length[Union[Flatten[Take[lo, All, 20]]]],"
                   "  Length[Union[Flatten[Take[lo, All, -20]]]]}]", "{1, 2, 2}", 0);
}


/* Volumetric morphology. The properties are algebraic identities, which is the strongest thing a
 * morphological operator offers: they hold exactly or the implementation is wrong. */
#define MVOL "Image3D[Table[N[Mod[z*13 + y*7 + x*3, 97]]/97, {z, 1, 12}, {y, 1, 14}, {x, 1, 16}]]"

static void test_volume_morphology_algebra(void) {
    assert_eval_eq("ImageDimensions[Dilation[" MVOL ", 2]]", "{16, 14, 12}", 0);
    /* DUALITY: eroding f is dilating -f and negating. Checked through 1 - f so the values stay in
     * the unit interval, which costs one rounding and is why this is 1e-16 rather than exact. */
    assert_eval_eq("Module[{v = " MVOL ", dv, neg},"
                   " dv = ImageData[v]; neg = Image3D[1. - dv];"
                   " Max[Abs[Flatten[ImageData[Erosion[v, 2]]"
                   "   - (1. - ImageData[Dilation[neg, 2]])]]] < 1.*^-15]", "True", 0);
    /* ORDERING: opening never brightens and closing never darkens, pointwise. */
    assert_eval_eq("Module[{v = " MVOL ", dv},"
                   " dv = ImageData[v];"
                   " {Max[Flatten[ImageData[Opening[v, 2]] - dv]] <= 0.,"
                   "  Min[Flatten[ImageData[Closing[v, 2]] - dv]] >= 0.}]", "{True, True}", 0);
    /* IDEMPOTENCE, exactly. This is what makes an opening an opening rather than merely a smoother,
     * and it holds only because both passes use the SAME element. */
    assert_eval_eq("Module[{v = " MVOL "},"
                   " {ImageData[Opening[Opening[v, 2], 2]] === ImageData[Opening[v, 2]],"
                   "  ImageData[Closing[Closing[v, 2], 2]] === ImageData[Closing[v, 2]]}]",
                   "{True, True}", 0);
    /* Monotone in the radius, and r = 0 is exactly the identity. */
    assert_eval_eq("Module[{v = " MVOL "},"
                   " {Min[Flatten[ImageData[Dilation[v, 3]] - ImageData[Dilation[v, 2]]]] >= 0.,"
                   "  Max[Flatten[ImageData[Erosion[v, 3]] - ImageData[Erosion[v, 2]]]] <= 0.}]",
                   "{True, True}", 0);
    assert_eval_eq("Module[{v = " MVOL ", dv}, dv = ImageData[v];"
                   " {ImageData[Dilation[v, 0]] === dv, ImageData[Erosion[v, 0]] === dv}]",
                   "{True, True}", 0);
    /* A uniform volume has nothing to spread, so both leave it alone. */
    assert_eval_eq("Module[{u = Image3D[Table[0.5, {8}, {8}, {8}]]},"
                   " {Union[Flatten[ImageData[Dilation[u, 3]]]],"
                   "  Union[Flatten[ImageData[Erosion[u, 3]]]]}]", "{{0.5}, {0.5}}", 0);
}

static void test_volume_morphology_matches_the_definition(void) {
    /* Three van Herk passes must equal the max over the whole cube, EXACTLY -- max and min are
     * order statistics, so unlike a sum there is no rounding to hide behind and `===` is the right
     * comparison. */
    assert_eval_eq("Module[{v = " MVOL ", dv, cl, ref}, dv = ImageData[v];"
                   " cl[q_, n_] := Max[1, Min[n, q]];"
                   " ref = Table[Max[Flatten[Table[dv[[cl[z + a, 12], cl[y + b, 14], cl[x + c, 16]]],"
                   "   {a, -2, 2}, {b, -2, 2}, {c, -2, 2}]]], {z, 1, 12}, {y, 1, 14}, {x, 1, 16}];"
                   " ImageData[Dilation[v, 2]] === ref]", "True", 0);
    assert_eval_eq("Module[{v = " MVOL ", dv, cl, ref}, dv = ImageData[v];"
                   " cl[q_, n_] := Max[1, Min[n, q]];"
                   " ref = Table[Min[Flatten[Table[dv[[cl[z + a, 12], cl[y + b, 14], cl[x + c, 16]]],"
                   "   {a, -2, 2}, {b, -2, 2}, {c, -2, 2}]]], {z, 1, 12}, {y, 1, 14}, {x, 1, 16}];"
                   " ImageData[Erosion[v, 2]] === ref]", "True", 0);
    /* An explicit structuring element DECLINES at rank 3. An arbitrary 3-D element is not separable,
     * so honouring it would mean a cubic walk costing hundreds of times more behind the same
     * spelling -- declining says so instead of hiding it. */
    assert_eval_eq("Head[Dilation[" MVOL ", {{{1.}}}]]", "Dilation", 0);
}


/* Volumetric MeanFilter and MedianFilter. */
#define FVOL "Image3D[Table[N[Mod[z*13 + y*7 + x*3, 97]]/97, {z, 1, 10}, {y, 1, 12}, {x, 1, 14}]]"
#define FCUBE \
  "cl[q_, n_] := Max[1, Min[n, q]];" \
  "cube[y_, x_, z_, r_] := Flatten[Table[dv[[cl[z + a, 10], cl[y + b, 12], cl[x + c, 14]]]," \
  "  {a, -r, r}, {b, -r, r}, {c, -r, r}]];"

static void test_volume_mean_and_median_match_the_definition(void) {
    assert_eval_eq("ImageDimensions[MeanFilter[" FVOL ", 2]]", "{14, 12, 10}", 0);
    /* The mean, to 1e-12. Not exact, and the reason is stated in the source: the window total comes
     * from differencing two prefix sums, whose error scales with the LINE length rather than the
     * window's. That buys radius-independence -- 0.55 ms at r = 1 and at r = 4 alike -- for rounding
     * no caller of an image mean can observe. */
    assert_eval_eq("Module[{v = " FVOL ", dv}, dv = ImageData[v];" FCUBE
                   " Max[Abs[Flatten[ImageData[MeanFilter[v, 1]]"
                   "   - Table[Mean[cube[y, x, z, 1]], {z, 1, 10}, {y, 1, 12}, {x, 1, 14}]]]]"
                   " < 1.*^-12]", "True", 0);
    /* The median EXACTLY, at two radii. A median is an order statistic, so unlike a sum there is no
     * rounding to hide behind and `===` is the right comparison -- which also pins the even-window
     * convention, since the reference takes the same lower middle. */
    assert_eval_eq("Module[{v = " FVOL ", dv}, dv = ImageData[v];" FCUBE
                   " And @@ Table[ImageData[MedianFilter[v, r]] ==="
                   "   Table[Module[{s = Sort[cube[y, x, z, r]]}, s[[(Length[s] + 1)/2]]],"
                   "     {z, 1, 10}, {y, 1, 12}, {x, 1, 14}], {r, {1, 2}}]]", "True", 0);
    /* r = 0 is EXACTLY the identity for both. For the mean this is short-circuited on purpose: the
     * prefix-sum path returns the input only to within the prefix's rounding, and this identity was
     * silently lost when that path replaced the separable convolution. */
    assert_eval_eq("Module[{v = " FVOL ", dv}, dv = ImageData[v];"
                   " {ImageData[MeanFilter[v, 0]] === dv, ImageData[MedianFilter[v, 0]] === dv}]",
                   "{True, True}", 0);
    assert_eval_eq("Module[{u = Image3D[Table[0.5, {8}, {8}, {8}]]},"
                   " {Union[Flatten[ImageData[MeanFilter[u, 2]]]],"
                   "  Union[Flatten[ImageData[MedianFilter[u, 2]]]]}]", "{{0.5}, {0.5}}", 0);
    /* (2r+1)^3 grows fast enough that a large radius is refused rather than silently taking minutes. */
    assert_eval_eq("Head[MedianFilter[" FVOL ", 20]]", "MedianFilter", 0);
}

static void test_volume_median_is_a_rank_filter(void) {
    /* AN ORDER STATISTIC: every output value is one of the input values. A mean invents values; a
     * median must not, and this is the property that says so about the whole output at once. */
    assert_eval_eq("Module[{v = " FVOL "},"
                   " Complement[Union[Flatten[ImageData[MedianFilter[v, 2]]]],"
                   "            Union[Flatten[ImageData[v]]]] === {}]", "True", 0);
    /* THE DISCRIMINATING PROPERTY, and the reason a median filter exists. A lone impulse in an
     * otherwise constant volume is removed COMPLETELY -- the output is the constant and nothing else
     * -- because the impulse is outvoted in every window it appears in. A mean cannot do this at any
     * radius: it spreads the impulse over the whole neighbourhood instead. Asserted as the exact set
     * of output values, which is absolute. */
    assert_eval_eq("Module[{imp},"
                   " imp = Image3D[Table[If[z == 5 && y == 6 && x == 7, 1., 0.25],"
                   "   {z, 1, 10}, {y, 1, 12}, {x, 1, 14}]];"
                   " {Union[Flatten[ImageData[MedianFilter[imp, 1]]]],"
                   "  Length[Union[Flatten[ImageData[MeanFilter[imp, 1]]]]] > 1}]",
                   "{{0.25}, True}", 0);
}


/* Volumetric binarization. */
#define BVOL "Image3D[Table[N[Mod[z*13 + y*7 + x*3, 97]]/97, {z, 1, 10}, {y, 1, 12}, {x, 1, 14}]]"
#define BCUBE \
  "cl[q_, n_] := Max[1, Min[n, q]];" \
  "cube[y_, x_, z_, r_] := Flatten[Table[dv[[cl[z + a, 10], cl[y + b, 12], cl[x + c, 14]]]," \
  "  {a, -r, r}, {b, -r, r}, {c, -r, r}]];"

static void test_volume_binarize(void) {
    /* Typed Bit and PACKED, which check-image-packing also enforces -- asserted here too because the
     * type and the storage are separate claims and a nested Bit volume would satisfy neither cheaply. */
    assert_eval_eq("{ImageDimensions[Binarize[" BVOL "]], Part[Binarize[" BVOL "], 2],"
                   " Head[Part[Binarize[" BVOL "], 1]],"
                   " Union[Flatten[ImageData[Binarize[" BVOL "]]]]}",
                   "{{14, 12, 10}, \"Bit\", NDArray, {0.0, 1.0}}", 0);
    /* A stated threshold is exactly a comparison -- no Otsu, no scaling, nothing to approximate. */
    assert_eval_eq("Module[{v = " BVOL ", dv}, dv = ImageData[v];"
                   " ImageData[Binarize[v, 0.5]]"
                   " === N[Table[If[dv[[z, y, x]] > 0.5, 1, 0], {z, 1, 10}, {y, 1, 12}, {x, 1, 14}]]]",
                   "True", 0);
}

static void test_volume_local_adaptive_binarize(void) {
    assert_eval_eq("{ImageDimensions[LocalAdaptiveBinarize[" BVOL ", 2]],"
                   " Part[LocalAdaptiveBinarize[" BVOL ", 2], 2],"
                   " Head[Part[LocalAdaptiveBinarize[" BVOL ", 2], 1]]}",
                   "{{14, 12, 10}, \"Bit\", NDArray}", 0);
    /* Against the longhand definition, EXACTLY, whenever the threshold is moved off the window mean.
     * The window statistics come from three separable prefix passes rather than a 3-D summed-volume
     * table: the table's box sum is an eight-term inclusion-exclusion whose sign pattern is easy to
     * get wrong and produces plausible output when it is, and mean3_boxsum was already written and
     * already checked. */
    assert_eval_eq("Module[{v = " BVOL ", dv}, dv = ImageData[v];" BCUBE
                   " ImageData[LocalAdaptiveBinarize[v, 2, {1, -0.3, 0.02}]]"
                   " === N[Table[Module[{wi = cube[y, x, z, 2], mu, sd},"
                   "   mu = Mean[wi]; sd = Sqrt[Max[0., Mean[wi^2] - mu^2]];"
                   "   If[dv[[z, y, x]] > mu - 0.3 sd + 0.02, 1, 0]],"
                   "   {z, 1, 10}, {y, 1, 12}, {x, 1, 14}]]]", "True", 0);
    assert_eval_eq("Module[{v = " BVOL ", dv}, dv = ImageData[v];" BCUBE
                   " ImageData[LocalAdaptiveBinarize[v, 2, {0.9, 0, 0.}]]"
                   " === N[Table[If[dv[[z, y, x]] > 0.9 Mean[cube[y, x, z, 2]], 1, 0],"
                   "   {z, 1, 10}, {y, 1, 12}, {x, 1, 14}]]]", "True", 0);
    /* Mean-only is the boundary case at rank 3 exactly as at rank 2: with the threshold equal to the
     * window mean a voxel can tie it, and the prefix-sum mean differs from a direct sum in the last
     * bit. Two voxels of 1680 here, and the assertion is that EVERY disagreement is a tie. */
    assert_eval_eq("Module[{v = " BVOL ", dv, mine, ref, bad}, dv = ImageData[v];" BCUBE
                   " mine = ImageData[LocalAdaptiveBinarize[v, 2]];"
                   " ref = N[Table[If[dv[[z, y, x]] > Mean[cube[y, x, z, 2]], 1, 0],"
                   "   {z, 1, 10}, {y, 1, 12}, {x, 1, 14}]];"
                   " bad = Position[mine - ref, _?(# != 0 &)];"
                   " And @@ Table[Abs[dv[[p[[1]], p[[2]], p[[3]]]]"
                   "   - Mean[cube[p[[2]], p[[3]], p[[1]], 2]]] < 1.*^-15, {p, bad}]]", "True", 0);
    /* A uniform volume: the window mean IS the value, so nothing exceeds it. */
    assert_eval_eq("Union[Flatten[ImageData[LocalAdaptiveBinarize["
                   "Image3D[Table[0.5, {8}, {8}, {8}]], 2]]]]", "{0.0}", 0);
    /* Raising the offset can only turn voxels off. */
    assert_eval_eq("Module[{v = " BVOL ", n},"
                   " n = Table[Total[Flatten[ImageData[LocalAdaptiveBinarize[v, 2, {1, 0, c}]]]],"
                   "   {c, {0., 0.05, 0.2}}];"
                   " And @@ Table[n[[i]] >= n[[i + 1]], {i, 1, 2}]]", "True", 0);
}

static void test_volume_local_adaptive_beats_global(void) {
    /* THE REASON IT EXISTS, at rank 3. A checkerboard volume under a lighting ramp along x: no single
     * number separates the two tones at both ends, so a global threshold must collapse one end to a
     * single value while the local form keeps the pattern at both. Set membership, not a percentage. */
    assert_eval_eq("Module[{ramp, gl, lo},"
                   " ramp = Image3D[Table[N[(0.1 + 0.8 (x - 1)/31)"
                   "   * If[Mod[Quotient[z - 1, 4] + Quotient[y - 1, 4] + Quotient[x - 1, 4], 2] == 0,"
                   "        0.55, 1.]], {z, 1, 16}, {y, 1, 16}, {x, 1, 32}]];"
                   " gl = ImageData[Binarize[ramp]]; lo = ImageData[LocalAdaptiveBinarize[ramp, 3]];"
                   " {Length[Union[Flatten[gl[[All, All, 1 ;; 10]]]]],"
                   "  Length[Union[Flatten[lo[[All, All, 1 ;; 10]]]]],"
                   "  Length[Union[Flatten[lo[[All, All, -10 ;; -1]]]]]}]", "{1, 2, 2}", 0);
}


/* The volumetric distance transform. This one admits the strongest tests in the file: squared
 * distances on an integer lattice are exact integers, so Pythagorean offsets give exact distances and
 * the comparison is `===` rather than a tolerance. */
static void test_volume_distance_transform_is_exactly_euclidean(void) {
    /* One background voxel at (6,6,6); the distance at an offset is sqrt of the sum of squares, and
     * these offsets are chosen so that root is an integer. A separable transform that added distances
     * per axis instead of SQUARED distances would fail every one of these while still producing a
     * plausible-looking gradient. */
    assert_eval_eq("Module[{v, r},"
                   " v = Image3D[Table[If[z == 6 && y == 6 && x == 6, 0., 1.],"
                   "   {z, 1, 12}, {y, 1, 12}, {x, 1, 12}]];"
                   " r = ImageData[DistanceTransform[v]];"
                   " {r[[6, 6, 6]], r[[6, 6, 9]] === 3., r[[6, 9, 10]] === 5.,"
                   "  r[[7, 8, 8]] === 3., r[[8, 9, 12]] === 7., r[[7, 10, 10]] === Sqrt[33.]}]",
                   "{0.0, True, True, True, True, True}", 0);
    /* Against the DEFINITION -- the minimum Euclidean distance to any background voxel, computed by
     * brute force over every one of them -- and the agreement is exact. */
    assert_eval_eq("Module[{u, du, bg, ref},"
                   " u = Image3D[Table[If[Mod[z*5 + y*3 + x*7, 11] == 0, 0., 1.],"
                   "   {z, 1, 7}, {y, 1, 8}, {x, 1, 9}]];"
                   " du = ImageData[u]; bg = Position[du, 0.];"
                   " ref = Table[Min[Table[Sqrt[N[(z - b[[1]])^2 + (y - b[[2]])^2"
                   "   + (x - b[[3]])^2]], {b, bg}]], {z, 1, 7}, {y, 1, 8}, {x, 1, 9}];"
                   " ImageData[DistanceTransform[u]] === ref]", "True", 0);
    /* An all-background volume is everywhere zero, exactly. */
    assert_eval_eq("Union[Flatten[ImageData[DistanceTransform["
                   "Image3D[Table[0., {6}, {6}, {6}]]]]]]", "{0.0}", 0);
    /* The threshold argument selects what counts as background, and 0.5 on a 0/1 volume is the same
     * partition as the default. */
    assert_eval_eq("Module[{v},"
                   " v = Image3D[Table[If[z == 6 && y == 6 && x == 6, 0., 1.],"
                   "   {z, 1, 12}, {y, 1, 12}, {x, 1, 12}]];"
                   " ImageData[DistanceTransform[v, 0.5]] === ImageData[DistanceTransform[v]]]",
                   "True", 0);
}


/* Volumetric GradientFilter and DerivativeFilter. A LINEAR RAMP is the right test object: a central
 * difference reproduces a linear function's slope exactly, and the smoothing along the other axes
 * preserves linearity, so the answers are exact equalities rather than tolerances. The slope is 1/32 --
 * a power of two, so it is representable and `===` means what it says. */
#define RAMPX "Image3D[Table[N[(x - 1)/32], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}]]"

static void test_volume_gradient_of_a_ramp_is_exactly_the_slope(void) {
    assert_eval_eq("Module[{g = ImageData[GradientFilter[" RAMPX "]]},"
                   " {g[[4, 5, 6]] === 1./32,"
                   "  Union[Flatten[g[[2 ;; 7, 2 ;; 9, 2 ;; 11]]]] === {1./32}}]",
                   "{True, True}", 0);
    /* A constant volume has no gradient at all -- exactly zero, not merely small. */
    assert_eval_eq("Union[Flatten[ImageData[GradientFilter["
                   "Image3D[Table[0.5, {8}, {8}, {8}]]]]]]", "{0.0}", 0);
    /* AXIS COVARIANCE: the magnitude does not care which axis the ramp runs along. This is what
     * catches a transposed stencil or a mis-indexed component, neither of which spoils the look of
     * the output. */
    assert_eval_eq("Module[{rz, ry},"
                   " rz = Image3D[Table[N[(z - 1)/32], {z, 1, 12}, {y, 1, 10}, {x, 1, 8}]];"
                   " ry = Image3D[Table[N[(y - 1)/32], {z, 1, 8}, {y, 1, 12}, {x, 1, 10}]];"
                   " {Union[Flatten[ImageData[GradientFilter[rz]][[2 ;; 11, 2 ;; 9, 2 ;; 7]]]]"
                   "    === {1./32},"
                   "  Union[Flatten[ImageData[GradientFilter[ry]][[2 ;; 7, 2 ;; 11, 2 ;; 9]]]]"
                   "    === {1./32}}]", "{True, True}", 0);
}

static void test_volume_derivative_filter_signs_and_orders(void) {
    /* THE SIGN, asserted as an exact signed value. The stencils are pre-flipped because ImageConvolve
     * reflects its kernel, and the rank-2 code records that this was caught only this way -- a gradient
     * magnitude squares the sign away, so a flipped derivative looks perfect in every other test. */
    assert_eval_eq("ImageData[DerivativeFilter[" RAMPX ", {0, 0, 1}]][[4, 5, 6]] === 1./32",
                   "True", 0);
    /* The derivative along an axis the ramp does not vary in is zero. */
    assert_eval_eq("Abs[ImageData[DerivativeFilter[" RAMPX ", {1, 0, 0}]][[4, 5, 6]]] < 1.*^-17",
                   "True", 0);
    /* Order 0 on every axis is a pure smoothing, which is a legitimate request and not a derivative:
     * on a linear ramp a symmetric smoothing returns the centre value untouched. */
    assert_eval_eq("ImageData[DerivativeFilter[" RAMPX ", {0, 0, 0}]][[4, 5, 6]] === N[5/32]",
                   "True", 0);
    /* The SECOND derivative of a linear function is zero -- an exact property of the function, not of
     * the filter, which is what makes it a good check on the order-2 stencil. */
    assert_eval_eq("Abs[ImageData[DerivativeFilter[" RAMPX ", {0, 0, 2}]][[4, 5, 6]]] < 1.*^-17",
                   "True", 0);
    assert_eval_eq("ImageDimensions[DerivativeFilter[" RAMPX ", {0, 0, 1}]]", "{12, 10, 8}", 0);
    /* A rank-2 order spec on a volume, and an order the stencil table does not define, both decline. */
    assert_eval_eq("{Head[DerivativeFilter[" RAMPX ", {0, 1}]],"
                   " Head[DerivativeFilter[" RAMPX ", {0, 0, 5}]]}",
                   "{DerivativeFilter, DerivativeFilter}", 0);
}


/* ImageReflect for a volume. A reflection is a pure index permutation, so every property here is an
 * EXACT identity -- nothing is interpolated, so nothing can drift. */
#define RVOL "Image3D[Table[N[Mod[z*13 + y*7 + x*3, 97]]/97, {z, 1, 6}, {y, 1, 8}, {x, 1, 10}]]"

static void test_volume_reflect_names_the_right_axes(void) {
    assert_eval_eq("ImageDimensions[ImageReflect[" RVOL "]]", "{10, 8, 6}", 0);
    /* The default is the HEIGHT axis, matching the planar default, and each named pair selects one
     * axis: Top/Bottom height, Left/Right width, Front/Back depth -- Mathematica's own vocabulary for
     * volumes. Checked against Reverse on the corresponding level of the data, which is what
     * "reflect about this axis" means and leaves no room for a transposition to hide. */
    assert_eval_eq("Module[{v = " RVOL ", dv}, dv = ImageData[v];"
                   " {ImageData[ImageReflect[v]] === Reverse[dv, 2],"
                   "  ImageData[ImageReflect[v, Top]] === Reverse[dv, 2],"
                   "  ImageData[ImageReflect[v, Left]] === Reverse[dv, 3],"
                   "  ImageData[ImageReflect[v, Front]] === Reverse[dv, 1]}]",
                   "{True, True, True, True}", 0);
    /* Either name of a pair is the same operation: reflecting to the top and to the bottom are one
     * thing, which is already true of the planar version. */
    assert_eval_eq("ImageReflect[" RVOL ", Bottom] === ImageReflect[" RVOL ", Top]", "True", 0);
    /* Front and Back DECLINE on a plane, which has no depth axis. Reinterpreting them as some other
     * axis would turn a caller's mistake into a wrong picture. */
    assert_eval_eq("Head[ImageReflect[Image[Table[N[i*j]/100, {i, 1, 6}, {j, 1, 8}]], Front]]",
                   "ImageReflect", 0);
    assert_eval_eq("Head[ImageReflect[" RVOL ", Sideways]]", "ImageReflect", 0);
}

static void test_volume_reflect_algebra(void) {
    /* SELF-INVERSE on every axis, bit for bit. */
    assert_eval_eq("Module[{v = " RVOL "},"
                   " And @@ Table[ImageReflect[ImageReflect[v, sd], sd] === v,"
                   "   {sd, {Top, Left, Front}}]]", "True", 0);
    /* Reflections about DIFFERENT axes COMMUTE, exactly. This is the property that catches an axis
     * confused with another: a swapped pair still reflects something, and still round-trips, but stops
     * commuting with the axis it was confused for. */
    assert_eval_eq("Module[{v = " RVOL "},"
                   " {ImageReflect[ImageReflect[v, Left], Front]"
                   "    === ImageReflect[ImageReflect[v, Front], Left],"
                   "  ImageReflect[ImageReflect[v, Top], Front]"
                   "    === ImageReflect[ImageReflect[v, Front], Top]}]", "{True, True}", 0);
    /* All three composed is a point inversion, and applying that twice is the identity. */
    assert_eval_eq("Module[{v = " RVOL ", dv, inv}, dv = ImageData[v];"
                   " inv = ImageReflect[ImageReflect[ImageReflect[v, Left], Top], Front];"
                   " {ImageData[inv] === Reverse[Reverse[Reverse[dv, 1], 2], 3],"
                   "  ImageReflect[ImageReflect[ImageReflect[inv, Left], Top], Front] === v}]",
                   "{True, True}", 0);
    /* A colour volume: the channels must ride along rather than be permuted with the axes. */
    assert_eval_eq("Module[{cv},"
                   " cv = Image3D[Table[N[Mod[z*13 + y*7 + x*3 + ch*29, 97]]/97,"
                   "   {z, 1, 4}, {y, 1, 5}, {x, 1, 6}, {ch, 1, 3}]];"
                   " {ImageReflect[ImageReflect[cv, Front], Front] === cv,"
                   "  Dimensions[ImageData[ImageReflect[cv, Left]]]}]", "{True, {4, 5, 6, 3}}", 0);
}


/* ColorConvert for a volume, and the floating-point facts that make its tests non-obvious. */
#define CCVOL "Image3D[Table[N[Mod[z*13 + y*7 + x*3 + ch*29, 97]]/97," \
              " {z, 1, 5}, {y, 1, 6}, {x, 1, 7}, {ch, 1, 3}]]"
#define CGVOL "Image3D[Table[N[Mod[z*13 + y*7 + x*3, 97]]/97, {z, 1, 5}, {y, 1, 6}, {x, 1, 7}]]"

static void test_volume_colorconvert(void) {
    assert_eval_eq("Module[{g = ColorConvert[" CCVOL ", \"Grayscale\"]},"
                   " {ImageDimensions[g], ImageChannels[g]}]", "{{7, 6, 5}, 1}", 0);
    /* THE ONE EXACT IDENTITY: an already-grey volume is copied, not weighted, so it comes back bit for
     * bit. */
    assert_eval_eq("ImageData[ColorConvert[" CGVOL ", \"Grayscale\"]] === ImageData[" CGVOL "]",
                   "True", 0);
    assert_eval_eq("ImageData[ColorConvert[" CCVOL ", \"Gray\"]]"
                   " === ImageData[ColorConvert[" CCVOL ", \"Grayscale\"]]", "True", 0);
    /* Against the Rec. 601 formula written out, to 1e-15 and NOT exactly. The difference is purely
     * summation order: Plus is Orderless, so a Mathilda-level reference adds the three weighted
     * channels in the canonical order while the implementation adds them in the order written. 27 of
     * 210 voxels differ here, every one by 1.11e-16. Asserting `===` would fail for a reason that has
     * nothing to do with the conversion. */
    assert_eval_eq("Module[{cv = " CCVOL ", d, ref},"
                   " d = ImageData[cv];"
                   " ref = Table[0.299 d[[z, y, x, 1]] + 0.587 d[[z, y, x, 2]]"
                   "   + 0.114 d[[z, y, x, 3]], {z, 1, 5}, {y, 1, 6}, {x, 1, 7}];"
                   " Max[Abs[Flatten[ImageData[ColorConvert[cv, \"Grayscale\"]] - ref]]] < 1.*^-15]",
                   "True", 0);
    assert_eval_eq("Head[ColorConvert[" CCVOL ", \"CMYK\"]]", "ColorConvert", 0);
}

static void test_colorconvert_equal_channels_is_value_dependent(void) {
    /* Three EQUAL channels do not reliably convert to that value, and this pins the ACTUAL behaviour
     * rather than a hope about it: the Rec. 601 weights sum to 0.9999999999999999 in the order they are
     * applied, though to exactly 1.0 in any order starting with 0.114, so the final rounding lands on
     * the input for some values and one ulp below for others.
     *
     * 0.75 is exact and 0.7 is not. Asserting BOTH is what documents the behaviour in the suite: if
     * someone "corrects" the weights to sum to 1.0 in double, this fails and points at why. */
    assert_eval_eq("Module[{f, a, b},"
                   " f[v_] := First[Union[Flatten[ImageData[ColorConvert["
                   "   Image3D[Table[v, {z, 1, 3}, {y, 1, 3}, {x, 1, 3}, {ch, 1, 3}]],"
                   "   \"Grayscale\"]]]]];"
                   " a = f[0.75]; b = f[0.7];"
                   " {a === 0.75, b === 0.7, Abs[b - 0.7] < 1.*^-15}]",
                   "{True, False, True}", 0);
    /* And the printed form cannot be trusted to reveal this: Mathilda prints that weight sum as 1.0
     * even through InputForm, so a test written from the printed value would assert the wrong thing
     * and pass for the wrong reason. The difference only shows by subtraction. */
    assert_eval_eq("1.0 - (0.299 + 0.587 + 0.114) != 0", "True", 0);
}


/* --------------------------------------------------------------- Import / Export
 *
 * The round trip is the only property worth asserting here, and it has to be asserted in the two
 * directions separately because they fail differently: a wrong Export writes a file no other program
 * can read (invisible to a test that reads it back with our own Import), and a wrong Import misreads
 * a correct file. Byte quantisation bounds the error at 1/510 per sample -- half a level -- so the
 * comparison is a bound, not an equality, and stating the bound is what makes the test able to catch
 * a channel swap: swapping red and blue stays well inside any loose tolerance but not inside this
 * one, for an image whose channels differ. */
#define IOIMG "Image[Table[{N[i/12], N[j/16], N[Mod[i + j, 8]]/8}, {i, 1, 12}, {j, 1, 16}], \"Real\"]"

static void test_png_round_trip_is_exact_to_the_quantisation(void) {
    assert_eval_eq("Module[{f = \"/tmp/mathilda_test_rt.png\", a = " IOIMG ", b},"
                   " b = Import[Export[f, a]];"
                   " {ImageDimensions[b] === ImageDimensions[a], ImageChannels[b],"
                   "  Max[Abs[Flatten[ImageData[b] - ImageData[a]]]] <= 1/510. + 1.*^-12}]",
                   "{True, 3, True}", 0);
    /* A grey file stays grey and an RGBA file keeps its alpha. Forcing three channels would invent
     * two for the first and silently discard transparency from the second, and both would still
     * round-trip a dimension check. */
    assert_eval_eq("Module[{g = Image[Table[N[i/8], {i, 1, 8}, {j, 1, 8}], \"Real\"],"
                   "        a = Image[Table[{0.2, 0.4, 0.6, 0.8}, {i, 1, 8}, {j, 1, 8}], \"Real\"]},"
                   " {ImageChannels[Import[Export[\"/tmp/mathilda_test_g.png\", g]]],"
                   "  ImageChannels[Import[Export[\"/tmp/mathilda_test_a.png\", a]]]}]",
                   "{1, 4}", 0);
    /* An imported image is the SAME representation a filter produces -- packed, canonical, "Real" --
     * so it is not a second-class citizen that every downstream head has to special-case. */
    assert_eval_eq("Module[{b = Import[Export[\"/tmp/mathilda_test_rt.png\", " IOIMG "]]},"
                   " {Head[Part[b, 1]], ImageType[b]}]", "{NDArray, \"Real\"}", 0);
}

static void test_export_clamps_rather_than_wraps(void) {
    /* An unsharp mask legitimately overshoots the unit interval and 8-bit output has nowhere to put
     * that. Wrapping would turn a bright highlight black, which reads as a bug in the filter rather
     * than in the writer. */
    assert_eval_eq("ImageData[Import[Export[\"/tmp/mathilda_test_clamp.png\","
                   " Image[{{2.0, -1.0}, {1.0, 0.0}}, \"Real\"]]]]",
                   "{{1.0, 0.0}, {1.0, 0.0}}", 0);
}

static void test_import_distinguishes_missing_from_unclaimed(void) {
    /* Two different failures, deliberately given two different answers: a file we can read and
     * could not is $Failed, and a format nothing here claims stays unevaluated so a later format
     * handler can pick it up. Collapsing both to $Failed would make Import[x] look implemented for
     * every format in existence. */
    assert_eval_eq("Import[\"/tmp/mathilda_no_such_file_2718.png\"]", "$Failed", 0);
    assert_eval_eq("Head[Import[\"/tmp/mathilda_whatever.xyz\"]]", "Import", 0);
    /* A volume has no single raster, and quietly writing its middle slice would be a lie about what
     * was exported. */
    assert_eval_eq("Head[Export[\"/tmp/mathilda_vol.png\","
                   " Image3D[Table[0.5, {z, 1, 2}, {y, 1, 2}, {x, 1, 2}], \"Real\"]]]",
                   "Export", 0);
}

static void test_lossy_and_lossless_formats_differ_as_expected(void) {
    /* JPEG is lossy at quality 90: bounded, not exact. Asserting a bound in BOTH directions is what
     * makes this a test rather than a tautology -- an upper bound alone would pass if Export silently
     * wrote a PNG, and a lower bound alone would pass on garbage. */
    assert_eval_eq("Module[{a = " IOIMG ", e},"
                   " e = Mean[Flatten[Abs[ImageData[Import[Export[\"/tmp/mathilda_test.jpg\", a]]]"
                   "                      - ImageData[a]]]];"
                   " {e > 0, e < 0.1}]", "{True, True}", 0);
    /* BMP and TGA are lossless too, so they carry the same bound as PNG. */
    assert_eval_eq("Module[{a = " IOIMG "},"
                   " Max[Table[Max[Abs[Flatten[ImageData[Import[Export[f, a]]] - ImageData[a]]]],"
                   "  {f, {\"/tmp/mathilda_t.bmp\", \"/tmp/mathilda_t.tga\"}}]] <= 1/510. + 1.*^-12]",
                   "True", 0);
    /* The format may also be stated instead of inferred, which is the only way to write a file whose
     * name does not end in one. */
    assert_eval_eq("ImageDimensions[Import[Export[\"/tmp/mathilda_noext\", " IOIMG ", \"PNG\"],"
                   " \"Image\"]]", "{16, 12}", 0);
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
    TEST(test_vanherk_agrees_with_an_independent_reference);
    TEST(test_correlation_is_convolution_reflected);
    TEST(test_ncc_scores_exactly_one_at_an_exact_match);
    TEST(test_right_angle_rotation_is_exact);
    TEST(test_reflection_is_exact_and_self_inverse);
    TEST(test_arbitrary_angle_rotation_round_trips_smooth_content);
    TEST(test_connectivity_is_the_discriminating_property);
    TEST(test_a_u_shape_needs_the_union_find);
    TEST(test_labels_are_contiguous_in_raster_order);
    TEST(test_components_interact_with_morphology);
    TEST(test_median_removes_an_outlier_exactly);
    TEST(test_median_is_not_separable);
    TEST(test_meanfilter_is_a_box_convolution);
    TEST(test_distance_transform_is_exactly_euclidean);
    TEST(test_distance_transform_degenerate_cases);
    TEST(test_distance_transform_relates_to_morphology);
    TEST(test_volume_downsampling_does_not_alias);
    TEST(test_volume_resize_respects_the_axis_order);
    TEST(test_levels_counts_sum_to_the_pixel_count);
    TEST(test_imageadjust_stretches_exactly_and_is_idempotent);
    TEST(test_imageadjust_parametric_curve);
    TEST(test_pad_then_crop_is_the_identity);
    TEST(test_pad_side_convention_and_modes);
    TEST(test_crop_trims_a_uniform_border);
    TEST(test_volume_pad_then_crop_is_the_identity);
    TEST(test_volume_pad_axis_conventions);
    TEST(test_volume_pad_fill_modes);
    TEST(test_ncc_finds_the_patch_it_was_given);
    TEST(test_ncc_agrees_with_its_definition);
    TEST(test_ncc_degenerate_and_invariant);
    TEST(test_fft_convolution_agrees_with_the_definition);
    TEST(test_fft_convolution_across_channels);
    TEST(test_volume_fft_convolution_agrees_with_the_definition);
    TEST(test_volume_fft_convolution_across_channels);
    TEST(test_corner_response_is_zero_where_there_is_no_corner);
    TEST(test_corner_response_finds_corners_and_rotates_with_them);
    TEST(test_corner_options_are_wired);
    TEST(test_corner_separation_is_respected);
    TEST(test_corner_feature_limit_keeps_the_strongest);
    TEST(test_volume_corner_eigenvalue_hierarchy);
    TEST(test_volume_corner_axis_covariance);
    TEST(test_volume_imagecorners);
    TEST(test_corner_named_options);
    TEST(test_local_adaptive_binarize_matches_the_definition);
    TEST(test_local_adaptive_binarize_properties);
    TEST(test_local_adaptive_beats_global_under_a_ramp);
    TEST(test_volume_morphology_algebra);
    TEST(test_volume_morphology_matches_the_definition);
    TEST(test_volume_mean_and_median_match_the_definition);
    TEST(test_volume_median_is_a_rank_filter);
    TEST(test_volume_binarize);
    TEST(test_volume_local_adaptive_binarize);
    TEST(test_volume_local_adaptive_beats_global);
    TEST(test_volume_distance_transform_is_exactly_euclidean);
    TEST(test_volume_gradient_of_a_ramp_is_exactly_the_slope);
    TEST(test_volume_derivative_filter_signs_and_orders);
    TEST(test_volume_reflect_names_the_right_axes);
    TEST(test_volume_reflect_algebra);
    TEST(test_volume_colorconvert);
    TEST(test_colorconvert_equal_channels_is_value_dependent);
    TEST(test_png_round_trip_is_exact_to_the_quantisation);
    TEST(test_export_clamps_rather_than_wraps);
    TEST(test_import_distinguishes_missing_from_unclaimed);
    TEST(test_lossy_and_lossless_formats_differ_as_expected);

    printf("All image tests passed.\n");
    return 0;
}
