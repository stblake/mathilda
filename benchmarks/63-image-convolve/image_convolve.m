(* 63-image-convolve -- ImageConvolve against scipy.ndimage.convolve.

   Both sides use the SAME kernel (a normalised Gaussian, sigma = r/2) and the SAME boundary
   rule (replicate the edge pixel: Mathilda's "Fixed", scipy's mode='nearest'), so the check
   values are directly comparable rather than approximately so.

   The image is built OUTSIDE the timed region on purpose. An early measurement that included
   Image[Table[...]] read as 72x slower than scipy; almost all of it was constructing 262144
   Expr reals, not convolving. Timing the construction here would measure expression building
   and call it convolution. *)

n = 512;
data = Table[N[Mod[x*7 + y*13, 251]/251], {y, n}, {x, n}];
img = Image[data];
small = Image[Table[N[Mod[x*7 + y*13, 251]/251], {y, 64}, {x, 64}]];

bench["IC gaussian r=1 512x512", ImageConvolve[img, GaussianMatrix[1]];];
bench["IC gaussian r=2 512x512", ImageConvolve[img, GaussianMatrix[2]];];
bench["IC gaussian r=4 512x512", ImageConvolve[img, GaussianMatrix[4]];];

(* r = 0 is a single tap, so it measures the marshalling floor -- loading the pixels into a
   buffer and rebuilding the result as expressions -- with essentially no arithmetic. Worth a
   row of its own, because at small radii that floor DOMINATES and it is the thing to optimise
   rather than the inner loop. *)
bench["IC marshalling floor r=0 512x512", ImageConvolve[img, GaussianMatrix[0]];];

bench["GaussianFilter r=2 512x512", GaussianFilter[img, 2];];

(* Checks on the 64x64 image, scaled to integers so both languages agree exactly on the text. *)
check["IC total r=2 64x64",
  Round[10^6 Total[Flatten[ImageData[ImageConvolve[small, GaussianMatrix[2]]]]]]];
check["IC corner r=2 64x64",
  Round[10^9 Part[ImageData[ImageConvolve[small, GaussianMatrix[2]]], 1, 1]]];
check["IC centre r=2 64x64",
  Round[10^9 Part[ImageData[ImageConvolve[small, GaussianMatrix[2]]], 32, 32]]];

(* A constant image through a kernel summing to 1 is unchanged EVERYWHERE including the border.
   This is what fails if the padding is ever changed to zero-fill, which darkens the edges. *)
check["IC constant preserved",
  Round[10^9 Max[Abs[Flatten[ImageData[ImageConvolve[Image[Table[0.25, {8}, {8}]],
                                                    GaussianMatrix[2]]] - 0.25]]]]];
