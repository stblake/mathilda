(* 64-image-resize -- ImageResize against skimage.

   Downsampling uses AREA AVERAGING on both sides (Mathilda's Automatic, skimage's
   downscale_local_mean), and upsampling uses bilinear with anti-aliasing off, so the two are
   doing the same arithmetic rather than approximately similar arithmetic.

   As in 63-image-convolve the image is built outside the timed region: constructing 262144 Expr
   reals is not resampling, and timing it would say so. *)

n = 512;
img = Image[Table[N[Mod[x*7 + y*13, 251]/251], {y, n}, {x, n}]];

bench["IR area 512->256", ImageResize[img, {256, 256}];];
bench["IR area 512->128", ImageResize[img, {128, 128}];];
bench["IR nearest 512->256", ImageResize[img, {256, 256}, Resampling -> "Nearest"];];
bench["IR bilinear 512->1024", ImageResize[img, {1024, 1024}];];

(* Area averaging PRESERVES THE MEAN exactly at an integer factor, because every source pixel
   contributes with equal total weight. A conservation law, so it is a check rather than a
   timing -- and it is what fails if the coverage weights are ever normalised wrongly. *)
check["IR area preserves mean",
  Round[10^9 (Mean[Flatten[ImageData[ImageResize[img, {256, 256}]]]]
              - Mean[Flatten[ImageData[img]]])]];

(* THE aliasing discriminator. A 4x4 checkerboard to 2x2: area averaging gives 0.5 everywhere,
   nearest gives a flat field with the pattern annihilated. *)
check["IR checkerboard area total",
  Round[10^6 Total[Flatten[ImageData[ImageResize[
    Image[{{0.,1.,0.,1.},{1.,0.,1.,0.},{0.,1.,0.,1.},{1.,0.,1.,0.}}], {2, 2}]]]]];

check["IR block mean exact",
  Round[10^9 Total[Flatten[ImageData[ImageResize[
    Image[{{0.,0.25,0.5,0.75},{0.,0.25,0.5,0.75},{1.,1.,1.,1.},{1.,1.,1.,1.}}], {2, 2}]]]]];
