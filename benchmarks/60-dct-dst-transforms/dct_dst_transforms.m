(* ==========================================================================
   Experiment 60 -- Discrete cosine/sine & 2-D transforms (group E)
   ==========================================================================
   WHAT IT MEASURES.  The real trigonometric transforms and the 2-D DFT in
   src/fourier.c: FourierDCT (types 1, 2, 4), FourierDST (types 1, 2), and a
   512x512 Fourier.  Experiment 17 covers the 1-D complex FFT and convolutions;
   this covers the cosine/sine transforms (JPEG/PDE workhorses) and the 2-D case.

   Baseline is scipy.fft.dct/dst and numpy.fft.fft2 (pocketfft, compiled), so
   this is a compiled-vs-compiled execution comparison at O(n log n).

   NORMALIZATION is convention-dependent and handled entirely in the check (as
   in experiment 17): Mathematica's FourierDCT is orthonormal-with-endpoint
   scaling and its Fourier divides by sqrt(n).  The Python side rescales scipy's
   result to match -- see the docstring in the .py file for the exact factor per
   type.  Only the types with a clean scipy equivalent are included (DCT 1/2/4,
   DST 1/2); DCT-3/DST-3 use a scaling with no single scipy setting and are left
   out rather than reconciled fragilely.

   TIMING runs on a 10^6 vector / 512x512 matrix; CHECKS use the deterministic
   {1..8} / {{1,2},{3,4}}, so the three systems pin the same number.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"FourierDCT", "FourierDST", "Fourier"}];

seed[];
n = 1000000;
v = rand01[{n}];
c8 = N[Range[8]];                       (* deterministic check input *)
mBig = rand01[{512, 512}];
m2 = N[{{1, 2}, {3, 4}}];               (* deterministic check input *)

bench["FourierDCT type 1 of 10^6", FourierDCT[v, 1];];
check["FourierDCT type 1 of 10^6", Round[10^6 Total[FourierDCT[c8, 1]]]];

bench["FourierDCT type 2 of 10^6", FourierDCT[v, 2];];
check["FourierDCT type 2 of 10^6", Round[10^6 Total[FourierDCT[c8, 2]]]];

bench["FourierDCT type 4 of 10^6", FourierDCT[v, 4];];
check["FourierDCT type 4 of 10^6", Round[10^6 Total[FourierDCT[c8, 4]]]];

bench["FourierDST type 1 of 10^6", FourierDST[v, 1];];
check["FourierDST type 1 of 10^6", Round[10^6 Total[FourierDST[c8, 1]]]];

bench["FourierDST type 2 of 10^6", FourierDST[v, 2];];
check["FourierDST type 2 of 10^6", Round[10^6 Total[FourierDST[c8, 2]]]];

bench["Fourier 2D 512x512", Fourier[mBig];];
check["Fourier 2D 512x512", Round[10^6 Total[Abs[Flatten[Fourier[m2]]]]]];
