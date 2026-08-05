(* Experiment 17 -- FFT, convolution, correlation.
   MEASURES src/fourier.c and src/convolutions.c.

   READ THE BUILD WARNING FIRST.  Without FFTW the makefile substitutes "a naive
   O(n^2) fallback" and every row here measures that instead: on this tree a
   32768-point transform was 1050 ms without FFTW and 1.0 ms with it -- a factor
   of a thousand that has nothing to do with src/fourier.c's own quality.

   Power-of-two AND non-power-of-two sizes are both here, because a mixed-radix
   size is exactly where a weak FFT falls back to something quadratic. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Fourier", "InverseFourier", "ListConvolve", "ListCorrelate",
         "FourierDCT", "PeriodogramArray"}];

(* Wolfram's Fourier normalises by 1/Sqrt[n]; the checks account for that. *)
small = {1., 2., 3., 4., 5., 6., 7., 8.};

d1 = rand01[{262144}];              (* 2^18: the friendly case *)
bench["Fourier 2^18 (262144)", Fourier[d1];];
check["Fourier 2^18 (262144)", Round[10^6 Re[First[Fourier[small]]]]];

d2 = rand01[{1200000}];             (* 2^7 * 3 * 5^5: mixed radix *)
bench["Fourier 1200000 (mixed radix)", Fourier[d2];];
check["Fourier 1200000 (mixed radix)", Round[10^6 Re[Last[Fourier[small]]]]];

d3 = rand01[{262143}];              (* 3^3 * 7 * ... : an awkward size *)
bench["Fourier 262143 (awkward size)", Fourier[d3];];
check["Fourier 262143 (awkward size)",
  Round[10^6 Abs[Im[Part[Fourier[small], 2]]]]];

bench["InverseFourier 2^18", InverseFourier[d1];];
check["InverseFourier 2^18",
  Round[10^6 Re[First[InverseFourier[Fourier[small]]]]]];

(* Convolution of two length-100000 sequences: the FFT path if there is one. *)
a = rand01[{100000}]; b = rand01[{2048}];
bench["ListConvolve 100000 x 2048", ListConvolve[b, a];];
check["ListConvolve 100000 x 2048",
  Round[10^6 First[ListConvolve[{1., 2.}, {3., 4., 5.}]]]];

bench["ListCorrelate 100000 x 2048", ListCorrelate[b, a];];
check["ListCorrelate 100000 x 2048",
  Round[10^6 First[ListCorrelate[{1., 2.}, {3., 4., 5.}]]]];
