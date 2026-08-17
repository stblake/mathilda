(* ==========================================================================
   Experiment 55 -- Special functions over machine arrays (group E)
   ==========================================================================
   WHAT IT MEASURES.  The elementwise special-function surface reached through a
   10^6 packed buffer -- src/ndkernels.c into src/special_functions/ -- the
   direct analogue of a scipy.special ufunc.  Experiment 11 covers
   Gamma/Erf/BesselJ/LogGamma/PolyGamma/AiryAi/Zeta; this covers the rest of the
   library with a scipy.special counterpart, and in doing so answers a sharper
   question: WHICH heads have a SIMD vector kernel and which fall back to
   scalar-by-scalar threading.

   THE FINDING.  The Fresnel integrals, the complementary/imaginary error
   functions, the exponential/trigonometric integrals and the two-argument Beta
   are vectorized -- tens of ms over 10^6, a fair overhead comparison against
   scipy's compiled ufuncs.  The modified Bessel functions BesselI/BesselK now
   also run their registered NDArray kernel (tens of ms over 10^6, faster than
   scipy.special.iv/kv): the transparency gate used to force any rule-bearing
   head -- and every Bessel head carries half-integer DownValues -- off the
   buffer path, so they threaded one element at a time (~20-40 s).  The gate now
   also admits a head with a matching-arity NDArray kernel (eval.c), with the
   .m rules guarded by !NDArrayQ so they decline on a buffer.  benchOnce is kept
   so the row stays comparable across the history.

   Checks are a rounded scalar at a fixed rational input, never a sum over the
   random timing data.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"FresnelC", "FresnelS", "Erfc", "Erfi", "ExpIntegralEi",
         "SinIntegral", "CosIntegral", "Beta", "BesselI", "BesselK"}];

seed[];
n = 1000000;
v = rand01[{n}] + 1/2;              (* in [0.5, 1.5]: safe for every function *)
w = v + 1/2;

(* ---- vectorized: buffer-speed, fair overhead comparison ---------------- *)

bench["FresnelC over 10^6", FresnelC[v];];
check["FresnelC over 10^6", Round[10^6 N[FresnelC[1/2]]]];

bench["FresnelS over 10^6", FresnelS[v];];
check["FresnelS over 10^6", Round[10^6 N[FresnelS[1/2]]]];

bench["Erfc over 10^6", Erfc[v];];
check["Erfc over 10^6", Round[10^6 N[Erfc[1/2]]]];

bench["Erfi over 10^6", Erfi[v];];
check["Erfi over 10^6", Round[10^6 N[Erfi[1/2]]]];

bench["ExpIntegralEi over 10^6", ExpIntegralEi[v];];
check["ExpIntegralEi over 10^6", Round[10^6 N[ExpIntegralEi[3/2]]]];

bench["SinIntegral over 10^6", SinIntegral[v];];
check["SinIntegral over 10^6", Round[10^6 N[SinIntegral[3/2]]]];

bench["CosIntegral over 10^6", CosIntegral[v];];
check["CosIntegral over 10^6", Round[10^6 N[CosIntegral[3/2]]]];

bench["Beta[.,.] over 10^6", Beta[v, w];];
check["Beta[.,.] over 10^6", Round[10^6 N[Beta[3/2, 2]]]];

(* ---- no vector kernel: scalar threading, benchOnce --------------------- *)

benchOnce["BesselI[0, .] over 10^6", BesselI[0, v];];
check["BesselI[0, .] over 10^6", Round[10^6 N[BesselI[0, 3/2]]]];

benchOnce["BesselK[0, .] over 10^6", BesselK[0, v];];
check["BesselK[0, .] over 10^6", Round[10^6 N[BesselK[0, 3/2]]]];
