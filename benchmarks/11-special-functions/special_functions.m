(* ==========================================================================
   Experiment 11 -- Special functions over machine arrays
   ==========================================================================
   WHAT IT MEASURES.  src/special_functions/ (41 files) reached through the
   packed-array path -- src/ndkernels.c and sf_machine.c.  Unlike group A this
   is an EXECUTION comparison: scipy.special is compiled C, so a gap here is
   overhead or a missing vector kernel, not a missing algorithm.

   Sizes are 10^6, comfortably above PACK_MIN_ELEMENTS (src/pack.h:84), so the
   buffer path is the one being measured.

   Checks are a rounded scalar from a small DETERMINISTIC input -- never a sum
   over the random timing data, which the three systems cannot align.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"Gamma", "LogGamma", "BesselJ", "BesselY", "Erf", "Erfc",
         "PolyGamma", "AiryAi", "Zeta", "ProductLog", "EllipticK", "Beta"}];

n = 1000000;
v = rand01[{n}] + 1/2;              (* in [0.5, 1.5]: safe for every function *)

bench["Gamma over 10^6", Gamma[v];];
check["Gamma over 10^6", Round[10^6 N[Gamma[3/2]]]];

bench["Erf over 10^6", Erf[v];];
check["Erf over 10^6", Round[10^6 N[Erf[1/2]]]];

bench["BesselJ[0, .] over 10^6", BesselJ[0, v];];
check["BesselJ[0, .] over 10^6", Round[10^6 N[BesselJ[0, 1/2]]]];

bench["LogGamma over 10^6", LogGamma[v];];
check["LogGamma over 10^6", Round[10^6 N[LogGamma[3/2]]]];

bench["PolyGamma[0, .] over 10^6", PolyGamma[0, v];];
check["PolyGamma[0, .] over 10^6", Round[10^6 N[PolyGamma[0, 3/2]]]];

benchIf["AiryAi over 10^6", "AiryAi", AiryAi[v];];
checkIf["AiryAi over 10^6", "AiryAi", Round[10^6 N[AiryAi[1/2]]]];

benchIf["Zeta over 10^6", "Zeta", Zeta[v + 1];];
checkIf["Zeta over 10^6", "Zeta", Round[10^6 N[Zeta[5/2]]]];
