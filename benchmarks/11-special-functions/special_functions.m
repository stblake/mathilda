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

(* The Legendre elliptic integrals need their own inputs: v reaches past 1,
   and m > 1 is where the complete integrals go complex and the machine kernels
   decline -- timing that would measure the DECLINE path, not the kernel.
   m in [0, 0.95], phi in [0, 1.5] radians keeps 1 - m Sin[phi]^2 >= 0.05, so
   the incomplete kernels stay inside the real principal domain too. *)
em = 95 rand01[{n}]/100;
ep = 3 rand01[{n}]/2;

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

(* --- Legendre elliptic integrals ---------------------------------------- *)
(* SciPy's parameter convention for ellipk/ellipe/ellipkinc/ellipeinc is m = k^2,
   the same as Mathilda's and Mathematica's, so these are direct counterparts.
   SciPy has no third kind; its column composes Carlson R_F and R_J, which is
   what a NumPy user would write (identity checked against mpmath at <= 1 ulp). *)
benchIf["EllipticK over 10^6", "EllipticK", EllipticK[em];];
checkIf["EllipticK over 10^6", "EllipticK", Round[10^6 N[EllipticK[1/2]]]];

benchIf["EllipticE over 10^6", "EllipticE", EllipticE[em];];
checkIf["EllipticE over 10^6", "EllipticE", Round[10^6 N[EllipticE[1/2]]]];

benchIf["EllipticF[., m] over 10^6", "EllipticF", EllipticF[ep, em];];
checkIf["EllipticF[., m] over 10^6", "EllipticF", Round[10^6 N[EllipticF[1, 1/2]]]];

benchIf["EllipticE[., m] over 10^6", "EllipticE", EllipticE[ep, em];];
checkIf["EllipticE[., m] over 10^6", "EllipticE", Round[10^6 N[EllipticE[1, 1/2]]]];

benchIf["EllipticPi[n, m] over 10^6", "EllipticPi", EllipticPi[1/2, em];];
checkIf["EllipticPi[n, m] over 10^6", "EllipticPi", Round[10^6 N[EllipticPi[1/3, 1/2]]]];

benchIf["EllipticPi[n, ., m] over 10^6", "EllipticPi", EllipticPi[1/2, ep, em];];
checkIf["EllipticPi[n, ., m] over 10^6", "EllipticPi",
        Round[10^6 N[EllipticPi[1/3, 1, 1/2]]]];
