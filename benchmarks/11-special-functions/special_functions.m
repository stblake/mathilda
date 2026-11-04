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

(* --- Mathieu functions ---------------------------------------------------- *)
(* MathieuC[a, q, z] at a generic (non-characteristic) a: SciPy has no counterpart
   (only the integer orders), so that row is Mathilda vs Mathematica. The plan
   (characteristic exponent + Floquet coefficients) is built once per (a, q),
   then the buffer is a short Clenshaw-like sum per element. a = 2, q = 1 is in a
   stability band, where every normalization convention coincides, so the check
   value is comparable across systems. *)
benchIf["MathieuC[a, q, .] over 10^6", "MathieuC", MathieuC[2, 1, v];];
checkIf["MathieuC[a, q, .] over 10^6", "MathieuC", Round[10^6 N[MathieuC[2, 1, 1/2]]]];

benchIf["MathieuCPrime[a, q, .] over 10^6", "MathieuCPrime", MathieuCPrime[2, 1, v];];
checkIf["MathieuCPrime[a, q, .] over 10^6", "MathieuCPrime",
        Round[10^6 N[MathieuCPrime[2, 1, 1/2]]]];

(* The characteristic order ce_2(z, 5): a = a_2(5) to 16 digits, which all three
   systems treat as the characteristic value (DLMF normalization). SciPy's
   mathieu_cem(2, 5, x) takes x in DEGREES. *)
benchIf["MathieuC ce_2(., 5) over 10^6", "MathieuC", MathieuC[7.449109739529178, 5., v];];
checkIf["MathieuC ce_2(., 5) over 10^6", "MathieuC",
        Round[10^6 MathieuC[7.449109739529178, 5., 0.5]]];

(* Arbitrary precision: one value at 10^4 digits (no SciPy counterpart). *)
benchIf["MathieuC at 10^4 digits", "MathieuC", MathieuC[2, 1, 3.7`10000];];
checkIf["MathieuC at 10^4 digits", "MathieuC", Round[10^6 MathieuC[2, 1, 3.7`30]]];

(* MathieuS[a, q, z], the odd Mathieu function: the same (a, q) plan as MathieuC
   (one cached set of Floquet coefficients; MathieuC + I MathieuS is the Floquet
   solution), a sine sum per element. Generic a again has no SciPy counterpart. *)
benchIf["MathieuS[a, q, .] over 10^6", "MathieuS", MathieuS[2, 1, v];];
checkIf["MathieuS[a, q, .] over 10^6", "MathieuS", Round[10^6 N[MathieuS[2, 1, 1/2]]]];

benchIf["MathieuSPrime[a, q, .] over 10^6", "MathieuSPrime", MathieuSPrime[2, 1, v];];
checkIf["MathieuSPrime[a, q, .] over 10^6", "MathieuSPrime",
        Round[10^6 N[MathieuSPrime[2, 1, 1/2]]]];

(* The characteristic order se_2(z, 5): a = b_2(5) to 16 digits (DLMF
   normalization in all three systems). SciPy's mathieu_sem(2, 5, x) takes x in
   DEGREES. *)
benchIf["MathieuS se_2(., 5) over 10^6", "MathieuS", MathieuS[2.0994604454866654, 5., v];];
checkIf["MathieuS se_2(., 5) over 10^6", "MathieuS",
        Round[10^6 MathieuS[2.0994604454866654, 5., 0.5]]];

(* Arbitrary precision: one value at 10^4 digits (no SciPy counterpart). *)
benchIf["MathieuS at 10^4 digits", "MathieuS", MathieuS[2, 1, 3.7`10000];];
checkIf["MathieuS at 10^4 digits", "MathieuS", Round[10^6 MathieuS[2, 1, 3.7`30]]];
