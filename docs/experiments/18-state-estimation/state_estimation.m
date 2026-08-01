(* ==========================================================================
   Experiment 18 -- State estimation: small dense matrices in a long loop
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.

       Mathilda      ./Mathilda -file state_estimation.m
       Mathematica   wolframscript -file state_estimation.m
       Python        python3 state_estimation.py

   WHAT IT MEASURES.  Every other benchmark in the suite is dominated by the
   cost of ARRAYS.  This one is dominated by the cost of a CALL.

   A 6-state constant-acceleration Kalman filter over 20000 measurements does,
   per step, two 6x6 products, a 6x6-by-6 matrix-vector, a 2x6-by-6x6, a 2x2
   inverse and half a dozen small adds.  The flop count for the whole run is
   about 2e7 -- a single dgemm in the linear-algebra suite does more than that
   in one call.  What is being measured is 20000 iterations x ~10 operations
   of per-operation dispatch on data too small for any kernel to matter.

   NumPy is not fast here either, and that is the point: the row shows where
   all three systems live when the data is small.  The ensemble filter beside
   it is the same estimator written the array way -- 4096 members at once --
   so the two regimes can be read together.

   DETERMINISM.  The measurement sequence is a fixed pair of sinusoids, so
   Tr[P] at the end is an exact cross-system check.  The filter is stable, so
   P converges to a steady state rather than drifting, which makes the check
   insensitive to the last-bit differences three systems will always have.
   ========================================================================== *)

(* ---- shared reporting helpers (identical in every experiment file) ------- *)

SetAttributes[bench, HoldRest];

bench[label_String, expr_] := Module[{ts},
  expr;
  ts = Table[First[AbsoluteTiming[expr]], {3}];
  Print[StringPadRight[label, 52], ToString[Round[1000. Min[ts], 0.001]], " ms"]
];

check[label_String, value_] :=
  Print[StringPadRight[label, 52], "check = ", value];

(* ---- model -------------------------------------------------------------- *)

kfn  = 20000;      (* measurements                                          *)
kfdt = 0.01;       (* sample interval                                       *)

(* Constant-acceleration transition in 2D: state is
   {x, y, vx, vy, ax, ay}, so position picks up v dt + a dt^2/2 and velocity
   picks up a dt. *)
kfF = Table[Which[i == j,                 1.,
                  j == i + 2 && i <= 4,   kfdt,
                  j == i + 4 && i <= 2,   0.5 kfdt^2,
                  True,                   0.], {i, 1, 6}, {j, 1, 6}];
kfFT = Transpose[kfF];

kfH  = Table[If[i == j, 1., 0.], {i, 1, 2}, {j, 1, 6}];   (* observe x, y   *)
kfHT = Transpose[kfH];
kfQ  = Table[If[i == j, 0.001, 0.], {i, 1, 6}, {j, 1, 6}]; (* process noise *)
kfR  = Table[If[i == j, 0.05,  0.], {i, 1, 2}, {j, 1, 2}]; (* meas. noise   *)
kfI  = Table[If[i == j, 1., 0.], {i, 1, 6}, {j, 1, 6}];

kfx0 = Table[0., {6}];
kfP0 = Table[If[i == j, 1., 0.], {i, 1, 6}, {j, 1, 6}];

(* A deterministic "trajectory plus sensor noise". *)
kfz[k_] := {Sin[0.01 k] + 0.02 Sin[7.1 k], Cos[0.013 k] + 0.02 Cos[5.3 k]};

(* ---- 1. the scalar filter ----------------------------------------------- *)

(* Textbook predict/update.  Every matrix here is 6x6 or smaller, so none of
   them is anywhere near the packing threshold -- which is the finding this
   experiment exists to make; see README.md. *)
kfstep[st_, k_] := Module[{xp, pp, yy, ss, kk},
  xp = kfF . st[[1]];                       (* predict state                *)
  pp = kfF . st[[2]] . kfFT + kfQ;          (* predict covariance           *)
  yy = kfz[k] - kfH . xp;                   (* innovation                   *)
  ss = kfH . pp . kfHT + kfR;               (* innovation covariance, 2x2   *)
  kk = pp . kfHT . Inverse[ss];             (* Kalman gain, 6x2             *)
  {xp + kk . yy, (kfI - kk . kfH) . pp}];

kalman[m_] := Module[{st, k},
  st = {kfx0, kfP0}; k = 1;
  While[k <= m, st = kfstep[st, k]; k = k + 1];
  Tr[st[[2]]]];

(* ---- 2. the ensemble filter --------------------------------------------- *)

enn = 4096;        (* ensemble members                                      *)
ens = 200;         (* steps                                                 *)

enE0 = Table[0.01 Sin[1. (i + 6 j)], {i, 1, enn}, {j, 1, 6}];

(* The same estimator, array-shaped: the covariance is never formed from a
   model, it is estimated from the spread of the ensemble.  One 4096x6 GEMM
   replaces the 6x6 algebra, so this row lives in the opposite regime to the
   one above. *)
enstep[e_, k_] := Module[{pr, mn, an, pc, kk, inn},
  pr  = e . kfFT;                           (* advance every member         *)
  mn  = Total[pr]/enn;                      (* ensemble mean                *)
  an  = pr - Table[mn, {enn}];              (* anomalies                    *)
  pc  = Transpose[an] . an/(enn - 1.);      (* sample covariance, 6x6       *)
  kk  = pc . kfHT . Inverse[kfH . pc . kfHT + kfR];
  inn = Table[kfz[k], {enn}] - pr . kfHT;
  pr + inn . Transpose[kk]];

enkf[m_] := Module[{e, k},
  e = enE0; k = 1;
  While[k <= m, e = enstep[e, k]; k = k + 1];
  Total[Total[e]/enn]];

(* ---- run ---------------------------------------------------------------- *)

Print["Experiment 18 -- state estimation"];
Print[""];

bench["Kalman filter, 6 states, 20000 steps",     kalman[kfn]];
check["Kalman filter (2000 steps)", kalman[2000]];

bench["ensemble Kalman, 4096 members, 200 steps", enkf[ens]];
check["ensemble Kalman (20 steps)", enkf[20]];

(* ---- what a small matrix costs ------------------------------------------ *)
(*
   A 6x6 matrix is 36 elements, far below the packing threshold, so it is
   correctly never packed -- and Inverse then cannot reach the numeric path at
   all and runs a fraction-free symbolic elimination instead.  These are the
   numbers behind that claim.  For scale: a 6x6 dgetrf+dgetri is about 1
   microsecond.
*)
Print[""];
Print["-- what one small dense operation costs --"];

sm6 = Table[If[i == j, 2., 0.3], {i, 1, 6}, {j, 1, 6}];
sv6 = Table[N[i], {i, 1, 6}];
sm2 = {{2., 0.3}, {0.3, 1.7}};

bench["Inverse[2x2] x 20000",   Do[Inverse[sm2], {20000}]];
bench["Inverse[6x6] x 2000",    Do[Inverse[sm6], {2000}]];
bench["m6 . v6      x 50000",   Do[sm6 . sv6, {50000}]];
bench["m6 . m6      x 20000",   Do[sm6 . sm6, {20000}]];

(* For contrast, the same operation on data ABOVE the packing threshold,
   where it reaches LAPACK.  Per element this is orders of magnitude cheaper. *)
sm60 = Table[If[i == j, 60., 0.3], {i, 1, 60}, {j, 1, 60}];
bench["Inverse[60x60] x 2000 (above threshold)", Do[Inverse[sm60], {2000}]];
