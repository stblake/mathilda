(* ================================================================= *)
(* FeynCalc microbenchmarks: the same file in Mathilda and in the    *)
(* reference kernel.                                                 *)
(*                                                                   *)
(*   ./Mathilda -file tools/feyncalc/bench_fc.m                      *)
(*   wolframscript -file tools/feyncalc/bench_fc.m                   *)
(*                                                                   *)
(* Each case prints   BENCH <TAB> name <TAB> cold <TAB> warm  with   *)
(* times in microseconds: cold is the first call (FeynCalc's MemSet  *)
(* caches start empty), warm the best of the next three.             *)
(* tools/feyncalc/bench_fc.py runs both kernels and prints ratios.   *)
(* Set FC_BENCH_REPS to repeat each warm call more times.            *)
(* ================================================================= *)

$FeynCalcStartupMessages = False;
$fcRoot = Environment["FC_ROOT"];
If[!StringQ[$fcRoot], $fcRoot = FileNameJoin[{Directory[], "packages", "feyncalc"}]];
$t0 = AbsoluteTime[];
Get[FileNameJoin[{$fcRoot, "FeynCalc", "FeynCalc.m"}]];
WriteString["stdout", "BENCH\tload\t", ToString[Round[10^6 (AbsoluteTime[] - $t0)]], "\t0\n"];

$reps = ToExpression[Environment["FC_BENCH_REPS"]];
If[!IntegerQ[$reps], $reps = 3];

SetAttributes[bench, HoldRest];
bench[name_String, expr_] := Module[{t, cold, warm = Infinity, r},
  t = AbsoluteTime[]; r = expr; cold = AbsoluteTime[] - t;
  Do[t = AbsoluteTime[]; r = expr; warm = Min[warm, AbsoluteTime[] - t], {$reps}];
  WriteString["stdout", "BENCH\t", name, "\t", ToString[Round[10^6 cold]], "\t",
    ToString[Round[10^6 warm]], "\t", ToString[LeafCount[r]], "\n"]];

(* Dirac algebra *)
bench["trace8", DiracSimplify[DiracTrace[GA[a1, a2, a3, a4, a5, a6, a7, a8]]]];
bench["traceSlashMass", DiracTrace[(GS[p] + m).GA[mu].(GS[q] + m).GA[nu].(GS[k] + m).GA[mu].GA[nu],
  DiracTraceEvaluate -> True]];
bench["traceG5", DiracTrace[GA[a, b, c, d, e, f, 5], DiracTraceEvaluate -> True]];
bench["diracSimplifyD", DiracSimplify[GAD[mu, nu, rho, sigma, tau, mu]]];
bench["diracChain", DiracSimplify[SpinorUBar[p1, m].GS[p2].GA[mu].GS[p1].GA[nu].SpinorU[p3, m]]];

(* Lorentz algebra *)
bench["contractLC", Contract[LC[mu, nu, rho, si] LC[mu, nu, al, be] FV[p, rho] FV[q, si] FV[k, al] FV[l, be]]];
bench["expandSP", ExpandScalarProduct[SP[p1 + p2 + p3 + p4 - k, p1 - p2 + p3 - p4 + k]]];
bench["contractChain", Contract[MT[mu, nu] MT[nu, rho] MT[rho, si] FV[p, mu] FV[q, si] +
  FV[p, mu] FV[k, mu] FV[q, nu] FV[l, nu]]];

(* Colour algebra *)
bench["sunTrace", SUNSimplify[SUNTrace[SUNT[a, b, c, a, b, c]]]];
bench["sunF", SUNSimplify[SUNF[a, b, c] SUNF[a, b, d] SUNF[c, e, f] SUNF[d, e, f]]];

(* Loop integrals *)
bench["tid2", TID[FVD[p, mu] FVD[p, nu] FAD[{p, m}, p - q], p]];
bench["tid3", TID[FVD[p, mu] FVD[p, nu] FVD[p, rho] FAD[{p, m0}, {p + q, m1}], p]];
bench["tdec", Tdec[{{q, mu}, {q, nu}, {q, rho}}, {p1, p2}]];
bench["apartFF", ApartFF[SPD[q, p] FAD[{q, m}, q - p, q + p], {q}]];
bench["fds", FDS[FAD[q, q - p, q + p, {q - 2 p, m}], q]];
bench["paveReduce", PaVeReduce[PaVe[1, 1, {SPD[p, p]}, {m0^2, m1^2}]]];
bench["toPakForm", FCLoopToPakForm[FAD[{q1, m}, {q2, m}, q1 - q2, q1 - p, q2 - p], {q1, q2}]];

WriteString["stdout", "BENCHEND\n"];
