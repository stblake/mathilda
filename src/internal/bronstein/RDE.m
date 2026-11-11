(* ::Package:: *)

(* ========================================================================== *)
(* RDE.m -- the tri-state Risch differential equation  D y + f y = g  over    *)
(*          the tower Q(x)(t_1..t_r), wrapping Risch`RischDECertified.        *)
(*                                                                            *)
(* Every certificate the package issues that rests on "this Risch DE has no   *)
(* solution" ([B] Prop 2.18 and the special-place reductions of [T] App A)    *)
(* must get that "no" from here, and only from the {"NoSolution"} answer:    *)
(* Risch`RischDE itself stays unevaluated both when no solution exists and    *)
(* when a sub-case is unimplemented, so it cannot back a certificate.  The    *)
(* C entry point runs inside a certification scope and answers "NoSolution"   *)
(* only when every decision on the path was a theorem of Bronstein Ch.6       *)
(* (src/calculus/integrate_risch_rde.h, "Certification").                     *)
(*                                                                            *)
(* The solution field is the WHOLE tower: the kernels are passed explicitly,  *)
(* since y' = 1/x has no solution in Q(x) but has Log[x] in Q(x, Log[x]).     *)
(* ========================================================================== *)

brRischDE::usage = "brRischDE[f, g, x, kernels] solves D[y] + f y == g for y in Q(x)(kernels) (kernels: the list of Log/Exp kernels of the tower, in surface form, e.g. {Log[x], E^x}; {} for Q(x)). It returns {\"Solution\", y} (y verified by the exact identity), {\"NoSolution\"} (a proof that none exists in that field), or {\"Unknown\", reason} (undecided: an unimplemented sub-case or an input outside rational functions of x and the kernels).";
brRischDE[f_, g_, x_Symbol, kernels_List] := Module[{r, why},
  r = If[kernels === {}, Risch`RischDECertified[f, g, x], Risch`RischDECertified[f, g, x, kernels]];
  Which[
    r === "NoSolution", {"NoSolution"},
    r === $Failed,
      why = Risch`RischDECertifiedReason[];
      {"Unknown", If[StringQ[why], why, "undecided"]},
    Head[r] === Risch`RischDECertified, {"Unknown", "unevaluated"},
    TrueQ[brRDEResidualZeroQ[f, g, x, r]], {"Solution", r},
    True, throwInternal["brRischDEWrongSolution", <|"f" -> f, "g" -> g, "y" -> r|>]]];

brRDEResidualZeroQ::usage = "brRDEResidualZeroQ[f, g, x, y] is True when D[y] + f y - g is exactly zero (Together of the expanded residual; f, g, y rational in x and Log/Exp kernels, whose derivatives D computes).";
brRDEResidualZeroQ[f_, g_, x_Symbol, y_] :=
  Together[Expand[Together[D[y, x] + f y - g]]] === 0 || PossibleZeroQ[D[y, x] + f y - g];

brRDENoSolutionQ::usage = "brRDENoSolutionQ[r] is True when the brRischDE answer r is the certified {\"NoSolution\"}.";
brRDENoSolutionQ[r_] := r === {"NoSolution"};

brRDESolution::usage = "brRDESolution[r] is the solution y of the brRischDE answer r = {\"Solution\", y}, or None.";
brRDESolution[{"Solution", y_}] := y;
brRDESolution[_] := None;
