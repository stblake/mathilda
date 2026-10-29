(* ::Package:: *)

(* ========================================================================== *)
(* PrincipalGen.m -- the principality test of [T] Ch 6 Sec 1 pp.64-65.        *)
(*                                                                            *)
(* Let delta be a divisor of degree 0 with no place at infinity, and I its    *)
(* ideal.  An element u of I has ord_P u >= delta_P at every finite place; if *)
(* moreover u has no pole at infinity then div(u) >= delta, and both having   *)
(* degree 0 forces div(u) = delta.  So delta is principal iff I contains a    *)
(* nonzero element regular at infinity, and that element is the generator     *)
(* (unique up to a constant).  With a basis of I normal at infinity such an   *)
(* element exists iff some basis row is itself regular at infinity (order     *)
(* k_i >= 0), because the order of Sum p_i row_i is the least order of its    *)
(* terms.  The test is exact: no heuristic verification is needed.            *)
(* ========================================================================== *)

principalGenerator::usage = "principalGenerator[F, B, I] decides whether the degree-0 divisor whose ideal is I (support at finite places only) is principal. It returns the generator u (coordinate list, div(u) equal to the divisor) or None.";
principalGenerator[F_Association, B_Association, I_Association] := Module[{x = B["x"], nrm, pos},
  nrm = normalizeAtInfinity[F, B, I["M"], I["den"], x];
  pos = SelectFirst[Range[Length[nrm["rows"]]], nrm["k"][[#]] >= 0 &, None];
  If[pos === None, None,
    rc[F, #/I["den"]] & /@ nrm["rows"][[pos]]]];

principalQ::usage = "principalQ[F, B, I] is True when the degree-0 divisor with ideal I is principal.";
principalQ[F_Association, B_Association, I_Association] := principalGenerator[F, B, I] =!= None;
