(* ========================================================================== *)
(* goursat/CubicField.m -- Goursat's Theorem 2 when the resolvent cubic is     *)
(* irreducible over Q.                                                         *)
(*                                                                            *)
(* goursat.py decides the criterion here by the trace (goursatF0) but leaves   *)
(* the elementary part F - F0 as an unevaluated integral, because sympy is too *)
(* slow in a cubic field (paper, Limitations (2)).  Here the reduction is      *)
(* carried out.  With F - F0 = (1/2) Sum_theta Fm_theta,                       *)
(*   Fm_theta = (F - F o S_theta)/2  anti-invariant under S_theta,              *)
(* the anti-invariant piece is integrated once, over Q(theta1) with theta1 =   *)
(* Root[res, 1], by the same conic reduction as the rational case               *)
(* (goursatPiece).  Every step of that reduction is an identity in the field   *)
(* Q(theta) (and the quadratic extension the conic may need), so the formula    *)
(* with theta1 replaced by a variable is the piece for every conjugate, and     *)
(*   Integral (F - F0) dt/y = RootSum[res, piece(#)/2 &].                       *)
(* The answer stays defined over Q.  If theta1 does not survive as a clean      *)
(* generator of the piece, the three conjugate pieces are computed explicitly.  *)
(* The front end verifies the final antiderivative either way.                 *)
(* ========================================================================== *)

(* {L, N, M} of the involution indexed by the root r of the cubic. *)
involutionAt[trip_, th_, r_] := Module[{L, N1, M1},
  {L, N1, M1} = gK /@ (trip /. th -> r);
  If[!isZero[L], {1, gK[N1/L], gK[M1/L]}, {0, 1, gK[M1/N1]}]];

(* The anti-invariant piece for the involution at the root r. *)
cubicPiece[F_, R_, t_, y_, trip_, th_, r_] := Module[{L, N1, M1, S, Fm},
  {L, N1, M1} = involutionAt[trip, th, r];
  S = (N1 t - M1)/(L t - N1);
  Fm = gK[(F - (F /. t -> S))/2];
  If[Fm === 0, {0, <||>}, goursatPiece[Fm, L, N1, M1, R, t, y]]];

(* Reduce a rational function of v modulo res(v): numerator and denominator    *)
(* separately (an identity at every root of res).                               *)
fieldRat[b_, v_, rv_] := Module[{nd = numDen[b]},
  gK[polyrem[nd[[1]], rv, v].(v^Range[Exponent[rv, v] - 1, 0, -1]) /
     polyrem[nd[[2]], rv, v].(v^Range[Exponent[rv, v] - 1, 0, -1])]];

(* Content of a polynomial over Q, made positive. *)
posContent[p_, vars_] := Abs[Module[{c = Flatten[CoefficientList[Expand[p], vars]]},
  If[AllTrue[c, ratQ], (GCD @@ Numerator[c])/(LCM @@ Denominator[c]), 1]]];

(* Tidy the generic piece: every radicand, logarithm and inverse-trigonometric  *)
(* argument that depends on v is reduced modulo res(v) (the conjugate formula   *)
(* is unchanged at every root), and the positive integer content of a square   *)
(* root's radicand is pulled out, Sqrt[c p] = Sqrt[c] Sqrt[p] for c > 0.        *)
(* Structural, like tidyArgs: a rule's right-hand side would stay unevaluated  *)
(* inside the held body of the inner RootSum.                                   *)
fieldTidy[e_, v_, rv_] := Which[
  AtomQ[e] || FreeQ[e, v], e,
  Head[e] === Power && MatchQ[e[[2]], _Rational] && Denominator[e[[2]]] > 1,
    With[{b2 = fieldRat[fieldTidy[e[[1]], v, rv], v, rv], q = e[[2]]},
      With[{c = posContent[Numerator[b2], v]/posContent[Denominator[b2], v]},
        Power[c, q] Power[gK[b2/c], q]]],
  MatchQ[Head[e], Log | ArcTan | ArcTanh] && Length[e] == 1 && FreeQ[e[[1]], _Slot] &&
    FreeQ[e[[1]], Power[_, q_Rational /; Denominator[q] > 1]],
    Head[e][fieldRat[e[[1]], v, rv]],
  Head[e] === RootSum && Length[e] == 2, RootSum @@ {e[[1]], fieldTidy[e[[2]], v, rv]},
  Head[e] === Function, Function @@ Append[Most[List @@ e], fieldTidy[Last[e], v, rv]],
  True, fieldTidy[#, v, rv] & /@ e];

cubicFieldReduce[F_, R_, t_, y_, th_, res_, trip_] := Module[
  {F0, Fel, resFn, r1, piece, info, v, generic, total},
  F0 = stage["cubic field: F0", goursatF0[F, R, t]];   (* the criterion, decided over Q *)
  Fel = gK[F - F0];
  If[Fel === 0, Return[{0, F0, {}}, Module]];
  resFn = With[{b = res /. th -> Slot[1]}, Function[b]];
  r1 = Root[resFn, 1];
  {piece, info} = Replace[stage["cubic field: piece", Catch[cubicPiece[F, R, t, y, trip, th, r1], "GoursatFail"]],
                    _ggFail -> {$Failed, <||>}];
  generic = stage["cubic field: generic", If[piece === $Failed, $Failed, fieldTidy[piece /. r1 -> v, v, res /. th -> v]]];
  (* The conjugate formula is checked, not assumed: a closed form may carry    *)
  (* simplifications that hold only at theta1 (Sqrt[r^2] -> r for r > 0), and   *)
  (* then its conjugates are not antiderivatives of theirs.  A generic RootSum *)
  (* piece is decided by one exact check of the sum; a closed form (explicit   *)
  (* Log/ArcTan, whose sum over theta is a costly nested-RootSum check) by      *)
  (* computing and checking the three conjugate pieces.                         *)
  total = If[generic =!= $Failed && FreeQ[generic, _Root] && !FreeQ[generic, _RootSum],
    (* built by Apply: RootSum holds its arguments, and With does not rename *)
    (* inside a Function's parameter list                                     *)
    RootSum @@ {resFn, Function @@ {{v}, generic/2}},
    $Failed];
  If[total === $Failed ||
     !stage["cubic field: check RootSum",
       derivZeroQ[total, -Fel/y, t, y, R, 2]],
    total = stage["cubic field: explicit conjugates", Catch[
      Sum[cubicPiece[F, R, t, y, trip, th, Root[resFn, j]][[1]]/2, {j, 3}], "GoursatFail"]];
    (* no verified reduction: goursat.py's unevaluated integral (UNEVAL) *)
    If[MatchQ[total, _ggFail] || !stage["cubic field: check conjugates", derivZeroQ[total, -Fel/y, t, y, R, 2]],
      total = ggIntegral[Fel/y, t]]];
  {total, F0,
   {<|"note" -> "cubic-field reduction", "elementaryPart" -> Fel, "piece" -> info,
      "generic" -> FreeQ[generic, _Root]|>}}];
