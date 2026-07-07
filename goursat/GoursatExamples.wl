(* ::Package:: *)

(* ===================================================================== *)
(*  GoursatExamples.wl                                                   *)
(*                                                                       *)
(*  Worked examples illustrating the GoursatPseudoElliptic package.      *)
(*  Run with                                                             *)
(*     Get["GoursatPseudoElliptic.wl"];                                  *)
(*     Get["GoursatExamples.wl"];                                        *)
(*  or paste cell-by-cell into a notebook.                               *)
(*                                                                       *)
(*  All examples below are pseudo-elliptic via the Klein 4-group V4 of   *)
(*  Theorem 2, which is what GoursatPseudoElliptic implements.  The      *)
(*  paper also discusses larger symmetry groups (period-3 in Sec. 4,     *)
(*  tetrahedral in Sec. 5, octahedral in Sec. 6); those would require    *)
(*  separate routines and are noted at the end.                          *)
(* ===================================================================== *)

Get["GoursatPseudoElliptic.wl"];

(* --------------------------------------------------------------------- *)
Print["================================================================="];
Print["  Example 1.  R(t) = (t^2-1)(t^2-4),  F(t) = t."];
Print["              Trivial Theorem 1 case: F is odd under S = -t."];
Print["================================================================="];

R1 = (t^2 - 1)(t^2 - 4);
F1 = t;

Print["Roots:                ", QuarticRoots[R1, t]];
Print["Three involutions:    ", Sigma = GoursatInvolutions[R1, t]];
Print["Goursat sum F+F(S1)+F(S2)+F(S3):"];
Print["  ", Together[F1 + Total[(F1 /. t -> #) & /@ Sigma]]];
Print["Pseudo-elliptic?      ", PseudoEllipticQ[F1, R1, t]];

ans1 = GoursatPseudoElliptic[F1, R1, t];
Print["Antiderivative:       ", ans1];
Print["Differentiation check (should be 0):"];
Print["  ", FullSimplify[D[ans1, t] - F1/Sqrt[R1]]];

(* --------------------------------------------------------------------- *)
Print[""];
Print["================================================================="];
Print["  Example 2.  Theorem 1 with FINITE fixed points."];
Print["              R(t) = (t^2-1)(t^2-4),  F(t) = (t^2-2)/t."];
Print["              F is anti-invariant under  S(t) = 2/t."];
Print["================================================================="];

R2 = (t^2 - 1)(t^2 - 4);
F2 = (t^2 - 2)/t;
S2 = 2/t;

Print["F + F(S):            ", Together[F2 + (F2 /. t -> S2)]];
Print["AntiInvariantQ:      ", AntiInvariantQ[F2, S2, t]];
Print["Fixed points of S:   ", InvolutionFixedPoints[S2, t]];

ans2 = GoursatReduceTheorem1[F2, R2, S2, t];
Print["Antiderivative:      ", ans2];
Print["Differentiation check:"];
Print["  ", FullSimplify[D[ans2, t] - F2/Sqrt[R2]]];

(* --------------------------------------------------------------------- *)
Print[""];
Print["================================================================="];
Print["  Example 3.  Symbolic modulus k.  R = (1-t^2)(1-k^2 t^2)."];
Print["              F(t) = t - 1/(k^2 t).                             "];
Print["              The V4 condition is f(t)+f(-t)+f(1/(kt))+f(-1/(kt))=0."];
Print["================================================================="];

R3 = (1 - t^2)(1 - k^2 t^2);
F3 = t - 1/(k^2 t);

Print["Pseudo-elliptic?     ", PseudoEllipticQ[F3, R3, t]];
ans3 = GoursatPseudoElliptic[F3, R3, t];
Print["Antiderivative:      ", ans3];
Print["Differentiation check:"];
Print["  ", FullSimplify[D[ans3, t] - F3/Sqrt[R3]]];

(* --------------------------------------------------------------------- *)
Print[""];
Print["================================================================="];
Print["  Example 4.  Theorem 2 with TWO components nonzero."];
Print["              R(t) = (t^2-1)(t^2-4),  F(t) = t + (t^2-2)/t."];
Print["              Linearity: each piece is pseudo-elliptic."];
Print["================================================================="];

R4 = (t^2 - 1)(t^2 - 4);
F4 = t + (t^2 - 2)/t;

Print["Pseudo-elliptic?     ", PseudoEllipticQ[F4, R4, t]];
ans4 = GoursatPseudoElliptic[F4, R4, t];
Print["Antiderivative:      ", ans4];
Print["Differentiation check:"];
Print["  ", FullSimplify[D[ans4, t] - F4/Sqrt[R4]]];

(* --------------------------------------------------------------------- *)
Print[""];
Print["================================================================="];
Print["  Example 5.  NEGATIVE example.  R = (t^2-1)(t^2-4),  F = t^2."];
Print["              t^2 is V4-invariant, so the sum is 4 t^2 != 0."];
Print["              The integral is genuinely elliptic, not elementary."];
Print["================================================================="];

R5 = (t^2 - 1)(t^2 - 4);
F5 = t^2;

Print["V4 sum F+F(S1)+F(S2)+F(S3):"];
Print["  ", Together[F5 + Total[(F5 /. t -> #) & /@ GoursatInvolutions[R5, t]]]];
Print["Pseudo-elliptic?     ", PseudoEllipticQ[F5, R5, t]];
Print["GoursatPseudoElliptic should refuse:"];
Quiet[Check[
  ans5 = GoursatPseudoElliptic[F5, R5, t];
  Print["  Result: ", ans5],
  Print["  (correctly returned $Failed)"]
]];

(* --------------------------------------------------------------------- *)
Print[""];
Print["================================================================="];
Print["  Example 6.  Cubic case.  R(t) = (t-1)(t-2)(t-3)."];
Print["              Treat 4th root as Infinity.                       "];
Print["              The V4 involutions pair these as                  "];
Print["              (1,2)|(3,Inf), (1,3)|(2,Inf), (1,Inf)|(2,3).      "];
Print["================================================================="];

R6 = (t - 1)(t - 2)(t - 3);
Print["Roots:               ", QuarticRoots[R6, t]];
Print["Involutions:         ", GoursatInvolutions[R6, t]];

(* Pick the involution swapping (1,Inf) and (2,3) *)
S6 = PairSwapInvolution[2, 3, 1, Infinity, t];
Print["Picked S6:           ", S6];
Print["Fixed points of S6:  ", InvolutionFixedPoints[S6, t]];

(* Construct F via the standard antisymmetrization F(t) = g(t) - g(S6(t))
   for any rational g.  Take g(t) = t. *)
F6 = Together[t - (t /. t -> S6)];
Print["F6 = t - S6(t):      ", F6];
Print["AntiInvariantQ:      ", AntiInvariantQ[F6, S6, t]];

ans6 = GoursatReduceTheorem1[F6, R6, S6, t];
Print["Antiderivative:      ", ans6];
Print["Differentiation check:"];
Print["  ", FullSimplify[D[ans6, t] - F6/Sqrt[R6]]];

(* --------------------------------------------------------------------- *)
Print[""];
Print["================================================================="];
Print["  Notes on what is NOT covered by this implementation."];
Print["================================================================="];
Print["The package implements Theorems 1 and 2 of Goursat (1887), which"];
Print["use only the Klein four-group V4 of root-pairing involutions."];
Print["Sections 4-6 of the paper extend the theory:"];
Print["  Sec. 4: Period-3 substitutions for R = t^3 - 1, e.g. the"];
Print["          Legendre/Clausen integral"];
Print["          Integrate[t/((t^3+8)Sqrt[t^3-1]), t]."];
Print["  Sec. 5: The full tetrahedral group A4 (12 elements) for t^3-1."];
Print["  Sec. 6: Octahedral cases for R = t(t^2+1) or R = t^4-1."];
Print["These higher-symmetry cases are NOT handled by GoursatPseudoElliptic."];
Print["For those, an extension implementing the period-3 character"];
Print["projections F(S) = alpha F  (alpha = e^(2 pi i/3))  is required;"];
Print["see Section 4 of the paper for the algorithm."];
