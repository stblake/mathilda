(* ::Package:: *)

(* ===================================================================== *)
(*  GoursatPseudoElliptic.wl                                             *)
(*                                                                       *)
(*  Implementation of E. Goursat's 1887 theory of pseudo-elliptic        *)
(*  integrals of the form  Integrate[F[t]/Sqrt[R[t]], t]                 *)
(*  where R is a polynomial of degree 3 or 4 with simple roots.          *)
(*                                                                       *)
(*  Reference:                                                           *)
(*    E. Goursat, "Note sur quelques integrales pseudo-elliptiques",     *)
(*    Bulletin de la S. M. F., t. 15 (1887), p. 106-120.                 *)
(*                                                                       *)
(*  Usage:                                                               *)
(*    Get["GoursatPseudoElliptic.wl"]                                    *)
(*    GoursatPseudoElliptic[t, (t^2-1)(t^2-4), t]                        *)
(* ===================================================================== *)

BeginPackage["GoursatPseudoElliptic`"];

PairSwapInvolution::usage =
  "PairSwapInvolution[a, b, c, d, t] returns the Mobius involution \
in the variable t that interchanges {a, b} <-> {c, d}, i.e. the unique \
non-trivial involution sending a<->b and c<->d.  Any one of a, b, c, d \
may be Infinity.";

QuarticRoots::usage =
  "QuarticRoots[R, t] returns the four roots of R(t).  If R(t) is cubic, \
Infinity is appended so the root list always has length 4.";

GoursatInvolutions::usage =
  "GoursatInvolutions[R, t] returns the list {S1, S2, S3} of three Mobius \
involutions of CP^1 that permute the four roots of R(t) in pairs \
(Goursat 1887, Section 3).  These three involutions, together with the \
identity, form a Klein four-group V4 acting on {a,b,c,d}.";

InvolutionFixedPoints::usage =
  "InvolutionFixedPoints[S, t] returns the two fixed points of the Mobius \
involution S(t), one of which may be Infinity.";

PseudoEllipticQ::usage =
  "PseudoEllipticQ[F, R, t] is True iff Integrate[F[t]/Sqrt[R[t]], t] \
satisfies Goursat's pseudo-elliptic criterion (Theorem 2): \
F + F[S1] + F[S2] + F[S3] == 0.";

AntiInvariantQ::usage =
  "AntiInvariantQ[F, S, t] is True iff F(t) + F(S(t)) == 0.  This is the \
hypothesis of Goursat's Theorem 1.";

GoursatReduceTheorem1::usage =
  "GoursatReduceTheorem1[F, R, S, t] applies Goursat's two-step \
substitution u = (t-alpha)/(t-beta) followed by x = u^2 to a Theorem-1 \
integrand (F + F[S] = 0) and returns an elementary antiderivative of \
F[t]/Sqrt[R[t]] in t.";

GoursatPseudoElliptic::usage =
  "GoursatPseudoElliptic[F, R, t] returns an elementary antiderivative of \
F[t]/Sqrt[R[t]] when the integrand satisfies Goursat's criterion \
(Theorem 2).  Internally, F is decomposed as F = F1 + F2 + F3 where each \
Fi is anti-invariant under at least one involution Si, and Theorem 1 is \
applied to each piece.  Returns $Failed if F is not pseudo-elliptic.";


Begin["`Private`"];


(* ====================================================================== *)
(*  Helpers                                                                *)
(* ====================================================================== *)

(* Robust zero test for rational/algebraic expressions *)
zeroQ[expr_] := Module[{e},
  e = Together[expr];
  If[PossibleZeroQ[e], Return[True]];
  e = Simplify[e];
  If[PossibleZeroQ[e], Return[True]];
  e = FullSimplify[e];
  PossibleZeroQ[e]
];

(* Substitute  u^2 -> x  in an expression that should depend only on u^2.
   Implemented as  u -> Sqrt[x] followed by PowerExpand and Together.    *)
substUSquared[expr_, u_, x_] :=
  Together[PowerExpand[expr /. u -> Sqrt[x]]];


(* ====================================================================== *)
(*  1.  Mobius involutions that swap a pair of pairs of points.           *)
(*      The general formula (Goursat, Section 3) is                       *)
(*         S(t) = ((ab-cd)t + (a+b)cd - (c+d)ab)                          *)
(*                / ((a+b-c-d)t - (ab-cd))                                *)
(*      which sends a<->b and c<->d.                                      *)
(* ====================================================================== *)

PairSwapInvolution[a_, b_, c_, d_, t_] /; FreeQ[{a, b, c, d}, Infinity] :=
  With[{
        LL = (a + b) - (c + d),
        NN = a*b - c*d,
        MM = a*b*(c + d) - c*d*(a + b)
      },
    Together[(NN*t - MM)/(LL*t - NN)]
  ];

(* Limit forms when one of the four points is at infinity.
   Take the formula above and let one root tend to Infinity.            *)
PairSwapInvolution[a_, b_, c_, Infinity, t_] :=
  Together[(c*t + a*b - c*(a + b))/(t - c)];
PairSwapInvolution[a_, b_, Infinity, c_, t_] :=
  PairSwapInvolution[a, b, c, Infinity, t];
PairSwapInvolution[Infinity, a_, b_, c_, t_] :=
  PairSwapInvolution[b, c, a, Infinity, t];
PairSwapInvolution[a_, Infinity, b_, c_, t_] :=
  PairSwapInvolution[b, c, a, Infinity, t];


(* ====================================================================== *)
(*  2.  Roots of R(t).                                                     *)
(* ====================================================================== *)

QuarticRoots[R_, t_Symbol] := Module[{rts, deg = Exponent[R, t]},
  If[!MemberQ[{3, 4}, deg], Return[$Failed]];
  rts = DeleteDuplicates[t /. Solve[R == 0, t]];
  Which[
    Length[rts] == 4, rts,
    Length[rts] == 3, Append[rts, Infinity],
    True, $Failed
  ]
];


(* ====================================================================== *)
(*  3.  The three involutions S1, S2, S3 of Goursat's Section 3.           *)
(*      They correspond to the three pairings of 4 points:                 *)
(*           {(a,b),(c,d)},  {(a,c),(b,d)},  {(a,d),(b,c)}.                *)
(* ====================================================================== *)

GoursatInvolutions[R_, t_Symbol] := Module[{rts},
  rts = QuarticRoots[R, t];
  If[rts === $Failed, Return[$Failed]];
  {
    PairSwapInvolution[rts[[1]], rts[[2]], rts[[3]], rts[[4]], t],
    PairSwapInvolution[rts[[1]], rts[[3]], rts[[2]], rts[[4]], t],
    PairSwapInvolution[rts[[1]], rts[[4]], rts[[2]], rts[[3]], t]
  }
];


(* ====================================================================== *)
(*  4.  Fixed points of an involution.                                     *)
(* ====================================================================== *)

InvolutionFixedPoints[S_, t_Symbol] := Module[
  {poly, sol, deg},
  poly = Expand[Numerator[Together[S - t]]];
  deg  = Exponent[poly, t];
  sol  = t /. Solve[poly == 0, t];
  Which[
    deg == 2 && Length[sol] == 2, sol,
    deg == 2 && Length[sol] == 1, {sol[[1]], sol[[1]]},  (* repeated -- not an involution *)
    deg == 1 && Length[sol] == 1, Append[sol, Infinity], (* one fixed point at infinity *)
    deg == 0,                     {Infinity, Infinity},
    True,                          $Failed
  ]
];


(* ====================================================================== *)
(*  5.  Tests for the Goursat conditions.                                  *)
(* ====================================================================== *)

PseudoEllipticQ[F_, R_, t_Symbol] := Module[{Sigma, sum},
  Sigma = GoursatInvolutions[R, t];
  If[Sigma === $Failed, Return[False]];
  sum = F + Total[(F /. t -> #) & /@ Sigma];
  zeroQ[sum]
];

AntiInvariantQ[F_, S_, t_Symbol] := zeroQ[F + (F /. t -> S)];


(* ====================================================================== *)
(*  6.  Theorem-1 reduction                                                *)
(*                                                                         *)
(*    Given an integrand F(t)/Sqrt[R(t)] satisfying F(t) + F(S(t)) = 0,    *)
(*    Goursat's reduction is:                                              *)
(*                                                                         *)
(*    (a) Find the two fixed points alpha, beta of S.  Substitute          *)
(*           u = (t - alpha)/(t - beta),                                   *)
(*           t = (u beta - alpha)/(u - 1).                                 *)
(*        Under this map, S becomes  u -> -u, so F(t(u)) becomes odd in u  *)
(*        and R(t(u)) (u-1)^4 becomes a polynomial in u^2.                 *)
(*                                                                         *)
(*    (b) Substitute  x = u^2.  Then  u du = dx/2, and the integrand       *)
(*        reduces to                                                       *)
(*           (1/2) G(x) / Sqrt[H(x)] dx                                    *)
(*        where H is at most quadratic in x.  This is elementary.          *)
(*                                                                         *)
(*    Special case beta = Infinity: S is then  t -> 2 alpha - t,           *)
(*        and we use u = t - alpha directly.                               *)
(* ====================================================================== *)

GoursatReduceTheorem1[F_, R_, S_, t_Symbol] := Module[
  {fps, alpha, beta, u, x, tu, FU, RUeven, GofX, HofX,
   integrand, antider, result},

  (* (1) Fixed points of S *)
  fps = InvolutionFixedPoints[S, t];
  If[fps === $Failed || Length[fps] != 2,
    Message[GoursatReduceTheorem1::nofp]; Return[$Failed]
  ];
  {alpha, beta} = fps;
  (* keep alpha finite when possible *)
  If[alpha === Infinity, {alpha, beta} = {beta, alpha}];

  (* (2) Substitute  t = t(u) *)
  Which[
    beta === Infinity,
      (* S is t -> 2 alpha - t.  Use u = t - alpha. *)
      tu     = u + alpha;
      FU     = Together[F /. t -> tu];           (* should be odd in u   *)
      RUeven = Together[R /. t -> tu]            (* should be even in u  *)
    ,
    True,
      (* Both fixed points finite. *)
      tu     = (u*beta - alpha)/(u - 1);
      (* dt/du = (alpha - beta)/(u-1)^2.  After clearing the common
         (u-1)^2 in numerator and the (u-1)^2 from Sqrt[R(t(u))],
         the effective integrand is (alpha-beta) F(t(u)) / Sqrt[R(t(u))(u-1)^4]. *)
      FU     = Together[(F /. t -> tu) (alpha - beta)];
      RUeven = Cancel[Together[(R /. t -> tu) (u - 1)^4]]
  ];

  (* (3) Substitute x = u^2.  FU/u should be a function of u^2 only,
         and so should RUeven. *)
  GofX = substUSquared[Together[FU/u], u, x];
  HofX = substUSquared[RUeven, u, x];

  If[!FreeQ[GofX, u] || !FreeQ[HofX, u],
    Message[GoursatReduceTheorem1::notpe]; Return[$Failed]
  ];

  (* (4) Reduced elementary integrand *)
  integrand = (1/2) GofX / Sqrt[HofX];
  antider   = Integrate[integrand, x];

  (* (5) Substitute back x = u^2 in terms of t *)
  result = Which[
    beta === Infinity, antider /. x -> (t - alpha)^2,
    True,              antider /. x -> ((t - alpha)/(t - beta))^2
  ];
  Together[result]
];

GoursatReduceTheorem1::nofp =
  "Could not determine the two fixed points of the involution.";
GoursatReduceTheorem1::notpe =
  "Anti-invariance check failed: the integrand is not Theorem-1 \
pseudo-elliptic for this involution.";


(* ====================================================================== *)
(*  7.  Theorem-2 (top-level) reduction.                                   *)
(*                                                                         *)
(*    Hypothesis:  F + F(S1) + F(S2) + F(S3) = 0.                          *)
(*                                                                         *)
(*    Decompose F into characters of the Klein four-group:                 *)
(*       F1 = (F - F(S1) - F(S2) + F(S3))/4   (character (+,-,-,+))        *)
(*       F2 = (F - F(S1) + F(S2) - F(S3))/4   (character (+,-,+,-))        *)
(*       F3 = (F + F(S1) - F(S2) - F(S3))/4   (character (+,+,-,-))        *)
(*    Then F = F1 + F2 + F3 (since the trivial character F+F(S1)+...=0).   *)
(*                                                                         *)
(*    Anti-invariance check (working out the V4 action carefully):         *)
(*       F1 is anti-inv under S1 and S2  (invariant under S3).             *)
(*       F2 is anti-inv under S1 and S3  (invariant under S2).             *)
(*       F3 is anti-inv under S2 and S3  (invariant under S1).             *)
(*    NOTE: This corrects an indexing slip in the proof of Theorem 2 in    *)
(*    the manuscript -- F2 is invariant under S2 (not anti-invariant).    *)
(*    The conclusion of the theorem is unaffected.                         *)
(*                                                                         *)
(*    For each component we apply Theorem 1 with a valid involution.       *)
(* ====================================================================== *)

GoursatPseudoElliptic[F_, R_, t_Symbol] := Module[
  {Sigma, S1, S2, S3, FS1, FS2, FS3, sum, F1, F2, F3, I1, I2, I3, total},

  Sigma = GoursatInvolutions[R, t];
  If[Sigma === $Failed,
    Message[GoursatPseudoElliptic::badpoly]; Return[$Failed]
  ];
  {S1, S2, S3} = Sigma;

  FS1 = Together[F /. t -> S1];
  FS2 = Together[F /. t -> S2];
  FS3 = Together[F /. t -> S3];
  sum = Together[F + FS1 + FS2 + FS3];

  If[!zeroQ[sum],
    Message[GoursatPseudoElliptic::nopseudo]; Return[$Failed]
  ];

  F1 = Together[(F - FS1 - FS2 + FS3)/4];   (* anti-inv under S1, S2 *)
  F2 = Together[(F - FS1 + FS2 - FS3)/4];   (* anti-inv under S1, S3 *)
  F3 = Together[(F + FS1 - FS2 - FS3)/4];   (* anti-inv under S2, S3 *)

  I1 = If[zeroQ[F1], 0, GoursatReduceTheorem1[F1, R, S1, t]];
  I2 = If[zeroQ[F2], 0, GoursatReduceTheorem1[F2, R, S1, t]];
  I3 = If[zeroQ[F3], 0, GoursatReduceTheorem1[F3, R, S2, t]];

  If[MemberQ[{I1, I2, I3}, $Failed], Return[$Failed]];

  total = Together[I1 + I2 + I3];
  total
];

GoursatPseudoElliptic::nopseudo =
  "F + F(S1) + F(S2) + F(S3) is not identically zero -- the integrand \
is not pseudo-elliptic by Goursat's criterion.";
GoursatPseudoElliptic::badpoly =
  "R(t) must be a polynomial of degree 3 or 4 in t with simple roots.";


End[]; (* `Private` *)
EndPackage[];
