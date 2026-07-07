(* ::Package:: *)

(* ===================================================================== *)
(*  CubicRootGoursat.wl                                                  *)
(*                                                                       *)
(*  A cube-root analog of Goursat's pseudo-elliptic reduction:           *)
(*  for integrals of the form  Integrate[F[t]/R[t]^(1/3), t]             *)
(*  where R is a cubic polynomial with simple roots.                     *)
(* ===================================================================== *)

BeginPackage["CubicRootGoursat`"];

CyclicMobius::usage =
  "CyclicMobius[a, b, c, t] returns the unique Mobius transformation \
of order 3 sending a -> b -> c -> a.";

CubicRootProject::usage =
  "CubicRootProject[F, S, t, k] returns the eigencomponent of F under \
the cyclic substitution S having eigenvalue omega^k (omega = e^(2 pi i/3)). \
The argument k must be 0, 1, or 2.";

CubicRootElementaryQ::usage =
  "CubicRootElementaryQ[F, S, t] is True iff the omega^1 eigencomponent \
of F under S vanishes -- equivalently, iff Integrate[F/R^(1/3), t] is \
elementary by the cube-root Goursat theorem (with R determined by S).";

CubicRootReduce::usage =
  "CubicRootReduce[F, R, S, t] reduces an integral whose F has no \
omega-eigencomponent under S to elementary form, returning an antiderivative \
of F[t]/R[t]^(1/3) in t.  Returns $Failed if the obstruction does not \
vanish.";


Begin["`Private`"];

omega = Exp[2 Pi I/3];

(* --- 1. Build the cyclic Mobius S with S(a)=b, S(b)=c, S(c)=a --- *)
CyclicMobius[a_, b_, c_, t_] := Module[{A, B, C, sol},
  sol = Solve[
    {(A a + B)/(C a + 1) == b,
     (A b + B)/(C b + 1) == c,
     (A c + B)/(C c + 1) == a},
    {A, B, C}
  ];
  If[sol === {}, Return[$Failed]];
  Together[(A t + B)/(C t + 1) /. sol[[1]]]
];

(* --- 2. eigencomponent of F under S, eigenvalue omega^k --- *)
CubicRootProject[F_, S_, t_, k_Integer] := Module[{FS, S2, FS2},
  FS = F /. t -> S;
  S2 = Together[S /. t -> S];
  FS2 = F /. t -> S2;
  Together[Simplify[(F + omega^(-k) FS + omega^(-2 k) FS2)/3]]
];

(* --- 3. Test whether F is "cube-pseudo-elementary" under S --- *)
CubicRootElementaryQ[F_, S_, t_] :=
  PossibleZeroQ[Together[Simplify[CubicRootProject[F, S, t, 1]]]];

(* --- 4. Full reduction.  S is assumed of order 3 permuting roots of R. --- *)
CubicRootReduce[F_, R_, S_, t_] := Module[
  {fps, alpha, beta, z, x, w, tz, dt, F0, F2, Fz, RzPoly, K, c,
   integrand0, integrand2, anti0, anti2, total},

  (* fixed points of S *)
  fps = t /. Solve[Together[Numerator[Together[S - t]]] == 0, t];
  If[Length[fps] != 2, Return[$Failed]];
  {alpha, beta} = fps;

  (* substitute t = (alpha - beta z)/(1 - z), so z = (t-alpha)/(t-beta) *)
  tz = If[beta === Infinity, z + alpha, (alpha - beta z)/(1 - z)];
  dt = D[tz, z];

  (* Check the omega^1 obstruction first *)
  If[!CubicRootElementaryQ[F, S, t],
    Message[CubicRootReduce::middle]; Return[$Failed]
  ];

  (* Get k=0 and k=2 components *)
  F0 = CubicRootProject[F, S, t, 0];
  F2 = CubicRootProject[F, S, t, 2];

  (* Transform R into z-coordinates: R(t(z)) (1-z)^3 = c (z^3 - K) *)
  Module[{Rz},
    Rz = Together[(R /. t -> tz)*If[beta === Infinity, 1, (1 - z)^3]];
    Rz = Cancel[Rz];
    (* Should be of the form c (z^3 - K) -- read off c, K *)
    {K, c} = If[beta === Infinity,
      {-Coefficient[Rz, z, 0]/Coefficient[Rz, z, 3], Coefficient[Rz, z, 3]},
      {-Coefficient[Expand[Rz], z, 0]/Coefficient[Expand[Rz], z, 3],
        Coefficient[Expand[Rz], z, 3]}
    ];
  ];

  (* k=0 piece: F0(t(z)) dt/R^(1/3) -> phi0(x) dx/(3 (x^2(x-K))^(1/3))
     then x = K w^3/(w^3-1) makes it -1/(w^3-1) d w. *)
  Module[{phi0, phi2, ans0, ans2},
    phi0 = Together[(F0 /. t -> tz) * dt /
                    If[beta === Infinity, 1, 1/(1-z)] * c^(1/3)];
    (* phi0 should be a function of z^3; extract phi0(z^3 -> x) *)
    phi0 = phi0 /. z -> x^(1/3);
    phi0 = Together[PowerExpand[phi0]];

    phi2 = Together[(F2 /. t -> tz) * dt /
                    If[beta === Infinity, 1, 1/(1-z)] * c^(1/3)];
    phi2 = phi2 /. z -> x^(1/3);
    phi2 = Together[PowerExpand[phi2]];

    (* Hand off to Mathematica's Integrate at the elementary level *)
    ans0 = Integrate[phi0/(x^2 (x - K))^(1/3), x];
    ans2 = Integrate[phi2/(x - K)^(1/3), x];

    total = ans0 + ans2;
    (* substitute back x = z^3 = ((t-alpha)/(t-beta))^3 *)
    total = total /. x -> If[beta === Infinity, (t - alpha)^3,
                                                 ((t - alpha)/(t - beta))^3];
    Together[total]
  ]
];

CubicRootReduce::middle = "The omega-eigencomponent (k=1) of F is nonzero; \
the integral is genuinely cube-root elliptic, not elementary.";


End[];
EndPackage[];


(* ===================================================================== *)
(*  Demonstration                                                         *)
(* ===================================================================== *)

(*
   Run with:
     << CubicRootGoursat.wl

   Then:

   (* The cyclic substitution for R = t^3-1 is just t -> omega t. *)
   S = Exp[2 Pi I/3] t;
   R = t^3 - 1;

   (* Test the three monomials F = 1, t, t^2: only F = t fails. *)
   CubicRootElementaryQ[1,   S, t]   (* True  *)
   CubicRootElementaryQ[t,   S, t]   (* False *)
   CubicRootElementaryQ[t^2, S, t]   (* True  *)

   (* Reduce F = 1 *)
   ans = CubicRootReduce[1, R, S, t];
   FullSimplify[D[ans, t] - 1/(t^3 - 1)^(1/3)]    (* should be 0 *)

   (* Reduce F = t^2 *)
   ans = CubicRootReduce[t^2, R, S, t];
   FullSimplify[D[ans, t] - t^2/(t^3 - 1)^(1/3)]  (* should be 0 *)

   (* Reduce F = t -- expected to fail with a specific error message *)
   CubicRootReduce[t, R, S, t]
*)
