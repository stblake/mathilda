(* ::Package:: *)

(* === GoursatAppendix.wl ===
   Companion implementation for the paper
       "A Generalisation of Goursat's Algorithm for Integration in Finite Terms"

   Diagnostic for pseudo-elementarity of integrals
       Integrate[F[t]/R[t]^p, t]
   for p in {1/2, 1/3, 2/3}:
     - p = 1/2: Goursat (Theorems 3.3, 3.5)
     - p = 1/3: cube-root analog (Theorem 4.3)
     - p = 2/3: dual cube-root analog (Corollary 4.8)

   Never calls Mathematica's Integrate during the diagnostic phase;
   reductions are constructive.  When the auto-extracting entry point
   IntegrateGoursat[integrand, t] is used, Integrate is invoked
   recursively on the rational reductions and the result is back-
   substituted to produce a closed-form antiderivative in t.

   Entry points:
     IntegrateGoursat[F, R, t, p, u]   -- explicit form (p optional, default 1/2)
     IntegrateGoursat[integrand, t]    -- auto-extracts F, R, p and integrates

   The substitution variable u in the explicit form is the symbol
   appearing in the OUTPUT integrands (G[u]/Sqrt[Q[u]] for p = 1/2,
   J0integrand[u] etc. for p = 1/3, 2/3).  It is independent of any
   internal Module locals and may be named freely by the user.

   Status strings:
     "elementary": the integral is elementary, and explicit rational
                   reductions are returned.
     "obstructed": the eigenspace criterion of the relevant theorem
                   FAILS (a non-zero obstructive eigencomponent is
                   identified and returned).  By Theorem 4.3(ii)/4.3(iii)
                   for p in {1/3, 2/3}, the obstructed verdict implies
                   non-elementarity WHENEVER the descended differential
                   on the obstruction curve is of the second kind, which
                   holds in particular when F is a polynomial in t.  For
                   F with rational poles at non-branch points, the
                   "obstructed" verdict signals failure of Goursat's
                   criterion but does not by itself prove non-elementarity:
                   a complete decision would require the full
                   Risch--Trager--Bronstein algorithm.
*)

(* ========================================================== *)
(* Helper: canonicalise across algebraic-number extensions     *)
(* ========================================================== *)

canonic[e_] :=
  Cancel[Together[e, Extension -> Automatic], Extension -> Automatic];

(* ========================================================== *)
(* Square-root case (Theorems 3.3 and 3.5)                     *)
(* ========================================================== *)

(* The Mobius involution swapping pairs {a,b} and {c,d}.  Handles  *)
(* Infinity in any slot by reducing to the canonical case d -> oo, *)
(* in which the limit of the general formula simplifies to:        *)
(*    S(t) = (c t + a b - c (a+b)) / (t - c).                      *)
MobiusInvolution[{a_, b_}, {c_, d_}, t_] :=
  Which[
    a === Infinity, MobiusInvolution[{c, d}, {b, a}, t],
    b === Infinity, MobiusInvolution[{c, d}, {a, b}, t],
    c === Infinity, MobiusInvolution[{a, b}, {d, c}, t],
    d === Infinity,
      canonic[(c t + a b - c (a + b))/(t - c)],
    True,
      canonic[((a b - c d) t + (a + b) c d - (c + d) a b)/
                 (((a + b) - (c + d)) t - (a b - c d))]
  ];

V4Projections[F_, {S1_, S2_, S3_}, t_] :=
  Module[{f, f1, f2, f3},
    {f, f1, f2, f3} =
      {F, F /. t -> S1, F /. t -> S2, F /. t -> S3};
    canonic /@
      {(f + f1 + f2 + f3)/4,
       (f + f1 - f2 - f3)/4,
       (f - f1 + f2 - f3)/4,
       (f - f1 - f2 + f3)/4}];

ToFunctionOfSquare[H_, u_, x_] :=
  If[PossibleZeroQ[Together[H]], 0,
    Module[{num, den, dN, dD},
      {num, den} = Through[{Numerator, Denominator}[canonic[H]]];
      dN = Exponent[num, u]; dD = Exponent[den, u];
      Sum[Coefficient[num, u, 2 k] x^k, {k, 0, dN/2}]/
      Sum[Coefficient[den, u, 2 k] x^k, {k, 0, dD/2}]]];

(* GoursatReduction: takes user-provided uOut as the output variable. *)
(* The intermediate Mobius substitution variable is a fresh local.    *)
GoursatReduction[Fj_, R_, t_, S_, uOut_] :=
  Module[{u, fps, alpha, beta, tu, Rfact, lc, Q, Fu, gu, gx, pre, back},
    fps = t /. Solve[S == t, t];
    {alpha, beta} = If[Length[fps] == 1, {fps[[1]], Infinity}, fps];
    tu = If[beta === Infinity, alpha + u, (alpha - beta u)/(1 - u)];
    Rfact = If[beta === Infinity,
               Expand[R /. t -> tu],
               Expand[canonic[(R /. t -> tu) (1 - u)^4]]];
    lc = Coefficient[Rfact, u, 4];
    Q = Sum[Coefficient[Rfact, u, 2 k]/lc uOut^k, {k, 0, 2}];
    Fu = canonic[Fj /. t -> tu];
    gu = canonic[Fu/u];
    gx = ToFunctionOfSquare[gu, u, uOut];
    pre = If[beta === Infinity, 1, alpha - beta]/(2 Sqrt[lc]);
    back = If[beta === Infinity,
              uOut -> (t - alpha)^2,
              uOut -> ((t - alpha)/(t - beta))^2];
    <|"S" -> S, "alpha" -> alpha, "beta" -> beta,
      "prefactor" -> pre, "G" -> gx, "Q" -> Q,
      "back" -> back|>];

PickAntiInvolution[j_, S_, t_] :=
  Module[{cands, withInf},
    cands = Complement[{1, 2, 3}, {j}];
    withInf = Select[cands,
      Length[t /. Solve[S[[#]] == t, t]] == 1 &];
    S[[If[Length[withInf] > 0, First[withInf], First[cands]]]]];

GoursatTest[F_, R_, t_, uOut_] :=
  Module[{roots, S, p, reds},
    roots = DeleteDuplicates[t /. Solve[R == 0, t]];
    (* The cubic case is identical with r_4 = Infinity. *)
    If[Length[roots] == 3, AppendTo[roots, Infinity]];
    If[Length[roots] != 4, Return[$Failed]];
    S = With[{r = roots},
      {MobiusInvolution[r[[{1, 2}]], r[[{3, 4}]], t],
       MobiusInvolution[r[[{1, 3}]], r[[{2, 4}]], t],
       MobiusInvolution[r[[{1, 4}]], r[[{2, 3}]], t]}];
    p = V4Projections[F, S, t];
    If[!PossibleZeroQ[p[[1]]],
      Return[<|"status" -> "obstructed",
              "involutions" -> S, "obstruction" -> p[[1]]|>]];
    reds = Association @ Table[
      If[!PossibleZeroQ[p[[j + 1]]],
         j -> GoursatReduction[p[[j + 1]], R, t,
                PickAntiInvolution[j, S, t], uOut],
         Nothing], {j, 1, 3}];
    <|"status" -> "elementary", "involutions" -> S,
      "F0" -> p[[1]], "F1" -> p[[2]],
      "F2" -> p[[3]], "F3" -> p[[4]],
      "reductions" -> reds|>];

(* ========================================================== *)
(* Cube-root case at exponent 1/3 (Theorem 4.3)                *)
(* ========================================================== *)

CyclicMobius[roots_List, t_] :=
  Module[{r, pos, A, B, C, sol},
    r = roots;
    pos = FirstPosition[r, Infinity, {0}, 1][[1]];
    If[pos > 0,
      r = RotateLeft[r, pos];
      Together[(r[[1]] t -
        (r[[1]]^2 - r[[1]] r[[2]] + r[[2]]^2))/(t - r[[2]])],
      sol = First[Solve[
        {(A r[[1]] + B)/(C r[[1]] + 1) == r[[2]],
         (A r[[2]] + B)/(C r[[2]] + 1) == r[[3]],
         (A r[[3]] + B)/(C r[[3]] + 1) == r[[1]]}, {A, B, C}]];
      canonic[(A t + B)/(C t + 1) /. sol]]];

EigenProjection[H_, z_, k_Integer] :=
  With[{w = Exp[2 Pi I/3]},
    canonic[Simplify[(H + w^(-k) (H /. z -> w z) +
      w^(-2 k) (H /. z -> w^2 z))/3]]];

ToFunctionOfCube[H_, z_, x_] :=
  If[PossibleZeroQ[Together[H]], 0,
    Module[{num, den, dN, dD},
      {num, den} = Through[{Numerator, Denominator}[canonic[H]]];
      dN = Exponent[num, z]; dD = Exponent[den, z];
      Sum[Coefficient[num, z, 3 k] x^k, {k, 0, dN/3}]/
      Sum[Coefficient[den, z, 3 k] x^k, {k, 0, dD/3}]]];

(* CubicTest: takes uOut for output variable.  Both J0integrand and *)
(* J2integrand use uOut, and the result includes "J0back" / "J2back" *)
(* substitution rules expressing uOut in terms of t.                  *)
CubicTest[F_, R_, t_, uOut_] :=
  Module[{roots, S, fps, alpha, beta, z, x, tz, Rz, cval, K,
          H, H0, H1, H2, phi0, psi2, J0, J2, J0back, J2back},
    roots = DeleteDuplicates[t /. Solve[R == 0, t]];
    If[Length[roots] == 2, AppendTo[roots, Infinity]];
    If[Length[roots] != 3, Return[$Failed]];
    S = CyclicMobius[roots, t];
    fps = t /. Solve[S == t, t];
    {alpha, beta} = If[Length[fps] == 1, {fps[[1]], Infinity}, fps];
    tz = If[beta === Infinity, alpha + z, (alpha - beta z)/(1 - z)];
    Rz = If[beta === Infinity,
            canonic[R /. t -> tz],
            canonic[(R /. t -> tz) (1 - z)^3]];
    cval = Coefficient[Expand[Rz], z, 3];
    K = -Coefficient[Expand[Rz], z, 0]/cval;
    H = Together[
      If[beta === Infinity,
         (F /. t -> tz)/cval^(1/3),
         (alpha - beta) (F /. t -> tz)/(cval^(1/3) (1 - z))]];
    {H0, H1, H2} = EigenProjection[H, z, #] & /@ {0, 1, 2};
    If[!PossibleZeroQ[Simplify[H1]],
      Return[<|"status" -> "obstructed", "S" -> S,
              "alpha" -> alpha, "beta" -> beta, "K" -> K,
              "obstruction" -> H1|>]];
    phi0 = ToFunctionOfCube[H0, z, x];
    psi2 = ToFunctionOfCube[Together[H2/z^2], z, x];
    J0 = Together[(phi0 /. x -> K/(1 - uOut^3)) uOut/(1 - uOut^3)];
    J2 = Together[(psi2 /. x -> uOut^3 + K) uOut];
    (* Back-substitutions u_J0 and u_J2 in terms of t,                *)
    (* derived from the canonical form R(t)(1-z)^3 = c(z^3-K).        *)
    J0back = If[beta === Infinity,
                uOut -> R^(1/3)/(cval^(1/3) (t - alpha)),
                uOut -> R^(1/3) (alpha - beta)/(cval^(1/3) (t - alpha))];
    J2back = If[beta === Infinity,
                uOut -> R^(1/3)/cval^(1/3),
                uOut -> R^(1/3) (alpha - beta)/(cval^(1/3) (t - beta))];
    <|"status" -> "elementary", "S" -> S,
      "alpha" -> alpha, "beta" -> beta, "K" -> K, "c" -> cval,
      "H0" -> H0, "H2" -> H2,
      "J0integrand" -> J0, "J2integrand" -> J2,
      "J0back" -> J0back, "J2back" -> J2back|>];

(* ========================================================== *)
(* Cube-root case at exponent 2/3 (Corollary 4.8)              *)
(* ========================================================== *)

CubicTest23[F_, R_, t_, uOut_] :=
  Module[{roots, S, fps, alpha, beta, z, x, tz, Rz, cval, K,
          Htilde, H0, H1, H2, phi1, phi2, J1, J2, J1back, J2back},
    roots = DeleteDuplicates[t /. Solve[R == 0, t]];
    If[Length[roots] == 2, AppendTo[roots, Infinity]];
    If[Length[roots] != 3, Return[$Failed]];
    S = CyclicMobius[roots, t];
    fps = t /. Solve[S == t, t];
    {alpha, beta} = If[Length[fps] == 1, {fps[[1]], Infinity}, fps];
    tz = If[beta === Infinity, alpha + z, (alpha - beta z)/(1 - z)];
    Rz = If[beta === Infinity,
            canonic[R /. t -> tz],
            canonic[(R /. t -> tz) (1 - z)^3]];
    cval = Coefficient[Expand[Rz], z, 3];
    K = -Coefficient[Expand[Rz], z, 0]/cval;
    Htilde = Together[
      If[beta === Infinity,
         (F /. t -> tz)/cval^(2/3),
         (alpha - beta) (F /. t -> tz)/cval^(2/3)]];
    {H0, H1, H2} = EigenProjection[Htilde, z, #] & /@ {0, 1, 2};
    If[!PossibleZeroQ[Simplify[H0]],
      Return[<|"status" -> "obstructed", "S" -> S,
              "alpha" -> alpha, "beta" -> beta, "K" -> K,
              "obstruction" -> H0|>]];
    phi1 = ToFunctionOfCube[Together[H1/z], z, x];
    phi2 = ToFunctionOfCube[Together[H2/z^2], z, x];
    J1 = Together[-(phi1 /. x -> K uOut^3/(uOut^3 - 1)) uOut/(uOut^3 - 1)];
    J2 = Together[phi2 /. x -> uOut^3 + K];
    J1back = If[beta === Infinity,
                uOut -> cval^(1/3) (t - alpha)/R^(1/3),
                uOut -> cval^(1/3) (t - alpha)/(R^(1/3) (alpha - beta))];
    J2back = If[beta === Infinity,
                uOut -> R^(1/3)/cval^(1/3),
                uOut -> R^(1/3) (alpha - beta)/(cval^(1/3) (t - beta))];
    <|"status" -> "elementary", "S" -> S,
      "alpha" -> alpha, "beta" -> beta, "K" -> K, "c" -> cval,
      "H1" -> H1, "H2" -> H2,
      "J1integrand" -> J1, "J2integrand" -> J2,
      "J1back" -> J1back, "J2back" -> J2back|>];

(* ========================================================== *)
(* Unified entry point with substitution variable             *)
(* ========================================================== *)

IntegrateGoursat[F_, R_, t_, p_ : 1/2, uOut_] :=
  Switch[p,
    1/2, GoursatTest[F, R, t, uOut],
    1/3, CubicTest[F, R, t, uOut],
    2/3, CubicTest23[F, R, t, uOut],
    _, "Only p = 1/2, 1/3, 2/3 are supported"];

(* ========================================================== *)
(* Auto-extracting entry point: parses F, R, p from the        *)
(* integrand and computes a closed-form antiderivative in t    *)
(* via a recursive call to Integrate on the rational pieces.   *)
(* ========================================================== *)

(* Parse expression of the form F[t]/R[t]^p, returning {F, R, p} *)
(* or $Failed if the form is not recognized.                       *)
ParseGoursatIntegrand[expr_, t_] :=
  Replace[expr, {
    f_. * Power[r_, q_]
      /; PolynomialQ[r, t] && MemberQ[{-1/2, -1/3, -2/3}, q] :>
        {canonic[f], r, -q},
    Power[r_, q_]
      /; PolynomialQ[r, t] && MemberQ[{-1/2, -1/3, -2/3}, q] :>
        {1, r, -q},
    _ :> $Failed
  }];

(* Auto-integrate p = 1/2 case: sum prefactor * Integrate[G/Sqrt[Q]] *)
ComputeAntiderivative12[result_, t_, uOut_] :=
  Module[{reductions, total = 0},
    reductions = result["reductions"];
    Do[
      With[{red = reductions[j]},
        total = total + red["prefactor"] *
          (Integrate[red["G"]/Sqrt[red["Q"]], uOut] /. red["back"])
      ],
      {j, Keys[reductions]}];
    total];

(* Auto-integrate p = 1/3 case: J0 in uOut + J2 in uOut, *)
(* each back-substituted to t.                            *)
ComputeAntiderivative13[result_, t_, uOut_] :=
  Module[{total = 0, J0, J2},
    J0 = result["J0integrand"];
    J2 = result["J2integrand"];
    If[!PossibleZeroQ[J0],
      total = total + (Integrate[J0, uOut] /. result["J0back"])];
    If[!PossibleZeroQ[J2],
      total = total + (Integrate[J2, uOut] /. result["J2back"])];
    total];

(* Auto-integrate p = 2/3 case: J1 in uOut + J2 in uOut. *)
ComputeAntiderivative23[result_, t_, uOut_] :=
  Module[{total = 0, J1, J2},
    J1 = result["J1integrand"];
    J2 = result["J2integrand"];
    If[!PossibleZeroQ[J1],
      total = total + (Integrate[J1, uOut] /. result["J1back"])];
    If[!PossibleZeroQ[J2],
      total = total + (Integrate[J2, uOut] /. result["J2back"])];
    total];

(* The 2-argument auto-extracting form *)
IntegrateGoursat[integrand_, t_] :=
  Module[{parsed, F, R, p, uOut, result},
    parsed = ParseGoursatIntegrand[integrand, t];
    If[parsed === $Failed, Return[$Failed]];
    {F, R, p} = parsed;
    uOut = Unique["u"];
    result = IntegrateGoursat[F, R, t, p, uOut];
    If[!AssociationQ[result] || result["status"] =!= "elementary",
      Return[result]];
    Switch[p,
      1/2, ComputeAntiderivative12[result, t, uOut],
      1/3, ComputeAntiderivative13[result, t, uOut],
      2/3, ComputeAntiderivative23[result, t, uOut]
    ]];


(* ========================================================== *)
(* Sample usage                                                *)
(* ========================================================== *)

(*
  Explicit form returns the diagnostic Association:

    IntegrateGoursat[t, (t^2-1)(t^2-4), t, 1/2, u]
      -> Example 5.1: status "elementary".  The "reductions" field
         has J=2,3 entries with G[u], Q[u], prefactor, and the
         back-substitution u -> t^2.

    IntegrateGoursat[1, t^3 - 1, t, 1/3, u]
      -> Example 5.2: status "elementary" with K=1, H_0=1, H_2=0,
         J0integrand = u/(1 - u^3),  J2integrand = 0,
         J0back = (u -> (t^3-1)^(1/3)/t).

    IntegrateGoursat[t, t^3 + 1, t, 1/2, u]
      -> The cubic-radicand sqrt case (Infinity treated as the 4th
         ramification point).  Returns status "obstructed" with
         obstruction F^(0) = (t^4 - 8 t)/(4 (t^3 + 1)).

  Auto-extracting form: parses F, R, p from the integrand and
  computes a closed-form antiderivative in t by recursively
  invoking Integrate.

    IntegrateGoursat[t/Sqrt[(t^2-1)(t^2-4)], t]
      -> elementary closed form (logs and arctangents/arctanhs).

    IntegrateGoursat[1/(t^3 - 1)^(1/3), t]
      -> the closed form of Example 5.2.

    IntegrateGoursat[t^2/(t^3 - 1)^(1/3), t]
      -> (t^3 - 1)^(2/3)/2 + C.

    IntegrateGoursat[t^2/(t^3 - 1)^(2/3), t]
      -> (t^3 - 1)^(1/3) + C.

    IntegrateGoursat[t/Sqrt[t^3 + 1], t]
      -> the diagnostic Association reporting non-elementarity
         (the integral is elliptic).
*)
