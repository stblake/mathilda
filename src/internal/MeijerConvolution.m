(* MeijerConvolution.m -- definite integrals over (0, Infinity) of
   C x^alpha K1(x) K2(x) by the Meijer G-function convolution theorem.

   Loaded lazily by the definite-integration stage "MeijerConvolution"
   (src/calculus/integrate.c).  All symbols are qualified in MeijerConvolution`.

   KERNELS.  Each kernel is written as  coef * G^{m,n}_{p,q}(w u | a; b)  in a
   common variable u = x^h.  With h = 2 every Bessel function is a G-function
   of x^2, and exp / cos / sin of x become G-functions of x^2 through the
   duplication formula:
     BesselJ[nu, b x] = G^{1,0}_{0,2}(b^2 u/4 | nu/2; -nu/2)
     BesselY[nu, b x] = G^{2,0}_{1,3}(b^2 u/4 | (-nu-1)/2; -nu/2, nu/2, (-nu-1)/2)
     BesselK[nu, b x] = 1/2 G^{2,0}_{0,2}(b^2 u/4 | nu/2, -nu/2)
     E^(-p x)         = Pi^(-1/2) G^{2,0}_{0,2}(p^2 u/4 | 0, 1/2)
     Cos[b x]         = Pi^(1/2) G^{1,0}_{0,2}(b^2 u/4 | 0; 1/2)
     Sin[b x]         = Pi^(1/2) G^{1,0}_{0,2}(b^2 u/4 | 1/2; 0)
     E^(-p x^2)       = G^{1,0}_{0,1}(p u | 0)
     (a^2 + x^2)^-r   = a^(-2 r)/Gamma[r] G^{1,1}_{1,1}(u/a^2 | 1 - r; 0)
   (each verified numerically against Mathematica).

   THEOREM.  With u^(s-1) absorbed into the first factor
   (z^sigma G(z | a; b) = G(z | a + sigma; b + sigma)),
     INT_0^Inf G1(w1 u) G2(w2 u) du = (1/w1) G^{n+mu, m+nu}_{q+r, p+t}( w2/w1 |
         -b_1..-b_m, c_1..c_nu ; -b_{m+1}..-b_q, c_{nu+1}..c_r ||
         -a_1..-a_n, d_1..d_mu ; -a_{n+1}..-a_p, d_{mu+1}..d_t )
   (Prudnikov-Brychkov-Marichev 2.24.1.1).  The resulting MeijerG is reduced to
   named functions by Mathilda's MeijerG evaluator where it can be.

   VALIDITY.  The theorem needs the integral to converge; the stage checks it
   from the kernels' behaviour at 0 (the power x^e0 each contributes) and at
   Infinity (exponential decay, or x^einf with or without oscillation).  A
   condition the assumptions do not prove is returned as a ConditionalExpression;
   a refuted one declines.  No numeric quadrature is used. *)

(* kernel -> {coef, {{an},{ap}}, {{bm},{bq}}, w, e0, einf, decay, freq}
   in u = x^2.  decay: "exp" | "osc" | "alg";  freq: oscillation frequency;
   an optional 9th element is a side condition the kernel form needs. *)
MeijerConvolution`Kernel2[k_, x_] := Module[{nu, b, p, a, r},
  Which[
    MatchQ[k, BesselJ[_, _]] && FreeQ[k[[1]], x] && MeijerConvolution`LinQ[k[[2]], x],
      nu = k[[1]]; b = Coefficient[k[[2]], x, 1];
      {1, {{}, {}}, {{nu/2}, {-nu/2}}, b^2/4, nu, -1/2, "osc", b},
    MatchQ[k, BesselY[_, _]] && FreeQ[k[[1]], x] && MeijerConvolution`LinQ[k[[2]], x],
      nu = k[[1]]; b = Coefficient[k[[2]], x, 1];
      {1, {{}, {(-nu - 1)/2}}, {{-nu/2, nu/2}, {(-nu - 1)/2}}, b^2/4, -Abs[nu], -1/2, "osc", b},
    MatchQ[k, BesselK[_, _]] && FreeQ[k[[1]], x] && MeijerConvolution`LinQ[k[[2]], x],
      nu = k[[1]]; b = Coefficient[k[[2]], x, 1];
      {1/2, {{}, {}}, {{nu/2, -nu/2}, {}}, b^2/4, -Abs[nu], 0, "exp", b},
    MatchQ[k, Power[E, _]] && MeijerConvolution`LinQ[k[[2]], x],
      p = -Coefficient[k[[2]], x, 1];
      {Pi^(-1/2), {{}, {}}, {{0, 1/2}, {}}, p^2/4, 0, 0, "exp", p},
    MatchQ[k, Power[E, _]] && PolynomialQ[k[[2]], x] && Exponent[k[[2]], x] == 2 &&
        Coefficient[k[[2]], x, 1] === 0 && Coefficient[k[[2]], x, 0] === 0,
      p = -Coefficient[k[[2]], x, 2];
      {1, {{}, {}}, {{0}, {}}, p, 0, 0, "exp", p},
    MatchQ[k, Cos[_]] && MeijerConvolution`LinQ[k[[1]], x],
      b = Coefficient[k[[1]], x, 1];
      {Pi^(1/2), {{}, {}}, {{0}, {1/2}}, b^2/4, 0, 0, "osc", b},
    MatchQ[k, Sin[_]] && MeijerConvolution`LinQ[k[[1]], x],
      b = Coefficient[k[[1]], x, 1];
      {Pi^(1/2), {{}, {}}, {{1/2}, {0}}, b^2/4, 1, 0, "osc", b},
    MatchQ[k, Power[_, _]] && FreeQ[k[[2]], x] && PolynomialQ[k[[1]], x] &&
        Exponent[k[[1]], x] == 2 && Coefficient[k[[1]], x, 1] === 0 &&
        Coefficient[k[[1]], x, 2] === 1 && FreeQ[Coefficient[k[[1]], x, 0], x],
      a = Coefficient[k[[1]], x, 0]; r = -k[[2]];
      (* Re[a] > 0 keeps the zeros +-I Sqrt[a] off the real axis: with
         a = -b^2 the integrand has real poles and the theorem does not apply *)
      {a^(-r)/Gamma[r], {{1 - r}, {}}, {{0}, {}}, 1/a, 0, -2 r, "alg", 0, Re[a] > 0},
    True, $Failed]];

MeijerConvolution`LinQ[e_, x_] := PolynomialQ[e, x] && Exponent[e, x] == 1 && Coefficient[e, x, 0] === 0;

(* Split f = C x^alpha K1 K2 (exactly two x-dependent kernels). *)
MeijerConvolution`Split[f_, x_] := Module[{fs, C = 1, alpha = 0, ks = {}},
  fs = If[Head[f] === Times, List @@ f, {f}];
  Do[
    Which[
      FreeQ[g, x], C = C g,
      g === x, alpha += 1,
      MatchQ[g, Power[x, _]] && FreeQ[g[[2]], x], alpha += g[[2]],
      MatchQ[g, Power[E, _]] && Head[g[[2]]] === Plus,
        (* E^(-p x - q x^2) etc. arrive as one Power: split the exponent *)
        ks = Join[ks, (E^# &) /@ (List @@ g[[2]])],
      True, AppendTo[ks, g]],
    {g, fs}];
  {C, alpha, ks}];

(* Prove a conjunction conjunct by conjunct, Refine first: Simplify leaves
   -b^2 > 0 under b > 0 as b^2 < 0 instead of refuting it. *)
MeijerConvolution`Prove[c_, assum_] := Module[{cs, r},
  cs = If[Head[c] === And, List @@ c, {c}];
  r = (With[{t = FixedPoint[Refine[#, assum] &, #, 4]}, If[BooleanQ[t], t, Simplify[t, assum]]] &) /@ cs;
  Which[MemberQ[r, False], False,
        AllTrue[r, TrueQ], True,
        True, And @@ DeleteCases[r, True]]];

MeijerConvolution`Run[f_, x_Symbol, assum_] := Catch[Module[
  {sp, C, alpha, ks, K1, K2, s, w1, w2, c0, cinf, cside, cond, g1a, g1b, g2a, g2b, res, an, ap, bm, bq, cn, cp, dm, dq},
  sp = MeijerConvolution`Split[f, x];
  {C, alpha, ks} = sp;
  If[Length[ks] != 2, Throw[$Failed, MeijerConvolution`fail]];
  K1 = MeijerConvolution`Kernel2[ks[[1]], x];
  K2 = MeijerConvolution`Kernel2[ks[[2]], x];
  If[K1 === $Failed || K2 === $Failed, Throw[$Failed, MeijerConvolution`fail]];
  (* convergence at 0: x^(alpha + e0_1 + e0_2), integrable iff exponent > -1 *)
  c0 = alpha + K1[[5]] + K2[[5]] > -1;
  (* convergence at Infinity *)
  cinf = Which[
    K1[[7]] === "exp" || K2[[7]] === "exp",
      And @@ Cases[{K1, K2}, k_ /; k[[7]] === "exp" :> k[[8]] > 0],
    K1[[7]] === "osc" && K2[[7]] === "osc",
      If[TrueQ[Simplify[K1[[8]]^2 == K2[[8]]^2, assum]],
        alpha + K1[[6]] + K2[[6]] < -1,          (* equal frequency: a non-oscillating part *)
        alpha + K1[[6]] + K2[[6]] < 0 && K1[[8]] > 0 && K2[[8]] > 0],
    K1[[7]] === "osc" || K2[[7]] === "osc",
      alpha + K1[[6]] + K2[[6]] < 0,
    True, alpha + K1[[6]] + K2[[6]] < -1];
  (* side condition of a kernel (element 9, when present) *)
  cside = And @@ Cases[{K1, K2}, k_ /; Length[k] >= 9 :> k[[9]]];
  cond = MeijerConvolution`Prove[c0 && cinf && cside, assum];
  If[cond === False, Throw[$Failed, MeijerConvolution`fail]];
  (* u = x^2: x^alpha dx = u^((alpha+1)/2 - 1) du / 2 ;  s = (alpha+1)/2 *)
  s = (alpha + 1)/2;
  (* Two oscillating kernels of different frequency make the integral
     discontinuous where the frequencies cross (Weber-Schafheitlin), and the
     theorem's G holds only on the side where its argument w2/w1 is below 1.
     Build both orientations and choose by that side; the order of the factors
     in the product must not decide the answer.  The crossing q == 1 itself is
     neither side: there the oscillations cancel to a non-oscillating tail
     x^(alpha+e1+e2), so the integral converges iff that exponent is < -1, and
     then (absolute convergence, uniform near q == 1) it is continuous in q
     with the value of the theorem at argument 1.  Where it diverges the
     convergence condition excludes q == 1.  No branch may fall through to
     Piecewise's implicit 0. *)
  If[K1[[7]] === "osc" && K2[[7]] === "osc" && !TrueQ[Simplify[K1[[8]]^2 == K2[[8]]^2, assum]],
    Module[{q = K2[[4]]/K1[[4]], lt, gt, le, ge, ceq, cq, rb, r12, r21, f12, f21},
      lt = MeijerConvolution`Prove[q < 1, assum];
      gt = MeijerConvolution`Prove[q > 1, assum];
      res = Which[
        lt === True, MeijerConvolution`Conv[K1, K2, C, s, assum],
        gt === True, MeijerConvolution`Conv[K2, K1, C, s, assum],
        True,
          ceq = MeijerConvolution`Prove[alpha + K1[[6]] + K2[[6]] < -1, assum];
          If[ceq =!= True,
            cq = MeijerConvolution`Prove[q != 1 || ceq, assum];
            cond = Which[cond === True, cq, cq === True, cond,
                         True, MeijerConvolution`Prove[cond && cq, assum]];
            If[cond === False, Throw[$Failed, MeijerConvolution`fail]]];
          rb = If[ceq === False, Indeterminate,
                 MeijerConvolution`Conv[K1, ReplacePart[K2, 4 -> K1[[4]]], C, s, assum]];
          le = MeijerConvolution`Prove[q <= 1, assum];
          ge = MeijerConvolution`Prove[q >= 1, assum];
          r12 = If[ge === True, None, MeijerConvolution`Conv[K1, K2, C, s, assum]];
          r21 = If[le === True, None, MeijerConvolution`Conv[K2, K1, C, s, assum]];
          f12 = ceq === True && MeijerConvolution`AtCrossingQ[r12, rb, K1[[8]], K2[[8]], assum];
          f21 = ceq === True && !f12 && MeijerConvolution`AtCrossingQ[r21, rb, K1[[8]], K2[[8]], assum];
          (* Piecewise is HoldAll: hand it evaluated pieces, not locals *)
          Which[
            le === True && ge === True, rb,
            le === True, If[f12, r12, Piecewise @@ {{{r12, Simplify[q < 1, assum]}}, rb}],
            ge === True, If[f21, r21, Piecewise @@ {{{r21, Simplify[q > 1, assum]}}, rb}],
            True, Piecewise @@ {{{r12, Simplify[If[f12, q <= 1, q < 1], assum]},
                                 {r21, Simplify[If[f21, q >= 1, q > 1], assum]}}, rb}]]],
    res = MeijerConvolution`Conv[K1, K2, C, s, assum]];
  If[cond === True, res, ConditionalExpression[res, cond]]],
  MeijerConvolution`fail];

(* Does the one-sided closed form r, continued to the crossing of the
   frequencies f1 == f2, equal the crossing value rb there?  Then r's branch
   owns the boundary and no separate boundary branch is needed. *)
MeijerConvolution`AtCrossingQ[r_, rb_, f1_, f2_, assum_] := Module[{vs, sub, d},
  If[r === None || rb === Indeterminate, Return[False]];
  vs = Variables[f2];
  If[vs === {}, vs = Variables[f1]];
  If[vs === {}, Return[False]];
  sub = Quiet[Solve[f2 == f1, First[vs]]];
  If[!MatchQ[sub, {{_Rule}, ___}], Return[False]];
  d = Quiet[Simplify[(r - rb) /. First[sub], assum]];
  d === 0];

(* The convolution theorem with u^(s-1) absorbed into the first kernel. *)
MeijerConvolution`Conv[K1_, K2_, C_, s_, assum_] := Module[
  {w1 = K1[[4]], w2 = K2[[4]], an, ap, bm, bq, cn, cp, dm, dq},
  an = K1[[2, 1]] + s - 1; ap = K1[[2, 2]] + s - 1;
  bm = K1[[3, 1]] + s - 1; bq = K1[[3, 2]] + s - 1;
  cn = K2[[2, 1]]; cp = K2[[2, 2]]; dm = K2[[3, 1]]; dq = K2[[3, 2]];
  Simplify[(C/2) K1[[1]] K2[[1]] w1^(1 - s) (1/w1) MeijerG[{Join[-bm, cn], Join[-bq, cp]},
        {Join[-an, dm], Join[-ap, dq]}, w2/w1], assum]];
