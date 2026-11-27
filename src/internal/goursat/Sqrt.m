(* ========================================================================== *)
(* goursat/Sqrt.m -- p = 1/2: the Klein four-group (Goursat's Theorems 1, 2),  *)
(* his order-3 substitutions, and the second level (multiplication by 2).     *)
(* goursat.py: involution_family, involutions, goursat_F0, conic_uniformizer,  *)
(* _conic_piece, _goursat_piece, quartic_I, order3_reduce, sqrt_reduce,        *)
(* quartic_covariants, sqrt_second_level.  No roots of R are ever computed:    *)
(* the involutions are indexed by the roots theta of the resolvent cubic       *)
(* (quartic R) or of R itself (cubic R).                                       *)
(* ========================================================================== *)

(* The three involutions S(t) = (N t - M)/(L t - N) pairing the four branch      *)
(* points of y^2 = R, as a family over the roots th of a cubic res(th): returns  *)
(* {res, {L, N, M}} with entries in Q[th] (th is the caller's symbol).          *)
involutionFamily[R_, t_, th_] := Module[
  {c, a, s1, s2, s3, s4, res, t0, t1, t2, trip, g},
  c = allCoeffs[R, t];
  a = c[[1]];
  If[Length[c] == 5,
    {s1, s2, s3, s4} = {-c[[2]]/a, c[[3]]/a, -c[[4]]/a, c[[5]]/a};
    res = th^3 - s2 th^2 + (s1 s3 - 4 s4) th - (s1^2 s4 + s3^2 - 4 s2 s4);
    (* with q1 q2 = R/a the pairing, l = s1'-s2', n = p1'-p2', m = p1's2'-p2's1': *)
    (* the three triples below are l (l,n,m), n (l,n,m), m (l,n,m).              *)
    t0 = {s1^2 - 4 s2 + 4 th, s1 th - 2 s3, 2 th (s2 - th) - s1 s3};
    t1 = {s1 th - 2 s3, th^2 - 4 s4, th s3 - 2 s1 s4};
    t2 = {2 th (s2 - th) - s1 s3, th s3 - 2 s1 s4, s3^2 - 4 s4 (s2 - th)};
    Do[trip = Expand[t0 + cc t1 + cc^2 t2];   (* (l + cc n + cc^2 m) != 0 at every root *)
       g = PolynomialGCD[trip[[1]], trip[[2]], trip[[3]], res, Extension -> Automatic];
       If[Exponent[g, th] == 0, Break[]],
       {cc, 0, 49}],
    (* cubic R: the fourth branch point is infinity *)
    {s1, s2} = {-c[[2]]/a, c[[3]]/a};
    res = Expand[(R /. t -> th)/a];
    trip = {1, th, 2 th s1 - 2 th^2 - s2}];
  {Expand[res], trip}];

(* Explicit involutions {L, N, M}, one for each root of the cubic (used only when *)
(* every root is rational or quadratic; CubicField.m handles an irreducible cubic). *)
involutions[R_, t_] := Module[{th, res, trip, out = {}, L, N1, M1},
  {res, trip} = involutionFamily[R, t, th];
  Do[{L, N1, M1} = gK /@ (trip /. th -> r);
     If[!isZero[L],
       {L, N1, M1} = {1, gK[N1/L], gK[M1/L]},          (* (L:N:M) is projective *)
       {L, N1, M1} = {0, 1, gK[M1/N1]}];
     AppendTo[out, {L, N1, M1}],
     {r, rootsOf[res, th]}];
  out];

(* F0 = (F + Sum_theta F o S_theta)/4 as a trace over the cubic, in rational  *)
(* arithmetic: "the sum is a symmetric function of the roots, so one can decide  *)
(* whether it vanishes without knowing them" (Goursat).  The trace of            *)
(* num(theta)/den(theta) is Tr(num(Theta) adj(den(Theta)))/det(den(Theta)),      *)
(* Theta the companion matrix of res.                                          *)
goursatF0[F_, R_, t_] := Module[
  {th, res, L, N1, M1, Fn, Fd, mdeg, Sn, Sd, hom, num, den, r2, r1, r0, Theta, pows,
   atTheta, Mn, Md, trace},
  {res, {L, N1, M1}} = involutionFamily[R, t, th];
  {Fn, Fd} = numDen[F];
  mdeg = Max[Exponent[Fn, t], Exponent[Fd, t]];
  Sn = N1 t - M1; Sd = L t - N1;
  (* F o S as a quotient of polynomials in (t, theta), built without gcds *)
  hom[P_] := Expand[Sum[Coefficient[P, t, i] Sn^i Sd^(mdeg - i), {i, 0, Exponent[P, t]}]];
  num = hom[Fn]; den = hom[Fd];
  {r2, r1, r0} = Rest[allCoeffs[res, th]];
  Theta = {{0, 0, -r0}, {1, 0, -r1}, {0, 1, -r2}};
  pows = {IdentityMatrix[3], Theta, Theta.Theta};
  atTheta[f_] := Module[{cs = CoefficientList[PolynomialRemainder[Expand[f], res, th], th]},
    Sum[Expand[cs[[i]]] pows[[i]], {i, Length[cs]}] + ConstantArray[0, {3, 3}]];
  Mn = atTheta[num];
  Md = Map[Expand, atTheta[den], {2}];
  trace = Expand[Tr[Mn.adj3[Md]]]/Expand[Det[Md]];
  gK[(F + trace)/4]];

(* A real-positive algebraic number? (decides the Euler substitution branch) *)
realPositiveQ[A_] := NumericQ[A] && Module[{v = N[A, 30]}, TrueQ[Im[v] == 0 && Re[v] > 0]];

(* Rational parametrisation of v^2 = Q(x).  Returns {x(m), v(m), kind, data}   *)
(* where m is recovered from (x, v) according to kind:                         *)
(*   "v": m = v;  "root": m = v/(x - x1);  "euler": m = v - sqrt(A) x.          *)
conicUniformizer[Q_, x_, m_] := Module[{cs, A, B, C, disc, sA, x1, x2, xm},
  cs = allCoeffs[Q, x];
  If[Length[cs] == 2,
    {B, C} = cs;
    Return[{(m^2 - C)/B, m, "v", None}, Module]];
  If[Length[cs] != 3, ggInternal["conic: Q not of degree 1 or 2"]];
  {A, B, C} = cs;
  disc = denest[sqrtPull[gK[B^2 - 4 A C]]];
  sA = denest[sqrtPull[A]];
  If[ratQ[disc] || (!ratQ[sA] && !realPositiveQ[A]),
    x1 = gK[(-B + disc)/(2 A)]; x2 = gK[(-B - disc)/(2 A)];
    xm = (x1 m^2 - A x2)/(m^2 - A);
    Return[{xm, gK[m (xm - x1)], "root", x1}, Module]];
  xm = (C - m^2)/(2 sA m - B);
  {xm, gK[m + sA xm], "euler", sA}];

(* Denest a square root of an algebraic number when it is cheap and exact   *)
(* (sympy: sqrtdenest).  Sqrt[a + b Sqrt[c]] = Sqrt[u] + Sqrt[v] when       *)
(* a^2 - b^2 c is a rational square.                                        *)
denest[e_] := Module[{r = e, a, b, c, s, u, v},
  If[MatchQ[e, Power[_Plus, 1/2]],
    With[{arg = Expand[e[[1]]]},
      If[MatchQ[arg, a0_?ratQ + b0_?ratQ Sqrt[c0_?ratQ] | a0_?ratQ + Sqrt[c0_?ratQ]],
        {a, b, c} = Replace[arg,
          {a0_ + b0_ Sqrt[c0_] :> {a0, b0, c0}, a0_ + Sqrt[c0_] :> {a0, 1, c0}}];
        s = Sqrt[a^2 - b^2 c];
        If[ratQ[s],
          u = (a + s)/2; v = (a - s)/2;
          r = Sqrt[u] + Sign[b] Sqrt[v];
          If[!isZero[r^2 - arg], r = e]]]]];
  r];

(* Integral Fm dt/y when it descends along x = X(t) to the conic v^2 = Q(x),    *)
(* v = vfac(t) y:   Fm dt/y = G(x) dx/v.  Returns {antiderivative, log}.        *)
conicPiece[Fm_, X_, vfac_, R_, t_, y_] := Module[
  {x, m, Q, G, xm, vm, kind, data, rho, v, mt},
  Q = descend[vfac^2 R, X, t, x];
  G = descend[Fm vfac/D[X, t], X, t, x];
  {xm, vm, kind, data} = conicUniformizer[Q, x, m];
  rho = gK[(G /. x -> xm) D[xm, m]/vm];
  v = vfac y;
  mt = gK[Switch[kind, "v", v, "root", v/(X - data), "euler", v - data X]];
  {ratInt[rho, m] /. m -> mt,
   <|"X" -> gK[X], "rho" -> rho, "m" -> mt, "kind" -> kind|>}];

(* Integral Fm dt/y for Fm anti-invariant under S = (N t - M)/(L t - N). *)
goursatPiece[Fm_, L_, N1_, M1_, R_, t_, y_] := Module[{X, Dn},
  If[!isZero[L],
    {X, Dn} = {(L t^2 - M1)/(L t - N1), L t - N1},      (* X = t + S(t) *)
    {X, Dn} = {(t - M1/(2 N1))^2, 1}];
  conicPiece[Fm, X, 1/Dn, R, t, y]];

(* The invariant I = 12ae - 3bd + c^2 of R as a binary quartic; I = 0 iff the  *)
(* four branch points are equianharmonic (j = 0).                              *)
quarticI[R_, t_] := Module[{c = allCoeffs[R, t]},
  c = Join[ConstantArray[0, 5 - Length[c]], c];
  Expand[12 c[[1]] c[[5]] - 3 c[[2]] c[[4]] + c[[3]]^2]];

(* Goursat's order-3 substitutions (sections 4-5 of his note), for j = 0.        *)
(* Call with F0 invariant under the Klein four-group.  Take g fixing a rational  *)
(* branch point r; its other fixed point r' is the second root of the Hessian of *)
(* R/(t - r).  With z = (t - r)/(t - r') and X = z^3,                           *)
(*   F0 = phi0(X) + phi1(X) z + phi2(X) z^2,                                    *)
(* and phi1(X) z dt/y is a genus-0 integral.  Returns {anti, residual, log} or   *)
(* None if no rational branch point exists.                                     *)
order3Reduce[F0_, R_, t_, y_] := Module[
  {x, r, rp, c, rat, H, others, z, X, phi, Fel, piece, info},
  If[Exponent[R, t] == 3,
    r = Infinity;
    c = allCoeffs[R, t];
    rp = -c[[2]]/(3 c[[1]]),
    rat = rationalRoots[R, t];
    If[rat === {}, Return[None, Module]];
    r = First[rat];
    H = First[cubicCovariants[gK[R/(t - r)], t]];
    others = Select[rootsOf[H, t], !isZero[# - r] &];
    rp = If[others === {}, Infinity, First[others]]];
  z = lin[r, t]/lin[rp, t];
  X = z^3;
  phi = coords[F0, {1, z, z^2}, X, t, x];
  Fel = gK[(phi[[2]] /. x -> X) z];
  If[Fel === 0, Return[{0, F0, {}}, Module]];
  {piece, info} = conicPiece[Fel, X, z/lin[rp, t]^2, R, t, y];
  {piece, gK[F0 - Fel], {Append[info, "order3FixedPoints" -> {r, rp}]}}];

cubicFieldBudget[] := 10;

(* Goursat's Theorem 2 without the roots of R.  The involutions are peeled off   *)
(* one at a time (greedy), each anti-invariant part integrated on a conic; an    *)
(* irreducible resolvent cubic is handed to CubicField.m.  Returns               *)
(* {antiderivative, residual F0, log}.                                           *)
sqrtReduce[F_, R_, t_, y_] := Module[
  {th, res, trip, invs, rank, total = 0, log = {}, rem = F, L, N1, M1, S, Fm, piece, info},
  {res, trip} = involutionFamily[R, t, th];
  If[AnyTrue[FactorList[res], Exponent[#[[1]], th] == 3 &],
    (* bounded: past cubicFieldBudget[] seconds the elementary part stays an    *)
    (* unevaluated integral, as in goursat.py                                     *)
    Return[TimeConstrained[cubicFieldReduce[F, R, t, y, th, res, trip], cubicFieldBudget[],
      With[{F0 = goursatF0[F, R, t]},
        {ggIntegral[gK[F - F0]/y, t], F0, {<|"note" -> "cubic-field reduction: budget exhausted"|>}}]],
      Module]];
  rank[tr_] := Count[tr, v_ /; !ratQ[v]];
  invs = involutions[R, t];
  invs = Join @@ Table[Select[invs, rank[#] == j &], {j, 0, 3}];   (* stable sort by rank *)
  Do[{L, N1, M1} = inv;
     S = (N1 t - M1)/(L t - N1);
     Fm = gK[(rem - (rem /. t -> S))/2];
     If[Fm =!= 0,
       {piece, info} = goursatPiece[Fm, L, N1, M1, R, t, y];
       AppendTo[log, Append[info, "S" -> gK[S]]];
       total += piece;
       rem = gK[rem - Fm]],
     {inv, invs}];
  {total, rem, log}];

(* Hessian H and sextic covariant T = R' H - R H' of R as a binary quartic (a   *)
(* cubic R is the quartic form with vanishing leading coefficient).             *)
quarticCovariants[R_, t_] := Module[{c = allCoeffs[R, t], a, b, cc, d, e, H, T},
  c = Join[ConstantArray[0, 5 - Length[c]], c];
  {a, b, cc, d, e} = {c[[1]], c[[2]]/4, c[[3]]/6, c[[4]]/4, c[[5]]};
  H = Expand[(a cc - b^2) t^4 + 2 (a d - b cc) t^3 + (a e + 2 b d - 3 cc^2) t^2 +
             2 (b e - cc d) t + (cc e - d^2)];
  T = Expand[D[R, t] H - R D[H, t]];
  {H, T}];

(* The residual of Goursat's Theorem 2 as a Goursat integral.  F0 is invariant   *)
(* under the Klein four-group, hence a function Phi of X = H/R; with Y = T/(R y)  *)
(* one has Y^2 = a3 X^3 + a1 X + a0 (the syzygy) and dt/y = -dX/Y:               *)
(*   F0 dt/y = -Phi(X) dX/Y,                                                    *)
(* a square-root integral on an isomorphic curve.  Returns {sub, back-rules}.    *)
sqrtSecondLevel[F0_, R_, t_, y_, depth_] := Module[
  {x, Xs, Ys, H, T, Phi, a3, a1, a0, sol, R2, sub},
  {H, T} = quarticCovariants[R, t];
  Phi = descend[F0, H/R, t, x];
  sol = solveIdentity[T^2 - a3 H^3 - a1 H R^2 - a0 R^3, {a3, a1, a0}, t];
  R2 = Expand[(a3 /. sol) Xs^3 + (a1 /. sol) Xs + (a0 /. sol)];
  sub = goursatIntegrate[gK[-(Phi /. x -> Xs)], R2, Xs, 1/2, Ys, depth + 1];
  {sub, {Xs -> gK[H/R], Ys -> T/(R y)}}];
