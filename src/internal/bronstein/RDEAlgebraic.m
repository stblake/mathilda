(* ::Package:: *)

(* ========================================================================== *)
(* RDEAlgebraic.m -- the Risch differential equation z' + f z = g over the    *)
(* algebraic function field K = Q(x)(Y), Y^m = G(x) ([B] Part 3 pp.160-170,   *)
(* with k = Q and the primitive monomial theta = x; ISSAC'91 for the shape).  *)
(*                                                                            *)
(* Decided by VALID BOUNDS and undetermined coefficients:                     *)
(*                                                                            *)
(*   poles.  Lemma 3.6: z has no pole at a finite place where f and g are     *)
(*     regular (there nu(z') = nu(z) - r would dominate).  The coordinates    *)
(*     on the integral basis w_i = Y^i/d_i are integral exactly when the      *)
(*     element is, so only the irreducible P dividing a coordinate            *)
(*     denominator of f or g carry poles.                                     *)
(*   order at a place p over P (ramification r = m/gcd(m, nu_P(G))).  The    *)
(*     three cases of Prop 3.3 give, together, the single valid bound         *)
(*         alpha_p >= min(0, gamma_p + r, -tau_p)                             *)
(*     (in case (i), beta < -r, the stated bound gamma - beta exceeds gamma + *)
(*     r), where -tau_p is used only when tau_p, the value at p of            *)
(*     r f P/P', is a positive integer.  gamma_p >= r min_i nu_P(g_i) since   *)
(*     the basis is integral.  The values tau_p over all places above P are   *)
(*     among the roots of Res_x(chi(x, z), P) with chi the characteristic     *)
(*     polynomial of r f P/P' over Q(x), its poles at P cleared (the values   *)
(*     at places where the element is regular are the roots of the limit).    *)
(*   order at infinity.  Prop 3.5: alpha >= min(0, gamma - r) or alpha =     *)
(*     r a, a a RATIONAL value of the x^-1 coefficient of f at a place over   *)
(*     infinity.  Those values are the roots of Res_v(w - A(v), v^m - lc G),  *)
(*     A(v) the x^-1 coefficient of the expansion with Y = v x^(d/m)(1 +      *)
(*     ...)^(1/m) on the branch v.                                            *)
(*   degrees.  Z = D z, D = Prod P^(n_P), n_P = ceil(-alpha_P/r), is integral *)
(*     at the finite places, so Z = Sum h_i w_i with h_i in Q[x]; at infinity *)
(*     the leading terms of the h_i w_i carry the distinct powers v^i of the  *)
(*     branch, so no cancellation across all branches: deg h_i <= deg D -     *)
(*     b/r_inf - (i d/m - deg d_i) (Prop 3.7), b the bound at infinity.       *)
(*   solve.  The ansatz in (R2): Z' + (f - D'/D) Z = D g, a linear system    *)
(*     over Q in the coefficients of the h_i.                                 *)
(*                                                                            *)
(* {"NoSolution"} is therefore a proof (valid bounds, an inconsistent exact   *)
(* system); {"Unknown", why} when a candidate polynomial degenerates.         *)
(* ========================================================================== *)

brAlgNuP::usage = "brAlgNuP[c, P, x] is the P-adic valuation of the rational function c in x (P irreducible): its multiplicity in the numerator minus that in the denominator; Infinity for c = 0.";
brAlgNuP[c_, P_, x_Symbol] := Module[{t = Together[c], mult},
  If[t === 0, Return[Infinity]];
  mult[p_] := Module[{k = 0, q = Expand[p]},
    While[Exponent[q, x] >= Exponent[P, x] && PolynomialRemainder[q, P, x] === 0,
      q = Expand[PolynomialQuotient[q, P, x]]; k++];
    k];
  mult[Numerator[t]] - mult[Denominator[t]]];

brAlgDegree::usage = "brAlgDegree[c, x] is the degree at infinity of the rational function c in x (deg numerator - deg denominator); -Infinity for 0.";
brAlgDegree[c_, x_Symbol] := With[{t = Together[c]},
  If[t === 0, -Infinity, Exponent[Numerator[t], x] - Exponent[Denominator[t], x]]];

brAlgPoleFactors::usage = "brAlgPoleFactors[els, x] is the list of the irreducible factors (over Q) of the coordinate denominators of the elements els (coordinate lists): the only finite places where a solution of the Risch DE can have a pole (Lemma 3.6).";
brAlgPoleFactors[els_List, x_Symbol] := Module[{dens},
  dens = Denominator[Together[#]] & /@ Flatten[els];
  DeleteDuplicates[Flatten[First /@ Select[FactorList[#], !FreeQ[#[[1]], x] &] & /@ dens]]];

brAlgRamification::usage = "brAlgRamification[m, G, P, x] is the ramification index of the places of Y^m = G above the irreducible P: m/gcd(m, nu_P(G)).";
brAlgRamification[m_Integer, G_, P_, x_Symbol] := m/GCD[m, brAlgNuP[G, P, x]];

brAlgTauCandidates::usage = "brAlgTauCandidates[B, f, P, r] is the list of the POSITIVE INTEGERS among the values at the places above P of u = r f P/P' (f a coordinate list on the basis B over Q): the integer roots of Res_x(chi, P), chi the characteristic polynomial of u with its poles at P cleared. It returns {\"Unknown\", why} when that resultant vanishes identically.";
brAlgTauCandidates[B_Association, f_List, P_, r_Integer] := Module[
  {x = B["x"], z = Unique["Bronstein`Private`brTz"], u, M, chi, cl, e, poly, L, T, roots},
  u = afScale[fieldQ[], Together[r P/D[P, x]], f];
  M = afMultMatrix[fieldQ[], B, u];
  chi = Together[Det[z IdentityMatrix[B["m"]] - M]];
  cl = CoefficientList[Expand[Numerator[chi]], z]/Denominator[chi];
  (* clear the poles at P: multiply by P^e, e the largest pole order *)
  e = Max[0, Max[-brAlgNuP[#, P, x] & /@ DeleteCases[cl, 0]]];
  cl = Together[# P^e] & /@ cl;
  L = PolynomialLCM @@ (Denominator /@ cl);
  poly = Expand[Total[Table[Together[cl[[k + 1]] L] z^k, {k, 0, Length[cl] - 1}]]];
  T = Expand[Resultant[poly, P, x]];
  If[T === 0, Return[{"Unknown", "the place-value polynomial vanishes identically"}]];
  If[FreeQ[T, z], Return[{}]];
  roots = -Coefficient[#, z, 0]/Coefficient[#, z, 1] & /@
    First /@ Select[FactorList[T], Exponent[#[[1]], z] === 1 &];
  Select[roots, IntegerQ[#] && # > 0 &]];

brAlgInfinityData::usage = "brAlgInfinityData[B] is <|\"r\" -> r, \"q\" -> d r/m, \"lc\" -> lc G|> for the places of Y^m = G above infinity (d = deg G): ramification r = m/gcd(m, d), and Y = v x^(d/m)(1 + ...)^(1/m), v^m = lc G, with x = s^-r and x^(d/m) = s^-q.";
brAlgInfinityData[B_Association] := Module[{m = B["m"], G = B["G"], x = B["x"], d},
  d = Exponent[G, x];
  <|"r" -> m/GCD[m, d], "q" -> (d/GCD[m, d]), "lc" -> Coefficient[G, x, d], "d" -> d|>];

brAlgInfinityCandidates::usage = "brAlgInfinityCandidates[B, f] is the list of the RATIONAL values a of the x^-1 coefficient of the expansion of f (a coordinate list on B) at the places above infinity: the rational roots of Res_v(w - A(v), v^m - lc G), A(v) the coefficient on the branch v (Prop 3.5).";
brAlgInfinityCandidates[B_Association, f_List] := Module[
  {x = B["x"], m = B["m"], G = B["G"], inf = brAlgInfinityData[B], s = Unique["Bronstein`Private`brSs"],
   v = Unique["Bronstein`Private`brVv"], w = Unique["Bronstein`Private`brWw"], r, q, lc, Gt, ex, ser, A, R},
  r = inf["r"]; q = inf["q"]; lc = inf["lc"];
  Gt = Together[(G /. x -> s^-r)/(lc s^(-r inf["d"])) - 1];           (* a polynomial in s^r *)
  ex = Sum[Together[(f[[i + 1]]/B["d"][[i + 1]]) /. x -> s^-r] (v s^-q)^i (1 + Gt)^(i/m), {i, 0, m - 1}];
  ser = Series[ex, {s, 0, r}];
  A = Expand[SeriesCoefficient[ser, r]];
  A = Expand[PolynomialRemainder[A, v^m - lc, v]];
  R = Expand[Resultant[w - A, v^m - lc, v]];
  If[FreeQ[R, w], Return[{}]];
  Select[-Coefficient[#, w, 0]/Coefficient[#, w, 1] & /@
    First /@ Select[FactorList[R], Exponent[#[[1]], w] === 1 &], Element[#, Rationals] === True &]];

brAlgGammaInfinity::usage = "brAlgGammaInfinity[B, g] is a lower bound for the order of g (coordinate list on B) at every place above infinity: min_i nu(g_i w_i) = -r (deg g_i + i d/m - deg d_i); Infinity for g = 0.";
brAlgGammaInfinity[B_Association, g_List] := Module[{inf = brAlgInfinityData[B], x = B["x"], ts},
  ts = Table[If[Together[g[[i + 1]]] === 0, Infinity,
      -inf["r"] (brAlgDegree[g[[i + 1]], x] + i inf["d"]/B["m"] - Exponent[B["d"][[i + 1]], x])], {i, 0, B["m"] - 1}];
  Min[ts]];

brAlgBounds::usage = "brAlgBounds[B, f, g] is <|\"D\" -> D, \"Bdeg\" -> {B_0, ..., B_(m-1)}|> for z' + f z = g over Q(x)(Y) (f, g coordinate lists on B): the denominator D = Prod P^(n_P) and the degree bounds of the numerator coordinates (see the header), or {\"Unknown\", why}.";
brAlgBounds[B_Association, f_List, g_List] := Catch[Module[
  {x = B["x"], m = B["m"], G = B["G"], Ps, Dd = 1, r, gam, tau, amin, nP, inf, rinf, ainf, b, Bdeg},
  Ps = brAlgPoleFactors[{f, g}, x];
  Do[
    r = brAlgRamification[m, G, pf, x];
    gam = r Min[brAlgNuP[#, pf, x] & /@ g];
    tau = brAlgTauCandidates[B, f, pf, r];
    If[MatchQ[tau, {"Unknown", _}], Throw[tau, "brAlgB"]];
    amin = Min[0, If[gam === Infinity, 0, gam + r], If[tau === {}, 0, -Max[tau]]];
    nP = Ceiling[-amin/r];
    Dd = Dd pf^nP,
    {pf, Ps}];
  inf = brAlgInfinityData[B]; rinf = inf["r"];
  ainf = rinf # & /@ brAlgInfinityCandidates[B, f];
  gam = brAlgGammaInfinity[B, g];
  b = Min[Join[{0, If[gam === Infinity, 0, gam - rinf]}, ainf]];
  Bdeg = Table[Floor[Exponent[Dd, x] - b/rinf - (i inf["d"]/m - Exponent[B["d"][[i + 1]], x])], {i, 0, m - 1}];
  <|"D" -> Expand[Dd], "Bdeg" -> Bdeg|>], "brAlgB"];

brAlgRDE::usage = "brAlgRDE[B, f, g] (or brAlgRDE[B, f, g, widen], every bound widened by widen: a test hook) solves z' + f z = g for z in Q(x)(Y) (Y^m = G the curve of the basis B; f, g coordinate lists over Q on B). It returns {\"Solution\", z} (a coordinate list, verified), {\"NoSolution\"} (a proof: valid bounds, an inconsistent system) or {\"Unknown\", why}.";
brAlgRDE[B_Association, f_List, g_List] := brAlgRDE[B, f, g, 0];
brAlgRDE[B_Association, f_List, g_List, widen_Integer] := Module[
  {x = B["x"], m = B["m"], F = fieldQ[], bd, Dd, Bdeg, cs, h, Z, fr, lhs, eqs, sol, z, unk, Pp},
  If[afZeroQ[F, g], Return[{"Solution", afZero[B]}]];
  bd = brAlgBounds[B, f, g];
  If[MatchQ[bd, {"Unknown", _}], Return[bd]];
  Dd = bd["D"]; Bdeg = bd["Bdeg"];
  (* widen > 0: every pole factor to a higher power and every degree bound up
     (a test of the bounds: the solvability must not change) *)
  If[widen > 0,
    Pp = Times @@ Join[{x + 7}, brAlgPoleFactors[{f, g}, x]];
    Dd = Expand[Dd Pp^widen];
    Bdeg = Bdeg + widen (1 + Exponent[Pp, x])];
  cs = Table[If[Bdeg[[i + 1]] < 0, {}, Table[Unique["Bronstein`Private`brCa"], {Bdeg[[i + 1]] + 1}]], {i, 0, m - 1}];
  unk = Flatten[cs];
  If[unk === {}, Return[{"NoSolution"}]];
  h = Table[If[cs[[i]] === {}, 0, cs[[i]] . (x^Range[0, Length[cs[[i]]] - 1])], {i, m}];
  (* (R2): Z' + (f - D'/D) Z - D g = 0, coordinatewise *)
  Z = h;
  fr = afPlus[F, f, afScale[F, -Together[D[Dd, x]/Dd], afOne[B]]];
  lhs = MapThread[Together[#1 + #2 - #3] &, {brAlgDLinear[B, Z], brAlgTimesLinear[B, fr, Z], (Together[Dd #] & /@ g)}];
  eqs = Flatten[CoefficientList[Expand[Numerator[#]], x] & /@ lhs];
  eqs = DeleteCases[Expand /@ eqs, 0];
  sol = If[eqs === {}, {Thread[unk -> 0]}, Quiet[Solve[Thread[eqs == 0], unk]]];
  If[sol === {} || !ListQ[sol], Return[{"NoSolution"}]];
  z = (Together[(# /. First[sol] /. Thread[unk -> 0])/Dd] & /@ h);
  (* verified exactly in K *)
  If[!afZeroQ[F, afMinus[F, afPlus[F, afD[F, B, z], afTimes[F, B, f, z]], g]],
    Return[{"Unknown", "the solution of the linear system does not verify"}]];
  {"Solution", z}];

brAlgDLinear::usage = "brAlgDLinear[B, a] is afD over Q for coordinates that may contain undetermined (linear) coefficients: (c_i' + c_i s_i), Together only (no canonical form over a field, which would treat the unknowns as transcendental but is unnecessary here).";
brAlgDLinear[B_Association, a_List] := Table[Together[D[a[[i]], B["x"]] + a[[i]] B["s"][[i]]], {i, B["m"]}];

brAlgTimesLinear::usage = "brAlgTimesLinear[B, a, b] is afTimes over Q with Together only (for coordinates linear in undetermined coefficients).";
brAlgTimesLinear[B_Association, a_List, b_List] := Module[{m = B["m"], out, kf},
  out = Table[0, {m}];
  Do[
    If[a[[i + 1]] =!= 0 && b[[j + 1]] =!= 0,
      kf = B["mult"][[i + 1, j + 1]];
      out[[kf[[1]] + 1]] += a[[i + 1]] b[[j + 1]] kf[[2]]],
    {i, 0, m - 1}, {j, 0, m - 1}];
  Together /@ out];

brRDEAlgebraic::usage = "brRDEAlgebraic[f, g, x, m, G, Y] solves z' + f z = g over Q(x)(Y), Y^m = G (reduced curve), with f, g expressions in x and the symbol Y (polynomial in Y, rational in x): {\"Solution\", z} (z an expression in x and Y), {\"NoSolution\"} (a proof) or {\"Unknown\", why}. The same interface as brRischDE.";
brRDEAlgebraic[f_, g_, x_Symbol, m_Integer, G_, Y_Symbol] := Module[{B = curveBasis[m, Expand[G], x], F = fieldQ[], fe, ge, r},
  fe = afFromY[F, B, brAlgReduceY[f, Y, m, G], Y];
  ge = afFromY[F, B, brAlgReduceY[g, Y, m, G], Y];
  r = brAlgRDE[B, fe, ge];
  If[First[r] === "Solution", {"Solution", Together[afToY[B, r[[2]], Y]]}, r]];

brAlgReduceY::usage = "brAlgReduceY[e, Y, m, G] is the expression e (rational in x and Y) written as a polynomial in Y of degree < m with coefficients rational in x: the denominator rationalised in Q(x)[Y]/(Y^m - G).";
brAlgReduceY[e_, Y_Symbol, m_Integer, G_] := Module[{t = Together[e], num, den, x, inv},
  num = Numerator[t]; den = Denominator[t];
  If[FreeQ[den, Y], Return[reduceYDegree[Expand[num], Y, m, G]/den]];
  x = First[Select[Variables[G], True &]];
  inv = powerBasisInverse[reduceYDegree[Expand[den], Y, m, G], Y, m, G, x];
  reduceYDegree[Expand[num (inv . (Y^Range[0, m - 1]))], Y, m, G]];
