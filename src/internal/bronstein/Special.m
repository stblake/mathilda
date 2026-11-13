(* ::Package:: *)

(* ========================================================================== *)
(* Special.m -- reduction at the special places of an EXPONENTIAL theta       *)
(*              ([B] "The algebraic part" and "Reduction at the special       *)
(*              places" pp.148-152, Prop 2.17-2.18, system (S_exp2)).         *)
(*                                                                            *)
(* f = c Y is REDUCED when it is normal and has no pole above theta = 0 or    *)
(* theta = infinity.  Any solution of (E1) has v = b Y with b in              *)
(* k[theta, 1/theta] at those places (Prop 2.17), and the principal part of   *)
(* f there equals that of v' (logarithmic derivatives have no pole at a       *)
(* special place, Lemma 1.8).  In the one live coordinate, with D theta =     *)
(* eta' theta and s = DG/(mG) regular at 0 and infinity:                      *)
(*   at 0:   the lowest term a theta^(-mu) of c, with -mu + e/m < 0 (a pole of *)
(*           c Y; e = mult of theta in G), must be the lowest term of          *)
(*           D(c1 theta^(-mu)) + s c1 theta^(-mu), i.e.                       *)
(*               c1' + (s(0) - mu eta') c1 = a;                               *)
(*   at oo:  the top term a theta^d of c with d + g/m > 0 likewise needs       *)
(*               c1' + (s(oo) + d eta') c1 = a,  s(oo) = lc'/(m lc) + (g/m) eta'.*)
(* Each is a Risch DE in k; its homogeneous solution lc^(-1/m) theta^-(d+g/m) *)
(* (resp. G0(0)^(-1/m) theta^(mu - e/m)) is not in k-bar (theta is            *)
(* transcendental and the exponent is nonzero), so the solution is unique,    *)
(* and lies in k when it lies in k-bar (Lemma 2.16).  "NoSolution" is         *)
(* therefore a certificate (Prop 2.18), the exponent being allowed to be      *)
(* fractional (ramified places): the leading-coefficient argument is the same *)
(* for Puiseux series.  The reduced integrand is c - (D b + s b), b the sum   *)
(* of the subtracted terms.                                                   *)
(* ========================================================================== *)

brSpecialPoleQ::usage = "brSpecialPoleQ[T, cv, c] is {atZero, atInfinity}: whether c Y has a pole above theta = 0 / theta = infinity (exponential theta).";
brSpecialPoleQ[T_Association, cv_Association, c_] := Module[{z = brPlaceDataZero[T, cv], i = brPlaceDataInf[T, cv]},
  {brSpecialOrderZero[T, c] + z["e"]/cv["m"] < 0, brSpecialDegreeInf[T, c] + i["g"]/cv["m"] > 0}];

brSpecialOrderZero::usage = "brSpecialOrderZero[T, c] is the order of c (in k(theta)) at theta = 0; Infinity for c = 0.";
brSpecialOrderZero[T_Association, c_] := Module[{t = Together[c], th = T["theta"]},
  If[t === 0, Infinity, brThetaOrder[Numerator[t], th] - brThetaOrder[Denominator[t], th]]];

brSpecialDegreeInf::usage = "brSpecialDegreeInf[T, c] is the degree of c (in k(theta)) at theta = infinity; -Infinity for c = 0.";
brSpecialDegreeInf[T_Association, c_] := Module[{t = Together[c], th = T["theta"]},
  If[t === 0, -Infinity, Exponent[Numerator[t], th] - Exponent[Denominator[t], th]]];

brSpecialReduce::usage = "brSpecialReduce[T, cv, c] reduces the component c Y at the special places of an exponential theta. It returns <|\"Algebraic\" -> b, \"Remainder\" -> r|> with c = D b + s b + r and r Y without poles above theta = 0 and theta = infinity; it throws NonElementary when a leading-coefficient Risch DE has no solution ([B] Prop 2.18) and Decline when one is undecided. A logarithmic theta is returned unchanged.";
brSpecialReduce[T_Association, cv_Association, c0_] := Module[
  {th = T["theta"], s = cv["s"], m = cv["m"], eta = brEta[T], TL, z, inf, c = Together[c0], b = 0, mu, d, a, sol, term, it = 0},
  If[T["kind"] =!= "Exp", Return[<|"Algebraic" -> 0, "Remainder" -> c|>]];
  TL = brLowerTower[T];
  z = brPlaceDataZero[T, cv]; inf = brPlaceDataInf[T, cv];
  (* at theta = 0 *)
  While[c =!= 0 && (mu = -brSpecialOrderZero[T, c]) - z["e"]/m > 0,
    If[++it > 500, throwInternal["brSpecialReduceIterations"]];
    a = brLaurentZero[T, c, -mu][[1, 2]];
    sol = brKSolveRDE[TL, Together[brSValueZero[T, cv] - mu eta], a];
    Switch[First[sol],
      "NoSolution", throwNonElementary["no reduction at theta = 0 ([B] Prop 2.18)", <|"mu" -> mu|>],
      "Unknown", throwUnimplemented["Special:rde", "reduction at theta = 0: " <> sol[[2]]]];
    term = Together[sol[[2]] th^(-mu)];
    b = Together[b + term];
    c = Together[c - (brD[T, term] + s term)]];
  (* at theta = infinity *)
  While[c =!= 0 && (d = brSpecialDegreeInf[T, c]) + inf["g"]/m > 0,
    If[++it > 500, throwInternal["brSpecialReduceIterations"]];
    a = brLaurentInf[T, c, d][[1, 2]];
    sol = brKSolveRDE[TL, Together[brSValueInf[T, cv] + d eta], a];
    Switch[First[sol],
      "NoSolution", throwNonElementary["no reduction at theta = infinity ([B] Prop 2.18)", <|"d" -> d|>],
      "Unknown", throwUnimplemented["Special:rde", "reduction at theta = infinity: " <> sol[[2]]]];
    term = Together[sol[[2]] th^d];
    b = Together[b + term];
    c = Together[c - (brD[T, term] + s term)]];
  <|"Algebraic" -> b, "Remainder" -> c|>];

brSValueZero::usage = "brSValueZero[T, cv] is s(0) for s = DG/(mG) and an exponential theta: (e eta' + D(G0(0))/G0(0))/m, G = theta^e G0.";
brSValueZero[T_Association, cv_Association] := Module[{z = brPlaceDataZero[T, cv], TL = brLowerTower[T]},
  Together[(z["e"] brEta[T] + brD[TL, z["G0"]]/z["G0"])/cv["m"]]];

brSValueInf::usage = "brSValueInf[T, cv] is s(infinity) for s = DG/(mG) and an exponential theta: (g eta' + lc'/lc)/m.";
brSValueInf[T_Association, cv_Association] := Module[{i = brPlaceDataInf[T, cv], TL = brLowerTower[T]},
  Together[(i["g"] brEta[T] + brD[TL, i["lc"]]/i["lc"])/cv["m"]]];
