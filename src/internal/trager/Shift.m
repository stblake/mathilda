(* ::Package:: *)

(* ========================================================================== *)
(* Shift.m -- behaviour at infinity and the change of variable x = a + 1/z,   *)
(*            [T] Ch 4 Sec 1 p.44, Sec 3 pp.48-50; [B] Sec 2.4.               *)
(*                                                                            *)
(* Hermite reduction is complete ([T] Ch 4 Sec 3, Theorem 1 p.50) when the    *)
(* differential f dx has no pole at infinity and the basis is normal there.   *)
(* The simple-radical basis is always normal at infinity, so a component      *)
(* whose differential is already regular at infinity is integrated in x as it *)
(* stands; otherwise infinity is moved to a finite point by x = a + 1/z, with *)
(* a chosen so that x = a is an ordinary point (G(a) != 0) where the          *)
(* integrand is finite.  In the z-frame the new infinity is x = a, so the      *)
(* differential is regular there.                                             *)
(*                                                                            *)
(* Order at infinity.  At a place P over x = Infinity with ramification index *)
(* r = m/gcd(m, M) (M = deg G):  ord_P(x) = -r, ord_P(Y) = -r M/m,            *)
(* ord_P(dx) = -r - 1.  For c(x) w_j,                                         *)
(*     ord_P(c w_j) = (r/m) (-m deg c + m deg d_j - j M) =: (r/m) lam(c, j),   *)
(* so c w_j dx is regular at infinity iff lam(c, j) >= m + gcd(m, M).        *)
(* ========================================================================== *)

infScaledOrder::usage = "infScaledOrder[F, B, c, j] is lam(c, j) = -m deg c + m deg d_j - j deg G: m/r times the order of c w_j at each place over infinity (r its ramification index). Infinity for c = 0.";
infScaledOrder[F_Association, B_Association, c_, j_Integer] := Module[{dc = rcDeg[F, c, B["x"]]},
  If[dc === -Infinity, Infinity,
    -B["m"] dc + B["m"] Exponent[B["d"][[j + 1]], B["x"]] - j B["deg"]]];

regularAtInfinityQ::usage = "regularAtInfinityQ[F, B, a] is True when the differential a dx (a an element on the basis B) has no pole at any place over x = Infinity.";
regularAtInfinityQ[F_Association, B_Association, a_List] := Module[
  {thr = B["m"] + GCD[B["m"], B["deg"]]},
  AllTrue[Range[0, B["m"] - 1], infScaledOrder[F, B, a[[# + 1]], #] >= thr &]];

integralAtInfinityQ::usage = "integralAtInfinityQ[F, B, a] is True when the function a has no pole at any place over x = Infinity.";
integralAtInfinityQ[F_Association, B_Association, a_List] :=
  AllTrue[Range[0, B["m"] - 1], infScaledOrder[F, B, a[[# + 1]], #] >= 0 &];

shiftCurve::usage = "shiftCurve[c, m, G, x, z] moves infinity to a finite point for the component c(x) Y dx on Y^m = G(x). It picks the first a in 0, 1, -1, 2, -2, ... (no bound: firstOrdinaryInteger) with G(a) != 0 and the denominator of c nonzero at a, substitutes x = a + 1/z and Y = Y1/z^k (k = Ceiling[M/m], M = deg G), and returns <|\"a\" -> a, \"k\" -> k, \"G\" -> Gz, \"c\" -> cz|>: the component is cz(z) Y1 dz on Y1^m = Gz(z) = z^(k m - M) z^M G(a + 1/z), again a reduced, absolutely irreducible curve of the same genus. The inverse substitution is z -> 1/(x - a), Y1 -> Y (x - a)^(-k).";
shiftCurve[c_, m_Integer, G_, x_Symbol, z_Symbol] := Module[
  {den = Denominator[Together[c]], M = Exponent[G, x], a, k, Gz, cz},
  (* G and den have finitely many roots: the search 0, 1, -1, 2, ... ends *)
  a = firstOrdinaryInteger[{G, den}, x];
  k = Ceiling[M/m];
  Gz = Expand[z^(k m - M) Expand[z^M (G /. x -> a + 1/z)]];
  (* dx = -dz/z^2 and Y = Y1 z^(-k) *)
  cz = Together[-(c /. x -> a + 1/z)/z^(k + 2)];
  <|"a" -> a, "k" -> k, "G" -> Gz, "c" -> cz|>];

firstOrdinaryInteger::usage = "firstOrdinaryInteger[polys, x] is the first integer a in 0, 1, -1, 2, -2, ... at which no polynomial of the list polys (each nonzero) vanishes. It exists since each has finitely many roots; no bound is imposed.";
firstOrdinaryInteger[polys_List, x_Symbol] := Module[{i = 0, a},
  If[MemberQ[Expand[polys], 0], throwInternal["firstOrdinaryIntegerZeroPolynomial"]];
  While[True,
    a = If[i === 0, 0, If[OddQ[i], (i + 1)/2, -i/2]];
    If[AllTrue[polys, (# /. x -> a) =!= 0 &], Return[a, Module]];
    i++]];

unshift::usage = "unshift[expr, sh, x, z, Y] maps an expression in z and the generator Y of the shifted curve back to x: z -> 1/(x - a), Y -> Y (x - a)^(-k), simultaneously.";
unshift[expr_, sh_Association, x_Symbol, z_Symbol, Y_Symbol] :=
  expr /. {z -> 1/(x - sh["a"]), Y -> Y (x - sh["a"])^(-sh["k"])};
