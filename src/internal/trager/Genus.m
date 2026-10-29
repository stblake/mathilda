(* ::Package:: *)

(* ========================================================================== *)
(* Genus.m -- the genus of the curve Y^m = G(x), [T] Ch 2 Sec 4 p.29.         *)
(*                                                                            *)
(* Riemann-Hurwitz for the cyclic cover (x, Y) -> x of degree m (tame, the    *)
(* characteristic being 0 or prime to m).  Over a root of an irreducible      *)
(* factor p^e of G there are gcd(m, e) places, each of ramification index     *)
(* m/gcd(m, e); over x = Infinity, the same with e replaced by deg G.  Hence  *)
(*     2 genus - 2 = -2 m + Sum_p deg(p) (m - gcd(m, e_p))                    *)
(*                         + (m - gcd(m, deg G)).                             *)
(* The same formula holds over GF(p), p not dividing m, which is what makes   *)
(* the good-reduction test of TorsionModP.m a comparison of genera.           *)
(* ========================================================================== *)

curveGenus::usage = "curveGenus[m, G, x] is the genus of the curve Y^m = G(x) (G a non-constant polynomial over Q, Y^m - G absolutely irreducible), by Riemann-Hurwitz. curveGenus[m, G, x, p] computes it for the reduction modulo the prime p (p not dividing m), from the factorisation of G over GF(p).";
curveGenus[m_Integer, G_, x_Symbol] := curveGenusFromFactors[m, Select[FactorList[G], !FreeQ[#[[1]], x] &], Exponent[G, x], x];
curveGenus[m_Integer, G_, x_Symbol, p_Integer] := Module[{Gp = gfMod[G, p]},
  curveGenusFromFactors[m, Select[FactorList[Gp, Modulus -> p], !FreeQ[#[[1]], x] &], Exponent[Gp, x], x]];

curveGenusFromFactors::usage = "curveGenusFromFactors[m, pairs, degG, x] evaluates the Riemann-Hurwitz formula from the list {{p_j, e_j}, ...} of irreducible factors of G and deg G.";
curveGenusFromFactors[m_Integer, pairs_List, degG_Integer, x_Symbol] := Module[{ram},
  ram = Total[(Exponent[#[[1]], x] (m - GCD[m, #[[2]]])) & /@ pairs] + (m - GCD[m, degG]);
  (-2 m + ram)/2 + 1];
