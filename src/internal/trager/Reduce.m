(* ::Package:: *)

(* ========================================================================== *)
(* Reduce.m -- the curve carrying one component a_i(x) yh^i.                  *)
(*                                                                            *)
(* yh^n = g with g = c Prod p_j^e_j (p_j distinct monic irreducibles over Q). *)
(* The element yh^i satisfies (yh^i)^n = c^i Prod p_j^(i e_j).  Write          *)
(*     i e_j = n q_j + r_j,   0 <= r_j < n,   s = Prod p_j^q_j,               *)
(* so Yt = yh^i / s satisfies Yt^n = c^i Prod p_j^r_j (every exponent < n).   *)
(* Let d = gcd(n, r_1, ..., r_k) ([T] Ch 3 Sec 3, p.39: the binomial Y^n - G   *)
(* is absolutely irreducible iff d = 1, the constants being n-th powers in an *)
(* algebraically closed field).                                               *)
(*                                                                            *)
(*   d = 1: the component lives on the curve Y^n = G, G = c^i Prod p_j^r_j,   *)
(*          with Y = yh^i / s exactly.                                        *)
(*   d > 1: Yt^n = c^i H^d with H = Prod p_j^(r_j/d).  Put m = n/d and         *)
(*          Yhat = H^(1/m) (an explicit radical).  Then (Yt/Yhat)^n is        *)
(*          constant, so rho = Yt/Yhat is LOCALLY CONSTANT, and               *)
(*              Int a Yt dx = rho Int a Yhat dx,                              *)
(*          an exact identity on every branch: the component is integrated    *)
(*          on the absolutely irreducible curve Y^m = H (m = 1 means it is    *)
(*          rational) and multiplied back by rho.  No algebraic constant      *)
(*          c^(1/d) is ever adjoined.                                         *)
(* ========================================================================== *)

componentCurve::usage = "componentCurve[n, g, i, x] describes the curve of the component yh^i, where yh^n = g in Q[x]. It returns <|\"m\" -> m, \"G\" -> G, \"s\" -> s, \"exact\" -> exact|>: the component is s Y (exact = True; Y = yh^i/s satisfies Y^m = G, m = n) or s rho Y (exact = False; Y = G^(1/m) is a new radical and rho = (yh^i/s)/G^(1/m) is locally constant). G is a polynomial over Q whose irreducible factors all have multiplicity < m, and Y^m - G is absolutely irreducible; m = 1 (G = 1) means the component is a rational function times rho.";
componentCurve[n_Integer, g_, i_Integer, x_Symbol] := Module[
  {fl, c, pairs, s = 1, red = {}, d, H, G},
  fl = FactorList[g];
  c = Times @@ (#[[1]]^#[[2]] & /@ Select[fl, FreeQ[#[[1]], x] &]);
  pairs = Select[fl, !FreeQ[#[[1]], x] &];
  Do[
    With[{p = pr[[1]], f = i pr[[2]]},
      s = s p^Quotient[f, n];
      If[Mod[f, n] > 0, AppendTo[red, {p, Mod[f, n]}]]],
    {pr, pairs}];
  d = GCD @@ Join[{n}, Last /@ red];
  If[d === 1,
    G = Expand[c^i Times @@ (#[[1]]^#[[2]] & /@ red)];
    <|"m" -> n, "G" -> G, "s" -> Expand[s], "exact" -> True|>,
    H = Expand[Times @@ (#[[1]]^(#[[2]]/d) & /@ red)];
    <|"m" -> n/d, "G" -> H, "s" -> Expand[s], "exact" -> False|>]];

reducedCurveQ::usage = "reducedCurveQ[m, G, x] is True when Y^m = G is in the form the per-component pipeline requires: m >= 2, G a non-constant polynomial over Q, every irreducible factor of multiplicity < m, and the gcd of m with the multiplicities equal to 1 (absolute irreducibility, [T] p.39).";
reducedCurveQ[m_Integer, G_, x_Symbol] := Module[{pairs},
  If[m < 2 || !rationalPolyQ[G, x] || FreeQ[G, x], Return[False]];
  pairs = Select[FactorList[G], !FreeQ[#[[1]], x] &];
  AllTrue[pairs, #[[2]] < m &] && GCD @@ Join[{m}, Last /@ pairs] === 1];
