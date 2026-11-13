(* ::Package:: *)

(* ========================================================================== *)
(* Curve.m -- the integrand on the curve y^n = g(theta) over k, and the       *)
(*            per-component curves Y^m = G ([T] App A Sec 1 p.74).            *)
(*                                                                            *)
(* brCanonicalForm generalises Trager's canonicalForm (trager/Validate.m)     *)
(* from Q(x) to k(theta): with radicand b = p/q (p, q in Q[x, lower][theta])  *)
(* the symbol yh = q b^(1/n) satisfies yh^n = g = q^(n-1) p exactly, and the  *)
(* integrand, rational in (k, theta, yh), is written f = Sum a_i yh^i with    *)
(* a_i in k(theta) (reduction modulo yh^n - g and inversion of the            *)
(* denominator in k(theta)[yh]/(yh^n - g), Trager's powerBasisInverse).       *)
(*                                                                            *)
(* Components: f is elementary iff each a_i yh^i is, and the algebraic part   *)
(* of the integral of a_i yh^i is v_i yh^i ([T] App A p.74).  Trager's        *)
(* componentCurve works over any UFD in the curve variable: with theta in     *)
(* place of x, the factors of g free of theta are the content c in k, and     *)
(* Y = yh^i/s satisfies Y^m = G with every theta-multiplicity < m.  Two cases *)
(* stay out of scope:                                                         *)
(*   - d = gcd(n, multiplicities) > 1 with a NON-constant content: then       *)
(*     rho = (yh^i/s)/G^(1/m) satisfies rho^n = c^i in k, not a constant, and  *)
(*     the constant field of the component is the radical extension           *)
(*     k(c^(i/d)) -- a coefficient field fieldFn does not represent;          *)
(*   - G free of theta after the reduction (the component's radical is a      *)
(*     radical of an element of k): the same constant-field extension, unless *)
(*     it is a rational constant, when the component is rational in theta.    *)
(* ========================================================================== *)

brSubstituteRadical::usage = "brSubstituteRadical[f, vars, b, n, y] replaces every power b^r (and (1/b)^r) of the radicand b in f by y^(r n), so that the result is rational in vars and y.";
brSubstituteRadical[f_, vars_List, b_, n_Integer, y_Symbol] :=
  f /. Power[c_, r_Rational] /; brDependsQ[c, vars] :>
    Which[Together[c - b] === 0, y^(r n),
          Together[c b - 1] === 0, y^(-r n),
          True, Power[c, r]];

brCanonicalForm::usage = "brCanonicalForm[T] puts the Shape A (or Algebraic/B) integrand of the tower record T in the form f = Sum_{i<n} a_i yh^i over k(theta). It returns <|\"n\" -> n, \"b\" -> b, \"q\" -> q, \"g\" -> g, \"y\" -> yh, \"coeffs\" -> {a_0..a_(n-1)}, \"yexpr\" -> q b^(1/n)|> with yh^n = g = q^(n-1) Numerator[b] a polynomial in the tower symbols, yh a fresh symbol standing for yexpr. Throws a Decline when the integrand is not rational in the tower and the radical, or the radical is reducible.";
brCanonicalForm[T_Association] := Module[{rad = T["radical"], vars, b, n, p, q, g, y, yh, fy, t, num, den, numC, denInv, coeffs},
  If[rad === None, throwInternal["brCanonicalFormNoRadical"]];
  vars = Prepend[T["syms"], T["x"]];
  b = rad["b"]; n = rad["n"];
  p = Numerator[Together[b]]; q = Denominator[Together[b]];
  g = Expand[q^(n - 1) p];
  y = Unique["Bronstein`Private`brY"]; yh = Unique["Bronstein`Private`brYh"];
  fy = brSubstituteRadical[T["f"], vars, b, n, y];
  If[!brRationalInQ[fy, Append[vars, y]],
    throwInternal["brCanonicalFormNotRational"]];
  t = Together[fy /. y -> yh/q];
  num = Numerator[t]; den = Denominator[t];
  numC = powerBasisCoords[num, yh, n, g];
  denInv = powerBasisInverse[den, yh, n, g, T["x"]];
  coeffs = Together /@ powerBasisCoords[(numC . (yh^Range[0, n - 1])) (denInv . (yh^Range[0, n - 1])), yh, n, g];
  <|"n" -> n, "b" -> b, "q" -> q, "g" -> g, "y" -> yh, "coeffs" -> coeffs, "yexpr" -> q b^(1/n)|>];

brComponentCurve::usage = "brComponentCurve[T, n, g, i] is the curve of the component yh^i (yh^n = g, g in Q[x, lower][theta]): Trager's componentCurve with theta as the curve variable, <|\"m\", \"G\", \"s\", \"exact\", \"rho\"|> where the component is s rho Y with Y^m = G, Y = yh^i/s (exact, rho = 1) or Y = G^(1/m) and rho a rational constant times the locally constant (yh^i/s)/G^(1/m). Throws a Decline when the content forces a radical extension of the constant field k (see the header).";
brComponentCurve[T_Association, n_Integer, g_, i_Integer] := Module[{th = T["theta"], cc, fl, c},
  cc = componentCurve[n, g, i, th];
  fl = FactorList[g];
  c = Times @@ (#[[1]]^#[[2]] & /@ Select[fl, FreeQ[#[[1]], th] &]);
  If[!cc["exact"] && brDependsQ[c^i, Prepend[T["lower"], T["x"]]],
    throwUnimplemented["Curve:radical-content", "component over a radical extension of the constant field (non-constant content)"]];
  If[cc["exact"] && FreeQ[cc["G"], th] && brDependsQ[cc["G"], Prepend[T["lower"], T["x"]]],
    throwUnimplemented["Curve:radical-of-k", "component radical of an element of the constant field (a radical below theta)"]];
  cc];

brCurve::usage = "brCurve[T, m, G] is the curve record of Y^m = G over k(theta): <|\"m\", \"G\", \"theta\", \"s\" -> (1/m) DG/G (so D Y = s Y, D the tower derivation), \"content\" -> d, \"primitive\" -> G0|> with G = d G0, d in k, G0 primitive in theta, and \"B\" -> Trager's basis descriptor curveBasis[m, G, theta] (the algebraic integral basis Y^i/d_i over k[theta]).";
brCurve[T_Association, m_Integer, G_] := Module[{th = T["theta"], dc},
  dc = brThetaContent[T, G];
  <|"m" -> m, "G" -> Expand[G], "theta" -> th,
    "s" -> Together[brD[T, G]/(m G)],
    "content" -> dc[[1]], "primitive" -> dc[[2]],
    "B" -> curveBasis[m, Expand[G], th]|>];

brReducedCurveQ::usage = "brReducedCurveQ[T, m, G] is True when Y^m = G is in the form the per-component pipeline requires: m >= 2, G a polynomial in theta over Q[x, lower] of positive degree, every theta-irreducible factor of multiplicity < m, and gcd of m with the multiplicities 1 (absolute irreducibility over the algebraic closure of k, [T] p.39).";
brReducedCurveQ[T_Association, m_Integer, G_] := Module[{th = T["theta"], pairs},
  If[m < 2 || !PolynomialQ[Expand[G], Join[{T["x"]}, T["lower"], {th}]] || FreeQ[G, th], Return[False]];
  pairs = Select[FactorList[G], !FreeQ[#[[1]], th] &];
  AllTrue[pairs, #[[2]] < m &] && GCD @@ Join[{m}, Last /@ pairs] === 1];
