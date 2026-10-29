(* ::Package:: *)

(* ========================================================================== *)
(* Validate.m -- recognise the radical and put the integrand in the canonical *)
(*               form f = Sum_{i=0}^{n-1} a_i(x) yh^i.                        *)
(*                                                                            *)
(* Input: an expression in x built from rational operations, rational numbers *)
(* and powers b^r (r a non-integer rational) of ONE radicand b in Q(x).       *)
(* Different powers of the same radicand are merged: b^(1/2) and b^(1/3) are  *)
(* y^3 and y^2 for y = b^(1/6).                                               *)
(*                                                                            *)
(* The radicand is then cleared of its denominator.  With b = p/q (p, q in    *)
(* Q[x]) the algebraic function yh = q y satisfies                            *)
(*        yh^n = q^n b = q^(n-1) p   in Q[x],                                 *)
(* an exact algebraic identity (no branch is chosen: yh is defined as q y).   *)
(* The integrand, rational in (x, y), is rational in (x, yh), and reducing     *)
(* modulo yh^n - g and inverting its denominator in the field                 *)
(* Q(x)[yh]/(yh^n - g) gives the unique coordinates a_i in the power basis.   *)
(* ========================================================================== *)

radicalPowers::usage = "radicalPowers[f, x] lists the distinct subexpressions b^r of f with r a non-integer rational and b depending on x.";
radicalPowers[f_, x_Symbol] := DeleteDuplicates[
  Cases[f, Power[b_, r_Rational] /; !FreeQ[b, x], {0, Infinity}]];

detectRadical::usage = "detectRadical[f, x] finds the single radical of f. It returns <|\"b\" -> b, \"n\" -> n|>, b in Q(x) and n >= 2, such that every radical subexpression of f is a power of b^(1/n); None when f has no radical; or throws a Decline when f is outside the scope (two different radicands, nested radicals, a transcendental function, a non-rational constant).";
detectRadical[f_, x_Symbol] := Module[{pows, bases, b0, exps, n},
  If[!FreeQ[f, (Sin | Cos | Tan | Cot | Sec | Csc | Exp | Log | ArcTan | ArcSin |
        ArcCos | Sinh | Cosh | Tanh | Erf | Gamma | Power[E, _])],
    throwDecline["transcendental function in the integrand"]];
  If[!FreeQ[f, Complex] || !FreeQ[f, _Real],
    throwDecline["coefficients outside Q"]];
  If[!FreeQ[f, Power[c_, _Rational] /; FreeQ[c, x]],
    throwDecline["algebraic constant in the integrand (coefficients must be rational)"]];
  pows = radicalPowers[f, x];
  If[pows === {}, Return[None]];
  If[!FreeQ[First /@ pows, Power[_, _Rational]],
    throwDecline["nested radicals"]];
  (* Normalise each base; a base that is the reciprocal of another base is    *)
  (* the same radicand with the exponent negated.                             *)
  bases = Together[First[#]] & /@ pows;
  b0 = First[bases];
  If[!rationalFunctionQ[b0, x], throwDecline["radicand is not a rational function over Q"]];
  exps = MapThread[
    Which[
      Together[#1 - b0] === 0, #2,
      Together[#1 b0 - 1] === 0, -#2,
      True, throwDecline["more than one radicand (a radical compositum)",
              <|"Radicands" -> DeleteDuplicates[bases]|>]] &,
    {bases, Last /@ pows}];
  n = LCM @@ (Denominator /@ exps);
  <|"b" -> b0, "n" -> n|>];

substituteRadical::usage = "substituteRadical[f, x, b, n, y] replaces every power b^r (and (1/b)^r) of the radicand in f by the corresponding power y^(r n) of the symbol y, so that the result is rational in x and y.";
substituteRadical[f_, x_Symbol, b_, n_Integer, y_Symbol] :=
  f /. Power[c_, r_Rational] /; !FreeQ[c, x] :>
    Which[Together[c - b] === 0, y^(r n),
          Together[c b - 1] === 0, y^(-r n),
          True, Power[c, r]];

reduceYDegree::usage = "reduceYDegree[p, y, n, g] rewrites the polynomial p in y modulo y^n - g, so that every power of y is below n.";
reduceYDegree[p_, y_Symbol, n_Integer, g_] := Module[{cl = CoefficientList[Expand[p], y], out},
  If[Length[cl] <= n, Return[Expand[p]]];
  out = Table[0, {n}];
  Do[out[[Mod[k, n] + 1]] += cl[[k + 1]] g^Quotient[k, n], {k, 0, Length[cl] - 1}];
  Expand[out . (y^Range[0, n - 1])]];

powerBasisCoords::usage = "powerBasisCoords[p, y, n, g] is the coordinate list {c_0, ..., c_(n-1)} of the polynomial p in y reduced modulo y^n - g.";
powerBasisCoords[p_, y_Symbol, n_Integer, g_] :=
  PadRight[CoefficientList[reduceYDegree[p, y, n, g], y], n];

powerBasisInverse::usage = "powerBasisInverse[q, y, n, g, x] is the coordinate list of 1/q in the power basis 1, y, ..., y^(n-1) of Q(x)[y]/(y^n - g), q a nonzero polynomial in y with coefficients in Q(x). It solves the linear system c.M = e_0 with M the multiplication-by-q matrix. Throws a Decline when q is a zero divisor (y^n - g reducible over Q(x)).";
powerBasisInverse[q_, y_Symbol, n_Integer, g_, x_Symbol] := Module[{M, sol},
  M = Table[powerBasisCoords[y^k q, y, n, g], {k, 0, n - 1}];
  If[Together[Det[M]] === 0,
    throwDecline["the denominator is a zero divisor: the radical is reducible"]];
  sol = LinearSolve[Transpose[M], UnitVector[n, 1]];
  Together /@ sol];

canonicalForm::usage = "canonicalForm[f, x] validates the integrand f and returns <|\"x\" -> x, \"b\" -> b, \"n\" -> n, \"q\" -> q, \"g\" -> g, \"y\" -> yh, \"coeffs\" -> {a_0, ..., a_(n-1)}|> with f = Sum a_i yh^i, yh = q b^(1/n) (q the denominator of b), yh^n = g = q^(n-1) Numerator[b] in Q[x], and yh a fresh symbol standing for q b^(1/n). For an integrand with no radical it returns n = 1 and coeffs = {f}. Throws a Decline outside the scope.";
canonicalForm[f_, x_Symbol] := Module[{rad, b, n, p, q, g, y, yh, fy, t, num, den, numC, denInv, coeffs},
  rad = detectRadical[f, x];
  If[rad === None,
    If[!rationalFunctionQ[f, x], throwDecline["integrand is not algebraic over Q(x)"]];
    Return[<|"x" -> x, "b" -> 1, "n" -> 1, "q" -> 1, "g" -> 1, "y" -> None, "coeffs" -> {Together[f]}|>]];
  b = rad["b"]; n = rad["n"];
  p = Numerator[Together[b]]; q = Denominator[Together[b]];
  g = Expand[q^(n - 1) p];
  y = Unique["y"]; yh = Unique["Y"];
  fy = substituteRadical[f, x, b, n, y];
  If[!rationalInXYQ[fy, x, y],
    throwDecline["integrand is not a rational function of x and the radical over Q"]];
  (* y = yh/q *)
  t = Together[fy /. y -> yh/q];
  num = Numerator[t]; den = Denominator[t];
  numC = powerBasisCoords[num, yh, n, g];
  denInv = powerBasisInverse[den, yh, n, g, x];
  coeffs = Together /@ powerBasisCoords[(numC . (yh^Range[0, n - 1])) (denInv . (yh^Range[0, n - 1])), yh, n, g];
  <|"x" -> x, "b" -> b, "n" -> n, "q" -> q, "g" -> g, "y" -> yh, "coeffs" -> coeffs|>];

radicalExpression::usage = "radicalExpression[cf] is the explicit expression for the symbol cf[\"y\"] of a canonical form: q b^(1/n).";
radicalExpression[cf_Association] := cf["q"] cf["b"]^(1/cf["n"]);
