(* ::Package:: *)

(* ========================================================================== *)
(* K0.m -- arithmetic in k0 = k(Y), Y^m = G(x), and in k0[theta]: Shape B's   *)
(*         algebraic base field (BRONSTEIN_PLAN.md "Shape B").                *)
(*                                                                            *)
(* k = Q(x)(t_1, ..., t_r) is the transcendental tower, t_r = theta at the    *)
(* top; Y is algebraic over Q(x), so k is NOT algebraically closed in k(Y)   *)
(* and Shape A's curve machinery (a curve over k(theta)) does not apply:      *)
(* here the log part is Rothstein-Trager over the FIELD k0, and its gcds need *)
(* arithmetic modulo Y^m - G.                                                 *)
(*                                                                            *)
(* An element is Trager's coordinate list {c_0, ..., c_(m-1)} on the integral *)
(* basis w_i = Y^i/d_i (curveBasis[m, G, x]); the coordinates are rational    *)
(* functions in x over the coefficient field F = fieldFn[{t_1..t_r}] (the     *)
(* tower symbols, theta included, are SCALARS of F: the curve variable is x), *)
(* with an algebraic extension Q(alpha) for a residue field.  Trager's af*    *)
(* arithmetic is field-agnostic, so this module is mostly naming plus the     *)
(* tower derivation and the Euclidean algorithm in theta over k0.             *)
(*                                                                            *)
(* Soundness of the zero test: {w_i} is a basis of k0 over k because Y^m - G *)
(* is irreducible over Q(x) (componentCurve / reducedCurveQ) and stays        *)
(* irreducible over k = Q(x)(t), the t_i being transcendental over Q(x)(Y).   *)
(* ========================================================================== *)

brK0Field::usage = "brK0Field[T] is the coefficient field of k0 for the tower record T: fieldFn[T[\"syms\"]] (the tower symbols are scalars; x is the curve variable). brK0Field[T, {t, mu, alpha}] adjoins the algebraic number alpha = t mod mu(t) (a residue field).";
brK0Field[T_Association] := fieldFn[T["syms"]];
brK0Field[T_Association, ext_List] := fieldFn[T["syms"], ext];

brK0Scalar::usage = "brK0Scalar[B, e] is the element e w_0 = e of k0 (e in k) as a coordinate list on the basis B.";
brK0Scalar[B_Association, e_] := PadRight[{e}, B["m"]];

brK0YMultiple::usage = "brK0YMultiple[B, e] is the element e Y of k0 (e in k): the coordinate d_1 e on w_1 = Y/d_1.";
brK0YMultiple[B_Association, e_] := If[B["m"] < 2, {e}, ReplacePart[Table[0, {B["m"]}], 2 -> Together[e B["d"][[2]]]]];

brK0FromY::usage = "brK0FromY[F, B, e, Y] is the element given by the expression e, a polynomial in the symbol Y of degree < m with coefficients in k (afFromY).";
brK0FromY[F_Association, B_Association, e_, Y_Symbol] := afFromY[F, B, e, Y];

brK0ToY::usage = "brK0ToY[B, a, Y] is the element a as the expression Sum_i c_i Y^i/d_i in the symbol Y (afToY).";
brK0ToY[B_Association, a_List, Y_Symbol] := afToY[B, a, Y];

brK0Times::usage = "brK0Times[F, B, a, b] is the product a b in k0 (afTimes).";
brK0Times[F_Association, B_Association, a_List, b_List] := afTimes[F, B, a, b];

brK0Inverse::usage = "brK0Inverse[F, B, a] is 1/a in k0 for a nonzero a (rationalElementInverse: the adjugate of the multiplication matrix over the norm).";
brK0Inverse[F_Association, B_Association, a_List] := (
  If[afZeroQ[F, a], throwInternal["brK0InverseZero"]];
  rationalElementInverse[F, B, a]);

brK0D::usage = "brK0D[T, F, B, a] is the tower derivative of the element a of k0: (D c_i + c_i s_i) on the basis, s_i = w_i'/w_i (curveBasis, a function of x), D the tower derivation of T (x' = 1, t_i' = D t_i).";
brK0D[T_Association, F_Association, B_Association, a_List] :=
  Table[rc[F, brTowerD[T, a[[i]]] + a[[i]] B["s"][[i]]], {i, B["m"]}];

brK0Degree::usage = "brK0Degree[a, th] is the degree in theta of the element a of k0[theta] (coordinates polynomial in theta), -Infinity for zero. A coordinate whose denominator involves theta is an internal error.";
brK0Degree[a_List, th_Symbol] := Module[{ds},
  ds = Function[c, With[{t = Together[c]},
      If[!FreeQ[Denominator[t], th], throwInternal["brK0NotPolynomial", <|"c" -> c|>]];
      If[t === 0, -Infinity, Exponent[Numerator[t], th]]]] /@ a;
  Max[ds]];

brK0Coefficient::usage = "brK0Coefficient[F, a, th, j] is the coefficient of theta^j of the element a of k0[theta], an element of k0 (coordinatewise).";
brK0Coefficient[F_Association, a_List, th_Symbol, j_Integer] :=
  rc[F, With[{t = Together[#]}, Coefficient[Expand[Numerator[t]], th, j]/Denominator[t]]] & /@ a;

brK0Monic::usage = "brK0Monic[F, B, a, th] is the element a of k0[theta] divided by its leading coefficient (an element of k0).";
brK0Monic[F_Association, B_Association, a_List, th_Symbol] := Module[{d = brK0Degree[a, th], lc},
  If[d === -Infinity, throwInternal["brK0MonicZero"]];
  lc = brK0Coefficient[F, a, th, d];
  brK0Times[F, B, a, brK0Inverse[F, B, lc]]];

brK0PolyRem::usage = "brK0PolyRem[F, B, a, b, th] is the remainder of a by b (b nonzero) in k0[theta]: schoolbook long division with leading coefficients inverted in k0.";
brK0PolyRem[F_Association, B_Association, a0_List, b_List, th_Symbol] := Module[
  {a = a0, db = brK0Degree[b, th], lcInv, da, q, guard = 0},
  If[db === -Infinity, throwInternal["brK0PolyRemByZero"]];
  lcInv = brK0Inverse[F, B, brK0Coefficient[F, b, th, db]];
  While[(da = brK0Degree[a, th]) >= db,
    If[++guard > 1000, throwInternal["brK0PolyRemNoProgress"]];
    q = afScale[F, th^(da - db), brK0Times[F, B, brK0Coefficient[F, a, th, da], lcInv]];
    a = afMinus[F, a, brK0Times[F, B, q, b]];
    (* the leading term cancels exactly; drop any rounding of the canonical form *)
    If[brK0Degree[a, th] >= da, throwInternal["brK0PolyRemLeadingTerm", <|"a" -> a|>]]];
  a];

brK0PolyGCD::usage = "brK0PolyGCD[F, B, a, b, th] is the monic gcd of a and b in k0[theta] (Euclid over the field k0); the gcd of 0 and 0 is 0.";
brK0PolyGCD[F_Association, B_Association, a0_List, b0_List, th_Symbol] := Module[{a = a0, b = b0, r, guard = 0},
  If[afZeroQ[F, a] && afZeroQ[F, b], Return[afZero[B]]];
  While[!afZeroQ[F, b],
    If[++guard > 1000, throwInternal["brK0PolyGCDNoProgress"]];
    r = brK0PolyRem[F, B, a, b, th];
    a = b; b = r];
  brK0Monic[F, B, a, th]];

brK0Trace::usage = "brK0Trace[Ft, e] is the trace from Q(t)(vars) down to Q(vars) of the rational function e over the extension field Ft (fieldFn with t mod mu): Sum_i of the t^i coefficient of e t^i, each in canonical form (rc), for coordinates that are rational in x (Trager's kTrace takes polynomials).";
brK0Trace[Ft_Association, e_] := Module[{t = Ft["t"], r},
  Together[Sum[r = Together[rc[Ft, e t^i]]; Coefficient[Expand[Numerator[r]], t, i]/Denominator[r],
    {i, 0, fieldDegree[Ft] - 1}]]];
