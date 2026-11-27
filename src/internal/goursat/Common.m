(* ========================================================================== *)
(* goursat/Common.m -- shared utilities of the GeneralisedGoursat package.     *)
(*                                                                            *)
(* Read INTO GeneralisedGoursat`Private` by GeneralisedGoursat.m.  Mirrors the  *)
(* "small utilities" section of the sympy reference goursat.py (function for   *)
(* function, camelCased): K, is_zero, lin, polyrem, descend, coords,           *)
(* normalize_algebraic, tidy.  Additions with no sympy counterpart: the        *)
(* failure protocol (ggFail / Throw), the linear solver used where sympy calls  *)
(* sp.solve / LUsolve, and numDen (sympy's fraction).                          *)
(* ========================================================================== *)

(* ---- failure protocol ---------------------------------------------------- *)
(* A stage that cannot continue throws ggFail[class, reason] with tag          *)
(* "GoursatFail".  class is "Decline" (out of scope, a resource bound, or a    *)
(* sub-step that sympy would leave unevaluated) or "Internal" (an invariant    *)
(* of the algorithm failed -- always a bug, reported by the test harness).     *)
ggThrow[class_String, reason_String] := Throw[ggFail[class, reason], "GoursatFail"];
ggDecline[reason_String] := ggThrow["Decline", reason];
ggInternal[reason_String] := ggThrow["Internal", reason];

(* ---- tracing -------------------------------------------------------------- *)
(* stage[label, e] evaluates e; with GeneralisedGoursat`$GoursatDebug True it   *)
(* also prints the label and the wall time to stderr.                           *)
SetAttributes[stage, HoldRest];
stage[label_String, e_] := If[TrueQ[GeneralisedGoursat`$GoursatDebug],
  Module[{tt, v},
    WriteString["stderr", "GeneralisedGoursat: ", label, " ...\n"];
    {tt, v} = AbsoluteTiming[e];
    WriteString["stderr", "GeneralisedGoursat: ", label, " ", ToString[NumberForm[tt, {6, 3}]], " s\n"];
    v],
  e];

(* ---- canonical forms ----------------------------------------------------- *)
(* K: canonical form of a rational function over Q(algebraic numbers, and     *)
(* square roots of parameter polynomials) -- sympy: cancel(together(e),       *)
(* extension=True).  Exact rationals pass through untouched: they dominate the *)
(* linear algebra of the rational cases.                                       *)
gK[e_Integer] := e;
gK[e_Rational] := e;
gK[e_] := With[{c = Cancel[Together[e], Extension -> Automatic]},
  (* Together / Cancel keep a numeric algebraic denominator as written (the  *)
  (* reference's form), so a chain of field operations grows its constants  *)
  (* without bound; sympy's field elements are always reduced.  The core's  *)
  (* canonical form over the number field reduces every constant.           *)
  If[FreeQ[c, Complex | Power[_?NumericQ, _Rational] | _Root], c,
    With[{r = Internal`NumberFieldCanonical[c]}, If[r === $Failed, c, r]]]];

(* is_zero.  A purely numeric (algebraic-number) residue that Together cannot  *)
(* collapse to 0 is decided by PossibleZeroQ, which is exact on algebraic       *)
(* numbers.                                                                    *)
isZero[e_] := Module[{r = gK[e]},
  Which[r === 0, True,
        NumericQ[r], numZeroQ[r],
        True, False]];

(* sympy's fraction(K(e)): {numerator, denominator}, both polynomials. *)
numDen[e_] := Module[{r = gK[e]}, {Numerator[r], Denominator[r]}];

(* The linear form vanishing at g (the constant 1 when g is Infinity). *)
lin[g_, t_] := If[g === Infinity, 1, t - g];

(* Descending coefficient list of a polynomial in t (sympy Poly.all_coeffs). *)
allCoeffs[f_, t_] := Module[{c = CoefficientList[Expand[f], t]},
  If[c === {}, {0}, Reverse[c]]];

(* Is v an exact rational number? *)
ratQ[v_] := MatchQ[v, _Integer | _Rational];

(* ---- polynomial remainder over an arbitrary coefficient field ------------ *)
(* Remainder of f modulo m as polynomials in t; the coefficients are rational  *)
(* functions of other symbols over Q(algebraics), canonicalised by gK at every *)
(* step.  Returns the dm coefficients of t^(dm-1), ..., t^0 (dm = deg m).      *)
polyrem[f_, m_, t_] := Module[{cm, cf, dm, q},
  cm = gK /@ allCoeffs[m, t];
  While[cm =!= {} && First[cm] === 0, cm = Rest[cm]];
  dm = Length[cm] - 1;
  cf = gK /@ allCoeffs[f, t];
  While[Length[cf] - 1 >= dm && cf =!= {},
    q = gK[cf[[1]]/cm[[1]]];
    Do[cf[[i]] = gK[cf[[i]] - q cm[[i]]], {i, 2, dm + 1}];
    cf = Rest[cf]];
  Join[ConstantArray[0, dm - Length[cf]], cf]];

(* Does the polynomial f divide g (in t)? *)
dividesQ[f_, g_, t_] := AllTrue[polyrem[g, f, t], # === 0 &];

(* ---- descent along a rational map ---------------------------------------- *)
(* g(t) is known to be a rational function of X(t): return G with             *)
(* G(X(t)) = g(t).  Uses num(g) = G den(g) modulo N(t) - x D(t), X = N/D; no   *)
(* inversion.  Throws Internal "not invariant" when g is not a function of X.  *)
descend[g0_, X_, t_, x_] := Module[{g, nn, dd, gn, gd, mm, rn, rd, G = None},
  g = gK[g0];
  If[g === 0, Return[0, Module]];
  {nn, dd} = numDen[X];
  {gn, gd} = numDen[g];
  mm = nn - x dd;
  rn = polyrem[gn, mm, t];
  rd = polyrem[gd, mm, t];
  Do[If[rd[[i]] =!= 0, G = gK[rn[[i]]/rd[[i]]]; Break[]], {i, Length[rd]}];
  If[G === None || !isZero[(G /. x -> X) - g], ggInternal["descend: not invariant"]];
  G];

(* Coordinates of F in a basis of Q(t) over Q(X(t)); Length[basis] = deg X. *)
coords[F_, basis_List, X_, t_, x_] := Module[
  {nn, dd, mm, fn, fd, bden = 1, bd, den, rhs, cols, sol},
  {nn, dd} = numDen[X];
  mm = nn - x dd;
  {fn, fd} = numDen[F];
  (* common denominator of the basis (no gcds needed) *)
  Do[bd = Denominator[gK[b]];
     If[!FreeQ[bd, t] && !dividesQ[bd, bden, t], bden = Expand[bden bd]],
     {b, basis}];
  den = Expand[fd bden];
  rhs = polyrem[Expand[fn bden], mm, t];
  cols = Table[polyrem[Expand[gK[b den]], mm, t], {b, basis}];
  sol = gK /@ linSolve[Transpose[cols], rhs];
  If[!isZero[Sum[(sol[[i]] /. x -> X) basis[[i]], {i, Length[basis]}] - F],
     ggInternal["coords: basis decomposition failed"]];
  sol];

(* ---- linear algebra -------------------------------------------------------- *)
(* A particular solution of A.u == b over Q(params)(algebraics), free unknowns  *)
(* set to 0 (as sympy's solve(dict=True) leaves them, and the Hermite driver    *)
(* zeroes them).  Throws Internal on an inconsistent system.                    *)
linSolve[A_List, b_List] := Module[{s},
  If[A === {} || A[[1]] === {}, Return[{}, Module]];
  s = Quiet[LinearSolve[A, b]];
  If[!ListQ[s], ggInternal["inconsistent linear system"]];
  s];

(* Solve the polynomial identity expr == 0 in t, expr linear in unknowns, for *)
(* the unknowns (sympy: solve(Poly(expr, t).all_coeffs(), unknowns)).  Returns *)
(* a list of rules, free unknowns -> 0.                                        *)
solveIdentity[expr_, unknowns_List, t_] := Module[{eqs, A, b, z, s},
  If[unknowns === {}, Return[{}, Module]];
  eqs = CoefficientList[Expand[expr], t];
  If[eqs === {}, Return[Thread[unknowns -> 0], Module]];
  z = Thread[unknowns -> 0];
  A = Table[D[eq, u], {eq, eqs}, {u, unknowns}];
  b = -(eqs /. z);
  s = linSolve[A, b];
  Thread[unknowns -> s]];

(* Adjugate of a 3 x 3 matrix (no division). *)
adj3[{{a_, b_, c_}, {d_, e_, f_}, {g_, h_, i_}}] :=
  {{e i - f h, c h - b i, b f - c e},
   {f g - d i, a i - c g, c d - a f},
   {d h - e g, b g - a h, a e - b d}};

(* The roots of a linear or quadratic polynomial, by the quadratic formula     *)
(* (sympy: roots), each in gK form, the square factors of the discriminant      *)
(* pulled out of the radical (sqrtPull).                                       *)
quadRoots[f_, t_] := Module[{c = CoefficientList[Expand[f], t], s},
  If[Length[c] == 2, Return[{gK[-c[[1]]/c[[2]]]}, Module]];
  s = sqrtPull[gK[c[[2]]^2 - 4 c[[1]] c[[3]]]];
  {gK[(-c[[2]] - s)/(2 c[[3]])], gK[(-c[[2]] + s)/(2 c[[3]])]}];

(* A square root of d with its square factors pulled out (either sign: callers *)
(* use both roots, or either branch of a conic): Sqrt[-3 (k - 1)^2] is          *)
(* (k - 1) Sqrt[-3], a numeric generator, not a symbolic radicand that may      *)
(* resurface in another spelling.                                               *)
sqrtPull[d_] := Module[{r = Together[d], f},
  If[NumericQ[r], Return[Sqrt[r], Module]];
  f[p_] := Module[{fl = FactorList[p]},
    {Times @@ (#[[1]]^Quotient[#[[2]], 2] & /@ fl), Times @@ (#[[1]]^Mod[#[[2]], 2] & /@ fl)}];
  With[{nu = f[Numerator[r]], de = f[Denominator[r]]},
    (nu[[1]]/de[[1]]) Sqrt[Expand[nu[[2]] de[[2]]]]/de[[2]]]];

rootsOf[p_, t_] := Module[{fl, out = {}, f, mult, d},
  d = Exponent[p, t];
  If[d < 1, Return[{}, Module]];
  (* algebraic-number coefficients: only the quadratic formula.  Coefficients   *)
  (* in Q(parameters) factor over Q(parameters) like rational ones (beyond      *)
  (* goursat.py, whose all_roots needs numbers).                                *)
  If[!FreeQ[CoefficientList[Expand[p], t], Complex | Power[_?NumericQ, _Rational] | _Root],
    If[d <= 2, Return[quadRoots[p, t], Module], ggDecline["roots: algebraic coefficients"]]];
  fl = FactorList[Together[p]];
  Do[{f, mult} = fm;
     Which[Exponent[f, t] == 0, Null,
           Exponent[f, t] <= 2, out = Join[out, Flatten[ConstantArray[quadRoots[f, t], mult]]],
           True, ggDecline["roots: irreducible factor of degree > 2"]],
     {fm, fl}];
  out];

(* The rational roots of a polynomial with rational coefficients, sorted. *)
rationalRoots[p_, t_] := Module[{fl},
  fl = FactorList[p];
  Sort[Flatten[Cases[fl, {f_, _} /; Exponent[f, t] == 1 :> (t /. Solve[f == 0, t])]]]];

(* ---- the algebraic function field y^n = R ---------------------------------- *)
(* Write a rational expression in (t, y), y^n = R, as Sum_j c_j(t) y^j. *)
normalizeAlgebraic[E0_, R_, t_, y_, n_] := Module[{c = normalizeCoeffs[E0, R, t, y, n]},
  Sum[Factor[c[[i + 1]]] y^i, {i, 0, n - 1}]];

(* The coefficient vector {c_0, ..., c_(n-1)} (gK form) of normalizeAlgebraic. *)
normalizeCoeffs[E0_, R_, t_, y_, n_] := Module[{E1, num, den, red, nv, dv, M, c},
  E1 = gK[E0];
  If[E1 === 0, Return[ConstantArray[0, n], Module]];
  {num, den} = {Numerator[E1], Denominator[E1]};
  red[f_] := Module[{v = ConstantArray[0, n], cl = CoefficientList[Expand[f], y]},
    Do[v[[Mod[j - 1, n] + 1]] += cl[[j]] R^Quotient[j - 1, n], {j, Length[cl]}];
    v];
  nv = red[num]; dv = red[den];
  (* multiplication by den on the basis 1, y, ..., y^(n-1) *)
  M = ConstantArray[0, {n, n}];
  Do[M[[Mod[i + j, n] + 1, i + 1]] += dv[[j + 1]] R^Quotient[i + j, n],
     {i, 0, n - 1}, {j, 0, n - 1}];
  If[den === 1 || FreeQ[den, y],
    gK /@ (nv/den),
    gK /@ linSolve[Map[gK, M, {2}], nv]]];

(* Heads that mark a transcendental (non-algebraic) term of an antiderivative. *)
transcendentalQ[e_] := !FreeQ[e, _Log | _ArcTan | _ArcTanh | _ArcCot | _ArcCoth |
                                 _RootSum | _Integrate | ggIntegral | ggDeferred];

(* Collect the purely algebraic terms of an antiderivative and reduce them;    *)
(* every argument of a logarithm or inverse trigonometric function is brought   *)
(* to its normal form in Q(t)[y]/(y^n - R) (tidyArgs).                          *)
tidy[anti_, R_, t_, y_, n_] := Module[{terms, alg = 0, rest = 0},
  terms = If[Head[anti] === Plus, List @@ anti, {anti}];
  Do[If[transcendentalQ[term], rest += term, alg += term], {term, terms}];
  normalizeAlgebraic[alg, R, t, y, n] + tidyArgs[rest, R, t, y, n]];

(* Normal form of an argument a(t, y); a slot-dependent argument (the Log[m - #] *)
(* of a RootSum) has its slot-free part normalised.  Back-substitution of a      *)
(* second-level answer (X -> H/R, Y -> T/(R y)) otherwise leaves nested          *)
(* fractions whose Together dominates verification.  Declines (keeps a) when    *)
(* the argument is not rational in t and y.                                     *)
normArg[a_, R_, t_, y_, n_] := Module[{p, q, r},
  If[FreeQ[a, _Slot],
    (* already a polynomial in y over Q(t) (every first-level uniformiser): kept *)
    (* as is -- its normal form over a field like Q(I, Sqrt[3], 2^(1/3)) costs   *)
    (* seconds of extension gcds and changes nothing                             *)
    If[FreeQ[a, Power[b_, e_?Negative] /; !FreeQ[b, y]] || FreeQ[Denominator[Together[a]], y],
      Return[a, Module]];
    r = Quiet[Catch[normalizeAlgebraic[a, R, t, y, n], "GoursatFail"]];
    Return[If[MatchQ[r, _ggFail], a, r], Module]];
  If[Head[a] =!= Plus, Return[a, Module]];
  p = Select[a, FreeQ[#, _Slot] &];
  With[{p2 = normArg[p, R, t, y, n]}, If[p2 === p, a, p2 + (a - p)]]];

(* Structurally, rebuilding RootSum and Function nodes by Apply so that each new *)
(* argument is evaluated before it goes back in: a rule's right-hand side is    *)
(* not evaluated inside a held RootSum body.                                    *)
tidyArgs[e_, R_, t_, y_, n_] := Which[
  AtomQ[e] || FreeQ[e, Log | ArcTan | ArcTanh | ArcCot | ArcCoth | RootSum], e,
  (* a node whose arguments are unchanged is returned as it is: rebuilding it    *)
  (* re-evaluates it (a RootSum over Q(k)(I Sqrt[3]) is renormalised), which was *)
  (* most of the cost of a tidy that changes nothing                            *)
  MatchQ[Head[e], Log | ArcTan | ArcTanh | ArcCot | ArcCoth] && Length[e] == 1 &&
    !FreeQ[e[[1]], t | y],
    With[{a2 = normArg[e[[1]], R, t, y, n]}, If[a2 === e[[1]], e, Head[e][a2]]],
  Head[e] === RootSum && Length[e] == 2,
    With[{b2 = tidyArgs[e[[2]], R, t, y, n]}, If[b2 === e[[2]], e, RootSum @@ {e[[1]], b2}]],
  Head[e] === Function,
    With[{b2 = tidyArgs[Last[e], R, t, y, n]},
      If[b2 === Last[e], e, Function @@ Append[Most[List @@ e], b2]]],
  True,
    With[{as = List @@ e},
      With[{as2 = tidyArgs[#, R, t, y, n] & /@ as}, If[as2 === as, e, Head[e] @@ as2]]]];
