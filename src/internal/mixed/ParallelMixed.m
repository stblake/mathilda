(* ::Package:: *)

(* ParallelMixed.wl -- a Mathilda port of parallel_mixed.py

   Parallel (Risch-Norman) integration over a simple radical y^2 = q sitting
   anywhere in a tower of transcendental generators, after

     S. Blake, Parallel Integration over Simple Radical Extensions II: Mixed Towers.

   Elements of L = K(y) are pairs {a, b} = a + b y for y^2 = q, and for a
   radical y^m = q of degree m >= 3 the m-tuples of coordinates on the Trager
   basis w_i = y^i/E_i (Part I, Prop. 3.2).  The tower is
     Tower[gens, derivs, q, m]   with derivs[[i]] = D(gens[[i]]) as a tuple,
   q = None for m = 1 (no radical).

   Stage 1 covers: pair arithmetic, the derivation, valuations and the
   classification (normal/special, delta_P), Hermite denominators, exact-order
   residues on every sheet with the exact constancy certificate, the quadratic
   y-split, tower specials, the ansatz and the linear system, and verification.
   Stage 2 adds the norm search, vanishing orders, deep residues (n = 1),
   elliptic torsion (division polynomials, Miller functions) and the v_oo
   residue; stage 3 the units (continued fractions over Q and GF(p)), the
   non-torsion certificates (Proposition 9.4: for the class at infinity by
   the continued fraction over GF(p), for a residue divisor at finite
   places by Cantor's algorithm in the Jacobian over GF(p), any genus), and
   the holomorphic-remainder certificate, with
   the special exponents and degree bounds of Algorithm 6 (Section 8.2:
   proved where Proposition 8.11 applies, the classical guess elsewhere).  Radicals of degree m >= 3 (y^m = q,
   q m-th-power-free) are handled on the Trager basis w_i = y^i/E_i:
   residue classes with the m x m norm, realisation by Hensel lifting and
   linear algebra with the poles distributed over the places at infinity,
   unit candidates from a bounded divisor search at infinity (UnitsGeneral).
   Not covered for m >= 3: S'-units over special primes, residue classes on
   several but not all sheets over one root, partially ramified branch
   primes, and the two m = 2 certificates (non-torsion, holomorphic
   remainder).  The port is complete.

   The function names are the Python names in CamelCase; the verbose trace
   lines are the same, so sessions can be compared between the two.          *)

BeginPackage["ParallelMixed`"];

Tower::usage = "Tower[gens, derivs, q, m] builds the differential tower L = F(gens)(y), y^m = q (m = 2 by default), over which ParallelIntegrateMixed works. For m >= 3 the elements are m-tuples of coordinates on the Trager basis w_i = y^i/E_i, q must be m-th-power-free, and the Association also carries \"m\", \"n\" (the tuple length), \"E\", the multiplication table \"mul\", the logarithmic derivatives \"Lam\" (D w_i = Lam_i w_i) and \"Dw\" = {D w_0, ..., D w_{m-1}}.
  gens   -- list of generator symbols {t1, ..., tn}; t1 = x is the base variable and must satisfy D(x) = 1.
  derivs -- list of pairs {a, b}, one per generator, giving D(ti) = a + b y with a, b rational functions of the generators. The pair for x is {1, 0}. Generators above the radical may have b != 0 (their derivatives involve y): this is the 'mixed' setting.
  q      -- the radicand, a squarefree polynomial in the generators (usually in x alone); omit it, or pass None, for m = 1 (a purely transcendental tower, e.g. after a pure root has been flattened into the basis).
The result is an Association with keys \"gens\", \"derivs\", \"q\" and, when q is present, \"Dy\" = D(y) as a pair, computed from Dq = Q0 + Q1 y as Dy = Dq/(2y) = Q1/2 + (Q0/(2q)) y.";

TowerD::usage = "TowerD[T, u] is the derivation of the tower T applied to the pair u = {a, b} = a + b y, returned as a pair. It is the sum of D applied coefficientwise to a and b (through every generator, using the chain rule with T[\"derivs\"]) and of b D(y), with y^2 = q used to keep the result on the basis {1, y}.";

ClassifyPrime::usage = "ClassifyPrime[T, p] classifies the height-one primes of the working ring lying over the irreducible polynomial p (Definition 4.3 of the paper) and returns {branch, eta, delta, special}:
  branch  -- True if p divides the radicand q (the primes over p are ramified, e_P = 2);
  eta     -- the order of the denominator divisor of D at P: eta = max(0, -v_P(D t_i), -v_P(D y));
  delta   -- the pole order of the derivation, delta = 1 + eta, so that v_P(D g) = v_P(g) - delta for v_P(g) != 0 at a normal prime (the valuation lemma);
  special -- True if P is special (v_P(D pi) > -eta for the uniformiser pi = p, or pi = y at a branch prime), False if P is normal. Special primes are carried as candidate logands and may appear in the denominator of the integral with a guessed exponent; normal primes obey the valuation lemma and carry residues.";

CertifyNonconstant::usage = "CertifyNonconstant[tau, gens] decides exactly whether a residue tau -- an element of the residue field kappa(P), given as an expression in the generators, algebraic numbers and at most one algebraic function Sqrt[A] -- lies outside the constant field Fbar. It returns
  True      -- tau is certified NOT to be a constant (this is the residue certificate of non-elementarity, Proposition 9.2(a));
  False     -- tau is a constant, possibly in disguise (e.g. (x + Sqrt[x^2+1]) (x - Sqrt[x^2+1]));
  Undecided -- the expression contains more than one algebraic function, or a non-square-root algebraic function, and no verdict is issued.
The decision reduces tau to the canonical form (a + b Sqrt[A]) / D0 with a, b, D0 rational in the generators; tau is constant iff b = 0 and a/D0 is free of the generators. A constant-in-disguise must never be mistaken for an obstruction, so a syntactic free-variable test is not used.";

ParallelIntegrateMixed::usage = "ParallelIntegrateMixed[integrand, x, opts] integrates an expression in x built from rational operations, one square root of a polynomial (the radical y^2 = q), pure roots that flatten (Lemma 3.2: a root whose base is linear in a generator with a constant coefficient, such as Sqrt[x], Sqrt[x + Log[x]], (x + Exp[x])^(1/3), Sqrt[Log[x]]), and logarithms, exponentials and tangents of expressions in the field generated so far. The tower is constructed automatically (BuildTower), the integral computed by the pair form below, and the result returned in terms of x. Option \"Verify\" -> True additionally checks D[result, x] against the integrand numerically at three points.
ParallelIntegrateMixed[{f0, f1}, T, opts] integrates f = f0 + f1 y over the tower T by the parallel (Risch-Norman) method of the paper (Algorithm 4), with y = Sqrt[q]. It returns one of
  an expression in the generators and Sqrt[q]   -- the integral, verified inside the routine by re-differentiation;
  {\"not elementary\", p, tau}                    -- certified: the residue tau at a normal prime over p is not a constant (Proposition 9.2(a));
  {\"failed\", reason, ...}                       -- no verdict: the reason names the guessed input in play (in stage 1 also the mechanisms not yet ported: norm search, torsion, deep residues, units, the place at infinity, the holomorphic-remainder certificate).
Options: \"SplitSpecials\" -> Automatic (default) retries a failed system with every special prime split into its linear factors over the algebraic closure, each with its own coefficient (needed e.g. for Sqrt[Tan[x]]); True always splits, False never; \"Bounds\" -> {b1, ..., bn} fixes the degree bounds of the numerator of the rational part in each generator (default: Algorithm 6 of the paper -- proved by Proposition 8.11 where it applies, else the maximal degree occurring + 2, raised on retry); \"Verbose\" -> True prints the classification of every prime, the residues on every sheet, the candidate logands and the special primes, in the same format as the SymPy implementation so sessions can be compared.
Steps (Algorithm 4): factor the denominator of f; for each prime, ClassifyPrime; specials become candidate logands and denominator factors; normal primes with a pole deeper than delta contribute Hermite factors to the denominator; normal primes with a pole of order exactly delta carry residues, computed on every sheet, certified constant or not, and realised as logands (equal residues -> the polynomial p; residues split by the curve -> the quadratic y-split; at a prime of degree > 4 the residues are computed in the residue field K[g]/(p) -- the residue polynomial of Rothstein and Trager -- and each residue class is realised from its sheet polynomial by Hensel lifting and linear algebra over K(c)); then the tower specials are added, the residual f - Sum c_i D u_i / u_i is formed, and the linear system for the rational part and the unknown logand coefficients is solved.";

BuildTower::usage = "BuildTower[integrand, x] constructs the differential tower for a surface expression. Returns {T, {f0, f1}, back, Y}: the Tower T, the integrand as a pair f0 + f1 y on T, the list of back-substitution rules mapping the generator symbols (and flattened roots) to their surface expressions in x, and the symbol Y used for y during the conversion. Processing is innermost-first: Log[a] -> new generator t with D t = D a / a; Exp[a] -> t with D t = t D a; Tan[a] -> t with D t = (1 + t^2) D a; a root base^(k/m) whose base is c g + r with g a generator, c a constant and r free of g is flattened, g -> (u^m - r)/c, D u = D(base)/(m u^(m-1)); any other root base^(k/m) with polynomial base becomes the radical y^m = base (pairs on {1, y} for m = 2, m-tuples on the Trager basis for m >= 3). Fails with a message if two distinct unflattenable radicals occur or if a radicand is not an m-th-power-free polynomial.";
BuildTower::radicals = "More than one unflattenable radical (`1` and `2`); the implementation handles one simple radical.";
BuildTower::mixedroots = "Roots of degree `1` and `2` of the same radicand `3`.";
BuildTower::powerfree = "The radicand `1` is not `2`-th-power-free.";
BuildTower::radical3 = "The derivative of `1` needs Sqrt[`2`], but the tower has the radical y^`3` = `2`.";
BuildTower::radicand = "The radicand `1` is not a squarefree polynomial in the generators.";
BuildTower::residual = "Could not reduce the integrand to a rational function of the generators and y: `1`.";
Undecided;

Begin["`Private`"];

(* ------------------------------------------------------------ pair arithmetic *)
(* Can[e]
   Canonical form of a rational function with algebraic-number coefficients:
   Together followed by Cancel with Extension -> Automatic, so that common
   factors are removed exactly even when the coefficients involve I, Sqrt[3],
   Root objects, etc.  Every pair operation runs its results through Can; it
   is the counterpart of SymPy's cancel(..., extension=True), the step that
   made the algebraic-number-heavy examples (Guenther) tractable.            *)
(* Field-aware canonicalisation.  When the coefficients involve algebraic
   constants (Root objects, numeric radicals, Gaussian Complex numbers),
   Cancel[Together[.], Extension -> Automatic] on the RAW radicals re-derives the
   algebraic tower every call and runs the polynomial GCD over it -- the dominant
   cost of the compositum cases (Charlwood P8 ~63% of its time, P4 ~33%).  Instead
   map the atoms to AlgebraicNumber[theta] in one common number field K = Q(theta)
   (FieldData = ToNumberField + a radical basis, MEMOISED by atom set so the field
   is built once, ~2 ms), run the native FLINT field Together/Cancel over K, and
   map the field coefficients back to radicals via FieldData's `back`.  Measured on
   P8's worst call: value-identical, 0.48 s -> 0.0025 s (~190x), and a more reduced
   result (LeafCount 390 -> 63) that also speeds every downstream op.  A rational
   input, a tiny one (below the size gate, where the field detour costs more than it
   saves), or one whose atoms do not generate a number field falls back to the
   original Extension -> Automatic path unchanged.  $CanFieldEnabled = False forces
   that path everywhere (A/B differential). *)
$CanFieldEnabled = True;
$CanFieldMemo = <||>;
CanRaw[e_] := Cancel[Together[e], Extension -> Automatic];
Can[e_] := Module[{atoms, key, fd},
  If[! TrueQ[$CanFieldEnabled], Return[CanRaw[e]]];
  atoms = DeleteDuplicates[Cases[e, _Root | _Complex | Power[_?NumericQ, _Rational], {0, Infinity}]];
  If[atoms === {} || LeafCount[e] < 40, Return[CanRaw[e]]];
  key = Sort[atoms];
  fd = If[KeyExistsQ[$CanFieldMemo, key], $CanFieldMemo[key], $CanFieldMemo[key] = FieldData[atoms]];
  If[! MatchQ[fd, {_, _}], Return[CanRaw[e]]];
  Cancel[Together[e /. fd[[1]]], Extension -> Automatic] /. a_AlgebraicNumber :> fd[[2]][a]];

(* $analyses: the once-per-integrand analyses of iPIM (Steps 1--14 of Algorithm 4), keyed by Hash[{f0, T}] *)
$analyses = <||>;
(* Padd[u, v], Pscale[s, u], Pmul[u, v, q], Pdiv[u, v, q]
   Arithmetic in L = K(y), y^2 = q, on pairs {a, b} = a + b y.
     Padd     componentwise sum;
     Pscale   multiplication by a scalar s of K;
     Pmul     (a0 + a1 y)(b0 + b1 y) = (a0 b0 + a1 b1 q) + (a0 b1 + a1 b0) y;
     Pdiv     u / v = u * conj(v) / N(v) with conj(v) = v0 - v1 y and the norm
              N(v) = v0^2 - q v1^2 in K (so the quotient is again a pair).
   For m = 1 (q = None) the second component is identically 0 and only Padd
   and Pscale are used.                                                      *)
Padd[u_, v_] := Can /@ (u + v);
Pscale[s_, u_] := Can /@ (s u);
Pmul[u_, v_, q_] := {Can[u[[1]] v[[1]] + u[[2]] v[[2]] q], Can[u[[1]] v[[2]] + u[[2]] v[[1]]]};
Pdiv[u_, v_, q_] := Module[{nrm = Can[v[[1]]^2 - q v[[2]]^2], w},
  w = Pmul[u, {v[[1]], -v[[2]]}, q];
  {Can[w[[1]]/nrm], Can[w[[2]]/nrm]}];

(* ------------------------------------------------------------ the tower *)
(* Tower[gens, derivs, q]  -- see the usage message.
   Internally: the derivative pairs are normalised with Can, and when q is
   present the pair Dy is stored so that TowerD can differentiate y without
   recomputing Dq each time.                                                 *)
Tower[gens_List, derivs_List, q_: None, m0_: 2] := Module[{T, Q, m, n, sfl, Qs = <||>, EE, mul, Lam, Qh, Dw},
  m = If[q === None, 1, m0]; n = If[m <= 2, 2, m];
  T = <|"gens" -> gens, "derivs" -> (Can /@ # & /@ (PadRight[#, n] & /@ derivs)), "q" -> q,
        "m" -> m, "n" -> n, "E" -> ConstantArray[1, n]|>;
  If[q =!= None && m == 2,
    Q = Dhat[T, q];                           (* Dq = Q0 + Q1 y *)
    T = Append[T, "Dy" -> {Can[Q[[2]]/2], Can[Q[[1]]/(2 q)]}];  (* Dy = Dq/(2y) = Q1/2 + Q0 y/(2q) *)
    T = Append[T, "Dw" -> {{0, 0}, T["Dy"]}]];
  If[q =!= None && m >= 3,
    (* the Trager basis w_i = y^i/E_i, E_i = prod_l Q_l^Floor[i l/m] for the
       squarefree decomposition q = lc prod_l Q_l^l, and the multiplication
       table w_i w_k = c_ik w_((i+k) mod m) with polynomial c_ik *)
    sfl = Select[FactorSquareFreeList[Expand[q]], ! FreeQ[#[[1]], Alternatives @@ gens] &];
    Do[Qs[fc[[2]]] = Lookup[Qs, fc[[2]], 1] fc[[1]], {fc, sfl}];
    EE = Table[Expand[Times @@ KeyValueMap[#2^Floor[i #1/m] &, Qs]], {i, 0, m - 1}];
    mul = Table[If[i + k < m,
        {i + k, Expand[Cancel[EE[[i + k + 1]]/(EE[[i + 1]] EE[[k + 1]])]]},
        {i + k - m, Expand[Cancel[q EE[[i + k - m + 1]]/(EE[[i + 1]] EE[[k + 1]])]]}],
      {i, 0, m - 1}, {k, 0, m - 1}];
    T = Join[T, <|"E" -> EE, "mul" -> mul|>];
    Qh = Dhat[T, q];
    Lam = Table[With[{lam = Pscale[(i/m)/q, Qh]},
        If[EE[[i + 1]] =!= 1, Padd[lam, Pscale[-1/EE[[i + 1]], Dhat[T, EE[[i + 1]]]]], lam]], {i, 0, m - 1}];
    T = Append[T, "Lam" -> Lam];
    Dw = Table[TowerD[T, TUnit[T, i]], {i, 0, m - 1}];
    T = Join[T, <|"Dw" -> Dw, "Dy" -> Dw[[2]]|>]];
  T];

(* tuples: padding, unit vectors, scalars *)
TPad[T_, u_] := PadRight[u, T["n"]];
TUnit[T_, i_] := UnitVector[T["n"], i + 1];
TScalar[T_, a_] := PadRight[{a}, T["n"]];
(* TimesW[T, t, i]: the tuple t times the basis element w_i (m >= 3) *)
TimesW[T_, t_, i_] := Module[{m = T["m"], r, idx, coef},
  r = ConstantArray[0, m];
  Do[If[t[[k + 1]] === 0, Continue[]];
    {idx, coef} = T["mul"][[k + 1, i + 1]];
    r[[idx + 1]] += t[[k + 1]] coef, {k, 0, m - 1}];
  Can /@ r];
(* TMul, TDiv, TNorm: the arithmetic of L on tuples for every m *)
TMul[T_, u_, v_] := Module[{m = T["m"], r, idx, coef, uu, vv},
  If[T["q"] === None, Return[{Can[u[[1]] v[[1]]], 0}]];
  If[m == 2, Return[Pmul[u, v, T["q"]]]];
  uu = TPad[T, u]; vv = TPad[T, v]; r = ConstantArray[0, m];
  Do[If[uu[[i + 1]] === 0 || vv[[k + 1]] === 0, Continue[]];
    {idx, coef} = T["mul"][[i + 1, k + 1]];
    r[[idx + 1]] += uu[[i + 1]] vv[[k + 1]] coef, {i, 0, m - 1}, {k, 0, m - 1}];
  Can /@ r];
(* TMatrix[T, v]: M[[j, k]] = coordinate j of v w_k, the multiplication matrix *)
TMatrix[T_, v_] := Transpose[Table[TimesW[T, TPad[T, v], k], {k, 0, T["m"] - 1}]];
TDiv[T_, u_, v_] := Module[{m = T["m"], uu, vv},
  If[T["q"] === None, Return[{Can[u[[1]]/v[[1]]], 0}]];
  If[m == 2, Return[Pdiv[u, v, T["q"]]]];
  uu = TPad[T, u]; vv = TPad[T, v];
  If[AllTrue[Rest[vv], # === 0 &], Return[Can /@ (uu/vv[[1]])]];
  Can /@ LinearSolve[TMatrix[T, vv], uu]];
TNorm[T_, u_] := Which[T["q"] === None, Can[u[[1]]], T["m"] == 2, Can[u[[1]]^2 - T["q"] u[[2]]^2],
  True, Can[Det[TMatrix[T, TPad[T, u]]]]];
(* FromY[T, expr, Y]: the coordinates of an expression rational in the
   generators and Y, Y^m = q (ToPair for m = 2; for m >= 3 numerator and
   denominator are reduced modulo Y^m - q, y^i = E_i w_i, and the quotient is
   taken in the basis) *)
FromY[T_, expr_, Y_] := Module[{m = T["m"], e, n, d, red, tup},
  If[T["q"] === None, Return[{Can[expr], 0}]];
  If[m == 2, Return[ToPair[expr, Y, T["q"]]]];
  e = Together[expr]; {n, d} = {Numerator[e], Denominator[e]};
  red[z_] := PolynomialRemainder[Expand[z], Y^m - T["q"], Y];
  tup[z_] := Table[Can[Coefficient[z, Y, i] T["E"][[i + 1]]], {i, 0, m - 1}];
  TDiv[T, tup[red[n]], tup[red[d]]]];
ToY[T_, u_, Y_] := Which[T["q"] === None, u[[1]], T["m"] == 2, u[[1]] + u[[2]] Y,
  True, Sum[TPad[T, u][[i + 1]] Y^i/T["E"][[i + 1]], {i, 0, T["m"] - 1}]];
BranchMult[T_, p_] := If[T["q"] === None, 0, Vp[T["q"], p, T["gens"]]];
(* RootOfUnity[m]: a primitive m-th root of unity in radicals (-1 for m = 2);
   MRoot[a, m]: the principal m-th root of a factored constant (Sqrt for
   m = 2), the real root for a negative rational and odd m *)
RootOfUnity[m_] := If[m == 2, -1, ComplexExpand[Exp[2 Pi I/m]]];
MRoot[a0_, m_] := Module[{a = Factor[a0]},
  If[m == 2, Return[PowerExpand[Sqrt[a]]]];
  If[MatchQ[a, _Integer | _Rational] && a < 0 && OddQ[m], -PowerExpand[(-a)^(1/m)], PowerExpand[a^(1/m)]]];
CompositionsOf[n_, 1] := {{n}};
CompositionsOf[n_, s_] := Join @@ Table[Prepend[#, first] & /@ CompositionsOf[n - first, s - 1], {first, 0, n}];

(* Dhat[T, c]
   The derivation applied to a scalar c of K (a rational function of the
   generators, no y): D c = Sum_i (d c / d t_i) D(t_i), a pair since the
   D(t_i) are pairs.  This is the chain rule over the tower; the result can
   have a y-component when some generator above the radical differentiates
   through y.                                                                *)
Dhat[T_, c_] := Module[{out = ConstantArray[0, T["n"]]},
  Do[out = Padd[out, Pscale[D[c, T["gens"][[i]]], T["derivs"][[i]]]], {i, Length[T["gens"]]}];
  out];

(* TowerD[T, u]  -- see the usage message.
   D(a + b y) = D a + (D b) y + b D y, assembled with Pmul so that the
   product (D b) y is reduced modulo y^2 = q.                                *)
TowerD[T_, u0_] := Module[{u = TPad[T, u0], out, h1, m = T["m"], t},
  out = Dhat[T, u[[1]]];
  If[T["q"] === None, Return[out]];
  If[m == 2,
    h1 = Dhat[T, u[[2]]];
    out = Padd[out, Pmul[h1, {0, 1}, T["q"]]];
    Return[Padd[out, Pscale[u[[2]], T["Dy"]]]]];
  (* m >= 3: D(a_i w_i) = (D a_i + a_i Lam_i) w_i, the product by the table *)
  Do[If[u[[i + 1]] === 0, Continue[]];
    t = Padd[Dhat[T, u[[i + 1]]], Pscale[u[[i + 1]], T["Lam"][[i + 1]]]];
    out = Padd[out, TimesW[T, t, i]], {i, 1, m - 1}];
  out];

(* ------------------------------------------------------------ valuations *)
(* Mult[poly, p, gens]
   Multiplicity of the irreducible polynomial p in the polynomial poly, by
   repeated exact division (PolynomialReduce with the single divisor p, whose
   remainder vanishes iff p divides).  Returns Infinity for poly = 0.       *)
Mult[poly_, p_, gens_] := Module[{m = 0, P = Expand[poly], qr},
  If[P === 0, Return[Infinity]];
  While[True,
    qr = PolynomialReduce[P, {p}, gens];
    If[qr[[2]] =!= 0, Return[m]];
    m++; P = qr[[1, 1]]]];

(* Vp[expr, p, gens]
   The p-adic valuation v_p on K = F(gens): multiplicity of p in the numerator
   minus multiplicity in the denominator of the cancelled expression;
   Infinity for 0.                                                           *)
Vp[expr_, p_, gens_] := Module[{e = Can[expr], n, d},
  If[e === 0, Return[Infinity]];
  {n, d} = {Numerator[e], Denominator[e]};
  Mult[n, p, gens] - Mult[d, p, gens]];

(* VP[u, p, T, branch]
   The valuation at the primes P of the working ring over p of the pair
   u = a + b y (Lemma 4.4 of the paper, m = 2, q squarefree):
     unramified (p does not divide q):  v_P(u) = min(v_p(a), v_p(b));
     branch     (p divides q):          v_P(u) = min(2 v_p(a), 2 v_p(b) + 1),
   since then v_P(p) = 2 and v_P(y) = 1, so a and b y have valuations of
   different parity and cannot cancel.  For m = 1 it is v_p(a).             *)
VP[u_, p_, T_, branch_] := Module[{g = T["gens"], m = T["m"], l0},
  If[T["q"] === None, Return[Vp[u[[1]], p, g]]];
  If[! branch, Return[Min[Vp[#, p, g] & /@ u]]];
  If[m == 2, Return[Min[2 Vp[u[[1]], p, g], 2 Vp[u[[2]], p, g] + 1]]];
  (* m >= 3, p | Q_l0 with gcd(m, l0) = 1: e_P = m, v_P(w_i) = i l0 mod m *)
  l0 = BranchMult[T, p];
  If[GCD[m, l0] =!= 1, Throw[{"failed", "partially ramified branch prime", p}, "PIM"]];
  Min[Table[m Vp[u[[i + 1]], p, g] + Mod[i l0, m], {i, 0, m - 1}]]];

(* Eta[T, p, branch]
   The order at P of the denominator divisor of the derivation:
   eta_P = max(0, -v_P(D t_i) over the generators, -v_P(D y)).  Together with
   the normality test this is Definition 4.3 in computable form; the pole
   order of D at a normal P is delta_P = 1 + eta_P.                          *)
Eta[T_, p_, branch_] := Module[{vals},
  vals = -VP[#, p, T, branch] & /@ T["derivs"];
  If[T["q"] =!= None, vals = Join[vals, -VP[#, p, T, branch] & /@ Rest[T["Dw"]]]];
  Max[Prepend[Select[vals, NumericQ[#] && # =!= -Infinity &], 0]]];

(* ClassifyPrime[T, p]  -- see the usage message.
   The uniformiser used in the normality test is pi = p at an unramified
   prime and pi = y at a branch prime (where v_P(y) = 1); P is normal iff
   v_P(D pi) = -eta_P, i.e. eta_P + v_P(D pi) = 0, and special iff the sum is
   >= 1 (D pi vanishes to higher order than the denominator divisor allows). *)
ClassifyPrime[T_, p_] := Module[{branch, eta, pi, special},
  branch = T["q"] =!= None && Vp[T["q"], p, T["gens"]] > 0;
  eta = Eta[T, p, branch];
  pi = If[branch, TUnit[T, If[T["m"] == 2, 1, PowerMod[BranchMult[T, p], -1, T["m"]]]], TScalar[T, p]];   (* uniformiser: w_{i*} (y for m = 2) or p *)
  special = eta + VP[TowerD[T, pi], p, T, branch] >= 1;
  {branch, eta, 1 + eta, special}];

(* ------------------------------------------------------------ residues *)
(* places over p: {g, rho, yval} with yval = None for m = 1 *)
(* PointsOver[T, p]
   The places of the curve over the prime p, for the residue computation.
   Looks for a generator g in which p is a polynomial of degree <= 4 with
   coefficients free of g (primes of higher degree go through the residue
   field: ResidueClasses), solves p = 0 for g (Solve may return Root objects), and
   for each root rho returns the two sheets {g, rho, +r} and {g, rho, -r}
   with r = Sqrt[q(rho)] simplified by PowerExpand (the sign convention of
   SymPy's positive symbols).  For m = 1 the third entry is None.  Returns
   None if p has no such generator direction.                               *)
PointsOver[T_, p_] := Module[{gens = T["gens"], g, cl, roots, qbar, r, pts},
 Catch[
  Do[g = gens[[i]];
    If[PolynomialQ[p, g] && 1 <= Exponent[p, g] <= 4 &&
       FreeQ[CoefficientList[p, g], g] &&
       (Exponent[p, g] <= 2 || FreeQ[CoefficientList[p, g], Alternatives @@ gens]),
      (* degree 3-4: only with constant coefficients (roots in radicals); a
         cubic/quartic prime with generator-dependent coefficients goes through
         the residue field instead (ResidueClasses), as in parallel_mixed.py *)
      roots = g /. Solve[p == 0, g];
      pts = {};
      Do[If[T["q"] === None, AppendTo[pts, {g, rho, None}],
          qbar = Can[T["q"] /. g -> rho];
          r = MRoot[qbar, T["m"]];
          If[NumericQ[r], r = RootReduce[r]];
          If[T["m"] == 2,
            AppendTo[pts, {g, rho, r}]; AppendTo[pts, {g, rho, -r}],
            Do[AppendTo[pts, {g, rho, If[k == 0, r, RRad[RootOfUnity[T["m"]]^k r]]}], {k, 0, T["m"] - 1}]]],
        {rho, roots}];
      Throw[pts, "pts"]],
    {i, Length[gens]}];
  None, "pts"]];

(* ReduceAt[expr, pt, Y]
   Evaluates an expression in the generators and the symbol Y (standing for
   y) at the place pt = {g, rho, yval}: substitutes Y -> yval, cancels, then
   substitutes g -> rho and simplifies.  The result is an element of the
   residue field kappa(P), expressed in the remaining generators.           *)
ReduceAt[expr_, pt_, Y_] := Module[{g, rho, yv, e},
  {g, rho, yv} = pt; e = expr;
  If[yv =!= None, e = e /. Y -> yv];
  Simplify[Can[e] /. g -> rho]];

(* ResFmt[v]
   Presentation of a residue: Simplify, and RootReduce when the value is a
   pure algebraic number (RootReduce is not applied to expressions containing
   the generators).                                                          *)
ResFmt[v_] := With[{s = Simplify[v]}, If[NumericQ[s], RootReduce[s], s]];

(* exact decision of tau \[Element] Fbar inside kappa(P): canonical form a + b Sqrt[A] *)
(* CertifyNonconstant[tau, gens]  -- see the usage message.
   Implementation: collect the powers with generator-dependent bases; refuse
   (Undecided) any exponent other than an integer or a half-integer, or more
   than one distinct square-root base A.  With no algebraic function the test
   is on the numerator and denominator of Can[tau].  With one, substitute
   Sqrt[A] -> S, reduce numerator and denominator modulo S^2 - A, multiply by
   the conjugate of the denominator so that D0 is free of S, and read off
   a + b S.  If A is itself a perfect square in F(gens) the expression was
   rational after all and the test recurses on (a + b r)/D0.               *)
CertifyNonconstant[tau0_, gens_] := Module[{tau, G, pows, algs, A, S, e, n, d, red, nn, dd, dbar, NN, D0, a, b, r, isSq},
  G = Alternatives @@ gens;
  tau = Together[tau0];
  If[FreeQ[tau, G], Return[False]];
  pows = Cases[{tau}, Power[base_, ex_] /; ! FreeQ[base, G], Infinity];
  If[AnyTrue[pows, ! (IntegerQ[#[[2]]] || Denominator[#[[2]]] == 2) &], Return[Undecided]];
  algs = DeleteDuplicates[Cases[pows, Power[base_, ex_] /; Denominator[ex] == 2 :> base]];
  If[algs === {},
    e = Can[tau];
    Return[! FreeQ[Numerator[e], G] || ! FreeQ[Denominator[e], G]]];
  If[Length[algs] =!= 1, Return[Undecided]];
  A = First[algs]; S = Unique["S"];
  e = tau /. Sqrt[A] -> S /. Power[A, -1/2] -> 1/S;
  e = Together[e]; {n, d} = {Numerator[e], Denominator[e]};
  red[z_] := PolynomialRemainder[Expand[z], S^2 - A, S];
  {nn, dd} = {red[n], red[d]};
  dbar = dd /. S -> -S;
  NN = red[Expand[nn dbar]]; D0 = red[Expand[dd dbar]];
  If[! FreeQ[D0, S], Return[Undecided]];
  a = NN /. S -> 0; b = Coefficient[Expand[NN], S, 1];
  r = PowerExpand[Sqrt[Factor[A]]];
  isSq = FreeQ[r, Power[base_, ex_] /; Denominator[ex] == 2 && ! FreeQ[base, G]];
  If[isSq, Return[CertifyNonconstant[Cancel[(a + b r)/D0], gens]]];
  If[Can[b/D0] =!= 0, Return[True]];
  CertifyNonconstant[Cancel[a/D0], gens]];

(* ------------------------------------------------------------ the algorithm *)


(* ============================================================ stage 2 ===== *)

(* IsZero[v]
   Exact zero test.  For algebraic constants RootReduce canonicalises
   (Sqrt[-I] - (1 - I)/Sqrt[2] -> 0, which Simplify alone misses); for
   expressions with free symbols Simplify is used.                           *)
IsZero[v_] := With[{w = Together[v]}, If[NumericQ[w], RootReduce[w] === 0, Simplify[w] === 0]];

(* RR[e]
   RootReduce for algebraic numbers, Together otherwise: the canonical form
   used for the coordinates in the group law and the Miller loop.          *)
RR[e_] := If[NumericQ[e], RootReduce[e], Together[e]];
(* RRad[e]: an algebraic number in radicals when it has them (ToRadicals,
   then ComplexExpand so that (-1)^(1/4) becomes (1 + I)/Sqrt[2]); Root
   objects of degree > 4 are kept as RootReduce leaves them.  Used for the
   residue classes and their logands, whose Root-object forms make the
   later Cancel/Together over the extension field intractable.            *)
RRad[e_] := If[NumericQ[e], With[{r = Quiet[ComplexExpand[ToRadicals[RootReduce[e]]]]}, If[FreeQ[r, Root], r, RootReduce[e]]], Together[e]];

(* NormSearchAll[q, p, g, kmax, dbmax]
   Every u = a + b y with a^2 - q b^2 = c p^k of small degree, b monic of
   degree db <= dbmax, 1 <= k <= kmax (Algorithm 3(c)).  Returns a list of
   {a, b, c, k}.  Different solutions vanish on different subsets of the
   places over p, so the realisation tries them all.                        *)
NormSearchAll[q_, p_, g_, kmax_: 2, dbmax_: 2] := Module[{dq, dp, out = {}, da, ac, bc, c, a, b, eq, sols},
  dq = Exponent[q, g]; dp = Exponent[p, g];
  Do[da = Quotient[Max[dq + 2 db, k dp] + 1, 2];
    ac = Table[na[i], {i, 0, da}]; bc = Table[nb[i], {i, 0, db - 1}]; c = nc;
    a = ac . g^Range[0, da];
    b = If[db == 0, 1, bc . g^Range[0, db - 1] + g^db];
    eq = Expand[a^2 - q b^2 - c p^k];
    sols = Quiet[Solve[Thread[CoefficientList[eq, g] == 0], Join[ac, bc, {c}]]];
    Do[With[{aa = Collect[Expand[a /. s], g, RR], bb = Collect[Expand[b /. s], g, RR], cc = RR[c /. s]},
        (* coefficients RootReduced: nested radicals from Solve would otherwise
           defeat the extension-field arithmetic of the linear system *)
        If[cc =!= 0 && FreeQ[aa, Alternatives @@ Join[ac, bc, {c}]] && NumericQ[cc],
          AppendTo[out, {aa, bb, cc, k}]]],
      {s, sols}],
    {k, 1, kmax}, {db, 0, dbmax}];
  out];

(* SheetSign[T, pt]
   For a place pt = {g, rho, yv} with constant coordinates, the sign eps
   with yv = eps Sqrt[q(rho)], decided exactly; None if neither sign matches.*)
SheetSign[T_, pt_] := Module[{g, rho, yv, r},
  {g, rho, yv} = pt; r = Sqrt[T["q"] /. g -> rho];
  Which[IsZero[yv - r], 1, IsZero[yv + r], -1, True, None]];

(* VanishOrder[T, u, pt, Y]
   Order of vanishing of the pair u = a + b y at a place with constant
   coordinates, by the series of u(rho + e, eps Sqrt[q(rho + e)]) in the
   base generator (1 if the place has non-constant coordinates).           *)
VanishOrder[T_, u_, pt_, Y_] := Module[{g, rho, yv, eps, e, ub, k = 0},
  {g, rho, yv} = pt;
  If[! FreeQ[rho, Alternatives @@ T["gens"]] || ! FreeQ[T["q"], Alternatives @@ DeleteCases[T["gens"], g]], Return[1]];
  eps = SheetSign[T, pt];
  If[eps === None, Return[1]];
  ub = (u[[1]] + u[[2]] eps Sqrt[T["q"] /. g -> rho + e]) /. g -> rho + e;
  While[k < 8 && IsZero[SeriesCoefficient[ub, {e, 0, k}]], k++];
  k];

(* DeepResidues[T, f, p, pts, Y]
   Laurent residues at a delta = 1 normal prime with a pole of order >= 2
   (Proposition 7.6 in the case where D acts trivially on the coefficient
   field: n = 1, D = d/dx).  The coefficient of e^-1 of f(rho + e, eps
   Sqrt[q(rho + e)]) at each place.                                         *)
DeepResidues[T_, f_, p_, pts_, Y_] := Module[{g = T["gens"][[1]], e, out = {}, eps, floc, m = T["m"], r, yser},
  Do[If[m == 2,
      eps = SheetSign[T, pt];
      If[eps === None, AppendTo[out, Undecided]; Continue[]];
      floc = (f[[1]] + f[[2]] eps Sqrt[T["q"] /. g -> pt[[2]] + e]) /. g -> pt[[2]] + e,
      (* m >= 3: the sheet y = eps q(rho)^(1/m), eps a root of unity *)
      r = MRoot[T["q"] /. g -> pt[[2]], m];
      eps = RRad[pt[[3]]/r];
      yser = eps (T["q"] /. g -> pt[[2]] + e)^(1/m);
      floc = Sum[(f[[i + 1]] /. g -> pt[[2]] + e) yser^i/(T["E"][[i + 1]] /. g -> pt[[2]] + e), {i, 0, T["n"] - 1}]];
    AppendTo[out, ResFmt[SeriesCoefficient[floc, {e, 0, -1}]]],
    {pt, pts}];
  out];

(* EllOps[q, g]
   Chord-and-tangent group law on y^2 = q(g), deg q = 3, after normalising
   to the monic model X = c3 g, Ytil = c3 y: Ytil^2 = X^3 + c2 X^2 + c1 c3 X
   + c0 c3^2.  Returns {toM, add, c}: toM maps a point (g, y) to the model,
   add is the group law with the identity represented by None, and c the
   coefficients {c3, c2, c1, c0} of q.                                       *)
EllOps[q_, g_] := Module[{c, c3, toM, add},
  c = PadLeft[Reverse[CoefficientList[q, g]], 4];       (* {c3, c2, c1, c0} *)
  c3 = c[[1]];
  toM = Function[P, {RR[c3 P[[1]]], RR[c3 P[[2]]]}];
  (* no Return here: inside a pure Function it would survive unevaluated
     as Return[...] instead of returning, so the identity None must be
     produced by the If/Which structure itself *)
  add = Function[{P, Q}, Which[
    P === None, Q,
    Q === None, P,
    True, Module[{x1, y1, x2, y2, lam, x3},
      {x1, y1} = P; {x2, y2} = Q;
      If[IsZero[x1 - x2] && IsZero[y1 + y2], None,
        lam = If[IsZero[x1 - x2],
          RR[(3 x1^2 + 2 c[[2]] x1 + c[[3]] c3)/(2 y1)],
          RR[(y2 - y1)/(x2 - x1)]];
        x3 = RR[lam^2 - c[[2]] - x1 - x2];
        {x3, RR[lam (x1 - x3) - y1]}]]]];
  {toM, add, c}];

(* DivisionPolyOrder[A, B, X0, Y0, Nmax]
   Order of the point (X0, Y0) on Ytil^2 = X^3 + A X + B if it is <= Nmax,
   else None: the division polynomials psi_n are built by the standard
   recursion with Ytil^2 reduced, and n P = O iff psi_n(P) = 0.            *)
DivisionPolyOrder[A_, B_, X0_, Y0_, Nmax_: 24] := Module[{X, yv, psi, red, get, n},
  red[e_] := PolynomialRemainder[Expand[e], yv^2 - (X^3 + A X + B), yv];
  psi[0] = 0; psi[1] = 1; psi[2] = 2 yv;
  psi[3] = 3 X^4 + 6 A X^2 + 12 B X - A^2;
  psi[4] = 4 yv (X^6 + 5 A X^4 + 20 B X^3 - 5 A^2 X^2 - 4 A B X - 8 B^2 - A^3);
  get[m_] := get[m] = With[{h = Quotient[m, 2]},
    If[OddQ[m],
      red[get[h + 2] get[h]^3 - get[h - 1] get[h + 1]^3],
      red[Cancel[get[h] (get[h + 2] get[h - 1]^2 - get[h - 2] get[h + 1]^2)/(2 yv)]]]];
  Do[get[k] = psi[k], {k, 0, 4}];
  n = SelectFirst[Range[2, Nmax], IsZero[get[#] /. {X -> X0, yv -> Y0}] &, None];
  n];

(* TorsionOrder[q, g, P, add, c, bound]
   The order m <= bound of the point P (model coordinates of EllOps): by
   the division polynomials when the model is depressed, otherwise by
   repeated addition; None if no order <= bound.                            *)
TorsionOrder[q_, g_, P_, add_, c_, bound_] := Module[{m = None, kP},
  If[c[[2]] === 0, m = DivisionPolyOrder[c[[3]] c[[1]], c[[4]] c[[1]]^2, P[[1]], P[[2]], bound]];
  If[m === None, kP = P; Do[kP = add[kP, P]; If[kP === None, m = k; Break[]], {k, 2, bound}]];
  m];

(* EllSum[q, g, Y, terms, add, c]
   The sum S of the points of terms = {{P, n}, ...}, n an integer, in the
   group law of EllOps, and the function h with div(h) = sum n (P - oo) -
   (S - oo): for S + P = R the line through S and P over the vertical
   through R has divisor S + P - R - oo, and -P = (x_P, -y_P) with
   div(x - x_P) = P + (-P) - 2 oo.  Accumulated in K[Y, g] modulo Y^2 - q;
   returns {S, {h0, h1}} with h = h0 + h1 y.  With terms = {{P, m}}, m the
   order of P, h is the Miller function with divisor m P - m oo.           *)
EllSum[q_, g_, Y_, terms_, add_, c_] := Module[{c3 = c[[1]], X0, S = None, num = 1, den = 1, P, n, x1, y1, lam, h0, h1},
  X0 = c3 g;
  Do[{P, n} = term;
    If[n < 0, P = {P[[1]], -P[[2]]}; n = -n; den *= (X0 - P[[1]])^n];
    Do[If[S === None, S = P; Continue[]];
      {x1, y1} = S;
      If[IsZero[x1 - P[[1]]] && IsZero[y1 + P[[2]]], num *= (X0 - x1); S = None; Continue[]];   (* S = -P: the vertical line *)
      lam = If[IsZero[x1 - P[[1]]], RR[(3 x1^2 + 2 c[[2]] x1 + c[[3]] c3)/(2 y1)], RR[(P[[2]] - y1)/(P[[1]] - x1)]];
      num = PolynomialRemainder[Expand[num (c3 Y - y1 - lam (X0 - x1))], Y^2 - q, Y];
      S = add[S, P];
      den *= (X0 - S[[1]]),
      {n}],
    {term, terms}];
  h0 = Can[(num + (num /. Y -> -Y))/2/den];
  h1 = Can[(num - (num /. Y -> -Y))/(2 Y)/den];
  {S, Collect[Numerator[#], g, RR]/Collect[Denominator[#], g, RR] & /@ {h0, h1}}];

(* TorsionOrderAndMiller[q, g, pt, Y, bound]
   On y^2 = q(g), deg q = 3: the order m of [P - oo] for the place
   pt = {g, rho, y0} with constant coordinates and the Miller function h
   with div(h) = m P - m oo, as a pair; None if no order <= bound is found.
   Used to generate S'-units over a special prime.                          *)
TorsionOrderAndMiller[q_, g_, pt_, Y_, bound_: 24] := Module[{toM, add, c, P, m},
  {toM, add, c} = EllOps[q, g];
  P = toM[{pt[[2]], pt[[3]]}];
  m = TorsionOrder[q, g, P, add, c, bound];
  If[m === None, None, {m, EllSum[q, g, Y, {{P, m}}, add, c][[2]]}]];

(* QCoords[consts]
   The coordinates over Q of algebraic constants in one common number
   field (ToNumberField), padded to a common length; a rational constant
   is its own coordinate.                                                   *)
QCoords[consts_List] := PadRows[If[Head[#] === AlgebraicNumber, #[[2]], {#}] & /@ ToNumberField[consts]];

(* PadRows[rows]: the rows padded with zeros to a common length. The one-argument
   PadRight now rectangularises a ragged list correctly (was A4, v0.170), so this
   delegates to it. *)
PadRows[rows_List] := PadRight[rows];

(* QBasisDivisors[taus]
   The residues taus (algebraic constants) over a Q-basis of their Q-span
   (Algorithm 3(d)): {{beta_j, ns_j}, ...} with ns_j integer vectors and
   tau_i = sum_j ns_j[[i]] beta_j, so that the divisor sum_i tau_i P_i is
   sum_j beta_j D_j with D_j = sum_i ns_j[[i]] P_i; beta_j is the residue at
   the j-th pivot divided by the common denominator of its row.  None if a
   coordinate is not rational.                                              *)
QBasisDivisors[taus_List] := Module[{A, Rr, piv, N0},
  A = Transpose[QCoords[taus]];                            (* column i: the residue tau_i *)
  If[! MatrixQ[A, MatchQ[#, _Integer | _Rational] &], Return[None]];
  Rr = DeleteCases[RowReduce[A], {0 ..}];                  (* A[[All, i]] = sum_j Rr[[j, i]] A[[All, piv[[j]]]] *)
  piv = FirstPosition[#, z_ /; z != 0][[1]] & /@ Rr;
  Table[N0 = LCM @@ Denominator[Rr[[j]]]; {taus[[piv[[j]]]]/N0, N0 Rr[[j]]}, {j, Length[piv]}]];

(* TorsionRealise[T, pending, Y, bound, verbose]
   Algorithm 3(d): the residue divisor of the pending primes {{p, pts,
   taus}, ...} on an elliptic curve with one place at infinity (deg q = 3,
   constant coordinates).  A place whose class [P - oo] is torsion of
   order mu is realised on its own by its Miller function; the residues of
   the remaining places are decomposed over a Q-basis (beta_j) of their
   Q-span, tau_P = sum_j n_Pj beta_j with n_Pj integers, so that the
   divisor is sum_j beta_j D_j with D_j = sum_P n_Pj (P - oo) of degree 0,
   and each D_j is realised by the group law: S_j = sum_P n_Pj P, mu_j the
   order of S_j (1 when S_j = oo), and EllSum on mu_j D_j gives u_j with
   div(u_j) = mu_j D_j, recorded with coefficient beta_j/mu_j.  Returns a
   list of {coefficient, pair} or None (a class of order > bound, no cubic
   model, or a place with non-constant coordinates).                        *)
TorsionRealise[T_, pending_, Y_, bound_: 24, verbose_: False] := Catch[Module[
  {gens = T["gens"], g, q = T["q"], places, toM, add, c, realise, out = {}, rest = {}, P, got, qb, ns, terms, beta},
  g = pending[[1, 2, 1, 1]];
  If[Exponent[q, g] =!= 3 || ! FreeQ[q, Alternatives @@ DeleteCases[gens, g]] || AnyTrue[pending, #[[2, 1, 1]] =!= g &], Throw[None, "tors"]];
  places = Cases[Flatten[Thread[{#[[2]], #[[3]]}] & /@ pending, 1], {_, tau_} /; tau =!= 0];
  If[! FreeQ[{#[[1, 2]], #[[1, 3]], #[[2]]} & /@ places, Alternatives @@ gens], Throw[None, "tors"]];
  {toM, add, c} = EllOps[q, g];
  (* {mu, u} with div(u) = mu * sum n (P - oo) for terms = {{P, n}, ...}, or None *)
  realise[terms0_] := Module[{S = None, mu},
    Do[Do[S = add[S, If[tm[[2]] > 0, tm[[1]], {tm[[1, 1]], -tm[[1, 2]]}]], {Abs[tm[[2]]]}], {tm, terms0}];
    mu = If[S === None, 1, TorsionOrder[q, g, S, add, c, bound]];
    If[mu === None, None, {mu, EllSum[q, g, Y, {#[[1]], mu #[[2]]} & /@ terms0, add, c][[2]]}]];
  Do[P = toM[{pl[[1, 2]], pl[[1, 3]]}];
    got = realise[{{P, 1}}];
    If[got === None, AppendTo[rest, {pl[[1]], P, pl[[2]]}],
      AppendTo[out, {RR[pl[[2]]/got[[1]]], got[[2]]}];
      If[verbose, Print["      torsion: [P - oo] of order ", got[[1]], " at (", pl[[1, 2]], ", ", pl[[1, 3]], "); Miller logand with coefficient ", RR[pl[[2]]/got[[1]]]]]],
    {pl, places}];
  If[rest =!= {},
    qb = QBasisDivisors[rest[[All, 3]]];
    If[qb === None, Throw[None, "tors"]];
    Do[{beta, ns} = bn;
      terms = Cases[Thread[{rest[[All, 2]], ns}], {_, n_} /; n != 0];
      got = realise[terms];
      If[got === None, Throw[None, "tors"]];
      AppendTo[out, {RR[beta/got[[1]]], got[[2]]}];
      If[verbose, Print["      torsion: the divisor ", StringRiffle[Cases[Transpose[{rest[[All, 1]], ns}], {pt_, n_} /; n != 0 :> ToString[n] <> "*(" <> ToString[pt[[2]], InputForm] <> ", " <> ToString[pt[[3]], InputForm] <> ")"], " + "],
        " - (", Total[ns], ") oo has order ", got[[1]], "; logand with coefficient ", RR[beta/got[[1]]]]],
      {bn, qb}]];
  out], "tors"];

(* NontorsionDivisor[T, pending, verbose]
   Proposition 9.4 for the residue divisor that Algorithm 3(a)-(d) leaves
   unrealised (pending = {{p, pts, taus}, ...} on y^2 = q with constant
   coordinates): each component over a Q-basis of the residues
   (QBasisDivisors) is reduced modulo good primes and the orders of its
   class in the Jacobians compared (NontorsionDivisorCertificate).  Returns
   {data, divisor} for the first component certified non-torsion, the
   integral then being non-elementary by Corollary 7.6; else None.       *)
NontorsionDivisor[T_, pending_, verbose_: False] := Catch[Module[{gens = T["gens"], g, q = T["q"], places, consts, nf, an, mu, coords, qb, pl, certd, dv, z},
  g = pending[[1, 2, 1, 1]];
  If[! FreeQ[q, Alternatives @@ DeleteCases[gens, g]] || AnyTrue[pending, #[[2, 1, 1]] =!= g &], Return[None]];
  places = Cases[Flatten[Thread[{#[[2]], #[[3]]}] & /@ pending, 1], {_, tau_} /; tau =!= 0];
  consts = Flatten[{#[[1, 2]], #[[1, 3]], #[[2]]} & /@ places];
  If[! FreeQ[consts, Alternatives @@ gens] || ! FreeQ[CoefficientList[q, g], Alternatives @@ gens], Return[None]];
  (* one number field Q(theta) for the coordinates and the residues: ONE
     ToNumberField call, whose theta the coordinates and mu share *)
  nf = ToNumberField[consts];
  an = FirstCase[nf, _AlgebraicNumber, None];
  mu = If[an === None, None, CoefficientList[MinimalPolynomial[an[[1]], z], z]];
  coords = PadRows[If[Head[#] === AlgebraicNumber, #[[2]], {#}] & /@ nf];
  If[! MatrixQ[coords, MatchQ[#, _Integer | _Rational] &], Return[None]];
  qb = QBasisDivisors[places[[All, 2]]];
  If[qb === None, Return[None]];
  Do[pl = Table[If[bn[[2, i]] == 0, Nothing, {{coords[[3 i - 2]], coords[[3 i - 1]]}, bn[[2, i]]}], {i, Length[places]}];
    certd = NontorsionDivisorCertificate[q, g, pl, mu];
    If[certd[[1]],
      dv = StringRiffle[Table[If[bn[[2, i]] == 0, Nothing, ToString[bn[[2, i]]] <> "*(" <> ToString[places[[i, 1, 2]], InputForm] <> ", " <> ToString[places[[i, 1, 3]], InputForm] <> ")"], {i, Length[places]}], " + "];
      If[verbose, Print["  the residue divisor ", dv, " (coefficient ", bn[[1]], ") is certified non-torsion by reduction mod p: ", certd[[2]]]];
      Throw[{certd[[2]], dv}, "ntdiv"]],
    {bn, qb}];
  None], "ntdiv"];

(* PairReduce[e, Ysym, q]
   The pair {a, b} = a + b y of an expression rational in the generators
   and Ysym with Ysym^2 = q (denominator rationalised); the general form of
   ToPair used by the change of model.                                     *)
PairReduce[e_, Ysym_, q_] := ToPair[e, Ysym, q];

(* QuarticToCubic[T, f]
   Change of model for a quartic radicand q(g) with a square leading
   coefficient (two rational places at infinity).  With p the polynomial
   part of Sqrt[q] and r = q - p^2 (deg r <= 2), the substitution
   s = y - p(g) gives 2 p(g) s = r(g) - s^2, a quadratic A' g^2 + B' g + C'
   over Q(s) whose discriminant Delta(s) is a cubic; with T = 2 A' g + B',
   T^2 = Delta(s), and (g, y) -> (s, T) is birational.  The whole tower is
   rewritten in the new generators (s replaces g at its position).  Returns
   {T', f', back}, back being the rules s -> y - p(g), Sqrt[Delta] ->
   2 A' g + B' in the old generators, applied to the result.  On the cubic
   model the torsion machinery (Miller functions) is available, which is
   what realises residue divisors and S'-units of large order on quartics. *)
QuarticToCubic[T_, f_] := Module[{gens = T["gens"], q = T["q"], g, p, r, s, TT, Y, pc, rc, Ap, Bp, Cp, Delta, gNew, yNew, conv,
    sNew, Dg, Ds, newDerivs, newGens, T2, f2, back},
  g = SelectFirst[gens, FreeQ[q, Alternatives @@ DeleteCases[gens, #]] &];
  p = SqrtPolyPart[q, g]; r = Expand[q - p^2];
  s = Unique["s"]; TT = Unique["T"]; Y = Unique["Y"];
  pc = PadRight[CoefficientList[p, g], 3]; rc = PadRight[CoefficientList[r, g], 3];
  Ap = Expand[2 pc[[3]] s - rc[[3]]]; Bp = Expand[2 pc[[2]] s - rc[[2]]]; Cp = Expand[2 pc[[1]] s - rc[[1]] + s^2];
  Delta = Expand[Bp^2 - 4 Ap Cp];
  gNew = Can[(TT - Bp)/(2 Ap)]; yNew = s + (p /. g -> gNew);
  conv[pair_] := PairReduce[(pair[[1]] + pair[[2]] Y) /. {g -> gNew, Y -> yNew}, TT, Delta];
  sNew = Unique["s"];
  Dg = T["derivs"][[Position[gens, g][[1, 1]]]];
  Ds = Padd[T["Dy"], Pscale[-D[p, g], Dg]];
  newDerivs = MapThread[If[#1 === g, conv[Ds], conv[#2]] &, {gens, T["derivs"]}] /. s -> sNew;
  newGens = gens /. g -> sNew;
  T2 = Tower[newGens, newDerivs, Delta /. s -> sNew];
  f2 = conv[f] /. s -> sNew;
  back = {sNew -> Sqrt[q] - p, Sqrt[Delta /. s -> sNew] -> (2 Ap g + Bp /. s -> sNew)};
  {T2, f2, back}];

(* RationalPoint[q, g]
   A rational point {g0, y0} on the conic y^2 = q(g), searching small
   integers and fractions; None if none is found in that range.            *)
RationalPoint[q_, g_] := Module[{v, y0},
  Do[v = q /. g -> g0; y0 = Sqrt[v];
    If[v >= 0 && MatchQ[y0, _Integer | _Rational], Return[{g0, y0}, Module]],
    {g0, DeleteDuplicates[Flatten[Table[{num/den, -num/den}, {num, 0, 12}, {den, 1, 3}]]]}];
  None];

(* ConicToLine[T, f]
   Rational parametrisation of a conic y^2 = q(g), deg q = 2, through a
   rational point (g0, y0): the pencil y = y0 + w (g - g0) meets the conic
   in the second point
     g = (2 y0 w - (q2 + w^2) g0 - q1)/(q2 - w^2),   y = y0 + w (g - g0),
   so the curve disappears and the whole tower is rewritten over Q(w) (w
   replaces g at its position; q becomes None).  Returns {T', f', back} with
   back the rule w -> (y - y0)/(g - g0), or None without a rational point.
   On a genus-0 curve every S'-unit and residue divisor is then a question
   about rational functions of w, which the flattened-tower machinery
   answers without any extension-field arithmetic over the curve.         *)
ConicToLine[T_, f_] := Module[{gens = T["gens"], q = T["q"], g, q2, q1, q0, pt, g0, y0, w, Y, gNew, yNew, conv, Dg, num, Dw, newDerivs, newGens, T2, f2},
  g = SelectFirst[gens, FreeQ[q, Alternatives @@ DeleteCases[gens, #]] &];
  {q0, q1, q2} = PadRight[CoefficientList[q, g], 3];
  pt = RationalPoint[q, g];
  If[pt === None, Return[None]];
  {g0, y0} = pt;
  w = Unique["w"]; Y = Unique["Y"];
  gNew = Cancel[(2 y0 w - (q2 + w^2) g0 - q1)/(q2 - w^2)];
  yNew = Cancel[y0 + w (gNew - g0)];
  conv[pair_] := {Cancel[(pair[[1]] + pair[[2]] Y) /. {g -> gNew, Y -> yNew}], 0};
  Dg = T["derivs"][[Position[gens, g][[1, 1]]]];
  (* D w = (Dy (g - g0) - (y - y0) Dg)/(g - g0)^2 as a pair over the conic *)
  num = Padd[Pscale[g - g0, T["Dy"]], Pscale[-1, Pmul[{-y0, 1}, Dg, q]]];
  Dw = Pscale[1/(g - g0)^2, num];
  newDerivs = MapThread[If[#1 === g, conv[Dw], conv[#2]] &, {gens, T["derivs"]}];
  newGens = gens /. g -> w;
  T2 = Tower[newGens, newDerivs, None];
  f2 = conv[f];
  {T2, f2, {w -> (Sqrt[q] - y0)/(g - g0)}}];

(* ------------------------------------------ residues in the residue field *)
(* ConstDir[T, p]
   The generator g in which p is a polynomial with coefficients free of every
   generator and, when there is a curve, in which q is a polynomial as well;
   None otherwise.                                                          *)
ConstDir[T_, p_] := Catch[Module[{gens = T["gens"], G},
  G = Alternatives @@ gens;
  Do[If[PolynomialQ[p, g] && Exponent[p, g] >= 1 && FreeQ[CoefficientList[p, g], G],
      Throw[If[T["q"] === None || FreeQ[T["q"], Alternatives @@ DeleteCases[gens, g]], g, None], "cd"]],
    {g, gens}];
  None], "cd"];

(* ResidueDir[T, p]
   A generator g in which p has positive degree, for the residue-field
   computation.  Prefers a constant-coefficient direction (ConstDir); for a
   transcendental tower (q = None) any direction with Exponent[p, g] >= 1 will
   do, the coefficients then involving the other generators, so that
   kappa(P) = Frac(F[gens]/(p)) is a field over F(other generators) -- the
   setting of Bronstein's parallel Risch, where every normal factor of the
   denominator carries a residue (Theorem 7.5).                              *)
ResidueDir[T_, p_] := Module[{g = ConstDir[T, p]},
  If[g =!= None, g,
    If[T["q"] === None, SelectFirst[T["gens"], PolynomialQ[p, #] && Exponent[p, #] >= 1 &, None], None]]];

(* CanonicalResidueField[T, f, p, g, delta, verbose]
   Algorithm 2 (the canonical residue, Proposition 7.7) at a normal prime P
   over the irreducible p of a transcendental tower (q = None; f the pair with
   scalar first coordinate), computed in the residue field
   kappa(P) = Frac(F[gens]/(p)) by reduction modulo p in the direction g --
   with NO substitution of a root, so the place may have non-constant
   coordinates and delta may exceed 1.  Runs the local Hermite reduction
   f <- f - D(c p^-k), c = -mu/(k lambda_P), lambda_P the class of
   p^(delta-1) D p, until the pole reaches order delta, then returns
   tau = (f p / Dp)|_P in kappa(P) (an expression reduced modulo p), or None
   when a denominator met on the way is not coprime to p.                     *)
CanonicalResidueField[T_, f_, p_, g_, delta_, verbose_: False] := Catch[Module[
  {gens = T["gens"], fc = Can[f[[1]]], Dp, lam, k, mu, c, i},
  Dp = TowerD[T, TScalar[T, p]][[1]];
  lam = ModP[Can[p^(delta - 1) Dp], p, g];
  If[lam === None || IsZero[lam], Throw[None, "crf"]];
  Do[
    k = -Vp[fc, p, gens] - delta;
    If[k < 0, Throw[0, "crf"]];
    If[k == 0, Throw[ModP[Can[fc p/Dp], p, g], "crf"]];
    mu = ModP[Can[p^(k + delta) fc], p, g];
    If[mu === None, Throw[None, "crf"]];
    c = ModP[Can[-mu/(k lam)], p, g];
    If[c === None, Throw[None, "crf"]];
    fc = Can[fc - TowerD[T, TScalar[T, Can[c/p^k]]][[1]]];
    If[verbose, Print["      Hermite step at (", p, "): order ", k + delta,
        " -> subtract D((", c, ") (", p, ")^(-", k, "))"]],
    {i, 16}];
  None], "crf"];

(* PolyZeroQ[P, g]: every coefficient of the polynomial P in g is zero (exact
   test on algebraic numbers by RootReduce, Together otherwise).            *)
PolyZeroQ[P_, g_] := AllTrue[CoefficientList[Expand[P], g], IsZero];

(* InvMod[d, p, g]
   The inverse of the polynomial d modulo p by the extended Euclidean
   algorithm, verified exactly; a linear solve when the Euclidean route does
   not return a certified inverse; None if d and p are not coprime.        *)
InvMod[d_, p_, g_] := Catch[Module[{eg, inv, n = Exponent[p, g], sv, s, sol},
  eg = Quiet[PolynomialExtendedGCD[d, p, g]];
  If[ListQ[eg] && FreeQ[eg[[1]], g] && ! IsZero[eg[[1]]],
    inv = Collect[Expand[eg[[2, 1]]/eg[[1]]], g, RR];
    If[PolyZeroQ[PolynomialRemainder[Expand[inv d], p, g] - 1, g], Throw[inv, "inv"]]];
  sv = Table[Unique["s"], {n}]; s = sv . g^Range[0, n - 1];
  sol = Quiet[Solve[Thread[CoefficientList[PolynomialRemainder[Expand[d s], p, g] - 1, g] == 0], sv]];
  If[sol === {}, Throw[None, "inv"]];
  Collect[Expand[s /. First[sol]], g, RR]], "inv"];

(* ModP[expr, p, g]
   A rational function of g whose denominator is coprime to p, reduced modulo
   p: a polynomial in g of degree < deg p.  Other generators and algebraic
   constants may occur in the coefficients.  None if not reducible.         *)
ModP[expr_, p_, g_] := Module[{e = Together[expr], nr, dr, di},
  nr = PolynomialRemainder[Numerator[e], p, g]; dr = PolynomialRemainder[Denominator[e], p, g];
  di = InvMod[dr, p, g];
  If[di === None, None, Collect[Expand[PolynomialRemainder[Expand[nr di], p, g]], g, Can]]];

Monic[P_, g_] := Collect[Expand[P/Coefficient[P, g, Exponent[P, g]]], g, RRad];

(* InFieldQ[a, gs]: is the algebraic number a in the number field Q(gs)?  By
   the degree of the field generated by gs with and without a.            *)
InFieldQ[a_, gs_] := Module[{g0 = Select[RR /@ gs, ! MatchQ[#, _Integer | _Rational] &], deg},
  deg[l_] := Length[First[ToNumberField[l]][[2]]];
  MatchQ[RR[a], _Integer | _Rational] || (g0 =!= {} && deg[Append[g0, a]] == deg[g0])];

(* ResidueClasses[T, p, tp, Y, verbose]
   Residues at the places over a prime p of ANY degree, computed in the
   residue field K[g]/(p) instead of at explicit roots (the residue-polynomial
   step of Rothstein and Trager inside Algorithm 2).  tp = {tau0, tau1} is
   the residue expression e (f h / Dh) as a pair; its classes modulo p are the
   residues tau0(alpha) +- tau1(alpha) y0 on the two sheets over each root
   alpha of p.  Their minimal polynomial over K is the residue polynomial
       R(z) = Res_g(p, (z - tau0)^2 - tau1^2 q)     (Res_g(p, z - tau0) without a curve);
   its irreducible factors give the residue values c and the fields K(c); the
   places with residue c lie over p_c = gcd(p, tau1^2 q - (c - tau0)^2) on
   the sheet y = Y0(g) := (c - tau0)/tau1 mod p_c (a polynomial with
   Y0^2 = q mod p_c); where tau1 vanishes both sheets carry c and the divisor
   is principal.  Returns {g, classes, principal} with classes = {{c, p_c, Y0}}
   and principal = {{c, p_c2}}; a status list {"not elementary", ...} or
   {"failed", ...}; or None when p has no usable direction.               *)
ResidueClasses[T_, p_, tp_, Y_, verbose_: False] := Catch[Module[
  {g, q = T["q"], m = T["m"], others, taus, tau0, tau1, z, R, t1q, facs, values = {}, classes = {}, principal = {},
   pc, EE, t1, p2, Y0, rts, c, MM, charp = None},
  g = ResidueDir[T, p];
  If[g === None, Throw[None, "rc"]];
  others = Alternatives @@ DeleteCases[T["gens"], g];
  taus = Table[ModP[tp[[i]], p, g], {i, If[q === None, 1, T["n"]]}];
  If[MemberQ[taus, None], Throw[{"failed", "residue not reducible modulo the prime", p}, "rc"]];
  tau0 = taus[[1]]; tau1 = If[q === None || m =!= 2, 0, taus[[2]]];
  If[Length[T["gens"]] > 1 && ! FreeQ[taus, others],
    Throw[{"not elementary", p, ToY[T, PadRight[taus, T["n"]], Y]}, "rc"]];          (* the residue depends on another generator *)
  z = Unique["z"];
  R = Which[q === None, Resultant[p, z - tau0, g],
        m == 2, t1q = ModP[Expand[tau1^2 q], p, g]; Resultant[p, Expand[(z - tau0)^2 - t1q], g],
        True,
        (* N(z - tau) = the characteristic polynomial of multiplication by tau *)
        MM = TMatrix[T, taus];
        charp = Expand[Det[z IdentityMatrix[m] - MM]];
        charp = Sum[ModP[Coefficient[charp, z, k], p, g] z^k, {k, 0, m}];
        Resultant[p, charp, g]];
  R = Expand[R];
  If[R === 0 || FreeQ[R, z], Throw[{"failed", "degenerate residue polynomial", p, R}, "rc"]];
  facs = Select[Rest[FactorList[R, Extension -> Automatic]], Exponent[#[[1]], z] >= 1 &];
  If[verbose, Print["      residue polynomial over (", p, "): ",
      StringRiffle[("(" <> ToString[#[[1]] /. z -> Global`z, InputForm] <> ")^" <> ToString[#[[2]]]) & /@ facs, " * "]]];
  Do[rts = z /. Solve[fac[[1]] == 0, z];
    Do[c = RRad[c0];
      If[IsZero[c], Continue[]];
      If[! AnyTrue[values, IsZero[# - c] &], AppendTo[values, c]],
      {c0, rts}],
    {fac, facs}];
  Do[
    If[q === None,
      pc = Monic[PolynomialGCD[p, tau0 - c, Extension -> Automatic], g];
      If[Exponent[pc, g] >= 1, AppendTo[principal, {c, pc}]];
      Continue[]];
    If[m >= 3,
      pc = Monic[PolynomialGCD[p, charp /. z -> c, Extension -> Automatic], g];
      If[Exponent[pc, g] < 1, Continue[]];
      Do[Which[sh[[2]] === None, Throw[{"failed", "residue class on several sheets over one root (m >= 3)", p, c}, "rc"],
               sh[[2]] === "all", AppendTo[principal, {c, sh[[1]]}],
               True, AppendTo[classes, {c, sh[[1]], sh[[2]]}]],
        {sh, SheetsOver[T, g, pc, taus, c]}];
      Continue[]];
    EE = Expand[tau1^2 q - (c - tau0)^2];
    pc = Monic[PolynomialGCD[p, EE, Extension -> Automatic], g];
    If[Exponent[pc, g] < 1, Continue[]];
    t1 = PolynomialRemainder[tau1, pc, g];
    p2 = Monic[PolynomialGCD[pc, t1, Extension -> Automatic], g];
    If[Exponent[p2, g] >= 1, AppendTo[principal, {c, p2}]; pc = Monic[PolynomialQuotient[pc, p2, g], g]];
    If[Exponent[pc, g] < 1, Continue[]];
    Y0 = Collect[Expand[PolynomialRemainder[Expand[(c - tau0) InvMod[PolynomialRemainder[t1, pc, g], pc, g]], pc, g]], g, RRad];
    If[! PolyZeroQ[PolynomialRemainder[Expand[Y0^2 - q], pc, g], g],
      Throw[{"failed", "internal: sheet polynomial does not square to q", p, c}, "rc"]];
    AppendTo[classes, {c, pc, Y0}],
    {c, values}];
  If[verbose,
    Do[Print["      residue class ", cl[[1]], ": places over ", cl[[2]], " on the sheet y = ", cl[[3]]], {cl, classes}];
    Do[Print["      residue ", pr[[1]], " on ", If[m == 2, "both", "all"], " sheets over ", pr[[2]], ": log(", pr[[2]], ")"], {pr, principal}]];
  {g, classes, principal}], "rc"];

(* PolyYStrip / PolyYRem / SheetsOver: polynomials in Y with coefficients in
   the residue ring K(c)[g]/(pc) (lists of polynomials in g, low degree first),
   and the sheets over pc on which the residue tau(y) = Sum taus_i w_i takes the
   value c: the gcd G(Y) of Y^m - q and c - tau(Y) over that ring.  Returns
   a list of {pc_j, sheet}, sheet the sheet polynomial Y0 (deg G = 1), "all"
   (deg G = m: every sheet, a principal divisor) or None (several but not all
   sheets, which the realisation does not represent).  A zero divisor in the
   Euclidean algorithm means pc is reducible over K(c): it is factored and
   the sheets computed over each factor. *)
PolyYStrip[A_, pc_, g_] := Module[{B = PolynomialRemainder[#, pc, g] & /@ A},
  While[B =!= {} && PolyZeroQ[Last[B], g], B = Most[B]]; B];
PolyYRem[A0_, B0_, pc_, g_] := Module[{A = PolyYStrip[A0, pc, g], B = PolyYStrip[B0, pc, g], inv, cf, sh},
  inv = InvMod[Last[B], pc, g];
  If[inv === None, Throw[$Failed, "zerodiv"]];
  While[Length[A] >= Length[B],
    cf = PolynomialRemainder[Expand[Last[A] inv], pc, g]; sh = Length[A] - Length[B];
    Do[A[[sh + i]] = PolynomialRemainder[Expand[A[[sh + i]] - cf B[[i]]], pc, g], {i, Length[B]}];
    A = PolyYStrip[A, pc, g]];
  A];
SheetsOver[T_, g_, pc_, taus_, c_] := Module[{m = T["m"], A, B, G, res, facs, Y0, inv},
  res = Catch[
    A = Join[{-T["q"]}, ConstantArray[0, m - 1], {1}];
    B = Table[With[{iv = InvMod[PolynomialRemainder[T["E"][[i + 1]], pc, g], pc, g]},
        If[iv === None, Throw[$Failed, "zerodiv"]];
        PolynomialRemainder[Expand[-taus[[i + 1]] iv], pc, g]], {i, 0, m - 1}];
    B[[1]] = PolynomialRemainder[Expand[B[[1]] + c], pc, g];
    A = PolyYStrip[A, pc, g]; B = PolyYStrip[B, pc, g];
    While[B =!= {}, {A, B} = {B, PolyYRem[A, B, pc, g]}];
    A, "zerodiv"];
  If[res === $Failed,
    facs = Select[Rest[FactorList[pc, Extension -> Automatic]], Exponent[#[[1]], g] >= 1 &][[All, 1]];
    If[Length[facs] <= 1, Return[{{pc, None}}]];
    Return[Join @@ (SheetsOver[T, g, Monic[#, g], taus, c] & /@ facs)]];
  G = res;
  Which[Length[G] - 1 == 1,
    inv = InvMod[G[[2]], pc, g];
    Y0 = Collect[PolynomialRemainder[Expand[-G[[1]] inv], pc, g], g, RRad];
    {{pc, Y0}},
    Length[G] - 1 == m, {{pc, "all"}},
    True, {{pc, None}}]];

(* RealiseClass[T, g, supp, c, mmax, verbose]
   Realisation of the residue class c supported on the places over the primes
   p_i of supp = {{p_i, Y0_i}}, sheet y = Y0_i(g): the logand u = a + b y with
   div(u) = m D_c - n_+ oo_+ - n_- oo_- for the least m <= mmax for which one
   exists, by LINEAR algebra over K(c).  Y0_i is lifted by Hensel's iteration
   to Y_i with Y_i^2 = q mod p_i^m, 'u vanishes to order >= m along D_c' is
   the linear condition a + b Y_i = 0 mod p_i^m, and the pole orders
   (n_+, n_-) at the two places at infinity (deg q even) are fixed in turn,
   so that N(u) = a^2 - q b^2 has degree exactly m deg D_c and is a constant
   times the product of the p_i^m.  Returns {c/m, {a, b}} or None.         *)
RealiseClass[T_, g_, supp_, c_, mmax_: 12, verbose_: False] := Catch[Module[
  {q = T["q"], dq, lc, sqrtlc, unbalancedOK, degSupp, lifts, n, configs, h1, Na, Nb, inf, basis, rows,
   Pc, M, Yk, dM, cols, nplus, nminus, eps, nsmall, w, Sw, sv, row, ns, vec, a, b, NN, Mall, quo, rmd, e, k},
  dq = Exponent[q, g]; lc = Coefficient[q, g, dq];
  sqrtlc = RR[Sqrt[lc]];
  unbalancedOK = EvenQ[dq] && InFieldQ[sqrtlc, Join[{c, lc}, Flatten[CoefficientList[#, g] & /@ Flatten[supp]], CoefficientList[q, g]]];
  degSupp = Total[Exponent[#[[1]], g] & /@ supp];
  lifts = {#[[1]], #[[1]], #[[2]]} & /@ supp;                (* {p_i, p_i^m, Y_i} *)
  Do[
    If[m > 1,
      lifts = Function[{l}, Module[{M2 = Expand[l[[2]] l[[1]]], Yn},
          Yn = Collect[PolynomialRemainder[Expand[(l[[3]] + q InvMod[l[[3]], M2, g])/2], M2, g], g, RR];
          {l[[1]], M2, Yn}]] /@ lifts];
    n = m degSupp;
    configs = If[OddQ[dq], {{Quotient[n, 2], Quotient[n - dq, 2], None}},
      h1 = dq/2;
      Select[Table[With[{Na0 = Max[np, n - np]}, {Na0, Na0 - h1, {np, n - np}}], {np, 0, n}],
             #[[2]] >= 0 && (#[[3, 1]] == #[[3, 2]] || unbalancedOK) &]];
    Do[{Na, Nb, inf} = cfg;
      If[Nb < 0, Continue[]];
      basis = Join[Table[{j, False}, {j, 0, Na}], Table[{j, True}, {j, 0, Nb}]];
      rows = {};
      Do[{Pc, M, Yk} = l; dM = Exponent[M, g];
        cols = Table[PadRight[CoefficientList[PolynomialRemainder[If[bs[[2]], Expand[g^bs[[1]] Yk], g^bs[[1]]], M, g], g], dM], {bs, basis}];
        rows = Join[rows, Transpose[cols]],
        {l, lifts}];
      If[inf =!= None && inf[[1]] != inf[[2]],
        {nplus, nminus} = inf; eps = If[nplus < nminus, 1, -1]; nsmall = Min[nplus, nminus];
        w = Unique["w"];
        Sw = Normal[Series[Sqrt[Expand[(q /. g -> 1/w) w^(2 h1)/lc]], {w, 0, Na - nsmall}]];
        sv = Table[RR[Coefficient[Sw, w, k]], {k, 0, Na - nsmall}];
        Do[row = Table[If[! bs[[2]], If[bs[[1]] == e, 1, 0],
                         With[{kk = bs[[1]] + h1 - e}, If[0 <= kk < Length[sv], eps sqrtlc sv[[kk + 1]], 0]]], {bs, basis}];
          AppendTo[rows, row], {e, Na, nsmall + 1, -1}]];
      If[rows === {}, Continue[]];
      ns = NullSpace[rows, ZeroTest -> (RootReduce[Together[#]] === 0 &)];
      If[ns === {}, Continue[]];
      vec = RRad /@ First[ns];
      a = Collect[Expand[vec[[1 ;; Na + 1]] . g^Range[0, Na]], g, RRad];
      b = Collect[Expand[vec[[Na + 2 ;;]] . g^Range[0, Nb]], g, RRad];
      NN = Expand[a^2 - q b^2];
      Mall = Expand[Times @@ (#[[2]] & /@ lifts)];
      rmd = PolynomialRemainder[NN, Mall, g]; quo = Expand[PolynomialQuotient[NN, Mall, g]];
      If[! PolyZeroQ[rmd, g] || Exponent[quo, g] > 0 || PolyZeroQ[NN, g], Continue[]];
      If[verbose, Print["      residue class ", c, ": div(u) = ", m, " D_c - ", If[inf === None, n, inf], " oo; u = ", a, " + (", b, ")*y, N(u) = ", RR[quo], " * (", Times @@ (#[[1]] & /@ lifts), ")^", m, "; coefficient ", RR[c/m]]];
      Throw[{RRad[c/m], {a, b}}, "rcls"],
      {cfg, configs}],
    {m, 1, mmax}];
  None], "rcls"];


(* RealisePoints[T, p, pts, taus, Y, verbose]
   Realisation of the residues taus at the places pts over a prime p that
   neither the polynomial p nor the y-split realises.  Without a radical the
   places are the principal primes g - rho and the logands are read off.
   Over a curve: try every norm-search
   solution u = a +- b y, accept the first whose vanishing set carries a
   constant ratio tau_P / ord_P(u) (the logand coefficient); the primes it
   leaves are realised jointly by TorsionRealise.  Returns a list of
   {coefficient, pair} or None.                                             *)
RealisePoints[T_, p_, pts_, taus_, Y_, verbose_: False] := Catch[Module[{g, out, ok, uu, hits, gammas, a, b, c, k},
  g = pts[[1, 1]];
  (* no curve: every place over p is the principal prime g - rho *)
  If[T["q"] === None,
    Throw[Cases[Thread[{pts, taus}], {{_, rho_, _}, tau_} /; tau =!= 0 :> {tau, {g - rho, 0}}], "realise"]];
  Do[{a, b, c, k} = sol; out = {}; ok = True;
    Do[uu = {a, sg b};
      hits = Pick[Thread[{pts, taus}], IsZero[ReduceAt[uu[[1]] + uu[[2]] Y, #, Y]] & /@ pts, True];
      gammas = Simplify[#[[2]]/VanishOrder[T, uu, #[[1]], Y]] & /@ hits;
      If[gammas === {} || ! AllTrue[gammas, IsZero[# - gammas[[1]]] &], ok = False; Break[]];
      If[gammas[[1]] =!= 0, AppendTo[out, {gammas[[1]], uu}]],
      {sg, {1, -1}}];
    If[ok,
      If[verbose, Do[Print["      norm factor ", o[[2, 1]], " + (", o[[2, 2]], ")*y (N = ", c, "*(", p, ")^", k, "): coefficient ", o[[1]]], {o, out}]];
      Throw[out, "realise"]],
    {sol, NormSearchAll[T["q"], p, g]}];
  None], "realise"];

(* VinftyResidue[T, f, Y]
   Residue at the place v_oo of a hypertangent top generator t with
   D t = eta (1 + t^2) (Lemma 8.1): delta = 1, uniformiser 1/t, and for
   v_oo(f) = -1 the residue is tau = -(lim f/t)/eta, an element of
   kappa(v_oo) = K_{n-1}(y) returned as a pair.  None if the top generator
   is not hypertangent or v_oo(f) != -1 (deeper poles need the reduction of
   Proposition 7.6 with D acting on the coefficients).                      *)
VinftyResidue[T_, f_, Y_] := Module[{t, et, vinf, lead, vs, ct, tau},
  t = Last[T["gens"]];
  (* Lemma 8.1 needs the hypertangent ABOVE the radical: q free of t, so that
     t = oo is one place with residue field K_{n-1}(y) and y a unit there;
     with q depending on t the curve has its own places over t = oo *)
  If[T["q"] === None || ! FreeQ[T["q"], t], Return[None]];
  et = Can[#/(1 + t^2)] & /@ Last[T["derivs"]];
  If[! FreeQ[et, t] || AllTrue[et, # === 0 &], Return[None]];
  vinf[cc_] := With[{c = Can[cc]}, If[c === 0, None, Exponent[Denominator[c], t] - Exponent[Numerator[c], t]]];
  vs = DeleteCases[vinf /@ f, None];
  If[vs === {} || Min[vs] =!= -1, Return[None]];
  lead[cc_] := With[{c = Can[cc]}, If[c === 0 || vinf[c] =!= -1, 0,
    Can[Coefficient[Numerator[c], t, Exponent[Numerator[c], t]]/Coefficient[Denominator[c], t, Exponent[Denominator[c], t]]]]];
  ct = lead /@ f;
  tau = TDiv[T, ct, et];
  Can /@ (-tau)];



(* ============================================================ stage 3 ===== *)

(* SqrtPolyPart[q, x]
   Polynomial part s of Sqrt[q] in Q((1/x)) for q of even degree with a
   square leading coefficient: the series of Sqrt[q(1/t) t^d] at t = 0 to
   order d/2, read back in x.                                                *)
SqrtPolyPart[q_, x_] := Module[{d = Exponent[q, x], t, qt, ser},
  qt = Expand[(q /. x -> 1/t) t^d];
  ser = Normal[Series[Sqrt[qt], {t, 0, d/2}]];
  Expand[(ser /. t -> 1/x) x^(d/2)]];

(* FundamentalUnit[q, x, maxsteps]
   A unit A + B y of Q[x][y]/(y^2 - q), q of even degree, from the
   continued fraction of Sqrt[q] (Abel-Chebyshev): the first convergent
   h + y k with h^2 - q k^2 a nonzero constant.  Returns {A, B, c} with
   A^2 - q B^2 = c and B monic, or None if no unit appears within maxsteps
   partial quotients or the coefficients grow past the height guard (the
   non-torsion signature: doubly exponential growth).                       *)
FundamentalUnit[q_, x_, maxsteps_: 60] := Catch[Module[{s, P = 0, Q = 1, hm2 = 0, hm1 = 1, km2 = 1, km1 = 0, a, h, k, c, lc, bits},
  If[OddQ[Exponent[q, x]], Throw[None, "unit"]];
  s = SqrtPolyPart[q, x];
  If[! FreeQ[s, Power[_, _Rational]], Throw[None, "unit"]];        (* lc not a square *)
  Do[a = PolynomialQuotient[Expand[P + s], Q, x];
    h = Expand[a hm1 + hm2]; k = Expand[a km1 + km2];
    c = Expand[h^2 - q k^2];
    If[FreeQ[c, x] && c =!= 0 && k =!= 0,
      lc = Coefficient[k, x, Exponent[k, x]];
      Throw[{Expand[h/lc], Expand[k/lc], Expand[c/lc^2]}, "unit"]];
    P = Expand[a Q - P];
    Q = Cancel[(q - P^2)/Q];
    bits = Max[0, Max[BitLength[Abs[Numerator[#]]] + BitLength[Abs[Denominator[#]]] & /@ CoefficientList[Q, x]]];
    If[bits > 8000, Throw[None, "unit"]];
    {hm2, hm1, km2, km1} = {hm1, h, km1, k},
    {maxsteps}];
  None], "unit"];

(* SqrtPolyPartModP[q, x, p]
   Polynomial part of Sqrt[q] in GF(p)[x]((1/x)) for q of even degree
   whose leading coefficient is a square mod p, by matching coefficients;
   None if lc is not a square.                                              *)
SqrtPolyPartModP[q_, x_, p_] := Module[{Q, d2, d, lc, r, coeffs, qc, inv2r, acc},
  Q = PolynomialMod[q, p]; d2 = Exponent[Q, x];
  If[OddQ[d2], Return[None]];
  d = d2/2; lc = Mod[Coefficient[Q, x, d2], p];
  r = SelectFirst[Range[0, p - 1], Mod[#^2 - lc, p] == 0 &, None];
  If[r === None, Return[None]];
  coeffs = ConstantArray[0, d + 1]; coeffs[[1]] = r;              (* coeffs[[k+1]] = coefficient of x^(d-k) *)
  qc = Table[Mod[Coefficient[Q, x, d2 - k], p], {k, 0, d2}];
  inv2r = PowerMod[2 r, -1, p];
  Do[acc = Mod[Sum[coeffs[[i + 1]] coeffs[[k - i + 1]], {i, 1, k - 1}], p];
    coeffs[[k + 1]] = Mod[(qc[[k + 1]] - acc) inv2r, p],
    {k, 1, d}];
  Sum[coeffs[[k + 1]] x^(d - k), {k, 0, d}]];

(* Modular polynomial arithmetic over GF(p). These now delegate directly to the
   Mathilda core, which gained a correct (and FLINT nmod_poly-accelerated)
   `Modulus -> p` for these heads (v0.170): PolynomialMod reduces rational and
   big-integer coefficients (was A6), PolynomialQuotient/Remainder honour
   `Modulus` (was A3), and PolynomialExtendedGCD is correct even when the inputs
   are coprime over Q (was A2). The former hand-rolled substitutes -- a
   CoefficientList/Sum round trip per step and a Euclidean loop over Q -- were
   both the workaround AND the section-C performance bottleneck; the core path is
   10-40x faster on the Jacobian arithmetic. `x` is retained in the signatures so
   call sites are unchanged. *)
ReduceModP[poly_, x_, p_] := PolynomialMod[poly, p];
PolyQuoModP[a_, b_, x_, p_] := PolynomialQuotient[a, b, x, Modulus -> p];
PolyRemModP[a_, b_, x_, p_] := PolynomialRemainder[a, b, x, Modulus -> p];
PolyExtGCDModP[a_, b_, x_, p_] := PolynomialExtendedGCD[a, b, x, Modulus -> p];

(* UnitDegreeModP[q, x, p, maxsteps]
   Order of the class [oo+ - oo-] on the reduction of y^2 = q modulo a
   good prime p: the degree of the fundamental unit of GF(p)[x][y], found
   from the always-periodic continued fraction of Sqrt[q] over GF(p).
   None if oo+- are not rational mod p.                                     *)
UnitDegreeModP[q_, x_, p_, maxsteps_: 10000] := Catch[Module[{s, Qp, P = 0, Q = 1, hm2 = 0, hm1 = 1, km2 = 1, km1 = 0, a, h, k, c},
  s = SqrtPolyPartModP[q, x, p];
  If[s === None, Throw[None, "udeg"]];
  Qp = PolynomialMod[q, p];
  Do[a = PolyQuoModP[P + s, Q, x, p];
    h = PolynomialMod[a hm1 + hm2, p]; k = PolynomialMod[a km1 + km2, p];
    c = PolynomialMod[h^2 - Qp k^2, p];
    If[Exponent[c, x] <= 0 && c =!= 0 && k =!= 0, Throw[Exponent[h, x], "udeg"]];
    P = PolynomialMod[a Q - P, p];
    Q = PolyQuoModP[Qp - P^2, Q, x, p];
    {hm2, hm1, km2, km1} = {hm1, h, km1, k},
    {maxsteps}];
  None], "udeg"];

(* NontorsionCertificate[q, x, nprimes, pmax]
   Proposition 9.4: if [oo+ - oo-] had finite order N, its reduction modulo
   a prime p of good reduction would have order N_p with N = N_p p^a, a >= 0
   (reduction is injective on prime-to-p torsion).  Two good primes whose
   N_p are incompatible with any common N certify non-torsion, hence the
   absence of non-constant units.  Returns {certified, {{p, N_p}, ...}}.  *)
NontorsionCertificate[q_, x_, nprimes_: 4, pmax_: 1000] := Catch[Module[{Q, lc, disc, bad, data = {}, Np, ps, p},
  Q = q; lc = Coefficient[Q, x, Exponent[Q, x]]; disc = Discriminant[Q, x];
  bad = 2 Numerator[lc] Denominator[lc] Numerator[disc] Times @@ (Denominator /@ CoefficientList[Q, x]);
  ps = Select[Range[3, pmax], PrimeQ];
  Do[If[Mod[bad, p] == 0, Continue[]];   (* Do over a packed list now iterates (was A1) *)
    Np = UnitDegreeModP[q, x, p];
    If[Np === None, Continue[]];
    AppendTo[data, {p, Np}];
    If[Length[data] >= nprimes, Break[]],
    {p, ps}];
  {IncompatibleQ[data], data}], "nt"];

(* IncompatibleQ[data]
   No finite N has N = N_p p^a for every {p, N_p} in data: for two primes
   N_{p1}/N_{p2} = p2^b/p1^a, so in lowest terms the numerator must be a
   power of p2 (or 1) and the denominator a power of p1 (or 1).          *)
IncompatibleQ[data_] := Or @@ Flatten[Table[
  With[{p1 = data[[i, 1]], p2 = data[[j, 1]], r = data[[i, 2]]/data[[j, 2]]},
    ! ((Numerator[r] == 1 || FactorInteger[Numerator[r]][[All, 1]] === {p2}) &&
       (Denominator[r] == 1 || FactorInteger[Denominator[r]][[All, 1]] === {p1}))],
  {i, Length[data]}, {j, i + 1, Length[data]}]];

(* ------------------------------------------- Jacobian arithmetic mod p
   Cantor's algorithm on w^2 = f(jv), deg f = 2g + 1 (one place at
   infinity), over GF(p): divisor classes as reduced Mumford pairs {u, v},
   u monic of degree <= g, v^2 = f mod u, polynomials in the symbol jv with
   integer coefficients in 0..p-1 (PolynomialMod); the identity is {1, 0}.
   Used by NontorsionDivisorCertificate (Proposition 9.4 for a residue
   divisor supported at finite places).
   The divisions are PolyQuoModP/PolyRemModP (division over Q, then
   reduction: Mathilda's PolynomialQuotient/PolynomialRemainder take no
   Modulus option) and the extended gcd is PolyExtGCDModP (Mathilda's
   PolynomialExtendedGCD ignores the Modulus for inputs coprime over Q). *)

CantorAdd[f_, D1_, D2_, p_] := Module[{u1, v1, u2, v2, d1, e1, e2, d, c1, c2, s1, s2, u, v, t, g, monic},
  monic[w_] := PolynomialMod[PowerMod[Coefficient[w, jv, Exponent[w, jv]], -1, p] w, p];
  {u1, v1} = D1; {u2, v2} = D2;
  {d1, {e1, e2}} = PolyExtGCDModP[u1, u2, jv, p];                                (* d1 = e1 u1 + e2 u2 *)
  {d, {c1, c2}} = PolyExtGCDModP[d1, v1 + v2, jv, p];                             (* d = c1 d1 + c2 (v1 + v2) *)
  s1 = PolynomialMod[c1 e1, p]; s2 = PolynomialMod[c1 e2, p];
  u = PolyQuoModP[u1 u2, d^2, jv, p];
  t = PolynomialMod[s1 u1 v2 + s2 u2 v1 + c2 (v1 v2 + f), p];
  v = PolyRemModP[PolyQuoModP[t, d, jv, p], u, jv, p];
  g = (Exponent[f, jv] - 1)/2;
  While[Exponent[u, jv] > g,
    u = monic[PolyQuoModP[f - v^2, u, jv, p]];
    v = PolyRemModP[-v, u, jv, p]];
  {monic[u], PolynomialMod[v, p]}];

CantorMul[f_, n_, D_, p_] := Module[{R = {1, 0}, Q = D, k = n},
  While[k > 0,
    If[OddQ[k], R = CantorAdd[f, R, Q, p]];
    Q = CantorAdd[f, Q, Q, p]; k = Quotient[k, 2]];
  R];

(* DivisorOrderModP[f, pts, p]
   The order of the class of sum_P n_P (P - oo) in the Jacobian of
   w^2 = f(jv) over GF(p), pts = {{{x0, y0}, n_P}, ...}: the Mumford
   representative of the divisor by Cantor's algorithm, a multiple M of the
   order by baby-step giant-step in the Weil interval [(Sqrt[p] - 1)^2g,
   (Sqrt[p] + 1)^2g] (which contains #J(GF(p))), and the least divisor of
   M annihilating the class; None if no multiple is found.                *)
DivisorOrderModP[f_, pts_, p_] := Module[{D = {1, 0}, P, g, lo, hi, m, baby, B, G, step, M = None, j},
  Do[P = {PolynomialMod[jv - pt[[1, 1]], p], Mod[Sign[pt[[2]]] pt[[1, 2]], p]};
    Do[D = CantorAdd[f, D, P, p], {Abs[pt[[2]]]}],
    {pt, pts}];
  If[D === {1, 0}, Return[1]];
  g = (Exponent[f, jv] - 1)/2;
  lo = Max[1, Floor[N[(Sqrt[p] - 1)^(2 g)]] - 1]; hi = Ceiling[N[(Sqrt[p] + 1)^(2 g)]] + 1;
  m = Floor[Sqrt[hi - lo]] + 1;
  baby = <||>; B = {1, 0};
  Do[If[! KeyExistsQ[baby, B], baby[B] = j]; B = CantorAdd[f, B, D, p], {j, 0, m - 1}];
  G = CantorMul[f, lo, D, p]; step = CantorMul[f, m, D, p];
  Do[j = Lookup[baby, Key[G], None];   (* Lookup with a list-valued key now works (was A5) *)
    If[j =!= None && lo + i m - j > 0, M = lo + i m - j; Break[]];        (* [lo + i m] D = [j] D *)
    G = CantorAdd[f, G, step, p],
    {i, 0, m}];
  If[M === None, Return[None]];
  SelectFirst[Divisors[M], CantorMul[f, #, D, p] === {1, 0} &, None]];

(* OddModel[f, r, p]
   For f in GF(p)[jv] of even degree d with a root r: the odd-degree model
   w^2 = f'(jv) of w^2 = f, f'(v) = v^d f(r + 1/v) of degree d - 1 (the
   coefficient of v^d is f(r) = 0), and the point map
   {x0, y0} -> {1/(x0 - r), y0/(x0 - r)^(d/2)}.                            *)
OddModel[f_, r_, p_] := With[{d = Exponent[f, jv]},
  {PolynomialMod[Expand[(f /. jv -> r + 1/jv) jv^d], p],
   Function[{x0, y0}, With[{s = PowerMod[Mod[x0 - r, p], -1, p]}, {s, Mod[y0 PowerMod[s, d/2, p], p]}]]}];

HornerModP[cs_, t_, p_] := Fold[Mod[#1 t + #2, p] &, 0, Reverse[cs]];      (* cs lowest first *)

(* NontorsionDivisorCertificate[q, x, places, mu, nprimes, pmax, budget]
   Certify that the class of a degree-0 divisor sum n_P P on y^2 = q with
   places P = (x0, y0) of constant coordinates is NOT torsion, by reduction
   modulo primes of good reduction (Proposition 9.4 for the residue divisor
   that Algorithm 3(d) leaves unrealised): if the class had finite order N,
   its reduction at a degree-one prime above p of the number field
   K = Q(theta) of the coordinates would have order N_p with N = N_p p^a,
   a >= 0 (reduction is injective on the prime-to-p torsion), so two primes
   with incompatible N_p exclude every finite order.  N_p is the order of
   the reduced class in the Jacobian over GF(p) (DivisorOrderModP), on an
   odd-degree model when deg q is even (a root of q mod p moved to
   infinity; primes without one are skipped); primes whose Weil bound
   exceeds budget are skipped.  places = {{{a, b}, n}, ...} with a, b the
   coordinates of x0, y0 over Q on the power basis of theta (lowest first),
   mu the minimal polynomial of theta as its CoefficientList (None for
   K = Q).  Returns {certified, {{p, N_p}, ...}}.                          *)
NontorsionDivisorCertificate[q_, x_, places_, mu_, nprimes_: 2, pmax_: 1000, budget_: 2 10^7] := Catch[Module[
  {d = Exponent[q, x], twoG, qc, muc, bad, data = {}, red, thetas, th, ev, pts, f, roots, mp, Np, z, ps, p},
  twoG = If[OddQ[d], d - 1, d - 2];
  qc = CoefficientList[q, x]; muc = If[mu === None, {0, 1}, mu];
  bad = 2 Times @@ (Numerator[#] Denominator[#] & /@ DeleteCases[Join[qc, muc,
      {Discriminant[q, x], Discriminant[Sum[muc[[i]] z^(i - 1), {i, Length[muc]}], z]}, Flatten[places[[All, 1]]]], 0]);
  ps = Select[Range[3, pmax], PrimeQ];
  Do[If[Mod[bad, p] == 0 || (Sqrt[p] + 1)^twoG > budget, Continue[]];   (* Do over a packed list now iterates (was A1) *)
    red = Function[c, Mod[Numerator[c] PowerMod[Denominator[c], -1, p], p]];
    thetas = Select[Range[0, p - 1], HornerModP[red /@ muc, #, p] == 0 &];
    If[Length[thetas] == 0, Continue[]];               (* no degree-one prime of K above p *)
    th = thetas[[1]];
    ev = Function[a, HornerModP[red /@ a, th, p]];
    pts = {{ev[#[[1, 1]]], ev[#[[1, 2]]]}, #[[2]]} & /@ places;
    f = Sum[red[qc[[i]]] jv^(i - 1), {i, Length[qc]}];
    If[EvenQ[d],
      roots = Select[Range[0, p - 1], Mod[f /. jv -> #, p] == 0 && ! MemberQ[pts[[All, 1, 1]], #] &];
      If[Length[roots] == 0, Continue[]];
      {f, mp} = OddModel[f, roots[[1]], p];
      pts = {mp[#[[1, 1]], #[[1, 2]]], #[[2]]} & /@ pts];
    If[AnyTrue[pts, Mod[#[[1, 2]]^2 - (f /. jv -> #[[1, 1]]), p] != 0 &], Continue[]];
    Np = DivisorOrderModP[f, pts, p];
    If[Np === None, Continue[]];
    AppendTo[data, {p, Np}];
    If[IncompatibleQ[data], Throw[{True, data}, "ntd"]];
    If[Length[data] >= nprimes, Break[]],
    {p, ps}];
  {False, data}], "ntd"];

(* ResidueFree[T, g, Y]
   For n = 1, D = d/dx: certify that the differential g dx (g a pair) has
   zero residue at every affine place.  Sub-critical poles carry no residue
   (Theorem 7.4(ii)); at the others the classical residues are computed by
   series, on both sheets at unramified places and with the local parameter
   x = rho + e^2 at branch places.  A False here means the realisation is
   not trusted and no certificate is issued.                                *)
ResidueFree[T_, g_, Y_] := Catch[Module[{x = T["gens"][[1]], q = T["q"], den, e, p, branch, rts, yser, loc, res},
  den = PolynomialLCM[Denominator[Can[g[[1]]]], Denominator[Can[g[[2]]]]];
  Do[p = fac[[1]];
    If[FreeQ[p, x], Continue[]];
    branch = Exponent[PolynomialGCD[p, q], x] > 0;
    If[VP[g, p, T, branch] > -(If[branch, 2, 1]), Continue[]];
    rts = x /. Solve[p == 0, x];
    Do[If[! IsZero[q /. x -> rho],
        yser = Sqrt[q /. x -> rho + e];
        Do[loc = (g[[1]] /. x -> rho + e) + sg (g[[2]] /. x -> rho + e) yser;
          res = SeriesCoefficient[loc, {e, 0, -1}];
          If[! IsZero[res], Throw[False, "rf"]],
          {sg, {1, -1}}],
        loc = 2 e ((g[[1]] /. x -> rho + e^2) + (g[[2]] /. x -> rho + e^2) e Sqrt[Cancel[(q /. x -> rho + e^2)/e^2]]);
        res = SeriesCoefficient[loc, {e, 0, -1}];
        If[! IsZero[res], Throw[False, "rf"]]],
      {rho, rts}],
    {fac, Rest[FactorList[den]]}];
  True], "rf"];


(* ------------------------------------------- bounds at the places (S8.2) *)

(* OrdInf[e, g]: v_oo(e) = deg_g(den) - deg_g(num) of a rational function;
   Infinity for 0.  LcInf[e, g]: the leading coefficient at g = oo, a
   function of the other generators.                                        *)
OrdInf[e_, g_] := With[{c = Can[e]}, If[c === 0, Infinity, Exponent[Denominator[c], g] - Exponent[Numerator[c], g]]];
LcInf[e_, g_] := With[{c = Can[e]}, Can[Coefficient[Numerator[c], g, Exponent[Numerator[c], g]]/
  Coefficient[Denominator[c], g, Exponent[Denominator[c], g]]]];

(* VInf[T, u, g]
   The valuations of the tuple u at the places over g = oo, one entry per
   place in the units of its uniformiser, and whether they are exact: for
   m <= 2 exact (a series at the two unramified places when the two
   coordinates could cancel), for m >= 3 with g in q the Gauss lower bound
   of Part I, S6.1.                                                          *)
VInf[T_, u0_, g_] := Module[{q = T["q"], m = T["m"], u = TPad[T, u0], d, a, b, va, vb, vy, uu, qs, lc, out, s, e, dp, vals, expr, ser, ord},
  If[q === None || FreeQ[q, g], Return[{{Min[OrdInf[#, g] & /@ u]}, True}]];
  d = Exponent[q, g];
  If[m >= 3,
    s = GCD[m, d]; e = m/s; dp = d/s;
    vals = Table[If[u[[i + 1]] === 0, Nothing,
      e OrdInf[u[[i + 1]], g] - (i dp - e If[T["E"][[i + 1]] === 1, 0, Exponent[T["E"][[i + 1]], g]])], {i, 0, m - 1}];
    Return[{ConstantArray[If[vals === {}, Infinity, Min[vals]], s], False}]];
  {a, b} = u[[1 ;; 2]];
  va = OrdInf[a, g]; vb = OrdInf[b, g];
  If[OddQ[d], Return[{{Min[2 va, 2 vb - d]}, True}]];           (* ramified: parities differ *)
  vy = -d/2;
  If[a === 0 || b === 0 || va =!= vb + vy, Return[{{Min[va, vb + vy], Min[va, vb + vy]}, True}]];
  uu = Unique["u"];
  qs = Expand[(q /. g -> 1/uu) uu^d]; lc = qs /. uu -> 0;          (* q(1/u) u^d, a unit at u = 0; Sqrt[lc] split off, or Series returns a Piecewise *)
  out = Table[
    expr = Together[(a /. g -> 1/uu) + sg Sqrt[lc] (b /. g -> 1/uu) Sqrt[Cancel[qs/lc]]/uu^(d/2)];
    ord = Infinity;
    Do[ser = Normal[Series[expr, {uu, 0, order}]];
      If[ser =!= 0, ord = Exponent[ser, uu, Min]; Break[]], {order, {va + 4, va + 12, va + 30}}];
    ord, {sg, {1, -1}}];
  {out, True}];

(* Kind[T, k]: "prim", "exp" or "tan" if t_k is a primitive, hyperexponential
   or hypertangent monomial over the field below it, else "other";
   IndependentQ[T, i]: no other generator's derivative, nor q, involves t_i. *)
Kind[T_, k_] := Module[{g = T["gens"][[k]], above = T["gens"][[k ;;]], d = T["derivs"][[k]]},
  Which[AllTrue[d, FreeQ[Can[#], Alternatives @@ above] &], "prim",
        AllTrue[d, FreeQ[Can[#/g], Alternatives @@ above] &], "exp",
        AllTrue[d, FreeQ[Can[#/(1 + g^2)], Alternatives @@ above] &], "tan",
        True, "other"]];
IndependentQ[T_, i_] := (T["q"] === None || FreeQ[T["q"], T["gens"][[i]]]) &&
  AllTrue[Delete[T["derivs"], i], FreeQ[#, T["gens"][[i]]] &];

(* LamCoeffs[T, k, g]: the leading coefficients at g = oo (g === None: the
   whole coordinates) of D t_k, of D t_k / t_k for a hyperexponential and of
   D t_k / (1 + t_k^2) for a hypertangent t_k.                             *)
LamCoeffs[T_, k_, g_] := Module[{gk = T["gens"][[k]], div, comps},
  div = Switch[Kind[T, k], "exp", gk, "tan", 1 + gk^2, _, 1];
  comps = Can[#/div] & /@ DeleteCases[T["derivs"][[k]], 0];
  If[g === None, comps, LcInf[#, g] & /@ comps]];

(* LcPlace[T, u, g]: the leading coefficients of the tuple u at the places
   over g = oo, one per place as in VInf, or None where the implementation
   does not read them (m >= 3, a ramified place, a coordinate carrying y off
   the curve variable).                                                    *)
LcPlace[T_, u0_, g_] := Module[{q = T["q"], m = T["m"], u = TPad[T, u0], d, a, b, va, uu, qs, lc, out, expr, ser, k},
  If[q === None || FreeQ[q, g],
    If[u[[1]] === 0 || AnyTrue[Rest[u], # =!= 0 &], Return[None]];
    Return[{LcInf[u[[1]], g]}]];
  d = Exponent[q, g];
  {a, b} = u[[1 ;; 2]];
  If[m >= 3 || OddQ[d] || (a === 0 && b === 0), Return[None]];
  va = Min[OrdInf[a, g], OrdInf[b, g] - d/2];
  uu = Unique["u"];
  qs = Expand[(q /. g -> 1/uu) uu^d]; lc = qs /. uu -> 0;
  out = Table[
    expr = Together[(a /. g -> 1/uu) + sg Sqrt[lc] (b /. g -> 1/uu) Sqrt[Cancel[qs/lc]]/uu^(d/2)];
    k = None;
    Do[ser = Normal[Series[expr, {uu, 0, order}]];
      If[ser =!= 0, k = Can[Coefficient[ser, uu, Exponent[ser, uu, Min]]]; Break[]], {order, {va + 4, va + 12, va + 30}}];
    k, {sg, {1, -1}}];
  If[MemberQ[out, None], None, out]];

(* ClassMod[e, pp, g, n, gens]: the class of e / pp^n modulo the prime pp of
   F[g] when it is a constant of Fbar (a number), else None.              *)
ClassMod[e_, pp_, g_, n_, gens_] := Module[{c = Can[e/pp^n], egcd, cls},
  egcd = Quiet[Check[PolynomialExtendedGCD[Denominator[c], pp, g], $Failed]];
  If[egcd === $Failed || ! FreeQ[egcd[[1]], g], Return[None]];
  cls = Quiet[Check[Can[PolynomialRemainder[Expand[Numerator[c] egcd[[2, 1]]/egcd[[1]]], pp, g]], $Failed]];
  If[cls === $Failed || ! FreeQ[cls, Alternatives @@ gens], None, cls]];

(* InQSpanQ[rho, mus]: rho lies in the Q-span of the algebraic numbers mus,
   the kernel-degree test of the extended (K1): N_v = {a : -a rho in
   sum_k Z mu_k}.  True on any failure, which leaves the place to the guess. *)
InQSpanQ[rho_, mus_] := Module[{nf, vecs, n},
  If[mus === {}, Return[False]];
  nf = Quiet[Check[ToNumberField[Prepend[mus, rho]], $Failed]];
  If[nf === $Failed || ! ListQ[nf] ||
     ! AllTrue[nf, Head[#] === AlgebraicNumber || Element[#, Rationals] === True &], Return[True]];
  vecs = If[Head[#] === AlgebraicNumber, Last[#], {#}] & /@ nf;
  n = Max[Length /@ vecs];
  vecs = PadRight[#, n] & /@ vecs;
  MatrixRank[vecs] == MatrixRank[Rest[vecs]]];

(* InfConsts[T, g, pl, names]: rho ("pi", the coefficient of pi^(s+1) in
   D pi = -D g pi^2) and the mu_k of the named hyperexponential generators
   at the pl-th place over g = oo, as numbers where LcPlace reads them.   *)
InfConsts[T_, g_, pl_, names_] := Module[{gens = T["gens"], i, lc},
  i = Position[gens, g][[1, 1]];
  Association @@ Table[
    lc = LcPlace[T, If[nm === "pi", -T["derivs"][[i]], Can[#/nm] & /@ T["derivs"][[Position[gens, nm][[1, 1]]]]], g];
    nm -> If[lc === None || Length[lc] < pl, None, lc[[pl]]], {nm, names}]];

(* SpecialConsts[T, pp, Dpi, s, g, names]: the same at the places over the
   special pp = pi of F[g]: the classes of D pi / pi^(s+1) and of
   (D t_k / t_k) / pi^s modulo pp, when they are constants of Fbar, the same
   at every place over pp.                                                 *)
SpecialConsts[T_, pp_, Dpi_, s_, g_, names_] := Module[{gens = T["gens"], pair, cls},
  Association @@ Table[
    pair = If[nm === "pi", Dpi, Can[#/nm] & /@ T["derivs"][[Position[gens, nm][[1, 1]]]]];
    cls = If[# === 0, 0, ClassMod[#, pp, g, If[nm === "pi", s + 1, s], gens]] & /@ pair;
    nm -> If[cls[[1]] === None || AnyTrue[Rest[cls], # =!= 0 &], None, cls[[1]]], {nm, names}]];

(* DecideBound: Algorithm 6, Step 3 -- the criterion (K1)-(K4) of Proposition
   8.11 that proves the bound at the place, or None.  own: the generators of
   the place; mono: the index of the monomial the place belongs to or None;
   A: the attaining generators; kinds: their kinds; rhoFree: the leading
   coefficient of D pi is free of A; lamK1/lamK4: the leading coefficients of
   the attaining generators are free of A / are constants; rhoR: the
   coefficient for (K3) when the residue field is F(x), else None; consts:
   names -> the values at the place of rho ("pi") and of the mu_k of the
   named hyperexponential generators, for (K1) with attaining
   hyperexponentials (Remark 8.14).                                       *)
DecideBound[T_, own_, mono_, r_, A_, kinds_, rhoFree_, lamK1_, lamK4_, rhoR_, upperAllAttain_, consts_] := Module[{x = T["gens"][[1]], R, exps, vals},
  Which[
    mono =!= None && Kind[T, mono] =!= "other" && IndependentQ[T, mono], "K2",
    r == 0 && rhoFree && lamK1 && AllTrue[kinds, MemberQ[{"prim", "exp"}, #] &] &&
      (exps = Pick[A, kinds, "exp"]; exps === {} ||
        (vals = consts[Prepend[exps, "pi"]]; vals =!= None && AllTrue[Values[vals], # =!= None && NumericQ[#] &] &&
          ! InQSpanQ[vals["pi"], vals[#] & /@ exps])), "K1",
    r >= 1 && rhoR =!= None && (R = Integrate[Can[rhoR], x]; ! FreeQ[R, Log | ArcTan | ArcTanh | RootSum]), "K3",
    r >= 1 && rhoFree && own =!= {} && T["q"] =!= None && upperAllAttain && lamK4 &&
      AllTrue[kinds, MemberQ[{"exp", "tan"}, #] &] && Count[kinds, "exp"] <= 1, "K4",
    True, None]];

(* InfData[key, T, rem]: Algorithm 6 at the places over t_i = oo for every
   generator, once per analysis: per generator {criterion or None, e_P, the
   list of {s, r, v_P(rem)} over the places, chain-possible}.              *)
InfData[key_, T_, rem_] := Module[{gens = T["gens"], q = T["q"], m = T["m"], curve, out, g, own, items, exact, vrem, d, eP, lam, rhoSyms, rhoR, crit, places, chain, shifts, s, r, A, c, i},
  If[KeyExistsQ[$analyses[key], "inf"], Return[$analyses[key]["inf"]]];
  curve = Select[gens, ! FreeQ[q, #] &];
  out = Table[
    g = gens[[i]];
    own = If[MemberQ[curve, g], Union[{g}, curve], {g}];
    items = Table[{gens[[k]], VInf[T, T["derivs"][[k]], g], VInf[T, TScalar[T, gens[[k]]], g][[1]]}, {k, Length[gens]}];
    If[q =!= None, AppendTo[items, {"y", VInf[T, T["Dw"][[2]], g], VInf[T, TUnit[T, 1], g][[1]]}]];
    exact = AllTrue[items, #[[2, 2]] &];
    vrem = VInf[T, rem, g][[1]];
    d = If[MemberQ[curve, g], Exponent[q, g], 0];
    eP = Which[! MemberQ[curve, g], 1, m == 2, If[OddQ[d], 2, 1], True, m/GCD[m, d]];
    lam = Table[LamCoeffs[T, k, g], {k, Length[gens]}];
    rhoSyms = LcInf[#, g] & /@ DeleteCases[T["derivs"][[i]], 0];
    rhoR = If[q === None && Length[gens] == 2 && i == 2 && T["derivs"][[1, 1]] === 1, -LcInf[T["derivs"][[2, 1]], g], None];
    crit = None; places = {}; chain = False;
    Do[shifts = Association @@ Table[If[it[[2, 1, pl]] === Infinity, Nothing, it[[1]] -> it[[2, 1, pl]] - it[[3, pl]]], {it, items}];
      s = Min[Values[shifts]];
      r = shifts[g] - s;
      A = Select[gens, Lookup[shifts, #, None] === s && ! MemberQ[own, #] &];
      c = If[exact,
        DecideBound[T, own, i, r, A, Kind[T, #] & /@ (Position[gens, #][[1, 1]] & /@ A),
          FreeQ[rhoSyms, Alternatives @@ A],
          AllTrue[A, FreeQ[lam[[Position[gens, #][[1, 1]]]], Alternatives @@ A] &],
          AllTrue[A, FreeQ[lam[[Position[gens, #][[1, 1]]]], Alternatives @@ Complement[gens, own]] &],
          rhoR, AllTrue[Complement[gens, own], MemberQ[A, #] &],
          Function[names, InfConsts[T, g, pl, names]]],
        None];
      If[c === None, crit = None;
        chain = r >= 1 && AnyTrue[A, Kind[T, Position[gens, #][[1, 1]]] === "prim" &]; Break[]];
      If[crit === None, crit = c];
      AppendTo[places, {s, r, vrem[[pl]]}],
      {pl, Length[vrem]}];
    {crit, eP, places, chain},
    {i, Length[gens]}];
  $analyses[key, "inf"] = out;
  out];

(* SpecialData[key, T, pp]: Algorithm 6 at the places over the special pp:
   {criterion or None, s, r, e_P, branch}; the derivatives being defined over
   F, a factor split off a constant-coefficient special shares the data of
   its parent.                                                              *)
SpecialData[key_, T_, pp_] := Module[{gens = T["gens"], q = T["q"], m = T["m"], parent, branch, own, vp, shifts, exact, vD, dk, pi, Dpi, sigPi, s, A, lam, rhoFree, mono, i, kd, crit, r},
  If[KeyExistsQ[$analyses[key]["spec"], pp], Return[$analyses[key]["spec"][pp]]];
  parent = Lookup[$analyses[key]["parent"], pp, None];
  If[parent =!= None, $analyses[key, "spec", pp] = SpecialData[key, T, parent]; Return[$analyses[key]["spec"][pp]]];
  branch = q =!= None && Vp[q, pp, gens] > 0;
  own = Select[gens, ! FreeQ[pp, #] &];
  vp[u_] := VP[u, pp, T, branch];
  shifts = <||>; exact = ! branch;
  Do[vD = vp[T["derivs"][[k]]];
    If[vD =!= Infinity, shifts[gens[[k]]] = vD - vp[TScalar[T, gens[[k]]]]];
    dk = T["derivs"][[k]];
    If[q =!= None && dk[[1]] =!= 0 && dk[[2]] =!= 0 && Vp[dk[[1]], pp, gens] == Vp[dk[[2]], pp, gens], exact = False],
    {k, Length[gens]}];
  If[q =!= None, vD = vp[T["Dw"][[2]]]; If[vD =!= Infinity, shifts["y"] = vD - vp[TUnit[T, 1]]]];
  pi = If[branch, TUnit[T, 1], TScalar[T, pp]];
  Dpi = TowerD[T, pi];
  sigPi = vp[Dpi] - vp[pi];
  s = Min[Append[Values[shifts], sigPi]];
  r = sigPi - s;
  A = Select[gens, Lookup[shifts, #, None] === s && ! MemberQ[own, #] &];
  lam = Table[LamCoeffs[T, k, None], {k, Length[gens]}];
  rhoFree = ! branch && FreeQ[DeleteCases[Can /@ Dpi, 0], Alternatives @@ A];
  mono = None;
  If[Length[own] == 1,
    i = Position[gens, own[[1]]][[1, 1]]; kd = Kind[T, i];
    If[(kd === "exp" && Expand[pp - gens[[i]]] === 0) || (kd === "tan" && Expand[(pp - gens[[i]])^2 + 1] === 0), mono = i]];
  crit = If[exact,
    DecideBound[T, own, mono, r, A, Kind[T, Position[gens, #][[1, 1]]] & /@ A, rhoFree,
      AllTrue[A, FreeQ[lam[[Position[gens, #][[1, 1]]]], Alternatives @@ A] &],
      AllTrue[A, FreeQ[lam[[Position[gens, #][[1, 1]]]], Alternatives @@ Complement[gens, own]] &],
      None, AllTrue[Complement[gens, own], MemberQ[A, #] &],
      If[! branch && Length[own] == 1, Function[names, SpecialConsts[T, pp, Dpi, s, own[[1]], names]], Function[names, None]]],
    None];
  $analyses[key, "spec", pp] = {crit, s, r, If[branch, m, 1], branch};
  $analyses[key]["spec"][pp]];

(* PlaceBounds[key, T, rem, denv, unkLogs, retry, verbose]
   Step 15 of Algorithm 4 by Algorithm 6: the exponent of every special
   prime and the degree bound in every generator, from the shift of the
   derivation at the places over sigma and over t_i = oo (Theorem 8.7,
   Corollary 8.9); the classical guess, raised by the retry count, at the
   places where Proposition 8.11 is silent.  Returns {bounds, exps (aligned
   with unkLogs), proved}.                                                  *)
PlaceBounds[key_, T_, rem_, denv0_, unkLogs_, retry_, verbose_] := Module[
  {gens = T["gens"], q = T["q"], nc, crits = <||>, exps, denv = denv0, nb, db, bounds, inf, crit, eP, places, chain, b, g, s, r, parent, vg, i, branch},
  nc = If[q === None, 1, T["n"]];
  exps = Table[{crit, s, r, eP, branch} = SpecialData[key, T, pl[[1]]];
    parent = Lookup[$analyses[key]["parent"], pl[[1]], None];
    vg = VP[rem, If[parent === None, pl[[1]], parent], T, branch];         (* a lower bound at the factor *)
    crits[{"s", pl[[1]]}] = crit;
    If[crit =!= None, Ceiling[Max[0, s - vg + r]/eP], Max[0, -vg] + retry],
    {pl, unkLogs}];
  Do[denv *= unkLogs[[i, 1]]^exps[[i]], {i, Length[unkLogs]}];
  nb = Numerator /@ rem[[1 ;; nc]]; db = (Denominator[#] denv) & /@ rem[[1 ;; nc]];
  inf = InfData[key, T, rem];
  bounds = ConstantArray[None, Length[gens]];
  Do[g = gens[[i]];
    {crit, eP, places, chain} = inf[[i]];
    b = If[crit =!= None,
      Exponent[denv, g] + Max[Ceiling[(Max[0, #[[1]] - #[[3]]] + #[[2]])/eP] & /@ places],
      Max[Exponent[#, g] & /@ Join[nb, db]] + 2 + If[chain, Max[Prepend[bounds[[i + 1 ;;]], 0]], 0]];
    bounds[[i]] = b; crits[{"g", g}] = crit,
    {i, Length[gens], 1, -1}];
  If[verbose, Print["  bounds: ", StringRiffle[Table[ToString[gens[[i]]] <> ": " <> ToString[bounds[[i]]] <> " [" <> ToString[crits[{"g", gens[[i]]}] /. None -> "guess"] <> "]", {i, Length[gens]}], ", "],
    If[unkLogs === {}, "", ";  specials: " <> StringRiffle[Table["(" <> ToString[unkLogs[[i, 1]]] <> ")^" <> ToString[exps[[i]]] <> " [" <> ToString[crits[{"s", unkLogs[[i, 1]]}] /. None -> "guess"] <> "]", {i, Length[unkLogs]}], ", "]],
    If[AllTrue[Values[crits], # =!= None &], "  (all proved)", "  (a guess is in play)"]]];
  {bounds, exps, AllTrue[Values[crits], # =!= None &]}];

(* VerifiedResidueFree[rem, T]: re-run the residue computation of Steps
   7--13 on the residual (through iPIM with a trivial ansatz): it is verified
   residue-free when the analysis finds nothing left to realise.           *)
VerifiedResidueFree[rem_, T_] := Module[{key = Hash[{rem, T}], gens = T["gens"], q = T["q"], comps, den, ok, p, branch, eta, delta, special},
  Catch[iPIM[rem, T, "Bounds" -> ConstantArray[0, Length[gens]], "Verbose" -> False], "PIM"];
  ok = KeyExistsQ[$analyses, key] && $analyses[key]["detLogs"] === {};
  (* the analysis evaluates a deeper-than-critical pole at a normal prime
     over the curve only for n = 1, D = d/dx (DeepResidues), and in a
     transcendental tower only along a residue direction; a pole it passes
     over in silence leaves the residual unverified *)
  If[ok,
    comps = DeleteCases[Can /@ TPad[T, rem], 0];
    den = If[comps === {}, 1, PolynomialLCM @@ (Denominator /@ comps)];
    Do[p = fac[[1]];
      If[FreeQ[p, Alternatives @@ gens], Continue[]];
      {branch, eta, delta, special} = ClassifyPrime[T, p];
      If[special || VP[rem, p, T, branch] >= -delta, Continue[]];
      If[q === None,
        If[ResidueDir[T, p] === None, ok = False; Break[]],
        If[! (Length[gens] == 1 && T["derivs"][[1]] === TUnit[T, 0] && delta == 1 && ! branch && PointsOver[T, p] =!= None), ok = False; Break[]]],
      {fac, Rest[FactorList[den]]}]];
  ok];


(* ------------------------------------------ m >= 3: places at infinity, realisation, units *)

(* InfinityData[T, g, lmax]
   The places over g = oo of y^m = q(g), deg q = d, lc(q) = c, for the curve
   variable g: s = gcd(m, d) places oo_j (j = 0..s-1), each of ramification
   index e = m/s, with local parameter t, g = t^-e and
       y = zeta^j c^(1/m) t^-d' S(t^e),  d' = d/s,  S(u) = (q(1/u) u^d / c)^(1/m),
   so that g^k w_i = zeta^(ij) c^(i/m) t^-v(k,i) F_i(t^e), F_i = S^i / Et_i,
   Et_i(u) = u^deg(E_i) E_i(1/u), v(k,i) = e k + i d' - e deg(E_i) the pole
   order.  InfCoef[inf, j, k, i, r] is the coefficient of t^r in the
   expansion of g^k w_i at oo_j (series to order u^lmax).                   *)
InfinityData[T_, g_, lmax_] := Module[{q = T["q"], m = T["m"], d, lc, s, e, dp, u, qt, S, F = {}, degE = {}, Spow = 1, Ei, dE, Et, Fi},
  d = Exponent[q, g]; lc = Coefficient[q, g, d];
  s = GCD[m, d]; e = m/s; dp = d/s;
  u = Unique["u"];
  qt = Expand[Expand[(q /. g -> 1/u) u^d]/lc];
  S = Normal[Series[qt^(1/m), {u, 0, lmax}]];
  Do[Ei = T["E"][[i + 1]]; dE = If[Ei === 1, 0, Exponent[Ei, g]];
    Et = Expand[u^dE (Ei /. g -> 1/u)];
    Fi = Normal[Series[Spow/Et, {u, 0, lmax}]];
    AppendTo[F, Table[RRad[Coefficient[Fi, u, l]], {l, 0, lmax}]];
    AppendTo[degE, dE];
    Spow = Normal[Series[Expand[Spow S], {u, 0, lmax}]],
    {i, 0, m - 1}];
  <|"s" -> s, "e" -> e, "dp" -> dp, "degE" -> degE, "lmax" -> lmax, "F" -> F, "zk" -> RootOfUnity[m], "ck" -> MRoot[lc, m]|>];
InfV[inf_, k_, i_] := inf["e"] k + i inf["dp"] - inf["e"] inf["degE"][[i + 1]];
InfCoef[inf_, j_, k_, i_, r_] := Module[{l, rr},
  {l, rr} = QuotientRemainder[r + InfV[inf, k, i], inf["e"]];
  If[rr != 0 || l < 0 || l > inf["lmax"], 0, RRad[inf["zk"]^(i j) inf["ck"]^i inf["F"][[i + 1, l + 1]]]]];

(* FindElement[T, g, affine, zorders, inf]
   An element u = Sum_i a_i(g) w_i of the integral closure O with u = 0
   modulo (M, sheet) for every affine condition {M, {B_0, ..., B_{m-1}}} of
   affine (B_i = Y^i/E_i mod M for the Hensel-lifted sheet Y: the linear
   conditions Sum_i a_i B_i = 0 mod M) and v_{oo_j}(u) >= zorders[[j]] at the
   s places at infinity (a negative entry allows a pole of that order, a
   positive one demands a zero), by linear algebra: the basis g^k w_i is
   bounded by the largest pole allowed (the leading terms cannot cancel at
   every place at once), the conditions at infinity are the vanishing of
   Laurent coefficients.  For a divisor of degree 0 the solution space has
   dimension <= 1.  Returns the tuple or None.                              *)
FindElement[T_, g_, affine_, zorders_, inf_] := Module[{m = T["m"], s = inf["s"], e = inf["e"], Dmax, kmax, basis, rows = {}, M, Bs, dM, cols, r, ns, vec, a},
  Dmax = Max[Join[-zorders, {0}]];
  kmax = Floor[Dmax/e] + Max[inf["degE"]] + 1;
  basis = Select[Flatten[Table[{k, i}, {i, 0, m - 1}, {k, 0, kmax}], 1], InfV[inf, #[[1]], #[[2]]] <= Dmax &];
  If[basis === {}, Return[None]];
  Do[{M, Bs} = af; dM = Exponent[M, g];
    cols = Table[r = PolynomialRemainder[Expand[g^bs[[1]] Bs[[bs[[2]] + 1]]], M, g];
                 PadRight[CoefficientList[r, g], dM], {bs, basis}];
    rows = Join[rows, Transpose[cols]], {af, affine}];
  Do[AppendTo[rows, Table[InfCoef[inf, j, bs[[1]], bs[[2]], r], {bs, basis}]],
    {j, 0, s - 1}, {r, -Dmax, zorders[[j + 1]] - 1}];
  If[rows === {}, Return[None]];
  ns = NullSpace[rows, ZeroTest -> (RootReduce[Together[#]] === 0 &)];
  If[ns === {}, Return[None]];
  vec = RRad /@ First[ns];
  a = ConstantArray[0, m];
  Do[a[[basis[[idx, 2]] + 1]] += vec[[idx]] g^basis[[idx, 1]], {idx, Length[basis]}];
  Collect[Expand[#], g, RRad] & /@ a];

(* RealiseClassGeneral[T, g, supp, c, mmax, verbose]
   Realisation of the residue class c for m >= 3 (the counterpart of
   RealiseClass): supp = {{p_i, Y0_i}} lists the primes and sheet polynomials
   (Y0_i^m = q mod p_i) of the places carrying c; the logand u in O with
   div(u) = N D_c - Sum_j n_j oo_j is sought for N = 1, ..., mmax and every
   distribution (n_j) of the N deg(D_c) poles over the s places at infinity,
   with the sheets lifted to Y_i^m = q mod p_i^N by Newton's iteration.  The
   norm N(u) = const prod p_i^N certifies the divisor.  Returns {c/N, u} or
   None.                                                                    *)
RealiseClassGeneral[T_, g_, supp_, c_, mmax_: 12, verbose_: False] := Catch[Module[
  {q = T["q"], m = T["m"], d, s, e, degSupp, lifts, inf = None, infL = -1, n, lmax, affine, u, Nu, Mall, quo, rmd},
  d = Exponent[q, g]; s = GCD[m, d]; e = m/s;
  degSupp = Total[Exponent[#[[1]], g] & /@ supp];
  lifts = {#[[1]], #[[1]], PolynomialRemainder[#[[2]], #[[1]], g]} & /@ supp;      (* {p_i, p_i^N, Y_i} *)
  Do[
    If[N0 > 1,
      lifts = Function[{l}, Module[{Pc = l[[1]], M2 = Expand[l[[2]] l[[1]]], Yk = l[[3]], inv},
          inv = InvMod[PolynomialRemainder[Expand[m Yk^(m - 1)], M2, g], M2, g];
          Yk = Collect[PolynomialRemainder[Expand[Yk - (Yk^m - q) inv], M2, g], g, RRad];
          {Pc, M2, Yk}]] /@ lifts];
    affine = Function[{l}, Module[{M = l[[2]], Yk = l[[3]]},
        {M, Table[Collect[PolynomialRemainder[Expand[Yk^i InvMod[PolynomialRemainder[T["E"][[i + 1]], M, g], M, g]], M, g], g, RRad], {i, 0, m - 1}]}]] /@ lifts;
    n = N0 degSupp;
    lmax = Floor[2 n/e] + 2;
    If[inf === None || lmax > infL, inf = InfinityData[T, g, lmax]; infL = lmax];
    Do[u = FindElement[T, g, affine, -cfg, inf];
      If[u === None, Continue[]];
      Nu = Expand[TNorm[T, u]];
      Mall = Expand[Times @@ (#[[2]] & /@ lifts)];
      rmd = PolynomialRemainder[Nu, Mall, g]; quo = Expand[PolynomialQuotient[Nu, Mall, g]];
      If[PolyZeroQ[Nu, g] || ! PolyZeroQ[rmd, g] || Exponent[quo, g] > 0, Continue[]];
      If[verbose, Print["      residue class ", c, ": div(u) = ", N0, " D_c - ", cfg, " oo; u = ", ToY[T, u, Global`y], ", N(u) = ", RR[quo], " * (", Times @@ (#[[1]] & /@ supp), ")^", N0, "; coefficient ", RRad[c/N0]]];
      Throw[{RRad[c/N0], u}, "rcg"],
      {cfg, SortBy[CompositionsOf[n, s], {Max[#], #} &]}],
    {N0, 1, mmax}];
  None], "rcg"];
RealiseClassAny[T_, g_, supp_, c_, mmax_: 12, verbose_: False] :=
  If[T["m"] == 2, RealiseClass[T, g, supp, c, mmax, verbose], RealiseClassGeneral[T, g, supp, c, mmax, verbose]];

(* UnitsGeneral[T, g, B, verbose]
   Unit candidates of O for m >= 3 (Algorithm 4, Step 6; the continued
   fraction of Part I is m = 2 only): the units are the elements whose
   divisor Sum_j z_j oo_j is supported at the s places at infinity, with
   Sum_j z_j = 0.  Every such divisor of degree <= B is tried by FindElement;
   a bounded height search, not complete (no non-torsion certificate for
   m >= 3).  For s = 1 there are no non-constant units.  Returns the list of
   unit tuples found (multiples of a found divisor are skipped).            *)
(* InLattice[z, basis, bound]: is the integer vector z a Z-combination of the
   basis vectors with coefficients of absolute value <= bound (a product of
   powers of units already found)? *)
InLattice[z_, basis_, bound_] := basis =!= {} &&
  AnyTrue[Tuples[Range[-bound, bound], Length[basis]], (# . basis) === z &];
UnitsGeneral[T_, g_, B0_: None, verbose_: False] := Module[{q = T["q"], m = T["m"], d, s, e, B = B0, inf, found = {}, z, u, Nu},
  d = Exponent[q, g]; s = GCD[m, d];
  If[s == 1, Return[{}]];
  e = m/s;
  If[B === None, B = If[s <= 3, 4, 3]];
  inf = InfinityData[T, g, Floor[2 B/e] + 2];
  Do[Do[If[AnyTrue[Transpose[{pos, neg}], #[[1]] > 0 && #[[2]] > 0 &], Continue[]];
      z = pos - neg;
      If[InLattice[z, found[[All, 1]], deg], Continue[]];
      u = FindElement[T, g, {}, z, inf];
      If[u === None, Continue[]];
      Nu = Expand[TNorm[T, u]];
      If[PolyZeroQ[Nu, g] || Exponent[Nu, g] != 0, Continue[]];
      AppendTo[found, {z, u}];
      If[verbose, Print["  unit candidate (m = ", m, "): ", ToY[T, u, Global`y], " with divisor ", z, " on the ", s, " places at infinity, N(u) = ", RR[Nu]]],
      {pos, CompositionsOf[deg, s]}, {neg, CompositionsOf[deg, s]}],
    {deg, 1, B}];
  found[[All, 2]]];

(* RealiseAtPoints[T, p, pts, taus, verbose]
   m >= 3, explicit places with constant coordinates over p: the residue
   classes are the distinct nonzero values c, each supported on the places
   (rho, y0) with residue c, i.e. the primes g - rho with the constant sheet
   y0; realised by RealiseClassGeneral.  Returns {{coefficient, u}, ...} or
   None.                                                                    *)
RealiseAtPoints[T_, p_, pts_, taus_, verbose_: False] := Module[{g = pts[[1, 1]], groups = {}, found, out = {}, supp, got},
  If[AnyTrue[pts, ! FreeQ[{#[[2]], #[[3]]}, Alternatives @@ T["gens"]] &], Return[None]];
  Do[If[IsZero[pr[[2]]], Continue[]];
    found = False;
    Do[If[IsZero[groups[[i, 1]] - pr[[2]]], AppendTo[groups[[i, 2]], pr[[1]]]; found = True; Break[]], {i, Length[groups]}];
    If[! found, AppendTo[groups, {pr[[2]], {pr[[1]]}}]],
    {pr, Thread[{pts, taus}]}];
  Do[supp = {Expand[g - #[[2]]], #[[3]]} & /@ grp[[2]];
    got = RealiseClassGeneral[T, g, supp, grp[[1]], 12, verbose];
    If[got === None, Return[None, Module]];
    AppendTo[out, got], {grp, groups}];
  out];

(* ---------------------------------------------------------- automatic towers *)

(* RadRules[q, Y]
   Rewriting rules sending every half-integral power of the radicand q to the
   basis element Y: q^(k + 1/2) -> q^k Y (so 1/Sqrt[q] -> Y/q).            *)
RadRules[q_, Y_, m_: 2] := If[m == 2, {Power[q, r_Rational] /; Denominator[r] == 2 :> q^(r - 1/2) Y},
  {Power[q, r_Rational] /; Divisible[m, Denominator[r]] :> Y^(r m)}];

(* ToPair[expr, Y, q]
   The pair {a, b} = a + b y of an expression rational in the generators and
   Y, Y^2 = q: numerator and denominator are reduced modulo Y^2 - q, the
   denominator is rationalised by its conjugate, and the two coefficients are
   read off.  For q = None the expression must be free of Y and the pair is
   {expr, 0}.                                                                *)
ToPair[expr_, Y_, q_] := Module[{e, n, d, dbar, nn, dd, red},
  If[q === None, Return[{Can[expr], 0}]];
  e = Together[expr];
  red[z_] := PolynomialRemainder[Expand[z], Y^2 - q, Y];
  n = red[Numerator[e]]; d = red[Denominator[e]];
  dbar = d /. Y -> -Y;
  nn = red[Expand[n dbar]]; dd = red[Expand[d dbar]];
  {Can[Coefficient[nn, Y, 0]/dd], Can[Coefficient[nn, Y, 1]/dd]}];

(* RationalFunctionQ[a, vars]
   True if a is a rational function of vars (numerator and denominator
   polynomials in vars, with any constant coefficients): the positive test
   for "innermost" subexpressions during the construction.  A negative test
   (free of Log, Tan, ...) would let an unknown function such as ArcTan[x]
   through as an opaque atom, which the linear solver would then treat as a
   constant and return a wrong integral.                                    *)
RationalFunctionQ[a_, vars_] := With[{e = Together[a]}, PolynomialQ[Numerator[e], vars] && PolynomialQ[Denominator[e], vars]];

(* Generator kinds handled by BuildTower: each h[a] becomes a generator t
   with the derivative pair listed, D a computed by the tower so far.
     Log      D a / a                    (primitive)
     ArcTan   D a / (1 + a^2)            (primitive)
     ArcTanh  D a / (1 - a^2)            (primitive)
     ArcSin   D a / Sqrt[1 - a^2]        (primitive through the radical y^2 = 1 - a^2)
     ArcCos  -D a / Sqrt[1 - a^2]
     ArcSinh  D a / Sqrt[1 + a^2]
     ArcCosh  D a / Sqrt[a^2 - 1]
     Exp      t D a                      (hyperexponential)
     Tan     (1 + t^2) D a               (hypertangent)
     Cot    -(1 + t^2) D a
     Tanh    (1 - t^2) D a
     Coth    (1 - t^2) D a
   Sin, Cos, Sinh, Cosh, Sec, Csc are not monomials and are rejected: rewrite
   them through Tan or Exp first.                                            *)
GenHeads = {Log, Tan, Cot, Tanh, Coth, ArcTan, ArcCot, ArcTanh, ArcCoth,
            ArcSin, ArcCos, ArcSinh, ArcCosh, ArcSec, ArcCsc, ArcSech, ArcCsch};
(* inverse functions with a rational derivative:  D t = coef(a) D a *)
InvRational = <|ArcTan -> Function[a, 1/(1 + a^2)], ArcCot -> Function[a, -1/(1 + a^2)],
                ArcTanh -> Function[a, 1/(1 - a^2)], ArcCoth -> Function[a, 1/(1 - a^2)]|>;
(* inverse functions with a derivative through a radical: D t = coef(a) D a / Sqrt[rad(a)];
   the third entry marks a derivative carrying |a|, whose sign the sample point decides *)
InvRadical = <|ArcSin -> {Function[a, 1 - a^2], Function[a, 1], False},
               ArcCos -> {Function[a, 1 - a^2], Function[a, -1], False},
               ArcSinh -> {Function[a, 1 + a^2], Function[a, 1], False},
               ArcCosh -> {Function[a, a^2 - 1], Function[a, 1], False},
               ArcSec -> {Function[a, a^2 - 1], Function[a, 1/a], True},
               ArcCsc -> {Function[a, a^2 - 1], Function[a, -1/a], True},
               ArcSech -> {Function[a, 1 - a^2], Function[a, -1/a], False},
               ArcCsch -> {Function[a, 1 + a^2], Function[a, -1/a], True}|>;
TrigHeads = {Sin, Cos, Sec, Csc}; HypHeads = {Sinh, Cosh, Sech, Csch};
BuildTower::function = "The function `1` is not a supported generator.";
BuildTower::radical2 = "The derivative of `1` needs the radical Sqrt[`2`], but the tower already has y^2 = `3`.";

(* SamplePoint[expr, x]
   A rational x0 at which every radicand and inverse-function argument of
   the integrand lies in its real domain: the point at which the signs of
   the square-root normalisations are decided.  None if no candidate works. *)
SamplePoint[expr_, x_] := Module[{rads, okay},
  rads = Cases[{expr}, Power[b_, e_Rational] :> {b, Denominator[e]}, Infinity];
  okay[x0_] := Quiet[Check[Module[{v},
    (* even roots need a positive radicand; odd roots a nonzero real one *)
    And @@ (With[{b = N[#[[1]] /. x -> x0], mr = #[[2]]}, NumericQ[b] && Abs[Im[b]] < 10^-12 && If[EvenQ[mr], Re[b] > 0, Re[b] != 0]] & /@ rads) &&
    And @@ Cases[{expr}, h_[a_] /; MemberQ[Join[GenHeads, TrigHeads, HypHeads], h] :>
      With[{av = N[a /. x -> x0]}, NumericQ[av] && Abs[Im[av]] < 10^-12 &&
        Switch[h, ArcSin | ArcCos | ArcTanh, -1 < Re[av] < 1, ArcCosh, Re[av] > 1, ArcCoth | ArcSec | ArcCsc, Abs[Re[av]] > 1,
                  ArcSech, 0 < Re[av] < 1, ArcCsch, Re[av] != 0, Log, Re[av] > 0, _, True]], Infinity] &&
    (v = N[expr /. x -> x0]; NumericQ[v] && Abs[Im[v]] < 10^-9 && Abs[v] < 10^12)], False]];
  SelectFirst[{1/2, 2, 1/3, 3, 3/2, 1/5, 5, -1/2, -2, 7/10, 3/10}, okay, None]];

(* BuildTower[integrand, x]  -- see the usage message. *)
BuildTower[integrand_, x_Symbol] := Module[
  {expr = integrand, gens = {x}, derivs = {{1, 0}}, q = None, m = 2, Y, back = {}, T, cands, c, kind, arg,
   tnew, ap, Da, base, r, mm, g, coef, rest, u, Dbase, pos, k = 0, fpair, rad, eulerStep, gq, a2, b2, c2, alpha, w, gw, yw, Dwp,
   trigStep, rootNormalize, trigC, wSurface, g0e, y0e, rcoef, sf, noRadical, sample, record, signFix, x0, icoef, toTuple, degs},
  Y = Unique["Y"];
  T = Tower[gens, derivs, q, m];
  toTuple[a_] := FromY[T, a /. RadRules[q, Y, m], Y];
  (* the sample point: every sign choice in a square-root normalisation is
     made so that the formal identity holds there, i.e. on the integrand's
     real domain (a monic or "positive generators" convention is wrong on one
     side of a branch cut, and flips the sign of the whole result) *)
  x0 = SamplePoint[integrand, x];
  sample = If[x0 === None, <||>, <|x -> x0|>];
  record[sym_, surface_] := If[Length[sample] > 0,
    With[{v = Quiet[N[surface /. sample]]},
      If[NumericQ[v] && Abs[Im[v]] < 10^-12, sample[sym] = Re[v]]];
    If[q =!= None && ! KeyExistsQ[sample, Y], With[{qv = Quiet[N[q /. sample]]}, If[NumericQ[qv], sample[Y] = qv^(1/m)]]]];
  signFix[sq_] := If[Length[sample] > 0,
    With[{v = Quiet[N[sq /. sample]]}, If[NumericQ[v] && Abs[Im[v]] < 10^-12 && Re[v] < 0, -sq, sq]], sq];
  (* Euler parametrisation: when a second radical is demanded of a conic
     y^2 = a g^2 + b g + c with a a square, substitute y = alpha g + w
     (alpha^2 = a), so that g = (w^2 - c)/(b - 2 alpha w) and y are rational
     in the new generator w and the curve disappears; the pending radical is
     then re-examined over Q(w).  The parametric form of Lemma 3.2. *)
  eulerStep[] := Module[{},
    If[m != 2, Return[False]];                       (* a conic y^2 = q only *)
    gq = SelectFirst[gens, PolynomialQ[q, #] && Exponent[q, #] == 2 && FreeQ[CoefficientList[q, #], Alternatives @@ gens] &, None];
    If[gq === None, Return[False]];
    {c2, b2, a2} = PadRight[CoefficientList[q, gq], 3];
    alpha = Sqrt[a2];
    w = Unique["w"];
    If[MatchQ[alpha, _Integer | _Rational],
      (* Euler's first substitution through the point at infinity, in the form
         y = -alpha g + w, i.e. w = y + alpha g: the LARGE root, positive for
         large positive g, so that the formal normalisations Sqrt[P/w^2] =
         Sqrt[P]/w hold on the natural real domain (w = y - alpha g would be
         negative there and flip the sign of the result) *)
      gw = Cancel[(w^2 - c2)/(b2 + 2 alpha w)]; yw = Cancel[w - alpha gw];
      Dwp = TowerD[T, {alpha gq, 1}];                             (* D(y + alpha g) as a pair *)
      wSurface = Sqrt[q] + alpha gq,
      (* otherwise through a rational point (g0, y0): y = y0 + w (g - g0) *)
      With[{pt = RationalPoint[q, gq]},
        If[pt === None, Return[False]];
        {g0e, y0e} = pt;
        gw = Cancel[(2 y0e w - (a2 + w^2) g0e - b2)/(a2 - w^2)]; yw = Cancel[y0e + w (gw - g0e)];
        Dwp = Pscale[1/(gq - g0e)^2, Padd[Pscale[gq - g0e, T["Dy"]], Pscale[-1, Pmul[{-y0e, 1}, derivs[[Position[gens, gq][[1, 1]]]], q]]]];
        wSurface = (Sqrt[q] - y0e)/(gq - g0e)]];
    Dwp = Cancel[(Dwp[[1]] + Dwp[[2]] Y) /. {Y -> yw, gq -> gw}];
    derivs = Cancel[(#[[1]] + #[[2]] Y) /. {Y -> yw, gq -> gw}] & /@ derivs;
    derivs = {#, 0} & /@ derivs;
    pos = Position[gens, gq][[1, 1]];
    gens = ReplacePart[gens, pos -> w]; derivs = ReplacePart[derivs, pos -> {Dwp, 0}];
    (* earlier back rules may mention this radical's symbol Y (e.g. a generator
       Log[g + Y]): record what Y was, then retire the symbol so that a later
       radical cannot be mistaken for it *)
    AppendTo[back, Y -> Sqrt[q]];
    AppendTo[back, w -> wSurface];
    record[w, wSurface]; sample = KeyDrop[sample, Y];
    expr = (expr /. RadRules[q, Y, m]) /. {Y -> yw, gq -> gw};
    q = None; Y = Unique["Y"]; T = Tower[gens, derivs, q, m];
    True];
  (* rootNormalize[base, k/mr]: a rational-function base num/den, or a
     polynomial base with mr-th-power factors: num den^(mr-1) = sq^mr sf with sf
     mr-th-power-free, base^(k/mr) = sq^k sf^(k/mr) / den^k; the signs of sq
     and den are chosen at the sample point so that the identity holds for the
     principal branches (for odd mr a negative sq is moved into sf as
     (-sq)^mr (-sf), and den is made positive).  The base and the exponent are
     passed separately: Power[base, 1/2] would be evaluated before the function
     sees it (Sqrt[a/4] -> Sqrt[a]/2 leaves a Times, not a Power). *)
  rootNormalize[bs_, ex_] := Module[{mr, num, den, sfl, sq = 1, sf = 1, pure, rest, const, sqf, dsgn},
    mr = Denominator[ex];
    {num, den} = {Numerator[Together[bs]], Denominator[Together[bs]]};
    If[OddQ[mr] && signFix[den] =!= den, {num, den} = {-num, -den}];
    sfl = FactorSquareFreeList[Expand[num den^(mr - 1)]];
    Do[sq *= fc[[1]]^Quotient[fc[[2]], mr]; sf *= fc[[1]]^Mod[fc[[2]], mr], {fc, sfl}];
    const = Cancel[Expand[num den^(mr - 1)]/(sq^mr sf)];
    If[! FreeQ[const, x] || ! FreeQ[const, Alternatives @@ Join[gens, {Y}]], Return[bs^ex]];
    sqf = signFix[sq];
    If[sqf =!= sq && OddQ[mr], const = -const];
    sq = sqf;
    dsgn = If[EvenQ[mr] && signFix[den] =!= den, (-1)^Numerator[ex], 1];
    (* pure-root factors: a factor of the squarefree radicand that is itself a
       generator g, positive at the sample point, is split off as g^ex so that
       it flattens (Lemma 3.2): Sqrt[2 t (1 + t^2)] = Sqrt[2] Sqrt[t] Sqrt[1 + t^2].
       This split is branch-unsafe (Sqrt[a b] = Sqrt[a] Sqrt[b] holds only for
       a, b > 0), so it is done ONLY when it is actually needed to reduce the
       radicand to a handled curve: either a simple radical y^m = q already exists
       (q =!= None -- keeping the factor fused would demand a second radical), or
       the fused squarefree radicand is not yet a genus-0 conic (total degree >= 3
       in the generators, so pulling a linear generator factor lowers the genus, as
       in the cubic example above).  When the fused radicand is already a genus-0
       conic (total degree <= 2) it can serve directly as the simple radical
       y^2 = q; keeping it fused then makes the back-substitution y -> Sqrt[q]
       globally faithful -- e.g. Sqrt[Sqrt[x] + x] stays fused instead of splitting
       into the (x > 0)-only x^(1/4) Sqrt[1 + Sqrt[x]]. *)
    pure = 1;
    If[q === None && With[{dm = Unique["d"]}, Exponent[Expand[sf /. Thread[gens -> gens dm]], dm]] <= 2,
      rest = const sf,                 (* genus-0 conic: keep the radical fused *)
      rest = const;
      Do[If[fc[[2]] == 1 && MemberQ[gens, fc[[1]]] && signFix[fc[[1]]] === fc[[1]],
          pure *= fc[[1]]^ex, rest *= fc[[1]]^fc[[2]]],
        {fc, FactorList[sf]}]];
    sf = Expand[rest];                 (* expanded, as before the split *)
    If[FreeQ[sf, Alternatives @@ Join[gens, {Y}]], Return[Together[dsgn sq^(mr ex) sf^ex/den^(mr ex)] pure]];
    Together[dsgn sq^(mr ex)/den^(mr ex)] pure Power[sf, ex]];
  (* trigStep[a, hyper]: rewrite every (hyperbolic) trigonometric function of
     the argument a through a tangent.  With s = sin a, c = cos a the
     integrand is a rational function of s, c and radicals of such; the
     substitution s -> t r, c -> r with r^2 = 1/(1 + t^2) (t = tan a) leaves
     an even integrand rational in t; otherwise the half-angle t = tan(a/2),
     s = 2t/(1 + t^2), c = (1 - t^2)/(1 + t^2), always does.  Hyperbolic:
     t = tanh a with r^2 = 1/(1 - t^2), or t = tanh(a/2), s = 2t/(1 - t^2),
     c = (1 + t^2)/(1 - t^2).  Tangents are preferred to exponentials: they
     keep the field real and their specials 1 + t^2, 1 - t^2 are the ones
     the theory handles. *)
  trigStep[a_, hyper_] := Module[{s, c, rules, e1, tn, rr, qr, bases, ph, evenQ, redB, e2, half, Da0, Dt,
      canChange, rr2, basesS, basesC, outerS, outerC, uNew},
    s = Unique["s"]; c = Unique["c"]; tn = Unique["t"]; rr = Unique["r"];
    rules = If[! hyper,
      {Sin[a] -> s, Cos[a] -> c, Tan[a] -> s/c, Cot[a] -> c/s, Sec[a] -> 1/c, Csc[a] -> 1/s},
      {Sinh[a] -> s, Cosh[a] -> c, Tanh[a] -> s/c, Coth[a] -> c/s, Sech[a] -> 1/c, Csch[a] -> 1/s}];
    e1 = expr /. rules;
    qr = If[! hyper, 1/(1 + tn^2), 1/(1 - tn^2)];
    (* radicals whose base involves s or c: reduce the bases separately *)
    bases = DeleteDuplicates[Cases[{e1}, Power[b_, _Rational] /; ! FreeQ[b, s | c] :> b, Infinity]];
    ph = Table[Unique["b"], {Length[bases]}];
    e2 = e1 /. Thread[bases -> ph];
    redB[z_] := ToPair[(z /. {s -> tn rr, c -> rr}), rr, qr];
    (* parity rules (change of the integration variable), when the trigonometric
       functions of x carry the whole x-dependence and no tower exists yet:
         odd in sin:  F = B(cos x) sin x,  int F dx = -int B(u) du,  u = cos x
         odd in cos:  F = B(sin x) cos x,  int F dx =  int B(u) du,  u = sin x
       (hyperbolic: u = cosh x with +, u = sinh x with +).  These keep the genus
       low where the half-angle substitution would inflate it. *)
    canChange = gens === {x} && a === x && FreeQ[e1, x] && q === None;
    If[canChange,
      rr2 = Unique["r"];
      basesS = ToPair[# /. s -> rr2, rr2, If[hyper, c^2 - 1, 1 - c^2]] & /@ bases;   (* s^2 = 1 - c^2, or c^2 - 1 *)
      outerS = ToPair[e2 /. s -> rr2, rr2, If[hyper, c^2 - 1, 1 - c^2]];
      If[AllTrue[basesS, #[[2]] === 0 &] && outerS[[1]] === 0,
        uNew = Unique["u"];
        expr = (If[hyper, 1, -1] outerS[[2]] /. c -> uNew) /. Thread[ph -> (basesS[[All, 1]] /. c -> uNew)];
        gens = {uNew}; derivs = {{1, 0}}; T = Tower[gens, derivs, None, m];
        AppendTo[back, uNew -> If[hyper, Cosh[a], Cos[a]]];
        record[uNew, If[hyper, Cosh[a], Cos[a]]];
        Return[Null, Module]];
      basesC = ToPair[# /. c -> rr2, rr2, 1 + If[hyper, 1, -1] s^2] & /@ bases;   (* c^2 = 1 - s^2, or 1 + s^2 *)
      outerC = ToPair[e2 /. c -> rr2, rr2, 1 + If[hyper, 1, -1] s^2];
      If[AllTrue[basesC, #[[2]] === 0 &] && outerC[[1]] === 0,
        uNew = Unique["u"];
        expr = (outerC[[2]] /. s -> uNew) /. Thread[ph -> (basesC[[All, 1]] /. s -> uNew)];
        gens = {uNew}; derivs = {{1, 0}}; T = Tower[gens, derivs, None, m];
        AppendTo[back, uNew -> If[hyper, Sinh[a], Sin[a]]];
        record[uNew, If[hyper, Sinh[a], Sin[a]]];
        Return[Null, Module]]];
    evenQ = AllTrue[Append[redB /@ bases, redB[e2]], #[[2]] === 0 &];
    Da0 = TowerD[T, toTuple[a]];
    If[evenQ,
      (* the outer expression, with the radical bases as placeholders, reduces
         to a rational function of t; each base reduces likewise and is put
         back under its root (normalised later by sqrtNormalize) *)
      expr = redB[e2][[1]] /. MapThread[#1 -> #2 &, {ph, redB[#][[1]] & /@ bases}];
      Dt = Pscale[If[! hyper, 1 + tn^2, 1 - tn^2], Da0];
      AppendTo[back, tn -> If[! hyper, Tan[a], Tanh[a]]]; record[tn, If[! hyper, Tan[a], Tanh[a]]],
      half = If[! hyper, {s -> 2 tn/(1 + tn^2), c -> (1 - tn^2)/(1 + tn^2)}, {s -> 2 tn/(1 - tn^2), c -> (1 + tn^2)/(1 - tn^2)}];
      expr = e1 /. half;
      Dt = Pscale[If[! hyper, (1 + tn^2)/2, (1 - tn^2)/2], Da0];
      AppendTo[back, tn -> If[! hyper, Tan[a/2], Tanh[a/2]]]; record[tn, If[! hyper, Tan[a/2], Tanh[a/2]]]];
    AppendTo[gens, tn]; AppendTo[derivs, Dt];
    T = Tower[gens, derivs, q, m]];
  (* general powers a^b with a non-rational exponent are Exp[b Log[a]] *)
  expr = expr /. Power[b_, e_] /; ! FreeQ[e, x] && b =!= E && ! MatchQ[e, _Rational | _Integer] :> Exp[e Log[b]];
  While[True,
    (* trigonometric families first: all Sin/Cos/Sec/Csc (or hyperbolic) of one argument together *)
    trigC = Cases[{expr}, (h_[a_] /; MemberQ[Join[TrigHeads, HypHeads], h]) /; RationalFunctionQ[a /. RadRules[q, Y, m], Join[gens, {Y}]] && ! FreeQ[a, Alternatives @@ Join[gens, {Y}]], Infinity];
    If[trigC =!= {},
      With[{c0 = First[SortBy[trigC, LeafCount]]}, trigStep[c0[[1]], MemberQ[HypHeads, Head[c0]]]];
      Continue[]];
    cands = Cases[{expr}, ((h_[a_] /; MemberQ[GenHeads, h]) | Power[E, a_] | Power[a_, _Rational]) /;
        RationalFunctionQ[a /. RadRules[q, Y, m], Join[gens, {Y}]] && ! FreeQ[a, Alternatives @@ Join[gens, {Y}]], Infinity];
    If[cands === {}, Break[]];
    c = First[SortBy[cands, LeafCount]];
    Switch[c,
      Power[E, _] | (h_[_] /; MemberQ[GenHeads, h]),
        arg = c[[If[Head[c] === Power, 2, 1]]];
        ap = toTuple[arg];
        Da = TowerD[T, ap];
        k++; tnew = Unique["t"];
        (* primitives through a radical: the radical must be the tower's *)
        rcoef = 1; noRadical = False; icoef = 1;   (* 1/Sqrt[rad] = (1/rcoef) y/q after normalisation *)
        If[KeyExistsQ[InvRadical, Head[c]],
          rad = Together[InvRadical[Head[c]][[1]][arg] /. RadRules[q, Y, m]];
          icoef = InvRadical[Head[c]][[2]][arg];
          If[InvRadical[Head[c]][[3]] && Length[sample] > 0,    (* the derivative carries |a| *)
            With[{av = Quiet[N[(arg /. Y -> If[q === None, Y, Sqrt[q]]) /. sample]]},
              If[NumericQ[av] && Re[av] < 0, icoef = -icoef]]];
          (* normalise Sqrt[rad] = rcoef Sqrt[sf] with sf a squarefree polynomial *)
          With[{nrm = rootNormalize[rad, 1/2]},
            With[{rts = Cases[{nrm}, Power[b_, 1/2] /; ! FreeQ[b, Alternatives @@ Join[gens, {Y}]] :> b, Infinity]},
              If[rts =!= {},
                sf = First[rts]; rcoef = Cancel[nrm/Sqrt[sf]];
                If[q === None,
                  q = Expand[sf]; m = 2; T = Tower[gens, derivs, q, m]; Da = TowerD[T, ap]; expr = expr /. RadRules[q, Y, m],
                  If[Expand[sf - q] =!= 0,
                    If[eulerStep[], Continue[], Message[BuildTower::radical2, c, sf, q]; Throw[$Failed, "build"]]];
                  If[m != 2, Message[BuildTower::radical3, c, q, m]; Throw[$Failed, "build"]]],
                rcoef = nrm; noRadical = True]]]];
        AppendTo[gens, tnew];
        AppendTo[derivs, Switch[Head[c],
          Log,     TDiv[T, Da, ap],
          ArcTan | ArcCot | ArcTanh | ArcCoth,
                   TMul[T, Da, toTuple[InvRational[Head[c]][arg]]],
          ArcSin | ArcCos | ArcSinh | ArcCosh | ArcSec | ArcCsc | ArcSech | ArcCsch,
                   With[{cp = toTuple[icoef/rcoef]},
                     If[noRadical, TMul[T, Da, cp], TMul[T, TMul[T, Da, cp], {0, 1/q}]]],
          Tan,     Pscale[1 + tnew^2, Da],
          Cot,     Pscale[-(1 + tnew^2), Da],
          Tanh,    Pscale[1 - tnew^2, Da],
          Coth,    Pscale[1 - tnew^2, Da],
          Power,   Pscale[tnew, Da]]];
        AppendTo[back, tnew -> c];
        record[tnew, c /. Y -> q^(1/m)];
        expr = expr /. c -> tnew,
      Power[_, _Rational],
        With[{nrm = rootNormalize[c[[1]], c[[2]]]},
          If[nrm =!= c, expr = expr /. c -> nrm; Continue[]]];
        {base, r} = List @@ c; mm = Denominator[r];          (* the degree of this root *)
        (* flattenable root: base = coef g + rest, coef constant, rest free of g  (Lemma 3.2) *)
        g = SelectFirst[gens, PolynomialQ[base, #] && Exponent[base, #] == 1 && FreeQ[base, Y] &&
                              FreeQ[Coefficient[base, #, 1], Alternatives @@ gens] &, None];
        If[g =!= None,
          coef = Coefficient[base, g, 1]; rest = base /. g -> 0;
          u = Unique["u"];
          Dbase = TowerD[T, toTuple[base]];
          pos = Position[gens, g][[1, 1]];
          gens = ReplacePart[gens, pos -> u];
          derivs = ReplacePart[derivs, pos -> Pscale[1/(mm u^(mm - 1)), Dbase]];
          derivs = derivs /. g -> (u^mm - rest)/coef;
          If[q =!= None, q = Expand[q /. g -> (u^mm - rest)/coef];
            If[! FreeQ[q, Y], Message[BuildTower::radicand, q]; Throw[$Failed, "build"]]];
          AppendTo[back, u -> base^(1/mm)];
          record[u, (base /. Y -> q^(1/m))^(1/mm)];
          expr = expr /. Power[base, rr_Rational] :> u^(mm rr);   (* the root itself, before g is eliminated *)
          expr = expr /. g -> (u^mm - rest)/coef,
          (* otherwise: the simple radical y^m = base *)
          If[q =!= None && Expand[base - q] =!= 0,
            If[eulerStep[], Continue[], Message[BuildTower::radicals, q, base]; Throw[$Failed, "build"]]];
          If[q =!= None,
            (* the radicand of the tower under another root base^(k/mm), mm | m *)
            If[! Divisible[m, mm], Message[BuildTower::mixedroots, m, mm, q]; Throw[$Failed, "build"]];
            expr = expr /. RadRules[q, Y, m]; Continue[]];
          If[! PolynomialQ[base, gens], Message[BuildTower::radicand, base]; Throw[$Failed, "build"]];
          (* every root of this radicand in the integrand: m = lcm of their degrees *)
          degs = DeleteDuplicates[Cases[{expr}, Power[b_, rr_Rational] /; Expand[b - base] === 0 :> Denominator[rr], Infinity]];
          m = If[Length[degs] > 1, LCM @@ degs, mm];
          If[AnyTrue[Select[FactorSquareFreeList[Expand[base]], ! FreeQ[#[[1]], Alternatives @@ gens] &], #[[2]] >= m &],
            Message[BuildTower::powerfree, base, m]; Throw[$Failed, "build"]];
          q = base;
          If[m > 2, derivs = PadRight[#, m] & /@ derivs];
          expr = expr /. RadRules[q, Y, m];
          If[Length[sample] > 0, With[{qv = Quiet[N[q /. sample]]}, If[NumericQ[qv], sample[Y] = qv^(1/m)]]]]];
    T = Tower[gens, derivs, q, m]];
  expr = expr /. RadRules[q, Y, m];
  If[! RationalFunctionQ[expr, Join[gens, {Y}]] || (q === None && ! FreeQ[expr, Y]),
    With[{fs = Cases[{expr}, h_[args___] /; ! MemberQ[{Plus, Times, Power, Rational, Complex, List}, h] && ! FreeQ[{args}, Alternatives @@ gens] :> h, Infinity]},
      If[fs =!= {}, Message[BuildTower::function, First[fs]], Message[BuildTower::residual, expr]]];
    Throw[$Failed, "build"]];
  fpair = FromY[T, expr, Y];
  If[q =!= None, AppendTo[back, Y -> q^(1/m)]];
  {T, fpair, back, Y}];

(* ParallelIntegrateMixed[integrand, x, opts]  -- the surface entry point.
   Builds the tower with BuildTower, runs the pair form, and maps the result
   back to x with ReplaceRepeated over the back rules (generator symbols to
   their defining expressions, flattened roots to base^(1/m), Y to Sqrt[q]).
   A certificate or failure tuple is returned as is, with the same
   back-substitution applied so the prime and residue read in surface terms.
   With "Verify" -> True the surface result is differentiated and compared
   with the integrand numerically at three points; a mismatch prints a
   warning and still returns the result.                                     *)
ParallelIntegrateMixed[integrand_, x_Symbol, opts : OptionsPattern[]] := TimeConstrained[Module[{bt, T, fpair, back, Y, res, surf, ok, sgn},
  bt = Catch[BuildTower[integrand, x], "build"];
  If[bt === $Failed, Return[$Failed]];
  {T, fpair, back, Y} = bt;
  If[OptionValue["Verbose"],
    Print["  tower: generators ", T["gens"], " with D = ", T["derivs"],
      If[T["q"] =!= None, ", y^" <> ToString[T["m"]] <> " = " <> ToString[T["q"], InputForm], ""]];
    Print["  integrand: ", fpair]];
  res = ParallelIntegrateMixed[fpair, T, "Bounds" -> OptionValue["Bounds"], "Verbose" -> OptionValue["Verbose"],
                                "SplitSpecials" -> OptionValue["SplitSpecials"]];
  If[ListQ[res], Return[res //. back]];
  surf = res //. back;
  (* Verify-or-decline soundness gate, with branch resolution.  A returned
     antiderivative MUST satisfy D[surf] == integrand on the integrand's real
     domain, checked numerically at REAL in-domain points (the discipline of the
     external corpus).  But the parallel method integrates over y^2 = q and maps
     the formal radical back as the PRINCIPAL Sqrt[q], so the surface can be the
     true antiderivative only up to the branch of that radical:
       - a global sign (D[surf] == -integrand everywhere): the decomposition took
         the other root of y^2 = q, and -surf is the antiderivative (A11); or
       - a sign that flips across a branch cut of the radicand on the real line
         (D[surf] == +integrand on one component, == -integrand on the next --
         e.g. Sqrt[1-Sin^6] == |Cos| Sqrt[...], A34): then no single elementary
         expression is a global antiderivative and surf is the conventional
         principal-branch answer, valid on each connected component.
     Pin the branch numerically: sgn = +1/-1 chooses surf/-surf; sgn = 0 declines.
     Soundness is preserved -- acceptance requires D[surf] to equal EXACTLY
     +integrand or -integrand (to 1e-12) at every finite real sample; a genuinely
     wrong surface gives a residual that is neither, and is still rejected.
     (verify-timeout: accept as principal branch rather than reject a
     likely-correct but large answer.) *)
  fac = TimeConstrained[
    Quiet[Module[{cand, pts, dv, fv, keep, tol = 10^-12, units, gu, u, dsurf},
      dsurf = D[surf, x];   (* differentiate once, reused at every sample below *)
      (* A dense fixed spread of positive rationals, not just the original 11
         anchors.  The sparse grid could be silently passed by a wrong,
         over-complete surface that happened to agree at those 11 points (a
         reverted A40 over-completion was wrong only off-grid, e.g. at x = 1/4;
         A19's non-deterministic S'-unit selection could produce an r wrong off
         the anchors).  The error of such a surface is a rational function that
         agrees with the integrand only at its finitely many zeros, so a denser
         grid of varied-denominator points rejects it.  The set is FIXED, so the
         gate stays deterministic and has no RNG side effect on the caller's
         random state (SeedRandom is global and BlockRandom is unavailable).
         Soundness is unchanged: more samples can only reject more wrong
         surfaces, never a correct one -- each kept point still tests
         D[surf] == u f EXACTLY (to tol).  Hoisting D[surf,x] out of the map
         makes this denser gate net-faster than the old 11-point one, which
         re-differentiated surf once per point. *)
      cand = {1/2, 1/3, 2/3, 1/5, 4/5, 2, 3, 3/2, 5/2, 7/3, 1/7,
              1/4, 3/4, 1/6, 5/6, 2/7, 5/7, 6/7, 3/8, 5/8, 7/8,
              2/9, 4/9, 7/9, 3/11, 7/11, 9/11, 5/13, 11/13,
              4/3, 5/4, 7/4, 9/7};
      pts = Select[cand,
        With[{v = N[integrand /. x -> #, 30]},
          NumericQ[v] && Abs[Im[v]] < 10^-10 && Abs[v] < 10^12] &];
      If[pts === {}, pts = {7/13 + 5 I/11, 4/9 - 3 I/7, 2/7 + 3 I/5}];   (* complex fallback *)
      dv = (N[dsurf /. x -> #, 30]) & /@ pts;
      fv = (N[integrand /. x -> #, 30]) & /@ pts;
      keep = Select[Range[Length[pts]],
        NumericQ[dv[[#]]] && NumericQ[fv[[#]]] && Abs[dv[[#]]] < 10^12 && Abs[fv[[#]]] < 10^12 &];
      (* The exact unit u with D[surf] == u f.  The formal radical y (y^2 = q) is
         mapped back as the PRINCIPAL Sqrt[q], so surf can be the true
         antiderivative up to a fourth root of unity: 1 or -1 for a real radical
         (the global-sign case A11, and the other root of y^2=q); and +-I when q
         is negative on the domain, where Sqrt[q] = I Sqrt[-q] leaves surf off by
         a factor of I (Charlwood A35).  Then D[surf/u] == f and surf/u =
         Conjugate[u] surf is the antiderivative.  Soundness: acceptance requires
         D[surf] to equal EXACTLY u f (to 1e-12) at every finite real sample for a
         single unit u; a genuinely wrong surface matches none. *)
      units = {1, -1, I, -I};
      gu = None;   (* the global unit u with D[surf] == u f, if any (# = keep index) *)
      Do[If[AllTrue[keep, Abs[dv[[#]] - u fv[[#]]] < tol &], gu = u; Break[]], {u, units}];
      Which[
        keep === {}, 0,                                     (* no usable sample *)
        gu =!= None, 1/gu,                                  (* global unit u: surf/u *)
        AllTrue[keep, Abs[dv[[#]] - fv[[#]]] < tol || Abs[dv[[#]] + fv[[#]]] < tol &],
          (* branch cut: D[surf] = +-f exactly at every point -> principal-branch
             answer, signed to match the sample nearest the origin *)
          With[{anchor = keep[[First[Ordering[Abs[N[pts[[keep]], 30]]]]]]},
           If[Abs[dv[[anchor]] - fv[[anchor]]] < tol, 1, -1]],
        True, 0]]],                                          (* not a unit multiple: wrong answer *)
    15, 1];
  If[fac === 0 || FreeQ[{1, -1, I, -I}, fac], Return[{"failed", "verification failed"}]];
  If[fac =!= 1, surf = fac surf];
  surf],
  $ParallelMixedTimeBudget, {"failed", "time budget exceeded"}];

(* Robustness backstop: a few integrands drive an internal stage (split-specials
   over Fbar + a number-field compositum, a two-generator ansatz blow-up, or a
   nested-radical BuildTower) into a very long run.  This coarse wall-clock
   budget converts such a run into a clean decline instead of an unbounded hang
   (interactive Integrate has no timeout of its own).  Set above the slowest
   legitimate solve observed on the stress corpus (~34 s) with headroom; it is a
   backstop, not a per-stage guard -- deterministic ansatz-size / number-field
   degree caps are follow-up work.  Overridable by the caller if needed. *)
$ParallelMixedTimeBudget = 45;

Options[ParallelIntegrateMixed] = {"Bounds" -> None, "Verbose" -> False, "Verify" -> False, "SplitSpecials" -> Automatic, "SpecialExponent" -> 0};
(* iPIM shares these options via OptionsPattern[ParallelIntegrateMixed]; an
   unpassed OptionValue now resolves its default against ParallelIntegrateMixed
   (the symbol named in OptionsPattern[...]), so a separate Options[iPIM]
   declaration is no longer needed (was A9, v0.170). *)

(* pair form: with "SplitSpecials" -> Automatic a failed system is retried with
   every special prime split into its linear factors over the algebraic
   closure (Theorem 6.1 describes the specials over Fbar; the integral of
   Sqrt[Tan[x]] needs distinct coefficients on the four factors of 1 + u^4). *)
ParallelIntegrateMixed[f00_List, T0_Association, opts : OptionsPattern[]] := Module[
  {f0 = f00, T = T0, r, split = OptionValue["SplitSpecials"], se, gq, isConic, isQuartic, res, T2, f2, back, r2, failedQ, g0, r0},
  failedQ[rr_] := ListQ[rr] && MemberQ[{"failed", "needs torsion realisation (milestone iii)"}, rr[[1]]];
  (* 0. a single generator g over the curve with Dg = r(g) != 1 -- a flattened
        root (Lemma 3.2) or the parameter of a parametrised conic -- is rescaled
        to d/dg (Lemma 3.4): int f dx = int (f/r) dg with the same antiderivative,
        and the tower is the setting of Part I (exact bounds, Prop. 9.2(b)) *)
  If[Length[T["gens"]] == 1 && T["q"] =!= None,
    g0 = T["gens"][[1]]; r0 = T["derivs"][[1, 1]];
    If[r0 =!= 1 && MatchQ[Rest[T["derivs"][[1]]], {0 ...}] && RationalFunctionQ[r0, {g0}] &&
        Complement[Cases[{r0}, _Symbol, Infinity], {g0}] === {},
      If[OptionValue["Verbose"], Print["  single generator ", g0, " with D", g0, " = ", r0, ": rescaled to d/d", g0, " (Lemma 3.4)"]];
      f0 = Can[#/r0] & /@ f0;
      T = Tower[{g0}, {{1}}, T["q"], T["m"]]]];
  gq = If[T["q"] === None, None, SelectFirst[T["gens"], FreeQ[T["q"], Alternatives @@ DeleteCases[T["gens"], #]] &, None]];
  isConic = T["m"] == 2 && gq =!= None && Exponent[T["q"], gq] == 2 && FreeQ[CoefficientList[T["q"], gq], Alternatives @@ T["gens"]];
  isQuartic = T["m"] == 2 && gq =!= None && Exponent[T["q"], gq] == 4 && FreeQ[CoefficientList[T["q"], gq], Alternatives @@ T["gens"]] &&
    MatchQ[Sqrt[Coefficient[T["q"], gq, 4]], _Integer | _Rational];
  (* 1. the base run, with the tower specials over Q *)
  r = Catch[iPIM[f0, T, "SplitSpecials" -> (split === True), opts], "PIM"];
  (* 2. the exponents of tower specials absent from the denominator are a
        guessed input (Remark 9.1): retry with exponents 1 and 2 *)
  Do[If[ListQ[r] && r[[1]] === "failed" && r[[2]] === "no solution within bounds" && MemberQ[r, "guess"] && OptionValue["Bounds"] === None,
      r = Catch[iPIM[f0, T, "SplitSpecials" -> (split === True), "SpecialExponent" -> se, opts], "PIM"]],
    {se, {1, 2}}];
  (* 3. a conic with a rational point: parametrise it away BEFORE any work
        over the algebraic closure -- on a genus-0 curve the S'-unit groups
        can have large rank, and the split-specials search over the curve
        produces nested-radical coefficients the linear algebra cannot digest *)
  If[failedQ[r] && isConic && OptionValue["Bounds"] === None,
    res = ConicToLine[T, f0];
    If[res =!= None,
      {T2, f2, back} = res;
      If[OptionValue["Verbose"], Print["  conic radicand: parametrised by ", back[[1, 1]], " = ", back[[1, 2]], "; the tower is now transcendental"]];
      r2 = ParallelIntegrateMixed[f2, T2, "Verbose" -> OptionValue["Verbose"], "SplitSpecials" -> split];
      If[! ListQ[r2], Return[r2 /. back, Module], If[r2[[1]] === "not elementary", Return[r2, Module]]]]];
  (* 4. specials split over the algebraic closure, with the exponent retries of 2 *)
  If[split === Automatic && ListQ[r] && r[[1]] === "failed" && r[[2]] === "no solution within bounds" && MemberQ[r, "splittable"],
    r = Catch[iPIM[f0, T, "SplitSpecials" -> True, opts], "PIM"];
    Do[If[ListQ[r] && r[[1]] === "failed" && r[[2]] === "no solution within bounds" && MemberQ[r, "guess"] && OptionValue["Bounds"] === None,
        r = Catch[iPIM[f0, T, "SplitSpecials" -> True, "SpecialExponent" -> se, opts], "PIM"]],
      {se, {1, 2}}]];
  r = If[ListQ[r] && r[[1]] === "failed", DeleteCases[r, "splittable" | "guess"], r];
  (* 5. a quartic radicand that blocks realisation: change to the cubic model *)
  If[failedQ[r] && isQuartic && OptionValue["Bounds"] === None,
    {T2, f2, back} = QuarticToCubic[T, f0];
    If[OptionValue["Verbose"], Print["  quartic radicand: changing to the cubic model y^2 = ", T2["q"]]];
    r2 = ParallelIntegrateMixed[f2, T2, "Verbose" -> OptionValue["Verbose"], "SplitSpecials" -> split];
    If[! ListQ[r2], r = (r2 /. Last[back]) /. First[back], If[r2[[1]] === "not elementary", r = r2]]];
  r];

(* ------------------------------------------------------------ the linear system over one number field
   The ansatz of iPIM is assembled and solved over ONE algebraic number field:
   the algebraic atoms of all the constants of the system (I, Sqrt[2], a nested
   radical, a Root object) are expressed in one field by a single ToNumberField
   call, the unknown-free columns of the system are polynomials with
   AlgebraicNumber coefficients (whose arithmetic is canonical, so that a zero
   is 0), the coefficient matrix is read off monomial by monomial and row
   reduced exactly, and the solution is mapped back to the radicals by linear
   algebra over Q.  Nothing here runs Together, Cancel or Solve on an
   expression that contains an unknown or a Root object; this is what SymPy's
   linsolve does with its QQ<alpha> domain and is the counterpart of
   _assemble_field / _FieldSystem of parallel_mixed.py.                      *)

NumDen[e_] := With[{t = Together[e]}, {Numerator[t], Denominator[t]}];

(* AlgAtoms[e]: the algebraic atoms of the constants of e (outermost first,
   so that a nested radical is one atom and is replaced as a whole) *)
AlgAtoms[e_] := Module[{c},
  c = Select[DeleteDuplicates[Cases[e, _Root | Power[_?NumericQ, _Rational], {0, Infinity}]],
    ! MatchQ[#, _Integer | _Rational] &];
  (* every Gaussian constant a + b I is I times rationals: one atom, I.
     Detect a Complex atom with Cases on the _Complex pattern, NOT
     FreeQ[e, Complex]: Mathilda's FreeQ with a bare type-symbol form does not
     match an atom's type-head (FreeQ[I + w, Complex] is True where WL is False),
     so a FreeQ[e, Complex] test misses I and the ansatz is assembled over Q with
     a bare Complex[0,1] and one spurious equation.

     Gated on $PMComplexField (only the split-specials retry sets it True): the
     assembly over Q(i) is 5-10x more FLINT multivariate-polynomial work than over
     Q, so an integrand whose Gaussian constants are incidental (a dsolve
     intermediate that solves over Q anyway -- e.g. y' == (Cos+1)/(2-Sin y))
     stays on the fast Q path, while the genuinely-algebraic split (Charlwood A1,
     A16, A19, A20, A37, whose specials factor over Q(i)) declines over Q, retries
     with the split, and there detects I and builds the field it needs.  The Q
     assembly is never WRONG when I is present -- it only DECLINES (the spurious
     equation) -- so trying it first is sound. *)
  If[Cases[e, _Complex, {0, Infinity}] =!= {}, AppendTo[c, I]];
  c];

(* FieldData[atoms]: None for no atoms; $Failed when the atoms do not generate
   a number field; otherwise {rules, back}: rules maps every atom to its
   AlgebraicNumber in one common field, back maps a field element to a linear
   combination over Q of monomials in the atoms (a basis of the field found by
   row reduction once), i.e. back to the original radicals *)
FieldData[atoms_] := Module[{nf, theta, n, zz, vec, degs, cands, vecs, red, piv, B, Binv, basisRad, back, rules},
  If[atoms === {}, Return[None]];
  nf = Quiet[Check[ToNumberField[atoms], $Failed]];
  If[! MatchQ[nf, {__AlgebraicNumber}], Return[$Failed]];
  theta = nf[[1, 1]];
  n = Exponent[MinimalPolynomial[theta, zz], zz];
  vec[e_] := Which[MatchQ[e, _AlgebraicNumber], PadRight[e[[2]], n],
                   MatchQ[e, _Integer | _Rational], PadRight[{e}, n], True, $Failed];
  (* the atom I stands for every Gaussian constant: Complex[a, b] -> a + b I *)
  rules = Join[MapThread[If[#1 === I, Complex[a_, b_] :> a + b #2, #1 -> #2] &, {atoms, nf}]];
  rules = Join[DeleteCases[rules, _RuleDelayed], Cases[rules, _RuleDelayed]];
  degs = Exponent[MinimalPolynomial[#, zz], zz] & /@ atoms;
  If[Times @@ degs > 8192,
    (* too many candidate monomials: the compact but opaque form *)
    back = Function[e, If[MatchQ[e, _AlgebraicNumber], RRad[RootReduce[e]], e]];
    Return[{rules, back}]];
  cands = Tuples[Range[0, # - 1] & /@ degs];
  vecs = vec[Times @@ (nf^#)] & /@ cands;
  red = RowReduce[Transpose[vecs]];
  piv = DeleteCases[LengthWhile[#, # === 0 &] + 1 & /@ red, _?(# > Length[cands] &)];
  B = vecs[[piv]]; Binv = Inverse[B];
  basisRad = (Times @@ (atoms^#)) & /@ cands[[piv]];
  back = Function[e, With[{v = vec[e]}, If[v === $Failed, e, (v . Binv) . basisRad]]];
  {rules, back}];

(* MonicPair[{num, den}, gens]: the fraction with den monic -- its leading
   coefficient (a rational or an AlgebraicNumber) moved into num -- so that
   denominators equal up to a constant factor coincide *)
MonicPair[{num_, den_}, gens_] := Module[{lc},
  If[num === 0, Return[{0, 1}]];
  lc = If[NumericQ[den], den, CoefficientRules[den, gens][[1, 2]]];
  If[lc === 1, {num, den}, {Expand[num/lc], Expand[den/lc]}]];

(* AnsatzSystem[T, rem, denv, units, unkLogs, css, monos, gammas, betas, unks]
   The linear system of the ansatz over one number field, solved: returns
   <|"sub" -> substitution for every unknown (the particular solution with the
   non-pivot unknowns 0 -- the unique such solution for the column order of
   unks, the same one parallel_mixed.py produces) or None when the system is
   inconsistent, "neq" -> number of equations|>, or $Failed when a constant of
   the system is not an algebraic number (the caller then takes the Expr
   route).  The columns: for m <= 2 every D(mono/D_v) and D(mono y/D_v) is a
   pair of polynomials over the one denominator Lfix = D_v^2 Dden 2q, by
   polynomial arithmetic alone (Dden a common multiple of the denominators of
   the D g); for m >= 3 the tower derivation on the Trager basis; the unit and
   special logands D u / u once each.                                          *)
AnsatzSystem[T_, rem_, denv_, units_, unkLogs_, css_, monos_, gammas_, betas_, unks_] := Module[
  {gens = T["gens"], q = T["q"], m = T["m"], nc = Length[css], cols, colFr, remFr, monoCols, Qh, atoms, fd, rules, back, alg,
   densG, Dden, Dg0s, Dg1s, denvF, ddenv, DdenF, qF, Qs0, Qs1, Lfix, parts, monoF, cg, n0, n1, s0, s1,
   colIdx, ncols, entries, rhs, nrows, rowOf, distinct, quo, dpos, poly, aug, red, xs, ok, sub, resid, partsM, remM},
  (* the unknown-free pieces of the system, in radicals *)
  cols = Join[
    Table[{gammas[[i]], TDiv[T, TowerD[T, units[[i, 1]]], units[[i, 1]]]}, {i, Length[units]}],
    Table[{betas[[i]], If[q =!= None,
        TDiv[T, TowerD[T, TScalar[T, unkLogs[[i, 1]]]], TScalar[T, unkLogs[[i, 1]]]],
        {Can[TowerD[T, {unkLogs[[i, 1]], 0}][[1]]/unkLogs[[i, 1]]], 0}]}, {i, Length[unkLogs]}]];
  colFr = {#[[1]], NumDen /@ Take[#[[2]], nc]} & /@ cols;
  remFr = NumDen /@ Take[rem, nc];
  monoCols = {};
  If[m >= 3,
    Do[With[{col = TowerD[T, (Times @@ (gens^monos[[a]]))/denv TUnit[T, i - 1]]},
        Do[AppendTo[monoCols, {css[[i, a]], j, NumDen[col[[j]]]}], {j, nc}]],
      {i, nc}, {a, Length[monos]}]];
  densG = DeleteDuplicates[Flatten[Table[Denominator[Together[c]], {d, T["derivs"]}, {c, d}]]];
  Qh = If[q =!= None, Dhat[T, q], None];
  (* the field *)
  atoms = AlgAtoms[{T["derivs"], q, denv, colFr, remFr, monoCols, Qh}];
  fd = FieldData[atoms];
  If[fd === $Failed, Return[$Failed]];
  {rules, back} = If[fd === None, {{}, Identity}, fd];
  alg = rules =!= {};
  (* the columns of the monomials *)
  parts = Table[{}, {nc}];                          (* per coordinate: {unknown, {num, den}} *)
  Dden = If[atoms === {} || FreeQ[densG, Alternatives @@ atoms], Fold[PolynomialLCM, 1, densG], Times @@ densG];
  denvF = Expand[denv /. rules];
  DdenF = Expand[Dden /. rules];
  If[m <= 2,
    Dg0s = Expand[Can[#[[1]] Dden] /. rules] & /@ T["derivs"];
    Dg1s = Expand[Can[#[[2]] Dden] /. rules] & /@ T["derivs"];
    ddenv = D[denvF, #] & /@ gens;
    If[q =!= None,
      qF = Expand[q /. rules];
      Qs0 = Expand[Can[Qh[[1]] Dden] /. rules]; Qs1 = Expand[Can[Qh[[2]] Dden] /. rules];
      Lfix = Expand[denvF^2 DdenF 2 qF],
      Lfix = Expand[denvF^2 DdenF]];
    Do[monoF = Times @@ (gens^monos[[a]]);
      cg = Table[Expand[D[monoF, gens[[k]]] denvF - monoF ddenv[[k]]], {k, Length[gens]}];
      n0 = Expand[cg . Dg0s]; n1 = Expand[cg . Dg1s];
      If[q =!= None, n0 = Expand[2 qF n0]; n1 = Expand[2 qF n1]];
      AppendTo[parts[[1]], {css[[1, a]], {n0, Lfix}}];
      If[nc > 1, AppendTo[parts[[2]], {css[[1, a]], {n1, Lfix}}]],
      {a, Length[monos]}];
    If[q =!= None,
      Do[monoF = Times @@ (gens^monos[[a]]);
        cg = Table[Expand[D[monoF, gens[[k]]] denvF - monoF ddenv[[k]]], {k, Length[gens]}];
        s0 = Expand[cg . Dg0s]; s1 = Expand[cg . Dg1s];
        (* D(m y/D_v) = c.(Dg0 + Dg1 y) y + (m/D_v)(Dy0 + Dy1 y),  y^2 = q *)
        n0 = Expand[2 qF qF s1 + monoF denvF Qs1 qF];
        n1 = Expand[2 qF s0 + monoF denvF Qs0];
        AppendTo[parts[[1]], {css[[2, a]], {n0, Lfix}}];
        AppendTo[parts[[2]], {css[[2, a]], {n1, Lfix}}],
        {a, Length[monos]}]],
    Do[AppendTo[parts[[mc[[2]]]], {mc[[1]], Expand[mc[[3]] /. rules]}], {mc, monoCols}]];
  Do[AppendTo[parts[[i]], {cf[[1]], Expand[cf[[2, i]] /. rules]}], {cf, colFr}, {i, nc}];
  remFr = Expand[remFr /. rules];
  (* the equations: per coordinate, every column times the quotient of the
     common denominator by its own (the lcm over Q, the product of the distinct
     denominators over an algebraic field), read off monomial by monomial *)
  colIdx = AssociationThread[unks -> Range[Length[unks]]];
  ncols = Length[unks];
  entries = {}; rhs = {}; nrows = 0; rowOf = <||>;
  Do[partsM = {#[[1]], MonicPair[#[[2]], gens]} & /@ parts[[i]];
    remM = MonicPair[remFr[[i]], gens];
    distinct = DeleteDuplicates[Join[partsM[[All, 2, 2]], {remM[[2]]}]];
    (* the common multiple of the distinct denominators and its quotients.
       Over Q (they are free of AlgebraicNumbers -- after the monic normalisation
       the usual case) use the minimal PolynomialLCM.  With AlgebraicNumber[theta]
       coefficients PolynomialLCM has NO number-field path -- it treats an
       AlgebraicNumber coefficient as opaque and returns a spurious factor -- and
       the old radical route (map back to radicals, Cancel[.., Extension ->
       Automatic]) could drop a factor of a special over the extension; both left
       an inconsistent, over-generated system (the Charlwood A2/A3/P8 divergence,
       MATHILDA_DIVERGENCES.md E.2.1).  Instead use the PRODUCT of the distinct
       denominators as the common multiple: any common multiple gives an
       equivalent system -- L is unknown-free, so multiplying the polynomial
       identity by an extra factor preserves its solution set -- and the quotient
       L/dd is the product of the OTHER denominators, exact by construction with
       no field gcd or division.  (The product is larger than the minimal LCM, so
       the system may carry redundant rows that row-reduce to 0 == 0; it stays
       consistent, which is what decides the solve.) *)
    quo = If[FreeQ[distinct, AlgebraicNumber],
      With[{L = Fold[PolynomialLCM, 1, distinct]}, Expand[Cancel[L/#]] & /@ distinct],
      With[{nd = Length[distinct]},
        Table[Expand[Times @@ Delete[distinct, k]], {k, nd}]]];
    dpos[dd_] := Position[distinct, _?(SameQ[#, dd] &), {1}, Heads -> False][[1, 1]];
    Do[If[part[[2, 1]] === 0, Continue[]];
      (* leave the product unexpanded: CoefficientRules expands it internally,
         and its native fmpq_mpoly read-off multiplies the factors natively and
         reads the terms straight out -- avoiding the big Expr Plus that the
         explicit Expand built and Orderless-sorted (the ans_eqn hot spot). *)
      poly = part[[2, 1]] quo[[dpos[part[[2, 2]]]]];
      Do[If[! KeyExistsQ[rowOf, {i, cr[[1]]}], nrows++; rowOf[{i, cr[[1]]}] = nrows];
        AppendTo[entries, {rowOf[{i, cr[[1]]}], colIdx[part[[1]]]} -> cr[[2]]],
        {cr, CoefficientRules[poly, gens]}],
      {part, partsM}];
    poly = remM[[1]] quo[[dpos[remM[[2]]]]];
    Do[If[! KeyExistsQ[rowOf, {i, cr[[1]]}], nrows++; rowOf[{i, cr[[1]]}] = nrows];
      AppendTo[rhs, {rowOf[{i, cr[[1]]}], ncols + 1} -> cr[[2]]],
      {cr, CoefficientRules[poly, gens]}],
    {i, nc}];
  If[nrows == 0, Return[<|"sub" -> Thread[unks -> 0], "neq" -> 0|>]];
  (* the solve: exact row reduction of the augmented matrix *)
  (* dense assembly from the {pos -> val} rules (SparseArray-free): each
     {row, col} is written once by the construction above, so a plain
     overwrite matches the reference (rows[r][c] = coeff in parallel_mixed.py). *)
  aug = ConstantArray[0, {nrows, ncols + 1}];
  Do[aug[[e[[1, 1]], e[[1, 2]]]] = e[[2]], {e, Join[entries, rhs]}];
  (* one-step row reduction: with AlgebraicNumber entries the default method
     is two orders of magnitude slower on the systems of the split specials *)
  red = RowReduce[aug, Method -> "OneStepRowReduction"];
  xs = ConstantArray[0, ncols]; ok = True;
  Do[With[{p = LengthWhile[red[[k]], # === 0 &] + 1},
      Which[p > ncols + 1, Null, p == ncols + 1, ok = False, True, xs[[p]] = red[[k, ncols + 1]]]],
    {k, nrows}];
  If[! ok, Return[<|"sub" -> None, "neq" -> nrows|>]];
  resid = aug[[All, 1 ;; ncols]] . xs - aug[[All, ncols + 1]];
  If[AnyTrue[resid, # =!= 0 &], Throw[{"failed", "internal: solution does not verify"}, "PIM"]];
  sub = Thread[unks -> (back /@ xs)];
  <|"sub" -> sub, "neq" -> nrows|>];


(* iPIM[f, T, opts]
   The body of ParallelIntegrateMixed; ParallelIntegrateMixed wraps it in
   Catch so that the early exits (certificates and failures, raised with
   Throw[..., "PIM"] from inside the prime loop) return from the routine and
   not merely from the loop.  Local variables mirror the Python names:
     detLogs -- realised logands {tau, {a, b}} with known coefficient tau;
     unkLogs -- special primes {p, 0} offered with an unknown coefficient;
     denv    -- the denominator D_v of the rational part (special s-parts and
                Hermite factors);
     rem     -- the residual f - Sum tau D u / u over the realised logands;
     EE      -- the pair D(V) - rem + Sum beta_i D s_i / s_i whose vanishing
                is the linear system.                                        *)
iPIM[f0_, T_, OptionsPattern[ParallelIntegrateMixed]] := Module[
  {gens = T["gens"], q = T["q"], m = T["m"], nc, css, verbose = OptionValue["Verbose"], bounds = OptionValue["Bounds"],
   Y, f, dlcm, detLogs = {}, unkLogs = {}, denv = 1, fl, p, mult, branch, eta, delta, special,
   eP, vP, Dp, tp, texpr, pts, taus, cert, done, seen, cand, rem, ld, nb, db, monos, cs0, cs1,
   V, EE, betas, eqs, unks, sol, sub, frees, y, surf, I0, split = OptionValue["SplitSpecials"], splittable, newLogs, PP, rts, torsion = {}, got, tinf, lower, nonconst, units = {}, unitsComplete = True, uu0, certd, gammas, B, r2, g0, sunits, sols, uuS, sexp = OptionValue["SpecialExponent"], gstar, a2, b2, c2, disc, s2, s, ok, pend, uu, tv,
   pendingClasses = {}, rc, gDir, classes, principal, groups, found, unrealised, sysK, neq, A, key, unitsBase, sunitsAll, spec, tau, cv, tauHi, exps, proved = False, curve, typeE},
  Y = Unique["y"];
  nc = If[q === None, 1, T["n"]];                      (* coordinates carrying the integrand *)
  f = TPad[T, f0]; f = Table[If[i > nc, 0, Can[f[[i]]]], {i, T["n"]}];
  (* Steps 1--14 of Algorithm 4 -- the classification, the residues and their
     realisation, the tower specials, the units and the residual -- are
     computed once per integrand and tower and reused by every rung of the
     retry ladder (special exponents, exact bounds, the split of the specials),
     which re-enters at Step 15; the counterpart of _Analysis in
     parallel_mixed.py.  The split variant of Steps 5--6 is memoised below. *)
  key = Hash[{f0, T}];
  A = Lookup[$analyses, key, None];
  unitsBase = {};
  If[A === None,
  (* the components must be rational functions of the generators: an opaque
     function would be treated as a constant by the solver *)
  If[! AllTrue[f, RationalFunctionQ[#, gens] &],
    Throw[{"failed", "integrand component is not a rational function of the generators", f}, "PIM"]];
  dlcm = PolynomialLCM @@ (Denominator /@ f);
  fl = Rest[FactorList[dlcm]];
  Do[{p, mult} = fac;
    If[FreeQ[p, Alternatives @@ gens], Continue[]];
    {branch, eta, delta, special} = ClassifyPrime[T, p];
    eP = If[branch, m, 1];
    vP = VP[f, p, T, branch];
    If[special,
      AppendTo[unkLogs, {p, 0}];
      If[verbose, Print["  (", p, "): special (multiplicity ", mult, "); candidate log(", p, ")"]];
      Continue[]];
    If[verbose, Print["  (", p, "): ", If[branch, "branch", "unramified"], ", delta = ", delta,
        ", v_P(f) = ", vP, If[vP > -delta, "  [sub-critical]", ""]]];
    tauHi = None;            (* canonical residue of a deeper pole over a constant-coefficient prime *)
    If[vP < -delta,
      denv *= p^Ceiling[(-vP - delta)/eP];
      If[q === None && ! branch,
        (* transcendental tower (pmint's setting): the canonical residue at any
           pole order, computed in kappa(P) with no root substitution, so the
           place may have non-constant coordinates and delta may exceed 1
           (Proposition 7.7, general case) *)
        gDir = ResidueDir[T, p];
        tau = If[gDir =!= None, CanonicalResidueField[T, f, p, gDir, delta, verbose], None];
        If[tau =!= None && ! IsZero[tau],
          If[verbose, Print["      canonical residue at order ", -vP, ": ", tau]];
          If[ConstDir[T, p] =!= None,
            tauHi = tau,       (* kappa(P) is algebraic over F: the residues at the places over p, below *)
            cv = CertifyNonconstant[tau, gens];
            Which[cv === True, Throw[{"not elementary", p, tau}, "PIM"],
                  cv === Undecided, Throw[{"failed", "residue constancy undecided", p, tau}, "PIM"],
                  True, AppendTo[detLogs, {tau, TScalar[T, p]}]]]],   (* log(p), coefficient tau *)
      (* deep residues (Proposition 7.6, constant-coefficient case: n = 1, D = d/dx) *)
      If[delta == 1 && ! branch && q =!= None && Length[gens] == 1 && T["derivs"][[1]] === TUnit[T, 0],
        pts = PointsOver[T, p];
        If[pts =!= None,
          taus = DeepResidues[T, f, p, pts, Y];
          If[verbose, Print["      deep residues at order ", -vP, ": ", taus]];
          If[MemberQ[taus, Undecided], Throw[{"failed", "deep residue sheet undecided", p}, "PIM"]];
          cert = CertifyNonconstant[#, gens] & /@ taus;
          If[MemberQ[cert, True], Throw[{"not elementary", p, First[Pick[taus, cert, True]]}, "PIM"]];
          If[MemberQ[cert, Undecided], Throw[{"failed", "residue constancy undecided", p, First[Pick[taus, cert, Undecided]]}, "PIM"]];
          If[AnyTrue[taus, # =!= 0 &],
            got = If[m == 2, RealisePoints[T, p, pts, taus, Y, verbose], RealiseAtPoints[T, p, pts, taus, verbose]];
            If[got === None, AppendTo[torsion, {p, pts, taus}], detLogs = Join[detLogs, got]]]]]]];
    If[vP == -delta || tauHi =!= None,
      (* tau_P = e (f h / Dh)|_P with h = p, or the canonical residue of the deeper pole *)
      If[tauHi =!= None, tp = {tauHi, 0}; texpr = tauHi,
        Dp = TowerD[T, TScalar[T, p]];
        tp = If[q =!= None, TDiv[T, Pscale[eP p, f], Dp], {Can[eP p f[[1]]/Dp[[1]]], 0}];
        texpr = ToY[T, tp, Y]];
      If[branch, Throw[{"failed", "critical branch residue", p}, "PIM"]];
      pts = PointsOver[T, p];
      (* a prime of degree > 4, and every constant-coefficient prime when
         m >= 3: residues in the residue field K[g]/(p) (no roots), the
         residue classes as (prime, sheet polynomial), each realised by
         Hensel lifting and linear algebra over K(c) *)
      If[pts === None || (m >= 3 && ConstDir[T, p] =!= None),
        rc = ResidueClasses[T, p, tp, Y, verbose];
        If[rc === None, Throw[{"failed", "residues over a prime with no constant-coefficient direction", p}, "PIM"]];
        If[StringQ[rc[[1]]] && ! (pts =!= None && rc[[1]] === "failed" && StringContainsQ[rc[[2]], "several sheets"]),
          Throw[rc, "PIM"]];
        If[! StringQ[rc[[1]]],
          {gDir, classes, principal} = rc;
          Do[AppendTo[detLogs, {pr[[1]], TScalar[T, pr[[2]]]}], {pr, principal}];
          Do[got = RealiseClassAny[T, gDir, {{cl[[2]], cl[[3]]}}, cl[[1]], 12, verbose];
            If[got === None, AppendTo[pendingClasses, {cl[[1]], gDir, cl[[2]], cl[[3]]}], AppendTo[detLogs, got]],
            {cl, classes}];
          Continue[]]];
      taus = ResFmt[ReduceAt[texpr, #, Y]] & /@ pts;
      If[verbose,
        If[q === None, Print["      residues ", taus],
          Print["      places ", ("y=" <> ToString[#[[3]]]) & /@ pts, ": residues ", taus]]];
      cert = CertifyNonconstant[#, gens] & /@ taus;
      If[MemberQ[cert, True],
        Throw[{"not elementary", p, First[Pick[taus, cert, True]]}, "PIM"]];
      If[MemberQ[cert, Undecided],
        Throw[{"failed", "residue constancy undecided", p, First[Pick[taus, cert, Undecided]]}, "PIM"]];
      If[AllTrue[taus, Simplify[# - taus[[1]]] === 0 &],
        If[taus[[1]] =!= 0, AppendTo[detLogs, {taus[[1]], TScalar[T, p]}]];
        Continue[]];
      If[q =!= None && m >= 3,
        got = RealiseAtPoints[T, p, pts, taus, verbose];
        If[got === None, AppendTo[torsion, {p, pts, taus}], detLogs = Join[detLogs, got]];
        Continue[]];
      (* quadratic y-split: disc = s^2 q  =>  factors 2 a g + b -+ s y, N = 4 a p *)
      done = False;
      If[q =!= None,
        Do[gstar = gg;
          If[! (PolynomialQ[p, gstar] && Exponent[p, gstar] == 2), Continue[]];
          {c2, b2, a2} = PadRight[CoefficientList[p, gstar], 3];
          If[! FreeQ[{a2, b2, c2}, gstar] || ! (Head[a2] === Integer || Head[a2] === Rational), Continue[]];
          disc = Can[b2^2 - 4 a2 c2];
          s2 = Can[disc/q];
          s = PowerExpand[Sqrt[Factor[s2]]];
          (* s must lie in K: no algebraic function of the generators (algebraic
             constants such as Sqrt[3] are allowed) *)
          If[Simplify[s^2 - s2] =!= 0 || ! FreeQ[s, Power[b_, e_Rational] /; ! FreeQ[b, Alternatives @@ gens]], Continue[]];
          ok = True; pend = {};
          Do[uu = {Expand[2 a2 gstar + b2], -sg s};
            tv = Pick[taus, Simplify[ReduceAt[uu[[1]] + uu[[2]] Y, #, Y]] === 0 & /@ pts, True];
            If[tv =!= {} && AllTrue[tv, # === tv[[1]] &], AppendTo[pend, {tv[[1]], uu}], ok = False],
            {sg, {1, -1}}];
          If[ok,
            Do[If[tv0[[1]] =!= 0, AppendTo[detLogs, tv0]];
              If[verbose, Print["      y-split factor ", tv0[[2, 1]], " + (", tv0[[2, 2]], ")*y: residue ", tv0[[1]]]],
              {tv0, pend}];
            done = True; Break[]],
          {gg, gens}]];
      If[done, Continue[]];
      got = RealisePoints[T, p, pts, taus, Y, verbose];
      If[got === None, AppendTo[torsion, {p, pts, taus}], detLogs = Join[detLogs, got]]],
    {fac, fl}];
  If[pendingClasses =!= {},
    (* a residue class may be supported over several primes: realise the
       classes with equal value jointly *)
    groups = {};
    Do[found = False;
      Do[If[groups[[i, 2]] === pcl[[2]] && IsZero[groups[[i, 1]] - pcl[[1]]],
          groups[[i, 3]] = Append[groups[[i, 3]], {pcl[[3]], pcl[[4]]}]; found = True; Break[]],
        {i, Length[groups]}];
      If[! found, AppendTo[groups, {pcl[[1]], pcl[[2]], {{pcl[[3]], pcl[[4]]}}}]],
      {pcl, pendingClasses}];
    unrealised = {};
    Do[got = If[Length[grp[[3]]] > 1, RealiseClassAny[T, grp[[2]], grp[[3]], grp[[1]], 12, verbose], None];
      If[got === None, AppendTo[unrealised, {Times @@ (#[[1]] & /@ grp[[3]]), {grp[[1]]}}], AppendTo[detLogs, got]],
      {grp, groups}];
    If[unrealised =!= {}, Throw[{"needs torsion realisation (milestone iii)", unrealised}, "PIM"]]];
  (* the divisor may be realisable only jointly over several primes: by
     torsion on a cubic model (Algorithm 3(d)) *)
  If[torsion =!= {} && m == 2,
    got = TorsionRealise[T, torsion, Y, 24, verbose];
    If[got =!= None, detLogs = Join[detLogs, got]; torsion = {},
      (* ... or not realisable at all: a component of the divisor that is
         provably non-torsion (Proposition 9.4 by reduction mod p) is a
         certificate of non-elementarity (Corollary 7.6) *)
      certd = NontorsionDivisor[T, torsion, verbose];
      If[certd =!= None, Throw[{"not elementary", "residue divisor not torsion: reduction mod p", certd[[1]], certd[[2]]}, "PIM"]]]];
  If[torsion =!= {}, Throw[{"needs torsion realisation (milestone iii)", {#[[1]], #[[3]]} & /@ torsion}, "PIM"]];

  (* residue at the hypertangent place at infinity (Lemma 8.1) *)
  If[q =!= None && Length[gens] >= 2,
    tinf = VinftyResidue[T, f, Y];
    If[tinf =!= None,
      If[verbose, Print["  v_oo (hypertangent top): delta = 1, v_oo(f) = -1, residue ", ToY[T, tinf, Global`y]]];
      lower = Most[gens];
      nonconst = AnyTrue[Rest[tinf], Can[#] =!= 0 &] || CertifyNonconstant[tinf[[1]], lower] === True;
      If[nonconst, Throw[{"not elementary", "v_oo", tinf}, "PIM"]]]];

  (* tower specials (Theorem 6.1): offered regardless of the integrand *)
  seen = unkLogs[[All, 1]];
  cand = {};
  Do[Do[If[comp === 0, Continue[]];
      Do[If[FreeQ[part, Alternatives @@ gens], Continue[]];
        Do[If[! FreeQ[fac2[[1]], Alternatives @@ gens], AppendTo[cand, fac2[[1]]]],
          {fac2, Rest[FactorList[part]]}],
        {part, {Numerator[Can[comp]], Denominator[Can[comp]]}}],
      {comp, dpair}],
    {dpair, Join[T["derivs"], If[q =!= None, Rest[T["Dw"]], {}]]}];
  Do[If[MemberQ[seen, pp], Continue[]];
    If[ClassifyPrime[T, pp][[4]],
      AppendTo[unkLogs, {pp, 0}]; AppendTo[seen, pp];
      If[verbose, Print["  tower special: candidate log(", pp, ")"]]],
    {pp, DeleteDuplicates[cand]}];

  (* specials over Fbar (Theorem 6.1): on request replace each special p(g)
     with constant coefficients and degree >= 2 by its linear factors g - r *)
  splittable = AnyTrue[unkLogs, Function[pl, AnyTrue[gens,
      PolynomialQ[pl[[1]], #] && Exponent[pl[[1]], #] >= 2 && FreeQ[CoefficientList[pl[[1]], #], Alternatives @@ gens] &]]];
  (* residue-invisible unit candidates (Remark 7.7): the fundamental unit of
     an even-degree radicand from the continued fraction of Sqrt[q]; when the
     search is inconclusive, Proposition 9.4 may certify that no unit exists *)
  If[q =!= None,
    Do[If[! FreeQ[q, Alternatives @@ DeleteCases[gens, g]], Continue[]];
      If[m == 2 && EvenQ[Exponent[q, g]],
        uu0 = FundamentalUnit[q, g];
        If[uu0 =!= None,
          AppendTo[unitsBase, {uu0[[1 ;; 2]], uu0[[3]]}];
          If[verbose, Print["  unit candidate: A + B*y with deg_", g, " B = ", Exponent[uu0[[2]], g]]],
          certd = NontorsionCertificate[q, g];
          If[certd[[1]],
            If[verbose, Print["  unit search inconclusive; [oo+ - oo-] certified non-torsion by reduction mod p: ", certd[[2]]]];
            unitsComplete = True,
            unitsComplete = False]]];
      If[m >= 3 && FreeQ[CoefficientList[q, g], Alternatives @@ gens],
        (* the bounded divisor search at the places at infinity *)
        unitsBase = Join[unitsBase, {#, None} & /@ UnitsGeneral[T, g, None, verbose]];
        unitsComplete = GCD[m, Exponent[q, g]] == 1];
      Break[],
      {g, gens}]];

  (* residual integrand *)
  rem = f;
  Do[ld = If[q =!= None, TDiv[T, TowerD[T, u[[2]]], u[[2]]], {Can[TowerD[T, u[[2]]][[1]]/u[[2, 1]]], 0}];
    rem = Padd[rem, Pscale[-u[[1]], ld]],
    {u, detLogs}];

    A = <|"f" -> f, "Y" -> Y, "detLogs" -> detLogs, "unkLogs" -> unkLogs, "denv" -> denv, "unitsBase" -> unitsBase,
          "unitsComplete" -> unitsComplete, "rem" -> rem, "splittable" -> splittable, "specials" -> <||>, "spec" -> <||>, "parent" -> <||>|>;
    $analyses[key] = A];
  {f, Y, detLogs, unkLogs, denv, unitsBase, unitsComplete, rem, splittable} =
    Lookup[A, {"f", "Y", "detLogs", "unkLogs", "denv", "unitsBase", "unitsComplete", "rem", "splittable"}];

  (* Steps 5--6 for the specials over Q or over Fbar: the special logands and
     the S'-units over them, once per variant *)
  spec = Lookup[A["specials"], split === True, None];
  If[spec === None,
    sunitsAll = {};
  If[split === True && splittable,
    newLogs = {};
    Do[PP = pl[[1]];
      gstar = SelectFirst[gens, PolynomialQ[PP, #] && Exponent[PP, #] >= 2 && FreeQ[CoefficientList[PP, #], Alternatives @@ gens] &, None];
      If[gstar === None, AppendTo[newLogs, pl],
        rts = RR /@ (gstar /. Solve[PP == 0, gstar]);          (* canonical algebraic numbers *)
        Do[AppendTo[newLogs, {gstar - r, 0}]; $analyses[key, "parent", gstar - r] = PP, {r, rts}]],
      {pl, unkLogs}];
    unkLogs = newLogs;
    If[verbose, Print["  specials split over Fbar: ", unkLogs[[All, 1]]]]];

  (* S'-units split over a special prime (Algorithm 4, Step 6: generators of
     O_{S'}^* ): u = a + b y with N(u) = c p^k for a special p, from the norm
     search -- e.g. 1 +- Sqrt[1 - x^2] over the special x of Log[x], needed
     for Log[x] ArcSin[x].  Generators only: constant b first, rational
     coefficients before algebraic ones, at most four per special. *)
     
  If[q =!= None && m == 2,
    Do[gstar = SelectFirst[gens, FreeQ[q, Alternatives @@ DeleteCases[gens, #]] && PolynomialQ[pl[[1]], #] &&
                                  Exponent[pl[[1]], #] >= 1 && FreeQ[CoefficientList[pl[[1]], #], Alternatives @@ gens] &, None];
      If[gstar === None, Continue[]];
      sols = Select[NormSearchAll[q, pl[[1]], gstar, 2, 1], #[[2]] =!= 0 &];
      (* Order the norm solutions so the four-generator cap covers conjugate
         orbits, not just one orbit's several generators.  NormSearchAll returns
         solutions of every norm class c (a^2 - q b^2 = c p^k); Galois-conjugate
         classes -- e.g. A40's c = 1/5 -+ 2I/5 over Q(Sqrt[5], I) -- are DISTINCT
         classes and BOTH are needed to span the log-part ansatz (Mathematica
         collects one generator, both y-signs, from EACH orbit).  The old SortBy
         keys were orbit-blind, so the cap filled all four slots from the first
         orbit's two generators and permanently missed the conjugate orbit ->
         A40 "no solution".  Fix: (1) a TOTAL canonical order (the trailing Re/Im
         keys break every tie deterministically, killing the Solve-branch-order
         non-determinism that made A19 wrong off-gate between runs); (2) a
         round-robin by norm class so the first generators taken are one per
         orbit.  A single-orbit case is left in its original within-orbit order
         (round-robin over one group is the identity), so this neither dedups
         within an orbit (A19 needs two independent same-class generators) nor
         raises the cap (over-completion let the exact solve pick a spurious
         particular solution). *)
      sols = SortBy[sols, {Exponent[#[[2]], gstar],
                           If[AllTrue[Join[CoefficientList[#[[1]], gstar], CoefficientList[#[[2]], gstar]], MatchQ[#, _Integer | _Rational] &], 0, 1],
                           #[[4]],
                           Re[N[#[[3]], 30]], Im[N[#[[3]], 30]],
                           Re[N[#[[1]] /. gstar -> 3/2, 30]], Im[N[#[[1]] /. gstar -> 3/2, 30]]} &];
      sols = With[{grps = GatherBy[sols, RR[#[[3]]] &]},
               Flatten[Table[Select[grps, Length[#] >= j &][[All, j]], {j, Max[Length /@ grps]}], 1]];
      sunits = {};
      Do[If[Length[sunits] >= 4, Break[]];
        Do[uuS = {Expand[sol[[1]]], Expand[sg sol[[2]]]};
          If[! MemberQ[sunits, uuS] && ! MemberQ[sunits, -uuS],
            AppendTo[sunits, uuS];
            If[verbose, Print["  S'-unit over special (", pl[[1]], "): ", uuS[[1]], " + (", uuS[[2]], ")*y  (N = ", sol[[3]], "*(", pl[[1]], ")^", sol[[4]], ")"]]],
          {sg, {1, -1}}],
        {sol, sols}];
      (* torsion classes of the places over the special (cubic models): the
         Miller functions are S'-units too *)
      If[Exponent[q, gstar] == 3,
        With[{ptsS = PointsOver[T, pl[[1]]]},
          If[verbose, Print["  torsion S'-units over (", pl[[1]], "): places ", ptsS]];
          If[ptsS =!= None && AllTrue[ptsS, FreeQ[{#[[2]], #[[3]]}, Alternatives @@ gens] &],
            Do[If[Length[sunits] >= 8, Break[]];
              With[{resS = TorsionOrderAndMiller[q, gstar, ptS, Y]},
                If[verbose, Print["      TorsionOrderAndMiller at ", ptS, " -> ", If[resS === None, None, resS[[1]]]]];
                If[resS =!= None && ! MemberQ[sunits, resS[[2]]],
                  AppendTo[sunits, resS[[2]]];
                  If[verbose, Print["  S'-unit over special (", pl[[1]], "): Miller function of [P - oo] of order ", resS[[1]], " at (", ptS[[2]], ", ", ptS[[3]], ")"]]]],
              {ptS, ptsS}]]]];
      sunitsAll = Join[sunitsAll, {#, None} & /@ sunits],
      {pl, unkLogs}]];
  If[q =!= None && m >= 3 && unkLogs =!= {} && verbose, Print["  S'-units over the specials are not searched for m = ", m]];

    $analyses[key, "specials", split === True] = {unkLogs, sunitsAll},
    {unkLogs, sunitsAll} = spec];
  units = Join[sunitsAll, unitsBase];

  (* Step 15 by Algorithm 6: the special exponents and the degree bounds from
     the shift of the derivation at every place, proved where Proposition
     8.11 applies, the classical guess raised by the retry count elsewhere *)
  If[bounds === None,
    {bounds, exps, proved} = PlaceBounds[key, T, rem, denv, unkLogs, sexp, verbose],
    proved = False;
    exps = Table[Max[0, -VP[rem, pl[[1]], T, q =!= None && Vp[q, pl[[1]], gens] > 0]] + sexp, {pl, unkLogs}]];
  Do[denv *= unkLogs[[i, 1]]^exps[[i]], {i, Length[unkLogs]}];
  monos = Tuples[Range[0, #] & /@ bounds];
  (* one family of coefficients per coordinate: c0[...], c1[...] for m <= 2 *)
  css = Table[cc[i - 1] @@@ monos, {i, nc}];
  cs0 = css[[1]]; cs1 = If[q =!= None, css[[2]], {}];
  V = PadRight[Table[Total[MapThread[#1 Times @@ (gens^#2) &, {css[[i]], monos}]]/denv, {i, nc}], T["n"]];
  gammas = Table[gamma[i], {i, Length[units]}];
  betas = Table[beta[i], {i, Length[unkLogs]}];
  unks = Join[Flatten[css], gammas, betas];
  (* the system is linear in the unknowns: every basis element of the ansatz
     (a monomial over D_v, a unit or S'-unit logand, a special logand)
     contributes an unknown-free column, and the assembly and the solve run
     over ONE number field (AnsatzSystem, the counterpart of _assemble_field
     in parallel_mixed.py).  The Expr route below -- Together on the
     unknown-bearing pairs, then Solve -- is taken only when a constant of the
     system is not an algebraic number. *)
  sysK = AnsatzSystem[T, rem, denv, units, unkLogs, css, monos, gammas, betas, unks];
  If[sysK === $Failed,
    (* assembled WITHOUT Cancel: with the unknown coefficients present a gcd
       would treat them as variables (catastrophic over an algebraic
       extension), and a spurious common factor only rescales an equation *)
    EE = Together /@ (TowerD[T, V] - rem);
    Do[EE = Together /@ (EE + gammas[[i]] TDiv[T, TowerD[T, units[[i, 1]]], units[[i, 1]]]), {i, Length[units]}];
    Do[EE = Together /@ (EE + betas[[i]] *
        If[q =!= None, TDiv[T, TowerD[T, TScalar[T, unkLogs[[i, 1]]]], TScalar[T, unkLogs[[i, 1]]]],
                       {Can[TowerD[T, {unkLogs[[i, 1]], 0}][[1]]/unkLogs[[i, 1]]], 0}]),
      {i, Length[unkLogs]}];
    eqs = Flatten[If[# === 0, {}, CoefficientRules[Numerator[Together[#]], gens][[All, 2]]] & /@ EE];
    neq = Length[eqs];
    sol = If[eqs === {}, {{}}, Quiet[Solve[Thread[eqs == 0], unks], Solve::svars]],
    neq = sysK["neq"];
    sol = If[sysK["sub"] === None, {}, {sysK["sub"]}]];
  If[verbose, Print["  ansatz: bounds ", bounds, ", ", Length[unks], " unknowns (", Length[Flatten[css]], " coefficients, ",
    Length[gammas], " unit/S'-unit logs, ", Length[betas], " special logs), ", neq, " equations; ",
    If[sol === {}, "no solution", "solved"]]];
  If[sol === {},
    (* Proposition 9.2(b): every bound in force is proved, the logand
       candidates are complete -- the unit group known, the specials split
       over Fbar, none of them over the curve variable, where the S'-units are
       only searched within a bound -- and the residual VERIFIED residue-free:
       the inconsistent system shows the residual to be a non-exact
       differential of the second kind (the holomorphic remainder) *)
    curve = Select[gens, ! FreeQ[q, #] &];
    typeE = q =!= None && AnyTrue[unkLogs, FreeQ[#[[1]], Alternatives @@ Complement[gens, curve]] &];
    (* the holomorphic-remainder (second-kind) certificate presupposes an
       algebraic curve of positive genus; on a no-curve tower (q === None) no
       such differential exists, so never certify there -- fall through to
       {"failed", ...} instead (this is the T2/T10 false-certificate guard). *)
    If[q =!= None && OptionValue["Bounds"] === None && proved && unitsComplete && ! typeE && ! (splittable && split =!= True) &&
        If[Length[gens] == 1 && q =!= None && m == 2 && T["derivs"][[1]] === {1, 0}, ResidueFree[T, rem, Y], VerifiedResidueFree[rem, T]],
      Throw[{"not elementary", "holomorphic remainder: residual second-kind differential is not exact (every bound in force is proved)", bounds}, "PIM"]];
    Throw[Join[{"failed", "no solution within bounds", bounds},
      If[splittable && split =!= True, {"splittable"}, {}], If[! proved, {"guess"}, {}]], "PIM"]];
  sub = First[sol];
  If[sysK === $Failed,
    frees = Complement[unks, sub[[All, 1]]];
    sub = Join[sub, Thread[frees -> 0]] /. (a_ -> b_) :> (a -> (b /. Thread[frees -> 0]));
    If[AnyTrue[EE, Can[# /. sub] =!= 0 &], Throw[{"failed", "internal: solution does not verify", bounds}, "PIM"]]];

  y = If[q =!= None, q^(1/m), None];
  surf[u_] := If[q === None, u[[1]], Sum[TPad[T, u][[i + 1]] y^i/T["E"][[i + 1]], {i, 0, T["n"] - 1}]];
  I0 = surf[Can /@ (V /. sub)]
       + Total[#[[1]] Log[surf[#[[2]]]] & /@ detLogs]
       + Total[Table[(gammas[[i]] /. sub) Log[surf[units[[i, 1]]]], {i, Length[units]}]]
       + Total[Table[(betas[[i]] /. sub) Log[unkLogs[[i, 1]]], {i, Length[unkLogs]}]];
  I0];

End[];
EndPackage[];

(* The Integrate`ParallelMixedTower method symbol is a C builtin (see
   src/calculus/integrate.c: builtin_integrate_pmt) that lazily loads this file
   on first use and delegates to ParallelMixed`ParallelIntegrateMixed, so the
   ~1s package load is paid only by sessions that actually reach the method. *)
