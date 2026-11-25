(* gen_compact.m -- benchmark 257: compact integrands with large integrals.

   Run:  SEED=1 CLASS=A QUOTA=300 ./Mathilda -file tools/integration_shootout/gen_compact.m > part.tsv
   (one process per class; benchmarks/257-compact-integrals/README.md has the
   seeds and quotas that built corpus.tsv).

   A SEARCH over compact integrands: each class draws f from a small grammar
   (shapes with the monomial's derivative as a factor, theta' R(theta, y), so
   that integrable ones are common; radicals with small exponents; irreducible
   denominators), puts it in compact form (compactForm), and asks Mathilda's
   Automatic Integrate for an antiderivative F.  A candidate is
     CASE  (accepted, expect = elem, F = Mathilda's verified antiderivative)
           when PossibleZeroQ[D[F] - f] and LeafCount[compactForm[F]] >=
           2 LeafCount[f];
     CAND  (handed to gen_compact.py, which asks FriCAS) when Mathilda neither
           closes nor certifies it -- so the corpus is not limited to what
           Mathilda can already do;
     dropped when Mathilda certifies it non-elementary, or closes it with a
           small integral.
   (Constructing F and differentiating, the first design, was tried: D[F]
   almost never collapses to half of F for algebraic and mixed F -- 1 case in
   600 attempts for class A.)

   Classes (the shootout's method map): R rational, T transcendental, A algebraic
   (Trager), M radical on top of a Log/Exp tower (Bronstein Shape A), B radical at
   the bottom (Bronstein Shape B, component and curve-monomial). *)
seed = ToExpression[Environment["SEED"]]; If[!IntegerQ[seed], seed = 1];
cls = Environment["CLASS"]; If[!StringQ[cls], cls = "A"];
quota = ToExpression[Environment["QUOTA"]]; If[!IntegerQ[quota], quota = 20];
gname = Environment["GRAMMAR"]; If[!StringQ[gname], gname = cls];
(* the leaf-count ratio an accepted F must reach over f: RATIO=2 (the default,
   benchmark 257) accepts LeafCount[Fc] >= 2 LeafCount[f]; any other RATIO is
   strict, LeafCount[Fc] > RATIO LeafCount[f] (benchmark 261: 4 LeafCount[D[F]] <
   LeafCount[F]).  ATTEMPTS is the attempt cap per accepted case (40). *)
ratio = ToExpression[Environment["RATIO"]]; If[!NumericQ[ratio], ratio = 2];
attFactor = ToExpression[Environment["ATTEMPTS"]]; If[!IntegerQ[attFactor], attFactor = 40];
ratioOK[lF_, lf_] := If[ratio === 2, lF >= 2 lf, lF > ratio lf];
fair = ratio =!= 2 || Environment["FAIR"] === "1";
SeedRandom[seed + 1000 First[ToCharacterCode[cls]]];
here = DirectoryName[$InputFileName];
Get[FileNameJoin[{here, "families.m"}]];
Get[FileNameJoin[{here, "common.m"}]];

(* ---------------- compact forms ---------------- *)
radicalOf[e_] := Module[{ps = DeleteDuplicates[Cases[{e}, Power[b_, q_Rational] /; !FreeQ[b, x] :> {Together[b], Denominator[q]}, {0, Infinity}]]},
  If[ps === {} || Length[DeleteDuplicates[First /@ ps]] > 1, None, {ps[[1, 1]], LCM @@ (Last /@ ps)}]];

(* e in canonical form over K(Y), Y^n = b: the denominator rationalised by the
   extended gcd with Y^n - b, every power of Y reduced below n *)
canonRad[e_] := Module[{rd = radicalOf[e], b, n, Y, ey, t, num, den, eg, inv, red, cl},
  If[rd === None, Return[e]];
  {b, n} = rd; Y = Unique["Y"];
  ey = e /. Power[c_, q_Rational] /; !FreeQ[c, x] && Together[c - b] === 0 :> Y^(q n);
  If[!FreeQ[ey, Power[c_, _Rational] /; !FreeQ[c, x]], Return[e]];
  t = Together[ey]; num = Numerator[t]; den = Denominator[t];
  red[p_] := PolynomialRemainder[Expand[p], Y^n - b, Y];
  If[!FreeQ[den, Y],
    eg = PolynomialExtendedGCD[red[den], Y^n - b, Y];
    If[!FreeQ[eg[[1]], Y], Return[e]];
    num = red[num eg[[2, 1]]]; den = eg[[1]]];
  cl = CoefficientList[red[num], Y];
  (* Y inside a non-polynomial head (Log[x + Y]) is not reduced: keep e *)
  If[!FreeQ[cl, Y] || !FreeQ[den, Y], Return[e]];
  Sum[Together[cl[[i]]/den] (b^(1/n))^(i - 1), {i, Length[cl]}]];

compactForm[e_] := Module[{cs},
  cs = {Together[e], Quiet[canonRad[e]]};
  cs = Join[cs, {Quiet[TimeConstrained[Factor[cs[[2]]], 2, $Aborted]]}];
  If[LeafCount[First[MinimalBy[cs, LeafCount]]] < 160,
    AppendTo[cs, Quiet[TimeConstrained[Simplify[First[MinimalBy[cs, LeafCount]]], 2, $Aborted]]]];
  cs = Select[cs, FreeQ[#, $Aborted | $Failed | Indeterminate | ComplexInfinity] &];
  First[MinimalBy[cs, LeafCount]]];

(* ---------------- grammar helpers ---------------- *)
sm[] := pick[{1, 2, 3, -1, -2, 1/2, 3/2}];
lpoly[d_] := Sum[r[-2, 2] x^k, {k, 0, d - 1}] + nz[1, 2] x^d;
upoly[u_, d_] := Sum[r[-2, 2] u^k, {k, 0, d - 1}] + nz[1, 2] u^d;
expo[] := pick[{1/2, -1/2, 1/3, -1/3, 2/3, 1/4, -1/4, 3/2, -3/2, 1/6}];
thQ[] := pick[{Log[x], Log[x + 1], Log[2 x + 1], E^x, E^(2 x), E^(-x), Log[x^2 + 1], E^(x^2), Log[1 + E^x]}];
binom[u_] := pick[{u^2 + nz[1, 3], u^3 + nz[-2, 2], u^4 + nz[1, 3], u^2 + u + nz[1, 2], u^3 - u + nz[1, 2], u + nz[1, 3], u^6 + 1}];

(* ---------------- the grammar per class ---------------- *)
grammar["R"] := Function[{}, Together[pick[{1, x, x^2, lpoly[1], lpoly[2]}]/
  pick[{lpoly[3], lpoly[4], lpoly[3] lpoly[1], (x^2 + nz[1, 3])^2 lpoly[1], lpoly[5], x^6 + nz[1, 3]}]]];

grammar["T"] := Function[{}, Module[{t = thQ[], dt},
  dt = Together[D[t, x]];
  pick[{dt upoly[t, pick[{0, 1}]]/upoly[t, pick[{2, 3, 4}]], dt Log[upoly[t, 2]], dt/(t upoly[t, 2]),
        t^pick[{1, 2}]/(x^2 + nz[1, 3]) + dt/(t^2 + nz[1, 3]), dt t/(t^4 + nz[1, 3]), x dt/(t^2 + nz[1, 2])}]]];

grammar["A"] := Function[{}, Module[{g = binom[x], e = expo[]},
  pick[{x^r[0, 3] g^e, x^r[0, 2] g^e/(x + nz[-2, 2]), g^e/(x^2 + nz[1, 3]), (x + nz[-2, 2]) g^e/lpoly[2],
        x^r[1, 3]/(lpoly[pick[{2, 3}]] g^Abs[e]), g^e/x, (x^2 + nz[-2, 2]) g^e/(x^2 - nz[1, 3])}]]];

grammar["M"] := Function[{}, Module[{t = thQ[], dt, g, e = expo[]},
  dt = Together[D[t, x]]; g = binom[t];
  pick[{dt t^r[0, 2] g^e, dt g^e/(t + nz[-2, 2]), dt g^e/(t^2 + nz[1, 3]), dt/(t g^Abs[e]),
        dt (t + nz[-2, 2]) g^e/upoly[t, 2], (x + t) dt g^e, dt g^e + x g^e/(x + 1)}]]];

grammar["B"] := Function[{}, Module[{t = thQ[], dt, g = binom[x], e = expo[]},
  dt = Together[D[t, x]];
  pick[{t g^e, t g^e/x, dt g^e, t^2 g^e/(x + nz[1, 3]), Log[x + g^(1/2)] g^e, E^(x g^(1/2)) g^e,
        dt g^e/(t + nz[1, 3]), Log[x + g^(1/2)]/x, (t + x) g^e/(x^2 + nz[1, 3])}]]];

(* B2: the radical-at-the-bottom grammar used for the later class-B seeds.  B's
   own grammar turned out to draw mostly kernel-free integrands (dt g^e with
   dt rational) and many non-elementary ones; B2 keeps to genus-0 radicals of
   x (linear or quadratic radicands) times a Log/Exp kernel, the shapes the
   component method and the genus-0 substitution decide, and its candidates
   must contain both a radical and a kernel (contentOK). *)
lin2[] := nz[1, 3] x + nz[-3, 3];
quad2[] := pick[{x^2 + nz[1, 4], x^2 - nz[1, 4], x^2 + x + nz[1, 3], x^2 + nz[-2, 2] x}];
grammar["B2"] := Function[{}, Module[{g = pick[{lin2[], lin2[], quad2[]}], e = expo[], L = Log[lin2[]], q = nz[1, 4]},
  pick[{L g^e, L g^e/x, x L g^e, L g^e/(x + nz[1, 3]), Log[x] g^e/(x + nz[-2, 2]),
        Log[x + Sqrt[x^2 + q]], x Log[x + Sqrt[x^2 + q]], Log[x + Sqrt[x^2 + q]]/x^2, Log[x + Sqrt[x + q]],
        E^Sqrt[lin2[]], x E^Sqrt[x + q], E^Sqrt[x + q] (x + q)^pick[{1/2, -1/2, 3/2}],
        Log[lin2[]]^2 g^e, (x + nz[-2, 2]) Log[x] g^e/(x^2 + q)}]]];
contentOK[f_] := If[MemberQ[{"B2", "B3", "B4", "M3", "M4"}, gname], !FreeQ[f, Log | Power[E, _]] && !FreeQ[f, Power[_, _Rational]], True];

(* EXCLUDE_EASY=1 (benchmark 262): an integrand in an easy class (easy.m) is
   dropped -- rational, linear-fractional radical, Euler-rationalisable, the
   "-sub" forms of the last two, and (after Integrate, so only for integrands
   that would otherwise be kept) any that Method -> "DerivativeDivides" closes. *)
excludeEasy = Environment["EXCLUDE_EASY"] === "1";
If[excludeEasy, Get[FileNameJoin[{DirectoryName[$InputFileName], "easy.m"}]]];
(* counts per class; a list, since Increment on an Association key is not
   supported (tasks/todo.md, ASSOC-INC-1) *)
easyNames = {"rational", "degenerate", "logradical", "linrad", "euler", "nested", "dd"}; nEasy = {0, 0, 0, 0, 0, 0, 0};
easyCount[c_] := With[{i = First[FirstPosition[easyNames, StringReplace[c, "-sub" -> ""]]]}, nEasy[[i]] = nEasy[[i]] + 1];
syntacticEasyQ[f_] := excludeEasy && Module[{c = easySyntacticClass[f]},
  If[c =!= "none", easyCount[c]; True, False]];
ddEasyQ[f_] := excludeEasy && If[ddClosesQ[f], easyCount["dd"]; True, False];

(* ---------------- benchmark 262: grammars away from the easy classes ----------------
   T3, A3, M3 and B3 draw integrands that are, by construction, mostly NOT
   rational, NOT in Q(x, ((a x + b)/(c x + d))^(1/m)), NOT in Q(x, Sqrt[quadratic])
   and NOT of the derivative-divides shape theta' R(theta): a rational cofactor in
   x that is not the kernel's derivative, two kernels, radicals of cubics and
   quartics, cube and fourth roots of quadratics, nested and paired radicals.  The
   grammar only biases the draw; EXCLUDE_EASY=1 (easy.m) is what keeps the easy
   ones out. *)
(* T3: D[Q] Log[S], Q rational with an irreducible or repeated denominator: by
   parts, Q Log[S] - Integrate[Q S'/S], always elementary, and its rational part
   brings logs and arctangents with algebraic constants.  Plus a few two-kernel
   and Log-of-rational shapes. *)
qrat[] := pick[{1/(x + nz[-3, 3]), x/(x^2 + nz[1, 3]), 1/(x^2 + nz[1, 3]), (x + nz[-2, 2])/(x^2 + x + nz[1, 3]),
  1/(x + nz[-3, 3])^2, x^2/(x^3 + nz[1, 3]), (x + nz[-2, 2])/((x + nz[-3, 3]) (x^2 + nz[1, 2])), 1/(x^3 + nz[-2, 2]), x/(x^2 - nz[1, 3])}];
srat[] := pick[{x + nz[-3, 3], x^2 + nz[1, 3], (x + nz[-3, 3])/(x + nz[-3, 3]), x^2 + x + nz[1, 3], x^3 + nz[-2, 2], x, x^2 - nz[1, 3]}];
tk3[] := pick[{Log[x], Log[x + nz[1, 3]], Log[x^2 + nz[1, 3]], E^x, E^(2 x), E^(-x)}];
grammar["T3"] := Function[{}, Module[{Q = qrat[], S = srat[]},
  pick[{D[Q, x] Log[S], D[Q, x] Log[S], D[Q, x] Log[S], D[Q, x] Log[S] + D[qrat[], x] Log[srat[]],
        Q Log[S] D[srat[], x]/srat[], D[Q, x] Log[S]^2 + 2 Q Log[S] D[S, x]/S,
        E^x (Q + D[Q, x]) + D[qrat[], x] Log[S], x E^x/(x + nz[1, 3])^2, Log[x]/(x + nz[1, 3])^r[2, 3]}]]];

(* A3: algebraic integrands with an elementary integral outside the easy classes.
   Chebyshev binomial differentials x^m (b x^n + a)^p with (m + 1)/n + p an
   integer and (m + 1)/n not (the substitution u = x^n, derivative-divides, would
   decide it); Goursat-type (x^n - c)/((x^n + c) (x^n + b)^(1/n)); pseudo-elliptic
   (x^2 - c)/((x^2 + c) Sqrt[x^4 + b x^2 + c^2]); and a wider search over radicals
   of cubics and quartics and cube and fourth roots of quadratics. *)
chebyshevBinom[] := Module[{n = pick[{3, 3, 4, 4, 6, 3}], p, ms},
  p = pick[Select[{1/3, -1/3, 2/3, -2/3, 1/4, -1/4, 3/4, -3/4, 1/6, -1/6, 5/6, -5/6, 4/3, -4/3, 5/4, -5/4}, Mod[n, Denominator[#]] == 0 &]];
  ms = Select[Range[-7, 7], IntegerQ[(# + 1)/n + p] && !IntegerQ[(# + 1)/n] &];
  If[ms === {}, Return[x^2 (x^3 + 1)^(1/3)]];
  x^pick[ms] (nz[1, 3] x^n + nz[-3, 3])^p];
g3[] := pick[{x^3 + nz[-2, 2], x^3 + x + nz[-2, 2], x^4 + nz[1, 3], x^4 - nz[1, 3], x^4 + nz[-2, 2] x^2 + nz[1, 3],
  x^6 + nz[1, 2], x (x^2 + nz[1, 3]), x^4 + x + nz[1, 2]}];
q3[] := pick[{x^2 + nz[1, 3], x^2 - nz[1, 3], x^2 + x + nz[1, 2], x^3 + nz[1, 2], x^2 + nz[-2, 2] x}];
grammar["A3"] := Function[{}, Module[{g = g3[], q = q3[], c = nz[1, 3], b = nz[-3, 3], n = pick[{3, 3, 4, 6}],
    es = pick[{1/2, -1/2, 3/2, -3/2}], eq = pick[{1/3, -1/3, 2/3, -2/3, 1/4, -1/4, 3/4}]},
  pick[{chebyshevBinom[], chebyshevBinom[], chebyshevBinom[], chebyshevBinom[] lpoly[1], chebyshevBinom[]/(x^n + nz[1, 3]),
        (x^n - c)/((x^n + c) (x^n + b)^(1/n)), x^r[0, 2] (x^n - c)/((x^n + c) (x^n + b)^(1/n)),
        (x^2 - c)/((x^2 + c) Sqrt[x^4 + b x^2 + c^2]), (x^2 + c)/((x^2 - c) Sqrt[x^4 + b x^2 + c^2]),
        lpoly[pick[{0, 1, 2}]] g^es/lpoly[pick[{1, 2}]], g^es/(x^2 + nz[1, 3]),
        q^eq/(x + nz[-2, 2]), x^r[0, 2] q^eq/(x^2 + nz[1, 3]), 1/((x + nz[-2, 2]) q^Abs[eq])}]]];

(* M3: a radical of a Log/Exp polynomial whose cofactor is not the kernel's
   derivative.  Elementary integrals of this shape are rare; the search keeps
   whatever it finds. *)
grammar["M3"] := Function[{}, Module[{t = thQ[], g, e = expo[]},
  g = binom[t];
  pick[{x g^e, g^e/(x + nz[-2, 2]), (x + t) g^e, x^r[1, 2] Together[D[t, x]] g^e, t g^e/(x^2 + nz[1, 3]),
        (x + nz[-2, 2]) g^e/(t + nz[-2, 2]), Together[D[t, x]] g^e + x g^e, x g^e/(t + nz[1, 3])}]]];

(* M4: a radical of an Exp/x polynomial built on a unit of the tower.  With
   g = P^2 - c E^h (P in x and E^x, h one of m x, x^2, ...), P^2 - g = c E^h is a
   unit, so D[Log[P + Sqrt[g]]] = D[h]/2 + A/Sqrt[g] with A = D[P] - D[h] P/2:
   A/Sqrt[g] integrates to Log[P + Sqrt[g]] - h/2 although no kernel divides it
   (x appears outside every kernel).  The grammar draws A/Sqrt[g] and its
   perturbations (a cofactor, a shift of A, a second radical term D[R Sqrt[g]]),
   and Integrate decides which are elementary.  Random M3 draws are almost
   never elementary; these are in about one case in five.  F of this shape is
   rarely twice f (Log[P + Sqrt[g]] holds P twice, as f does), so M runs at
   ratio 1 (gen_compact.py, class_ratio). *)
m4P[] := pick[{x + nz[-3, 3], x^2 + nz[-3, 3], x + nz[-2, 2] E^x, x^2 + nz[-2, 2] x, 2 x + nz[-3, 3],
  x^2 + x + nz[-2, 2], 1/x + nz[-2, 2], x + 1/x, x^3 + nz[-2, 2], x E^x + nz[-2, 2]}];
m4h[] := pick[{x, x, 2 x, -x, -2 x, 3 x, x^2, -x^2, x^2 + x}];
grammar["M4"] := Function[{}, Module[{P = m4P[], h = m4h[], c = nz[-3, 3], g, A, R},
  g = Together[Expand[P^2 - c E^h]];
  A = Together[D[P, x] - D[h, x] P/2];
  R = pick[{x, 1/x, P, 1, E^x}];
  pick[{A/Sqrt[g], x A/Sqrt[g], A/(x Sqrt[g]), (A + nz[-2, 2])/Sqrt[g], (A + nz[-2, 2] x)/Sqrt[g],
        A/(P Sqrt[g]), A/((P + nz[-2, 2]) Sqrt[g]), A/g^(3/2), A Sqrt[g]/P^2, A E^x/Sqrt[g],
        A/Sqrt[g] + Together[D[R Sqrt[g], x]], nz[1, 3] A/Sqrt[g] + Together[D[R, x]] Sqrt[g],
        A Log[x]/Sqrt[g], (x + nz[-2, 2]) A/((x + nz[-2, 2]) Sqrt[g])}]]];

(* B3: a radical of x that no linear-fractional or Euler substitution removes,
   with a Log/Exp kernel: D[Q] Log[S] with Q an algebraic function whose
   derivative is compact (a Chebyshev binomial or a radical of a cubic or quartic
   over x^k), and Log/Exp times the A3 shapes. *)
qalg[] := pick[{(x^3 + nz[-2, 2])^pick[{1/3, 2/3}]/x^r[1, 2], Sqrt[g3[]]/x^r[1, 2], (x^4 + nz[1, 3])^(1/4)/x, x (x^3 + nz[1, 3])^(1/3),
  Sqrt[g3[]] x, (x^2 + nz[1, 3])^(1/3)}];
grammar["B3"] := Function[{}, Module[{Q = qalg[], S = srat[], t = tk3[]},
  pick[{Together[D[Q, x]] Log[S], Together[D[Q, x]] Log[x], Together[D[Q, x]] Log[S] + Q D[S, x]/S, t chebyshevBinom[],
        Log[x] chebyshevBinom[], E^x (Q + Together[D[Q, x]]), Log[x + Sqrt[g3[]]]}]]];

(* B4: B3 drew few distinct integrands (110 unique for 118 cases in 262).  B4
   draws the derivative of a product that holds a radical y of a cubic or a
   quartic (no linear-fractional or Euler substitution removes it) at the bottom
   of a Log/Exp kernel, so the unperturbed draw is elementary by construction:
     D[R E^y]            E^y (R' + R y')          the radical inside Exp
     D[R Log[S + y]]     R' Log[S + y] + R (S' + y')/(S + y)
     D[P y E^(k x)]      E^(k x) (D[P y] + k P y)
     D[y Log[S]]         y' Log[S] + y S'/S
     D[R y Log[x]]       a radical, a rational and a Log factor at once
   with R, S rational, P polynomial.  None is theta' h(theta) for one kernel
   theta (x stays outside the kernel in R and S); a perturbation (one term
   dropped, the cofactor shifted) is elementary or not as Integrate decides. *)
y4[] := pick[{Sqrt[g3[]], Sqrt[g3[]], (x^3 + nz[-2, 2])^(1/3), (x^3 + x + nz[1, 2])^(1/3), (x^4 + nz[1, 3])^(1/4),
  Sqrt[x^3 + nz[-2, 2] x^2 + nz[1, 2]]}];
r4[] := pick[{x, 1/x, x + nz[-2, 2], 1/(x + nz[-2, 2]), x^2, 1/x^2, nz[1, 3], x/(x^2 + nz[1, 2])}];
s4[] := pick[{x, x^2, x + nz[-2, 2], nz[1, 3], x^2 + nz[1, 2], 1/x}];
grammar["B4"] := Function[{}, Module[{y = y4[], R = r4[], S = s4[], P = pick[{1, x, x^2, x + nz[-2, 2]}],
    k = pick[{1, -1, 2}], base, pert},
  base = pick[{D[R E^y, x], D[R E^y, x], D[R Log[S + y], x], D[R Log[S + y], x], D[P y E^(k x), x],
    D[y Log[S + nz[1, 2]], x], D[R y Log[x], x], D[x E^y, x] + D[R Log[x + y], x]}];
  base = Together[base];
  pert = pick[{0, 0, 0, 1, 2}];
  Which[pert === 0, base,
    pert === 1 && Head[Expand[Numerator[base]]] === Plus, Together[(Expand[Numerator[base]] /. Plus[a_, b__] :> Plus[b])/Denominator[base]],
    True, Together[base + nz[-2, 2] R base/(x + nz[1, 3])]]]];

(* an accepted reference must be elementary: an ALLOWLIST of heads, since a
   blocklist of special functions misses the ones nobody listed (Erfc and
   HypergeometricPFQ references once slipped into the corpus labelled elem) *)
elementaryHeads = {Plus, Times, Power, Log, Sqrt, Exp, Sin, Cos, Tan, Cot, Sec, Csc, Sinh, Cosh, Tanh, Coth, Sech, Csch,
  ArcTan, ArcTanh, ArcSin, ArcCos, ArcCot, ArcCoth, ArcSinh, ArcCosh, ArcSec, ArcCsc, ArcSech, ArcCsch, Root, RootSum, Function, Slot};
elementaryQ[F_] := Complement[Union[Cases[{F}, h_[___] :> h, {0, Infinity}]], elementaryHeads, {List, Rational, Complex}] === {};  (* Rational, Complex: atoms that Mathilda's h_[___] matches *)

(* ---------------- search ---------------- *)
emitCand[id_, f0_] := Module[{f = stripConstants[f0], pts, fv},   (* no constant factors: common.m *)
  If[!KeyExistsQ[emitted, f],
    pts = TimeConstrained[points[f], 30, {}];
    If[Length[pts] === 3,
      emitted[f] = True;
      fv = N[f /. x -> #, 40] & /@ pts;
      Print["CAND\t", id, "\t", cls, "\t", "c" <> cls, "\t", "elem", "\t", ToString[InputForm[f]], "\tNone\t",
        ToString[InputForm[pts]], "\t", ToString[InputForm[fv]], "\t", ToString[FullForm[f]], "\t", LeafCount[f]]]]];

accepted = 0; attempts = 0; ncand = 0; nescaped = 0; maxAttempts = attFactor quota;
While[accepted < quota && attempts < maxAttempts,
  attempts++;
  (* An abort that escapes its TimeConstrained (an Integrate that cached a
     value with Abort[] inside during an inner timeout and reused it later, a
     known evaluator bug, tasks/todo.md) costs this attempt, not the run. *)
  CheckAbort[Module[{f, F, Fc, nonelem = False, id},
    f = Quiet[TimeConstrained[compactForm[grammar[gname][]], 10, $Aborted]];
    If[f === $Aborted || FreeQ[f, x] || !contentOK[f] || LeafCount[f] > 120 || KeyExistsQ[emitted, f], Continue[]];
    If[syntacticEasyQ[f], Continue[]];
    F = Quiet[Check[TimeConstrained[Integrate[f, x], 30, $Aborted], nonelem = True; $Failed, {Integrate::nonelem}]];
    If[nonelem, Continue[]];
    id = StringJoin["c", cls, "-", ToString[seed], "-", ToString[attempts]];
    If[F === $Aborted || !FreeQ[F, Integrate],
      If[!ddEasyQ[f], emitCand[id, f]; ncand++]; Continue[]];
    (* fair mode (every RATIO but 2): an integrand Mathilda answers wrongly or
       only with special functions goes to FriCAS like an unanswered one, so a
       Mathilda miss is not silently left out of the corpus *)
    If[!elementaryQ[F], If[fair && !ddEasyQ[f], emitCand[id, f]; ncand++]; Continue[]];
    Fc = Quiet[TimeConstrained[compactForm[F], 10, F]];
    If[!ratioOK[LeafCount[Fc], LeafCount[f]], Continue[]];
    If[!TrueQ[Quiet[TimeConstrained[PossibleZeroQ[D[F, x] - f], 20, False]]],
      If[fair && !ddEasyQ[f], emitCand[id, f]; ncand++]; Continue[]];
    If[ddEasyQ[f], Continue[]];
    With[{before = Length[emitted]},
      emit[id, cls, "c" <> cls <> "-mathilda", "elem", f, F];
      If[Length[emitted] > before, accepted++]]],
    nescaped++; WriteString["stderr", "escaped abort at attempt ", attempts, "\n"]]];
WriteString["stderr", "class ", cls, ": accepted ", accepted, ", FriCAS candidates ", ncand, ", attempts ", attempts,
  ", escaped aborts ", nescaped, If[excludeEasy, ", easy excluded " <> ToString[InputForm[Thread[easyNames -> nEasy]]], ""], "\n"];
