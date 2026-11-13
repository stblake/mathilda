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
contentOK[f_] := If[gname === "B2", !FreeQ[f, Log | Power[E, _]] && !FreeQ[f, Power[_, _Rational]], True];

(* an accepted reference must be elementary: an ALLOWLIST of heads, since a
   blocklist of special functions misses the ones nobody listed (Erfc and
   HypergeometricPFQ references once slipped into the corpus labelled elem) *)
elementaryHeads = {Plus, Times, Power, Log, Sqrt, Exp, Sin, Cos, Tan, Cot, Sec, Csc, Sinh, Cosh, Tanh, Coth, Sech, Csch,
  ArcTan, ArcTanh, ArcSin, ArcCos, ArcCot, ArcCoth, ArcSinh, ArcCosh, ArcSec, ArcCsc, ArcSech, ArcCsch, Root, RootSum, Function, Slot};
elementaryQ[F_] := Complement[Union[Cases[{F}, h_[___] :> h, {0, Infinity}]], elementaryHeads, {List, Rational, Complex}] === {};  (* Rational, Complex: atoms that Mathilda's h_[___] matches *)

(* ---------------- search ---------------- *)
emitCand[id_, f_] := Module[{pts, fv},
  If[!KeyExistsQ[emitted, f],
    pts = TimeConstrained[points[f], 30, {}];
    If[Length[pts] === 3,
      emitted[f] = True;
      fv = N[f /. x -> #, 40] & /@ pts;
      Print["CAND\t", id, "\t", cls, "\t", "c" <> cls, "\t", "elem", "\t", ToString[InputForm[f]], "\tNone\t",
        ToString[InputForm[pts]], "\t", ToString[InputForm[fv]], "\t", ToString[FullForm[f]], "\t", LeafCount[f]]]]];

accepted = 0; attempts = 0; ncand = 0; maxAttempts = 40 quota;
While[accepted < quota && attempts < maxAttempts,
  attempts++;
  Module[{f, F, Fc, nonelem = False, id},
    f = Quiet[TimeConstrained[compactForm[grammar[gname][]], 10, $Aborted]];
    If[f === $Aborted || FreeQ[f, x] || !contentOK[f] || LeafCount[f] > 120 || KeyExistsQ[emitted, f], Continue[]];
    F = Quiet[Check[TimeConstrained[Integrate[f, x], 30, $Aborted], nonelem = True; $Failed, {Integrate::nonelem}]];
    If[nonelem, Continue[]];
    id = StringJoin["c", cls, "-", ToString[seed], "-", ToString[attempts]];
    If[F === $Aborted || !FreeQ[F, Integrate],
      emitCand[id, f]; ncand++; Continue[]];
    If[!elementaryQ[F], Continue[]];
    Fc = Quiet[TimeConstrained[compactForm[F], 10, F]];
    If[LeafCount[Fc] < 2 LeafCount[f], Continue[]];
    If[!TrueQ[Quiet[TimeConstrained[PossibleZeroQ[D[F, x] - f], 20, False]]], Continue[]];
    With[{before = Length[emitted]},
      emit[id, cls, "c" <> cls <> "-mathilda", "elem", f, F];
      If[Length[emitted] > before, accepted++]]]];
WriteString["stderr", "class ", cls, ": accepted ", accepted, ", FriCAS candidates ", ncand, ", attempts ", attempts, "\n"];
