(* gen_bottom.m -- benchmark 260, a radical at the bottom: integrands in
   Q(x)(y_1..y_m)(t_1..t_k), y_i simple radicals of x and the t_i logarithms and
   exponentials of elements of Q(x) AND of Q(x)(y) (a radical inside a kernel).
   This is class B of benchmark 256, the one class where FriCAS led at v0.6646;
   256 has only six B families and no ground truth for its perturbed ones
   (Campaign D0, src/internal/bronstein/BRONSTEIN_PLAN.md).

   Run:  COUNT=60 ./Mathilda -file tools/integration_shootout/gen_bottom.m \
           > benchmarks/260-radical-bottom/corpus.tsv

   Families (COUNT cases per half):
   - "<fam>-elem": f = D[F], F random in the family's field, so f has an
     elementary integral: a "not elementary" from any system is a false one;
   - "<fam>-rand": random integrands of the same field, expect "unknown";
   - "pert-nonelem": f = g + D[F] with g a classically non-elementary integrand
     over Sqrt[x] (E^(k x)/Sqrt[x + c] and E^(k x) Sqrt[x + c], erfi type;
     Sqrt[x]/Log[x] and 1/(Sqrt[x] Log[x]), li type) and F over a second,
     independent radical: f is not elementary (expect "nonelem"), and an answer
     from any system is a wrong one.
   Points as in 256 (gen.m): three rationals where every radicand and log
   argument is positive and f is real. *)
x =.;
count = ToExpression[Environment["COUNT"] /. $Failed -> "60"];
r[a_, b_] := RandomInteger[{a, b}];
nz[a_, b_] := RandomChoice[{-1, 1}] RandomInteger[{a, b}];
Get["tools/integration_shootout/common.m"];

(* the radicals: genus 0 (one square root of a quadratic, a square root of x),
   genus >= 1 (square roots of cubics and quartics), and a cube root *)
rad0 = {Sqrt[1 + x^2], Sqrt[2 + x^2], Sqrt[x], Sqrt[1 + x]};
rad1 = {Sqrt[1 + x^3], Sqrt[2 + x^3], Sqrt[1 + x^4], Sqrt[3 + x^4], (1 + x^3)^(1/3)};
rads = Join[rad0, rad1];
cf := RandomChoice[{-1, 1, 2, 1/2, 3}];
(* kernels over Q(x), and over Q(x)(y) for the radical y *)
kQ = {E^x, E^(2 x), Log[x], Log[1 + x], E^(x^2)};
kY[y_] := {Log[x + y], Log[1 + y], E^y, E^(x y)};
co[y_] := RandomChoice[{1, y, 1/y, x y, y/(1 + x), x/y, y^2}];
term[y_, ks_] := cf x^r[0, 1] co[y] RandomChoice[Append[ks, 1]] RandomChoice[Append[ks, 1]];

fams = {
  (* one radical, kernels over Q(x) only *)
  {"one", Function[{}, With[{y = RandomChoice[rads]}, Sum[term[y, kQ], {r[1, 2]}]]],
          Function[{}, With[{y = RandomChoice[rads]}, RandomChoice[{1, x, 1/x, 1/(1 + x)}] co[y] RandomChoice[kQ]]]},
  (* one radical inside a kernel, further monomials above it (Phase 5) *)
  {"mid", Function[{}, With[{y = RandomChoice[rads]}, Sum[term[y, Join[kQ, kY[y]]], {r[1, 2]}]]],
          Function[{}, With[{y = RandomChoice[rads]}, RandomChoice[{1, x, 1/x}] co[y] RandomChoice[kY[y]] RandomChoice[kQ]]]},
  (* two independent radicals, a term over each and a mixed one *)
  {"two", Function[{}, With[{ys = RandomSample[rads, 2]},
            term[ys[[1]], kQ] + term[ys[[2]], Join[kQ, kY[ys[[2]]]]] + RandomChoice[{0, 1}] cf ys[[1]] ys[[2]] RandomChoice[kQ]]],
          Function[{}, With[{ys = RandomSample[rads, 2]},
            co[ys[[1]]] RandomChoice[kQ]/RandomChoice[{1, x, 1 + x}] + co[ys[[2]]] RandomChoice[Join[kQ, kY[ys[[2]]]]]/RandomChoice[{1, x, 1 + x}]]]}
};

nonelemG := With[{k = RandomChoice[{1, 2, -1}], c = RandomChoice[{0, 1, 2}]},
  RandomChoice[{E^(k x)/Sqrt[x + c], E^(k x) Sqrt[x + c], Sqrt[x]/Log[x], 1/(Sqrt[x] Log[x])}]];
pertF := With[{y = RandomChoice[{Sqrt[1 + x^2], Sqrt[2 + x^2], Sqrt[1 + x^3], Sqrt[2 + x^3], Sqrt[1 + x^4]}]},
  Sum[term[y, Join[kQ, kY[y]]], {r[1, 2]}]];

Do[
  Module[{nm = fam[[1]], k, tries, f, F, id, before},
    SeedRandom[26000 + StringLength[nm] 7 + First[ToCharacterCode[nm]]];
    {k, tries} = {0, 0};
    While[k < count && tries < 8 count, tries++;
      F = fam[[2]][]; f = Quiet[D[F, x]];
      If[f =!= 0 && Head[f] =!= D, id = nm <> "-elem-" <> ToString[k + 1];
        before = Length[emitted]; emit[id, "B", nm <> "-elem", "elem", f, F];
        If[Length[emitted] > before, k++]]];
    {k, tries} = {0, 0};
    While[k < count && tries < 8 count, tries++;
      f = Quiet[fam[[3]][]];
      If[!FreeQ[f, x], id = nm <> "-rand-" <> ToString[k + 1];
        before = Length[emitted]; emit[id, "B", nm <> "-rand", "unknown", f, None];
        If[Length[emitted] > before, k++]]]],
  {fam, fams}];

Module[{k = 0, tries = 0, f, id, before},
  SeedRandom[26099];
  While[k < count && tries < 8 count, tries++;
    f = nonelemG + Quiet[D[pertF, x]];
    id = "pert-nonelem-" <> ToString[k + 1];
    before = Length[emitted]; emit[id, "B", "pert-nonelem", "nonelem", f, None];
    If[Length[emitted] > before, k++]]];
