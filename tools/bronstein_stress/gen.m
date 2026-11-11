(* gen.m -- seeded generator of Shape A stress integrands for the Bronstein
   package (BRONSTEIN_PLAN.md Phase 9).  Every case is f = D[F, x] for a random
   elementary F, so the correct verdict is "elementary": a "not elementary"
   answer is a FALSE CERTIFICATE and a wrong antiderivative is WRONG.

   Run:  SEED=1 COUNT=40 ./Mathilda -file tools/bronstein_stress/gen.m
   Prints one line per case:  CASE<TAB>id<TAB>family<TAB>InputForm[f]          *)
seed = ToExpression[Environment["SEED"]]; count = ToExpression[Environment["COUNT"]];
If[!IntegerQ[seed], seed = 1]; If[!IntegerQ[count], count = 40];
SeedRandom[seed];
r[a_, b_] := RandomInteger[{a, b}];
nz[a_, b_] := RandomChoice[DeleteCases[Range[a, b], 0]];
theta[] := RandomChoice[{Log[x], Log[x + r[1, 3]], Log[nz[1, 3] x + r[1, 3]], E^x, E^(nz[-2, 2] x),
   E^(x^2), Log[1 + E^x], Log[x^2 + 1]}];
lin[th_] := nz[1, 3] th + r[-2, 2] x + nz[-3, 3];
families = {
  (* algebraic part only *)
  "alg" -> Function[{th, n}, With[{Y = lin[th]^(1/n)}, (r[1, 3] x + r[-2, 2] th + r[-1, 1]) Y^RandomChoice[{1, -1, n - 1}]]],
  (* log of x + Y: residues at finite places *)
  "logxy" -> Function[{th, n}, With[{Y = lin[th]^(1/n)}, r[0, 2] x Y + nz[-2, 2] Log[x + r[1, 2] + Y]]],
  (* arcsinh: residues at the special places *)
  "asinh" -> Function[{th, n}, With[{Y = Sqrt[th^2 + nz[1, 3]]}, nz[-2, 2] Log[th + Y] + r[0, 1] Y]],
  (* genus 1: Log[p + Y] on Y^2 = th^3 + c *)
  "genus1" -> Function[{th, n}, With[{Y = Sqrt[th^3 + nz[1, 3]]}, Log[r[1, 2] + Y] + r[0, 1] th Y]],
  (* algebraic residues: Log of a quadratic in Y over Q(sqrt) *)
  "algres" -> Function[{th, n}, With[{Y = Sqrt[th^2 + nz[1, 2]]}, Log[(th + 1 + Y)/(th - 1 + Y)] + r[0, 1] Y/th]],
  (* rational multiple of a component: b(th) Y with b rational in th *)
  "ratb" -> Function[{th, n}, With[{Y = lin[th]^(1/n)}, Y/(th + nz[1, 3]) + r[0, 1] Y^(n - 1)]]
};
(* FAMILY=random: random R(x, theta) Y^k, mostly NOT elementary -- the verdicts
   are cross-checked against FriCAS (fricas_check.py), not assumed *)
If[Environment["FAMILY"] === "random",
  Do[Module[{th = theta[], n = RandomChoice[{2, 2, 3}], Y, f},
     Y = lin[th]^(1/n);
     f = (r[-2, 2] + r[-1, 1] x + nz[-2, 2] th)/(r[1, 2] x + RandomChoice[{0, 1}] th^r[1, 2] + nz[1, 3]) Y^RandomChoice[{1, -1}];
     Print["CASE\t", seed, "-r", i, "\trandom\t", ToString[InputForm[f]]]],
    {i, count}];
  Quit[]];
Do[Module[{fam = RandomChoice[families], th = theta[], n = RandomChoice[{2, 2, 2, 3, 3, 4}], F, f},
   F = fam[[2]][th, n];
   f = Together[D[F, x]];
   Print["CASE\t", seed, "-", i, "\t", fam[[1]], "\t", ToString[InputForm[f]]]],
  {i, count}];
