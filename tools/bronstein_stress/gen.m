(* gen.m -- seeded generator of Shape A stress integrands for the Bronstein
   package (BRONSTEIN_PLAN.md Phase 9).  Every case is f = D[F, x] for a random
   elementary F, so the correct verdict is "elementary": a "not elementary"
   answer is a FALSE CERTIFICATE and a wrong antiderivative is WRONG.

   Run:  SEED=1 COUNT=40 ./Mathilda -file tools/bronstein_stress/gen.m
   Prints one line per case:  CASE<TAB>id<TAB>family<TAB>InputForm[f]          *)
seed = ToExpression[Environment["SEED"]]; count = ToExpression[Environment["COUNT"]];
If[!IntegerQ[seed], seed = 1]; If[!IntegerQ[count], count = 40];
SeedRandom[seed];
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "integration_shootout", "families.m"}]];
families = mixedFamilies;
(* FAMILY=random: random R(x, theta) Y^k, mostly NOT elementary -- the verdicts
   are cross-checked against FriCAS (fricas_check.py), not assumed *)
(* Every other case is PERTURBED: D[F] of a known-elementary F plus a small
   component (x Y, Y/x, theta Y), so the integrand looks like an elementary one
   and the verdict is usually "not elementary" -- the near-miss class a sloppy
   certificate gets wrong. *)
If[Environment["FAMILY"] === "random",
  Do[Module[{th = theta[], n = RandomChoice[{2, 2, 3}], Y, f, fam},
     Y = lin[th]^(1/n);
     If[OddQ[i],
       f = (r[-2, 2] + r[-1, 1] x + nz[-2, 2] th)/(r[1, 2] x + RandomChoice[{0, 1}] th^r[1, 2] + nz[1, 3]) Y^RandomChoice[{1, -1}],
       fam = RandomChoice[Take[families, 6]];
       With[{F = fam[[2]][th, n]},
         (* the perturbation on F's own radical (a second radical is out of scope) *)
         Y = First[Cases[F, Power[b_, q_Rational] :> b^(1/Denominator[q]), {0, Infinity}], Y];
         f = Together[D[F, x] + RandomChoice[{x Y, Y/x, th Y, Y/(x + 1)}]]]];
     Print["CASE\t", seed, "-r", i, "\t", If[OddQ[i], "random", "perturb"], "\t", ToString[InputForm[f]]]],
    {i, count}];
  Quit[]];
Do[Module[{fam = RandomChoice[families], th = theta[], n = RandomChoice[{2, 2, 2, 3, 3, 4}], F, f},
   F = fam[[2]][th, n];
   f = Together[D[F, x]];
   Print["CASE\t", seed, "-", i, "\t", fam[[1]], "\t", ToString[InputForm[f]], "\t", ToString[InputForm[F]]]],
  {i, count}];
