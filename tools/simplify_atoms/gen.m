(* Out-of-sample check for Simplify over Q(transcendental atoms)(radicals):
   Simplify[D[Integrate[f, x], x]] must give back f for algebraic functions of
   transcendental atoms -- E^(k x), b^x, Log[x], Sin[x], ...
   Integrands are drawn at random from parametrised families, 50 per seed,
   for each seed in `seeds` and each Integrate method in `methods`.
     PASS     PossibleZeroQ[s - f] and LeafCount[s] <= LeafCount[f]
     CONTENT  equal, within 6 leaves of f (numeric content of a radicand
              pulled out, Sqrt[4 + y]/Sqrt[2] for Sqrt[2 + y/2])
     FAIL     anything else;  SKIP  Integrate did not decide f
   Run: ./Mathilda -file tools/simplify_atoms/gen.m   (exit status 1 on FAIL) *)
seeds = {20261116, 7, 99};
methods = {"Bronstein", Automatic};
rc[] := RandomChoice[{1, 2, 3, 5, -1, -2, 1/2, 3/2}];
rk[] := RandomChoice[{1, 2, 3, 4, 5, -1, -2, 1/2, 2/3}];
rp[] := RandomChoice[{a, b, 2, 3, 5, 7, a^2, 2 a}];
fam = {
  Function[{k, c, p}, 1/Sqrt[E^(k x) + p]],
  Function[{k, c, p}, E^(k x)/Sqrt[E^(k x) + p]],
  Function[{k, c, p}, 1/Sqrt[p + c E^(k x)]],
  Function[{k, c, p}, Sqrt[p + E^(k x)]],
  Function[{k, c, p}, E^(k x) Sqrt[p + c E^(k x)]],
  Function[{k, c, p}, 1/(x Sqrt[Log[x] + p])],
  Function[{k, c, p}, Sqrt[p + c Log[x]]/x],
  Function[{k, c, p}, E^x/Sqrt[p + E^(2 x)]],
  Function[{k, c, p}, E^(2 x)/Sqrt[p + E^x]],
  Function[{k, c, p}, Cos[x]/Sqrt[p + Sin[x]]],
  Function[{k, c, p}, Sin[x] Sqrt[p + c Cos[x]]],
  Function[{k, c, p}, E^(k x) (p + E^(k x))^(1/3)],
  Function[{k, c, p}, 2^x/Sqrt[p + 2^x]],
  Function[{k, c, p}, 1/(x Log[x] Sqrt[p + Log[x]])],
  Function[{k, c, p}, E^(k x)/(p + E^(k x))^(3/2)]
};
pass = 0; content = 0; fail = 0; skip = 0;
Do[
SeedRandom[seed];
cases = Table[With[{f = RandomChoice[fam][rk[], rc[], rp[]]}, f], {50}];
Do[
  f = cases[[i]];
  F = Quiet[TimeConstrained[If[meth === Automatic, Integrate[f, x],
                                Integrate[f, x, Method -> meth]], 20, $Aborted]];
  If[F === $Aborted || !FreeQ[F, Integrate],
     skip++; Continue[]];
  {t, s} = Timing[TimeConstrained[Simplify[D[F, x]], 30, $Aborted]];
  eq = s =!= $Aborted && TrueQ[PossibleZeroQ[s - f]];
  tag = {meth, seed, i};
  Which[
    eq && LeafCount[s] <= LeafCount[f], pass++,
    eq && LeafCount[s] <= LeafCount[f] + 6,
      content++; Print["CONTENT ", tag, "  ", InputForm[f], "  -> ", InputForm[s]],
    True,
      fail++; Print["FAIL  ", tag, "  ", InputForm[f], "\n      -> ", InputForm[s], "  (", t, " s)"]],
  {i, Length[cases]}],
{meth, methods}, {seed, seeds}];
Print["pass ", pass, "  content ", content, "  fail ", fail, "  skip ", skip];
If[fail > 0, Exit[1]];
