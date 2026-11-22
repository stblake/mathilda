(* m2_tower_stress.m -- the false-certificate stress of Phase 5 M2 (the general
   tower over the curve, src/internal/bronstein/BaseTower.m).

   Every integrand is f = D[F] with F random in Q(x)(y)(t1)(t2)(t3): y a simple
   radical, the t_i logarithms and exponentials of elements of Q(x), of Q(x)(y),
   and of the field below them (Log[x + E^(x y)], E^(Log[x + y]^2), ...), at
   least one kernel with a radical argument.  So f HAS an elementary integral and
   any "not elementary" is a FALSE CERTIFICATE.  Each case runs
   BronsteinIntegrate (Method -> "Bronstein") under TimeConstrained; an answer is
   verified with PossibleZeroQ.

   Run (sharded):  N=330 SHARDS=6 SHARD=k TC=30 ./Mathilda -file \
                     tools/bronstein_stress/m2_tower_stress.m
   One TSV line per case on stdout: id, verdict, monomials, seconds, f.
   verdict: answer | WRONG | FALSECERT | incomplete | timeout | other. *)
LoadModule["bronstein/Bronstein.m"];
$ContextPath = Join[$ContextPath, {"Bronstein`Private`", "Trager`Private`"}];
env[k_, d_] := ToExpression[Environment[k] /. $Failed -> d];
n = env["N", "330"]; shards = env["SHARDS", "1"]; shard = env["SHARD", "0"]; tc = env["TC", "30"];
SeedRandom[20261122];
rads = {Sqrt[1 + x^2], Sqrt[1 + x^3], Sqrt[2 + x^3], Sqrt[1 + x^4], (1 + x^3)^(1/3), Sqrt[x], Sqrt[1 + x]};
cf := RandomChoice[{-1, 1, 2, 1/2, 3, -2/3}];
co[y_] := RandomChoice[{1, y, 1/y, x y, y/(1 + x), x/y, y^2, x, 1/x}];
(* kernels in three tiers: over Q(x), over Q(x)(y), over a kernel below *)
k0[y_] := RandomChoice[{Log[x], Log[1 + x], E^x, E^(x^2)}];
k1[y_] := RandomChoice[{Log[x + y], Log[1 + y], E^y, E^(x y), Log[1 + x y]}];
(* an exponential linear in a logarithm collapses to a power (E^(x Log[a]) = a^x), so
   over a logarithm only non-linear exponents *)
k2[y_, k_] := RandomChoice[Join[{Log[x + k], Log[1 + x k], E^(k^2), Log[1 + k^2]},
  If[Head[k] === Log, {E^(1/(1 + k))}, {E^(x k), E^(y k)}]]];
mkF[] := Module[{y = RandomChoice[rads], a, b, c, ks, mono},
  a = k1[y]; b = RandomChoice[{k0[y], k1[y]}]; c = k2[y, RandomChoice[{a, b}]];
  ks = {a, b, c};
  mono := cf co[y] Times @@ RandomSample[Append[ks, 1], RandomInteger[{1, 3}]];
  Total[Table[mono, {RandomInteger[{1, 3}]}]] + cf a b c/RandomChoice[{1, x, 1 + x}]];
cases = Table[mkF[], {n}];
Do[If[Mod[i - 1, shards] === shard,
    Module[{F = cases[[i]], f, t0, r, v, nm, T},
      f = Quiet[D[F, x]];
      T = Quiet[Check[Catch[brTower[f, x], _], $Failed]];
      nm = If[AssociationQ[T], ToString[T["shape"]] <> ":" <> ToString[Length[T["monos"]]], "decline"];
      t0 = AbsoluteTime[];
      r = Quiet[TimeConstrained[BronsteinIntegrate[f, x], tc, "TIMEOUT"]];
      v = Which[
        r === "TIMEOUT", "timeout",
        ListQ[r] && r =!= {} && First[r] === "not elementary", "FALSECERT",
        ListQ[r] && r =!= {} && StringQ[First[r]] && StringContainsQ[First[r], "incomplete"], "incomplete",
        ListQ[r] && r =!= {} && StringQ[First[r]], "other:" <> First[r],
        TrueQ[Quiet[TimeConstrained[PossibleZeroQ[D[r, x] - f], 60, False]]], "answer",
        True, "WRONG"];
      Print[StringRiffle[{"m2-" <> ToString[i], v, nm, ToString[NumberForm[AbsoluteTime[] - t0, {6, 2}]],
        ToString[f, InputForm]}, "\t"]]]],
  {i, n}];
