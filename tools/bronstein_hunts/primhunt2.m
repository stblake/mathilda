(* seeded Liouvillian-primitive hunt: f = D[F]; a False is a false certificate
   only when every Liouvillian kernel of F also occurs in f *)
SeedRandom[3];
rr[] := RandomInteger[{-3, 3}];
rq[] := RandomChoice[{1, 2, -1, 1/2, 3}];
rp[d_] := Sum[rr[] x^k, {k, 0, d}];
rat[] := rp[RandomInteger[{0, 2}]]/(x + RandomInteger[{1, 3}])^RandomInteger[{0, 1}];
lk = (Erf | Erfi | Erfc | ExpIntegralEi | LogIntegral)[_];
kern[] := RandomChoice[{Erf[x], Erfi[x], Erf[2 x], ExpIntegralEi[x], ExpIntegralEi[2 x], Erfc[x], Erf[x + 1]}];
mkF["poly"] := rat[] kern[] + rat[];
mkF["sq"] := Module[{k = kern[]}, rq[] k^2 + rat[] k];
mkF["exp"] := rat[] E^(rr[] x) kern[] + rat[] E^x;
mkF["log"] := rat[] Log[x] kern[] + rat[] kern[];
mkF["two"] := rat[] kern[] kern[] + rat[];
fams = {"poly", "sq", "exp", "log", "two"};
stats = <||>;
Do[Do[
  Module[{F, f, r, key, kF, kf},
    F = mkF[fam]; f = Together[D[F, x]];
    If[FreeQ[F, lk] || PossibleZeroQ[f], Continue[]];
    kF = Union[Cases[{F}, lk, Infinity]]; kf = Union[Cases[{f}, lk, Infinity]];
    r = TimeConstrained[Quiet[Risch`ElementaryIntegralQ[f, x]], 60, $Aborted];
    key = Which[r === True, "True", r === False,
        If[SubsetQ[kf, kF], "FALSE-CERT", "False(ok: kernel not in f)"],
      r === $Aborted, "timeout", True, "undecided"];
    If[key === "FALSE-CERT" || key === "undecided", Print[key, " ", fam, ": F = ", InputForm[F]]];
    stats[{fam, key}] = Lookup[stats, Key[{fam, key}], 0] + 1],
  {trial, 24}], {fam, fams}];
Scan[Print[#[[1]], " -> ", #[[2]]] &, SortBy[Normal[stats], First]];
