(* tests/bronstein_bench/capture.m -- freeze the inputs of the ideal and field
   operations (Phase 4.5, step 1).  Runs the integrands below with every
   operation of opNames wrapped: each call is timed and its arguments kept;
   the 3 slowest calls per (operation, field type) are written to
   tests/bronstein_bench/ideal_inputs.m for ideal_ops.m.
   ./Mathilda -file tests/bronstein_bench/capture.m   (about 10 minutes) *)
Integrate[Sqrt[x], x, Method -> "Bronstein"];
capOps = {"idealTimes", "idealTimesGenerators", "idealPower", "idealInverse", "idealReduce", "hnfModDet", "hnf",
  "principalGenerator", "normalizeAtInfinity", "kNorm", "kExtGCD", "residueIdeal", "afMultMatrix", "logDerivativeNF"};
capKeep = 3; capCalls = <||>; capIn = <||>;
capFieldType[F_Association] := StringJoin[F["Type"], If[KeyExistsQ[F, "t"], "+ext", ""], If[KeyExistsQ[F, "consts"], "+consts", ""], If[KeyExistsQ[F, "I"], "+I", ""]];
capFieldType[_] := "?";
Do[With[{s = Symbol["Trager`Private`" <> nm], name = nm},
    Unprotect[s];
    s[args___] /; !TrueQ[capIn[name]] := Module[{t0 = AbsoluteTime[], r, dt, key},
      capIn[name] = True;
      r = s[args];
      capIn[name] = False;
      dt = AbsoluteTime[] - t0;
      key = {name, capFieldType[First[{args}]]};
      If[dt > 0.005,
        capCalls[key] = Take[Reverse[SortBy[Append[Lookup[capCalls, Key[key], {}], {dt, {args}}], First]], UpTo[capKeep]]];
      r];
    DownValues[s] = RotateRight[DownValues[s], 1]],
  {nm, capOps}];
capCases = {
  {"G-perf-2", "Bronstein", Together[D[Log[Log[x] + (1 + Log[x] + Log[x]^3)^(1/3)] + 2 Log[Log[x]^2 + Log[x] (1 + Log[x] + Log[x]^3)^(1/3) + (1 + Log[x] + Log[x]^3)^(2/3)], x]], 200},
  {"G-perf-3", "Bronstein", Tan[x]^(1/4), 200},
  {"G-perf-4", "Trager", x/((2^(1/3) x^3 - 1) Sqrt[x^2 + 2]), 60},
  {"E3-1", "Bronstein", 1/(x (1 + Log[x]^3)^(1/3)), 60},
  {"E3-2", "Bronstein", E^x/(1 + E^(3 x))^(1/3), 60},
  {"E3-5", "Bronstein", (1 + Log[x])/(x (1 + Log[x] + Log[x]^3)^(1/3)), 60},
  {"T-1", "Trager", 1/((x - Sqrt[2]) Sqrt[x^3 + 1]), 60},
  {"T-2", "Trager", 1/((x^2 + I) Sqrt[x + 1]), 60},
  {"T-3", "Trager", (1 + (1 + I) x^4)^(3/4), 60},
  {"T-4", "Trager", x/Sqrt[x^4 + 10 x^2 - 96 x - 71], 60},
  {"T-5", "Trager", 1/((x - 2^(1/3)) Sqrt[x^2 + 1]), 60}};
Do[Module[{t = AbsoluteTiming[TimeConstrained[
      If[c[[2]] === "Trager", Trager`TragerIntegrate[c[[3]], x], Bronstein`BronsteinIntegrate[c[[3]], x]], c[[4]], "TIMEOUT"]]},
    capIn = <||>;
    Print[c[[1]], "  ", N[t[[1]], 3], "  ", If[ListQ[t[[2]]] || t[[2]] === "TIMEOUT", t[[2]], "answer"]]],
  {c, capCases}];
capOut = Flatten[KeyValueMap[Function[{k, v}, <|"op" -> k[[1]], "field" -> k[[2]], "time" -> #[[1]], "args" -> #[[2]]|> & /@ v], capCalls]];
Put[capOut, "tests/bronstein_bench/ideal_inputs.m"];
Print["CAPTURED ", Length[capOut], " calls"];
Do[Print[r["op"], "  ", r["field"], "  ", N[r["time"], 3]], {r, SortBy[capOut, {#["op"] &, -#["time"] &}]}];
