(* association_bench.m — timing harness for the hash-backed Association ops.
 *
 * Run inside the REPL with:   Get["examples/association_bench.m"]
 * or from a shell:            ./Mathilda < examples/association_bench.m
 *
 * Companion to examples/association_bench.py and examples/association-benchmarks.md.
 * Counts / GroupBy / bulk Lookup are backed by an open-addressing hash index
 * (src/assoc.c), so each is amortised O(n). naiveCounts rebuilds the whole
 * association on every insert, so it is O(n^2) — included to show the gap the
 * bulk builtins close. *)

(* ---- Hash-backed frequency table: Counts, at increasing N ---- *)
Print["== Counts (hash-backed, O(n)) =="];
Do[
  SeedRandom[1];
  data = RandomInteger[1000, n];
  Print["  Counts   N=", n, ": ", 1000. Timing[Counts[data]][[1]], " ms"],
  {n, {20000, 40000, 80000, 160000}}];

(* ---- GroupBy: cost is dominated by evaluating the classifier f per element ---- *)
Print["== GroupBy (O(n) + n * cost of f) =="];
Print["  GroupBy  N=80000: ",
  1000. Timing[GroupBy[Range[80000], Mod[#, 100] &]][[1]], " ms"];

(* ---- Bulk Lookup: one index build, then O(1) per key ---- *)
SeedRandom[1];
big = Counts[RandomInteger[1000, 160000]];
keys = Keys[big];
Print["== Lookup of ", Length[keys], " keys (single index build, O(n+m)) =="];
Print["  Lookup   : ", 1000. Timing[Lookup[big, keys]][[1]], " ms"];

(* ---- Naive O(n^2) counting: rebuild the association on every insert ---- *)
(* Note: Lookup[a, x] (not a[[x]]) is used for key access — for integer keys
 * a[[x]] is *positional*, not a key lookup, in Mathilda as in Wolfram. *)
naiveCounts[list_] := Module[{a = <||>, x, i},
  Do[
    x = list[[i]];
    If[KeyExistsQ[a, x],
       AssociateTo[a, x -> Lookup[a, x] + 1],
       AssociateTo[a, x -> 1]],
    {i, Length[list]}];
  a];

Print["== naiveCounts (O(n^2), no bulk hashing) =="];
Do[
  SeedRandom[1];
  s = RandomInteger[200, n];
  Print["  naive    N=", n, ": ", 1000. Timing[naiveCounts[s]][[1]], " ms"],
  {n, {1000, 2000, 4000}}];
