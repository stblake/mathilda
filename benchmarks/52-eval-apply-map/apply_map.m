(* Apply / MapApply / Map fast-path benchmark -- runs UNCHANGED in Mathilda and
   WolframKernel. Symbolic data throughout so both stay on their general
   evaluators (ff/gg/hh/aa/... are all undefined). Lists are built ONCE, outside
   timing, so each case measures only the Apply/Map restructuring. Emits bare
   integer-microsecond timings (in order), then bare check values (in order). *)

$Reps = 5;
SetAttributes[bench, HoldAll];
bench[e_] := Module[{ts}, e; ts = Table[First[AbsoluteTiming[e]], {$Reps}];
  Round[1.*^6 Min[ts]]];

symP = Table[c[k] x, {k, 5000}];                                 (* 5000 symbolic terms c[k] x *)
tri  = Table[{aa[k], bb[k], cc[k]}, {k, 5000}];                  (* 5000 triples *)
comp = Table[gg[aa[k], bb[k], cc[k]], {k, 5000}];               (* compound elements (depth 2) *)
deep = Table[gg[hh[aa[k], bb[k]], hh[cc[k], dd[k]]], {k, 5000}]; (* deep elements (depth 3) *)

(* ---- timings (integer microseconds, in order) ---- *)
Print[bench[ff @@ symP]];              (* 1. Apply @@  (level 0)            *)
Print[bench[ff @@@ tri]];              (* 2. MapApply @@@ (level 1)         *)
Print[bench[ff /@ comp]];              (* 3. Map /@ compound (level 1)      *)
Print[bench[ff /@ deep]];              (* 4. Map /@ deep (level 1)          *)
Print[bench[Map[ff, comp, {0}]]];      (* 5. Map at level 0                 *)
Print[bench[Scan[ff, comp]]];          (* 6. Scan compound (level 1)        *)
Print[bench[Scan[ff, deep]]];          (* 7. Scan deep (level 1)            *)
Print[bench[MapIndexed[ff, comp]]];    (* 8. MapIndexed compound (level 1)  *)
Print[bench[MapIndexed[ff, deep]]];    (* 9. MapIndexed deep (level 1)      *)
Print[bench[Apply[Plus, tri, {1}]]];   (* 10. MapApply Plus (level 1)       *)
Print[bench[Length[Plus @@ symP]]];    (* 11. Plus @@ grouping (5000 terms) *)

(* ---- checks (bare values, in order; must agree across engines) ---- *)
Print[Length[ff @@ symP]];
Print[LeafCount[ff @@@ tri]];
Print[LeafCount[ff /@ comp]];
Print[LeafCount[ff /@ deep]];
Print[LeafCount[Map[ff, comp, {0}]]];
Print[Module[{s = 0}, Scan[(s = s + 1) &, comp]; s]];
Print[Module[{s = 0}, Scan[(s = s + 1) &, deep]; s]];
Print[LeafCount[MapIndexed[ff, comp]]];
Print[LeafCount[MapIndexed[ff, deep]]];
Print[LeafCount[Apply[Plus, tri, {1}]]];
Print[Length[Plus @@ symP]];
