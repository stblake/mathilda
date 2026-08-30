(* geometry_e2e.m -- GEO-1 end-to-end acceptance script.
   Run: ./Mathilda -file tests/scripts/geometry_e2e.m
   Evaluates every GEO-1 acceptance input through the full REPL pipeline
   (parser -> evaluator -> printer) and compares the printed form against
   the Wolfram-kernel-verified oracle strings (mathilda printer forms).
   Prints PASS/FAIL per row and a final summary line.
   CONTRACT: success is the literal marker line "geometry_e2e: ALL PASS".
   Callers MUST assert the marker (e.g. `./Mathilda -file ... | grep -q
   "geometry_e2e: ALL PASS"`): Quit[1] does NOT propagate a nonzero exit
   through -file mode (fault-injection receipt, 2026-08-27), so the process
   exit code alone is NOT a verdict. *)

errs = 0;
check[input_, expected_] := Module[{got = ToString[input]},
  If[got === expected,
    Print["PASS  ", expected],
    errs = errs + 1; Print["FAIL  expected: ", expected, "  got: ", got]]];

SetAttributes[check, HoldFirst];

check[Area[Polygon[{{0,0},{1,0},{1/2,1/2}}]], "1/4"];
check[Area[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}]], "10"];
check[Area[Polygon[{{0,0},{1.5,0},{1.5,1},{0,1}}]], "1.5"];
check[Area[Polygon[{{0,0},{1,0},{1,1},{0,1},{0,0}}]], "1"];
check[Area[Polygon[{{0,0},{1,0}}]], "Undefined"];
check[Perimeter[Polygon[{{0,0},{1,0},{0,1}}]], "2 + Sqrt[2]"];
check[Perimeter[Polygon[{{0,0},{3.,0},{3.,4.}}]], "12.0"];
check[Perimeter[Polygon[{{0,0},{1,0}}]], "Undefined"];
check[Perimeter[Polygon[{{0,0},{1/2,0},{0,1}}]], "3/2 + 1/2 Sqrt[5]"];
check[RegionCentroid[Polygon[{{0,0},{1,0},{0,1}}]], "{1/3, 1/3}"];
check[RegionCentroid[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}]], "{2, 7/5}"];
check[RegionMember[Polygon[{{0,0},{2,0},{2,2},{0,2}}], {1,1}], "True"];
check[RegionMember[Polygon[{{0,0},{2,0},{2,2},{0,2}}], {2,1}], "True"];
check[RegionMember[Polygon[{{0,0},{2,0},{2,2},{0,2}}], {0,0}], "True"];
check[RegionMember[Polygon[{{0,0},{2,0},{2,2},{0,2}}], {3,1}], "False"];
check[RegionMember[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}], {2,3}], "False"];
check[RegionMember[Polygon[{{0,0},{2,0},{2,2},{0,2}}], {1/2,1/2}], "True"];
check[ConvexHullRegion[{{0,0},{2,0},{1,0},{2,2},{0,2},{1,1}}], "Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}]"];
check[ConvexHullRegion[{{0,0},{1,0},{0,1},{2,2},{1/2,1/2}}], "Polygon[{{0, 0}, {1, 0}, {2, 2}, {0, 1}}]"];
check[ConvexHullRegion[{{0,0},{1,1},{2,2},{3,3}}], "Line[{{0, 0}, {3, 3}}]"];
check[ConvexHullRegion[{{1,2}}], "Point[{1, 2}]"];
check[Area[Polygon[{{0,0},{a,0},{0,1}}]], "Area[Polygon[{{0, 0}, {a, 0}, {0, 1}}]]"];
check[Area[ConvexHullRegion[{{0,0},{2,0},{1,0},{2,2},{0,2},{1,1}}]], "4"];
check[Attributes[Area], "{Protected, ReadProtected}"];
check[ConvexHullRegion[NDArray[{{0.,0.},{2.,0.},{2.,2.},{0.,2.},{1.,1.}}]], "Polygon[{{0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}}]"];
check[Area[Polygon[NDArray[{{0.,0.},{2.,0.},{2.,2.},{0.,2.}}]]], "4.0"];

(* Adversarial-review findings (2026-08-27): each was a wrong answer before. *)
check[Area[Polygon[{{0,0},{0,0},{1,1}}]], "Undefined"];
check[Perimeter[Polygon[{{0,0},{1,0},{1,0}}]], "Undefined"];
check[RegionMember[Polygon[{{0,0},{10^400,0},{0,1}}], {1,1}], "False"];
check[Area[Polygon[Table[{i, i^2}, {i, 0, 5}]]], "20"];
check[Area[Polygon[NDArray[{{0,0},{2,0},{2,2},{0,2}}]]], "4.0"];

If[errs === 0,
  Print["geometry_e2e: ALL PASS"],
  Print["geometry_e2e: ", errs, " FAILURES"]; Quit[1]];
