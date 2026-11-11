(* valgrind_sample.m -- run the PDE corpus check on a fixed cross-section of
 * records, $pdReps times, for the leak audit:
 *
 *   valgrind --leak-check=full ./Mathilda -file DSolve_test_status/pde/valgrind_sample.m
 *
 * Mathilda's startup has a fixed one-time "definitely lost" baseline (dyld /
 * Accelerate / FLINT init, ~12.8 KB in 400 blocks on macOS).  A per-call leak
 * makes the total GROW with the number of calls, so compare a run with
 * $pdReps = 1 against $pdReps = 3 (and a trivial Print[1] script): equal totals
 * mean the PDE methods and the verifier are leak-free.  The sample is the
 * labels in valgrind_sample.txt (one per line), chosen to touch every PDE method
 * on a PASS record plus a spread of declines. *)
Get["DSolve_test_status/pde/dsolve_pde_corpus_prelude.m"];
If[!IntegerQ[$pdReps], $pdReps = 1];
pdAll = Join @@ (Get /@ FileNames["DSolve_test_status/pde/DE_PDE_*.m"]);
pdLabels = ReadList["DSolve_test_status/pde/valgrind_sample.txt", String];
pdSample = Select[pdAll, MemberQ[pdLabels, #[[1]]] &];
Print["records: ", Length[pdSample], "  reps: ", $pdReps];
pdTally = Table[pdeCheckCode @@ pdRec[[1 ;; 6]], {$pdReps}, {pdRec, pdSample}];
Print["verdicts: ", Tally[Flatten[pdTally]]];
