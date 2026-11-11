(* valgrind_dsolve_only.m -- the valgrind_sample.txt records through DSolve
 * ALONE (no verifier), $pdReps times: isolates leaks in the PDE methods from
 * leaks in the verifier's PossibleZeroQ / Simplify / TimeConstrained probes.
 * Compare totals for $pdReps = 1 and 3, as in valgrind_sample.m. *)
If[!IntegerQ[$pdReps], $pdReps = 1];
pdAll = Join @@ (Get /@ FileNames["DSolve_test_status/pde/DE_PDE_*.m"]);
pdLabels = ReadList["DSolve_test_status/pde/valgrind_sample.txt", String];
pdSample = Select[pdAll, MemberQ[pdLabels, #[[1]]] &];
(* NO TimeConstrained: a timeout's siglongjmp abandons the expression nodes of
 * the aborted computation by design (time-and-date.md), and under valgrind's
 * 30-50x slowdown an 8 s budget fires on most records -- the totals would then
 * measure the documented abort leak, not the methods.  The whole sample runs
 * in ~8 s natively. *)
Do[Quiet[DSolve[pdRec[[3]], pdRec[[4]], pdRec[[5]], Sequence @@ pdRec[[6]]]],
   {$pdReps}, {pdRec, pdSample}];
Print["records: ", Length[pdSample], "  reps: ", $pdReps];
