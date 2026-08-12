(* ==========================================================================
   benchmarks/harness.m -- the ONLY copy of the reporting helpers.
   ==========================================================================

   Loaded by every experiment as the first line:

       Get["../harness.m"];

   `Get` resolves relative to the CURRENT DIRECTORY, not the script's directory
   ($InputFileName and DirectoryName do not exist in Mathilda).  run_all.py
   therefore sets cwd to the experiment's folder before launching either CAS.
   To run one by hand, cd into the folder first:

       cd benchmarks/01-integrate-rational
       ../../Mathilda    -file integrate_rational.m
       wolframscript     -file integrate_rational.m

   This file is loaded UNMODIFIED by both Mathilda and Mathematica, so the two
   CAS columns cannot be timing two different reporting methods.

   ---------------------------------------------------------------------------
   THE EMIT CONTRACT.  One tagged tab-separated line per CASE -- a case being one
   named kernel, timed in every system and joined across them by its label.
   run_all.py joins
   the `.m` and `.py` halves of an experiment on <label>, so the two files must
   spell each label identically -- the runner fails the experiment if they do
   not, because silent case-drift is how a table starts comparing two different
   programs.

       BENCH   <label>  <milliseconds>
       CHECK   <label>  <value>
       REQUIRE <head>   present|absent
       SKIP    <label>  <head>            a case not run, because <head> is absent

   Tags rather than the older StringPadRight[label, 52] column format, which is
   parseable only by splitting on whitespace and hoping no label contains a
   number.  run_all.py renders the readable table from these tags, so nothing is
   lost by emitting them.
   ========================================================================== *)

(* ---- timing method ------------------------------------------------------ *)

(* One untimed warm-up, then the MINIMUM of $BenchReps timed runs.  The
   minimum, not the mean: we are measuring the cost of the work, and every
   source of noise on a loaded machine can only add.

   AbsoluteTiming, never Timing[] -- Timing sums CPU time across threads and so
   overstates every threaded NDArray path and every BLAS call by roughly the
   core count. *)

$BenchReps = 3;

SetAttributes[bench, HoldRest];

(* TWO SHARP EDGES, both found by this harness emitting garbage cases, both
   worth reading before editing anything below.

   1. The two arities do NOT delegate to each other.  `HoldRest` holds every
      argument after the first, so `bench[label, expr] := bench[label, expr,
      $BenchReps]` would pass the *symbol* `$BenchReps` rather than 3, the
      `reps_Integer` pattern would not match, and the 2-argument form would
      silently emit no case at all.

   2. The repetition count must be a GLOBAL symbol or a literal -- never a
      Module local.  `Table[e, {k}]` with `k` local to an enclosing Module
      returns UNEVALUATED in Mathilda (a literal, a global, and a
      pattern-substituted integer all work).  Writing the obvious
      `Module[{reps = $BenchReps}, Table[..., {reps}]]` produced cases reading
      `0.001 Round[1e+06 Table[0.0, {reps}]]`.  That is a Mathilda bug, tracked
      in this week's changelog; `$BenchReps` stays global to route around it. *)

bench[label_String, expr_] := Module[{ts},
  expr;                                          (* warm-up, untimed *)
  ts = Table[First[AbsoluteTiming[expr]], {$BenchReps}];
  Print["BENCH\t", label, "\t", ToString[Round[1000. Min[ts], 0.001]]];
];

bench[label_String, expr_, reps_Integer] := Module[{ts},
  expr;                                          (* warm-up, untimed *)
  ts = Table[First[AbsoluteTiming[expr]], {reps}];
  Print["BENCH\t", label, "\t", ToString[Round[1000. Min[ts], 0.001]]];
];

(* A case whose single run is already expensive (a 10^6-digit N[Pi], a 1050^2
   matmul x12) gets reps = 1 and NO warm-up: the warm-up would double a
   multi-second case for no information.  Used by 31-wolframmark. *)
SetAttributes[benchOnce, HoldRest];
benchOnce[label_String, expr_] :=
  Print["BENCH\t", label, "\t",
        ToString[Round[1000. First[AbsoluteTiming[expr]], 0.001]]];

(* ---- the value gate ---------------------------------------------------- *)

(* A cheap scalar recorded from every system and compared by run_all.py.  A
   timing case is MEANINGLESS until the systems agree on the check: the check is
   what pins the algorithm, so the columns cannot silently be timing three
   different computations.  A case whose check disagrees is reported CHECK-FAIL
   and its timing is DISCARDED, not shown with a caveat. *)

(* InputForm IS LOAD-BEARING, not decoration.  Bare `ToString[1/2]` renders in
   `wolframscript` as Wolfram's 2-D typeset fraction --

       1
       -
       2

   -- three lines, and run_all.py parses line by line, so the check value
   arrived as "1".  Mathilda prints "1/2" on one line, so the two systems
   "disagreed" on a value they both computed as Rational[1, 2], and a correct
   `Limit` was reported CHECK-FAIL with its timing discarded.  InputForm is
   guaranteed single-line ASCII in both systems and round-trips exactly. *)

check[label_String, value_] :=
  Print["CHECK\t", label, "\t", ToString[value, InputForm]];

(* ---- the coverage / absence register ----------------------------------- *)

(* `require` declares a head the experiment depends on, and reports whether
   this system has it.  run_all.py turns these into the coverage percentage and
   into ABSENT.md.  An ABSENT head NEVER gets a ratio attached: "we do not have
   this function" is feature work, and pooling it into a speed number would
   turn a missing algorithm into a fake 40x. *)

require[head_String] :=
  Print["REQUIRE\t", head, "\t", If[Names[head] === {}, "absent", "present"]];

require[heads_List] := Scan[require, heads];

haveQ[head_String] := Names[head] =!= {};

(* Run a case only if the head exists; otherwise emit SKIP so the runner can
   tell "absent" apart from "silently never measured".  This is the difference
   between an honest absence register and a table with holes in it. *)
SetAttributes[benchIf, HoldRest];
benchIf[label_String, head_String, expr_] :=
  If[haveQ[head],
     bench[label, expr],
     Print["SKIP\t", label, "\t", head]];

benchIf[label_String, head_String, expr_, reps_Integer] :=
  If[haveQ[head],
     bench[label, expr, reps],
     Print["SKIP\t", label, "\t", head]];

SetAttributes[checkIf, HoldRest];
checkIf[label_String, head_String, value_] :=
  If[haveQ[head], check[label, value], Null];

(* ---- portability shim, and why it exists ------------------------------- *)

(* WolframMark spells its random data `RandomReal[{}, dims]` -- an EMPTY range,
   which Mathematica reads as the default {0, 1}.  Mathilda returns that form
   UNEVALUATED, so every downstream Fourier / Eigenvalues / Dot / Transpose /
   SVD / LinearSolve in a verbatim port silently operates on a symbolic head
   and returns in ~0.1 ms -- a fictitious 10-100x "win" over Mathematica that
   measures nothing at all.  Found by the value gate during research.
   Tracked in ABSENT.md as a real compatibility gap.

   rand01 is the one spelling both systems agree on, so the benchmark measures
   arithmetic instead of measuring that bug.  The bug itself is still reported:
   31-wolframmark calls require["RandomRealEmptyRange"] via a probe case. *)

rand01[dims_] := RandomReal[{0, 1}, dims];

(* ---- deterministic seeding --------------------------------------------- *)

(* Every experiment seeds before generating, so a rerun on the same host
   reproduces.  The three systems cannot be made to draw the SAME stream, which
   is why value checks either use exactly-reproducible integer data or a
   separate small deterministic case -- never a float sum over random data. *)

seed[] := SeedRandom[1];
