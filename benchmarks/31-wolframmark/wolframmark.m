(* ==========================================================================
   Experiment 31 -- WolframMark
   ==========================================================================

   The fifteen tests of Wolfram Research's own WolframMark suite, at their
   official sizes and repetition counts, transcribed from the `Benchmarking`
   package's test list.  Runs unmodified in Mathilda and in Mathematica.

       cd benchmarks/31-wolframmark
       ../../Mathilda -file wolframmark.m
       wolframscript  -file wolframmark.m
       python3 wolframmark.py

   WHAT THIS SUITE ACTUALLY IS -- read before quoting any number from it.
   WolframMark is Wolfram's HARDWARE benchmark.  It holds Mathematica constant
   and varies the machine, reporting a score relative to a reference system:
   the question it was built to answer is "how fast is this laptop at running
   Mathematica", NOT "how good is this CAS".  Using it cross-implementation is a
   repurposing, and three things follow:

     1. It is NOT a coverage measure.  The fifteen tests were chosen to exercise
        hardware dimensions -- FPU throughput, memory bandwidth, cache, bignum --
        not to span a CAS's functionality.  Group A/B/C do that job here.
     2. The sizes encode Mathematica's performance profile.  Each test is tuned
        to take ~1 s in Mathematica on reference hardware, which says nothing
        about how long the same work takes in another implementation.  Concrete
        case: the 1200000-point x 11 DFT below is a ~1 s test in Mathematica and
        would have needed roughly FORTY HOURS against Mathilda's pre-FFTW
        O(n^2) Fourier fallback.  A test being "one second" is a fact about
        Mathematica-on-hardware, not about the kernel.
     3. The aggregate WolframMark SCORE is meaningless across implementations
        and is deliberately not computed here.  These are 15 more cases, joined
        and classified like every other case; nothing is summed into a score.

   WHAT IT IS STILL GOOD FOR.  The workload SELECTION is a third party's and
   predates this project, so unlike groups A-C it cannot be accused of having
   been picked to flatter us.  That is the whole of its value here.

   TWO DELIBERATE DEVIATIONS, both forced, both documented:

   1. `Benchmarking`Private`BenchmarkTiming[...]` is replaced by this corpus's
      `benchOnce[label, ...]`.  Same measurement -- one timed run of the whole
      body, wall clock -- but it emits the tagged line run_all.py joins on.  The
      benchmark *code* inside is character-for-character the official test.

   2. THE `RandomReal[{}, dims]` BUG.  Every official test spells its random
      data `RandomReal[{}, dims]` -- an EMPTY range, which Mathematica reads as
      the default {0, 1}.  Mathilda returns that form UNEVALUATED.  A verbatim
      port therefore hands a symbolic `RandomReal[{}, {1050, 1050}]` to Dot,
      Fourier, Eigenvalues, Transpose, SingularValueDecomposition and
      LinearSolve, all of which return in ~0.1 ms having computed nothing -- a
      fictitious 10-100x "win" over Mathematica.  So the timed rows use
      `rand01[dims]` from the harness (`RandomReal[{0, 1}, dims]`), the one
      spelling both systems agree on, and the bug itself is measured as its own
      row below so it appears in the report rather than being quietly worked
      around.

   VALUE CHECKS.  The timed data is random and the three systems cannot be made
   to draw the same stream, so no check is ever a float sum over the timed
   input.  Each check instead runs the SAME operation over a small exactly
   reproducible case from data.m.  A timing is meaningless until the systems
   agree on that.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"FindFit", "AccuracyGoal", "Fourier", "Eigenvalues", "DiagonalMatrix",
         "SingularValueDecomposition", "LinearSolve", "NIntegrate", "Expand",
         "Sort", "Transpose", "Gamma", "ArcTan", "SeedRandom", "RandomReal",
         "RandomInteger", "Do", "N"}];

(* ---- 0. The compatibility bug, as its own row --------------------------- *)

(* Not an official test.  It is here because the port could not be verbatim,
   and a deviation that is not measured is a deviation that gets forgotten.
   Mathilda leaves this unevaluated, so its "timing" is ~0 and its check
   disagrees with the other two systems -- run_all.py reports CHECK-FAIL and
   DISCARDS the timing, which is the correct outcome. *)
bench["RandomReal empty-range form", RandomReal[{}, {200000}];];
check["RandomReal empty-range form", Head[RandomReal[{}, {4}]] === List];

(* ---- 1. Data Fitting --------------------------------------------------- *)

(* Nonlinear parameter fit.  `Fit` is linear-basis fitting and cannot stand in
   for it, so on a system without FindFit this emits SKIP and lands in the
   absence register with no ratio attached. *)
benchIf["Data Fitting", "FindFit",
  Module[{data},
    data = Flatten[Table[{x, y, z,
        Log[120 x] - Abs[Cos[z/300]/(140 y)]},
        {x, 0.2, 10, 0.22}, {y, 0.2, 10, 0.22}, {z, 0.2, 10, 0.22}], 2];
    FindFit[data, Log[a x] - Abs[Cos[b z]/(c y)], {a, b, c}, {x, y, z},
            AccuracyGoal -> 6];], 1];

(* ---- 2. Digits of Pi --------------------------------------------------- *)

benchOnce["Digits of Pi", N[Pi, 1000000];];
check["Digits of Pi", Round[10^20 (N[Pi, 40] - 3)]];

(* ---- 3. Discrete Fourier Transform ------------------------------------- *)

benchOnce["Discrete Fourier Transform",
  Module[{data}, seed[]; data = rand01[{1200000}];
    Do[Fourier[data], {11}]]];
(* Deterministic check: the DFT of a small exact sequence, first component. *)
check["Discrete Fourier Transform",
  Round[10^6 Re[First[Fourier[{1., 2., 3., 4., 5., 6., 7., 8.}]]]]];

(* ---- 4. Eigenvalues of a Matrix ---------------------------------------- *)

(* Note a^-1 is the ELEMENTWISE reciprocal, not Inverse -- Power is Listable.
   That is what the official test does; it is a similarity transform only by
   accident of notation, and copying it exactly is the point. *)
benchOnce["Eigenvalues of a Matrix",
  Module[{a, b, m}, seed[];
    a = rand01[{420, 420}];
    b = DiagonalMatrix[rand01[{420}]];
    m = a.b.a^-1;
    Do[Eigenvalues[m], {6}]]];
check["Eigenvalues of a Matrix",
  Round[10^6 Total[Sort[Re[Eigenvalues[N[hilbertLike[6]]]]]]]];

(* ---- 5. Elementary Functions ------------------------------------------- *)

benchOnce["Elementary Functions",
  Module[{m1, m2}, seed[];
    m1 = rand01[{2200000}]; m2 = rand01[{2200000}];
    Do[Exp[m1]; Sin[m1]; ArcTan[m1, m2], {30}]]];
check["Elementary Functions",
  Round[10^6 (Exp[0.5] + Sin[0.5] + ArcTan[0.5, 1.5])]];

(* ---- 6. Gamma Function ------------------------------------------------- *)

(* 55 gammas of ~85000: exact integers of ~400000 digits each.  This is a
   bignum test wearing a special-function hat. *)
benchOnce["Gamma Function",
  Module[{a}, seed[]; a = RandomInteger[{80000, 90000}, {55}]; Gamma[a];]];
check["Gamma Function", Gamma[21]];

(* ---- 7. Large Integer Multiplication ----------------------------------- *)

benchOnce["Large Integer Multiplication",
  Module[{a}, seed[];
    a = RandomInteger[{10^1100000, 10^(1100000 + 1)}];
    Do[a (a + 1), {20}]]];
check["Large Integer Multiplication", Mod[(10^50 + 1) (10^50 + 2), 10^12]];

(* ---- 8. Matrix Arithmetic ---------------------------------------------- *)

benchOnce["Matrix Arithmetic",
  Module[{m}, seed[]; m = rand01[{840, 840}];
    Do[(1. + 0.5 m)^127, {50}];]];
check["Matrix Arithmetic", Round[10^6 ((1. + 0.5 0.25)^127)]];

(* ---- 9. Matrix Multiplication ------------------------------------------ *)

benchOnce["Matrix Multiplication",
  Module[{m1, m2}, seed[];
    m1 = rand01[{1050, 1050}]; m2 = rand01[{1050, 1050}];
    Do[m1.m2, {12}]]];
check["Matrix Multiplication",
  Round[10^6 Total[Flatten[N[spdIntMatrix[8]].N[spdIntMatrix[8]]]]]];

(* ---- 10. Matrix Transpose ---------------------------------------------- *)

benchOnce["Matrix Transpose",
  Module[{m}, seed[]; m = rand01[{2070, 2070}];
    Do[Transpose[m], {40}]]];
check["Matrix Transpose", Total[Flatten[Transpose[{{1, 2, 3}, {4, 5, 6}}]]]];

(* ---- 11. Numerical Integration ----------------------------------------- *)

benchOnce["Numerical Integration",
  NIntegrate[Sin[x^2 + y^2], {x, -2.6 Pi, 2.6 Pi}, {y, -2.6 Pi, 2.6 Pi}];];
(* A 1-D integral with an exact answer, so the check does not inherit the 2-D
   cubature's accuracy warning. *)
check["Numerical Integration", Round[10^4 NIntegrate[Sin[x], {x, 0, Pi}]]];

(* ---- 12. Polynomial Expansion ------------------------------------------ *)

benchOnce["Polynomial Expansion",
  Expand[Times @@ Table[(c + x)^3, {c, 350}]];];
check["Polynomial Expansion",
  Length[Expand[Times @@ Table[(c + x)^3, {c, 12}]]]];

(* ---- 13. Random Number Sort -------------------------------------------- *)

benchOnce["Random Number Sort",
  Module[{a}, seed[]; a = RandomInteger[{1, 50000}, {520000}];
    Do[Sort[a], {15}]]];
check["Random Number Sort", Total[Take[Sort[lcgInts[2000, 50000]], 10]]];

(* ---- 14. Singular Value Decomposition ---------------------------------- *)

benchOnce["Singular Value Decomposition",
  Module[{m}, seed[]; m = rand01[{860, 860}];
    Do[SingularValueDecomposition[m], {2}]]];
check["Singular Value Decomposition",
  Round[10^4 Total[Flatten[N[Part[SingularValueDecomposition[
    N[spdIntMatrix[6]]], 2]]]]]];

(* ---- 15. Solving a Linear System --------------------------------------- *)

benchOnce["Solving a Linear System",
  Module[{m, v}, seed[];
    m = rand01[{1150, 1150}]; v = rand01[{1150}];
    Do[LinearSolve[m, v], {16}]]];
check["Solving a Linear System",
  Round[10^6 Total[LinearSolve[N[spdIntMatrix[8]],
                               N[Table[i, {i, 8}]]]]]];
