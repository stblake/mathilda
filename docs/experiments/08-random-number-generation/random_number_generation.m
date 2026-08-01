(* ==========================================================================
   Experiment 8 -- A machine-precision random number generator
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.

       Mathilda      ./Mathilda -file random_number_generation.m
       Mathematica   wolframscript -file random_number_generation.m
       Python        python3 random_number_generation.py

   WHAT IT MEASURES.  Bulk draws, and one Monte Carlo that consumes them.
   Random number generation is invisible in a benchmark suite until something
   draws 10^7 of them, at which point a generator that allocates a GMP integer
   per draw becomes the whole cost of every stochastic kernel in the system --
   which is exactly what the Black-Scholes and Monte-Carlo rows of the later
   sweeps do.

   NOTE ON THE MATHEMATICA COLUMN.  Wolfram's generator is a different
   algorithm with different statistical guarantees, so the comparison is of
   THROUGHPUT only.  The Monte Carlo row is the meaningful one: it consumes
   the draws and its answer must converge to Pi in every system.
   ========================================================================== *)

SetAttributes[bench, HoldRest];

(* bench[label, expr] -- one untimed warm-up, then the MINIMUM of three timed
   runs.  The minimum, not the mean: we are measuring the cost of the work, and
   every source of noise on a loaded machine can only add. *)
bench[label_String, expr_] := Module[{ts},
  expr;
  ts = Table[First[AbsoluteTiming[expr]], {3}];
  Print[StringPadRight[label, 52], ToString[Round[1000. Min[ts], 0.001]], " ms"]
];

check[label_String, value_] :=
  Print[StringPadRight[label, 52], "check = ", value];


Print["Experiment 8 -- random number generation"];
Print[""];
Print["-- bulk draws, 10^7 --"];

bench["RandomReal[{0,1}, 10^7]",        RandomReal[{0, 1}, 10^7]];
bench["RandomInteger[{0,100}, 10^7]",   RandomInteger[{0, 100}, 10^7]];
bench["RandomInteger[{0,255}, 10^7]",   RandomInteger[{0, 255}, 10^7]];
bench["RandomReal[{0,1}, {1000, 1000}]", RandomReal[{0, 1}, {1000, 1000}]];

Print[""];
Print["-- a Monte Carlo that consumes them --"];

mcpi[m_] := Module[{u, v},
  u = RandomReal[{0, 1}, m];
  v = RandomReal[{0, 1}, m];
  4. Total[UnitStep[1. - (u u + v v)]]/m];

bench["Monte Carlo pi, 10^7 samples", mcpi[10^7]];
Print[StringPadRight["Monte Carlo pi, 10^7 samples", 52],
      "value = ", mcpi[10^7], "   (converges to Pi, not exact)"];

Print[""];
Print["-- the draw is only useful if the distribution is right --"];
sq = RandomReal[{0, 1}, 10^6];
Print["mean     (want 0.5)       = ", Round[Total[sq]/10^6, 1.*^-4]];
Print["variance (want 1/12)      = ",
      Round[Total[(sq - Total[sq]/10^6)^2]/10^6, 1.*^-4]];
si = RandomInteger[{0, 9}, 10^6];
Print["integer mean (want 4.5)   = ", Round[N[Total[si]]/10^6, 1.*^-4]];
