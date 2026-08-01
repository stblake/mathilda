(* ==========================================================================
   Experiment 4 -- Auto-compilation: the compiler runs without Compile[]
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.

       Mathilda      ./Mathilda -file auto_compilation.m
       Mathematica   wolframscript -file auto_compilation.m
       Python        python3 auto_compilation.py

   WHAT IT MEASURES.  Ordinary user code -- no Compile[] anywhere -- that the
   system recognises as numeric and compiles on its own.  Two families:

     * the numeric BUILTINS, which evaluate a user function many times at
       machine precision (NIntegrate, Plot, FindRoot, NSum, NProduct, NDSolve)
     * the FUNCTIONAL heads, where the function is applied element by element
       (Nest, NestList, FoldList, Scan, Accumulate)

   In both cases the win is the same one experiment 1 measured, and the point
   of THIS experiment is that the user never asked for it.

   Every timing below can be compared against the un-compiled path by setting
   the environment variable MATHILDA_NO_AUTOCOMPILE=1 and re-running, which is
   how the README's before/after column was produced.
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


Print["Experiment 4 -- auto-compilation"];
Print[""];
Print["-- numeric builtins over a user function --"];

f[x_] := Sin[x] Exp[-x/10.] + Sqrt[x + 1.];

bench["NIntegrate[f, {x, 0, 50}]",        NIntegrate[f[x], {x, 0, 50}]];
check["NIntegrate[f, {x, 0, 50}]",        Round[NIntegrate[f[x], {x, 0, 50}], 1.*^-6]];

bench["NSum[f[k]/k^2, {k, 1, 2000}]",     NSum[f[k]/k^2, {k, 1, 2000}]];
bench["NProduct[1 + 1/k^2, {k, 1, 2000}]", NProduct[1 + 1./k^2, {k, 1, 2000}]];
bench["FindRoot[f[x] == 2, {x, 1}]",      FindRoot[f[x] == 2, {x, 1}]];

Print[""];
Print["-- functional heads over a machine-number list --"];

lst = RandomReal[{0, 1}, 10^6];

bench["Nest[3.5 # (1-#) &, 0.31, 10^6]",   Nest[3.5 # (1. - #) &, 0.31, 10^6]];
bench["NestList[3.5 # (1-#) &, 0.31, 10^6]", NestList[3.5 # (1. - #) &, 0.31, 10^6]];
bench["FoldList[#1 + Sin[#2] &, 0., list]", FoldList[#1 + Sin[#2] &, 0., lst]];
bench["Scan[Sin[#]^2 + 1. &, list]",        Scan[Sin[#]^2 + 1. &, lst]];
bench["Accumulate[list]",                   Accumulate[lst]];
bench["NestWhileList[# + 1. &, 0., # < 10^6 &]",
      NestWhileList[# + 1. &, 0., # < 1.*^6 &]];

Print[""];
Print["-- an ODE and an interpolation, both of which evaluate a user body --"];

Clear[ly];
bench["NDSolve, Lorenz to t = 200",
      NDSolve[{lx'[t] == 10.(ly[t] - lx[t]),
               ly'[t] == lx[t](28. - lz[t]) - ly[t],
               lz'[t] == lx[t] ly[t] - (8./3.) lz[t],
               lx[0] == 1., ly[0] == 1., lz[0] == 1.},
              {lx, ly, lz}, {t, 0, 200}]];

ni = 10^4;
tbl = Table[{xx, Sin[xx]}, {xx, 0., 100., 100./(ni - 1)}];
ifn = Interpolation[tbl];
pts = RandomReal[{0, 100}, ni];
bench["Interpolation, 10^4 nodes, 10^4 evaluations", ifn /@ pts];
check["Interpolation at 3.3", Round[ifn[3.3], 1.*^-6]];

Print[""];
Print["re-run with MATHILDA_NO_AUTOCOMPILE=1 for the interpreted column"];
