(* ==========================================================================
   Experiment 1 -- Compile[]: a typed bytecode VM for numeric code
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.

       Mathilda      ./Mathilda -file compile_engine.m
       Mathematica   wolframscript -file compile_engine.m
       Python        python3 compile_engine.py

   WHAT IT MEASURES.  The same four kernels twice: once as ordinary
   interpreted expressions, once through Compile[].  Both CAS have a Compile[]
   with the same surface syntax, so the two columns compare two
   implementations of the same idea rather than two different ideas.

   WHY THESE FOUR.  Each is a tight loop with no array in it -- a scalar
   recurrence, a 2-D domain with an early exit, an all-pairs sum, and a
   vectorised Monte Carlo.  That is the regime a bytecode VM exists for: when
   the work per evaluation step is one arithmetic operation, the interpreter's
   per-step cost is the whole cost.
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

(* ---- 1. logistic map: a scalar recurrence, 10^7 iterations -------------- *)
(* Nothing to vectorise: each iterate depends on the last. *)
logi = Compile[{{n, _Integer}},
  Module[{x = 0.31, k = 0},
    While[k < n, x = 3.9 x (1. - x); k = k + 1]; x]];

logiInterp[n_] := Module[{x = 0.31, k = 0},
  While[k < n, x = 3.9 x (1. - x); k = k + 1]; x];

(* ---- 2. Mandelbrot: a 2-D domain with a data-dependent early exit ------- *)
mand = Compile[{{w, _Integer}, {maxit, _Integer}},
  Module[{c = 0, i = 0, j = 0, zr = 0., zi = 0., cr = 0., ci = 0., k = 0, t = 0.},
    i = 0;
    While[i < w,
      j = 0;
      While[j < w,
        cr = -2. + 3. i/w; ci = -1.5 + 3. j/w;
        zr = 0.; zi = 0.; k = 0;
        While[k < maxit && zr zr + zi zi < 4.,
          t = zr zr - zi zi + cr; zi = 2. zr zi + ci; zr = t; k = k + 1];
        c = c + k; j = j + 1];
      i = i + 1];
    c]];

(* ---- 3. Lennard-Jones energy: an all-pairs scalar sum ------------------- *)
ljc = Compile[{{n, _Integer}},
  Module[{e = 0., i = 1, j = 1, dx = 0., dy = 0., dz = 0., r2 = 0., ir6 = 0.},
    i = 1;
    While[i <= n,
      j = i + 1;
      While[j <= n,
        dx = Sin[1. i] - Sin[1. j];
        dy = Cos[1.3 i] - Cos[1.3 j];
        dz = Sin[0.7 i] - Sin[0.7 j];
        r2 = dx dx + dy dy + dz dz + 0.5;
        ir6 = 1./(r2 r2 r2);
        e = e + 4.(ir6 ir6 - ir6);
        j = j + 1];
      i = i + 1];
    e]];

(* ---- 4. Monte Carlo pi, vectorised ------------------------------------- *)
(* Included because it is the one of the four that is NOT a scalar loop: it
   is a compiled expression over whole arrays, which is experiment 2's
   subject and the reason the two experiments sit next to each other. *)
mcpi[m_] := Module[{u, v}, u = RandomReal[{0, 1}, m]; v = RandomReal[{0, 1}, m];
  4. Total[UnitStep[1. - (u u + v v)]]/m];

(* ---- run ---------------------------------------------------------------- *)

Print["Experiment 1 -- Compile[]: a typed bytecode VM"];
Print[""];

bench["logistic map, 10^7 iterations (compiled)",   logi[10^7]];
check["logistic map, 10^7 iterations",              logi[10^7]];
bench["logistic map, 10^5 iterations (interpreted)", logiInterp[10^5]];

bench["Mandelbrot, 800^2, 100 iterations",          mand[800, 100]];
check["Mandelbrot, 800^2, 100 iterations",          mand[800, 100]];

bench["Lennard-Jones energy, 1452 bodies",          ljc[1452]];
check["Lennard-Jones energy, 1452 bodies",          ljc[1452]];

bench["Monte Carlo pi, 10^7 samples (vectorised)",  mcpi[10^7]];

Print[""];
Print["the interpreted logistic map runs 100x fewer iterations, so divide the"];
Print["compiled number by 100 before comparing -- that ratio is the VM's worth"];
