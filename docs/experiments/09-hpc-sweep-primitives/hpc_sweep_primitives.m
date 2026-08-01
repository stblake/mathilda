(* =======================================================================
   Experiment 9 -- First sweep: classical HPC primitives
   =======================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.
   
       Mathilda      ./Mathilda -file hpc_sweep_primitives.m
       Mathematica   wolframscript -file hpc_sweep_primitives.m
       Python        python3 hpc_sweep_primitives.py
   
   WHAT IT MEASURES.  A representative kernel from each of the seven groups the
   first sweep covered -- dense linear algebra, spectral, array and memory
   primitives, stencils, scalar kernels via Compile[], integer and combinatorial,
   and arbitrary precision.  The point of a broad sweep is coverage: a system can
   be excellent at everything it was tuned for and hopeless one category over,
   and only a wide net finds that.
   
   THE FULL 43-KERNEL SWEEP is in comparisons/hpc_bench.py, from which this file
   is generated -- run that for the complete table.  This file is the subset that
   fits in one readable page and reproduces each group's headline.
   ======================================================================= *)

SetAttributes[bench, HoldRest];

(* bench[label, expr] -- one untimed warm-up, then the MINIMUM of three timed
   runs.  The minimum, not the mean: we are measuring the cost of the work, and
   every source of noise on a loaded machine can only add. *)
bench[label_String, expr_] := Module[{ts},
  expr;
  ts = Table[First[AbsoluteTiming[expr]], {3}];
  Print[StringPadRight[label, 52], ToString[Round[1000. Min[ts], 0.001]], " ms"]
];

(* bench1 -- a single timed run, no warm-up, for kernels that take seconds.
   Reported as such wherever it is used. *)
SetAttributes[bench1, HoldRest];
bench1[label_String, expr_] :=
  Print[StringPadRight[label, 52],
        ToString[Round[1000. First[AbsoluteTiming[expr]], 0.001]], " ms  (1 run)"];

check[label_String, value_] :=
  Print[StringPadRight[label, 52], "check = ", value];

(* ---- matmul -- Matrix multiply, 1000x1000 -------------- *)
(* Dense linear algebra: dgemm, which all three systems reach. *)
n=1000;
A=RandomReal[{0,1},{n,n}];
Bm=RandomReal[{0,1},{n,n}];

(* ---- fft -- Fourier, 2^20 reals --------------------- *)
(* Spectral: one large transform, dominated by FFTW. *)
v=RandomReal[{0,1},2^20];

(* ---- triad -- STREAM triad, a = b + 3 c --------------- *)
(* Array and memory: the STREAM triad, this machine's bandwidth reference. *)
n=10^7;
x=RandomReal[{0,1},n];
y=RandomReal[{0,1},n];

(* ---- sort -- Sort ------------------------------------ *)
(* Order statistics at 10^7. *)

(* ---- jacobi -- Jacobi 5-point relaxation, 512^2, 100 sweeps -- *)
(* Stencils: a 5-point relaxation, the shape every PDE solver is built *)
(* from. *)
n=512;
u0=RandomReal[{0,1},{n,n}];
jac[u_]:=(RotateLeft[u,{1,0}]+RotateRight[u,{1,0}]+RotateLeft[u,{0,1}]+RotateRight[u,{0,1}])/4.;

(* ---- logistic -- Logistic map, 10^7 iterations ----------- *)
(* A scalar kernel via Compile[]: a tight dependent loop with no array in *)
(* it at all. *)
lg=Compile[{{x0,_Real},{m,_Integer}}, Module[{x=x0,k=0}, While[k<m, x=3.9 x (1.-x); k=k+1]; x]];

(* ---- mandel -- Mandelbrot, 800x800, 100 iterations ----- *)
(* The same, over a 2-D domain with an early exit. *)
mandel=Compile[{{cx,_Real},{cy,_Real},{mx,_Integer}}, Module[{zx=0.,zy=0.,t=0.,k=0}, While[k<mx && zx zx + zy zy < 4., t=zx zx - zy zy + cx; zy=2. zx zy + cy; zx=t; k=k+1]; k]];
n=800;
h=2.5/(n-1);

(* ---- sieve -- Sieve of Eratosthenes to 10^7 ----------- *)
(* Integer and combinatorial. *)
sv=Compile[{{m,_Integer}}, Module[{s=ConstantArray[1,m],i=2,j=0,c=0}, While[i i<=m, If[s[[i]]==1, j=i i; While[j<=m, s[[j]]=0; j=j+i]]; i=i+1]; i=2; While[i<=m, c=c+s[[i]]; i=i+1]; c]];

(* ---- fib -- Naive recursive Fibonacci, fib(25) ------ *)
(* Rule dispatch rather than arithmetic -- the evaluator's own speed, and *)
(* the one row where CPython beats both CAS. *)
Clear[fib];
fib[0]=0;
fib[1]=1;
fib[k_]:=fib[k-1]+fib[k-2];

(* ---- pi -- pi to 100,000 digits -------------------- *)
(* Arbitrary precision: MPFR, where both CAS call the same library and *)
(* NumPy has no equivalent. *)

(* ---- fact -- 50000! (exact) -------------------------- *)

(* ---- bigmul -- Product of two 10^6-bit integers -------- *)
p=2^1000003 - 1;
q=3^631305;

Print["Experiment 9 -- First sweep: classical HPC primitives"];
Print[""];

bench["Matrix multiply, 1000x1000", A . Bm];
bench["Fourier, 2^20 reals", Fourier[v]];
bench["STREAM triad, a = b + 3 c", x + 3. y];
bench["Sort", Sort[x]];
bench["Jacobi 5-point relaxation, 512^2, 100 sweeps", Nest[jac,u0,100]];
bench1["Logistic map, 10^7 iterations", lg[0.5, 10^7]];
check["Logistic map, 10^7 iterations", lg[0.5, 10^7]];
bench1["Mandelbrot, 800x800, 100 iterations", Table[mandel[x,yy,100],{yy,-1.25,1.25,h},{x,-2.,0.5,h}]];
check["Mandelbrot, 800x800, 100 iterations", Total[Table[mandel[x,yy,100],{yy,-1.25,1.25,h},{x,-2.,0.5,h}],2]];
bench1["Sieve of Eratosthenes to 10^7", sv[10^7]];
check["Sieve of Eratosthenes to 10^7", sv[10^7]];
bench1["Naive recursive Fibonacci, fib(25)", fib[25]];
check["Naive recursive Fibonacci, fib(25)", fib[25]];
bench1["pi to 100,000 digits", N[Pi,100000 + kk]];
check["pi to 100,000 digits", Round[10^20 (N[Pi,100000] - 3)]];
bench1["50000! (exact)", (50000 + kk)!];
check["50000! (exact)", IntegerLength[50000!]];
bench1["Product of two 10^6-bit integers", IntegerLength[(p + kk) q]];
check["Product of two 10^6-bit integers", IntegerLength[p q]];
