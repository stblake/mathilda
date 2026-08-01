(* =======================================================================
   Experiment 5 -- Machine integers: int64 as a peer of float64
   =======================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.
   
       Mathilda      ./Mathilda -file machine_integers.m
       Mathematica   wolframscript -file machine_integers.m
       Python        python3 machine_integers.py
   
   WHAT IT MEASURES.  An array system that can hold only machine REALS must
   choose, for every integer computation, between an exact answer and a fast one.
   These four kernels are integer end to end and none can accept a rounded
   answer, so they are the test of whether an int64 buffer type earns its place.
   
   THE EXACTNESS CONTRACT, which is what makes this hard: an int64 buffer must
   ABANDON the whole array on overflow rather than wrap, so the ordinary List
   path re-runs the operation and GMP answers exactly.  A wrapped sum that is
   merely fast is not a correct answer at any speed.
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

(* ---- sieve -- Sieve of Eratosthenes to 10^7 ----------- *)
(* A pure integer-array workload inside Compile[]: a flag vector, strided *)
(* writes, and a count. The compiled form is the one that was measured; an *)
(* interpreted strided Part assignment is a different and far slower *)
(* program. *)
sv=Compile[{{m,_Integer}}, Module[{s=ConstantArray[1,m],i=2,j=0,c=0}, While[i i<=m, If[s[[i]]==1, j=i i; While[j<=m, s[[j]]=0; j=j+i]]; i=i+1]; i=2; While[i<=m, c=c+s[[i]]; i=i+1]; c]];

(* ---- collatz -- Collatz longest chain below 10^6 -------- *)
(* A data-dependent scalar loop -- the case Compile[] exists for, and the *)
(* one where an integer type must not silently become a float: 3n+1 leaves *)
(* a double's exact range long before it leaves int64. *)
cz=Compile[{{m,_Integer}}, Module[{bl=0,k=1,q=0,len=0}, While[k<=m, q=k; len=0; While[q>1, If[Mod[q,2]==0, q=Quotient[q,2], q=3q+1]; len=len+1]; If[len>bl, bl=len]; k=k+1]; bl]];

(* ---- life -- Game of Life, 256^2, 100 generations ---- *)
(* Eight rotations, an integer sum and two comparisons. The neighbour count *)
(* is an int64 grid, so its UnitStep sees integers and not reals -- which *)
(* is what the narrowing kernel category of this experiment is for. *)
n=256;
g0=RandomInteger[1,{n,n}];
nb[g_]:=Sum[RotateLeft[g,{i,j}],{i,-1,1},{j,-1,1}]-g;
life[g_]:=With[{k=nb[g]}, UnitStep[k-3]UnitStep[3-k] + UnitStep[k-2]UnitStep[2-k] g];

(* ---- fib -- Naive recursive Fibonacci, fib(25) ------ *)
(* Not an array kernel at all: pure rule dispatch, and the control that *)
(* says how much of the above is the buffer and how much is the evaluator. *)
Clear[fib];
fib[0]=0;
fib[1]=1;
fib[k_]:=fib[k-1]+fib[k-2];

Print["Experiment 5 -- Machine integers: int64 as a peer of float64"];
Print[""];

bench1["Sieve of Eratosthenes to 10^7", sv[10^7]];
check["Sieve of Eratosthenes to 10^7", sv[10^7]];
bench1["Collatz longest chain below 10^6", cz[10^6]];
check["Collatz longest chain below 10^6", cz[10^6]];
bench["Game of Life, 256^2, 100 generations", Nest[life,g0,100]];
bench1["Naive recursive Fibonacci, fib(25)", fib[25]];
check["Naive recursive Fibonacci, fib(25)", fib[25]];
