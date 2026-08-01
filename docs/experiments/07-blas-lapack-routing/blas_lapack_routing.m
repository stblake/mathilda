(* =======================================================================
   Experiment 7 -- Reaching the vendor BLAS and LAPACK kernels
   =======================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.
   
       Mathilda      ./Mathilda -file blas_lapack_routing.m
       Mathematica   wolframscript -file blas_lapack_routing.m
       Python        python3 blas_lapack_routing.py
   
   WHAT IT MEASURES.  Seven dense linear-algebra calls.  On this host all three
   systems link the SAME Apple Accelerate BLAS, so on the rows that reach it the
   three columns are running byte-identical kernels and any spread is pure
   overhead -- which makes this the cleanest possible test of whether a fast path
   is actually being reached.
   
   Where a row is far from the other two, the kernel is NOT being reached: that
   is the finding, and README.md names which ones and why.
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
(* dgemm -- the row every system reaches, and therefore the calibration for *)
(* the rest. *)
n=1000;
A=RandomReal[{0,1},{n,n}];
Bm=RandomReal[{0,1},{n,n}];

(* ---- solve -- LinearSolve, 1000x1000 ------------------ *)
(* dgesv. *)
n=1000;
A=RandomReal[{0,1},{n,n}];
bv=RandomReal[{0,1},n];

(* ---- inverse -- Inverse, 500x500 ------------------------ *)
(* dgetrf + dgetri. *)
n=500;
A=RandomReal[{0,1},{n,n}];

(* ---- det -- Det, 500x500 ---------------------------- *)
(* An LU factorisation and a product of the diagonal. *)

(* ---- qr -- QRDecomposition, 500x500 ---------------- *)
(* dgeqrf + dorgqr. *)

(* ---- eigen -- Eigenvalues, 300x300 symmetric ---------- *)
(* dsyevd for a symmetric matrix. The ordering convention is the reason *)
(* this one is hard to route: LAPACK's order is not Wolfram's, and parity *)
(* must be preserved. *)
n=300;
A=RandomReal[{0,1},{n,n}];
S=(A+Transpose[A])/2.;

(* ---- svd -- SingularValueDecomposition, 300x300 ----- *)
(* dgesdd. *)
n=300;
A=RandomReal[{0,1},{n,n}];

Print["Experiment 7 -- Reaching the vendor BLAS and LAPACK kernels"];
Print[""];

bench["Matrix multiply, 1000x1000", A . Bm];
bench["LinearSolve, 1000x1000", LinearSolve[A,bv]];
bench["Inverse, 500x500", Inverse[A]];
bench["Det, 500x500", Det[A]];
bench["QRDecomposition, 500x500", QRDecomposition[A]];
bench["Eigenvalues, 300x300 symmetric", Eigenvalues[S]];
bench["SingularValueDecomposition, 300x300", SingularValueDecomposition[A]];
