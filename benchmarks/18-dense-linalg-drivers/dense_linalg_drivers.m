(* Experiment 18 -- Dense linear algebra drivers.
   MEASURES src/linalg/: lapack_bridge.c, blas_bridge.c and the eigen/QR/SVD
   kernels.  On this host NumPy links the SAME Accelerate BLAS Mathilda does, so
   on these rows all three systems call byte-identical kernels and any spread is
   PURE OVERHEAD -- the cleanest measurement in the whole suite.
   Roadmap item 5 ("reach the remaining LAPACK drivers") lives here.
   Checks use an exactly reproducible SPD integer matrix from data.m. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Dot", "Inverse", "Det", "LinearSolve", "Eigenvalues", "Eigenvectors",
         "SingularValueDecomposition", "QRDecomposition", "LUDecomposition",
         "MatrixRank", "PseudoInverse", "Norm"}];

n = 600;
a = rand01[{n, n}];
b = rand01[{n, n}];
v = rand01[{n}];
s8 = N[spdIntMatrix[8]];

bench["Dot 600x600", a.b;];
check["Dot 600x600", Round[10^6 Total[Flatten[s8.s8]]]];

bench["LinearSolve 600x600", LinearSolve[a, v];];
check["LinearSolve 600x600",
  Round[10^6 Total[LinearSolve[s8, N[Table[i, {i, 8}]]]]]];

bench["Inverse 600x600", Inverse[a];];
check["Inverse 600x600", Round[10^6 Total[Flatten[Inverse[s8]]]]];

bench["Det 600x600", Det[a];];
check["Det 600x600", Round[10^4 Det[s8]]];

bench["Eigenvalues 600x600 (general)", Eigenvalues[a], 1];
check["Eigenvalues 600x600 (general)",
  Round[10^4 Total[Sort[Re[Eigenvalues[s8]]]]]];

bench["SingularValueDecomposition 400x400",
  SingularValueDecomposition[rand01[{400, 400}]], 1];
check["SingularValueDecomposition 400x400",
  Round[10^4 Total[Flatten[Part[SingularValueDecomposition[s8], 2]]]]];

bench["QRDecomposition 400x400", QRDecomposition[rand01[{400, 400}]], 1];
check["QRDecomposition 400x400",
  Round[10^4 Abs[Total[Flatten[Part[QRDecomposition[s8], 2]]]]]];
