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

(* ---- the decompositions and matrix functions added in v0.6180-v0.6188 ------

   Timing on random matrices (shifted to be SPD / Hermitian positive definite
   / spectrum in the right half-plane where the head needs it); checks on the
   small reproducible matrices below, as order- and sign-invariant scalars.
   The matrix functions (MatrixLog, MatrixFunction, non-integer MatrixPower)
   are timed at 100 x 100: the reference column costs ~10 s at 200 x 200 for
   each of them, which would spend the file's budget on four cases. *)

c4 = {{4., 1. + I, 0., 0.5 I}, {1. - I, 5., 2., 0.}, {0., 2., 6., 1. - 0.5 I},
      {-0.5 I, 0., 1. + 0.5 I, 7.}};                     (* Hermitian PD *)
g4 = N[{{3, 1, 0, 2}, {2, 4, 1, 0}, {0, 1, 5, 1}, {1, 0, 1, 6}}];

spd6 = a . Transpose[a] + n IdentityMatrix[n];
h4 = With[{c = rand01[{400, 400}] + I rand01[{400, 400}]},
          c . ConjugateTranspose[c] + 400 IdentityMatrix[400]];
g3 = rand01[{300, 300}];
c3 = rand01[{300, 300}] + I rand01[{300, 300}];
g2 = rand01[{200, 200}];
p1 = rand01[{100, 100}] + 100 IdentityMatrix[100];      (* spectrum near 100 *)
r1 = rand01[{100, 100}]/10;
m10 = rand01[{10, 10}];

require[{"CholeskyDecomposition", "HessenbergDecomposition", "Eigensystem",
         "SingularValueList", "Minors", "Adjugate", "Orthogonalize",
         "MatrixExp", "MatrixLog", "MatrixFunction", "MatrixPower"}];

benchIf["CholeskyDecomposition 600x600", "CholeskyDecomposition",
  CholeskyDecomposition[spd6];];
checkIf["CholeskyDecomposition 600x600", "CholeskyDecomposition",
  Round[10^6 Total[Flatten[CholeskyDecomposition[s8]]]]];

benchIf["CholeskyDecomposition 400x400 complex Hermitian", "CholeskyDecomposition",
  CholeskyDecomposition[h4];];
checkIf["CholeskyDecomposition 400x400 complex Hermitian", "CholeskyDecomposition",
  Round[10^6 Re[Total[Flatten[CholeskyDecomposition[c4]]]]]];

(* H's off-diagonal signs follow the Householder convention; |H| does not. *)
benchIf["HessenbergDecomposition 400x400", "HessenbergDecomposition",
  HessenbergDecomposition[rand01[{400, 400}]];, 1];
checkIf["HessenbergDecomposition 400x400", "HessenbergDecomposition",
  Round[10^6 Total[Flatten[Abs[Last[HessenbergDecomposition[s8]]]]]]];

(* Eigenvectors are unit-normalised in every column: {sum |lambda|, sum |v|^2}. *)
benchIf["Eigensystem 300x300 (general)", "Eigensystem", Eigensystem[g3];, 1];
checkIf["Eigensystem 300x300 (general)", "Eigensystem",
  With[{es = Eigensystem[g4]},
    {Round[10^6 Max[Abs[First[es]]]], Round[10^6 Total[Flatten[Abs[Last[es]]^2]]]}]];

benchIf["Eigensystem 300x300 complex", "Eigensystem", Eigensystem[c3];, 1];
checkIf["Eigensystem 300x300 complex", "Eigensystem",
  Round[10^6 Total[Abs[First[Eigensystem[c4]]]]]];

benchIf["SingularValueList 400x400", "SingularValueList",
  SingularValueList[rand01[{400, 400}]];];
checkIf["SingularValueList 400x400", "SingularValueList",
  Round[10^6 Total[SingularValueList[s8]]]];

benchIf["SingularValueList 300x300 complex", "SingularValueList",
  SingularValueList[c3];];
checkIf["SingularValueList 300x300 complex", "SingularValueList",
  Round[10^6 Total[SingularValueList[c4]]]];

(* 252^2 = 63504 determinants of 5 x 5 submatrices. *)
benchIf["Minors[m, 5] 10x10", "Minors", Minors[m10, 5];];
checkIf["Minors[m, 5] 10x10", "Minors", Round[10^3 Total[Flatten[Minors[s8, 4]]]]];

benchIf["Adjugate 200x200", "Adjugate", Adjugate[g2];];
checkIf["Adjugate 200x200", "Adjugate", Round[10^3 Total[Flatten[Adjugate[s8]]]]];

benchIf["Orthogonalize 200x200", "Orthogonalize", Orthogonalize[g2];];
checkIf["Orthogonalize 200x200", "Orthogonalize",
  Round[10^6 Total[Flatten[Orthogonalize[s8]]]]];

benchIf["MatrixExp 300x300", "MatrixExp", MatrixExp[g3/10];];
checkIf["MatrixExp 300x300", "MatrixExp", Round[10^6 Total[Flatten[MatrixExp[s8/10]]]]];

benchIf["MatrixExp 300x300 complex", "MatrixExp", MatrixExp[c3/10];];
checkIf["MatrixExp 300x300 complex", "MatrixExp",
  Round[10^6 Re[Total[Flatten[MatrixExp[c4/10]]]]]];

benchIf["MatrixLog 100x100", "MatrixLog", MatrixLog[p1];, 1];
checkIf["MatrixLog 100x100", "MatrixLog", Round[10^6 Total[Flatten[MatrixLog[s8]]]]];

benchIf["MatrixFunction[Sqrt] 100x100", "MatrixFunction", MatrixFunction[Sqrt, p1];, 1];
checkIf["MatrixFunction[Sqrt] 100x100", "MatrixFunction",
  Round[10^6 Total[Flatten[MatrixFunction[Sqrt, s8]]]]];

benchIf["MatrixFunction[Cos] 100x100", "MatrixFunction", MatrixFunction[Cos, r1];, 1];
checkIf["MatrixFunction[Cos] 100x100", "MatrixFunction",
  Round[10^6 Total[Flatten[MatrixFunction[Cos, s8/10]]]]];

benchIf["MatrixPower[m, 0.3] 100x100", "MatrixPower", MatrixPower[p1, 0.3];, 1];
checkIf["MatrixPower[m, 0.3] 100x100", "MatrixPower",
  {Round[10^6 Total[Flatten[MatrixPower[s8, 0.3]]]],
   Round[10^6 Total[Flatten[MatrixPower[g4, 0.3]]]]}];
