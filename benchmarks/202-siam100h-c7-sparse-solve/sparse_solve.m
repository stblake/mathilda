(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 202 -- SIAM 100-Digit Challenge, HOLDOUT A for Problem 7
   ("Too Large to Invert"): ten new entries / functionals of the inverse of a
   large sparse integer matrix.  As in experiment 172 the matrix is the
   problem's data -- a SparseArray of Band specifications, built (and N[]-ed)
   outside the timed expression -- and the case is LinearSolve alone.
   P(n) = the first n primes; "+-o" = 1 at every a_ij with |i - j| = o.

     hc7.1  n = 20000, diag P, +-1 (tridiagonal):          (A^-1)_11.
     hc7.2  n = 10000, diag P, +-1, +-2 (pentadiagonal):    (A^-1)_11.
     hc7.3  n = 5000,  diag P, +-3^k (1, 3, 9, ..., 2187):  (A^-1)_11.
     hc7.4  n = 3000,  diag P, +-Fibonacci (1,2,3,5,...,2584): (A^-1)_11
            (15 band pairs; not diagonally dominant near the top).
     hc7.5  n = 1000,  diag 1, 2, ..., n (NOT primes), +-2^k: (A^-1)_11
            (an indefinite matrix).
     hc7.6  n = 3000,  diag P, +-2^k (C7's matrix, smaller): (A^-1)_nn.
     hc7.7  n = 1000,  diag P, +-2^k: the off-diagonal entry (A^-1)_12.
     hc7.8  5-point Laplacian (4 on the diagonal, -1 to grid neighbours,
            Dirichlet) on a 60 x 60 grid, n = 3600: (A^-1)_11 (corner node).
     hc7.9  n = 1000,  diag P, +-2^k: 1^T A^-1 1 (Total of the solution of
            A x = (1,...,1)).
     hc7.10 7-point Laplacian on a 20 x 20 x 20 grid, n = 8000: (A^-1)_11
            (3-D fill-in: the harder one).
   Rows: "<id> machine" = machine reals; "<id> wp50" = N[A, 60].  Every case
   is capped at 60 s, and the machine rows run first: Mathilda's
   arbitrary-precision LinearSolve densifies the SparseArray and its parallel
   MPFR LU did not honour the TimeConstrained cap in a full-file run (v0.6377:
   hc7.2 wp50, run right after the machine rows of hc7.1/hc7.2, overran its
   60 s cap by > 12 min; in a fresh process the same case stops at 60 s), so a
   wp50 row can stall the rest of the file.  True values: truth.tsv (siam100_refs/hc7_ref.py:
   mixed-precision iterative refinement with exact integer residuals vs
   wolframscript LinearSolve at WorkingPrecision 80 / eigen-sum closed forms
   for the Laplacians / the continued fraction for the tridiagonal one).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"SparseArray", "Band", "LinearSolve", "Prime"}];

mat[n_, offs_, d_] := SparseArray[Join[{Band[{1, 1}] -> d},
   Flatten[Table[{Band[{1, 1 + o}] -> 1, Band[{1 + o, 1}] -> 1}, {o, offs}]]], {n, n}];
pow2[n_] := 2^Range[0, Floor[Log2[n - 1]]];
(* grid Laplacians: -1 to each neighbour inside the grid *)
lap2[m_] := With[{n = m^2, r = Table[If[Mod[i, m] == 0, 0, -1], {i, m^2 - 1}]},
   SparseArray[{Band[{1, 1}] -> 4, Band[{1, 2}] -> r, Band[{2, 1}] -> r,
                Band[{1, 1 + m}] -> -1, Band[{1 + m, 1}] -> -1}, {n, n}]];
lap3[m_] := With[{n = m^3,
    r1 = Table[If[Mod[i, m] == 0, 0, -1], {i, m^3 - 1}],
    rm = Table[If[Mod[i, m^2] > m^2 - m || Mod[i, m^2] == 0, 0, -1], {i, m^3 - m}]},
   SparseArray[{Band[{1, 1}] -> 6, Band[{1, 2}] -> r1, Band[{2, 1}] -> r1,
                Band[{1, 1 + m}] -> rm, Band[{1 + m, 1}] -> rm,
                Band[{1, 1 + m^2}] -> -1, Band[{1 + m^2, 1}] -> -1}, {n, n}]];

fib[n_] := Select[Union[Fibonacci[Range[2, 30]]], # < n &];
mats = {
  mat[20000, {1}, Prime[Range[20000]]],
  mat[10000, {1, 2}, Prime[Range[10000]]],
  mat[5000, 3^Range[0, 7], Prime[Range[5000]]],
  mat[3000, fib[3000], Prime[Range[3000]]],
  mat[1000, pow2[1000], Range[1000]],
  mat[3000, pow2[3000], Prime[Range[3000]]],
  mat[1000, pow2[1000], Prime[Range[1000]]],
  lap2[60],
  mat[1000, pow2[1000], Prime[Range[1000]]],
  lap3[20]};
Am = N /@ mats; Aw = N[#, 60] & /@ mats;
dim[i_] := Length[mats[[i]]];
e1[i_] := UnitVector[dim[i], 1];

vrun["hc7.1 machine", First[LinearSolve[Am[[1]], e1[1]]], 10, 3, 60];
vrun["hc7.2 machine", First[LinearSolve[Am[[2]], e1[2]]], 10, 3, 60];
vrun["hc7.3 machine", First[LinearSolve[Am[[3]], e1[3]]], 10, 1, 60];
vrun["hc7.4 machine", First[LinearSolve[Am[[4]], e1[4]]], 10, 1, 60];
vrun["hc7.5 machine", First[LinearSolve[Am[[5]], e1[5]]], 10, 3, 60];
vrun["hc7.6 machine", Last[LinearSolve[Am[[6]], UnitVector[3000, 3000]]], 10, 1, 60];
vrun["hc7.7 machine", LinearSolve[Am[[7]], e1[7]][[2]], 10, 3, 60];
vrun["hc7.8 machine", First[LinearSolve[Am[[8]], e1[8]]], 10, 3, 60];
vrun["hc7.9 machine", Total[LinearSolve[Am[[9]], ConstantArray[1, 1000]]], 10, 3, 60];
vrun["hc7.10 machine", First[LinearSolve[Am[[10]], e1[10]]], 10, 1, 60];
vrun["hc7.1 wp50", First[LinearSolve[Aw[[1]], e1[1]]], 48, 1, 60];
vrun["hc7.2 wp50", First[LinearSolve[Aw[[2]], e1[2]]], 48, 1, 60];
vrun["hc7.3 wp50", First[LinearSolve[Aw[[3]], e1[3]]], 48, 1, 60];
vrun["hc7.4 wp50", First[LinearSolve[Aw[[4]], e1[4]]], 48, 1, 60];
vrun["hc7.5 wp50", First[LinearSolve[Aw[[5]], e1[5]]], 48, 1, 60];
vrun["hc7.6 wp50", Last[LinearSolve[Aw[[6]], UnitVector[3000, 3000]]], 48, 1, 60];
vrun["hc7.7 wp50", LinearSolve[Aw[[7]], e1[7]][[2]], 48, 1, 60];
vrun["hc7.8 wp50", First[LinearSolve[Aw[[8]], e1[8]]], 48, 1, 60];
vrun["hc7.9 wp50", Total[LinearSolve[Aw[[9]], ConstantArray[1, 1000]]], 48, 1, 60];
vrun["hc7.10 wp50", First[LinearSolve[Aw[[10]], e1[10]]], 48, 1, 60];
