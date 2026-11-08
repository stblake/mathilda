(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 232 -- SIAM 100-Digit Challenge, HOLDOUT B for Problem 7
   ("Too Large to Invert").  Ten new problems that need the same capability:
   build a large SPARSE matrix and read one component of the solution of a
   linear system with it, here (A^-1)_jj = (LinearSolve[A, e_j])_j unless
   stated otherwise.  p_i is the i-th prime.

     bc7.1  n = 1000, tridiagonal: p_i on the diagonal, 1 on the first sub-
            and super-diagonal.  (A^-1)_11.                         (easy)
     bc7.2  n = 10^4: the 5-point Laplacian on a 100 x 100 grid (4 on the
            diagonal, -1 to each grid neighbour, Dirichlet) plus I.  The
            diagonal entry of A^-1 at grid node (50,50), index 49*100+50.
     bc7.3  n = 20000: p_i on the diagonal, 1 at every |i-j| a Fibonacci
            number (1, 2, 3, 5, ..., 17711).  (A^-1)_11.
     bc7.4  n = 2000: p_i on the diagonal, 1 at every |i-j| a perfect square
            (1, 4, ..., 1936).  (A^-1)_11.
     bc7.5  n = 2000, NONSYMMETRIC: p_i on the diagonal, +1 at j - i = 2^k,
            -1 at i - j = 2^k (k = 0..10).  (A^-1)_11.
     bc7.6  n = 5000, the original matrix (1 at |i-j| = 2^k): the FIRST
            component of the solution of A x = (1, 1, ..., 1).
     bc7.7  n = p = 2003, an unstructured sparse matrix: vertices 0..p-1
            (index v+1), 7 on the diagonal, 1 joining v to 2v, 3v and 5v
            mod p (v = 1..p-1) symmetrically.  The diagonal entry of A^-1
            at vertex 1 (index 2).
     bc7.8  n = 27000: the 7-point Laplacian on a 30 x 30 x 30 grid (6 on the
            diagonal, Dirichlet).  (A^-1)_11 (a corner node).
     bc7.9  n = 22500, INDEFINITE: the 5-point Laplacian on a 150 x 150 grid
            minus I/2 (a shifted Helmholtz operator).  The diagonal entry of
            A^-1 at grid node (75,75), index 74*150+75.
     bc7.10 n = 2047: graph Laplacian of the complete binary tree with 11
            levels (node i has children 2i, 2i+1) plus I.  (A^-1)_11 (root).

   Every matrix is a SparseArray built outside the timed expression (it is
   the problem's data), converted to machine reals (N) or 60-digit reals
   (N[., 60]); the timing is the LinearSolve alone.
     "<id> machine"  LinearSolve on the machine-real SparseArray
     "<id> wp50"     LinearSolve on the 60-digit SparseArray (an arbitrary-
                     precision sparse solve; expected to exceed the cap in
                     most systems for the large n)
   True values (two independent methods each): truth.tsv,
   siam100_refs/bc7_ref.py.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"SparseArray", "Band", "LinearSolve", "Prime", "KroneckerProduct"}];

primeBands[n_, offs_, sup_, sub_] := SparseArray[Join[{Band[{1, 1}] -> Prime[Range[n]]},
   Flatten[Table[{Band[{1, 1 + o}] -> sup, Band[{1 + o, 1}] -> sub}, {o, offs}]]], {n, n}];
lap1[m_] := SparseArray[{Band[{1, 1}] -> 2, Band[{1, 2}] -> -1, Band[{2, 1}] -> -1}, {m, m}];
eye[m_] := SparseArray[Band[{1, 1}] -> 1, {m, m}];
lap2[m_] := KroneckerProduct[lap1[m], eye[m]] + KroneckerProduct[eye[m], lap1[m]];
lap3[m_] := KroneckerProduct[lap1[m], eye[m^2]] + KroneckerProduct[eye[m], lap2[m]];

a[1] = primeBands[1000, {1}, 1, 1];
a[2] = lap2[100] + eye[10^4];
a[3] = primeBands[20000, Select[Union[Fibonacci[Range[2, 30]]], # < 20000 &], 1, 1];
a[4] = primeBands[2000, Range[44]^2, 1, 1];
a[5] = primeBands[2000, 2^Range[0, 10], 1, -1];
a[6] = primeBands[5000, 2^Range[0, 12], 1, 1];
a[7] = SparseArray[Flatten[Table[{{v + 1, Mod[g v, 2003] + 1} -> 1, {Mod[g v, 2003] + 1, v + 1} -> 1},
          {v, 1, 2002}, {g, {2, 3, 5}}], 2], {2003, 2003}] + 7 eye[2003];
a[8] = lap3[30];
a[9] = lap2[150] - eye[22500]/2;
a[10] = With[{n = 2047}, SparseArray[Join[
     Table[{i, i} -> Which[i == 1, 3, 2 i <= n, 4, True, 2], {i, n}],
     Flatten[Table[{{i, c} -> -1, {c, i} -> -1}, {i, 1, 1023}, {c, {2 i, 2 i + 1}}], 2]], {n, n}]];
idx = <|1 -> 1, 2 -> 4950, 3 -> 1, 4 -> 1, 5 -> 1, 6 -> 1, 7 -> 2, 8 -> 1, 9 -> 11175, 10 -> 1|>;
rhs[k_, n_] := If[k == 6, ConstantArray[1, n], UnitVector[n, idx[k]]];

Do[With[{n = Length[a[k]]},
    am[k] = N[a[k]]; aw[k] = N[a[k], 60];
    bm[k] = N[rhs[k, n]]; bw[k] = N[rhs[k, n], 60]], {k, 10}];

Do[With[{k = k, j = idx[k]},
    vrun["bc7." <> ToString[k] <> " machine", LinearSolve[am[k], bm[k]][[j]], 10, 1, 90];
    vrun["bc7." <> ToString[k] <> " wp50", LinearSolve[aw[k], bw[k]][[j]], 48, 1, 75]],
  {k, 10}];
