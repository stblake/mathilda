(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 172 -- SIAM 100-Digit Challenge, Problem 7 ("Too Large to Invert").

       A is the 20000 x 20000 matrix with the primes 2, 3, 5, ..., 224737 on
       the diagonal and 1 at every a_ij with |i - j| = 1, 2, 4, ..., 16384.
       What is (A^-1)_11 ?                                  = 0.7250783462...

   Idiomatic statement only: A as a SparseArray (Band specifications), and
   (A^-1)_11 as the first component of LinearSolve[A, e1].  The matrix is
   built outside the timed expression (it is the problem's data); the timing
   is the solve alone.
     c7.machine         n = 20000, machine reals
     c7.wp50            n = 20000, 60-digit reals (expected to exceed the cap
                        in most systems: an arbitrary-precision sparse solve)
     c7.machine n=2000  the same construction at n = 2000 (primes up to 17389,
                        offsets 1..1024): a size every system can attempt
     c7.wp50 n=2000     ditto at 60 digits
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"SparseArray", "Band", "LinearSolve", "Prime"}];

mat[n_] := SparseArray[Join[{Band[{1, 1}] -> Prime[Range[n]]},
   Flatten[Table[{Band[{1, 1 + 2^k}] -> 1, Band[{1 + 2^k, 1}] -> 1},
                 {k, 0, Floor[Log2[n - 1]]}]]], {n, n}];

a20 = mat[20000]; e20 = UnitVector[20000, 1];
a2 = mat[2000];   e2 = UnitVector[2000, 1];
a20m = N[a20]; a20w = N[a20, 60];
a2m = N[a2];   a2w = N[a2, 60];

vrun["c7.machine", First[LinearSolve[a20m, e20]], 10];
vrun["c7.wp50", First[LinearSolve[a20w, e20]], 48];
vrun["c7.machine n=2000", First[LinearSolve[a2m, e2]], 10, 3];
vrun["c7.wp50 n=2000", First[LinearSolve[a2w, e2]], 48];
