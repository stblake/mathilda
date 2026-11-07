(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 184 -- SIAM 100-Digit Challenge, Appendix D Problem 9
   ("A triangular walk"; problem 11 on the book's web page).

     A particle starts at the top vertex of a triangular lattice with 30
     points on each side (465 points) and takes 60 random steps, each to a
     neighbour of its CURRENT point chosen uniformly (6 neighbours inside,
     4 on an edge, 2 at a corner).  What is the probability that it ends in
     the bottom row?        Answer: 9.5123435020743320973... e-6 (exact
     rational, numerator 58 digits / denominator 63 digits).

   The idiomatic statement is the 465 x 465 row-stochastic transition matrix P
   and the first row of P^60, summed over the bottom-row columns:
     machine : MatrixPower of the machine-real matrix;
     wp50    : exact rational arithmetic -- the start vector pushed through
               60 exact Dot products (v . P), then N[., 60].  (An exact
               MatrixPower of P costs ~6 dense rational 465^3 products; the
               60 vector products are what a user does for the exact answer.)
   Point (r, c), 0 <= c <= r <= 29, has index r(r+1)/2 + c + 1; its
   neighbours are (r, c+-1), (r-1, c-1), (r-1, c), (r+1, c), (r+1, c+1)
   that lie in the triangle.   True value and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"MatrixPower", "SparseArray"}];

nT = 30;
idx[{r_, c_}] := r (r + 1)/2 + c + 1;
nbrs[{r_, c_}] := Select[{{r, c - 1}, {r, c + 1}, {r - 1, c - 1}, {r - 1, c},
                          {r + 1, c}, {r + 1, c + 1}}, 0 <= #[[2]] <= #[[1]] < nT &];
pts = Flatten[Table[{r, c}, {r, 0, nT - 1}, {c, 0, r}], 1];
P = Normal[SparseArray[Flatten[Table[With[{l = nbrs[p]},
        ({idx[p], idx[#]} -> 1/Length[l]) & /@ l], {p, pts}]], {Length[pts], Length[pts]}]];
bottom = Table[idx[{nT - 1, c}], {c, 0, nT - 1}];
e1 = UnitVector[Length[pts], 1];

vrun["d9.machine", Total[MatrixPower[N[P], 60][[1, bottom]]], 10, 3];
vrun["d9.wp50", N[Total[Nest[#.P &, e1, 60][[bottom]]], 60], 48];
