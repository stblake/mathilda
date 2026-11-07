(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 178 -- SIAM 100-Digit Challenge, Appendix D Problem 3
   (Trefethen): the Thomson problem for N = 100.

       E_N = min over N points x_j on the unit sphere of
             Sum_{j<k} |x_j - x_k|^-1,       E_100 = 4448.3506343313...

   Idiomatic statement: FindMinimum of the energy over spherical coordinates
   from a deterministic start, the Fibonacci (golden-angle) lattice
       theta_k = ArcCos[1 - 2k/N],  phi_k = k Pi (3 - Sqrt[5]),  k = 0..N-1.
   Point 0 is then the north pole and is held there; phi_1 is held at its
   start value (the two fix the rotation group), leaving 2N - 3 = 197
   unknowns.  From this start a plain quasi-Newton descent reaches the global
   minimum (scipy BFGS does; see siam100_refs/d3_ref.py: random starts reach
   the global minimum only about 40% of the time, and the classical Rakhmanov-
   Saff-Zhou spiral start lands in the 4448.5069 local minimum).  A system
   that lands in a different local minimum has given a WRONG answer.
   Rows: machine precision; WorkingPrecision -> 60 (the minimum VALUE is
   accurate to about twice the digits of the minimiser, so the default goals
   at WP 60 suffice for 50 digits of E).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindMinimum"}];

npts = 100;
ph1 = Pi (3 - Sqrt[5]);
pt[0] = {0, 0, 1};
pt[1] = {Sin[th[1]] Cos[ph1], Sin[th[1]] Sin[ph1], Cos[th[1]]};
pt[k_Integer] := {Sin[th[k]] Cos[ph[k]], Sin[th[k]] Sin[ph[k]], Cos[th[k]]};
energy = Sum[1/Sqrt[(pt[j] - pt[k]).(pt[j] - pt[k])], {j, 0, npts - 2}, {k, j + 1, npts - 1}];
start = Join[Table[{th[k], N[ArcCos[1 - 2 k/npts]]}, {k, 1, npts - 1}],
             Table[{ph[k], N[k ph1]}, {k, 2, npts - 1}]];
start60 = Join[Table[{th[k], N[ArcCos[1 - 2 k/npts], 60]}, {k, 1, npts - 1}],
               Table[{ph[k], N[k ph1, 60]}, {k, 2, npts - 1}]];

vrun["d3.machine", First[FindMinimum[energy, start]], 10];
vrun["d3.wp50", First[FindMinimum[energy, start60, WorkingPrecision -> 60]], 48, 1, 600];
