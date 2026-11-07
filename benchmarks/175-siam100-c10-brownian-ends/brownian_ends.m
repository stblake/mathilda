(* bench-timeout: 400 *)
(* ==========================================================================
   Experiment 175 -- SIAM 100-Digit Challenge, Problem 10 ("Hitting the Ends").

       A particle at the centre of a 10 x 1 rectangle does Brownian motion
       until it hits the boundary.  What is the probability that it hits one
       of the ends rather than one of the sides?           = 3.837587979...e-7

   The probability is u(0,0), u harmonic on [-5,5] x [-1/2,1/2], u = 1 on the
   ends x = +-5 and u = 0 on the sides y = +-1/2.  Two statements:
     c10.machine / c10.wp50  the separated (Fourier) solution of that problem
         evaluated at the centre, summed with NSum:
             u(0,0) = (4/Pi) Sum_k (-1)^k/(2k+1) sech(5 (2k+1) Pi)
         -- the direct statement every system can make.  (The book's closed
         form (2/Pi) arcsin(...) via singular moduli is NOT used as a case.)
     c10.machine pde  the boundary-value problem itself, NDSolveValue on the
         Rectangle with DirichletCondition, evaluated at the centre (default
         mesh).  Only WL has a PDE solver among the systems; the others SKIP.
         need = 3: the digits a default-mesh FEM solve reaches on a value
         whose relative size is 4e-7 of the boundary data.
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSum", "NDSolveValue", "DirichletCondition"}];

vrun["c10.machine", 4/Pi NSum[(-1)^k/(2 k + 1) Sech[5 (2 k + 1) Pi], {k, 0, Infinity}], 10, 3];
vrun["c10.wp50", 4/Pi NSum[(-1)^k/(2 k + 1) Sech[5 (2 k + 1) Pi], {k, 0, Infinity},
                           WorkingPrecision -> 60], 48];
vrun["c10.machine pde",
  NDSolveValue[{Laplacian[u[x, y], {x, y}] == 0,
      DirichletCondition[u[x, y] == 1, Abs[x] == 5],
      DirichletCondition[u[x, y] == 0, Abs[y] == 1/2]},
    u, Element[{x, y}, Rectangle[{-5, -1/2}, {5, 1/2}]]][0, 0], 3];
