(* bench-timeout: 500 *)
(* ==========================================================================
   Experiment 173 -- SIAM 100-Digit Challenge, Problem 8 ("Hot Plate").

       The square plate [-1,1]^2 is at u = 0.  At t = 0 the side x = 1 is
       raised to u = 5 and held there, the other three stay at 0, and
       u_t = Laplacian u.  When does the centre reach u = 1?
                                                          = 0.4240113870...

   Two statements:
     c8.machine / c8.wp50  the separated (Fourier) solution at the centre,
         u(0,0,t) = (5/4) (1 - w(t)^2),
         w(t) = (4/Pi) Sum_k (-1)^k/(2k+1) Exp[-(2k+1)^2 Pi^2 t/4]
         (w = centre value of the 1-D problem on [-1,1], zero ends, initial 1;
         the 2-D solution with all four sides at 5 is 5 (1 - w^2), and one
         side is a quarter of it by symmetry), summed with NSum inside a
         FindRoot for u(0,0,t) == 1.  The book's theta-function root
         (Appendix B) is NOT used as a case.
     c8.machine pde  the initial-boundary-value problem itself: NDSolveValue
         (method of lines) on the square, then FindRoot on the interpolated
         centre value.  The boundary data (5 on x = 1) is inconsistent with
         the initial data (0) -- that is the problem as posed.  need = 4: what
         a default method-of-lines grid can deliver on a discontinuous corner.
         Only WL has a PDE solver among the systems; the others SKIP.
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSum", "FindRoot", "NDSolveValue"}];

w[t_?NumericQ] := 4/Pi NSum[(-1)^k/(2 k + 1) Exp[-(2 k + 1)^2 Pi^2 t/4], {k, 0, Infinity}];
w60[t_?NumericQ] := 4/Pi NSum[(-1)^k/(2 k + 1) Exp[-(2 k + 1)^2 Pi^2 t/4], {k, 0, Infinity},
                              WorkingPrecision -> 60];

vrun["c8.machine", t /. FindRoot[5/4 (1 - w[t]^2) == 1, {t, 0.4}], 10, 3];
vrun["c8.wp50", t /. FindRoot[5/4 (1 - w60[t]^2) == 1, {t, 4/10}, WorkingPrecision -> 60], 48];
vrun["c8.machine pde",
  Module[{s},
    s = NDSolveValue[{D[u[t, x, y], t] == D[u[t, x, y], x, x] + D[u[t, x, y], y, y],
          u[0, x, y] == 0, u[t, -1, y] == 0, u[t, 1, y] == 5,
          u[t, x, -1] == 0, u[t, x, 1] == 0},
        u, {t, 0, 1}, {x, -1, 1}, {y, -1, 1}];
    t /. FindRoot[s[t, 0, 0] == 1, {t, 0.4}]], 4];
