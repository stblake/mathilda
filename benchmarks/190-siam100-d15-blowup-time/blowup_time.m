(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 190 -- SIAM 100-Digit Challenge, Appendix D, Problem 15
   (Trefethen): blow-up time of a reaction-diffusion equation.

       u_t = Laplacian u + e^u on the square [0, 3]^2, u = 0 on the boundary,
       u(x, y, 0) = 0.  At what time t_inf does u blow up?
                                                     answer = 2.1197288...

   The square is supercritical (Bratu: a steady state exists only for side
   < Sqrt[6.808...] = 2.609), so the solution blows up, at the centre, at
   t_inf ~ 2.12, with u(centre) ~ -log(t_inf - t) and a profile that shrinks
   like Sqrt[(t_inf - t) |log(t_inf - t)|].  A fixed grid therefore stops
   resolving the solution shortly before t_inf, which is what limits every
   method-of-lines answer (see truth.tsv: only ~8 digits are established).

   Asked of each tier: the idiomatic statement is the PDE handed to NDSolve
   (method of lines, its default spatial discretisation) integrated past the
   blow-up; the blow-up time is where NDSolve stops ("step size is
   effectively zero; singularity suspected"), i.e. the end of the solution's
   domain.
     machine: default options (no grid or tolerance tuning).
     wp50:    WorkingPrecision -> 60.  No system can approach 50 digits for a
              2-D blow-up problem (and the reference has ~8), so this tier is
              N/A; it is written so that it scores need = 6 against the 8 known
              digits and so that its cost is recorded.
   True value and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve"}];

tinf[wp_] := (u /. First[NDSolve[{
    D[u[t, x, y], t] == D[u[t, x, y], x, x] + D[u[t, x, y], y, y] + Exp[u[t, x, y]],
    u[0, x, y] == 0, u[t, 0, y] == 0, u[t, 3, y] == 0, u[t, x, 0] == 0, u[t, x, 3] == 0},
    u, {t, 0, 3}, {x, 0, 3}, {y, 0, 3}, WorkingPrecision -> wp]])["Domain"][[1, 2]];

vrun["d15.machine", tinf[MachinePrecision], 6];
vrun["d15.wp50", tinf[60], 6];
