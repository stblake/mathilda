(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 167 -- SIAM 100-Digit Challenge, Problem 2 ("Chaotic Scattering").

       A photon moving at speed 1 in the x-y plane starts at t = 0 at
       (1/2, 1/10) heading due east.  Around every integer lattice point
       (i, j) there is a circular mirror of radius 1/3.  How far from (0, 0)
       is the photon at t = 10?                         =  0.9952629194...

   Idiomatic statement (the same in every system): straight-line flight,
   the next mirror hit from the system's polynomial root finder on the exact
   circle-intersection equation |p + t v - c|^2 == 1/9 (WL: NSolve over the
   Reals), specular reflection v' = v - 2 (v.n) n.  The only glue is the
   candidate filter: lattice points c ahead of the photon ((c - p).v > 0)
   within distance 1/3 of its line (|det(v, c - p)| < 1/3), in the box around
   the remaining segment.

   The flight is chaotic: each reflection amplifies an input error, ~13
   digits are lost over t = 10 (measured: a 16-digit run keeps ~3, a
   60-digit run ~47).  So:
     c2.machine     t = 10 in doubles: a perfect double computation has ~3
                    correct digits (need 3 -- the finding is how close each
                    system's machine arithmetic comes to that);
     c2.machine-t3  t = 3 in doubles: the reduced instance (|p(3)| =
                    1.3753586978...), where a double keeps ~13 digits;
     c2.wp50        t = 10 at 60-digit working precision (need 45).
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSolve", "Det"}];

photon[T_, prec_] := Module[{p = N[{1/2, 1/10}, prec], v = N[{1, 0}, prec],
    left = N[T, prec], last = None, done = False, q, cands, hits, t, c, n},
  While[!done,
    q = {p, p + left v};
    cands = Select[Tuples[{Range[Floor[Min[q[[All, 1]]]] - 1, Ceiling[Max[q[[All, 1]]]] + 1],
                           Range[Floor[Min[q[[All, 2]]]] - 1, Ceiling[Max[q[[All, 2]]]] + 1]}],
      # =!= last && (# - p).v > 0 && Abs[Det[{v, # - p}]] < 1/3 &];
    hits = Flatten[Function[cc, {#, cc} & /@ Select[
        tt /. NSolve[(p + tt v - cc).(p + tt v - cc) == 1/9, tt, Reals,
                     WorkingPrecision -> prec], # > 0 &]] /@ cands, 1];
    If[hits === {} || Min[hits[[All, 1]]] >= left,
      p = p + left v; done = True,
      {t, c} = First[SortBy[hits, First]];
      p = p + t v; left = left - t; last = c; n = 3 (p - c);
      v = v - 2 (v.n) n]];
  Norm[p]];

vrun["c2.machine", photon[10, MachinePrecision], 3, 3];
vrun["c2.machine-t3", photon[3, MachinePrecision], 10, 3];
vrun["c2.wp50", photon[10, 60], 45];
