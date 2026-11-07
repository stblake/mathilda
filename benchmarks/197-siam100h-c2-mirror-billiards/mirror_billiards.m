(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 197 -- SIAM 100-Digit Challenge, holdout A for Problem 2
   (family hc2: a photon reflecting off circular mirrors centred on lattice
   points).  Speed 1, start p0, unit heading v0, mirrors of radius r on the
   lattice points (i, j) selected by `sel`; report a coordinate of the photon
   at time T:

     k   p0          v0           r     T  mirrors          quantity
     1   (1/2,1/10)  (1,0)        1/3   1  all              |p(T)|          (easy)
     2   (1/2,1/10)  (1,0)        1/4   4  all              |p(T)|
     3   (1/2,1/2)   (3/5,4/5)    1/3   3  all              x(T)
     4   (0,1/2)     (5/13,12/13) 2/5   3  all              |p(T)|
     5   (1/2,1/10)  (1,0)        1/3   2  all              |p(T)|          (easy)
     6   (1/2,1/2)   (4/5,3/5)    1/2   5  i+j even         |p(T)|          (checkerboard)
     7   (1/2,1/3)   (8/17,15/17) 1/5   8  all              y(T)            (harder)
     8   (1/2,1/7)   (-7/25,24/25) 2/5  3  all              |p(T) - p0|
     9   (0,1/2)     (12/13,5/13) 1/3   6  i odd            |p(T)|          (mirror columns)
     10  (1/2,1/10)  (1,0)        1/3   5  all              |p(T)|          (C2 at T = 5)

   The statement is C2's (experiment 167), unchanged: straight flight, next
   hit from NSolve over the Reals on |p + t v - c|^2 == r^2, reflection
   v' = v - 2 (v.n) n, over the candidate filter (mirrors ahead of the photon
   within distance r of its line, in the box around the remaining segment).
   Chaotic loss of digits, measured (siam100_refs/hc2_ref.py: the digits a
   16-digit / 60-digit run keeps): 16.8/62, 12.3/58, 12.5/57.8, 12.8/56.5,
   15.3/59.9, 12.1/56.1, 6.7/50.8, 14.2/58.7, 11.6/56.4, 11.0/55.4.  So every
   machine row needs 10 except hc2.7 (need 6: a perfect double computation
   keeps 6.7 digits); every wp50 row needs 48.
   True values: truth.tsv; references siam100_refs/hc2_ref.{py,wl}.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSolve", "Det"}];

sel["all"] = (True &); sel["even"] = (EvenQ[#[[1]] + #[[2]]] &); sel["oddcol"] = (OddQ[#[[1]]] &);

photon[{p0_, v0_, r_, T_, s_, qty_}, prec_] := Module[{p = N[p0, prec], v = N[v0, prec],
    left = N[T, prec], last = None, done = False, q, cands, hits, t, c, n},
  While[!done,
    q = {p, p + left v};
    cands = Select[Tuples[{Range[Floor[Min[q[[All, 1]]]] - 1, Ceiling[Max[q[[All, 1]]]] + 1],
                           Range[Floor[Min[q[[All, 2]]]] - 1, Ceiling[Max[q[[All, 2]]]] + 1]}],
      sel[s][#] && # =!= last && (# - p).v > 0 && Abs[Det[{v, # - p}]] < r &];
    hits = Flatten[Function[cc, {#, cc} & /@ Select[
        tt /. NSolve[(p + tt v - cc).(p + tt v - cc) == r^2, tt, Reals,
                     WorkingPrecision -> prec], # > 0 &]] /@ cands, 1];
    If[hits === {} || Min[hits[[All, 1]]] >= left,
      p = p + left v; done = True,
      {t, c} = First[SortBy[hits, First]];
      p = p + t v; left = left - t; last = c; n = (p - c)/r;
      v = v - 2 (v.n) n]];
  Switch[qty, "norm", Norm[p], "x", p[[1]], "y", p[[2]], "fromstart", Norm[p - N[p0, prec]]]];

hc2 = {
  {{1/2, 1/10}, {1, 0}, 1/3, 1, "all", "norm"},
  {{1/2, 1/10}, {1, 0}, 1/4, 4, "all", "norm"},
  {{1/2, 1/2}, {3/5, 4/5}, 1/3, 3, "all", "x"},
  {{0, 1/2}, {5/13, 12/13}, 2/5, 3, "all", "norm"},
  {{1/2, 1/10}, {1, 0}, 1/3, 2, "all", "norm"},
  {{1/2, 1/2}, {4/5, 3/5}, 1/2, 5, "even", "norm"},
  {{1/2, 1/3}, {8/17, 15/17}, 1/5, 8, "all", "y"},
  {{1/2, 1/7}, {-7/25, 24/25}, 2/5, 3, "all", "fromstart"},
  {{0, 1/2}, {12/13, 5/13}, 1/3, 6, "oddcol", "norm"},
  {{1/2, 1/10}, {1, 0}, 1/3, 5, "all", "norm"}};
needM = {10, 10, 10, 10, 10, 10, 6, 10, 10, 10};

Do[
  With[{pr = hc2[[k]], nd = needM[[k]]},
    vrun["hc2." <> ToString[k] <> " machine", photon[pr, MachinePrecision], nd, 1, 60];
    vrun["hc2." <> ToString[k] <> " wp50", photon[pr, 60], 48, 1, 60]],
  {k, 10}];
