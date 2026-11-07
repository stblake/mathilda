(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 227 -- SIAM 100-Digit Challenge, sealed holdout B for Problem 2
   ("Chaotic Scattering"): ten more unit-speed photons among fixed mirrors --
   disks seen from outside and straight walls -- with specular reflection.

     bc2.1   square lattice of disks r = 1/3; start (1/2,1/10), v = (3,4)/5;
             |p(T)| at T = 6
     bc2.2   square lattice r = 1/4; start (1/2,1/2), v = (4,3)/5; |p(5)|
     bc2.3   TRIANGULAR lattice (centres (i + j/2, j sqrt3/2)), r = 1/3;
             start (1/2, 1/5), v = (1,0); |p(4)|
     bc2.4   Sinai billiard: walls of [-1,1]^2 + disk r = 1/2 at 0;
             start (3/4,0), v = (3,4)/5; |p(12)|
     bc2.5   box [0,2]x[0,1] + disks r = 1/4 at (1/2,1/2), (3/2,1/2);
             start (1,1/10), v = (3,4)/5; |p(6)|
     bc2.6   OPEN three-disk scatterer: r = 1 at (0,6/5), (+-3sqrt3/5,-3/5);
             start (0,0), v = (12,-5)/13; |p(20)| (after escape)
     bc2.7   Sinai billiard as bc2.4, start (-1/2,-3/4), v = (5,12)/13;
             the TIME of the 16th reflection
     bc2.8   box [0,3]x[0,2] only (integrable); start (1/3,1/7),
             v = (12,5)/13; |p(100)|
     bc2.9   square lattice r = 9/20 (gaps 1/10); start (1/2,1/2),
             v = (3,4)/5; |p(3)|
     bc2.10  Sinai billiard as bc2.4 for a longer time, |p(25)| (harder:
             a double keeps fewer digits -- see truth.tsv)

   Idiomatic statement (the same in every system, as in experiment 167):
   straight flight; the next hit from the system's polynomial root finder on
   the exact intersection equation |p + t v - c|^2 == r^2 (disk) or
   n.(p + t v) == c (wall) -- WL: NSolve over the Reals; reflection
   v' = v - 2 (v.n) n with the unit normal n.  Glue: the candidate filter --
   the mirror just hit is excluded; a disk must be ahead ((c - p).v > 0) and
   within r of the line (|det(v, c - p)| < r); a wall must not be parallel;
   lattice disks are taken from the box around the remaining segment (capped
   at length 50).  Rows: machine precision ("bc2.k machine") and
   WorkingPrecision 60 ("bc2.k wp50").  True values: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSolve", "Det"}];

(* {lattice, finite mirrors, p0, v0, T, number of reflections} (exact);
   lattice = {"sq"|"tri", r} or None; mirror = {"disk", cx, cy, r} | {"line", nx, ny, c} *)
box2[x0_, x1_, y0_, y1_] := {{"line", 1, 0, x1}, {"line", -1, 0, -x0},
                             {"line", 0, 1, y1}, {"line", 0, -1, -y0}};
sinai = Append[box2[-1, 1, -1, 1], {"disk", 0, 0, 1/2}];
scene[1] = {{"sq", 1/3}, {}, {1/2, 1/10}, {3/5, 4/5}, 6, Infinity};
scene[2] = {{"sq", 1/4}, {}, {1/2, 1/2}, {4/5, 3/5}, 5, Infinity};
scene[3] = {{"tri", 1/3}, {}, {1/2, 1/5}, {1, 0}, 4, Infinity};
scene[4] = {None, sinai, {3/4, 0}, {3/5, 4/5}, 12, Infinity};
scene[5] = {None, Join[box2[0, 2, 0, 1], {{"disk", 1/2, 1/2, 1/4}, {"disk", 3/2, 1/2, 1/4}}],
            {1, 1/10}, {3/5, 4/5}, 6, Infinity};
scene[6] = {None, {{"disk", 0, 6/5, 1}, {"disk", -3 Sqrt[3]/5, -3/5, 1}, {"disk", 3 Sqrt[3]/5, -3/5, 1}},
            {0, 0}, {12/13, -5/13}, 20, Infinity};
scene[7] = {None, sinai, {-1/2, -3/4}, {5/13, 12/13}, Infinity, 16};
scene[8] = {None, box2[0, 3, 0, 2], {1/3, 1/7}, {12/13, 5/13}, 100, Infinity};
scene[9] = {{"sq", 9/20}, {}, {1/2, 1/2}, {3/5, 4/5}, 3, Infinity};
scene[10] = {None, sinai, {3/4, 0}, {3/5, 4/5}, 25, Infinity};

(* lattice disks around the segment from p to q, each tagged with its id *)
latDisks[None, _, _, _] := {};
latDisks[{"sq", r_}, p_, q_, prec_] := Module[{xs = {p[[1]], q[[1]]}, ys = {p[[2]], q[[2]]}},
  Flatten[Table[{{"disk", N[i, prec], N[j, prec], N[r, prec]}, {i, j}},
    {i, Floor[Min[xs]] - 1, Floor[Max[xs]] + 2}, {j, Floor[Min[ys]] - 1, Floor[Max[ys]] + 2}], 1]];
latDisks[{"tri", r_}, p_, q_, prec_] := Module[{xs = {p[[1]], q[[1]]}, ys = {p[[2]], q[[2]]},
    hh = N[Sqrt[3]/2, prec]},
  Flatten[Table[{{"disk", N[i + j/2, prec], N[j Sqrt[3]/2, prec], N[r, prec]}, {i, j}},
    {j, Floor[Min[ys]/hh] - 1, Floor[Max[ys]/hh] + 2},
    {i, Floor[Min[xs] - j/2] - 1, Floor[Max[xs] - j/2] + 2}], 1]];

(* the positive real roots t of the hit equation of mirror m from p along v *)
hitTimes[{"disk", cx_, cy_, r_}, p_, v_, prec_] :=
  Select[tt /. NSolve[(p + tt v - {cx, cy}).(p + tt v - {cx, cy}) == r^2, tt, Reals,
                      WorkingPrecision -> prec], # > 0 &];
hitTimes[{"line", nx_, ny_, c_}, p_, v_, prec_] :=
  Select[tt /. NSolve[{nx, ny}.(p + tt v) == c, tt, Reals, WorkingPrecision -> prec], # > 0 &];
okQ[{"disk", cx_, cy_, r_}, p_, v_] := ({cx, cy} - p).v > 0 && Abs[Det[{v, {cx, cy} - p}]] < r;
okQ[{"line", nx_, ny_, c_}, p_, v_] := {nx, ny}.v != 0;
normal[{"disk", cx_, cy_, r_}, p_] := (p - {cx, cy})/r;
normal[{"line", nx_, ny_, c_}, p_] := {nx, ny};

photon[k_, prec_] := Module[{sc = scene[k], fin, p, v, left, nmax, last = None, clock = 0,
    nb = 0, done = False, res, cands, hits, t, m, n},
  fin = MapIndexed[{N[#1, prec], #2} &, sc[[2]]];
  p = N[sc[[3]], prec]; v = N[sc[[4]], prec];
  left = If[sc[[5]] === Infinity, Infinity, N[sc[[5]], prec]]; nmax = sc[[6]];
  clock = N[0, prec];
  While[!done,
    cands = Select[Join[fin, latDisks[sc[[1]], p, p + Min[left, 50] v, prec]],
                   #[[2]] =!= last && okQ[#[[1]], p, v] &];
    hits = Flatten[Function[cd, {#, cd} & /@ hitTimes[cd[[1]], p, v, prec]] /@ cands, 1];
    If[hits === {} || Min[hits[[All, 1]]] >= left,
      p = p + left v; res = Norm[p]; done = True,
      {t, m} = First[SortBy[hits, First]];
      p = p + t v; left = left - t; clock = clock + t; nb++; last = m[[2]];
      n = normal[m[[1]], p]; v = v - 2 (v.n) n;
      If[nb == nmax, res = clock; done = True]]];
  res];

Do[vrun["bc2." <> ToString[k] <> " machine", photon[k, MachinePrecision], 10, 3, 60];
   vrun["bc2." <> ToString[k] <> " wp50", photon[k, 60], 48, 1, 120],
  {k, 10}];
