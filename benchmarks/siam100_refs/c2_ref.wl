(* c2 reference (method 2): independent implementation in WL arbitrary
   precision (precision 200, explicitly tracked): the next hit is found with
   NSolve on the quadratic |p + t v - c|^2 == 1/9 for every lattice point c in
   the box around the remaining segment; reflection about the normal. *)
dist[T_, prec_] := Module[{p = SetPrecision[{1/2, 1/10}, prec], v = SetPrecision[{1, 0}, prec],
    left = SetPrecision[T, prec], last = {}, hits, best, t, c, n, q, box},
  While[True,
    q = {p, p + left v};
    box = Tuples[{Range[Floor[Min[q[[All, 1]]]] - 1, Ceiling[Max[q[[All, 1]]]] + 1],
                  Range[Floor[Min[q[[All, 2]]]] - 1, Ceiling[Max[q[[All, 2]]]] + 1]}];
    hits = Flatten[Table[If[c === last, {},
        {#, c} & /@ Select[tt /. NSolve[(p + tt v - c).(p + tt v - c) == 1/9, tt, Reals, WorkingPrecision -> prec],
          # > 0 &]], {c, box}], 1];
    If[hits === {} || Min[hits[[All, 1]]] >= left, Return[Norm[p + left v]]];
    {t, c} = First[SortBy[hits, First]];
    p = p + t v; left = left - t; last = c; n = (p - c) 3;
    v = v - 2 (v.n) n]];
Print[N[dist[10, 200], 70]];
Print[N[dist[3, 200], 70]];
