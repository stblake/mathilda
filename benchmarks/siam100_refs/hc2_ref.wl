(* hc2 (holdout A, family C2) reference, method B: independent WL
   implementation at precision 150: every mirror in the box around the
   remaining segment, hit times from NSolve on |p + t v - c|^2 == r^2 over the
   Reals, reflection about the normal.  Method A: hc2_ref.py (mpmath). *)
sel["all"] = (True &); sel["even"] = (EvenQ[#[[1]] + #[[2]]] &); sel["oddcol"] = (OddQ[#[[1]]] &);
prob = <|
  1 -> {{1/2, 1/10}, {1, 0}, 1/3, 1, "all", "norm"},
  2 -> {{1/2, 1/10}, {1, 0}, 1/4, 4, "all", "norm"},
  3 -> {{1/2, 1/2}, {3/5, 4/5}, 1/3, 3, "all", "x"},
  4 -> {{0, 1/2}, {5/13, 12/13}, 2/5, 3, "all", "norm"},
  5 -> {{1/2, 1/10}, {1, 0}, 1/3, 2, "all", "norm"},
  6 -> {{1/2, 1/2}, {4/5, 3/5}, 1/2, 5, "even", "norm"},
  7 -> {{1/2, 1/3}, {8/17, 15/17}, 1/5, 8, "all", "y"},
  8 -> {{1/2, 1/7}, {-7/25, 24/25}, 2/5, 3, "all", "fromstart"},
  9 -> {{0, 1/2}, {12/13, 5/13}, 1/3, 6, "oddcol", "norm"},
  10 -> {{1/2, 1/10}, {1, 0}, 1/3, 5, "all", "norm"}|>;
fly[{p0_, v0_, r_, T_, s_, q_}, prec_] := Module[{p = SetPrecision[p0, prec],
    v = SetPrecision[v0, prec], left = SetPrecision[T, prec], last = {}, hits, t, c, n, qq, box},
  While[True,
    qq = {p, p + left v};
    box = Select[Tuples[{Range[Floor[Min[qq[[All, 1]]]] - 1, Ceiling[Max[qq[[All, 1]]]] + 1],
                         Range[Floor[Min[qq[[All, 2]]]] - 1, Ceiling[Max[qq[[All, 2]]]] + 1]}], sel[s]];
    hits = Flatten[Table[If[c === last, {},
        {#, c} & /@ Cases[tt /. NSolve[(p + tt v - c).(p + tt v - c) == r^2, tt, Reals,
                                       WorkingPrecision -> prec], w_ /; w > 0]], {c, box}], 1];
    If[hits === {} || Min[hits[[All, 1]]] >= left, p = p + left v; Break[]];
    {t, c} = First[SortBy[hits, First]];
    p = p + t v; left = left - t; last = c; n = (p - c)/r;
    v = v - 2 (v.n) n];
  Switch[q, "norm", Norm[p], "x", p[[1]], "y", p[[2]], "fromstart", Norm[p - p0]]];
Do[Print[k, " ", ToString[NumberForm[N[fly[prob[k], 150], 70], 70]]], {k, 10}];
