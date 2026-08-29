(* Prototype of rt_cherry_dilog_exp in Mathilda's own language. *)

ClearAll[dilogExpTry];
dilogExpTry[f0_, x_] := Module[
  {f = f0, exps, rates, c, th, F, outer, rhoList, roots, uvars, urules,
   normLogRule, ekRule, norm, cands, gcandidates, basis, dbasis, syms,
   ci, resid, rnum, cl, eqs, unk, sol, ans, Q, chk, matchvars, Ftest},

  (* ---- 1. exponential rate c ---- *)
  exps = DeleteDuplicates[Cases[f, Power[E, e_] :> e, {0, Infinity}]];
  exps = Select[exps, ! FreeQ[#, x] &];
  If[exps === {}, Return[$Failed]];
  (* each exponent must be an exact rate*x (linear, zero constant term) *)
  If[! (And @@ (PolynomialQ[#, x] && Exponent[#, x] == 1 &&
       Coefficient[#, x, 0] === 0 & /@ exps)), Return[$Failed]];
  rates = Coefficient[#, x, 1] & /@ exps;
  If[! (And @@ (IntegerQ /@ rates)), Return[$Failed]];  (* proto: integer rates *)
  c = GCD @@ rates;                                       (* fundamental rate    *)
  th = rmT;

  (* ---- kernelize E^(k c x) -> th^k ---- *)
  ekRule = Power[E, e_] :> th^(Coefficient[e, x, 1]/c);
  F = f /. ekRule;

  (* ---- 2. outer logs Log[P(th)] and their roots ---- *)
  outer = DeleteDuplicates[Cases[F, Log[p_] :> p, {0, Infinity}]];
  outer = Select[outer, ! FreeQ[#, th] &];

  (* roots from denominator of F (in th) and from outer-log arguments *)
  roots = {};
  Module[{den = Denominator[Together[F /. Log[_] -> 0]]},
    roots = Join[roots,
      Cases[Solve[den == 0, th], {th -> r_} :> r]]];
  Do[roots = Join[roots, Cases[Solve[p == 0, th], {th -> r_} :> r]], {p, outer}];
  roots = DeleteDuplicates[Select[roots, (Element[#, Rationals] && # =!= 0) &]];
  rhoList = roots;

  (* ---- kernel logs u_k = Log[th - rho_k]; and Log[th] = c x ---- *)
  uvars = Table[Unique["uu"], {Length[rhoList]}];
  urules = Join[
    Table[Log[th - rhoList[[k]]] -> uvars[[k]], {k, Length[rhoList]}],
    {Log[th] -> c x}];

  normLogRule = Log[z_] :> PowerExpand[Log[Factor[Together[z]]]];
  norm[e_] := (((e /. ekRule) /. normLogRule) //. urules);
  (* reject a candidate derivative that, after norm, still carries a complex/Pi
     (reversed pair Log[-1]=I Pi) or an unmatched th/x-dependent Log. *)
  badLog[e_] := (! FreeQ[e, Complex]) || (! FreeQ[e, Pi]) ||
                (! FreeQ[e, Log[a_] /; (! FreeQ[a, th]) || (! FreeQ[a, x])]);

  (* gate: F, after logs->u, must be rational in th and total deg<=1 in {x,u} *)
  Ftest = norm[F];
  (* (soft gate for prototype; rely on diff-back) *)

  (* ---- 3. Moebius dilog-argument candidates ---- *)
  gcandidates = {};
  Do[(* (A) rho/th ; (B) th/rho *)
     AppendTo[gcandidates, rhoList[[i]]/th];
     AppendTo[gcandidates, th/rhoList[[i]]],
     {i, Length[rhoList]}];
  Do[If[i != j,
       AppendTo[gcandidates, (th - rhoList[[i]])/(rhoList[[j]] - rhoList[[i]])]; (* C *)
       AppendTo[gcandidates, (th - rhoList[[i]])/(th - rhoList[[j]])]],          (* D *)
     {i, Length[rhoList]}, {j, Length[rhoList]}];
  gcandidates = DeleteDuplicates[gcandidates];

  (* ---- ansatz basis (answer terms) ---- *)
  basis = {};
  Do[AppendTo[basis, PolyLog[2, g]], {g, gcandidates}];
  Do[AppendTo[basis, x Log[th - rhoList[[k]]]];
     AppendTo[basis, Log[th - rhoList[[k]]]], {k, Length[rhoList]}];
  Do[If[i <= j,
       AppendTo[basis, Log[th - rhoList[[i]]] Log[th - rhoList[[j]]]]],
     {i, Length[rhoList]}, {j, Length[rhoList]}];
  AppendTo[basis, x^2/2];
  AppendTo[basis, x];

  (* substitute th -> E^(c x) so D knows the derivative *)
  basis = basis /. th -> E^(c x);

  (* keep only admissible candidates (drop reversed pairs / unmatched logs) *)
  Module[{keep, db, bb},
    db = {}; bb = {};
    Do[Module[{d = norm[D[basis[[i]], x]]},
         If[! badLog[d], AppendTo[db, d]; AppendTo[bb, basis[[i]]]]],
       {i, Length[basis]}];
    basis = bb;
  ];
  syms = Table[Unique["cc"], {Length[basis]}];
  dbasis = Table[norm[D[syms[[i]] basis[[i]], x]], {i, Length[basis]}];

  matchvars = Join[{th, x}, uvars];
  resid = norm[F] - Total[dbasis];
  rnum = Numerator[Together[resid]];
  cl = DeleteCases[Flatten[CoefficientList[rnum, matchvars]], 0];
  eqs = Thread[cl == 0];
  sol = Solve[eqs, syms];
  If[! (MatchQ[sol, {{__Rule}}] || MatchQ[sol, {{}}]), Return[$Failed]];

  ans = Total[Table[syms[[i]] basis[[i]], {i, Length[basis]}]];
  Q = ans /. First[sol];
  Q = Q /. Table[s -> 0, {s, syms}];       (* pin free unknowns *)
  Q = Q /. th -> E^(c x);

  chk = Simplify[PowerExpand[(D[Q, x] - f) /. Log[a_] :> Log[Factor[Together[a]]]]];
  If[chk === 0, Q, $Failed]
];

test[f_] := Module[{r},
  r = dilogExpTry[f, x];
  Print["f = ", f];
  Print["  Q = ", r];
  If[r =!= $Failed,
    Print["  diff-back = ", Simplify[PowerExpand[(D[r, x] - f) /. Log[a_] :> Log[Factor[a]]]]]];
  Print["----"]];

test[x/(-1 + E^x)];
test[x/(1 + E^x)];
test[x E^x/(-1 + E^x)];
test[x/(-1 + E^(2 x))];
test[Log[1 + E^x]];
test[Log[1 - E^x]];
test[Log[1 + E^-x]];
test[2 x/(-1 + E^x) + 3 Log[1 + E^x]];
Print["=== declines (should be $Failed) ==="];
test[x Log[1 + E^x]];
test[1/(E^x - 1)];
