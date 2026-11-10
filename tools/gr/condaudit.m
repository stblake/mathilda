(* condaudit.m -- audit the book's CONDITIONS of the definite table records.

   Driven by tools/gr/condaudit.py.  The oracle (verify.m) checked each G&R
   formula at two parameter instantiations; the table (src/internal/GRTable.m)
   serves a symbolic query whenever the caller's Assumptions imply the book's
   conditions, without a numeric check (there are symbols left).  A condition
   the book states too weakly therefore becomes a wrong answer: G&R 3.197.7
   says Re mu > 0 where the integral converges only for mu > 1/2.

   For every definite record this samples further parameter points that
   satisfy the book's conditions (exact rationals, deterministic per record),
   evaluates the book's value there and compares it with NIntegrate.

   A record that disagrees somewhere is wrong as stated on part of its stated
   domain.  One correction is derived from first principles and applied
   automatically: when the integrand and the limits are EVEN in a parameter p
   (unchanged by p -> -p), the integral is even in p, so a value verified for
   p > 0 is the integral for p > 0 and the formula's domain is p > 0 (the book
   writes a^2 > b^2 where it means a > 0; G&R 3.165.2's value is odd in a).
   The record is re-audited with p > 0 added for one such parameter at a time,
   then for all of them; the first set under which it agrees everywhere is the
   correction, and without one the record is left out.

     @@A  id  agree  disagree  nint-failed  value-nonfinite  points  first-disagreement  fix  fix-agree  fix-disagree *)

grAuditVals = {1/7, 1/3, 3/5, 5/6, 1, 7/5, 9/4, 3, 11/2, -1/4, -2/3, -3/2, -7/3};

grAuditSeed[id_String] := Mod[Total[ToCharacterCode[id] Range[StringLength[id]]], 2^31];

(* {agree, disagree, nint-failed, value-nonfinite, points, first disagreement} *)
grAuditRun[id_String, f_, x_, lim_, v_, ps_List, conds_List] :=
  Module[{pts = {}, tries = 0, pt, ok, nv, ni, msg, agree = 0, dis = 0, nfail = 0, vbad = 0, first = ""},
    SeedRandom[grAuditSeed[id]];
    While[Length[pts] < 4 && tries < 300,
      tries++;
      pt = Thread[ps -> RandomChoice[grAuditVals, Length[ps]]];
      If[MemberQ[pts, pt], Continue[]];
      ok = Quiet[TimeConstrained[And @@ (TrueQ[Refine[#]] & /@ (conds /. pt)), 2, False]];
      If[ok, AppendTo[pts, pt]]];
    Do[
      nv = Quiet[TimeConstrained[N[v /. pt, 20], 10, $Failed]];
      If[!NumberQ[nv] || Abs[nv] > 10^200, vbad++; Continue[]];
      (* a quadrature that did not converge is no verdict -- unless its
         estimate is 10^6 times the book's finite value: the integral
         diverges at a point the conditions admit (G&R 3.197.7 at mu = 1/7) *)
      msg = True; ni = $Failed;
      Quiet[TimeConstrained[msg = Check[ni = NIntegrate[f /. pt, {x, lim[[1]] /. pt, lim[[2]] /. pt},
                                                        MaxRecursion -> 20]; False, True], 20, Null]];
      If[!NumberQ[ni], nfail++; Continue[]];
      If[msg && Abs[ni] <= 10^6 Max[1, Abs[nv]], nfail++; Continue[]];
      If[Abs[ni - nv] <= 10^-5 Max[1, Abs[nv]], agree++,
        dis++; If[first === "", first = ToString[{pt, N[nv, 10], ni}, InputForm]]],
      {pt, pts}];
    {agree, dis, nfail, vbad, Length[pts], first}];

(* the parameters the integrand and the limits are even in *)
grEvenParams[f_, lim_, ps_List] := Select[ps,
  FreeQ[lim, #] && TrueQ[Quiet[PossibleZeroQ[(f /. # -> -#) - f]]] &];

grAuditOne[id_String, f_, x_, lim_, v_, ps_List, conds_List] :=
  Module[{r, ev, fix = "", rf = {0, 0}},
    r = grAuditRun[id, f, x, lim, v, ps, conds];
    If[r[[2]] > 0,
      ev = grEvenParams[f, lim, ps];
      (* the smallest correction: one even parameter at a time, then all *)
      Do[If[fix =!= "", Break[]];
         rf = grAuditRun[id <> "+fix", f, x, lim, v, ps, Join[conds, (# > 0) & /@ set]][[{1, 2}]];
         If[rf[[1]] > 0 && rf[[2]] == 0,
           fix = StringRiffle[ToString[# > 0, InputForm] & /@ set, " && "]],
         {set, If[Length[ev] > 1, Append[List /@ ev, ev], {ev}]}];
      If[fix === "", rf = {0, 0}]];
    Print["@@A\t", id, "\t", StringRiffle[ToString /@ r[[1 ;; 5]], "\t"], "\t", r[[6]], "\t",
          fix, "\t", rf[[1]], "\t", rf[[2]]]];

GRTable`Private`add[key_, id_String, f_, x_, lim_, Hold[v_], ps_List, Hold[conds_List]] :=
  If[lim =!= None && If[gAuditIds === All, Mod[grAuditSeed[id], gAuditN] === gAuditK,
                         MemberQ[gAuditIds, id]],
    grAuditOne[id, f, x, lim, v, ps, conds]];
