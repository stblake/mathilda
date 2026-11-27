(* ========================================================================== *)
(* goursat/Hermite.m -- Hermite reduction on y^n = R and the kind of the       *)
(* reduced differential (paper Proposition "Hermite reduction"; goursat.py:    *)
(* hermite, classify_reduced).                                                 *)
(* ========================================================================== *)

(* F dt/y^k = d(A y^(n-k)) + (P + N/Ds) dt/y^k  with deg P <= deg R - 2, Ds     *)
(* squarefree and coprime to R, deg N < deg Ds.  The reduced form is unique, so *)
(* it vanishes iff F dt/y^k is exact.  Returns {A, P, N, Ds}.  One linear solve: *)
(* the poles of A are read off the factorisation of den(F).                     *)
hermite[F0_, R_, t_, n_, k_] := Module[
  {F, e, d, Fn, Fd, DA = 1, Ds = 1, fl, f, mult, degF, dPA, dN, PA, P, Nn, num,
   ca, cp, cn, unknowns, sol},
  F = gK[F0];
  If[F === 0, Return[{0, 0, 0, 1}, Module]];
  e = (n - k)/n;
  d = Exponent[R, t];
  {Fn, Fd} = numDen[F];
  fl = FactorList[Fd, Extension -> Automatic];
  Do[{f, mult} = fac;
     If[!FreeQ[f, t],
       If[dividesQ[f, R, t],
         DA = DA f^mult,                       (* a pole at a branch point: removed entirely *)
         DA = DA f^(mult - 1); Ds = Ds f]],
     {fac, fl}];
  DA = Expand[DA]; Ds = Expand[Ds];
  degF = Exponent[Fn, t] - Exponent[Fd, t];
  dPA = Exponent[DA, t] + Max[degF - d + 1, -1];
  dN = Exponent[Ds, t];
  PA = Sum[ca[i] t^i, {i, 0, dPA}];
  P = Sum[cp[i] t^i, {i, 0, d - 2}];
  Nn = Sum[cn[i] t^i, {i, 0, dN - 1}];
  (* (A' R + e A R' + P + N/Ds - F) over the common denominator DA^2 Ds Fd, A = PA/DA *)
  num = ((D[PA, t] DA - PA D[DA, t]) R + e PA DA D[R, t]) Ds Fd +
        (P Ds + Nn) DA^2 Fd - Fn DA^2 Ds;
  unknowns = Join[Table[ca[i], {i, 0, dPA}], Table[cp[i], {i, 0, d - 2}],
                  Table[cn[i], {i, 0, dN - 1}]];
  sol = solveIdentity[num, unknowns, t];
  {gK[(PA /. sol)/DA], Expand[gK[P /. sol]], Expand[gK[Nn /. sol]], Ds}];

(* Kind of the Hermite-reduced differential (P + N/Ds) dt/y^k:                  *)
(* "zero", "second-kind" (not exact, no residues: a proof of non-elementarity),  *)
(* or "third-kind" (residues; elementarity is a torsion question).              *)
classifyReduced[P_, Nn_, R_, t_, n_, k_] := Module[{d, s, tau, a, ser, res = 0, i},
  If[P === 0 && Nn === 0, Return["zero", Module]];
  If[Nn =!= 0, Return["third-kind", Module]];          (* simple poles at non-branch points *)
  d = Exponent[R, t];
  If[Mod[d k, n] != 0, Return["second-kind", Module]]; (* infinity is a branch point *)
  s = d k/n;
  a = Coefficient[R, t, d];
  ser = Normal[Series[(Expand[(R /. t -> 1/tau) tau^d]/a)^(-k/n), {tau, 0, d + 1}]];
  (* coefficient of 1/t in P R^(-k/n) at infinity *)
  Do[i = j - s + 1;
     If[i >= 0, res += Coefficient[P, t, j] If[i == 0, ser /. tau -> 0, Coefficient[ser, tau, i]]],
     {j, 0, Exponent[P, t]}];
  If[isZero[res], "second-kind", "third-kind"]];
