(* Experiment 30 -- Schur decomposition.
   MEASURES src/linalg/schurdecomp*.c: the standard form (dgees/zgees) and the
   generalized QZ form (dgges/zgges), plus the packed/NDArray fast path.
   NumPy/SciPy link the SAME Accelerate LAPACK Mathilda does, so any spread on
   these rows is pure Expr<->buffer marshalling -- and the "packed" row shows it
   gone (a visible NDArray input stays on the buffer end to end, no boxing).
   Checks reconstruct a reproducible symmetric matrix; a correct decomposition
   gives residual 0 regardless of the (non-unique) Schur basis. *)

Get["../harness.m"];
Get["../data.m"];

require[{"SchurDecomposition"}];

s8  = N[spdIntMatrix[8]];
s8b = N[spdIntMatrix[8] + IdentityMatrix[8]];
c8  = s8 + I s8b;   (* reproducible complex matrix for the zgees check *)

schurResid[m_] := Module[{q, t}, {q, t} = SchurDecomposition[m];
  Round[10^6 Max[Abs[Flatten[m - q . t . ConjugateTranspose[q]]]]]];
genResid[m_, a_] := Module[{q, s, p, t}, {q, s, p, t} = SchurDecomposition[{m, a}];
  Round[10^6 Max[Abs[Flatten[Join[m - q . s . ConjugateTranspose[p],
                                   a - q . t . ConjugateTranspose[p]]]]]]];

(* rand01 auto-packs to a transparent packed-list, so this IS the packed fast
   path (data stays on the buffer, factors come back packed, no Expr boxing). *)
bench["SchurDecomposition 300x300", SchurDecomposition[rand01[{300, 300}]], 1];
check["SchurDecomposition 300x300", schurResid[s8]];

bench["SchurDecomposition complex 200x200",
  SchurDecomposition[rand01[{200, 200}] + I rand01[{200, 200}]], 1];
check["SchurDecomposition complex 200x200", schurResid[c8]];

bench["SchurDecomposition generalized 200x200",
  SchurDecomposition[{rand01[{200, 200}], rand01[{200, 200}]}], 1];
check["SchurDecomposition generalized 200x200", genResid[s8, s8b]];
