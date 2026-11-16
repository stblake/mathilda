(* Reproducer: the C transcendental Risch cannot certify non-elementary
   polynomials in ArcTan[u] of degree >= 2 (README.md). Run from the repo root:
   ./Mathilda -file tasks/risch_arctan_nonelem/repro.m
   Columns: the real form through Integrate`RischTranscendental, the TrigToExp
   (two-log, Q(i)) form through it, Risch`ElementaryIntegralQ on the real form,
   and ParallelMixedTower (the oracle). Times in seconds. *)
st[r_, h_] := Which[Head[r] === h, "UNDEC", ListQ[r], If[First[r] === "not elementary", "CERT", "FAIL"], r === False, "CERT", r === True, "ELEMQ", FreeQ[r, Integrate | Integrate`RischTranscendental], "ELEM", True, "UNDEC"];
Do[With[{f = u^m ArcTan[u]^k}, With[{g = TrigToExp[f]},
  With[{a = AbsoluteTiming[TimeConstrained[Quiet[Integrate`RischTranscendental[f, u]], 30, $TO]],
        b = AbsoluteTiming[TimeConstrained[Quiet[Integrate`RischTranscendental[g, u]], 30, $TO]],
        c = AbsoluteTiming[TimeConstrained[Quiet[Risch`ElementaryIntegralQ[f, u]], 30, $TO]],
        d = AbsoluteTiming[TimeConstrained[Quiet[Integrate`ParallelMixedTower[f, u]], 30, $TO]]},
   Print[InputForm[f], " | real Risch ", st[a[[2]], Integrate`RischTranscendental], " ", Round[a[[1]], 0.01],
     " | 2-log Risch ", If[b[[2]] === $TO, "TIMEOUT", st[b[[2]], Integrate`RischTranscendental]], " ", Round[b[[1]], 0.01],
     " | EIQ(real) ", c[[2]], " ", Round[c[[1]], 0.01],
     " | PMT ", st[d[[2]], Integrate`ParallelMixedTower], " ", Round[d[[1]], 0.01]]]]],
  {k, 1, 3}, {m, {0, 1, 5}}];
