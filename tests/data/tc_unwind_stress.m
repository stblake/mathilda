(* TimeConstrained timeouts landing at random points of PossibleZeroQ's sampler
   (tools/tc_unwind_stress.py runs this in parallel copies).  Until v0.6664 the
   SIGPROF backstop siglongjmp'd out of evaluator-heavy C code that discards the
   Abort[] it gets back, leaving shared state torn: ~6% of runs segfaulted under
   load (expr_alloc_node, expr_hash, eval_valueless_symbol). *)
e = (Sin[x]^2 + Cos[x]^2 - 1) (x^7 + 3 x - Sqrt[2])^5 + Exp[x] Log[1 + x^2] - Exp[x] Log[1 + x^2];
n = 0;
Do[TimeConstrained[PossibleZeroQ[e + k (Sin[x+k]^2 + Cos[x+k]^2 - 1)], 10^-4 (1 + Mod[k, 37])]; n++, {k, 1, 3000}];
SeedRandom[5]; a = RandomInteger[10^30, 3];
SeedRandom[5]; PossibleZeroQ[Sin[x]^2 + Cos[x]^2 - 1]; b = RandomInteger[10^30, 3];
Print[n, " ", a === b];
