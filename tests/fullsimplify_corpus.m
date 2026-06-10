(* FullSimplify per-rule corpus.

   A single List of {input, expected} pairs. The runner (test_fullsimplify_-
   corpus.c) evaluates FullSimplify[input] and checks it is structurally equal
   to expected. Each entry exercises an identity that the relevance engine
   selects and that genuinely lowers complexity (so it is actually taken),
   plus a couple of >=-Simplify sanity cases. *)
{
    (* --- gamma family: recurrences (genuine gaps vs Simplify) --- *)
    {Gamma[x + 1]/Gamma[x],                 x},
    {Gamma[x + 1] - x Gamma[x],             0},
    {LogGamma[x + 1] - LogGamma[x],         Log[x]},
    {PolyGamma[0, x + 1] - PolyGamma[0, x], 1/x},

    (* --- error function: complementary identity --- *)
    {Erf[x] + Erfc[x],                      1},

    (* --- polylogarithm: dilogarithm duplication --- *)
    {PolyLog[2, z] + PolyLog[2, -z],        PolyLog[2, z^2]/2},

    (* --- real radical: Surd[x,n]^n == x --- *)
    {Surd[x, 3]^3,                          x},

    (* --- >= Simplify sanity: FullSimplify never does worse --- *)
    {(x - 1) (x + 1) (x^2 + 1) + 1,         x^4},
    {Sin[x]^2 + Cos[x]^2,                   1}
}
