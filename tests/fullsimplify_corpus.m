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

    (* --- log of trig/hyperbolic combinations (logexp.m): every rule in every
           configuration -- both signs, both argument orders, renamed variable,
           integer and rational coefficients (factored '+' and distributed '-'),
           half-angle forms, and non-firing guards that must stay unchanged. --- *)

    (* Sec/Tan -> ArcTanh[Sin] *)
    {Log[Sec[x] + Tan[x]],                  ArcTanh[Sin[x]]},
    {Log[Tan[x] + Sec[x]],                  ArcTanh[Sin[x]]},
    {Log[Sec[t] + Tan[t]],                  ArcTanh[Sin[t]]},
    {Log[2 Sec[x] + 2 Tan[x]],              Log[2] + ArcTanh[Sin[x]]},
    {Log[Sec[x]/2 + Tan[x]/2],              Log[1/2] + ArcTanh[Sin[x]]},
    {Log[Sec[x] - Tan[x]],                  -ArcTanh[Sin[x]]},
    {Log[2 Sec[x] - 2 Tan[x]],              Log[2] - ArcTanh[Sin[x]]},

    (* Csc/Cot -> ArcTanh[Cos] *)
    {Log[Csc[x] + Cot[x]],                  ArcTanh[Cos[x]]},
    {Log[Cot[x] + Csc[x]],                  ArcTanh[Cos[x]]},
    {Log[2 Csc[x] + 2 Cot[x]],              Log[2] + ArcTanh[Cos[x]]},
    {Log[Csc[x] - Cot[x]],                  -ArcTanh[Cos[x]]},
    {Log[2 Csc[x] - 2 Cot[x]],              Log[2] - ArcTanh[Cos[x]]},

    (* Cosh/Sinh -> linear *)
    {Log[Cosh[x] + Sinh[x]],                x},
    {Log[2 Cosh[x] + 2 Sinh[x]],            Log[2] + x},
    {Log[Cosh[x] - Sinh[x]],                -x},
    {Log[2 Cosh[x] - 2 Sinh[x]],            Log[2] - x},

    (* Coth/Csch -> ArcCoth[Cosh] *)
    {Log[Coth[x] + Csch[x]],                ArcCoth[Cosh[x]]},
    {Log[Csch[x] + Coth[x]],                ArcCoth[Cosh[x]]},
    {Log[2 Coth[x] + 2 Csch[x]],            Log[2] + ArcCoth[Cosh[x]]},
    {Log[Coth[x] - Csch[x]],                -ArcCoth[Cosh[x]]},
    {Log[2 Coth[x] - 2 Csch[x]],            Log[2] - ArcCoth[Cosh[x]]},

    (* Half-angle tangents/cotangents *)
    {Log[Tan[x/2]],                         -ArcTanh[Cos[x]]},
    {Log[Cot[x/2]],                         ArcTanh[Cos[x]]},
    {Log[Tanh[x/2]],                        -ArcCoth[Cosh[x]]},
    {Log[Coth[x/2]],                        ArcCoth[Cosh[x]]},
    {Log[3 Tan[x/2]],                       Log[3] - ArcTanh[Cos[x]]},
    {Log[Cot[t/2]],                         ArcTanh[Cos[t]]},

    (* Non-firing guards: unequal coefficients, mismatched angles, wrong
       fractional angle -- must be returned unchanged. *)
    {Log[2 Sec[x] + 3 Tan[x]],              Log[2 Sec[x] + 3 Tan[x]]},
    {Log[Sec[x] + Tan[y]],                  Log[Sec[x] + Tan[y]]},
    {Log[Tan[x/3]],                         Log[Tan[x/3]]},

    (* --- >= Simplify sanity: FullSimplify never does worse --- *)
    {(x - 1) (x + 1) (x^2 + 1) + 1,         x^4},
    {Sin[x]^2 + Cos[x]^2,                   1}
}
