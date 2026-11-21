/* Known gaps for test_numeric_stress.c — combinations where NEITHER the
 * machine nor the high-precision path produces a finite number.
 *
 * Included into a `static const KnownGap kKnownGaps[]` initializer. Each
 * entry is "<function label>@<argument label>" plus the reason it is not a
 * bug. The sweep fails on an unlisted gap AND on a listed gap that no
 * longer occurs, so this file is a live inventory of what large-argument
 * numerics does not cover — not a suppression list. Fix one of these and
 * the suite will tell you to delete the line.
 *
 * Note what is NOT here: a case where one path answers and the other does
 * not is never a gap, it is a divergence between two spellings of the same
 * computation, and the sweep reports it as a failure.
 */


/* Sinh — e^|x| past MPFR's exponent range (~10^323228458) on both paths */
{ "Sinh@2^53", "e^|x| past MPFR's exponent range (~10^323228458) on both paths" },
{ "Sinh@10^20", "e^|x| past MPFR's exponent range (~10^323228458) on both paths" },
{ "Sinh@2^53+1", "e^|x| past MPFR's exponent range (~10^323228458) on both paths" },
{ "Sinh@pi-digits", "e^|x| past MPFR's exponent range (~10^323228458) on both paths" },
{ "Sinh@10^25", "e^|x| past MPFR's exponent range (~10^323228458) on both paths" },
{ "Sinh@10^40", "e^|x| past MPFR's exponent range (~10^323228458) on both paths" },
{ "Sinh@-10^25", "e^|x| past MPFR's exponent range (~10^323228458) on both paths" },
{ "Sinh@10^20/3", "e^|x| past MPFR's exponent range (~10^323228458) on both paths" },

/* Cosh — e^|x| past MPFR's exponent range on both paths */
{ "Cosh@2^53", "e^|x| past MPFR's exponent range on both paths" },
{ "Cosh@10^20", "e^|x| past MPFR's exponent range on both paths" },
{ "Cosh@2^53+1", "e^|x| past MPFR's exponent range on both paths" },
{ "Cosh@pi-digits", "e^|x| past MPFR's exponent range on both paths" },
{ "Cosh@10^25", "e^|x| past MPFR's exponent range on both paths" },
{ "Cosh@10^40", "e^|x| past MPFR's exponent range on both paths" },
{ "Cosh@-10^25", "e^|x| past MPFR's exponent range on both paths" },
{ "Cosh@10^20/3", "e^|x| past MPFR's exponent range on both paths" },

/* Exp — e^x past MPFR's exponent range on both paths */
{ "Exp@2^53", "e^x past MPFR's exponent range on both paths" },
{ "Exp@10^20", "e^x past MPFR's exponent range on both paths" },
{ "Exp@2^53+1", "e^x past MPFR's exponent range on both paths" },
{ "Exp@pi-digits", "e^x past MPFR's exponent range on both paths" },
{ "Exp@10^25", "e^x past MPFR's exponent range on both paths" },
{ "Exp@10^40", "e^x past MPFR's exponent range on both paths" },
{ "Exp@10^20/3", "e^x past MPFR's exponent range on both paths" },

/* SinhIntegral — grows like e^x/x; past MPFR's exponent range */
{ "SinhIntegral@2^53", "grows like e^x/x; past MPFR's exponent range" },
{ "SinhIntegral@10^20", "grows like e^x/x; past MPFR's exponent range" },
{ "SinhIntegral@2^53+1", "grows like e^x/x; past MPFR's exponent range" },
{ "SinhIntegral@pi-digits", "grows like e^x/x; past MPFR's exponent range" },
{ "SinhIntegral@10^25", "grows like e^x/x; past MPFR's exponent range" },
{ "SinhIntegral@10^40", "grows like e^x/x; past MPFR's exponent range" },
{ "SinhIntegral@-10^25", "grows like e^x/x; past MPFR's exponent range" },
{ "SinhIntegral@10^20/3", "grows like e^x/x; past MPFR's exponent range" },

/* CoshIntegral — grows like e^x/x; past MPFR's exponent range */
{ "CoshIntegral@2^53", "grows like e^x/x; past MPFR's exponent range" },
{ "CoshIntegral@10^20", "grows like e^x/x; past MPFR's exponent range" },
{ "CoshIntegral@2^53+1", "grows like e^x/x; past MPFR's exponent range" },
{ "CoshIntegral@pi-digits", "grows like e^x/x; past MPFR's exponent range" },
{ "CoshIntegral@10^25", "grows like e^x/x; past MPFR's exponent range" },
{ "CoshIntegral@10^40", "grows like e^x/x; past MPFR's exponent range" },
{ "CoshIntegral@10^20/3", "grows like e^x/x; past MPFR's exponent range" },

/* ExpIntegralEi — grows like e^x/x; past MPFR's exponent range */
{ "ExpIntegralEi@2^53", "grows like e^x/x; past MPFR's exponent range" },
{ "ExpIntegralEi@10^20", "grows like e^x/x; past MPFR's exponent range" },
{ "ExpIntegralEi@2^53+1", "grows like e^x/x; past MPFR's exponent range" },
{ "ExpIntegralEi@pi-digits", "grows like e^x/x; past MPFR's exponent range" },
{ "ExpIntegralEi@10^25", "grows like e^x/x; past MPFR's exponent range" },
{ "ExpIntegralEi@10^40", "grows like e^x/x; past MPFR's exponent range" },
{ "ExpIntegralEi@10^20/3", "grows like e^x/x; past MPFR's exponent range" },

/* Gamma — overflows; declines as unevaluated on one path and ComplexInfinity on the other */
{ "Gamma@2^53", "overflows; declines as unevaluated on one path and ComplexInfinity on the other" },
{ "Gamma@10^20", "overflows; declines as unevaluated on one path and ComplexInfinity on the other" },
{ "Gamma@2^53+1", "overflows; declines as unevaluated on one path and ComplexInfinity on the other" },
{ "Gamma@pi-digits", "overflows; declines as unevaluated on one path and ComplexInfinity on the other" },
{ "Gamma@10^25", "overflows; declines as unevaluated on one path and ComplexInfinity on the other" },
{ "Gamma@10^40", "overflows; declines as unevaluated on one path and ComplexInfinity on the other" },
{ "Gamma@-10^25", "overflows; declines as unevaluated on one path and ComplexInfinity on the other" },
{ "Gamma@10^20/3", "overflows; declines as unevaluated on one path and ComplexInfinity on the other" },

/* Beta — Beta[n, 2] = 1/(n(n+1)) via Pochhammer, which overflows first */
{ "Beta@2^53", "Beta[n, 2] = 1/(n(n+1)) via Pochhammer, which overflows first" },
{ "Beta@10^20", "Beta[n, 2] = 1/(n(n+1)) via Pochhammer, which overflows first" },
{ "Beta@2^53+1", "Beta[n, 2] = 1/(n(n+1)) via Pochhammer, which overflows first" },
{ "Beta@pi-digits", "Beta[n, 2] = 1/(n(n+1)) via Pochhammer, which overflows first" },
{ "Beta@10^25", "Beta[n, 2] = 1/(n(n+1)) via Pochhammer, which overflows first" },
{ "Beta@10^40", "Beta[n, 2] = 1/(n(n+1)) via Pochhammer, which overflows first" },

/* PolyLog — no large-|z| continuation implemented; stays unevaluated */

/* HarmonicNumber — reduces to PolyGamma[0, z] at large negative z, which is unimplemented there */
{ "HarmonicNumber@-10^25", "reduces to PolyGamma[0, z] at large negative z, which is unimplemented there" },

/* Erfi — series does not terminate at large argument: times out on both paths */
{ "Erfi@2^53", "series does not terminate at large argument: times out on both paths" },
{ "Erfi@10^20", "series does not terminate at large argument: times out on both paths" },
{ "Erfi@2^53+1", "series does not terminate at large argument: times out on both paths" },
{ "Erfi@pi-digits", "series does not terminate at large argument: times out on both paths" },
{ "Erfi@10^25", "series does not terminate at large argument: times out on both paths" },
{ "Erfi@10^40", "series does not terminate at large argument: times out on both paths" },
{ "Erfi@-10^25", "series does not terminate at large argument: times out on both paths" },
{ "Erfi@10^20/3", "series does not terminate at large argument: times out on both paths" },

/* AiryBi — grows like e^(2/3 x^(3/2)); no large-argument path, stays unevaluated */
{ "AiryBi@2^53", "grows like e^(2/3 x^(3/2)); no large-argument path, stays unevaluated" },
{ "AiryBi@10^20", "grows like e^(2/3 x^(3/2)); no large-argument path, stays unevaluated" },
{ "AiryBi@2^53+1", "grows like e^(2/3 x^(3/2)); no large-argument path, stays unevaluated" },
{ "AiryBi@pi-digits", "grows like e^(2/3 x^(3/2)); no large-argument path, stays unevaluated" },
{ "AiryBi@10^25", "grows like e^(2/3 x^(3/2)); no large-argument path, stays unevaluated" },
{ "AiryBi@10^40", "grows like e^(2/3 x^(3/2)); no large-argument path, stays unevaluated" },
{ "AiryBi@10^20/3", "grows like e^(2/3 x^(3/2)); no large-argument path, stays unevaluated" },

/* BesselI0 — grows like e^x/sqrt(x); no large-argument path, stays unevaluated */
{ "BesselI0@2^53", "grows like e^x/sqrt(x); no large-argument path, stays unevaluated" },
{ "BesselI0@10^20", "grows like e^x/sqrt(x); no large-argument path, stays unevaluated" },
{ "BesselI0@2^53+1", "grows like e^x/sqrt(x); no large-argument path, stays unevaluated" },
{ "BesselI0@pi-digits", "grows like e^x/sqrt(x); no large-argument path, stays unevaluated" },
{ "BesselI0@10^25", "grows like e^x/sqrt(x); no large-argument path, stays unevaluated" },
{ "BesselI0@10^40", "grows like e^x/sqrt(x); no large-argument path, stays unevaluated" },
{ "BesselI0@-10^25", "grows like e^x/sqrt(x); no large-argument path, stays unevaluated" },
{ "BesselI0@10^20/3", "grows like e^x/sqrt(x); no large-argument path, stays unevaluated" },

/* BernoulliB — defined for integer index only; a non-integer index stays symbolic */
{ "BernoulliB@22/7", "defined for integer index only; a non-integer index stays symbolic" },
{ "BernoulliB@-3/2", "defined for integer index only; a non-integer index stays symbolic" },

/* EulerE — defined for integer index only; a non-integer index stays symbolic */
{ "EulerE@22/7", "defined for integer index only; a non-integer index stays symbolic" },
{ "EulerE@-3/2", "defined for integer index only; a non-integer index stays symbolic" },

/* Underflow (v0.66, 1ff092b98): a value below MPFR's exponent range is
 * Underflow[] with General::unfl at N[.., p], and machine N[] reports the machine
 * 0. (N[Underflow[]] == 0., v0.6660) -- the mirror image of the
 * Sinh/Cosh/Exp overflow gaps above, and the reference's answer too. */

/* Sech */
{ "Sech@2^53", "e^-|x| below MPFR's exponent range (~10^-323228458): machine 0., Underflow[] at 30 digits, as the reference" },
{ "Sech@10^20", "e^-|x| below MPFR's exponent range (~10^-323228458): machine 0., Underflow[] at 30 digits, as the reference" },
{ "Sech@2^53+1", "e^-|x| below MPFR's exponent range (~10^-323228458): machine 0., Underflow[] at 30 digits, as the reference" },
{ "Sech@pi-digits", "e^-|x| below MPFR's exponent range (~10^-323228458): machine 0., Underflow[] at 30 digits, as the reference" },
{ "Sech@10^25", "e^-|x| below MPFR's exponent range (~10^-323228458): machine 0., Underflow[] at 30 digits, as the reference" },
{ "Sech@10^40", "e^-|x| below MPFR's exponent range (~10^-323228458): machine 0., Underflow[] at 30 digits, as the reference" },
{ "Sech@-10^25", "e^-|x| below MPFR's exponent range (~10^-323228458): machine 0., Underflow[] at 30 digits, as the reference" },
{ "Sech@10^20/3", "e^-|x| below MPFR's exponent range (~10^-323228458): machine 0., Underflow[] at 30 digits, as the reference" },

/* Csch */
{ "Csch@2^53", "e^-|x| below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "Csch@10^20", "e^-|x| below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "Csch@2^53+1", "e^-|x| below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "Csch@pi-digits", "e^-|x| below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "Csch@10^25", "e^-|x| below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "Csch@10^40", "e^-|x| below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "Csch@-10^25", "e^-|x| below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "Csch@10^20/3", "e^-|x| below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },

/* Erfc */
{ "Erfc@2^53", "e^-x^2 below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "Erfc@10^20", "e^-x^2 below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "Erfc@2^53+1", "e^-x^2 below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "Erfc@pi-digits", "e^-x^2 below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "Erfc@10^25", "e^-x^2 below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "Erfc@10^20/3", "e^-x^2 below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },

/* ExpIntegralEi */
{ "ExpIntegralEi@-10^25", "e^x/x at x = -10^25 below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },

/* AiryAi */
{ "AiryAi@2^53", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },
{ "AiryAi@10^20", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },
{ "AiryAi@2^53+1", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },
{ "AiryAi@pi-digits", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },
{ "AiryAi@10^25", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },
{ "AiryAi@10^20/3", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },

/* AiryAiPrime */
{ "AiryAiPrime@2^53", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },
{ "AiryAiPrime@10^20", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },
{ "AiryAiPrime@2^53+1", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },
{ "AiryAiPrime@pi-digits", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },
{ "AiryAiPrime@10^25", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },
{ "AiryAiPrime@10^20/3", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },

/* BesselK0 */
{ "BesselK0@2^53", "e^-x below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "BesselK0@10^20", "e^-x below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "BesselK0@2^53+1", "e^-x below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "BesselK0@pi-digits", "e^-x below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "BesselK0@10^25", "e^-x below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "BesselK0@10^40", "e^-x below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "BesselK0@10^20/3", "e^-x below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },

/* The same underflow at 10^40, where machine N[] already answered 0. through
 * the raised MPFR pass; since v0.6660 both magnitudes give machine 0. and
 * N[.., 30] Underflow[], exactly as the reference. */
{ "Erfc@10^40", "e^-x^2 below MPFR's exponent range: machine 0., Underflow[] at 30 digits, as the reference" },
{ "AiryAi@10^40", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., Underflow[] at 30 digits" },
{ "AiryAiPrime@10^40", "e^-(2/3)x^(3/2) below MPFR's exponent range: machine 0., -Underflow[] at 30 digits" },
