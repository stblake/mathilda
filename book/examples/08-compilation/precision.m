# Section 8.7: machine integers, overflow, and arbitrary precision.
# By default a compiled integer that overflows the machine word is detected and
# the call is abandoned, so the interpreter re-runs it and promotes to a bignum:
# the compiled answer is always the interpreter's.
Compile[{{n, _Integer}}, n^3][3000000]
3000000^3
# WorkingPrecision -> n compiles the arithmetic in MPFR at n digits; the result
# equals what the interpreter would compute at the same precision.
Compile[{{x, _Real}}, Sqrt[x], WorkingPrecision -> 40][2]
N[Sqrt[2], 40]
CompileDiagnostics[{{x, _Real}}, Sqrt[x], WorkingPrecision -> 40]
# "BigIntegers" -> True keeps integer arithmetic exact (GMP) inside the compiled
# code, so it stays compiled AND correct rather than bailing on overflow.
Compile[{{n, _Integer}}, n^3, "BigIntegers" -> True][3000000]
CompileDiagnostics[{{n, _Integer}}, n^3, "BigIntegers" -> True]
# RuntimeOptions -> "Speed" skips the overflow check: an overflowing result is
# left wrapped modulo 2^64 -- a different, generally wrong answer -- instead of
# promoting. The wrapped value is exactly the modular reduction.
Compile[{{n, _Integer}}, n^3, RuntimeOptions -> "Speed"][3000000]
Mod[3000000^3, 2^64]
