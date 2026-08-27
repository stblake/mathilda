# Section 4.5.6: NLimit -- limits recovered by sequence acceleration.
NLimit[(1 - Cos[x])/x^2, x -> 0, AccuracyGoal -> 10]
NLimit[Zeta[s] - 1/(s - 1), s -> 1, AccuracyGoal -> 10]
NLimit[n (2^(1/n) - 1), n -> Infinity, Method -> "SequenceLimit", AccuracyGoal -> 10]
