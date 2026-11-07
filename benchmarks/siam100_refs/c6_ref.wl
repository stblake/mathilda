(* c6 reference (method 2): the book's closed form (6.1):
   eps = Sqrt[1 - eta^2]/4 with AGM(Sqrt[4-(eta-1)^2], Sqrt[4-(eta+1)^2]) == 1 *)
eta = eta /. FindRoot[ArithmeticGeometricMean[Sqrt[4 - (eta - 1)^2], Sqrt[4 - (eta + 1)^2]] == 1,
   {eta, 99/100}, WorkingPrecision -> 120, PrecisionGoal -> 100, AccuracyGoal -> 100];
Print[N[Sqrt[1 - eta^2]/4, 70]];
