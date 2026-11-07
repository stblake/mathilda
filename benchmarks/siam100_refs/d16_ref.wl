(* Reference B for SIAM100 App. D problem 16 (exp 191) -- Wolfram Language,
   independent implementation of the sine-series solution:
   b_k = -(2/3) Im phihat(k Pi/3), phihat by the two-scale infinite product
   (J factors, plain truncation, no tail correction -- J chosen so the
   neglected factors are 1 + O(2^-J w) < 10^-80), u_x on the trapezoid grid by
   FourierDCT (type I) at precision 80, trapezoid sum, FindRoot (secant) in t.
   Run: wolframscript -file d16_ref.wl            (length 5.4)
        wolframscript -file d16_ref.wl reduced    (length 4.5) *)
prec = 80;
(* target length and secant start: 5.4 (the problem) or 4.5 (reduced row) *)
{target, tApprox, M} = If[MemberQ[$ScriptCommandLine, "reduced"],
  {45/10, 53943397/10^10, 2^12}, {54/10, 5212782547/10^16, 2^19}];
p = {1 + Sqrt[3], 3 + Sqrt[3], 3 - Sqrt[3], 1 - Sqrt[3]}/8;
K = Ceiling[Sqrt[(prec Log[10.] + 25)/(0.9 tApprox)] 3/Pi];
J = Ceiling[Log2[K 3.2]] + Ceiling[prec Log2[10.]] + 10;
ks = Range[K];
w = N[ks Pi/3, prec + 20];
ph = ConstantArray[N[1, prec + 20], K];
Do[z = Exp[-I w/2^j]; ph = ph (p[[1]] + p[[2]] z + p[[3]] z^2 + p[[4]] z^3), {j, 1, J}];
b = -2/3 Im[ph];
Print["K=", K, " J=", J, " b1=", N[b[[1]], 30]];
L[t_?NumericQ] := Module[{c, u, f},
  c = Join[{0}, b w Exp[-w^2 t], ConstantArray[0, M - K]];
  u = FourierDCT[c, 1] Sqrt[M/2];
  f = Sqrt[1 + u^2];
  3/M (Total[f] - (First[f] + Last[f])/2)];
sol = FindRoot[L[t] - target, {t, SetPrecision[tApprox, prec], SetPrecision[tApprox (1 + 10^-8), prec]},
  WorkingPrecision -> prec, AccuracyGoal -> Infinity, PrecisionGoal -> 70, MaxIterations -> 30];
Print["t* = ", N[t /. sol, 75]];
