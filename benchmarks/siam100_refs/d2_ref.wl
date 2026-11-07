(* D2 reference, method B (independent of d2_ref.py's digamma closed form):
   the first 2^16 - 1 terms are summed DIRECTLY (term by term at 110 digits,
   no grouping); the tail n >= 2^16 is grouped by bit length k >= 17 with block
   sums B(k) = PolyGamma[0, 2^k] - PolyGamma[0, 2^(k-1)] (numeric argument:
   an exact integer argument expands to a giant HarmonicNumber), and the self-similar identity
   Sum_{k>=17} 1/f(k) = S - Sum_{n=1}^{16} 1/f(n) turns the tail into
   Log2 (S - P16) + Sum_{k>=17} (B(k) - Log2)/f(k); S then solves a linear
   equation (FindRoot at 100 digits). *)
f[n_Integer] := If[n <= 2, n, n f[BitLength[n]]];
direct = Total[Table[N[1/f[n], 110], {n, 1, 2^16 - 1}]];
p16 = Sum[1/f[n], {n, 1, 16}];
em = Total[Table[(PolyGamma[0, N[2^k, 120]] - PolyGamma[0, N[2^(k - 1), 120]] - Log[N[2, 120]])/f[k], {k, 17, 400}]];
sol = s /. FindRoot[s == direct + Log[2] (s - p16) + em, {s, 2}, WorkingPrecision -> 100];
Print[NumberForm[sol, 80]];
