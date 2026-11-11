(* Deficiencies in PolynomialGCD / Together over a nested-radical field,
 * surfaced by issue #133.  Run: ./Mathilda -file tasks/polygcd_together_deficiencies/repro.m
 * Each line prints the CURRENT answer; see README.md for the correct one. *)

(* D1: PolynomialGCD misses a non-trivial common factor -> 1 (wrong).
 *     Correct: d1 up to a unit. *)
q = Sqrt[3/4 - 1/Sqrt[2]];
d1 = 7 + 5 Sqrt[2] - 4 x q;
d2 = 3 - Sqrt[2] x + q^3;
Print["D1  ", PolynomialGCD[Expand[d1 d2], Expand[d1^2]]];
Print["D1x ", PolynomialGCD[Expand[d1 d2], Expand[d1^2], Extension -> Automatic]];

(* D2: Together does not finish on a 4639-leaf rational function of x whose
 *     coefficients lie in Q(I, Sqrt[2], Sqrt[3/4 - 1/Sqrt[2]]).  Its value at
 *     x = 1/3 is 1.4714045207910316829. *)
e = Get[FileNameJoin[{DirectoryName[$InputFileName], "d2_together_input.m"}]];
Print["D2  ", AbsoluteTiming[TimeConstrained[Together[e], 20, "TIMEOUT"]]];
Print["D2x ", AbsoluteTiming[TimeConstrained[Together[e, Extension -> Automatic], 20, "TIMEOUT"]]];
