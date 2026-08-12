(* ==========================================================================
   benchmarks/data.m -- shared, deterministic input generators.
   ==========================================================================

   Loaded by an experiment right after the harness:

       Get["../harness.m"]; Get["../data.m"];

   WHY THIS IS SHARED.  The value gate is only meaningful if both halves of an
   experiment compute over the SAME inputs.  Putting the generators in one place
   per language makes that true by construction rather than by two authors
   agreeing on a formula; `data.py` builds the identical objects.

   Every generator here is EXACTLY reproducible -- integer coefficients, closed
   forms, no floating-point randomness -- so a check value can be compared
   across three systems.  Where an experiment needs random float data for its
   TIMING (a 10^6-element buffer), it draws that locally with `rand01` from the
   harness and checks a separate small deterministic case instead.  The three
   systems cannot be made to draw the same float stream, so a float sum over
   random data is never a check.
   ========================================================================== *)

(* ---- univariate polynomials -------------------------------------------- *)

(* NEVER `Product[body, {i, 1, k}]` OR `Sum[body, {k, 0, n}]` FOR THESE.
   Both are symbolic-summation calls, and `Product` legitimately returns a
   CLOSED FORM rather than an expanded polynomial:

       Product[x^2 + i x + (i+1), {i, 1, 2}]
         -> (1 + x)^2 (1 + (1 + x^2)/(1 + x)) (2 + (1 + x^2)/(1 + x))

   which is the right value in the wrong shape.  `Exponent` of that reads 30
   for k = 10 instead of 20, and `PolynomialGCD` handed two of them returns
   unevaluated -- which showed up as Mathilda "finding no common factor" where
   sympy found a degree-20 one.  That looked exactly like a PolynomialGCD bug
   and was not one.

   `Times @@ Table[...]` and `Total[Table[...]]` build the object directly, and
   the wrapping `Expand` pins a canonical form, so both systems hold the same
   polynomial and `Exponent` means what it says. *)

(* A dense degree-n polynomial with small, reproducible integer coefficients.
   Deliberately not random: mild coefficient growth keeps Factor's cost about
   degree, not about coefficient size. *)
denseUPoly[n_Integer, x_] :=
  Total[Table[Mod[k^2 + 3 k + 1, 17] x^k, {k, 0, n}]] + 1;

(* A product of distinct irreducible-ish quadratics -- a factoring target whose
   answer is known by construction. Expanded, so its degree is 2k as written. *)
factorableUPoly[k_Integer, x_] :=
  Expand[Times @@ Table[x^2 + i x + (i + 1), {i, 1, k}]];

(* A sparse high-degree polynomial: cheap to store, expensive to factor. *)
sparseUPoly[n_Integer, x_] := x^n + x^(Quotient[n, 2]) + 1;

(* ---- rational functions ------------------------------------------------ *)

(* A proper rational function with a squarefree denominator of degree n: the
   canonical Hermite-reduction / partial-fraction target. *)
ratFn[n_Integer, x_] :=
  denseUPoly[n - 1, x] / (Times @@ Table[x - i, {i, 1, n}]);

(* A rational function whose denominator has repeated factors, which is the
   branch Hermite reduction actually exists for. *)
ratFnRepeated[n_Integer, x_] :=
  denseUPoly[n, x] / (Times @@ Table[(x - i)^2, {i, 1, n}]);

(* ---- multivariate ------------------------------------------------------- *)

multiPoly[k_Integer, vars_List] :=
  Expand[Times @@ Table[Total[vars] + i, {i, 1, k}]];

(* ---- polynomial systems (Groebner targets) ----------------------------- *)

(* cyclic-n: the standard Groebner basis stress family. *)
cyclicSystem[vars_List] := Module[{n = Length[vars]},
  Append[
    Table[Total[Table[Times @@ Table[vars[[Mod[i + k, n, 1]]], {k, 0, d - 1}],
                      {i, 1, n}]], {d, 1, n - 1}],
    (Times @@ vars) - 1]];

(* ---- matrices ---------------------------------------------------------- *)

(* A symmetric positive-definite integer matrix: exactly reproducible, so
   LinearSolve / Eigenvalues / Det checks are comparable across systems. *)
spdIntMatrix[n_Integer] :=
  Table[If[i === j, n + 2, 1/(1 + Abs[i - j])], {i, n}, {j, n}];

(* An exactly-representable ill-conditioned classic. *)
hilbertLike[n_Integer] := Table[1/(i + j - 1), {i, n}, {j, n}];

(* ---- integer data ------------------------------------------------------ *)

(* Reproducible pseudo-random integers WITHOUT depending on either system's
   RNG stream: a linear congruential sequence, identical in data.py.  This is
   what makes a Sort or a Total check comparable across three systems. *)
lcgInts[n_Integer, modulus_Integer] :=
  Module[{s = 12345},
    Table[s = Mod[1103515245 s + 12345, 2147483648]; Mod[s, modulus] + 1,
          {n}]];
