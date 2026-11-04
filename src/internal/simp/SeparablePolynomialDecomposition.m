(* ====================================================================== *)
(* SeparablePolynomialDecomposition.m -- minimal sum of separable products. *)
(*                                                                          *)
(* Loaded lazily by Simplify's polynomial pipeline (try_separable_decomp in *)
(* src/simp/simp_search.c) the first time it sees a polynomial in two or    *)
(* more variables; it can also be loaded by hand with                       *)
(* LoadModule["simp/SeparablePolynomialDecomposition.m"].  Simplify tries    *)
(* every bipartition {S, T} of the variables and keeps the decomposition     *)
(* only when it is STRICTLY simpler than Expand / Factor / FactorTerms.      *)
(* ====================================================================== *)

BeginPackage["SeparablePolynomialDecomposition`"];

SeparablePolynomialDecomposition::usage =
  "SeparablePolynomialDecomposition[p, S, T] rewrites the polynomial p as a minimal-length sum of \
products B[S] Q[T], with each factor returned in factored form.";

Begin["`Private`"];

(*
SeparablePolynomialDecomposition[p, S, T]

Decompose a polynomial into a minimal sum of separable products.

Write ``p`` as ``Sum[B[[j]] Q[[j]], {j, r}]`` where each ``B[[j]]`` is a
polynomial in the variables ``S`` only and each ``Q[[j]]`` a polynomial in
the variables ``T`` only, with ``r`` as small as possible for that split.
Every factor is then passed through ``Factor``.

Parameters
----------
p : polynomial
    Expanded polynomial in the variables ``Join[S, T]`` with rational
    coefficients.  Non-polynomial atoms are treated as constants.
S : list of symbols
    First block of variables.
T : list of symbols
    Second block of variables.  ``S`` and ``T`` should be disjoint and
    together cover the variables of ``p``.

Returns
-------
expr : polynomial
    An expression equal to ``p``, of the form
    ``Factor[B1] Factor[Q1] + ... + Factor[Br] Factor[Qr]``.
    ``0`` when ``p`` is zero.  $Failed when ``p`` is not a polynomial in
    ``Join[S, T]`` or when ``S`` and ``T`` are not disjoint, non-empty lists
    (the Simplify stage then declines).

Inexact input
-------------
A ``p`` with inexact (machine or arbitrary-precision) numbers is first
rationalised with ``Rationalize[p, 0]`` and the exact form decomposed (the
result is then exact).  If inexact numbers survive ``Rationalize``, ``p``
itself is returned unchanged.  Simplify does not use this path: its stage
gates on exact input, so Simplify never changes the precision of an
inexact expression.

Notes
-----
Let ``M`` be the coefficient matrix of ``p`` with rows indexed by
``T``-monomials and columns by ``S``-monomials (the ``(S, T)``-flattening).
Then ``r = MatrixRank[M]`` is the minimal number of separable terms for
this split, and the decomposition is unique up to a change of basis in
``GL(r, Q)``.  The basis is chosen greedily: coefficient polynomials of
``p`` with respect to ``T`` are sorted by number of terms, and a maximal
independent subset is taken in that order via the pivot columns of the
reduced row echelon form of ``Transpose[M]``.  This favours bases that
``Factor`` recognises, but no optimality of the final expression size is
guaranteed (minimal formula size is NP-hard in general).

Complexity is dominated by ``RowReduce`` on an ``m x n`` rational matrix,
where ``m`` and ``n`` are the numbers of distinct ``T``- and
``S``-monomials.

If the greedy choice lands on a scalar multiple of a natural basis
element, the corresponding ``Q[[j]]`` acquires rational coefficients;
normalising ``cs`` with ``FactorTermsList`` beforehand removes this.

The pivot test is ``# === 0 &`` (an exact zero entry of the RREF row);
EqualTo[0] is not available in Mathilda.

See Also
--------
Factor, Collect, CoefficientRules, RowReduce

Examples
--------
>>> p = Expand[(3 a b - b^2)^3 (x^2 - y + 1)^2 - (a - b)^2 (x - y)];
>>> SeparablePolynomialDecomposition[p, {a, b}, {x, y}]
-(a - b)^2 (x - y) + b^3 (3 a - b)^3 (1 + x^2 - y)^2
*)

(* Inexact input: decompose the exact form Rationalize[p, 0]; if any inexact
   number survives rationalisation, answer p unchanged. *)
SeparablePolynomialDecomposition[p_, S_List, T_List] /;
   !FreeQ[p, _?InexactNumberQ] :=
 With[{q = Rationalize[p, 0]},
  If[FreeQ[q, _?InexactNumberQ], SeparablePolynomialDecomposition[q, S, T], p]]

SeparablePolynomialDecomposition[p_, S_List, T_List] /;
   S =!= {} && T =!= {} && Intersection[S, T] === {} &&
   FreeQ[p, _?InexactNumberQ] && PolynomialQ[p, Join[S, T]] :=
 Module[{r, mons, cs, crs, smons, M, R, piv, B, Q},
  (* Coefficient rules w.r.t. T, sparsest coefficient polynomial first. *)
  r = SortBy[CoefficientRules[p, T], Length@MonomialList@Last@# &];
  mons = Times @@ (T^First@#) & /@ r;   (* T-monomials *)
  cs = Last /@ r;                         (* Their coefficients, polys in S. *)
  (* (S,T)-flattening: rows = T-monomials, columns = S-monomials.
     CoefficientRules[c, S] is computed ONCE per row and every column read
     from its Association in one Lookup: recomputing it per (row, column)
     entry made the matrix O(m n) CoefficientRules calls (3.8 s of a 4.8 s
     run on the 27225-term (1+a+b+c)^8 (1+x+y+z)^8; 0.04 s this way). *)
  crs = Association /@ (CoefficientRules[#, S] & /@ cs);
  smons = Union @@ (Keys /@ crs);
  M = Lookup[#, smons, 0] & /@ crs;
  (* RREF of the transpose expresses each cs[[i]] in terms of pivot ones. *)
  R = DeleteCases[RowReduce@Transpose@M, {0 ..}];
  If[R === {}, Return[0, Module]];               (* rank 0: p == 0 *)
  piv = LengthWhile[#, # === 0 &] + 1 & /@ R;    (* pivot column indices *)
  B = cs[[piv]];                                 (* basis of Span[cs] *)
  Q = R . mons;                                  (* companion polys in T *)
  (Factor /@ B) . (Factor /@ Q)]

(* Anything outside the class above is a decline. *)
SeparablePolynomialDecomposition[_, _, _] := $Failed

End[];

SetAttributes[SeparablePolynomialDecomposition, Protected];

EndPackage[];
