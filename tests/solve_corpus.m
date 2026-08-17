(* solve_corpus.m
 *
 * Form-invariant stress corpus for Solve[].  A single List literal of
 * records
 *
 *     {label, equation(s), vars, domain-or-option, expectedCount, refNote}
 *
 * consumed held (unevaluated) by the C runner tests/test_solve_corpus.c
 * and by the standalone driver tests/solve_corpus_run.m.  Each case is
 * checked by BACK-SUBSTITUTION (see solve_check_prelude.m): the residual
 * lhs-rhs of every equation must vanish under each returned solution,
 * the solution count must equal expectedCount, and concrete values must
 * lie in the requested domain.  Verification is independent of solution
 * ordering and radical/formatting spelling.
 *
 * expectedCount is the MATHEMATICALLY CORRECT number of solution
 * rule-lists (cross-checked against Mathematica's Solve; sympy/maxima
 * agree on the sets).  Cases marked [Dn] are the empirically-confirmed
 * defects the corpus drives to zero — they FAIL/UNEVAL on the current
 * build and pass once their phase lands.
 *
 * Domain slot: Automatic (Complexes), Reals, Integers, Rationals, or an
 * option rule (Modulus->p, Cubics->True, VerifySolutions->True,
 * GeneratedParameters->K) forwarded as Solve's third argument.
 *
 * Lexer note: write coefficients on transcendental heads with a space
 * (`3 Exp[x]`, not `3Exp[x]`) or use E^x form — a digit fused to `Exp`
 * mis-lexes the digit as a float exponent.
 *)

{
  (* ===== A. Univariate polynomial (rational + radical roots) ===== *)
  {"A-quad-rat",        x^2 - 5 x + 6 == 0,        x, Automatic, 2, "MMA/sympy/maxima: {2,3}"},
  {"A-quad-mult2",      (x - 1)^2 == 0,            x, Automatic, 2, "double root, listed twice"},
  {"A-cubic-rat",       x^3 - 6 x^2 + 11 x - 6 == 0, x, Automatic, 3, "{1,2,3}"},
  {"A-biquadratic",     x^4 - 10 x^2 + 1 == 0,     x, Automatic, 4, "+-Sqrt[5+-2Sqrt6]"},
  {"A-cubic-binom",     x^3 - 2 == 0,              x, Automatic, 3, "1 real + 2 complex cube roots"},
  {"A-quintic-root",    x^5 - x - 1 == 0,          x, Automatic, 5, "Root[] objects, verified numerically"},
  {"A-quad-imag",       x^2 + 1 == 0,              x, Automatic, 2, "{-I,I}"},
  {"A-quad-nonmonic",   2 x^2 + 3 x + 1 == 0,      x, Automatic, 2, "{-1,-1/2}"},
  {"A-sextic-unity",    x^6 - 1 == 0,              x, Automatic, 6, "6th roots of unity"},
  {"A-cubic-factored",  x^3 + x^2 + x + 1 == 0,    x, Automatic, 3, "{-1,-I,I}"},
  {"A-quartic-unity",   x^4 - 1 == 0,              x, Automatic, 4, "{-1,1,-I,I}"},
  {"A-quad-radical",    x^2 - 2 == 0,              x, Automatic, 2, "+-Sqrt[2]"},

  (* ===== B. Domains (membership + count) ===== *)
  {"B-reals-none",      x^2 + 1 == 0,              x, Reals,    0, "no real roots"},
  {"B-reals-radical",   x^2 - 2 == 0,              x, Reals,    2, "+-Sqrt[2], Im=0"},
  {"B-int-none",        x^2 - 2 == 0,              x, Integers, 0, "no integer roots"},
  {"B-int-pair",        x^2 - 4 == 0,              x, Integers, 2, "{-2,2}"},
  {"B-reals-unity",     x^4 - 1 == 0,              x, Reals,    2, "{-1,1}"},
  {"B-rat-pair",        x^2 - 4 == 0,              x, Rationals, 2, "[D3] {-2,2}; unevaluated today"},
  {"B-rat-none",        x^2 - 2 == 0,              x, Rationals, 0, "[D3] no rational roots -> {}"},
  {"B-mod7-quad",       x^2 - 2 == 0,              x, Modulus -> 7, 2, "[D1] {3,4}; Modulus ignored today"},
  {"B-mod5-nonresidue", x^2 - 2 == 0,              x, Modulus -> 5, 0, "[D1] 2 is a non-residue mod 5 -> {}"},
  {"B-mod7-linear",     3 x - 1 == 0,              x, Modulus -> 7, 1, "[D1] {5}: 3*5=15=1 mod 7"},
  {"B-mod7-irred",      x^2 + x + 1 == 0,          x, Modulus -> 7, 2, "[D1] {2,4} mod 7"},
  {"B-reals-linear",    x + 1 == 0,                x, Reals,    1, "{-1}"},

  (* ===== C. Parametric / symbolic-coefficient (explicit rules only) ===== *)
  {"C-param-linear",    a x + b == 0,              x, Automatic, 1, "x->-b/a"},
  {"C-param-quad",      a x^2 + b x + c == 0,      x, Automatic, 2, "quadratic formula, verified symbolically"},
  {"C-param-sqrt",      x^2 == a,                  x, Automatic, 2, "x->+-Sqrt[a]"},
  {"C-param-factored",  (x - a) (x - b) == 0,      x, Automatic, 2, "{a,b}"},
  {"C-xy-hyperbola",    x y == 1,                  {x, y}, Automatic, 1, "[D4] x->1/y; nsdim today"},
  {"C-xy-circle",       x^2 + y^2 == 1,            {x, y}, Automatic, 2, "[D4] x->+-Sqrt[1-y^2]; nsdim today"},
  {"C-solve-one-param", x/y == 2,                  x, Automatic, 1, "x->2y"},

  (* ===== D. Transcendental / inverse-function ===== *)
  {"D-log1",            Log[x] == 1,               x, Automatic, 1, "x->E"},
  {"D-exp2",            E^x == 2,                   x, Automatic, 1, "periodic: Log[2]+2 I C[1] Pi"},
  {"D-exp1",            E^x == 1,                   x, Automatic, 1, "periodic: 2 I C[1] Pi"},
  {"D-sin0",            Sin[x] == 0,               x, Automatic, 2, "two periodic families"},
  {"D-sin-half",        Sin[x] == 1/2,             x, Automatic, 2, "two periodic families"},
  {"D-arcsin",          ArcSin[x] == Pi/6,         x, Automatic, 1, "x->1/2"},
  {"D-log-poly",        Log[x]^2 - 3 Log[x] + 2 == 0, x, Automatic, 2, "[D2] {E,E^2}; unevaluated today"},
  {"D-exp-poly",        E^(2 x) - 3 E^x + 2 == 0,  x, Automatic, 2, "[D2] u=E^x -> x=0,Log[2]; unevaluated today"},

  (* ===== E. Radical / nested-radical ===== *)
  {"E-sqrt-simple",     Sqrt[x] == 3,              x, Automatic, 1, "x->9"},
  {"E-sqrt-linear",     Sqrt[x + 5] == x - 1,      x, Automatic, 1, "x->4; x=-1 extraneous"},
  {"E-cuberoot",        x^(1/3) == 2,              x, Automatic, 1, "x->8"},
  {"E-sqrt-plus-x",     Sqrt[x] + x == 6,          x, Automatic, 1, "x->4; u=Sqrt[x] gives u=2 only"},
  {"E-sqrt-diff",       Sqrt[x] - Sqrt[x - 1] == 1, x, Automatic, 1, "x->1"},
  {"E-sqrt-quad",       Sqrt[x + 1] == x - 1,      x, Automatic, 1, "x->3; x=0 extraneous"},

  (* ===== F. Systems (zero-dimensional) ===== *)
  {"F-lin-2x2",         {x + y == 3, x - y == 1},          {x, y}, Automatic, 1, "{2,1}"},
  {"F-lin-3x3",         {x + y + z == 6, x - y == 1, y - z == 1}, {x, y, z}, Automatic, 1, "one solution"},
  {"F-circle-line",     {x^2 + y^2 == 25, x - y == 1},     {x, y}, Automatic, 2, "two intersections"},
  {"F-circle-line2",    {x^2 + y^2 == 1, x + y == 1},      {x, y}, Automatic, 2, "{1,0},{0,1}"},
  {"F-sum-prod",        {x y == 6, x + y == 5},            {x, y}, Automatic, 2, "{2,3},{3,2}"},
  {"F-inconsistent",    {x == 1, x == 2},                  x,      Automatic, 0, "{}"},
  {"F-dependent-lin",   {x + y == 2, 2 x + 2 y == 4},      {x, y}, Automatic, 1, "rank-1: x->2-y (svars)"},
  {"F-overdet-consist", {x + y == 3, x - y == 1, 2 x == 4}, {x, y}, Automatic, 1, "{2,1} satisfies all"},

  (* ===== G. Degenerate ===== *)
  {"G-tauto-zero",      0 == 0,                    x, Automatic, 1, "{{}}"},
  {"G-false-num",       1 == 0,                    x, Automatic, 0, "{}"},
  {"G-tauto-var",       x == x,                    x, Automatic, 1, "{{}}"},
  {"G-false-num2",      2 == 3,                    x, Automatic, 0, "{}"},
  {"G-identity",        x^2 == x^2,                x, Automatic, 1, "{{}}"},

  (* ===== H. Options / behaviour pinning ===== *)
  {"H-cubics-radical",  x^3 - 2 == 0,              x, Cubics -> True,  3, "radical form; verified numerically"},
  {"H-list-form",       {x^2 - 5 x + 6 == 0},      x, Automatic, 2, "single-equation list == bare equation"},
  {"H-genparam-K",      Sin[x] == 0,               x, GeneratedParameters -> K, 2, "families use K[k] not C[k]"},
  {"H-verify-radical",  Sqrt[x + 5] == x - 1,      x, VerifySolutions -> True, 1, "[D5] extraneous x=-1 dropped"},

  (* ===== I. Systems in unexpanded / product / power form ===== *)
  (* [D6] gb_from_expr's single-term parser cannot ingest Power[Plus,k] /
   * Times[Plus,...]; solvenlsys now runs internal_expand first (as the
   * GroebnerBasis builtin already did).  Unevaluated on builds before the
   * expand fix; correct after. *)
  {"I-two-circles",     {x^2 + y^2 == 1, (x - 1)^2 + y^2 == 1}, {x, y}, Automatic, 2, "[D6] {1/2,+-Sqrt3/2}; unexpanded (x-1)^2"},
  {"I-circles-radaxis", {x^2 + y^2 == 5, (x - 3)^2 + y^2 == 2}, {x, y}, Automatic, 2, "[D6] radical-axis intersection {2,+-1}"},
  {"I-product-line",    {(x + y) (x - y) == 0, x + y == 2},     {x, y}, Automatic, 1, "[D6] x+y=0 branch contradicts x+y=2 -> only {1,1}"},
  {"I-cubic-shift",     {(x + 1)^3 == y, y == x + 1},           {x, y}, Automatic, 3, "[D6] u=x+1: u^3=u -> x=-1,0,-2"},
  {"I-expanded-pin",    {x^2 + y^2 == 25, x - y == 1},          {x, y}, Automatic, 2, "regression pin: expanded form always worked"},

  (* ===== J. Rational-function systems (denominator clearing) ===== *)
  (* [D7] solvenlsys now puts each equation over a common denominator
   * (Together), solves Numerator == 0, and prunes any tuple on a non-constant
   * denominator's zero locus.  Unevaluated before the denominator-clearing
   * pass; correct after.  (Together cancels removable factors, matching
   * Solve's treatment of (x^2-1)/(x-1) as x+1.) *)
  {"J-reciprocal-sum",  {1/x + 1/y == 1, x + y == 4},          {x, y}, Automatic, 1, "[D7] xy=4 & x+y=4 -> {2,2}"},
  {"J-denom-prune",     {(x - y)/(x + y - 2) == 0, x y == 1},  {x, y}, Automatic, 1, "[D7] {1,1} zeroes x+y-2 -> pruned; only {-1,-1}"},
  {"J-all-denom-zero",  {1/(x - 1) == 1/(y - 1), x + y == 2},  {x, y}, Automatic, 0, "[D7] sole candidate x=y=1 is a pole -> {}"},
  {"J-const-denom-pin", {x/2 + y == 2, x - y == 1},            {x, y}, Automatic, 1, "constant denom (rational coeff) -> {2,1}; pin"},

  (* ===== K. Parametric nonlinear systems (coefficients in Q(params)) ===== *)
  (* [D8] free symbols that are not solve variables become coefficients in
   * Q(params); solvenlsys routes these through the RationalFunctions-coefficient
   * Gröbner engine and back-substitutes symbolically.  Unevaluated before the
   * parametric path (gb_from_expr rejected the free symbols); correct after. *)
  {"K-vieta",           {x + y == a, x y == b},                {x, y}, Automatic, 2, "[D8] Vieta: (a+-Sqrt[a^2-4b])/2"},
  {"K-quad-param",      {x^2 == a, x + y == 0},                {x, y}, Automatic, 2, "[D8] x->+-Sqrt[a], y->-/+Sqrt[a]"},
  {"K-circle-radius",   {x^2 + y^2 == r^2, y == 0},            {x, y}, Automatic, 2, "[D8] x->+-Sqrt[r^2], y->0"},

  (* ===== N. Systems over a finite field GF(p) (Modulus -> p, p prime) ===== *)
  (* [D9] Solve[{system}, vars, Modulus -> p] over GF(p): finite-field Gröbner
   * basis (gbmod.c) + triangular residue enumeration.  Refused before (systems
   * declined by solvemod.c); solved after.  Requires p prime. *)
  {"N-mod7-circ-diag",  {x^2 + y^2 == 1, x == y},              {x, y}, Modulus -> 7, 2, "[D9] {2,2},{5,5} mod 7"},
  {"N-mod7-sum-prod",   {x + y == 5, x y == 6},                {x, y}, Modulus -> 7, 2, "[D9] {2,3},{3,2} mod 7"},
  {"N-mod7-no-residue", {x^2 + y^2 == 1, (x - 1)^2 + y^2 == 1}, {x, y}, Modulus -> 7, 0, "[D9] y^2=6 non-residue mod 7 -> {}"},
  {"N-mod5-inconsist",  {x + y == 0, x + y == 1},              {x, y}, Modulus -> 5, 0, "[D9] unit ideal -> {}"},
  {"N-mod3-posdim",     {x + y == 1},                          {x, y}, Modulus -> 3, 3, "[D9] free var: all of GF(3) -> 3 pts"},
  {"N-mod5-three-var",  {x + y + z == 1, x y + y z + z x == 0, x y z == 0}, {x, y, z}, Modulus -> 5, 3, "[D9] perms of {0,0,1} mod 5"},

  (* ===== M. System regression pins (already correct; guard against drift) ===== *)
  {"M-cyclic3",         {x + y + z == 0, x y + y z + z x == -1, x y z == 0}, {x, y, z}, Automatic, 6, "six permutations of {-1,0,1}"},
  {"M-overdet-consist", {x^2 + y^2 == 1, x == 0, y == 1},       {x, y}, Automatic, 1, "overdetermined, consistent -> {0,1}"},
  {"M-inconsistent",    {x^2 + y^2 == 1, x^2 + y^2 == 2},       {x, y}, Automatic, 0, "unit ideal -> {}"},
  {"M-int-system",      {x^2 + y^2 == 25, x - y == 1},          {x, y}, Integers, 2, "integer lattice pts {4,3},{-3,-4}"},
  {"M-root-system",     {x^2 + y^2 == 1, x^3 + y == 1/2},       {x, y}, Automatic, 6, "Root[] tuples, verified numerically"}
}
