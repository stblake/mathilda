(* reduce_corpus.m
 *
 * Semantic stress corpus for Reduce[].  Each record is
 *     {label, expr, vars, dom, expected}
 * where `dom` is Automatic for the default (2-argument) form and `expected` is
 *   True | False        the statement fully decides (soundness cross-checked);
 *   "solved"            a formula EQUIVALENT to the input (grid + witnesses);
 *   "decline"           Reduce must bubble back unevaluated (the soundness
 *                       invariant: never invent a formula it cannot justify).
 *
 * The file is read as a HELD List literal by the C runner
 * (tests/test_reduce_corpus.c); each case is verified by back-sampling in
 * tests/reduce_check_prelude.m.  Every record is constructed to PASS today, so
 * the runner's baseline is 0 and any non-PASS is a real regression.  A genuine
 * EQUIVALENCE failure (not a mere non-minimal form) is a soundness bug.
 *)
{
  (* ---- equations over Complexes (default) : grid + witness ---- *)
  {"eq-quad-int",        x^2 == 4,                         x,       Automatic, "solved"},
  {"eq-quad-neg",        x^2 == -1,                        x,       Automatic, "solved"},
  {"eq-quad-surd",       x^2 == 2,                         x,       Automatic, "solved"},
  {"eq-quad-mid",        x^2 - 5 x + 6 == 0,               x,       Automatic, "solved"},
  {"eq-factored",        (x - 1)(x - 2) == 0,              x,       Automatic, "solved"},
  {"eq-linear-num",      2 x == 6,                         x,       Automatic, "solved"},
  {"eq-quartic-unity",   x^4 == 1,                         x,       Automatic, "solved"},
  {"eq-param-axb",       a x == b,                         x,       Automatic, "solved"},
  {"eq-param-axplusb",   a x + b == 0,                     x,       Automatic, "solved"},
  {"eq-param-ax0",       a x == 0,                         x,       Automatic, "solved"},

  (* ---- univariate real inequalities : grid-strong ---- *)
  {"re-gt1",             x^2 > 1,                          x,       Reals,     "solved"},
  {"re-ge1",             x^2 >= 1,                         x,       Reals,     "solved"},
  {"re-lt1",             x^2 < 1,                          x,       Reals,     "solved"},
  {"re-le1",             x^2 <= 1,                         x,       Reals,     "solved"},
  {"re-ne1",             x^2 != 1,                         x,       Reals,     "solved"},
  {"re-cube-pos",        x^3 > 0,                          x,       Reals,     "solved"},
  {"re-surd-bound",      x^2 - 2 > 0,                      x,       Reals,     "solved"},
  {"re-cubic-sign",      (x - 1)(x - 2)(x - 3) > 0,        x,       Reals,     "solved"},
  {"re-quad-factor",     x^2 - 5 x + 6 > 0,                x,       Reals,     "solved"},
  {"re-band",            x^2 > 1 && x < 5,                 x,       Reals,     "solved"},
  {"re-interval",        x > 0 && x < 1,                   x,       Reals,     "solved"},
  {"re-annulus",         x^2 > 1 && x^2 < 9,               x,       Reals,     "solved"},
  {"re-disj",            x^2 <= 1 || x >= 3,               x,       Reals,     "solved"},
  {"re-mixed-eq-ineq",   x^2 == 4 && x > 0,                x,       Reals,     "solved"},
  {"re-eq-surd",         x^2 - 2 == 0,                     x,       Reals,     "solved"},
  {"re-neq",             x != 0,                           x,       Reals,     "solved"},
  {"re-zero-point",      x^2 <= 0,                         x,       Reals,     "solved"},

  (* ---- univariate real : fully decided ---- *)
  {"re-true-sq",         x^2 >= 0,                         x,       Reals,     True},
  {"re-true-sq1",        x^2 + 1 > 0,                      x,       Reals,     True},
  {"re-false-sq",        x^2 < 0,                          x,       Reals,     False},
  {"re-false-sq1",       x^2 + 1 < 0,                      x,       Reals,     False},

  (* ---- linear real systems (Fourier-Motzkin) ---- *)
  {"fm-triangle",        x + y < 1 && x > 0 && y > 0,      {x, y},  Reals,     "solved"},
  {"fm-eq-ineq",         x + y == 1 && x > 0,              {x, y},  Reals,     "solved"},
  {"fm-eq-determined",   x + y == 1 && x - y == 3,         {x, y},  Reals,     "solved"},
  {"fm-rational",        2 x + 3 y <= 6 && x >= 0 && y >= 0, {x, y}, Reals,    "solved"},
  {"fm-infeasible",      x > 1 && x < 0 && y > 0,          {x, y},  Reals,     False},
  {"fm-disj",            x < 0 || x > 1,                   {x, y},  Reals,     "solved"},

  (* ---- parametric linear systems over Complexes ---- *)
  {"sys-determined",     x + y == 3 && x - y == 1,         {x, y},  Automatic, "solved"},
  {"sys-underdet",       x + y == 1,                       {x, y},  Automatic, "solved"},
  {"sys-param-cramer",   a x + y == 1 && x + y == 0,       {x, y},  Automatic, "solved"},
  {"sys-overdet",        a x == 1 && x == 2,               x,       Automatic, "solved"},
  {"sys-three-var",      x + y + z == 6 && x - y == 0 && z == 2, {x, y, z}, Automatic, "solved"},

  (* ---- Integers / Rationals ---- *)
  {"int-quad",           x^2 == 4,                         x,       Integers,  "solved"},
  {"int-quad-none",      x^2 == 2,                         x,       Integers,  False},
  {"int-lin-none",       2 x == 5,                         x,       Integers,  False},
  {"int-lin",            2 x + 1 == 5,                     x,       Integers,  "solved"},
  {"int-band",           x^2 <= 4,                         x,       Integers,  "solved"},
  {"int-band2",          x^2 < 5 && x > -3,                x,       Integers,  "solved"},
  {"int-range",          1 <= x <= 3,                      x,       Integers,  "solved"},
  {"int-bounded",        x^2 < 10 && x > 0,                x,       Integers,  "solved"},
  {"int-sys-bounded",    x + y == 5 && x > 0 && y > 0,     {x, y},  Integers,  "solved"},
  {"int-diophantine",    x + y == 3,                       {x, y},  Integers,  "solved"},
  {"int-diophantine2",   2 x + 3 y == 1,                   {x, y},  Integers,  "solved"},
  {"rat-quad",           x^2 == 4,                         x,       Rationals, "solved"},
  {"rat-quad-none",      x^2 == 2,                         x,       Rationals, False},

  (* ---- phase 6: multivariate nonlinear inequalities over the reals (CAD) ---- *)
  {"cad-disk-closed",    x^2 + y^2 <= 1,                   {x, y},  Reals,     "solved"},
  {"cad-disk-open",      x^2 + y^2 < 1,                    {x, y},  Reals,     "solved"},
  {"cad-disk-empty",     x^2 + y^2 < 0,                    {x, y},  Reals,     False},
  {"cad-plane",          x^2 + y^2 >= 0,                   {x, y},  Reals,     True},
  {"cad-product",        x y > 0,                          {x, y},  Reals,     "solved"},
  {"cad-lines",          (x - y) (x + y) > 0,              {x, y},  Reals,     "solved"},
  {"cad-ellipse",        x^2 + 4 y^2 < 4,                  {x, y},  Reals,     "solved"},
  {"cad-hyperbola",      x^2 - y^2 > 1,                    {x, y},  Reals,     "solved"},
  {"cad-half-disk",      x^2 + y^2 <= 1 && x >= 0,         {x, y},  Reals,     "solved"},
  {"cad-param-lin",      a x + y < 1,                      {x, y},  Reals,     "solved"},
  {"cad-absent-var",     x^2 > 1,                          {x, y},  Reals,     "solved"},
  {"cad-eq-curve",       x^2 + y^2 == 1,                   {x, y},  Reals,     "solved"},
  {"cad-punctured",      x^2 + y^2 <= 1 && x != 0,         {x, y},  Reals,     "solved"},
  {"cad-half-upper",     x^2 + y^2 <= 1 && y >= 0,         {x, y},  Reals,     "solved"},
  {"cad-point",          x^2 + y^2 <= 0,                   {x, y},  Reals,     "solved"},

  (* ---- soundness invariant: these MUST decline (never a guessed formula) ---- *)
  {"dec-default-ineq",   x^2 > 1,                          x,       Automatic, "decline"},
  {"dec-int-unbounded",  x > 0,                            x,       Integers,  "decline"},
  {"dec-int-unbounded2", x < 5,                            x,       Integers,  "decline"},
  {"dec-nonconst-denom", 1/x < 1,                          x,       Reals,     "decline"},
  {"dec-abs",            Abs[x] < 1,                       x,       Reals,     "decline"},
  {"dec-bad-domain",     x^2 == 4,                         x,       Booleans,  "decline"},
  {"dec-nl-multivar3",   x^2 + y^2 + z^2 < 1,              {x, y, z}, Reals,   "decline"},
  {"dec-nl-eq-complex",  x^2 + y^2 == 1,                   {x, y},  Automatic, "decline"},
  {"dec-complexes-neq",  x != 0,                           x,       Automatic, "decline"},
  {"dec-complexes-conj", x^2 == 4 && x^3 == 8,             x,       Automatic, "decline"},
  {"dec-nl-system",      x y == 1 && x + y == 3,           {x, y},  Automatic, "decline"}
}
