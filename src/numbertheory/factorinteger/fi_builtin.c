/*
 * fi_builtin.c — FactorInteger[...] and the FactorInteger`Method[...] entry
 * points: argument and option parsing, rational and Gaussian inputs, partial
 * factorisation, messages, output, registration and docstrings.
 *
 *   FactorInteger[n]                       complete factorisation
 *   FactorInteger[n, k]                    at most k entries; the last may be
 *                                          an unfactored remainder
 *   FactorInteger[n, Automatic]            same as FactorInteger[n]
 *   FactorInteger[n, Method -> m]          m a method name, Automatic, or
 *                                          {"Name", "Sub" -> v, ...}
 *   FactorInteger[n, GaussianIntegers -> True]   factor over Z[i]
 *   FactorInteger`Name[n, ...]             = FactorInteger[n, Method -> {"Name", ...}]
 *
 * n may be an integer, a rational (the denominator's primes get negative
 * exponents), or a Gaussian integer / Gaussian rational (factored over Z[i]
 * automatically, unit first, primes in the first quadrant sorted by real
 * then imaginary part — Mathematica's normalisation).
 */
#include "fi_internal.h"

#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "message.h"
#include "print.h"
#include "sym_names.h"
#include "symtab.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The parent docstring of Blake's rational base descent: deliberately long,
 * since the method is new (S. Blake, "Rational Base Descent: A Deterministic
 * Algorithm for Factoring Structured Semiprimes" and "Lattice-Augmented
 * Rational Base Descent"; both in factorisation_temp/). */
static const char RBD_DOC[] =
    "FactorInteger`BlakeRationalBaseDescent[n] factors n by Blake's rational base\n"
    "\tdescent (RBD), a special-purpose method for integers with a factor that is\n"
    "\tclose to a power of a rational number.\n"
    "FactorInteger`BlakeRationalBaseDescent[n, \"Sub\" -> v, ...] sets suboptions; the\n"
    "\tsame method is FactorInteger[n, Method -> {\"BlakeRationalBaseDescent\", ...}].\n"
    "\n"
    "Idea.  Suppose n = P M where the larger factor is P ~ c (a/b)^k for a base\n"
    "\ta/b > 1 (coprime a > b >= 1), a coefficient c and an exponent k.  Starting\n"
    "\tfrom Q = n, the descent repeats Q <- floor(b Q / a); each step divides by the\n"
    "\tbase and peels off one power.  After k steps Q has divided P down to about c,\n"
    "\tso Q lies within a small constant of the integer c M, and gcd(Q - j, n) for\n"
    "\ta small offset j exposes M (or M gcd(c, P)).  A single descent tests every\n"
    "\texponent k and every coefficient c at once, at a cost polynomial in log n.\n"
    "\n"
    "When it works.  The theorem: if (a/b)^k > M and P differs from floor(c (a/b)^k)\n"
    "\tby less than (a/b)^k / M, a factor is found at step k.  The structure must\n"
    "\tlie on the LARGER factor; the cofactor M is unconstrained apart from its\n"
    "\tsize.  The base must be known or found by search.  The success region in\n"
    "\tbase space is a comb of isolated teeth of relative width ~1/(k M) (with\n"
    "\tharmonic teeth at a/d for divisors d of a): nothing signals a nearby base,\n"
    "\tso bases are enumerated, not optimised.  For an integer base the descent\n"
    "\ttruncates base-a digits, and truncations inherit the algebraic divisors of\n"
    "\tCunningham numbers and repunits (truncation divisibility).  The method is no\n"
    "\tthreat to balanced RSA moduli: against an unstructured n the expected search\n"
    "\tis ~sqrt(n) bases (~n^(1/4) with the lattice detector).\n"
    "\n"
    "Variants (\"Variant\" -> name; each is also the head\n"
    "\tFactorInteger`BlakeRationalBaseDescent<name>):\n"
    "\t\"Exact\"      default.  The descent runs on S n with S = 2^s >= 2a/(a-b), so\n"
    "\t             the target is always one of four consecutive values, for any\n"
    "\t             base.  Nothing is tested above sqrt(cmax n), since no c M can\n"
    "\t             lie there; that depth is reached with one division.  Gcds are\n"
    "\t             batched: one per block of steps, replayed only on a hit.\n"
    "\t\"Classic\"    the published reference algorithm: Q <- floor(b Q/a), testing\n"
    "\t             gcd(Q - j, n) for |j| <= ceil(2 + a/(a-b)) at every depth.\n"
    "\t\"AllHits\"    runs the reference descent to the end, records every distinct\n"
    "\t             non-trivial gcd, and refines them into pairwise-coprime\n"
    "\t             divisors of n.  This is the mode for Cunningham numbers and\n"
    "\t             repunits with an integer base.\n"
    "\t\"Lattice\"    replaces the gcd window by Howgrave-Graham's approximate-divisor\n"
    "\t             detector (LLL on a d-dimensional lattice).  It finds c M from Q\n"
    "\t             within X <= n^(beta^2) when M >= n^beta, so P may be much\n"
    "\t             further from c (a/b)^k.  It runs only at depths X <= Q <=\n"
    "\t             sqrt(cmax n) + X, one LLL reduction per depth.\n"
    "\t\"TwoPass\"    Exact over the bases, then Lattice over a coarse grid of\n"
    "\t             \"BaseRange\": spacing X/(2 c k n^beta), each grid point realised\n"
    "\t             by the simplest rational near it.  Without a range it uses the\n"
    "\t             Stern-Brocot bases (height 16 unless \"Height\" is given).\n"
    "\t\"Coppersmith\" a known base, coefficient and exponent: one detector call at\n"
    "\t             floor(c (a/b)^k), which is factoring with high bits known.\n"
    "\n"
    "Suboptions:\n"
    "\t\"Base\" -> a/b or {a1/b1, ...}: the bases to use.  A base is reduced, and\n"
    "\t\tinverted if below 1.  Default: enumerate.\n"
    "\t\"Order\" -> \"Weighted\" (default: ascending a / log2(a/b), height times\n"
    "\t\tdescent length, so cheap simple bases go first), \"HeightDoubling\"\n"
    "\t\t(Stern-Brocot in height bands 16, 32, 64, ...), \"SternBrocot\"\n"
    "\t\t(breadth-first over the Stern-Brocot tree, coarse to fine) or \"Sum\"\n"
    "\t\t(a + b ascending).  All skip perfect-power bases (x/y)^m.\n"
    "\t\"Height\" -> H: largest numerator enumerated (default 1000).\n"
    "\t\"Bases\" -> B: enumeration budget (default 200000).\n"
    "\t\"MaxCoefficient\" -> cmax: largest coefficient c searched for (default\n"
    "\t\t2^64).  It sets the deferred-testing ceiling; Infinity tests every depth.\n"
    "\t\"Multiplier\" -> m: descend from m n, with gcds still taken with n.  This\n"
    "\t\thelps when m P, rather than P, is close to a power of the base.\n"
    "\t\"Tolerance\" -> X, \"Beta\" -> beta (default 1/2), \"Dimension\" -> d (default\n"
    "\t\t12) and \"Multiplicity\" -> t set the lattice detector.  When omitted, X is\n"
    "\t\tthe largest tolerance the LLL bound guarantees for d.\n"
    "\t\"BaseRange\" -> {lo, hi}: the TwoPass lattice grid.\n"
    "\t\"Coefficient\" -> c, \"Exponent\" -> k: required by \"Coppersmith\".\n"
    "\n"
    "Behaviour.  The suboptions describe the input.  A composite cofactor that the\n"
    "\tconfigured search cannot split is retried with a light default search\n"
    "\t(Exact, Stern-Brocot height 64), except under AllHits, whose pieces are\n"
    "\tre-run with AllHits instead.  Anything still composite is returned with\n"
    "\tFactorInteger::nofac.  Many-base searches, and a single base's lattice depths,\n"
    "\trun on all available threads.  The result does not depend on the thread\n"
    "\tcount: the lowest-indexed success wins.  Lattice reduction uses FLINT's LLL\n"
    "\twhen available and a native exact LLL otherwise.";

/* ------------------------------------------------------------------ */
/* Method table                                                        */
/* ------------------------------------------------------------------ */

enum { SUB_BOUND = 1, SUB_ITER = 2, SUB_STARTS = 4, SUB_B1 = 8, SUB_B2 = 16,
       SUB_CURVES = 32, SUB_SIGMA = 64, SUB_BASE = 128, SUB_RBD = 256 };

typedef struct {
    const char* name;       /* canonical Method string = FactorInteger`name */
    FiMethod    m;
    unsigned    subs;       /* accepted suboptions */
    const char* doc;        /* FactorInteger`name docstring */
} MethodInfo;

static const MethodInfo METHODS[] = {
    { "TrialDivision", FI_TRIAL, SUB_BOUND,
      "FactorInteger`TrialDivision[n] factors n by trial division only.\n"
      "\tPrimes are tested in word-sized batches (one bignum remainder per\n"
      "\tbatch).  Default bound: sqrt(n) when that is below 10^7, else 10^7;\n"
      "\ta cofactor left above the bound is returned with FactorInteger::nofac.\n"
      "\tSuboption: \"Bound\" -> B (at most 2^40)." },
    { "PollardRho", FI_RHO, SUB_ITER | SUB_STARTS,
      "FactorInteger`PollardRho[n] factors n with Brent's variant of Pollard rho\n"
      "\t(batched gcds, backtracking; Montgomery arithmetic on one and two words).\n"
      "\tFinds a factor p in about sqrt(p) steps.\n"
      "\tSuboptions: \"Iterations\" -> steps per start (default 2^24 on bignums,\n"
      "\t2^26 on words), \"Starts\" -> number of polynomials y^2 + c (default 4/16)." },
    { "Fermat", FI_FERMAT, SUB_ITER,
      "FactorInteger`Fermat[n] factors n by Fermat's method: it finds n = a^2 - b^2\n"
      "\tquickly when n has two factors close to sqrt(n).  Candidates are\n"
      "\tfiltered by quadratic-residue tables modulo nine small moduli.\n"
      "\tSuboption: \"Iterations\" -> number of a values tried (default 2^26)." },
    { "Lehman", FI_LEHMAN, SUB_ITER,
      "FactorInteger`Lehman[n] factors n by Lehman's O(n^(1/3)) method: trial\n"
      "\tdivision to n^(1/3), then Fermat on 4 k n for k <= n^(1/3).  Complete\n"
      "\tbelow 2^42; on larger n it runs within \"Iterations\" (default 2^24)." },
    { "ShanksSquareForms", FI_SQUFOF, SUB_ITER,
      "FactorInteger`ShanksSquareForms[n] factors n by Shanks' square forms\n"
      "\t(SQUFOF), trying the 16 Gower-Wagstaff multipliers in turn.  Expected\n"
      "\twork O(n^(1/4)); runs in word arithmetic for k n < 2^124 and declines\n"
      "\tbeyond.  Suboption: \"Iterations\" -> forward steps per multiplier." },
    { "CFRAC", FI_CFRAC, SUB_ITER,
      "FactorInteger`CFRAC[n] factors n by the continued-fraction method of\n"
      "\tMorrison and Brillhart: Knuth-Schroeppel multiplier, factor base sized\n"
      "\tby digit count, early-abort trial division, single large primes, and\n"
      "\tGF(2) elimination.  Suboption: \"Iterations\" -> continued-fraction\n"
      "\tsteps allowed (default 2^32)." },
    { "Dixon", FI_DIXON, SUB_ITER,
      "FactorInteger`Dixon[n] factors n by Dixon's random-squares method, with\n"
      "\tthe squares drawn near sqrt(k n) so the residues are square-root\n"
      "\tsized, early-abort trial division, single large primes and GF(2)\n"
      "\telimination.  Suboption: \"Iterations\" -> squares tried (default 2^30)." },
    { "QuadraticSieve", FI_QS, 0,
      "FactorInteger`QuadraticSieve[n] factors n by the self-initialising\n"
      "\tquadratic sieve, the fastest general method for 30-100 digit\n"
      "\tcomposites with no small factors.  Below 2^64 the word engine is used." },
    { "PollardPMinus1", FI_PM1, SUB_B1 | SUB_B2,
      "FactorInteger`PollardPMinus1[n] factors n by Pollard's p-1 method: it\n"
      "\tfinds a prime p of n when p - 1 is B1-smooth apart from one prime up to\n"
      "\tB2.  Uses GMP-ECM's fast stage 2 when available.\n"
      "\tSuboptions: \"B1\" (default 10^6), \"B2\" (default: GMP-ECM's choice)." },
    { "WilliamsPPlus1", FI_PP1, SUB_B1 | SUB_B2,
      "FactorInteger`WilliamsPPlus1[n] factors n by Williams' p+1 method: it\n"
      "\tfinds p when p + 1 is smooth (three seeds are run, since one seed\n"
      "\tworks only when its discriminant is a non-residue mod p).\n"
      "\tSuboptions: \"B1\" (default 10^6), \"B2\"." },
    { "ECM", FI_ECM, SUB_B1 | SUB_B2 | SUB_CURVES | SUB_SIGMA,
      "FactorInteger`ECM[n] factors n by Lenstra's elliptic-curve method (GMP-ECM),\n"
      "\twhose running time depends on the size of the factor found, not of n.\n"
      "\tBy default it climbs GMP-ECM's optimal levels (B1 = 2000, 11000, 5e4,\n"
      "\t2.5e5, 1e6, ...) with the expected curve count at each.\n"
      "\tSuboptions: \"B1\", \"B2\", \"Curves\", \"Sigma\" (first curve parameter)." },
    { "BlakeRationalBaseDescent", FI_RBD, SUB_BASE | SUB_RBD, RBD_DOC },
};
#define NMETHODS (sizeof METHODS / sizeof METHODS[0])

static const struct { const char* alias; FiMethod m; } ALIASES[] = {
    { "Automatic", FI_AUTOMATIC },
    { "PollardP-1", FI_PM1 }, { "PollardPM1", FI_PM1 }, { "P-1", FI_PM1 },
    { "WilliamsP+1", FI_PP1 }, { "WilliamsPP1", FI_PP1 }, { "P+1", FI_PP1 },
    { "SQUFOF", FI_SQUFOF },
    { "QS", FI_QS }, { "SIQS", FI_QS }, { "MPQS", FI_QS },
    { "ContinuedFraction", FI_CFRAC },
    { "EllipticCurve", FI_ECM },
};

/* Symbol aliases registered in addition to the canonical names. */
static const struct { const char* sym; FiMethod m; } SYM_ALIASES[] = {
    { "SQUFOF", FI_SQUFOF }, { "QS", FI_QS }, { "SIQS", FI_QS },
};

/* BlakeRationalBaseDescent variants: "Variant" -> name, and the alias heads
 * FactorInteger`BlakeRationalBaseDescent<name>. */
static const struct { const char* name; FiRbdVariant v; const char* doc; } RBD_VARIANTS[] = {
    { "Exact", FI_RBD_EXACT,
      "FactorInteger`BlakeRationalBaseDescentExact[n] is rational base descent by\n"
      "\tthe scaled exact descent: the descent runs on S n so that the target lies\n"
      "\tin a base-independent window of four values; nothing is tested above\n"
      "\tsqrt(\"MaxCoefficient\" n), which is reached by one division, and the gcds\n"
      "\tare batched.  Suboptions as for FactorInteger`BlakeRationalBaseDescent." },
    { "Classic", FI_RBD_CLASSIC,
      "FactorInteger`BlakeRationalBaseDescentClassic[n] is the reference rational\n"
      "\tbase descent: Q <- floor(b Q/a) from n, testing gcd(Q - j, n) for\n"
      "\t|j| <= ceil(2 + a/(a-b)) at every depth.  Suboptions as for\n"
      "\tFactorInteger`BlakeRationalBaseDescent." },
    { "AllHits", FI_RBD_ALLHITS,
      "FactorInteger`BlakeRationalBaseDescentAllHits[n] runs the reference descent\n"
      "\tto the end, records every non-trivial gcd, and refines them into pairwise\n"
      "\tcoprime divisors of n; this recovers the algebraic factors of Cunningham\n"
      "\tnumbers and repunits with an integer \"Base\" (truncation divisibility)." },
    { "Lattice", FI_RBD_LATTICE,
      "FactorInteger`BlakeRationalBaseDescentLattice[n] is lattice-augmented rational\n"
      "\tbase descent: Howgrave-Graham's approximate-divisor detector replaces the\n"
      "\tgcd window at each depth with Q in [X, sqrt(\"MaxCoefficient\" n) + X],\n"
      "\twidening the tolerance from O(a/(a-b)) to \"Tolerance\" X <= n^(beta^2).\n"
      "\tSuboptions \"Tolerance\", \"Beta\" (divisor >= n^beta, default 1/2),\n"
      "\t\"Multiplicity\", \"Dimension\" (default 12), plus the base options." },
    { "TwoPass", FI_RBD_TWOPASS,
      "FactorInteger`BlakeRationalBaseDescentTwoPass[n] runs the exact descent over\n"
      "\tthe bases, then the lattice descent over a coarse grid of \"BaseRange\" ->\n"
      "\t{lo, hi} (spacing X/(2 c k n^beta), each point realised by its simplest\n"
      "\trational), or over the same bases when no range is given." },
    { "Coppersmith", FI_RBD_COPPERSMITH,
      "FactorInteger`BlakeRationalBaseDescentCoppersmith[n, \"Base\" -> a/b,\n"
      "\t\"Coefficient\" -> c, \"Exponent\" -> k] finds a factor of n within the\n"
      "\tlattice tolerance of floor(c (a/b)^k) with one Howgrave-Graham detector\n"
      "\tcall (factoring with high bits known)." },
};
#define NRBD_VARIANTS (sizeof RBD_VARIANTS / sizeof RBD_VARIANTS[0])
static const char RBD_PREFIX[] = "BlakeRationalBaseDescent";

static const MethodInfo* method_info(FiMethod m) {
    for (size_t i = 0; i < NMETHODS; i++) if (METHODS[i].m == m) return &METHODS[i];
    return NULL;
}

/* Method name -> method; *variant is set (>= 0) for the RBD variant
 * aliases "BlakeRationalBaseDescent<Variant>", else left alone. */
static int lookup_method(const char* s, FiMethod* out, int* variant) {
    for (size_t i = 0; i < NMETHODS; i++)
        if (strcmp(METHODS[i].name, s) == 0) { *out = METHODS[i].m; return 1; }
    for (size_t i = 0; i < sizeof ALIASES / sizeof ALIASES[0]; i++)
        if (strcmp(ALIASES[i].alias, s) == 0) { *out = ALIASES[i].m; return 1; }
    size_t np = sizeof RBD_PREFIX - 1;
    if (strncmp(s, RBD_PREFIX, np) == 0)
        for (size_t i = 0; i < NRBD_VARIANTS; i++)
            if (strcmp(s + np, RBD_VARIANTS[i].name) == 0) {
                *out = FI_RBD;
                if (variant) *variant = (int)RBD_VARIANTS[i].v;
                return 1;
            }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Messages                                                            */
/* ------------------------------------------------------------------ */

/* One message per unfactored composite (a call can leave several). */
void fi_warn_incomplete(FiCtx* c, const mpz_t m) {
    if (c->warned) return;              /* internal callers: silent */
    char* s = mpz_get_str(NULL, 10, m);
    mth_message_gated(g_arith_warnings_muted, "FactorInteger", "nofac",
        "%s is composite but no factor was found within the search bounds; "
        "it is returned unfactored with exponent 1.", s ? s : "The argument");
    free(s);
}

/* ------------------------------------------------------------------ */
/* Expression helpers                                                  */
/* ------------------------------------------------------------------ */

static int is_head(const Expr* e, const char* sym, size_t argc) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym && e->data.function.arg_count == argc;
}

static int is_rule(const Expr* e) {
    return is_head(e, SYM_Rule, 2) || is_head(e, SYM_RuleDelayed, 2);
}

/* n as a rational p/q (q > 0).  Returns 0 if n is not an exact rational. */
static int get_rational(const Expr* e, mpz_t p, mpz_t q) {
    if (expr_is_integer_like(e)) { expr_to_mpz((Expr*)e, p); mpz_init_set_ui(q, 1); return 1; }
    if (is_head(e, SYM_Rational, 2) && expr_is_integer_like(e->data.function.args[0]) &&
        expr_is_integer_like(e->data.function.args[1])) {
        expr_to_mpz(e->data.function.args[0], p);
        expr_to_mpz(e->data.function.args[1], q);
        if (mpz_sgn(q) < 0) { mpz_neg(p, p); mpz_neg(q, q); }
        return 1;
    }
    return 0;
}

static Expr* mpz_expr(const mpz_t z) { return expr_bigint_normalize(expr_new_bigint_from_mpz(z)); }

static Expr* pair_expr(Expr* a, int64_t e) {
    Expr* args[2] = { a, expr_new_integer(e) };
    return expr_new_function(expr_new_symbol(SYM_List), args, 2);
}

/* ------------------------------------------------------------------ */
/* Option parsing                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    int gaussian;
    int have_k;
    int64_t k;
} CallOpts;

/* Read a numeric suboption value: positive integer or real. */
static int sub_number(const Expr* v, double* out) {
    if (v->type == EXPR_INTEGER) { *out = (double)v->data.integer; return v->data.integer > 0; }
    if (v->type == EXPR_REAL) { *out = v->data.real; return v->data.real > 0; }
    if (v->type == EXPR_BIGINT) { mpz_t z; expr_to_mpz((Expr*)v, z); *out = mpz_get_d(z); mpz_clear(z); return *out > 0; }
    return 0;
}

/* A positive exact integer into z. */
static int pos_integer(const Expr* v, mpz_t z) {
    if (!expr_is_integer_like(v)) return 0;
    expr_to_mpz((Expr*)v, z);
    return mpz_sgn(z) > 0;
}

/* A positive real number given exactly or approximately. */
static int pos_real(const Expr* v, double* out) {
    if (is_head(v, SYM_Rational, 2) && expr_is_integer_like(v->data.function.args[0]) &&
        expr_is_integer_like(v->data.function.args[1])) {
        mpq_t q;
        mpq_init(q);
        expr_to_mpz(v->data.function.args[0], mpq_numref(q));
        expr_to_mpz(v->data.function.args[1], mpq_denref(q));
        mpq_canonicalize(q);
        *out = mpq_get_d(q);
        mpq_clear(q);
        return *out > 0;
    }
    return sub_number(v, out);
}

/* A positive rational (exact, or a machine real converted exactly). */
static int pos_rational(const Expr* v, mpq_t q) {
    if (v->type == EXPR_REAL) { if (!(v->data.real > 0)) return 0; mpq_set_d(q, v->data.real); return 1; }
    mpz_t p, d;
    if (!get_rational(v, p, d)) return 0;
    mpz_set(mpq_numref(q), p);
    mpz_set(mpq_denref(q), d);
    mpq_canonicalize(q);
    mpz_clears(p, d, NULL);
    return mpq_sgn(q) > 0;
}

/* The BlakeRationalBaseDescent suboptions other than "Base".  Returns 1
 * when applied, 0 when k is one of them but the value is malformed, -1 when
 * k is not an RBD suboption. */
static int apply_rbd_sub(FiRbdOpts* o, const char* k, const Expr* val) {
    double x = 0;
    if (strcmp(k, "Variant") == 0) {
        if (val->type != EXPR_STRING) return 0;
        for (size_t i = 0; i < NRBD_VARIANTS; i++)
            if (strcmp(val->data.string, RBD_VARIANTS[i].name) == 0) { o->variant = (int)RBD_VARIANTS[i].v; return 1; }
        return 0;
    }
    if (strcmp(k, "Order") == 0) {
        if (val->type != EXPR_STRING) return 0;
        if (strcmp(val->data.string, "Weighted") == 0) { o->order = FI_RBD_ORDER_WEIGHTED; return 1; }
        if (strcmp(val->data.string, "HeightDoubling") == 0) { o->order = FI_RBD_ORDER_HD; return 1; }
        if (strcmp(val->data.string, "SternBrocot") == 0) { o->order = FI_RBD_ORDER_SB; return 1; }
        if (strcmp(val->data.string, "Sum") == 0) { o->order = FI_RBD_ORDER_SUM; return 1; }
        return 0;
    }
    if (strcmp(k, "MaxCoefficient") == 0) {
        if ((val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_Infinity) ||
            (is_head(val, SYM_DirectedInfinity, 1) && val->data.function.args[0]->type == EXPR_INTEGER &&
             val->data.function.args[0]->data.integer == 1)) {
            o->cmax_set = o->cmax_inf = 1;
            return 1;
        }
        if (!pos_integer(val, o->cmax)) return 0;
        o->cmax_set = 1;
        o->cmax_inf = 0;
        return 1;
    }
    if (strcmp(k, "Multiplier") == 0) return (o->mult_set = pos_integer(val, o->mult));
    if (strcmp(k, "Tolerance") == 0) return (o->tol_set = pos_integer(val, o->tol));
    if (strcmp(k, "Coefficient") == 0) return (o->coef_set = pos_integer(val, o->coef));
    if (strcmp(k, "Beta") == 0) {
        if (!pos_real(val, &x) || x > 1) return 0;
        o->beta = x;
        return 1;
    }
    if (strcmp(k, "BaseRange") == 0) {
        if (!is_head(val, SYM_List, 2)) return 0;
        if (!pos_rational(val->data.function.args[0], o->lo) ||
            !pos_rational(val->data.function.args[1], o->hi)) return 0;
        if (mpq_cmp_ui(o->lo, 1, 1) <= 0 || mpq_cmp(o->lo, o->hi) >= 0) return 0;
        o->range_set = 1;
        return 1;
    }
    int which = strcmp(k, "Height") == 0 ? 1 : strcmp(k, "Bases") == 0 ? 2 :
                strcmp(k, "Multiplicity") == 0 ? 3 : strcmp(k, "Dimension") == 0 ? 4 :
                strcmp(k, "Exponent") == 0 ? 5 : 0;
    if (!which) return -1;
    if (val->type != EXPR_INTEGER || val->data.integer < 1) return 0;
    int64_t v = val->data.integer;
    switch (which) {
    case 1: o->height = (uint64_t)v; break;
    case 2: o->max_bases = (uint64_t)v; break;
    case 3: if (v > 47) return 0; o->t = (int)v; break;
    case 4: if (v < 2 || v > 48) return 0; o->d = (int)v; break;
    default: o->expo = (uint64_t)v; break;
    }
    return 1;
}

/* Apply one method suboption ("Name" -> value).  Returns 0 (after a message)
 * on an unknown or malformed suboption. */
static int apply_sub(FiCtx* c, const Expr* key, const Expr* val) {
    const char* k = key->type == EXPR_STRING ? key->data.string :
                    key->type == EXPR_SYMBOL ? key->data.symbol.name : NULL;
    const MethodInfo* mi = method_info(c->method);
    unsigned allowed = mi ? mi->subs : 0;
    const char* mname = mi ? mi->name : "Automatic";
    double x = 0;
    if (!k) goto bad;
    if (strcmp(k, "Base") == 0 && (allowed & SUB_BASE)) {
        /* a/b, or a list of them; each must be a positive rational != 1. */
        const Expr* const* items = &val;
        size_t ni = 1;
        if (is_head(val, SYM_List, val->type == EXPR_FUNCTION ? val->data.function.arg_count : 0)) {
            items = (const Expr* const*)val->data.function.args;
            ni = val->data.function.arg_count;
            if (ni == 0) goto bad;
        }
        fi_rbd_opts_reset_bases(&c->rbd);
        for (size_t i = 0; i < ni; i++) {
            mpz_t p, q;
            if (!get_rational(items[i], p, q)) goto bad;
            int ok = mpz_sgn(p) > 0 && mpz_cmp(p, q) != 0;
            if (ok) fi_rbd_opts_add_base(&c->rbd, p, q);
            mpz_clears(p, q, NULL);
            if (!ok) goto bad;
        }
        return 1;
    }
    if (allowed & SUB_RBD) {
        int r = apply_rbd_sub(&c->rbd, k, val);
        if (r == 1) return 1;
        if (r == 0) goto bad;
    }
    if (!sub_number(val, &x)) goto bad;
    if      (strcmp(k, "Bound") == 0 && (allowed & SUB_BOUND))       c->bound = (uint64_t)x;
    else if (strcmp(k, "Iterations") == 0 && (allowed & SUB_ITER))   c->iterations = (uint64_t)x;
    else if (strcmp(k, "Starts") == 0 && (allowed & SUB_STARTS))     c->starts = (uint64_t)x;
    else if (strcmp(k, "Curves") == 0 && (allowed & SUB_CURVES))     c->curves = (uint64_t)x;
    else if (strcmp(k, "B1") == 0 && (allowed & SUB_B1))             c->B1 = x;
    else if (strcmp(k, "B2") == 0 && (allowed & SUB_B2))             c->B2 = x;
    else if (strcmp(k, "Sigma") == 0 && (allowed & SUB_SIGMA))       c->sigma = (uint64_t)x;
    else goto bad;
    return 1;
bad:
    mth_message("FactorInteger", "mopt",
        "%s is not a valid suboption setting for Method -> \"%s\".",
        k ? k : "The suboption", mname);
    return 0;
}

/* Method -> spec.  Returns 0 (after a message) when spec is not recognised. */
static int apply_method(FiCtx* c, const Expr* spec) {
    if (spec->type == EXPR_SYMBOL && spec->data.symbol.name == SYM_Automatic) {
        c->method = FI_AUTOMATIC;
        return 1;
    }
    if (spec->type == EXPR_STRING) {
        int variant = -1;
        if (lookup_method(spec->data.string, &c->method, &variant)) {
            if (variant >= 0) c->rbd.variant = variant;
            return 1;
        }
        mth_message("FactorInteger", "method",
            "Method -> \"%s\" is not a known factoring method; see ?FactorInteger.",
            spec->data.string);
        return 0;
    }
    if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL &&
        spec->data.function.head->data.symbol.name == SYM_List && spec->data.function.arg_count >= 1) {
        if (!apply_method(c, spec->data.function.args[0])) return 0;
        for (size_t i = 1; i < spec->data.function.arg_count; i++) {
            const Expr* r = spec->data.function.args[i];
            if (!is_rule(r) || !apply_sub(c, r->data.function.args[0], r->data.function.args[1]))
                return 0;
        }
        return 1;
    }
    mth_message("FactorInteger", "method",
        "The Method setting is not a method name, Automatic, or {name, suboptions...}.");
    return 0;
}

/* Parse args[first..] of FactorInteger (or of a FactorInteger`X call when
 * `named` is set: bare suboption rules allowed, Method not).  Returns 0 to
 * leave the call unevaluated. */
static int parse_args(const Expr* res, size_t first, FiCtx* c, CallOpts* o, int named) {
    size_t argc = res->data.function.arg_count;
    for (size_t i = first; i < argc; i++) {
        const Expr* a = res->data.function.args[i];
        if (a->type == EXPR_INTEGER && i == first) {
            if (a->data.integer < 1) return 0;          /* FactorInteger[n, 0]: unevaluated */
            o->have_k = 1;
            o->k = a->data.integer;
            continue;
        }
        /* A positive BigInt k (10^30) asks for more factors than any n has:
         * the complete factorisation, as Infinity does. */
        if (a->type == EXPR_BIGINT && i == first && mpz_sgn(a->data.bigint) > 0) continue;
        /* Automatic, All and Infinity as k: the complete factorisation. */
        if (a->type == EXPR_SYMBOL && i == first &&
            (a->data.symbol.name == SYM_Automatic || a->data.symbol.name == SYM_All ||
             a->data.symbol.name == SYM_Infinity)) continue;
        if (!is_rule(a)) {
            /* A non-rule past the k slot: FactorInteger::nonopt.
             * A non-integer k itself stays unevaluated silently. */
            if (i > first || (i == first && first == 1 && argc > 2)) {
                char* bad = expr_to_string((Expr*)a);
                char* call = expr_to_string((Expr*)res);
                mth_message("FactorInteger", "nonopt",
                    "Options expected (instead of %s) beyond position %zu in %s. "
                    "An option must be a rule or a list of rules.",
                    bad ? bad : "?", i, call ? call : "?");
                free(bad); free(call);
            }
            return 0;
        }
        const Expr* lhs = a->data.function.args[0];
        const Expr* rhs = a->data.function.args[1];
        if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == SYM_GaussianIntegers) {
            if (rhs->type != EXPR_SYMBOL) return 0;
            if (rhs->data.symbol.name == SYM_True) o->gaussian = 1;
            else if (rhs->data.symbol.name == SYM_False) o->gaussian = 0;
            else return 0;
            continue;
        }
        if (!named && lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == SYM_Method) {
            if (!apply_method(c, rhs)) return 0;
            continue;
        }
        if (named && lhs->type == EXPR_STRING) {
            if (!apply_sub(c, lhs, rhs)) return 0;
            continue;
        }
        if (lhs->type == EXPR_SYMBOL) {
            mth_message("FactorInteger", "optx", "Unknown option %s in FactorInteger.",
                        lhs->data.symbol.name);
        }
        return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* Rational factorisation                                              */
/* ------------------------------------------------------------------ */

/* Factor p/q into out (q's primes negative). */
static void factor_rational(const mpz_t p, const mpz_t q, FiFactors* out, FiCtx* c) {
    mpz_t a;
    mpz_init(a);
    mpz_abs(a, p);
    fi_factor(a, out, c);
    if (mpz_cmp_ui(q, 1) > 0) {
        FiFactors den;
        fi_factors_init(&den);
        int64_t k = c->k_limit;
        c->k_limit = 0;                 /* partial applies to the whole answer */
        fi_factor(q, &den, c);
        c->k_limit = k;
        for (size_t i = 0; i < den.n; i++) fi_factors_add(out, den.v[i].p, -den.v[i].e);
        fi_factors_clear(&den);
    }
    mpz_clear(a);
}

/* Build the output list for a rational input. */
static Expr* rational_result(const mpz_t p, const mpz_t q, FiFactors* f) {
    fi_factors_sort(f);
    int neg = mpz_sgn(p) < 0;
    size_t total = f->n + (size_t)neg;
    Expr** items = malloc((total + 1) * sizeof(Expr*));
    size_t ni = 0;
    if (neg) items[ni++] = pair_expr(expr_new_integer(-1), 1);
    for (size_t i = 0; i < f->n; i++) items[ni++] = pair_expr(mpz_expr(f->v[i].p), f->v[i].e);
    (void)q;
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, ni);
    free(items);
    return out;
}

/* ---- FactorInteger[n, k]: the partial-factorisation form ------------------ *
 * Fixed by some 250 verified cases, numerator and denominator each on
 * their own:
 *   1. Trial division by the primes below 2^16 (the bound: 65521 is
 *      found, 65537 is not) gives the small prime powers T, ascending, and a
 *      cofactor C.
 *   2. C = 1: a part with at most k primes is listed in full; otherwise its
 *      k - 1 LARGEST prime powers are kept and the rest multiplied into one
 *      entry {rest, 1} -- FactorInteger[360, 2] is {{5, 1}, {72, 1}}.
 *   3. C > 1 and |T| > k - 1: C stays one unsplit entry, the k - 2 largest of
 *      T are kept, the rest of T merged -- FactorInteger[10^30 + 1, 2] is
 *      {{165573604901641, 1}, {6039610000603961, 1}}: the five trial primes
 *      merged, the two large ones left as their product.
 *   4. C > 1 with room r = k - |T| > 1: C is split into at most r entries,
 *      smallest factors first, the remainder left unsplit. Which factors are
 *      split off first is a matter of search order, so another splitting
 *      would be equally valid (C = 1000003 1000033 1000037 with r = 2 is
 *      {1000003, 1000033 1000037}; {1000033, 1000003 1000037} is as good).
 *   5. A negative n has no separate {-1, 1}; the sign goes into the
 *      denominator's part if there is one, else the numerator's, onto: the
 *      merged trial entry if any; else an unsplit C; else (C split into primes)
 *      the smallest entry; else (C = 1) the largest prime power -- as a base
 *      with exponent 1 (FactorInteger[-36, 2] is {{-9, 1}, {2, 2}}).
 *   6. Entries sorted by base value; k = 1 is {{n, 1}}; k = All / Infinity /
 *      Automatic is the complete factorisation.
 * Steps 1 and 3 mean a large hard cofactor is never factored when k is small. */
typedef struct { mpq_t base; int64_t e; } FiPEntry;

static int fi_pentry_cmp(const void* a, const void* b) {
    return mpq_cmp(((const FiPEntry*)a)->base, ((const FiPEntry*)b)->base);
}

#define FI_PARTIAL_TRIAL_BOUND 65535u

/* A growable entry list. */
typedef struct { FiPEntry* v; size_t n, cap; } FiPList;

static void fi_pentry_push_z(FiPList* L, const mpz_t b, int64_t e) {
    if (L->n == L->cap) {
        L->cap = L->cap ? 2 * L->cap : 16;
        L->v = realloc(L->v, L->cap * sizeof(FiPEntry));
    }
    mpq_init(L->v[L->n].base);
    mpq_set_z(L->v[L->n].base, b);
    L->v[L->n].e = e;
    L->n++;
}

/* Append the entries of one part (m > 1; dir = +1 numerator, -1 denominator)
 * to L. Returns the index of the entry that carries a sign (rule 5). */
static size_t fi_partial_part(const mpz_t m_in, int dir, int64_t k, FiCtx* c, FiPList* L) {
    mpz_t m, rest, pw;
    mpz_inits(m, rest, pw, NULL);
    mpz_set(m, m_in);
    FiFactors T;
    fi_factors_init(&T);
    fi_trial_divide(m, FI_PARTIAL_TRIAL_BOUND, 1, &T);
    fi_factors_sort(&T);
    size_t t = T.n;
    size_t carrier;

    /* merge T[0 .. upto) into one {rest, dir} entry */
    #define FI_MERGE_TRIAL(upto) do {                                         \
        mpz_set_ui(rest, 1);                                                  \
        for (size_t i_ = 0; i_ < (upto); i_++) {                              \
            mpz_pow_ui(pw, T.v[i_].p, (unsigned long)T.v[i_].e);              \
            mpz_mul(rest, rest, pw);                                          \
        }                                                                     \
        fi_pentry_push_z(L, rest, dir);                               \
    } while (0)

    if (mpz_cmp_ui(m, 1) == 0) {                                   /* rule 2 */
        if (t <= (size_t)k) {
            for (size_t i = 0; i < t; i++) fi_pentry_push_z(L, T.v[i].p, dir * T.v[i].e);
            carrier = L->n - 1;                                    /* largest */
        } else {
            size_t keep_from = t - ((size_t)k - 1);
            carrier = L->n;
            FI_MERGE_TRIAL(keep_from);
            for (size_t i = keep_from; i < t; i++) fi_pentry_push_z(L, T.v[i].p, dir * T.v[i].e);
        }
    } else if (t > (size_t)k - 1) {                                /* rule 3 */
        size_t keep_from = t - ((size_t)k - 2);
        carrier = L->n;
        FI_MERGE_TRIAL(keep_from);
        for (size_t i = keep_from; i < t; i++) fi_pentry_push_z(L, T.v[i].p, dir * T.v[i].e);
        fi_pentry_push_z(L, m, dir);
    } else {                                                       /* rule 4 */
        size_t first = L->n;
        for (size_t i = 0; i < t; i++) fi_pentry_push_z(L, T.v[i].p, dir * T.v[i].e);
        int64_t room = k - (int64_t)t;
        if (room <= 1 || fi_is_prime(m)) {
            carrier = L->n;
            fi_pentry_push_z(L, m, dir);
        } else {
            FiFactors cf;
            fi_factors_init(&cf);
            c->k_limit = room;                    /* smallest first, remainder unsplit */
            fi_factor(m, &cf, c);
            c->k_limit = 0;
            fi_factors_sort(&cf);
            /* The early stop may overshoot by one (a last prime is emitted as
             * itself): keep the room - 1 smallest, multiply the rest into one
             * unsplit entry. */
            size_t keep = cf.n;
            if (cf.n > (size_t)room) {
                keep = (size_t)room - 1;
                mpz_set_ui(rest, 1);
                for (size_t i = keep; i < cf.n; i++) {
                    mpz_pow_ui(pw, cf.v[i].p, (unsigned long)cf.v[i].e);
                    mpz_mul(rest, rest, pw);
                }
            }
            carrier = (size_t)-1;
            for (size_t i = 0; i < keep; i++) {
                if (carrier == (size_t)-1 && !fi_is_prime(cf.v[i].p)) carrier = L->n;
                fi_pentry_push_z(L, cf.v[i].p, dir * cf.v[i].e);
            }
            if (keep < cf.n) {
                carrier = L->n;
                fi_pentry_push_z(L, rest, dir);
            }
            if (carrier == (size_t)-1) {                 /* all prime: the smallest */
                carrier = first;
                for (size_t i = first + 1; i < L->n; i++)
                    if (mpq_cmp(L->v[i].base, L->v[carrier].base) < 0) carrier = i;
            }
            fi_factors_clear(&cf);
        }
    }
    #undef FI_MERGE_TRIAL
    fi_factors_clear(&T);
    mpz_clears(m, rest, pw, NULL);
    return carrier;
}

static Expr* partial_result(const mpz_t p, const mpz_t q, int64_t k, FiCtx* c) {
    if (k == 1) {                                         /* {{n, 1}} */
        Expr* items[1] = { pair_expr(make_rational_mpz(p, q), 1) };
        return expr_new_function(expr_new_symbol(SYM_List), items, 1);
    }
    mpz_t a;
    mpz_init(a);
    mpz_abs(a, p);
    FiPList L = { NULL, 0, 0 };
    c->k_limit = 0;
    size_t s_num = (mpz_cmp_ui(a, 1) > 0) ? fi_partial_part(a, 1, k, c, &L) : (size_t)-1;
    size_t s_den = (mpz_cmp_ui(q, 1) > 0) ? fi_partial_part(q, -1, k, c, &L) : (size_t)-1;
    FiPEntry* ent = L.v;
    size_t n_ent = L.n;

    if (mpz_sgn(p) < 0) {
        size_t at = (s_den != (size_t)-1) ? s_den : s_num;
        if (at != (size_t)-1) {
            /* base^e -> (-(base^|e|))^(sign of e) as a single base, exponent 1 */
            FiPEntry* t = &ent[at];
            mpz_t pw;
            mpz_init(pw);
            int64_t e = t->e < 0 ? -t->e : t->e;
            mpz_pow_ui(pw, mpq_numref(t->base), (unsigned long)e);
            mpq_set_z(t->base, pw);
            if (t->e < 0) mpq_inv(t->base, t->base);
            mpq_neg(t->base, t->base);
            t->e = 1;
            mpz_clear(pw);
        }
    }
    qsort(ent, n_ent, sizeof(FiPEntry), fi_pentry_cmp);
    Expr** items = malloc((n_ent + 1) * sizeof(Expr*));
    for (size_t i = 0; i < n_ent; i++) {
        items[i] = pair_expr(make_rational_mpz(mpq_numref(ent[i].base), mpq_denref(ent[i].base)), ent[i].e);
        mpq_clear(ent[i].base);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, n_ent);
    free(items);
    free(ent);
    mpz_clear(a);
    return out;
}

/* ------------------------------------------------------------------ */
/* Gaussian integers                                                   */
/* ------------------------------------------------------------------ */

typedef struct { mpz_t re, im; int64_t e; } GEntry;
typedef struct { GEntry* v; size_t n, cap; } GList;

static void gl_add(GList* g, const mpz_t re, const mpz_t im, int64_t e) {
    for (size_t i = 0; i < g->n; i++)
        if (mpz_cmp(g->v[i].re, re) == 0 && mpz_cmp(g->v[i].im, im) == 0) { g->v[i].e += e; return; }
    if (g->n == g->cap) { g->cap = g->cap ? 2 * g->cap : 8; g->v = realloc(g->v, g->cap * sizeof(GEntry)); }
    mpz_init_set(g->v[g->n].re, re);
    mpz_init_set(g->v[g->n].im, im);
    g->v[g->n].e = e;
    g->n++;
}
static void gl_clear(GList* g) {
    for (size_t i = 0; i < g->n; i++) mpz_clears(g->v[i].re, g->v[i].im, NULL);
    free(g->v);
}

/* Rotate (re, im) by powers of i into the first quadrant: re > 0, im >= 0. */
static void normalize_quadrant(mpz_t re, mpz_t im) {
    mpz_t t;
    mpz_init(t);
    for (int k = 0; k < 4 && !(mpz_sgn(re) > 0 && mpz_sgn(im) >= 0); k++) {
        mpz_neg(t, im); mpz_set(im, re); mpz_set(re, t);   /* times i */
    }
    mpz_clear(t);
}

/* p == 1 (mod 4) prime: write p = c^2 + d^2 (Cornacchia / Hermite–Serret). */
static void two_squares(const mpz_t p, mpz_t c, mpz_t d) {
    mpz_t t, x, e, a, b, r;
    mpz_inits(t, x, e, a, b, r, NULL);
    for (unsigned long tv = 2;; tv++) {
        mpz_set_ui(t, tv);
        if (mpz_legendre(t, p) == -1) break;
    }
    mpz_sub_ui(e, p, 1);
    mpz_fdiv_q_2exp(e, e, 2);
    mpz_powm(x, t, e, p);                   /* x^2 == -1 (mod p) */
    mpz_set(a, p);
    mpz_set(b, x);
    mpz_mul(r, b, b);
    while (mpz_cmp(r, p) > 0) {
        mpz_mod(r, a, b);
        mpz_set(a, b);
        mpz_set(b, r);
        mpz_mul(r, b, b);
    }
    mpz_set(c, b);
    mpz_sub(r, p, r);
    mpz_sqrt(d, r);
    mpz_clears(t, x, e, a, b, r, NULL);
}

/* If (u + v i) divides (a + b i), replace (a, b) by the quotient. */
static int gdiv(mpz_t a, mpz_t b, const mpz_t u, const mpz_t v) {
    mpz_t N, re, im, t;
    mpz_inits(N, re, im, t, NULL);
    mpz_mul(N, u, u); mpz_addmul(N, v, v);
    mpz_mul(re, a, u); mpz_addmul(re, b, v);
    mpz_mul(im, b, u); mpz_submul(im, a, v);
    int ok = mpz_divisible_p(re, N) && mpz_divisible_p(im, N);
    if (ok) { mpz_divexact(a, re, N); mpz_divexact(b, im, N); }
    mpz_clears(N, re, im, t, NULL);
    return ok;
}

/* Factor the nonzero Gaussian integer a + b i into first-quadrant primes;
 * (ua, ub) receives the unit.  Returns 0 if the norm could not be fully
 * factored (the call is then left unevaluated). */
static int gaussian_factor(const mpz_t a0, const mpz_t b0, GList* out, mpz_t ua, mpz_t ub, FiCtx* c) {
    mpz_t N, a, b, cc, dd, cr, ci;
    mpz_inits(N, a, b, cc, dd, cr, ci, NULL);
    mpz_set(a, a0);
    mpz_set(b, b0);
    mpz_mul(N, a, a);
    mpz_addmul(N, b, b);
    FiFactors nf;
    fi_factors_init(&nf);
    int64_t kl = c->k_limit;
    c->k_limit = 0;
    fi_factor(N, &nf, c);
    c->k_limit = kl;
    int ok = !c->incomplete;
    for (size_t i = 0; i < nf.n && ok; i++) {
        const mpz_t* p = (const mpz_t*)&nf.v[i].p;
        int64_t e = nf.v[i].e;
        if (mpz_cmp_ui(*p, 2) == 0) {
            mpz_set_ui(cr, 1); mpz_set_ui(ci, 1);
            for (int64_t k = 0; k < e; k++) gdiv(a, b, cr, ci);
            gl_add(out, cr, ci, e);
        } else if (mpz_fdiv_ui(*p, 4) == 3) {
            mpz_set(cr, *p); mpz_set_ui(ci, 0);
            for (int64_t k = 0; k < e / 2; k++) gdiv(a, b, cr, ci);
            gl_add(out, cr, ci, e / 2);
        } else {
            two_squares(*p, cc, dd);
            /* pi = c + d i and its conjugate, each taken in its first-quadrant
             * form BEFORE dividing, so the leftover unit is relative to the
             * primes actually reported. */
            int64_t ep = 0, eq = 0;
            mpz_set(cr, cc); mpz_set(ci, dd);
            normalize_quadrant(cr, ci);
            while (ep < e && gdiv(a, b, cr, ci)) ep++;
            if (ep) gl_add(out, cr, ci, ep);
            mpz_set(cr, cc); mpz_neg(ci, dd);
            normalize_quadrant(cr, ci);
            while (eq < e - ep && gdiv(a, b, cr, ci)) eq++;
            if (eq) gl_add(out, cr, ci, eq);
        }
    }
    mpz_set(ua, a);                         /* what remains is a unit */
    mpz_set(ub, b);
    fi_factors_clear(&nf);
    mpz_clears(N, a, b, cc, dd, cr, ci, NULL);
    return ok;
}

static int gcmp(const void* x, const void* y) {
    const GEntry* p = x;
    const GEntry* q = y;
    int c = mpz_cmp(p->re, q->re);
    return c ? c : mpz_cmp(p->im, q->im);
}

static Expr* gaussian_expr(const mpz_t re, const mpz_t im) {
    if (mpz_sgn(im) == 0) return mpz_expr(re);
    return make_complex(mpz_expr(re), mpz_expr(im));
}

/* Factor z = (pr/qr) + (pi/qi) i over Z[i]. */
static Expr* gaussian_result(const mpz_t pr, const mpz_t qr, const mpz_t pi, const mpz_t qi, FiCtx* c) {
    mpz_t m, a, b, ua, ub, va, vb, t;
    mpz_inits(m, a, b, ua, ub, va, vb, t, NULL);
    mpz_lcm(m, qr, qi);
    mpz_divexact(t, m, qr); mpz_mul(a, pr, t);
    mpz_divexact(t, m, qi); mpz_mul(b, pi, t);
    Expr* out = NULL;
    if (mpz_sgn(a) == 0 && mpz_sgn(b) == 0) {
        Expr* items[1] = { pair_expr(expr_new_integer(0), 1) };
        out = expr_new_function(expr_new_symbol(SYM_List), items, 1);
        goto done;
    }
    GList g = { NULL, 0, 0 };
    int ok = gaussian_factor(a, b, &g, ua, ub, c);
    mpz_set_ui(va, 1); mpz_set_ui(vb, 0);
    if (ok && mpz_cmp_ui(m, 1) > 0) {
        GList h = { NULL, 0, 0 };
        mpz_set_ui(t, 0);
        ok = gaussian_factor(m, t, &h, va, vb, c);
        for (size_t i = 0; i < h.n; i++) gl_add(&g, h.v[i].re, h.v[i].im, -h.v[i].e);
        gl_clear(&h);
    }
    if (ok) {
        /* unit = u / v = u * conj(v) for units. */
        mpz_t wr, wi;
        mpz_inits(wr, wi, NULL);
        mpz_mul(wr, ua, va); mpz_addmul(wr, ub, vb);
        mpz_mul(wi, ub, va); mpz_submul(wi, ua, vb);
        size_t kept = 0;
        for (size_t i = 0; i < g.n; i++) if (g.v[i].e) g.v[kept++] = g.v[i]; else mpz_clears(g.v[i].re, g.v[i].im, NULL);
        g.n = kept;
        qsort(g.v, g.n, sizeof(GEntry), gcmp);
        int unit_one = mpz_cmp_ui(wr, 1) == 0 && mpz_sgn(wi) == 0;
        size_t ni = 0;
        Expr** items = malloc((g.n + 1) * sizeof(Expr*));
        if (!unit_one || g.n == 0) items[ni++] = pair_expr(gaussian_expr(wr, wi), 1);
        for (size_t i = 0; i < g.n; i++) items[ni++] = pair_expr(gaussian_expr(g.v[i].re, g.v[i].im), g.v[i].e);
        out = expr_new_function(expr_new_symbol(SYM_List), items, ni);
        free(items);
        mpz_clears(wr, wi, NULL);
    }
    gl_clear(&g);
done:
    mpz_clears(m, a, b, ua, ub, va, vb, t, NULL);
    return out;
}

/* ------------------------------------------------------------------ */
/* Entry points                                                        */
/* ------------------------------------------------------------------ */

static Expr* factorinteger_core(Expr* res, FiCtx* c, int named) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1) {
        mth_message("FactorInteger", "argt",
                    "FactorInteger called with 0 arguments; 1 or 2 arguments are expected.");
        return NULL;
    }
    Expr* n = res->data.function.args[0];
    CallOpts o = { 0, 0, 0 };
    if (!parse_args(res, 1, c, &o, named)) return NULL;
    c->k_limit = o.have_k ? o.k : 0;
    if (c->method == FI_RBD && c->rbd.variant == FI_RBD_COPPERSMITH &&
        (c->rbd.nbases != 1 || !c->rbd.coef_set || c->rbd.expo == 0)) {
        mth_message("FactorInteger", "rbdcop",
            "The Coppersmith variant of BlakeRationalBaseDescent needs one \"Base\", "
            "a \"Coefficient\" and an \"Exponent\".");
        return NULL;
    }

    mpz_t p, q;
    if (get_rational(n, p, q)) {
        Expr* out;
        if (o.gaussian && o.have_k && o.k == 1) {
            Expr* items[1] = { pair_expr(make_rational_mpz(p, q), 1) };
            out = expr_new_function(expr_new_symbol(SYM_List), items, 1);
        } else if (o.gaussian) {
            mpz_t zi, oi;
            mpz_init_set_ui(zi, 0);
            mpz_init_set_ui(oi, 1);
            out = gaussian_result(p, q, zi, oi, c);
            mpz_clears(zi, oi, NULL);
        } else if (mpz_sgn(p) == 0) {
            Expr* items[1] = { pair_expr(expr_new_integer(0), 1) };
            out = expr_new_function(expr_new_symbol(SYM_List), items, 1);
        } else if (mpz_cmpabs_ui(p, 1) == 0 && mpz_cmp_ui(q, 1) == 0) {
            Expr* items[1] = { pair_expr(mpz_expr(p), 1) };
            out = expr_new_function(expr_new_symbol(SYM_List), items, 1);
        } else if (o.have_k) {
            out = partial_result(p, q, o.k, c);
        } else {
            FiFactors f;
            fi_factors_init(&f);
            factor_rational(p, q, &f, c);
            out = rational_result(p, q, &f);
            fi_factors_clear(&f);
        }
        mpz_clears(p, q, NULL);
        return out;
    }
    Expr *re, *im;
    if (is_complex(n, &re, &im) && o.have_k && o.k == 1) {
        /* k = 1 is {{n, 1}} in Z[i] too. */
        Expr* items[1] = { pair_expr(expr_copy(n), 1) };
        return expr_new_function(expr_new_symbol(SYM_List), items, 1);
    }
    if (is_complex(n, &re, &im)) {
        mpz_t pr, qr, pi, qi;
        if (!get_rational(re, pr, qr)) goto inexact;
        if (!get_rational(im, pi, qi)) { mpz_clears(pr, qr, NULL); goto inexact; }
        Expr* out = gaussian_result(pr, qr, pi, qi, c);
        mpz_clears(pr, qr, pi, qi, NULL);
        return out;
    }
inexact:
    /* Mathematica: FactorInteger::exact for an inexact or symbolic argument
     * (lists thread through Listable before reaching here). */
    if (n->type == EXPR_REAL || n->type == EXPR_SYMBOL || is_complex(n, &re, &im)) {
        if (!(n->type == EXPR_SYMBOL && n->data.symbol.name == SYM_List)) {
            char* s = expr_to_string(n);
            mth_message("FactorInteger", "exact", "Argument %s in FactorInteger is not an exact number.",
                        s ? s : "");
            free(s);
        }
    }
    return NULL;
}

Expr* builtin_factorinteger(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* Fast path: FactorInteger[n] for a machine integer n >= 2 — no context,
     * no bignums, no heap until the result list (most calls are this). */
    if (res->data.function.arg_count == 1 && res->data.function.args[0]->type == EXPR_INTEGER &&
        res->data.function.args[0]->data.integer >= 2) {
        uint64_t p[64];
        int64_t e[64];
        int k = fi_u64_factor_array((uint64_t)res->data.function.args[0]->data.integer, p, e);
        Expr* items[64];
        for (int i = 0; i < k; i++) items[i] = pair_expr(expr_new_integer((int64_t)p[i]), e[i]);
        return expr_new_function(expr_new_symbol(SYM_List), items, (size_t)k);
    }
    FiCtx c;
    fi_ctx_init(&c, FI_AUTOMATIC);
    Expr* out = factorinteger_core(res, &c, 0);
    fi_ctx_clear(&c);
    return out;
}

/* FactorInteger`Name[n, ...]: the method comes from the head's name. */
static Expr* builtin_factorinteger_method(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* full = res->data.function.head->data.symbol.name;
    const char* tick = strrchr(full, '`');
    FiMethod m;
    int variant = -1;
    if (!tick || !lookup_method(tick + 1, &m, &variant)) return NULL;
    FiCtx c;
    fi_ctx_init(&c, m);
    if (variant >= 0) c.rbd.variant = variant;
    /* Fast path: FactorInteger`Method[n] on a machine integer, for the
     * methods with a word splitter. */
    if (res->data.function.arg_count == 1 && res->data.function.args[0]->type == EXPR_INTEGER &&
        res->data.function.args[0]->data.integer >= 2) {
        uint64_t p[64];
        int64_t e[64];
        int k = fi_word_method_array((uint64_t)res->data.function.args[0]->data.integer, &c, p, e);
        if (k >= 0) {
            Expr* items[64];
            for (int i = 0; i < k; i++) items[i] = pair_expr(expr_new_integer((int64_t)p[i]), e[i]);
            fi_ctx_clear(&c);
            return expr_new_function(expr_new_symbol(SYM_List), items, (size_t)k);
        }
    }
    Expr* out = factorinteger_core(res, &c, 1);
    fi_ctx_clear(&c);
    return out;
}

/* ------------------------------------------------------------------ */
/* Public helpers for other modules                                    */
/* ------------------------------------------------------------------ */

/* Complete factorisation of |n| with the Automatic method (see facint.h). */
int fi_factor_complete_mpz(const mpz_t n, FiFactors* out) {
    FiCtx c;
    fi_ctx_init(&c, FI_AUTOMATIC);
    c.warned = 1;                         /* callers decline on incompleteness */
    mpz_t a;
    mpz_init(a);
    mpz_abs(a, n);
    fi_factor(a, out, &c);
    mpz_clear(a);
    int ok = !c.incomplete;
    fi_ctx_clear(&c);
    fi_factors_sort(out);
    return ok;
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

static const char FACTORINTEGER_DOC[] =
    "FactorInteger[n] gives a list of the prime factors of the integer n, together\n"
    "\twith their exponents: {{p1, e1}, {p2, e2}, ...}.  For a negative n the first\n"
    "\tentry is {-1, 1}; for a rational n the denominator's primes have negative\n"
    "\texponents; a Gaussian integer or Gaussian rational is factored over Z[i]\n"
    "\t(unit first, primes in the first quadrant).\n"
    "FactorInteger[n, k] does a partial factorisation into at most k entries:\n"
    "\tthe small primes found by trial division (below 2^16) keep\n"
    "\ttheir k - 1 largest and merge the rest, a large cofactor stays unsplit when there\n"
    "\tis no room (FactorInteger[360, 2] is {{5, 1}, {72, 1}}). k = All or Infinity is\n"
    "\tthe complete factorisation; k = 1 gives {{n, 1}}.\n"
    "FactorInteger[n, Automatic] is FactorInteger[n].\n"
    "Options:\n"
    "\tGaussianIntegers -> True factors over the Gaussian integers.\n"
    "\tMethod -> m chooses the algorithm; m is a name below, Automatic, or\n"
    "\t{name, \"Sub\" -> v, ...} with suboptions.  Each method is also callable\n"
    "\tdirectly as FactorInteger`name[n, \"Sub\" -> v, ...] (see ?FactorInteger`name).\n"
    "Methods:\n"
    "\tAutomatic                  word engine below 2^64 (micro-ECM, rho, SQUFOF);\n"
    "\t                           above: trial division to 2^16, rho, micro-ECM,\n"
    "\t                           Fermat probe, p-1, p+1, then ECM levels, switching\n"
    "\t                           to the quadratic sieve for composites up to 100 digits\n"
    "\t\"TrialDivision\"            batched trial division        (\"Bound\")\n"
    "\t\"PollardRho\"               Pollard-Brent rho             (\"Iterations\", \"Starts\")\n"
    "\t\"Fermat\"                   Fermat, residue-sieved        (\"Iterations\")\n"
    "\t\"Lehman\"                   Lehman, O(n^(1/3))            (\"Iterations\")\n"
    "\t\"ShanksSquareForms\"        SQUFOF; alias \"SQUFOF\"        (\"Iterations\")\n"
    "\t\"CFRAC\"                    continued fractions           (\"Iterations\")\n"
    "\t\"Dixon\"                    random squares                (\"Iterations\")\n"
    "\t\"QuadraticSieve\"           SIQS; aliases \"QS\", \"SIQS\"\n"
    "\t\"PollardPMinus1\"           p-1; alias \"PollardP-1\"       (\"B1\", \"B2\")\n"
    "\t\"WilliamsPPlus1\"           p+1; alias \"WilliamsP+1\"      (\"B1\", \"B2\")\n"
    "\t\"ECM\"                      elliptic curves (GMP-ECM)     (\"B1\", \"B2\", \"Curves\", \"Sigma\")\n"
    "\t\"BlakeRationalBaseDescent\" rational base descent         (\"Base\", \"Variant\", ...)\n"
    "A composite no factor of which is found within a method's bounds is returned\n"
    "\twith exponent 1 and FactorInteger::nofac.";

static void reg(const char* name, const char* doc) {
    symtab_add_builtin(name, builtin_factorinteger_method);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_set_docstring(name, doc);
}

void fi_init(void) {
    symtab_add_builtin("FactorInteger", builtin_factorinteger);
    symtab_get_def("FactorInteger")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_set_docstring("FactorInteger", FACTORINTEGER_DOC);
    char buf[96];
    for (size_t i = 0; i < NMETHODS; i++) {
        snprintf(buf, sizeof buf, "FactorInteger`%s", METHODS[i].name);
        reg(buf, METHODS[i].doc);
    }
    for (size_t i = 0; i < sizeof SYM_ALIASES / sizeof SYM_ALIASES[0]; i++) {
        snprintf(buf, sizeof buf, "FactorInteger`%s", SYM_ALIASES[i].sym);
        reg(buf, method_info(SYM_ALIASES[i].m)->doc);
    }
    for (size_t i = 0; i < NRBD_VARIANTS; i++) {
        snprintf(buf, sizeof buf, "FactorInteger`%s%s", RBD_PREFIX, RBD_VARIANTS[i].name);
        reg(buf, RBD_VARIANTS[i].doc);
    }
}
