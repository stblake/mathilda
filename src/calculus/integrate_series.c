/* integrate_series.c
 *
 * Definite integration by series expansion + term-by-term integration +
 * resummation.  See integrate_series.h for the method and its correctness
 * argument.  This is the alternating sibling of integrate_ratlog01_try
 * (integrate_beta.c): the 1 + x^d family it closes resums to Catalan /
 * DirichletBeta / eta via the existing Sum[] engine, which integrate_ratlog01's
 * HurwitzZeta form cannot represent.
 */

#include "integrate_series.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "arithmetic.h"      /* arith_warnings_mute_push/pop */

#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* ---- construction / evaluation helpers (mirrors integrate_beta.c) -------- */

static Expr* cp(const Expr* e) { return expr_copy((Expr*)e); }
static Expr* mk_int(long v) { return expr_new_integer((int64_t)v); }

static Expr* mk_fn(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
static Expr* mk_fn1(const char* h, Expr* a) { Expr* v[1]={a}; return mk_fn(h,v,1); }
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) { Expr* v[2]={a,b}; return mk_fn(h,v,2); }
static Expr* mk_fn3b(const char* h, Expr* a, Expr* b, Expr* c) { Expr* v[3]={a,b,c}; return mk_fn(h,v,3); }
static Expr* Pow_(Expr* a, Expr* b)   { return mk_fn2("Power", a, b); }
static Expr* Times_(Expr* a, Expr* b) { return mk_fn2("Times", a, b); }
static Expr* Plus_(Expr* a, Expr* b)  { return mk_fn2("Plus", a, b); }

static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);   /* evaluate does not free its input */
    expr_free(call);
    return r;
}
static Expr* ev1(const char* name, Expr* a) { return eval_take(mk_fn1(name, a)); }
static Expr* simp(Expr* e) { return ev1("Simplify", e); }
static Expr* simp2(Expr* e, Expr* as) {
    if (!as) return simp(e);
    return eval_take(mk_fn2("Simplify", e, cp(as)));
}

static bool head_name_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}
static bool sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}
static bool is_symbol(const Expr* e, const Expr* x) {
    return e->type == EXPR_SYMBOL && x->type == EXPR_SYMBOL &&
           e->data.symbol.name == x->data.symbol.name;
}
static bool contains_symbol(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return is_symbol(e, x);
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol(e->data.function.args[i], x)) return true;
    return false;
}

/* True iff `e` mentions a divergence / indeterminacy symbol or an unevaluated
 * Integrate/Sum -- i.e. it is NOT a usable finite closed form. */
static bool mentions_nonfinite(const Expr* e) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL)
        return sym_is(e, "Infinity") || sym_is(e, "ComplexInfinity") ||
               sym_is(e, "Indeterminate") || sym_is(e, "Undefined");
    if (head_name_is(e, "DirectedInfinity") || head_name_is(e, "Integrate") ||
        head_name_is(e, "Sum"))
        return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (mentions_nonfinite(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (mentions_nonfinite(e->data.function.args[i])) return true;
    return false;
}

static Expr* coeff_(const Expr* e, const Expr* x, long k) {
    return eval_take(mk_fn3b("Coefficient", cp(e), cp((Expr*)x), mk_int(k)));
}
static bool is_zero_simp_(Expr* e) {   /* consumes e */
    Expr* s = simp(e);
    bool z = (s->type == EXPR_INTEGER && s->data.integer == 0);
    expr_free(s);
    return z;
}
/* long-valued Exponent[e, x]; returns -1 for the zero polynomial (Exponent 0 is
 * -Infinity) or any non-integer result. */
static long exponent_(const Expr* e, const Expr* x) {
    Expr* ex = eval_take(mk_fn2("Exponent", cp(e), cp((Expr*)x)));
    long d = (ex->type == EXPR_INTEGER) ? (long) ex->data.integer : -1;
    expr_free(ex);
    return d;
}
static bool is_true_(Expr* e) { bool t = sym_is(e, "True"); expr_free(e); return t; }
static bool poly_in_x(const Expr* e, const Expr* x) {
    return is_true_(eval_take(mk_fn2("PolynomialQ", cp(e), cp((Expr*)x))));
}

/* Evaluate N[e]; extract a finite machine real.  Returns 1 on success. */
static int series_real_double(Expr* e, double* out) {
    Expr* n = ev1("N", cp(e));
    if (!n) return 0;
    int ok = 0;
    if (n->type == EXPR_REAL)         { *out = n->data.real;            ok = 1; }
    else if (n->type == EXPR_INTEGER) { *out = (double) n->data.integer; ok = 1; }
    expr_free(n);
    return ok;
}

/* Numeric backstop: the closed form V must agree with NIntegrate to a loose
 * relative tolerance.  Rejects (returns false) if either side is not a finite
 * real -- a free symbolic parameter, a complex value, or an NIntegrate that
 * cannot produce a number all decline, keeping the zero-wrong-answer invariant.
 * (Every in-scope integrand on [0,1] is fully numeric.) */
static bool series_crosscheck_ok(Expr* f, Expr* x, Expr* a, Expr* b, Expr* V) {
    double vd;
    if (!series_real_double(V, &vd)) return false;

    arith_warnings_mute_push();
    Expr* spec = mk_fn3b("List", cp(x), cp(a), cp(b));
    Expr* ni   = eval_take(mk_fn2("NIntegrate", cp(f), spec));   /* spec consumed */
    arith_warnings_mute_pop();

    double nd;
    int k = ni ? series_real_double(ni, &nd) : 0;
    if (ni) expr_free(ni);
    if (k != 1) return false;
    return fabs(vd - nd) <= 1e-3 * (1.0 + fabs(vd));
}

/* The exact weight Integrate[x^m Log[x]^k, {x,0,1}] = (-1)^k k! / (m+1)^(k+1),
 * returned as sign*k! times (m+1)^-(k+1).  k >= 1, m >= 0. */
static Expr* logmoment_weight(long k, long m) {
    long fact = 1;
    for (long i = 2; i <= k; i++) fact *= i;
    long sign = (k % 2 == 0) ? 1 : -1;
    return Times_(mk_int(sign * fact), Pow_(mk_int(m + 1), mk_int(-(k + 1))));
}

#define SERIES_DMAX 24

/* ---- the method --------------------------------------------------------- */

Expr* integrate_series_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    /* Interval exactly [0,1] in v1. */
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    if (!(b->type == EXPR_INTEGER && b->data.integer == 1)) return NULL;
    if (!contains_symbol(f, x)) return NULL;

    /* 1. Strip Log[x]^k (k >= 1); collect the remaining factors as `rest`. */
    int k = 0;
    Expr* restparts[64]; size_t nrest = 0;
    Expr** fac; size_t nf; Expr* one[1];
    if (head_name_is(f, "Times")) { fac = f->data.function.args; nf = f->data.function.arg_count; }
    else { one[0] = f; fac = one; nf = 1; }
    for (size_t i = 0; i < nf; i++) {
        Expr* fe = fac[i];
        if (head_name_is(fe, "Log") && fe->data.function.arg_count == 1 &&
            expr_eq(fe->data.function.args[0], (Expr*)x)) { k += 1; continue; }
        if (head_name_is(fe, "Power") && fe->data.function.arg_count == 2 &&
            head_name_is(fe->data.function.args[0], "Log") &&
            fe->data.function.args[0]->data.function.arg_count == 1 &&
            expr_eq(fe->data.function.args[0]->data.function.args[0], (Expr*)x) &&
            fe->data.function.args[1]->type == EXPR_INTEGER &&
            fe->data.function.args[1]->data.integer >= 1) {
            k += (int) fe->data.function.args[1]->data.integer; continue;
        }
        if (nrest < 64) restparts[nrest++] = fe; else return NULL;
    }
    if (k < 1 || k > 20) return NULL;   /* k>=1 is the Log weight; k=0 is N-L's */

    Expr* rest;
    if (nrest == 0) rest = mk_int(1);
    else if (nrest == 1) rest = cp(restparts[0]);
    else {
        Expr** c2 = malloc(nrest * sizeof(Expr*));
        for (size_t i = 0; i < nrest; i++) c2[i] = cp(restparts[i]);
        rest = mk_fn("Times", c2, nrest);
        free(c2);
    }

    /* 2. rest rational in x; split N(x)/Q(x), Q non-constant, Q(0) != 0. */
    Expr* tog = ev1("Together", rest);          /* consumes rest */
    Expr* N   = ev1("Numerator", cp(tog));
    Expr* Q   = ev1("Denominator", tog);        /* consumes tog */
    if (!poly_in_x(N, x) || !poly_in_x(Q, x) || !contains_symbol(Q, x)) {
        expr_free(N); expr_free(Q); return NULL;
    }
    Expr* q0 = eval_take(mk_fn2("ReplaceAll", cp(Q), mk_fn2("Rule", cp(x), mk_int(0))));
    bool q0_zero = is_zero_simp_(q0);           /* consumes q0 */
    if (q0_zero) { expr_free(N); expr_free(Q); return NULL; }   /* pole at x=0 */

    /* 3. Smallest d, sign eps with Q | (1 - eps x^d).  denom_poly = 1 - eps x^d,
     *    so eps=-1 is 1+x^d and eps=+1 is 1-x^d.  Test eps=-1 first so a factor
     *    with no real pole at x=1 (e.g. 1+x^2) takes its natural alternating
     *    form rather than a larger 1-x^d multiple. */
    long d_found = 0; int eps_found = 0; Expr* cof = NULL;
    for (long d = 1; d <= SERIES_DMAX && !d_found; d++) {
        int signs[2] = { -1, +1 };
        for (int si = 0; si < 2 && !d_found; si++) {
            int eps = signs[si];
            /* denom_poly = 1 - eps*x^d */
            Expr* denom = Plus_(mk_int(1), Times_(mk_int(-eps), Pow_(cp(x), mk_int(d))));
            Expr* rem = eval_take(mk_fn3b("PolynomialRemainder", cp(denom), cp(Q), cp(x)));
            bool divides = (rem->type == EXPR_INTEGER && rem->data.integer == 0);
            expr_free(rem);
            if (divides) {
                cof = eval_take(mk_fn3b("PolynomialQuotient", cp(denom), cp(Q), cp(x)));
                d_found = d; eps_found = eps;
            }
            expr_free(denom);
        }
    }
    if (!d_found) { expr_free(N); expr_free(Q); return NULL; }
    long d = d_found; int eps = eps_found;

    /* N := Expand[N * cof]; integrand (sans Log) = N / (1 - eps x^d). */
    Expr* Nn = ev1("Expand", Times_(N, cof));   /* consumes N and cof */
    expr_free(Q);

    /* 4. Split N = Quot (1 - eps x^d) + S, deg S < d. */
    Expr* denom_poly = Plus_(mk_int(1), Times_(mk_int(-eps), Pow_(cp(x), mk_int(d))));
    Expr* Quot = eval_take(mk_fn3b("PolynomialQuotient", cp(Nn), cp(denom_poly), cp(x)));
    Expr* S    = eval_take(mk_fn3b("PolynomialRemainder", cp(Nn), denom_poly, cp(x))); /* denom_poly consumed */
    expr_free(Nn);

    Expr* acc = mk_int(0);

    /* 4a. Polynomial part Quot(x): finite exact sum of the log-moment identity. */
    long degQ = exponent_(Quot, x);             /* -1 if Quot == 0 */
    for (long m = 0; m <= degQ; m++) {
        Expr* cm = coeff_(Quot, x, m);
        if (cm->type == EXPR_INTEGER && cm->data.integer == 0) { expr_free(cm); continue; }
        acc = Plus_(acc, Times_(cm, logmoment_weight(k, m)));
    }
    expr_free(Quot);

    /* 4b. Proper part S/(1 - eps x^d): one Sum per residue class j (0 <= j < d).
     *     s_j (-1)^k k! Sum[eps^n / (d n + j + 1)^(k+1), {n,0,Infinity}]. */
    long degS = exponent_(S, x);                /* -1 if S == 0 */
    bool ok = true;
    for (long j = 0; j <= degS && ok; j++) {
        Expr* sj = coeff_(S, x, j);
        if (sj->type == EXPR_INTEGER && sj->data.integer == 0) { expr_free(sj); continue; }

        Expr* n = expr_new_symbol("Integrate`Series`n");
        /* base = d*n + (j+1) */
        Expr* base = Plus_(Times_(mk_int(d), cp(n)), mk_int(j + 1));
        Expr* summand = (eps == 1)
            ? Pow_(base, mk_int(-(k + 1)))
            : Times_(Pow_(mk_int(-1), cp(n)), Pow_(base, mk_int(-(k + 1))));
        Expr* sumexpr = mk_fn2("Sum", summand,
                               mk_fn3b("List", cp(n), mk_int(0), expr_new_symbol("Infinity")));
        Expr* inner = eval_take(sumexpr);

        /* Reject an unclosed / symbolic-in-n / nonfinite resummation. */
        bool bad = !inner || head_name_is(inner, "Sum") || mentions_nonfinite(inner) ||
                   (inner && contains_symbol(inner, n));
        expr_free(n);
        if (bad) {
            if (inner) expr_free(inner);
            expr_free(sj);
            ok = false;
            break;
        }
        acc = Plus_(acc, Times_(sj, Times_(logmoment_weight(k, 0), inner)));
    }
    expr_free(S);
    if (!ok) { expr_free(acc); return NULL; }

    Expr* val = assumptions ? simp2(acc, assumptions) : simp(acc);
    if (!val || contains_symbol(val, x) || mentions_nonfinite(val)) {
        if (val) expr_free(val);
        return NULL;
    }

    /* 5. Numeric backstop against NIntegrate. */
    if (!series_crosscheck_ok(f, x, a, b, val)) { expr_free(val); return NULL; }
    return val;
}

Expr* builtin_integrate_series(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!head_name_is(spec, "List") || spec->data.function.arg_count != 3) return NULL;
    Expr* x = spec->data.function.args[0];
    Expr* a = spec->data.function.args[1];
    Expr* b = spec->data.function.args[2];
    if (x->type != EXPR_SYMBOL) return NULL;
    Expr* as = NULL;
    for (size_t t = 2; t < res->data.function.arg_count; t++) {
        Expr* opt = res->data.function.args[t];
        if (opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2 &&
            opt->data.function.head->type == EXPR_SYMBOL &&
            (opt->data.function.head->data.symbol.name == SYM_Rule ||
             opt->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
            Expr* lhs = opt->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == SYM_Assumptions) {
                as = opt->data.function.args[1];
                continue;
            }
        }
        return NULL;
    }
    return integrate_series_try(f, x, a, b, as);
}

void integrate_series_init(void) {
    symtab_add_builtin("Integrate`SeriesExpansion", builtin_integrate_series);
    symtab_get_def("Integrate`SeriesExpansion")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`SeriesExpansion",
        "Integrate`SeriesExpansion[f, {x, 0, 1}] evaluates a log-weighted "
        "rational integral R(x) Log[x]^k over [0,1] (k >= 1) whose denominator "
        "is a root-of-unity cyclotomic factor 1 +/- x^d, by expanding the "
        "geometric series, integrating term-by-term via "
        "Integrate[x^m Log[x]^k, {x,0,1}] = (-1)^k k!/(m+1)^(k+1), and resumming "
        "with Sum[] (-> Catalan / DirichletBeta / Zeta).  The alternating "
        "1 + x^d case that Integrate`RationalLog (HurwitzZeta) cannot represent "
        "is this mechanism's.  Numerically cross-checked against NIntegrate; "
        "returns unevaluated when the integrand is not of this form, the "
        "interval is not [0,1], or the resummation does not close.");
}
