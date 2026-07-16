/* risch_hermite.c — the Hermite reduction over a monomial extension (§5.3).
 *
 * See risch_hermite.h.  Literal transcription of Bronstein's quadratic
 * HermiteReduce (Symbolic Integration I, 2nd ed., p.139):
 *
 *   (f_p, f_s, f_n) <- CanonicalRepresentation(f, D)
 *   (a, d) <- (numerator(f_n), denominator(f_n))        (d monic)
 *   (d_1, ..., d_m) <- SquareFree(d)
 *   g <- 0
 *   for i <- 2 to m such that deg(d_i) > 0:
 *       v <- d_i;  u <- d / v^i
 *       for j <- i-1 downto 1:
 *           (b, c) <- ExtendedEuclidean(u D[v], v, -a/j)   (b u D[v] + c v = -a/j)
 *           g <- g + b / v^j
 *           a <- -j c - u D[b]
 *       d <- u v
 *   (q, r) <- PolyDivide(a, d)
 *   return (g, r/d, q + f_p + f_s)
 *
 * D[v], D[b] are the MONOMIAL DERIVATION (not d/dt); the squarefree factoring
 * inside SquareFree(d) is w.r.t. t (repeated t-roots).  The Diophantine solve
 * is gcd-based: gcd(u D[v], v) = 1 because v is normal (gcd(D[v], v) = 1) and
 * gcd(u, v) = 1.
 */

#include "risch_hermite.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"
#include "risch_field.h"
#include "risch_canonical.h"

#include <stdlib.h>

/* ------------------------------------------------------------------ */
/* Evaluation helpers (mirroring the sibling risch_ modules).          */
/* ------------------------------------------------------------------ */

static Expr* rh_cp(const Expr* e) { return expr_copy((Expr*)e); }
static Expr* rh_eval_adopt(Expr* call) { Expr* r = evaluate(call); expr_free(call); return r; }
static Expr* rh_fn(const char* head, Expr** args, size_t n) {
    return rh_eval_adopt(expr_new_function(expr_new_symbol(head), args, n));
}
static Expr* rh_call1(const char* head, Expr* a) { return rh_fn(head, (Expr*[]){ a }, 1); }
/* Unevaluated builders. */
static Expr* rh_times(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){ a, b }, 2);
}
static Expr* rh_pow(Expr* a, long n) {
    return expr_new_function(expr_new_symbol("Power"), (Expr*[]){ a, expr_new_integer(n) }, 2);
}
/* Expand[a] (polynomial normal form) and Cancel[Together[a]] (reduced fraction);
 * each adopts a. */
static Expr* rh_expand(Expr* a) { return rh_call1("Expand", a); }
static Expr* rh_norm(Expr* a)   { return rh_call1("Cancel", rh_call1("Together", a)); }

/* Expanded polynomial product a*b (a, b borrowed). */
static Expr* rh_mul(const Expr* a, const Expr* b) {
    return rh_expand(rh_times(rh_cp(a), rh_cp(b)));
}
/* Expanded power base^n, n >= 0 (base borrowed). */
static Expr* rh_powexp(const Expr* base, long n) {
    return rh_expand(rh_pow(rh_cp(base), n));
}

/* ------------------------------------------------------------------ */
/* HermiteReduce (§5.3, p.139).                                        */
/* ------------------------------------------------------------------ */

bool risch_hermite_reduce(const Expr* f, const Expr* t, const RischDeriv* d,
                          Expr** go, Expr** ho, Expr** ro) {
    if (t->type != EXPR_SYMBOL || !risch_deriv_lookup(d, t)) return false;

    /* (f_p, f_s, f_n) <- CanonicalRepresentation(f). */
    Expr* fp = NULL; Expr* fs = NULL; Expr* fn = NULL;
    risch_canonical_representation(f, t, d, &fp, &fs, &fn);

    /* (a, dcur) <- num/den of f_n, dcur monic in t. */
    Expr* a = NULL; Expr* dcur = NULL;
    risch_field_num_den_t(fn, t, &a, &dcur);
    expr_free(fn);

    /* (d_1, ..., d_m) <- SquareFree(dcur). */
    size_t m = 0;
    Expr** facs = risch_squarefree_t(dcur, t, &m);

    Expr* g = expr_new_integer(0);

    for (size_t i = 2; i <= m; i++) {
        const Expr* v = facs[i - 1];               /* v = d_i (borrowed) */
        if (risch_field_degree_t(v, t) <= 0) continue;
        Expr* vi = rh_powexp(v, (long)i);          /* v^i */
        Expr* u  = risch_field_divexact_t(dcur, vi, t);   /* u = dcur / v^i */
        expr_free(vi);
        if (!u) continue;                          /* defensive: not a factor */
        Expr* Dv  = risch_field_deriv(v, d);       /* monomial derivation */
        Expr* uDv = rh_mul(u, Dv);                 /* u D[v] */
        expr_free(Dv);

        for (long j = (long)i - 1; j >= 1; j--) {
            /* rhs = -a / j. */
            Expr* rhs = rh_expand(rh_times(
                expr_new_function(expr_new_symbol("Rational"),
                    (Expr*[]){ expr_new_integer(-1), expr_new_integer(j) }, 2), rh_cp(a)));
            /* (b, c): b (u D[v]) + c v = -a/j,  deg_t(b) < deg_t(v). */
            Expr* b = NULL; Expr* c = NULL;
            risch_field_diophantine_t(uDv, v, rhs, t, &b, &c);
            expr_free(rhs);

            /* g <- g + b / v^j. */
            Expr* vj = rh_powexp(v, j);
            Expr* bvj = rh_norm(rh_times(rh_cp(b), rh_pow(vj, -1)));   /* adopts vj */
            g = rh_norm(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ g, bvj }, 2));

            /* a <- -j c - u D[b]. */
            Expr* Db   = risch_field_deriv(b, d);
            Expr* uDb  = rh_mul(u, Db);
            Expr* jint = expr_new_integer(j);
            Expr* jc   = rh_mul(jint, c);          /* rh_mul copies -> free jint */
            expr_free(Db); expr_free(jint);
            Expr* anew = rh_expand(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ rh_times(expr_new_integer(-1), jc),      /* -j c */
                           rh_times(expr_new_integer(-1), uDb) }, 2));  /* - u D[b] */
            expr_free(a); a = anew;
            expr_free(b); expr_free(c);
        }
        /* dcur <- u v. */
        Expr* uv = rh_mul(u, v);
        expr_free(dcur); dcur = uv;
        expr_free(u); expr_free(uDv);
    }

    for (size_t i = 0; i < m; i++) expr_free(facs[i]);
    free(facs);

    /* (q, r) <- PolyDivide(a, dcur);  h = r/dcur (simple);  r_out = q + f_p + f_s. */
    Expr* q = NULL; Expr* rr = NULL;
    if (!risch_field_divmod_t(a, dcur, t, &q, &rr)) {   /* dcur == 0 impossible (monic) */
        q = expr_new_integer(0); rr = rh_cp(a);
    }
    expr_free(a);
    Expr* h = rh_norm(rh_times(rr, rh_pow(dcur, -1)));  /* adopts rr, dcur */
    Expr* r = rh_norm(expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ q, fp, fs }, 3));                    /* adopts q, fp, fs */

    *go = g; *ho = h; *ro = r;
    return true;
}

/* ------------------------------------------------------------------ */
/* Builtin.                                                            */
/* ------------------------------------------------------------------ */

/* Risch`HermiteReduce[f, t, deriv] -> {g, h, r}. */
static Expr* builtin_risch_hermite_reduce(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    const Expr* f = res->data.function.args[0];
    const Expr* t = res->data.function.args[1];
    if (t->type != EXPR_SYMBOL) return NULL;
    RischDeriv d;
    if (!risch_deriv_from_rules(res->data.function.args[2], &d)) return NULL;
    Expr* g = NULL; Expr* h = NULL; Expr* r = NULL;
    bool ok = risch_hermite_reduce(f, t, &d, &g, &h, &r);
    risch_deriv_free(&d);
    if (!ok) return NULL;
    return expr_new_function(expr_new_symbol("List"), (Expr*[]){ g, h, r }, 3);
}

/* ------------------------------------------------------------------ */
/* Registration.                                                       */
/* ------------------------------------------------------------------ */

void integrate_risch_hermite_init(void) {
    symtab_add_builtin("Risch`HermiteReduce", builtin_risch_hermite_reduce);
    symtab_get_def("Risch`HermiteReduce")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Risch`HermiteReduce",
        "Risch`HermiteReduce[f, t, deriv] is the Hermite reduction of f in k(t)\n"
        "over the monomial derivation deriv (Bronstein 5.3): returns {g, h, r}\n"
        "with f == D[g] + h + r, h simple (squarefree normal denominator) and r\n"
        "reduced.  The simple part h is then handled by the residue/log part and\n"
        "r by the polynomial reduction.");
}
