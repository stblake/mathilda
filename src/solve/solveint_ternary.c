/*
 * solveint_ternary.c
 *
 * Part of the Solve[..., Integers] engine; split out of solveint.c.
 * See solveint_internal.h for the shared SICtx/SearchState substrate.
 *
 * Tier 2 E of SOLVE_INTEGERS.md: the homogeneous ternary quadratic /
 * Legendre solver.  For a single homogeneous diagonal degree-2 equation in
 * exactly 3 variables over an unbounded domain, decide solvability RIGOROUSLY
 * (Legendre's theorem, a proof -- not a bounded search) and either
 *   - return the complete 2-parameter integer family (a conic with one rational
 *     point has a rational parametrisation), or
 *   - return the trivial-only list {{x->0, y->0, z->0}} when Legendre proves no
 *     nontrivial solution exists (this is what Mathematica returns, NOT {}).
 *
 * Phase 1 scope: the symmetric pattern  x^2 + y^2 == k z^2  (after clearing gcd
 * and normalising: two coefficients equal -- hence +-1 once pairwise coprime --
 * and the third of opposite sign).  This covers the Wikipedia targets
 *   x^2 + y^2 == z^2   (Pythagorean, k = 1),
 *   x^2 + y^2 == 3 z^2 (trivial-only, Legendre fails),
 *   x^2 + y^2 == 2 z^2 (solvable, k = 2).
 * Cross-term (non-diagonal) forms and the a != b general ternary decline (a
 * documented follow-up); a decline is always safe.
 */
#include "solveint.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "eval.h"
#include "expr.h"
#include "poly/mpoly.h"
#include "numbertheory/numbertheory_internal.h"
#include "solveint_internal.h"


/* Read the diagonal coefficients (a, b, c) of a homogeneous ternary quadratic
 * over the three active variables vs[0..2].  Returns false on any term that is
 * not a pure square of one of the three variables (a cross term x*y, a linear
 * term, or a constant) -- those shapes are the general-BQF follow-up. */
static bool ternary_extract(const MPoly* eq, int n, const int* vs,
                            mpz_t a, mpz_t b, mpz_t c) {
    mpz_set_ui(a, 0); mpz_set_ui(b, 0); mpz_set_ui(c, 0);
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        /* No variable outside the three actives may appear. */
        for (int v = 0; v < n; v++)
            if (v != vs[0] && v != vs[1] && v != vs[2] && ex[v] != 0) return false;
        int d0 = ex[vs[0]], d1 = ex[vs[1]], d2 = ex[vs[2]];
        if      (d0 == 2 && d1 == 0 && d2 == 0) mpz_set(a, eq->coefs[t]);
        else if (d0 == 0 && d1 == 2 && d2 == 0) mpz_set(b, eq->coefs[t]);
        else if (d0 == 0 && d1 == 0 && d2 == 2) mpz_set(c, eq->coefs[t]);
        else return false;                       /* cross / linear / constant */
    }
    return (mpz_sgn(a) != 0 && mpz_sgn(b) != 0 && mpz_sgn(c) != 0);
}


/* True iff m (nonzero) is squarefree.  Via the trial/df factoriser; a factor
 * failure (astronomically large m) conservatively reports "not squarefree" so
 * the solver declines rather than guess. */
static bool ternary_squarefree(const mpz_t m) {
    mpz_t am; mpz_init(am); mpz_abs(am, m);
    if (mpz_cmp_ui(am, 1) <= 0) { mpz_clear(am); return true; }
    mpz_t* primes = NULL; unsigned long* exps = NULL; size_t np = 0;
    bool sf = true;
    if (df_factor_mpz(am, &primes, &exps, &np)) {
        for (size_t i = 0; i < np; i++) if (exps[i] > 1) { sf = false; break; }
        for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
        free(primes); free(exps);
    } else sf = false;
    mpz_clear(am);
    return sf;
}


/* Legendre solvability of  x^2 + y^2 == k z^2  (k > 0 squarefree): a nontrivial
 * integer solution exists iff -1 is a quadratic residue modulo every odd prime
 * factor of k (the residue condition mod 2 is vacuous).  Equivalently, no prime
 * factor of k is == 3 (mod 4).  This is the Hasse-Minkowski / Legendre content
 * specialised to this form. */
static bool ternary_solvable(const mpz_t k) {
    mpz_t primes[PR_MAX_DISTINCT_PRIMES]; size_t np = 0;
    for (size_t i = 0; i < PR_MAX_DISTINCT_PRIMES; i++) mpz_init(primes[i]);
    bool ok = pr_collect_distinct_primes(k, primes, &np);
    mpz_t neg1; mpz_init_set_si(neg1, -1);
    for (size_t i = 0; i < np && ok; i++) {
        if (mpz_cmp_ui(primes[i], 2) == 0) continue;          /* mod 2 vacuous */
        if (mpz_legendre(neg1, primes[i]) == -1) ok = false;  /* -1 non-residue */
    }
    mpz_clear(neg1);
    for (size_t i = 0; i < PR_MAX_DISTINCT_PRIMES; i++) mpz_clear(primes[i]);
    return ok;
}


/* C[k]. */
static Expr* mk_Ck(int k) { return mk_fn1("C", mk_int(k)); }

/* A fresh quadratic form  css*C[1]^2 + ctt*C[2]^2 + cst*C[1]*C[2]. */
static Expr* mk_qform(int64_t css, int64_t ctt, int64_t cst) {
    Expr* terms[3]; int nt = 0;
    if (css) terms[nt++] = mk_fn2("Times", mk_int(css),
                                  mk_fn2("Power", mk_Ck(1), mk_int(2)));
    if (ctt) terms[nt++] = mk_fn2("Times", mk_int(ctt),
                                  mk_fn2("Power", mk_Ck(2), mk_int(2)));
    if (cst) terms[nt++] = expr_new_function(mk_sym("Times"),
                                  (Expr*[]){ mk_int(cst), mk_Ck(1), mk_Ck(2) }, 3);
    if (nt == 0) return mk_int(0);
    return expr_new_function(mk_sym("Plus"), terms, nt);
}

/* C[3] * (css*C1^2 + ctt*C2^2 + cst*C1*C2). */
static Expr* mk_scaled_q(int64_t css, int64_t ctt, int64_t cst) {
    return mk_fn2("Times", mk_Ck(3), mk_qform(css, ctt, cst));
}

/* Assemble a solution tuple {v_pa->va, v_pb->vb, v_po->vo} with the rules placed
 * in ascending variable-index order (a canonical tuple). */
static Expr* ternary_tuple(SICtx* c, int pa, int pb, int po,
                           Expr* va, Expr* vb, Expr* vo) {
    Expr* byidx[3] = {NULL, NULL, NULL};
    byidx[pa] = mk_rule(expr_copy(c->var[pa]), va);
    byidx[pb] = mk_rule(expr_copy(c->var[pb]), vb);
    byidx[po] = mk_rule(expr_copy(c->var[po]), vo);
    Expr* rr[3]; int m = 0; for (int i = 0; i < 3; i++) if (byidx[i]) rr[m++] = byidx[i];
    return mk_list(rr, 3);
}


/* The trivial-only answer  {{v0->0, v1->0, v2->0}}. */
static Expr* ternary_trivial(SICtx* c, const int* vs) {
    Expr* rules[3];
    for (int i = 0; i < 3; i++) rules[i] = mk_rule(expr_copy(c->var[vs[i]]), mk_int(0));
    Expr* tuple = mk_list(rules, 3);
    return mk_list(&tuple, 1);
}


/* Homogeneous ternary quadratic / Legendre solver.  Returns the owned answer
 * (a parametric family List, or the trivial-only list), or NULL to decline. */
Expr* si_solve_ternary_quadratic(SICtx* c) {
    if (c->neq != 1 || c->n != 3) return NULL;

    /* Fully unbounded, no orderings / disequations / abs-orderings, everything
     * captured -- defer any boxed or constrained instance to the bounded engine
     * (which already solves them completely). */
    for (int i = 0; i < 3; i++) if (c->has_lo[i] || c->has_hi[i]) return NULL;
    if (c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured)
        return NULL;

    int vs[3] = {0, 1, 2};
    for (int i = 0; i < 3; i++) if (mpoly_deg_var(c->eq[0], vs[i]) < 1) return NULL;

    mpz_t a, b, cc, g;
    mpz_inits(a, b, cc, g, NULL);
    Expr* result = NULL;
    bool declined = false;

    if (!ternary_extract(c->eq[0], c->n, vs, a, b, cc)) { declined = true; goto done; }

    /* Primitive: divide out gcd(a, b, c) (signs preserved). */
    mpz_gcd(g, a, b); mpz_gcd(g, g, cc);
    if (mpz_sgn(g) != 0 && mpz_cmp_ui(g, 1) != 0) {
        mpz_divexact(a, a, g); mpz_divexact(b, b, g); mpz_divexact(cc, cc, g);
    }

    /* Phase 1: squarefree and pairwise coprime.  (No-op for the targets; a
     * non-squarefree / non-coprime form is a documented follow-up -> decline.) */
    if (!ternary_squarefree(a) || !ternary_squarefree(b) || !ternary_squarefree(cc)) {
        declined = true; goto done;
    }
    { mpz_t gab, gac, gbc; mpz_inits(gab, gac, gbc, NULL);
      mpz_gcd(gab, a, b); mpz_gcd(gac, a, cc); mpz_gcd(gbc, b, cc);
      bool coprime = (mpz_cmp_ui(gab, 1) == 0 && mpz_cmp_ui(gac, 1) == 0 &&
                      mpz_cmp_ui(gbc, 1) == 0);
      mpz_clears(gab, gac, gbc, NULL);
      if (!coprime) { declined = true; goto done; }
    }

    /* Real condition: all same sign => only the trivial real point. */
    if (mpz_sgn(a) == mpz_sgn(b) && mpz_sgn(b) == mpz_sgn(cc)) {
        result = ternary_trivial(c, vs); goto done;
    }

    /* Find the symmetric pair (two equal coefficients -> both +-1 by coprimality)
     * and the odd-one-out.  pa,pb = the pair's variable indices; po = the "z". */
    int pa, pb, po;
    if      (mpz_cmp(a, b) == 0) { pa = vs[0]; pb = vs[1]; po = vs[2]; }
    else if (mpz_cmp(a, cc) == 0) { pa = vs[0]; pb = vs[2]; po = vs[1]; }
    else if (mpz_cmp(b, cc) == 0) { pa = vs[1]; pb = vs[2]; po = vs[0]; }
    else { declined = true; goto done; }         /* a != b != c: follow-up */

    /* k = |coeff of the odd var|; the form is  x^2 + y^2 == k z^2. */
    mpz_t k; mpz_init(k);
    if (po == vs[0]) mpz_abs(k, a);
    else if (po == vs[1]) mpz_abs(k, b);
    else mpz_abs(k, cc);
    if (!mpz_fits_slong_p(k)) { mpz_clear(k); declined = true; goto done; }
    int64_t kk = mpz_get_si(k);

    /* Legendre solvability. */
    if (!ternary_solvable(k)) { mpz_clear(k); result = ternary_trivial(c, vs); goto done; }

    /* Single-representation gate.  A k with two or more distinct prime factors
     * == 1 (mod 4) has several inequivalent sum-of-two-squares representations,
     * hence several solution CLASSES that one witness's family cannot all reach
     * (e.g. 65 = 5*13: both (1,8) and (4,7)).  Emitting one witness there would
     * be an incomplete -- i.e. wrong -- answer, so decline (a follow-up: one
     * family per representation, mirroring genpell's per-class families). */
    { int64_t m = kk, primes1mod4 = 0;
      for (int64_t d = 2; d * d <= m; d++)
          if (m % d == 0) { if (d % 4 == 1) primes1mod4++; while (m % d == 0) m /= d; }
      if (m > 1 && m % 4 == 1) primes1mod4++;
      if (primes1mod4 > 1) { mpz_clear(k); declined = true; goto done; }
    }
    mpz_clear(k);

    /* Witness on  A^2 + B^2 == k  (guaranteed once Legendre passes: k is then a
     * sum of two squares).  P0 = (A0, B0, 1). */
    int64_t A0 = -1, B0 = -1;
    { int64_t amax = si_isqrt_i64(kk);
      for (int64_t A = 0; A <= amax; A++) {
          int64_t rem = kk - A * A;
          int64_t r = si_isqrt_i64(rem);
          if (r * r == rem) { A0 = A; B0 = r; break; }
      } }
    if (A0 < 0) { declined = true; goto done; }   /* unreachable if Legendre holds */

    /* Complete parametric family.  Chord/tangent from P0 (drop the z axis;
     * V = s e_A + t e_B, s=C[1], t=C[2], overall scale C[3]):
     *   phi_A = -A0 s^2 + A0 t^2 - 2 B0 s t
     *   phi_B =  B0 s^2 - B0 t^2 - 2 A0 s t
     *   phi_C =  s^2 + t^2                (unit coeff -> content 1, no reduction)
     * Only squares appear in x^2+y^2==k z^2, so every sign/swap image of a
     * solution is a solution; the z-sign is the C[3] (resp. C[1]) sign.  The
     * complete set is the union of two orbits, each closed under the group
     * <x<->y, x->-x, y->-y> (order 8, but z-sign folds into the scale so 8
     * tuples suffice per orbit): the phi-family, plus the tangent family C[1]*P0
     * -- which supplies the odd multiples of P0 that phi (giving 2*P0 at the
     * tangent direction) misses.  Verified complete + sound against brute force
     * for every single-representation solvable k up to 200. */
    Expr* raw[16]; int nraw = 0;
    const int64_t pAs = -A0, pAt = A0, pAst = -2 * B0;   /* phi_A coeffs */
    const int64_t pBs =  B0, pBt = -B0, pBst = -2 * A0;  /* phi_B coeffs */
    for (int sx = 1; sx >= -1; sx -= 2)
        for (int sy = 1; sy >= -1; sy -= 2) {
            /* (sx phi_A, sy phi_B, phi_C) and its A<->B swap */
            raw[nraw++] = ternary_tuple(c, pa, pb, po,
                mk_scaled_q(sx * pAs, sx * pAt, sx * pAst),
                mk_scaled_q(sy * pBs, sy * pBt, sy * pBst), mk_scaled_q(1, 1, 0));
            raw[nraw++] = ternary_tuple(c, pa, pb, po,
                mk_scaled_q(sx * pBs, sx * pBt, sx * pBst),
                mk_scaled_q(sy * pAs, sy * pAt, sy * pAst), mk_scaled_q(1, 1, 0));
        }
    for (int sx = 1; sx >= -1; sx -= 2)
        for (int sy = 1; sy >= -1; sy -= 2) {
            /* tangent family C[1]*(sx A0, sy B0, 1) and its A0<->B0 swap */
            raw[nraw++] = ternary_tuple(c, pa, pb, po,
                mk_fn2("Times", mk_int(sx * A0), mk_Ck(1)),
                mk_fn2("Times", mk_int(sy * B0), mk_Ck(1)), mk_Ck(1));
            raw[nraw++] = ternary_tuple(c, pa, pb, po,
                mk_fn2("Times", mk_int(sx * B0), mk_Ck(1)),
                mk_fn2("Times", mk_int(sy * A0), mk_Ck(1)), mk_Ck(1));
        }
    /* Evaluate and drop duplicate tuples (A0==B0 or a zero coord collapses some
     * sign/swap images to the same family). */
    Expr* uniq[16]; int nu = 0;
    for (int i = 0; i < nraw; i++) {
        Expr* e = eval_and_free(raw[i]);
        bool dup = false;
        for (int j = 0; j < nu; j++) if (expr_eq(e, uniq[j])) { dup = true; break; }
        if (dup) expr_free(e); else uniq[nu++] = e;
    }
    result = mk_list(uniq, (size_t)nu);

done:
    mpz_clears(a, b, cc, g, NULL);
    if (declined) { if (result) expr_free(result); return NULL; }
    return result;
}
