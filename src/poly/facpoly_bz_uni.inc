/* ===================================================================== */
/*  Univariate factorisation over Z: the classical Berlekamp-Zassenhaus  */
/*  pipeline (Cantor-Zassenhaus mod p, quadratic Hensel lifting, and     */
/*  Zassenhaus recombination).                                           */
/*                                                                       */
/*  This is the factoriser every build uses for a univariate integer     */
/*  polynomial that does not go through FLINT: all of Factor in a        */
/*  USE_FLINT=0 build, and -- in either build -- the univariate images   */
/*  that the multivariate pipeline (bivariate Hensel, the irreducibility */
/*  probe) and the rational-function path hand to bz_factor_to_expr.    */
/*                                                                       */
/*  Pipeline for a polynomial P in Z[x]:                                 */
/*    1. P = c * x^v * f with c the signed content (sign of lc(P)), f    */
/*       primitive with positive leading coefficient and f(0) != 0.      */
/*    2. Squarefree decomposition of f (Yun, over Z via the primitive    */
/*       subresultant gcd zupoly_gcd).  A mod-p squarefree certificate   */
/*       skips Yun for the common squarefree input.                      */
/*       x^n - 1 and x^n + 1 are answered by their cyclotomic factors   */
/*       directly (the inputs on which recombination is exponential).    */
/*    3. Each squarefree part g of degree 2 .. BZ_MAX_DEG:               */
/*       a. Choose a good prime: several small odd primes p with         */
/*          p | lc(g) false and g mod p squarefree; the distinct-degree  */
/*          factorisation of each gives its modular factor count and     */
/*          degree pattern.  The prime with the fewest modular factors   */
/*          wins, and the subset-sum sets of all patterns are intersected */
/*          into the set of degrees a true factor can have.  One modular */
/*          factor (or no admissible proper degree) proves irreducibility.*/
/*       b. Equal-degree factorisation (Cantor-Zassenhaus) mod the       */
/*          chosen p, in word arithmetic: p < 2^16, so every product of  */
/*          two residues is < 2^32 and a dot product of up to 2^31 terms */
/*          fits in a uint64_t -- no 128-bit type, no overflow.          */
/*       c. Quadratic Hensel lifting of the monic modular factors along  */
/*          a balanced factor tree to M = p^(2^j) > 2 B, with B the      */
/*          Landau-Mignotte bound |lc| 2^(n-1) ||g||_2 on lc(g) times    */
/*          the coefficients of any proper factor (all in GMP).          */
/*       d. Zassenhaus recombination over subsets of increasing size,    */
/*          enumerated combinatorially (no bitmask), pruned by the       */
/*          admissible-degree set and a constant-term divisibility test, */
/*          confirmed by EXACT division over Z, and with every found     */
/*          factor's modular factors removed before the search continues.*/
/*                                                                       */
/*  Correctness is unconditional: every reported factor is confirmed by  */
/*  exact division, so a run that stops early (subset budget exhausted,  */
/*  TimeConstrained deadline) returns a coarser but still exact product. */
/*  The bound only governs completeness.                                 */
/* ===================================================================== */

/* ---------------------------------------------------------------------- */
/*  Word-size polynomials over F_p, p an odd prime < 2^16                 */
/* ---------------------------------------------------------------------- */

typedef struct {
    int deg;          /* -1 for the zero polynomial */
    uint64_t* c;      /* c[0..deg], each in [0, p) */
} NPoly;

static NPoly* np_alloc(int deg) {
    NPoly* a = (NPoly*)malloc(sizeof(NPoly));
    a->deg = deg;
    a->c = (uint64_t*)calloc((size_t)(deg >= 0 ? deg + 1 : 1), sizeof(uint64_t));
    return a;
}

static void np_free(NPoly* a) {
    if (!a) return;
    free(a->c);
    free(a);
}

static void np_trim(NPoly* a) {
    while (a->deg >= 0 && a->c[a->deg] == 0) a->deg--;
}

static NPoly* np_copy(const NPoly* a) {
    NPoly* r = np_alloc(a->deg);
    if (a->deg >= 0) memcpy(r->c, a->c, sizeof(uint64_t) * (size_t)(a->deg + 1));
    return r;
}

static NPoly* np_const(uint64_t v) {
    NPoly* r = np_alloc(0);
    r->c[0] = v;
    np_trim(r);
    return r;
}

/* a^-1 mod p for a != 0 mod p (p prime). */
static uint64_t nm_inv(uint64_t a, uint64_t p) {
    int64_t t = 0, newt = 1;
    int64_t r = (int64_t)p, newr = (int64_t)(a % p);
    while (newr != 0) {
        int64_t q = r / newr, tmp;
        tmp = t - q * newt; t = newt; newt = tmp;
        tmp = r - q * newr; r = newr; newr = tmp;
    }
    if (t < 0) t += (int64_t)p;
    return (uint64_t)t;
}

static void np_make_monic(NPoly* a, uint64_t p) {
    if (a->deg < 0 || a->c[a->deg] == 1) return;
    uint64_t inv = nm_inv(a->c[a->deg], p);
    for (int i = 0; i <= a->deg; i++) a->c[i] = a->c[i] * inv % p;
}

static NPoly* np_sub(const NPoly* a, const NPoly* b, uint64_t p) {
    int d = a->deg > b->deg ? a->deg : b->deg;
    NPoly* r = np_alloc(d);
    for (int i = 0; i <= d; i++) {
        uint64_t x = (i <= a->deg ? a->c[i] : 0) + p - (i <= b->deg ? b->c[i] : 0);
        r->c[i] = x >= p ? x - p : x;
    }
    np_trim(r);
    return r;
}

/* Schoolbook product.  Residues are < 2^16, so each partial product is
 * < 2^32 and a coefficient accumulates at most min(deg)+1 < 2^31 of them
 * before the single reduction at the end: no overflow in a uint64_t. */
static NPoly* np_mul(const NPoly* a, const NPoly* b, uint64_t p) {
    if (a->deg < 0 || b->deg < 0) return np_alloc(-1);
    NPoly* r = np_alloc(a->deg + b->deg);
    for (int i = 0; i <= a->deg; i++) {
        uint64_t ai = a->c[i];
        if (!ai) continue;
        uint64_t* rc = r->c + i;
        for (int j = 0; j <= b->deg; j++) rc[j] += ai * b->c[j];
    }
    for (int i = 0; i <= r->deg; i++) r->c[i] %= p;
    np_trim(r);
    return r;
}

/* Long division a = q b + r (b != 0).  Returns q (if q_out) and r (if
 * r_out); both fresh.  Reduction is delayed: each step adds (p - coef) b[j]
 * (< 2^32) to a slot, and a slot takes at most deg(b) + 1 such additions
 * before it is reduced (the leading slot when it becomes the pivot, the
 * rest at the end), so with deg(b) < 2^31 nothing overflows a uint64_t. */
static void np_divrem(const NPoly* a, const NPoly* b, uint64_t p,
                      NPoly** q_out, NPoly** r_out) {
    NPoly* r = np_copy(a);
    int db = b->deg;
    NPoly* q = np_alloc(a->deg >= db ? a->deg - db : -1);
    if (a->deg >= db) {
        uint64_t inv = nm_inv(b->c[db], p);
        for (int i = r->deg; i >= db; i--) {
            uint64_t lead = r->c[i] % p;
            if (!lead) continue;
            uint64_t coef = lead * inv % p;
            q->c[i - db] = coef;
            uint64_t neg = p - coef;
            uint64_t* rc = r->c + (i - db);
            for (int j = 0; j < db; j++) rc[j] += neg * b->c[j];
        }
        r->deg = db - 1;
        for (int j = 0; j <= r->deg; j++) r->c[j] %= p;
        np_trim(r);
        np_trim(q);
    }
    if (q_out) *q_out = q; else np_free(q);
    if (r_out) *r_out = r; else np_free(r);
}

static NPoly* np_rem(const NPoly* a, const NPoly* b, uint64_t p) {
    NPoly* r;
    np_divrem(a, b, p, NULL, &r);
    return r;
}

static NPoly* np_mulmod(const NPoly* a, const NPoly* b, const NPoly* f, uint64_t p) {
    NPoly* t = np_mul(a, b, p);
    NPoly* r = np_rem(t, f, p);
    np_free(t);
    return r;
}

static NPoly* np_powmod(const NPoly* a, uint64_t e, const NPoly* f, uint64_t p) {
    NPoly* res = np_const(1);
    NPoly* base = np_rem(a, f, p);
    while (e > 0) {
        if (e & 1) {
            NPoly* t = np_mulmod(res, base, f, p);
            np_free(res); res = t;
        }
        e >>= 1;
        if (e) {
            NPoly* t = np_mulmod(base, base, f, p);
            np_free(base); base = t;
        }
    }
    np_free(base);
    return res;
}

/* Monic gcd (zero when both inputs are zero). */
static NPoly* np_gcd(const NPoly* a, const NPoly* b, uint64_t p) {
    NPoly* x = np_copy(a);
    NPoly* y = np_copy(b);
    while (y->deg >= 0) {
        NPoly* r = np_rem(x, y, p);
        np_free(x);
        x = y;
        y = r;
    }
    np_free(y);
    np_make_monic(x, p);
    return x;
}

/* s a + t b = g (monic gcd), deg s < deg b, deg t < deg a. */
static void np_xgcd(const NPoly* a, const NPoly* b, uint64_t p,
                    NPoly** g_out, NPoly** s_out, NPoly** t_out) {
    NPoly *r0 = np_copy(a), *r1 = np_copy(b);
    NPoly *s0 = np_const(1), *s1 = np_alloc(-1);
    NPoly *t0 = np_alloc(-1), *t1 = np_const(1);
    while (r1->deg >= 0) {
        NPoly *q, *r2;
        np_divrem(r0, r1, p, &q, &r2);
        NPoly* qs = np_mul(q, s1, p);
        NPoly* s2 = np_sub(s0, qs, p);
        NPoly* qt = np_mul(q, t1, p);
        NPoly* t2 = np_sub(t0, qt, p);
        np_free(q); np_free(qs); np_free(qt);
        np_free(r0); r0 = r1; r1 = r2;
        np_free(s0); s0 = s1; s1 = s2;
        np_free(t0); t0 = t1; t1 = t2;
    }
    if (r0->deg >= 0) {
        uint64_t inv = nm_inv(r0->c[r0->deg], p);
        for (int i = 0; i <= r0->deg; i++) r0->c[i] = r0->c[i] * inv % p;
        for (int i = 0; i <= s0->deg; i++) s0->c[i] = s0->c[i] * inv % p;
        for (int i = 0; i <= t0->deg; i++) t0->c[i] = t0->c[i] * inv % p;
    }
    np_free(r1); np_free(s1); np_free(t1);
    *g_out = r0; *s_out = s0; *t_out = t0;
}

static NPoly* np_deriv(const NPoly* a, uint64_t p) {
    if (a->deg <= 0) return np_alloc(-1);
    NPoly* r = np_alloc(a->deg - 1);
    for (int i = 1; i <= a->deg; i++) r->c[i - 1] = a->c[i] * ((uint64_t)i % p) % p;
    np_trim(r);
    return r;
}

static NPoly* np_from_z(const ZUPoly* f, uint64_t p) {
    NPoly* r = np_alloc(f->deg);
    for (int i = 0; i <= f->deg; i++) r->c[i] = mpz_fdiv_ui(f->c[i], (unsigned long)p);
    np_trim(r);
    return r;
}

static ZUPoly* np_to_z(const NPoly* a) {
    ZUPoly* r = zupoly_new(a->deg >= 0 ? a->deg + 1 : 1);
    for (int i = 0; i <= a->deg; i++) mpz_set_ui(r->c[i], (unsigned long)a->c[i]);
    zupoly_normalize(r);
    return r;
}

/* Deterministic xorshift64* stream for the EDF random splitters. */
static uint64_t bz_rng_next(uint64_t* s) {
    uint64_t x = *s;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    *s = x;
    return x * UINT64_C(2685821657736338717);
}

/* ---------------------------------------------------------------------- */
/*  Distinct- and equal-degree factorisation mod p                        */
/* ---------------------------------------------------------------------- */

typedef struct {
    NPoly** g;   /* g[k] = product of the irreducible factors of degree d[k] */
    int* d;
    int n, cap;
} BzDdf;

static void bz_ddf_push(BzDdf* L, NPoly* g, int d) {
    if (L->n == L->cap) {
        L->cap = L->cap ? 2 * L->cap : 8;
        L->g = (NPoly**)realloc(L->g, sizeof(NPoly*) * (size_t)L->cap);
        L->d = (int*)realloc(L->d, sizeof(int) * (size_t)L->cap);
    }
    L->g[L->n] = g;
    L->d[L->n] = d;
    L->n++;
}

static void bz_ddf_clear(BzDdf* L) {
    for (int i = 0; i < L->n; i++) np_free(L->g[i]);
    free(L->g); free(L->d);
    L->g = NULL; L->d = NULL; L->n = L->cap = 0;
}

/* Number of irreducible factors recorded by a DDF. */
static int bz_ddf_count(const BzDdf* L) {
    int r = 0;
    for (int i = 0; i < L->n; i++) r += L->g[i]->deg / L->d[i];
    return r;
}

/* Distinct-degree factorisation of the monic squarefree f mod p.
 * Returns false (with L cleared) if the TimeConstrained deadline passed. */
static bool bz_ddf(const NPoly* f, uint64_t p, BzDdf* L) {
    NPoly* fs = np_copy(f);
    NPoly* X = np_alloc(1); X->c[1] = 1;
    NPoly* h = np_powmod(X, p, fs, p);          /* x^(p^i) mod fs, i = 1 */
    int i = 1;
    bool ok = true;
    while (fs->deg >= 2 * i) {
        if (tc_deadline_passed()) { ok = false; break; }
        NPoly* hx = np_sub(h, X, p);
        NPoly* g = np_gcd(fs, hx, p);
        np_free(hx);
        if (g->deg > 0) {
            NPoly* q;
            np_divrem(fs, g, p, &q, NULL);
            np_free(fs); fs = q;
            bz_ddf_push(L, g, i);
            NPoly* hr = np_rem(h, fs, p);
            np_free(h); h = hr;
        } else {
            np_free(g);
        }
        if (fs->deg < 2 * (i + 1)) break;
        NPoly* hn = np_powmod(h, p, fs, p);
        np_free(h); h = hn;
        i++;
    }
    if (ok && fs->deg > 0) bz_ddf_push(L, np_copy(fs), fs->deg);
    np_free(fs); np_free(X); np_free(h);
    if (!ok) bz_ddf_clear(L);
    return ok;
}

/* Equal-degree factorisation: f monic, product of deg(f)/d distinct
 * irreducibles of degree d mod the odd prime p.  Appends the monic
 * irreducible factors to out[*nout]. */
static void bz_edf(const NPoly* f, int d, uint64_t p, uint64_t* rng,
                   NPoly*** out, int* nout, int* cap) {
    if (f->deg <= d) {
        if (*nout == *cap) {
            *cap = *cap ? 2 * *cap : 16;
            *out = (NPoly**)realloc(*out, sizeof(NPoly*) * (size_t)*cap);
        }
        (*out)[(*nout)++] = np_copy(f);
        return;
    }
    for (;;) {
        /* Random a with deg a < deg f. */
        NPoly* a = np_alloc(f->deg - 1);
        for (int i = 0; i < f->deg; i++) a->c[i] = bz_rng_next(rng) % p;
        np_trim(a);
        if (a->deg <= 0) { np_free(a); continue; }
        /* b = a^((p^d - 1)/2) = (a^(1 + p + ... + p^(d-1)))^((p-1)/2):
         * the norm form keeps every exponent below p, so nothing like
         * p^d is ever formed. */
        NPoly* acc = np_copy(a);
        NPoly* t = np_copy(a);
        for (int k = 1; k < d; k++) {
            NPoly* tn = np_powmod(t, p, f, p);
            np_free(t); t = tn;
            NPoly* an = np_mulmod(acc, t, f, p);
            np_free(acc); acc = an;
        }
        np_free(t);
        NPoly* b = np_powmod(acc, (p - 1) / 2, f, p);
        np_free(acc);
        np_free(a);
        NPoly* one = np_const(1);
        NPoly* bm1 = np_sub(b, one, p);
        np_free(one); np_free(b);
        NPoly* g = np_gcd(f, bm1, p);
        np_free(bm1);
        if (g->deg > 0 && g->deg < f->deg) {
            NPoly* q;
            np_divrem(f, g, p, &q, NULL);
            np_make_monic(q, p);
            bz_edf(g, d, p, rng, out, nout, cap);
            bz_edf(q, d, p, rng, out, nout, cap);
            np_free(g); np_free(q);
            return;
        }
        np_free(g);
    }
}

/* ---------------------------------------------------------------------- */
/*  Arithmetic in (Z/mZ)[x] on GMP polynomials                            */
/* ---------------------------------------------------------------------- */

/* Reduce every coefficient into [0, m). */
static void zp_mod(ZUPoly* a, const mpz_t m) {
    for (int i = 0; i <= a->deg; i++) mpz_fdiv_r(a->c[i], a->c[i], m);
    zupoly_normalize(a);
}

static ZUPoly* zp_mul(const ZUPoly* a, const ZUPoly* b, const mpz_t m) {
    ZUPoly* r = zupoly_mul(a, b);
    zp_mod(r, m);
    return r;
}

static ZUPoly* zp_sub(const ZUPoly* a, const ZUPoly* b, const mpz_t m) {
    ZUPoly* r = zupoly_sub(a, b);
    zp_mod(r, m);
    return r;
}

static ZUPoly* zp_add(const ZUPoly* a, const ZUPoly* b, const mpz_t m) {
    ZUPoly* r = zupoly_add(a, b);
    zp_mod(r, m);
    return r;
}

/* a = q h + r mod m with h MONIC (lc exactly 1). */
static void zp_divrem_monic(const ZUPoly* a, const ZUPoly* h, const mpz_t m,
                            ZUPoly** q_out, ZUPoly** r_out) {
    ZUPoly* r = zupoly_copy(a);
    zp_mod(r, m);
    int dh = h->deg;
    int dq = r->deg - dh;
    ZUPoly* q = zupoly_new(dq >= 0 ? dq + 1 : 1);
    if (dq >= 0) {
        mpz_t coef; mpz_init(coef);
        for (int i = r->deg; i >= dh; i--) {
            mpz_fdiv_r(coef, r->c[i], m);
            if (mpz_sgn(coef) == 0) continue;
            mpz_set(q->c[i - dh], coef);
            for (int j = 0; j <= dh; j++) mpz_submul(r->c[i - dh + j], coef, h->c[j]);
        }
        mpz_clear(coef);
        for (int i = dh; i < r->cap; i++) mpz_set_ui(r->c[i], 0);
        zp_mod(r, m);
        zupoly_normalize(q);
    }
    if (q_out) *q_out = q; else zupoly_free(q);
    if (r_out) *r_out = r; else zupoly_free(r);
}

/* One quadratic Hensel step (von zur Gathen-Gerhard Alg. 15.10):
 * from f = g h mod m, s g + t h = 1 mod m (h monic, deg s < deg h,
 * deg t < deg g) produce the same relations mod m2 = m^2.  s, t are
 * updated only when `want_st`. */
static void zp_hensel_step(const ZUPoly* f, ZUPoly** g, ZUPoly** h,
                           ZUPoly** s, ZUPoly** t, const mpz_t m2, bool want_st) {
    ZUPoly* gh = zp_mul(*g, *h, m2);
    ZUPoly* e = zp_sub(f, gh, m2);
    zupoly_free(gh);

    ZUPoly* se = zp_mul(*s, e, m2);
    ZUPoly *q, *r;
    zp_divrem_monic(se, *h, m2, &q, &r);
    zupoly_free(se);

    ZUPoly* te = zp_mul(*t, e, m2);
    ZUPoly* qg = zp_mul(q, *g, m2);
    ZUPoly* g1 = zp_add(*g, te, m2);
    ZUPoly* gs = zp_add(g1, qg, m2);
    zupoly_free(g1); zupoly_free(te); zupoly_free(qg); zupoly_free(q);
    ZUPoly* hs = zp_add(*h, r, m2);
    zupoly_free(r); zupoly_free(e);

    if (want_st) {
        ZUPoly* sg = zp_mul(*s, gs, m2);
        ZUPoly* th = zp_mul(*t, hs, m2);
        ZUPoly* b0 = zp_add(sg, th, m2);
        zupoly_free(sg); zupoly_free(th);
        ZUPoly* one = zupoly_from_int(1);
        ZUPoly* b = zp_sub(b0, one, m2);
        zupoly_free(one); zupoly_free(b0);
        ZUPoly* sb = zp_mul(*s, b, m2);
        ZUPoly *c, *d;
        zp_divrem_monic(sb, hs, m2, &c, &d);
        zupoly_free(sb);
        ZUPoly* s_new = zp_sub(*s, d, m2);
        ZUPoly* tb = zp_mul(*t, b, m2);
        ZUPoly* cg = zp_mul(c, gs, m2);
        ZUPoly* t1 = zp_sub(*t, tb, m2);
        ZUPoly* t_new = zp_sub(t1, cg, m2);
        zupoly_free(tb); zupoly_free(cg); zupoly_free(t1);
        zupoly_free(c); zupoly_free(d); zupoly_free(b);
        zupoly_free(*s); *s = s_new;
        zupoly_free(*t); *t = t_new;
    }
    zupoly_free(*g); *g = gs;
    zupoly_free(*h); *h = hs;
}

/* Lift the factorisation F = lc(F) * u[0] * ... * u[r-1] mod p (u monic,
 * pairwise coprime) to monic factors mod M = p^(2^steps), written into
 * out[0..r-1].  F is any integer polynomial congruent to the lifted
 * product mod M (the root call passes the input, the children their
 * lifted cofactor), with lc(F) a unit mod p.  Balanced factor tree: each
 * level lifts a total degree of deg F, so the cost is O(log r) lifts. */
static void bz_hensel_tree(const ZUPoly* F, NPoly** u, int r, uint64_t p,
                           int steps, const mpz_t M, ZUPoly** out) {
    if (r == 1) {
        ZUPoly* m = zupoly_copy(F);
        zp_mod(m, M);
        mpz_t inv; mpz_init(inv);
        mpz_invert(inv, m->c[m->deg], M);
        for (int i = 0; i <= m->deg; i++) {
            mpz_mul(m->c[i], m->c[i], inv);
            mpz_fdiv_r(m->c[i], m->c[i], M);
        }
        mpz_clear(inv);
        zupoly_normalize(m);
        out[0] = m;
        return;
    }
    int k = r / 2;
    uint64_t lcp = mpz_fdiv_ui(F->c[F->deg], (unsigned long)p);
    NPoly* g0 = np_const(lcp);
    for (int i = 0; i < k; i++) { NPoly* t = np_mul(g0, u[i], p); np_free(g0); g0 = t; }
    NPoly* h0 = np_const(1);
    for (int i = k; i < r; i++) { NPoly* t = np_mul(h0, u[i], p); np_free(h0); h0 = t; }
    NPoly *gg, *s0, *t0;
    np_xgcd(g0, h0, p, &gg, &s0, &t0);    /* s0 g0 + t0 h0 = 1 */
    np_free(gg);

    ZUPoly* g = np_to_z(g0);
    ZUPoly* h = np_to_z(h0);
    ZUPoly* s = np_to_z(s0);
    ZUPoly* t = np_to_z(t0);
    np_free(g0); np_free(h0); np_free(s0); np_free(t0);

    mpz_t m; mpz_init_set_ui(m, (unsigned long)p);
    for (int j = 0; j < steps; j++) {
        mpz_mul(m, m, m);
        zp_hensel_step(F, &g, &h, &s, &t, m, j + 1 < steps);
    }
    mpz_clear(m);
    zupoly_free(s); zupoly_free(t);

    bz_hensel_tree(g, u, k, p, steps, M, out);
    bz_hensel_tree(h, u + k, r - k, p, steps, M, out + k);
    zupoly_free(g); zupoly_free(h);
}

/* ---------------------------------------------------------------------- */
/*  Zassenhaus over Z for one primitive squarefree polynomial             */
/* ---------------------------------------------------------------------- */

/* Upper limit on the subsets the recombination examines (after the
 * degree filter) before it stops and reports the unresolved cofactor as
 * one factor.  The result stays exact; only its completeness is capped. */
#define BZ_RECOMBINE_BUDGET 1000000L
/* ... and on the subsets merely enumerated (most rejected by the degree
 * filter at a cost of s additions each). */
#define BZ_VISIT_BUDGET 50000000L

/* Largest squarefree degree sent through the modular pipeline. */
#define BZ_MAX_DEG 1000

static void bz_out_push(ZUPoly*** out, int* n, int* cap, ZUPoly* g) {
    if (*n == *cap) {
        *cap = *cap ? 2 * *cap : 8;
        *out = (ZUPoly**)realloc(*out, sizeof(ZUPoly*) * (size_t)*cap);
    }
    (*out)[(*n)++] = g;
}

/* Subset-sum set of a DDF degree pattern: bits[k] = 1 iff some sub-product
 * of the modular factors has degree k. */
static void bz_degree_set(const BzDdf* L, int n, unsigned char* bits) {
    memset(bits, 0, (size_t)n + 1);
    bits[0] = 1;
    for (int i = 0; i < L->n; i++) {
        int cnt = L->g[i]->deg / L->d[i], d = L->d[i];
        for (int c = 0; c < cnt; c++)
            for (int k = n - d; k >= 0; k--)
                if (bits[k]) bits[k + d] = 1;
    }
}

/* Necessary condition for G | F over Z, in word arithmetic: G divides F
 * modulo two primes q near 2^16 that do not divide lc(G).  Rejects a wrong
 * candidate in O(n^2) word operations instead of a long division over Z
 * whose remainders grow without bound. */
static bool bz_divides_mod_q(const ZUPoly* F, const ZUPoly* G) {
    static const uint64_t qs[] = { 65521, 65519, 65497, 65479 };
    int used = 0;
    for (int k = 0; k < 4 && used < 2; k++) {
        uint64_t q = qs[k];
        if (mpz_divisible_ui_p(G->c[G->deg], (unsigned long)q)) continue;
        used++;
        NPoly* fq = np_from_z(F, q);
        NPoly* gq = np_from_z(G, q);
        NPoly* rq = np_rem(fq, gq, q);
        bool ok = (rq->deg < 0);
        np_free(fq); np_free(gq); np_free(rq);
        if (!ok) return false;
    }
    return true;
}

/* The next odd prime after p (trial division; p stays below 2^16). */
static int bz_next_prime(int p) {
    for (int q = p + 2;; q += 2) {
        bool pr = true;
        for (int d = 3; d * d <= q; d += 2) if (q % d == 0) { pr = false; break; }
        if (pr) return q;
    }
}

/* Factor the primitive, squarefree f (deg >= 1, lc > 0, f(0) != 0) into
 * irreducibles over Z, appended (primitive, positive lc) to out.  Returns
 * false only if the TimeConstrained deadline passed (out then holds an
 * exact but possibly coarser factorisation, still multiplying to f). */
static bool bz_zassenhaus(const ZUPoly* f, ZUPoly*** out, int* nout, int* cap) {
    int n = f->deg;
    /* Degree 1 is irreducible.  Above BZ_MAX_DEG the quadratic-arithmetic
     * distinct-degree factorisation (O(n^3 log p) word operations per
     * prime) is too slow to be useful, so the squarefree part is reported
     * whole: exact, possibly unsplit. */
    if (n <= 1 || n > BZ_MAX_DEG) { bz_out_push(out, nout, cap, zupoly_copy(f)); return true; }

    /* --- prime selection -------------------------------------------- */
    /* More primes buy a smaller modular factor count and a sharper degree
     * filter; DDF is cheap next to a recombination over too many factors,
     * until the degree makes each DDF itself expensive. */
    int want = n <= 300 ? 7 : (n <= 600 ? 4 : 2);
    unsigned char* allowed = (unsigned char*)malloc((size_t)n + 1);
    unsigned char* bits = (unsigned char*)malloc((size_t)n + 1);
    memset(allowed, 1, (size_t)n + 1);
    BzDdf best = {NULL, NULL, 0, 0};
    uint64_t best_p = 0;
    NPoly* best_f = NULL;
    int best_r = INT_MAX, good = 0;
    bool aborted = false;
    for (int p = 3; p < 65536 && good < want; p = bz_next_prime(p)) {
        if (mpz_divisible_ui_p(f->c[n], (unsigned long)p)) continue;
        NPoly* fb = np_from_z(f, (uint64_t)p);
        np_make_monic(fb, (uint64_t)p);
        NPoly* dfb = np_deriv(fb, (uint64_t)p);
        NPoly* gg = np_gcd(fb, dfb, (uint64_t)p);
        bool sqf = (gg->deg == 0);
        np_free(dfb); np_free(gg);
        if (!sqf) { np_free(fb); continue; }
        BzDdf L = {NULL, NULL, 0, 0};
        if (!bz_ddf(fb, (uint64_t)p, &L)) { np_free(fb); aborted = true; break; }
        good++;
        int r = bz_ddf_count(&L);
        bz_degree_set(&L, n, bits);
        for (int k = 0; k <= n; k++) allowed[k] &= bits[k];
        if (r < best_r) {
            bz_ddf_clear(&best); np_free(best_f);
            best = L; best_f = fb; best_r = r; best_p = (uint64_t)p;
        } else {
            bz_ddf_clear(&L); np_free(fb);
        }
        if (best_r <= 2) break;      /* 1: irreducible; 2: recombination is a single test */
    }
    free(bits);
    bool proper_degree = false;
    for (int k = 1; k < n; k++) if (allowed[k]) { proper_degree = true; break; }
    if (aborted || good == 0 || best_r <= 1 || !proper_degree) {
        /* Irreducible (one modular factor, or no admissible proper degree),
         * or no usable prime / deadline: report f whole -- exact. */
        bz_ddf_clear(&best); np_free(best_f); free(allowed);
        bz_out_push(out, nout, cap, zupoly_copy(f));
        return !aborted;
    }
    uint64_t p = best_p;

    /* --- equal-degree split mod p ------------------------------------ */
    NPoly** u = NULL; int r = 0, ucap = 0;
    uint64_t rng = UINT64_C(0x9E3779B97F4A7C15) ^ (uint64_t)n;
    for (int i = 0; i < best.n; i++) bz_edf(best.g[i], best.d[i], p, &rng, &u, &r, &ucap);
    bz_ddf_clear(&best); np_free(best_f);

    /* --- bound and Hensel lift --------------------------------------- */
    /* B = |lc| 2^(n-1) ceil(||f||_2); M = p^(2^steps) > 2B. */
    mpz_t B, M;
    mpz_init(B); mpz_init(M);
    for (int i = 0; i <= n; i++) mpz_addmul(B, f->c[i], f->c[i]);
    mpz_sqrt(B, B);
    mpz_add_ui(B, B, 1);
    mpz_mul_2exp(B, B, (mp_bitcnt_t)n);              /* 2 * 2^(n-1) */
    mpz_mul(B, B, f->c[n]);
    mpz_abs(B, B);
    int steps = 0;
    mpz_set_ui(M, (unsigned long)p);
    while (mpz_cmp(M, B) <= 0) { mpz_mul(M, M, M); steps++; }

    ZUPoly** L = (ZUPoly**)malloc(sizeof(ZUPoly*) * (size_t)r);
    bz_hensel_tree(f, u, r, p, steps, M, L);
    int* ldeg = (int*)malloc(sizeof(int) * (size_t)r);
    for (int i = 0; i < r; i++) { ldeg[i] = u[i]->deg; np_free(u[i]); }
    free(u);

    /* --- recombination ------------------------------------------------ */
    /* Cheap necessary conditions, all exact because M > 2B bounds every
     * candidate's 1-norm: for a = 0, 1, -1 the symmetric residue of
     * lc * prod u_i(a) mod M IS (lc/lc(G)) G(a), which must divide
     * lc * f*(a) whenever f*(a) != 0.  The a = 0 test alone is weak when
     * the roots are roots of unity (every constant term is then a unit
     * mod M, and +-1 turns up for a fraction 2/(p-1) of the subsets). */
    mpz_t halfM, lc, v;
    mpz_t tgt[3];
    mpz_init(halfM); mpz_init(lc); mpz_init(v);
    for (int a = 0; a < 3; a++) mpz_init(tgt[a]);
    mpz_fdiv_q_2exp(halfM, M, 1);
    mpz_t* ev = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)r * 3);  /* u_i(0), u_i(1), u_i(-1) mod M */
    for (int i = 0; i < r; i++) {
        mpz_init(ev[3 * i]); mpz_init(ev[3 * i + 1]); mpz_init(ev[3 * i + 2]);
        if (L[i]->deg >= 0) mpz_set(ev[3 * i], L[i]->c[0]);
        zupoly_eval_si(L[i], 1, ev[3 * i + 1]);
        mpz_fdiv_r(ev[3 * i + 1], ev[3 * i + 1], M);
        zupoly_eval_si(L[i], -1, ev[3 * i + 2]);
        mpz_fdiv_r(ev[3 * i + 2], ev[3 * i + 2], M);
    }
    ZUPoly* fs = zupoly_copy(f);
    int* act = (int*)malloc(sizeof(int) * (size_t)r);   /* live modular factors */
    int ra = r;
    for (int i = 0; i < r; i++) act[i] = i;
    int* comb = (int*)malloc(sizeof(int) * (size_t)(r + 1));
    long tested = 0, visited = 0;
    bool stop = false;

    for (int s = 1; 2 * s <= ra && !stop; ) {
        bool found = false;
        mpz_set(lc, fs->c[fs->deg]);
        mpz_mul(tgt[0], lc, fs->c[0]);                   /* lc * f*(a) */
        zupoly_eval_si(fs, 1, tgt[1]);  mpz_mul(tgt[1], tgt[1], lc);
        zupoly_eval_si(fs, -1, tgt[2]); mpz_mul(tgt[2], tgt[2], lc);
        for (int i = 0; i < s; i++) comb[i] = i;
        for (;;) {
            int dsum = 0;
            for (int i = 0; i < s; i++) dsum += ldeg[act[comb[i]]];
            if ((++visited & 4095) == 0 && tc_deadline_passed()) { stop = true; aborted = true; break; }
            if (visited > BZ_VISIT_BUDGET) { stop = true; break; }
            if (dsum < fs->deg && allowed[dsum]) {
                if (++tested > BZ_RECOMBINE_BUDGET) { stop = true; break; }
                bool pass = true;
                for (int a = 0; a < 3 && pass; a++) {
                    if (mpz_sgn(tgt[a]) == 0) continue;   /* f*(a) = 0: no information */
                    mpz_set(v, lc);
                    for (int i = 0; i < s; i++) {
                        mpz_mul(v, v, ev[3 * act[comb[i]] + a]);
                        mpz_fdiv_r(v, v, M);
                    }
                    if (mpz_cmp(v, halfM) > 0) mpz_sub(v, v, M);
                    pass = (mpz_sgn(v) != 0 && mpz_divisible_p(tgt[a], v));
                }
                if (pass) {
                    ZUPoly* G = zupoly_from_int(1);
                    mpz_set(G->c[0], lc);
                    for (int i = 0; i < s; i++) {
                        ZUPoly* t2 = zp_mul(G, L[act[comb[i]]], M);
                        zupoly_free(G); G = t2;
                    }
                    for (int i = 0; i <= G->deg; i++)
                        if (mpz_cmp(G->c[i], halfM) > 0) mpz_sub(G->c[i], G->c[i], M);
                    zupoly_normalize(G);
                    ZUPoly* Gp = zupoly_primitive_part(G);
                    zupoly_free(G);
                    ZUPoly* Q = (Gp->deg >= 1 && bz_divides_mod_q(fs, Gp))
                              ? zupoly_divexact(fs, Gp) : NULL;
                    if (Q) {
                        bz_out_push(out, nout, cap, Gp);
                        zupoly_free(fs); fs = Q;
                        /* drop the used modular factors from act[] */
                        int w = 0, ci = 0;
                        for (int i = 0; i < ra; i++) {
                            if (ci < s && comb[ci] == i) { ci++; continue; }
                            act[w++] = act[i];
                        }
                        ra = w;
                        found = true;
                        break;
                    }
                    zupoly_free(Gp);
                }
            }
            /* next combination of s out of ra */
            int i = s - 1;
            while (i >= 0 && comb[i] == ra - s + i) i--;
            if (i < 0) break;
            comb[i]++;
            for (int j = i + 1; j < s; j++) comb[j] = comb[j - 1] + 1;
        }
        if (!found) s++;
    }
    if (fs->deg >= 1) bz_out_push(out, nout, cap, fs); else zupoly_free(fs);

    for (int i = 0; i < r; i++) {
        zupoly_free(L[i]);
        mpz_clear(ev[3 * i]); mpz_clear(ev[3 * i + 1]); mpz_clear(ev[3 * i + 2]);
    }
    free(ev);
    for (int a = 0; a < 3; a++) mpz_clear(tgt[a]);
    free(L); free(ldeg); free(act); free(comb); free(allowed);
    mpz_clear(B); mpz_clear(M);
    mpz_clear(halfM); mpz_clear(lc); mpz_clear(v);
    return !aborted;
}

/* ---------------------------------------------------------------------- */
/*  Squarefree decomposition over Z (Yun)                                 */
/* ---------------------------------------------------------------------- */

static ZUPoly* bz_deriv(const ZUPoly* f) {
    ZUPoly* d = zupoly_new(f->deg > 0 ? f->deg : 1);
    for (int i = 1; i <= f->deg; i++) mpz_mul_ui(d->c[i - 1], f->c[i], (unsigned long)i);
    zupoly_normalize(d);
    return d;
}

/* True when f (primitive) is certified squarefree by a prime p that keeps
 * its degree and leaves f mod p squarefree.  False means "unknown". */
static bool bz_certify_squarefree(const ZUPoly* f) {
    int tried = 0;
    for (int p = 3; p < 65536 && tried < 3; p = bz_next_prime(p)) {
        if (mpz_divisible_ui_p(f->c[f->deg], (unsigned long)p)) continue;
        tried++;
        NPoly* fb = np_from_z(f, (uint64_t)p);
        NPoly* d = np_deriv(fb, (uint64_t)p);
        NPoly* g = np_gcd(fb, d, (uint64_t)p);
        bool ok = (d->deg >= 0 && g->deg == 0);
        np_free(fb); np_free(d); np_free(g);
        if (ok) return true;
    }
    return false;
}

/* Yun's algorithm on the primitive f with positive lc: f = prod parts[i]^mult[i],
 * every part primitive, squarefree, positive lc, pairwise coprime. */
static void bz_squarefree(const ZUPoly* f, ZUPoly*** parts, int** mults, int* n) {
    int cap = 4;
    *parts = (ZUPoly**)malloc(sizeof(ZUPoly*) * (size_t)cap);
    *mults = (int*)malloc(sizeof(int) * (size_t)cap);
    *n = 0;
    if (f->deg <= 1 || bz_certify_squarefree(f)) {
        (*parts)[0] = zupoly_copy(f); (*mults)[0] = 1; *n = 1;
        return;
    }
    ZUPoly* fp = bz_deriv(f);
    ZUPoly* a = zupoly_gcd(f, fp);
    if (a->deg <= 0) {
        zupoly_free(a); zupoly_free(fp);
        (*parts)[0] = zupoly_copy(f); (*mults)[0] = 1; *n = 1;
        return;
    }
    ZUPoly* b = zupoly_divexact(f, a);
    ZUPoly* c = zupoly_divexact(fp, a);
    zupoly_free(a); zupoly_free(fp);
    if (!b || !c) {            /* cannot happen for a true gcd; stay exact */
        zupoly_free(b); zupoly_free(c);
        (*parts)[0] = zupoly_copy(f); (*mults)[0] = 1; *n = 1;
        return;
    }
    ZUPoly* bd = bz_deriv(b);
    ZUPoly* d = zupoly_sub(c, bd);
    zupoly_free(bd); zupoly_free(c);
    for (int i = 1; b->deg >= 1; i++) {
        ZUPoly* ai = zupoly_gcd(b, d);
        ZUPoly* b2 = zupoly_divexact(b, ai);
        ZUPoly* c2 = zupoly_divexact(d, ai);
        if (!b2 || !c2) {      /* defensive: keep the product exact */
            zupoly_free(b2); zupoly_free(c2); zupoly_free(ai);
            break;
        }
        if (ai->deg >= 1) {
            if (*n == cap) {
                cap *= 2;
                *parts = (ZUPoly**)realloc(*parts, sizeof(ZUPoly*) * (size_t)cap);
                *mults = (int*)realloc(*mults, sizeof(int) * (size_t)cap);
            }
            (*parts)[*n] = ai; (*mults)[*n] = i; (*n)++;
        } else {
            zupoly_free(ai);
        }
        zupoly_free(b); zupoly_free(d);
        b = b2;
        ZUPoly* b2d = bz_deriv(b);
        d = zupoly_sub(c2, b2d);
        zupoly_free(b2d); zupoly_free(c2);
    }
    zupoly_free(b); zupoly_free(d);

    /* Verify prod parts^mult == f; on any mismatch fall back to f itself
     * (exact, merely unsplit). */
    ZUPoly* prod = zupoly_from_int(1);
    for (int k = 0; k < *n; k++)
        for (int e = 0; e < (*mults)[k]; e++) {
            ZUPoly* t = zupoly_mul(prod, (*parts)[k]);
            zupoly_free(prod); prod = t;
        }
    bool ok = zupoly_eq(prod, f);
    zupoly_free(prod);
    if (!ok) {
        for (int k = 0; k < *n; k++) zupoly_free((*parts)[k]);
        (*parts)[0] = zupoly_copy(f); (*mults)[0] = 1; *n = 1;
    }
}

/* ---------------------------------------------------------------------- */
/*  x^n - 1 and x^n + 1: cyclotomic factorisation                         */
/* ---------------------------------------------------------------------- */

/* x^n - 1 = prod_{d | n} Phi_d and x^n + 1 = prod_{d | 2n, d not | n} Phi_d,
 * every Phi_d irreducible over Q.  These are the inputs on which Zassenhaus
 * recombination is exponential (Phi_d splits into phi(d)/ord_d(p) factors
 * mod every p), so they are answered directly.  f must be primitive with
 * lc > 0; returns false (out untouched) unless f is exactly x^n -+ 1. */
#define BZ_CYCLO_MAX_DEG 20000
static bool bz_cyclotomic_split(const ZUPoly* f, ZUPoly*** out, int* nout, int* cap) {
    int n = f->deg;
    if (n < 2 || n > BZ_CYCLO_MAX_DEG || mpz_cmp_ui(f->c[n], 1) != 0) return false;
    for (int i = 1; i < n; i++) if (mpz_sgn(f->c[i]) != 0) return false;
    int sgn;
    if (mpz_cmp_si(f->c[0], -1) == 0) sgn = -1;          /* x^n - 1 */
    else if (mpz_cmp_ui(f->c[0], 1) == 0) sgn = 1;       /* x^n + 1 */
    else return false;
    int N = (sgn < 0) ? n : 2 * n;
    /* phi[d] for d | N, built bottom-up: Phi_d = (x^d - 1) / prod_{e | d, e < d} Phi_e */
    ZUPoly** phi = (ZUPoly**)calloc((size_t)N + 1, sizeof(ZUPoly*));
    for (int d = 1; d <= N; d++) {
        if (N % d) continue;
        ZUPoly* q = zupoly_new(d + 1);
        mpz_set_ui(q->c[d], 1);
        mpz_set_si(q->c[0], -1);
        zupoly_normalize(q);
        for (int e = 1; e < d && q; e++) {
            if (d % e || !phi[e]) continue;
            ZUPoly* t = zupoly_divexact(q, phi[e]);
            zupoly_free(q);
            q = t;
        }
        if (!q) {                                       /* impossible; stay exact */
            for (int k = 0; k <= N; k++) zupoly_free(phi[k]);
            free(phi);
            return false;
        }
        phi[d] = q;
    }
    for (int d = 1; d <= N; d++) {
        if (!phi[d]) continue;
        if (sgn > 0 && n % d == 0) continue;            /* x^n + 1: d | 2n, d not | n */
        bz_out_push(out, nout, cap, phi[d]);
        phi[d] = NULL;
    }
    for (int k = 0; k <= N; k++) zupoly_free(phi[k]);
    free(phi);
    return true;
}

/* ---------------------------------------------------------------------- */
/*  Public entry                                                          */
/* ---------------------------------------------------------------------- */

/* Factor the univariate polynomial P in `var` over Z.  Returns a fresh
 * Expr: Times[content, var^v, f_1^e_1, ...] (evaluated), each f_i
 * primitive and irreducible with positive leading coefficient, and the
 * signed integer content in front.  When P is irreducible and primitive
 * with positive leading coefficient -- or not a polynomial in `var` with
 * integer coefficients -- a copy of P is returned unchanged. */
Expr* bz_factor_to_expr(Expr* P, Expr* var) {
    if (!P || P->type != EXPR_FUNCTION) return expr_copy(P);
    ZUPoly* zp = expr_to_zupoly(P, var);
    if (!zp) return expr_copy(P);
    if (zp->deg <= 0) { zupoly_free(zp); return expr_copy(P); }

    /* Signed content: P = cont * f, f primitive with lc > 0. */
    mpz_t cont; mpz_init(cont);
    zupoly_content(zp, cont);
    if (mpz_sgn(zp->c[zp->deg]) < 0) mpz_neg(cont, cont);
    ZUPoly* f = zupoly_new(zp->deg + 1);
    for (int i = 0; i <= zp->deg; i++) mpz_divexact(f->c[i], zp->c[i], cont);
    zupoly_normalize(f);
    zupoly_free(zp);

    /* Strip var^v. */
    int v = 0;
    while (v <= f->deg && mpz_sgn(f->c[v]) == 0) v++;
    if (v > 0) {
        ZUPoly* g = zupoly_new(f->deg - v + 1);
        for (int i = v; i <= f->deg; i++) mpz_set(g->c[i - v], f->c[i]);
        zupoly_normalize(g);
        zupoly_free(f);
        f = g;
    }

    /* Squarefree parts, then each factored. */
    ZUPoly** fac = NULL; int* fmul = NULL; int nf = 0, fcap = 0;
    {
        ZUPoly** cyc = NULL; int nc = 0, ccap = 0;
        if (bz_cyclotomic_split(f, &cyc, &nc, &ccap)) {
            fac = cyc; nf = fcap = nc;
            fmul = (int*)malloc(sizeof(int) * (size_t)(nc > 0 ? nc : 1));
            for (int i = 0; i < nc; i++) fmul[i] = 1;
            f->deg = 0;                                  /* consumed */
        }
    }
    if (f->deg >= 1) {
        ZUPoly** parts; int* mults; int np;
        bz_squarefree(f, &parts, &mults, &np);
        for (int k = 0; k < np; k++) {
            ZUPoly** sub = NULL; int ns = 0, scap = 0;
            bz_zassenhaus(parts[k], &sub, &ns, &scap);
            for (int j = 0; j < ns; j++) {
                if (nf == fcap) {
                    fcap = fcap ? 2 * fcap : 8;
                    fac = (ZUPoly**)realloc(fac, sizeof(ZUPoly*) * (size_t)fcap);
                    fmul = (int*)realloc(fmul, sizeof(int) * (size_t)fcap);
                }
                fac[nf] = sub[j]; fmul[nf] = mults[k]; nf++;
            }
            free(sub);
            zupoly_free(parts[k]);
        }
        free(parts); free(mults);
    }
    zupoly_free(f);

    bool trivial = (mpz_cmp_ui(cont, 1) == 0 && v == 0 && nf == 1 && fmul[0] == 1);
    if (trivial) {
        zupoly_free(fac[0]); free(fac); free(fmul); mpz_clear(cont);
        return expr_copy(P);
    }

    size_t nargs = (size_t)nf + 2, k = 0;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * nargs);
    if (mpz_cmp_ui(cont, 1) != 0) {
        args[k++] = mpz_fits_slong_p(cont) ? expr_new_integer((int64_t)mpz_get_si(cont))
                                           : expr_new_bigint_from_mpz(cont);
    }
    if (v > 0) {
        args[k++] = (v == 1) ? expr_copy(var)
            : expr_new_function(expr_new_symbol(SYM_Power),
                                (Expr*[]){ expr_copy(var), expr_new_integer(v) }, 2);
    }
    for (int i = 0; i < nf; i++) {
        Expr* base = zupoly_to_expr(fac[i], var);
        zupoly_free(fac[i]);
        args[k++] = (fmul[i] == 1) ? base
            : expr_new_function(expr_new_symbol(SYM_Power),
                                (Expr*[]){ base, expr_new_integer(fmul[i]) }, 2);
    }
    free(fac); free(fmul); mpz_clear(cont);
    Expr* res = (k == 1) ? args[0]
        : expr_new_function(expr_new_symbol(SYM_Times), args, k);
    free(args);
    return eval_and_free(res);
}
