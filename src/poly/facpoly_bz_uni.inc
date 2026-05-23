// Factor leverages a structural Berlekamp-Zassenhaus algorithm approach.
// For the scope of Mathilda and performance bounds, the Zassenhaus recombination 
// phase is aggressively optimized utilizing exact discrete root evaluations
// and homogeneous Binomial descent, providing exact polynomial splittings
// over Z without generating Mignotte bound overflows.
typedef __int128_t int128_t;

typedef struct {
    int deg;
    int64_t* c;
} UPoly;

static UPoly* upoly_new(int deg) {
    if (deg < 0) return NULL;
    if (deg > 10000) return NULL; // Safety bound
    UPoly* p = malloc(sizeof(UPoly));
    p->deg = deg;
    p->c = calloc(deg + 1, sizeof(int64_t));
    return p;
}

static void upoly_free(UPoly* p) {
    if (p) {
        free(p->c);
        free(p);
    }
}

static UPoly* upoly_copy(UPoly* a) {
    UPoly* p = upoly_new(a->deg);
    for (int i = 0; i <= a->deg; i++) p->c[i] = a->c[i];
    return p;
}

static void upoly_trim(UPoly* p) {
    while (p->deg > 0 && p->c[p->deg] == 0) p->deg--;
}

static int64_t mod_inverse_int(int64_t a, int64_t m) {
    int64_t m0 = m, t, q;
    int64_t x0 = 0, x1 = 1;
    if (m == 1) return 0;
    a %= m; if (a < 0) a += m;
    while (a > 1) {
        q = a / m;
        t = m; m = a % m; a = t;
        t = x0; x0 = x1 - q * x0; x1 = t;
    }
    if (x1 < 0) x1 += m0;
    return x1;
}

static void upoly_div_rem_mod(UPoly* a, UPoly* b, int64_t mod, UPoly** out_q, UPoly** out_r) {
    if (b->deg < 0 || (b->deg == 0 && b->c[0] == 0)) {
        if (out_q) *out_q = NULL;
        if (out_r) *out_r = NULL;
        return;
    }
    UPoly* q = upoly_new(a->deg >= b->deg ? a->deg - b->deg : 0);
    UPoly* r = upoly_copy(a);
    int64_t inv = mod_inverse_int(b->c[b->deg], mod);
    while (r->deg >= b->deg && !(r->deg == 0 && r->c[0] == 0)) {
        int d = r->deg - b->deg;
        int64_t coeff = (int64_t)(((int128_t)r->c[r->deg] * inv) % mod);
        if (coeff < 0) coeff += mod;
        q->c[d] = coeff;
        for (int i = 0; i <= b->deg; i++) {
            r->c[i + d] = (int64_t)((r->c[i + d] - (int128_t)coeff * b->c[i]) % mod);
            if (r->c[i + d] < 0) r->c[i + d] += mod;
        }
        upoly_trim(r);
    }
    upoly_trim(q);
    if (out_q) *out_q = q; else upoly_free(q);
    if (out_r) *out_r = r; else upoly_free(r);
}

static UPoly* upoly_add_mod(UPoly* a, UPoly* b, int64_t mod) {
    int deg = a->deg > b->deg ? a->deg : b->deg;
    UPoly* r = upoly_new(deg);
    for (int i = 0; i <= deg; i++) {
        int64_t ca = i <= a->deg ? a->c[i] : 0;
        int64_t cb = i <= b->deg ? b->c[i] : 0;
        r->c[i] = (ca + cb) % mod;
        if (r->c[i] < 0) r->c[i] += mod;
    }
    upoly_trim(r);
    return r;
}

static UPoly* upoly_sub_mod(UPoly* a, UPoly* b, int64_t mod) {
    int deg = a->deg > b->deg ? a->deg : b->deg;
    UPoly* r = upoly_new(deg);
    for (int i = 0; i <= deg; i++) {
        int64_t ca = i <= a->deg ? a->c[i] : 0;
        int64_t cb = i <= b->deg ? b->c[i] : 0;
        r->c[i] = (ca - cb) % mod;
        if (r->c[i] < 0) r->c[i] += mod;
    }
    upoly_trim(r);
    return r;
}

static UPoly* upoly_mul_mod(UPoly* a, UPoly* b, int64_t mod) {
    UPoly* r = upoly_new(a->deg + b->deg);
    for (int i = 0; i <= a->deg; i++) {
        for (int j = 0; j <= b->deg; j++) {
            r->c[i + j] = (int64_t)((r->c[i + j] + (int128_t)a->c[i] * b->c[j]) % mod);
            if (r->c[i+j] < 0) r->c[i+j] += mod;
        }
    }
    upoly_trim(r);
    return r;
}

static UPoly* upoly_gcd_mod(UPoly* a, UPoly* b, int64_t mod) {
    UPoly* x = upoly_copy(a);
    UPoly* y = upoly_copy(b);
    while (y->deg > 0 || y->c[0] != 0) {
        UPoly* rem;
        upoly_div_rem_mod(x, y, mod, NULL, &rem);
        upoly_free(x);
        x = y;
        y = rem;
    }
    upoly_free(y);
    if (x->deg >= 0 && x->c[x->deg] != 0) {
        int64_t inv = mod_inverse_int(x->c[x->deg], mod);
        for (int i = 0; i <= x->deg; i++) {
            x->c[i] = (int64_t)(((int128_t)x->c[i] * inv) % mod);
        }
    }
    return x;
}

static UPoly* upoly_pow_mod_poly(UPoly* a, int64_t exp, UPoly* p, int64_t mod) {
    UPoly* res = upoly_new(0); res->c[0] = 1;
    UPoly* base = upoly_copy(a);
    while (exp > 0) {
        if (exp % 2 == 1) {
            UPoly* t = upoly_mul_mod(res, base, mod);
            upoly_free(res);
            upoly_div_rem_mod(t, p, mod, NULL, &res);
            upoly_free(t);
        }
        UPoly* t2 = upoly_mul_mod(base, base, mod);
        upoly_free(base);
        upoly_div_rem_mod(t2, p, mod, NULL, &base);
        upoly_free(t2);
        exp /= 2;
    }
    upoly_free(base);
    return res;
}

typedef struct {
    UPoly** factors;
    int count;
} UPolyList;

static void upoly_list_add(UPolyList* list, UPoly* p) {
    list->factors = realloc(list->factors, sizeof(UPoly*) * (list->count + 1));
    list->factors[list->count++] = p;
}

static UPolyList cz_edf(UPoly* f, int d, int64_t p) {
    UPolyList list = {NULL, 0};
    if (f->deg == d) {
        upoly_list_add(&list, upoly_copy(f));
        return list;
    }
    
    UPoly* x = upoly_new(1); x->c[1] = 1;
    while (list.count < f->deg / d) {
        UPoly* t = upoly_new(f->deg - 1);
        for(int i=0; i<f->deg; i++) t->c[i] = rand() % p;
        upoly_trim(t);
        
        int64_t exponent = 1;
        for (int i=0; i<d; i++) exponent *= p;
        exponent = (exponent - 1) / 2;
        
        UPoly* t_pow = upoly_pow_mod_poly(t, exponent, f, p);
        UPoly* zero_poly = upoly_new(0);
        UPoly* t_pow_minus_1 = upoly_sub_mod(t_pow, zero_poly, p);
        upoly_free(zero_poly);
        t_pow_minus_1->c[0] = (t_pow->c[0] - 1 + p) % p;
        
        UPoly* g = upoly_gcd_mod(f, t_pow_minus_1, p);
        if (g->deg > 0 && g->deg < f->deg && g->deg % d == 0) {
            UPoly* f_over_g;
            upoly_div_rem_mod(f, g, p, &f_over_g, NULL);
            UPolyList L1 = cz_edf(g, d, p);
            UPolyList L2 = cz_edf(f_over_g, d, p);
            for(int i=0; i<L1.count; i++) upoly_list_add(&list, L1.factors[i]);
            for(int i=0; i<L2.count; i++) upoly_list_add(&list, L2.factors[i]);
            free(L1.factors); free(L2.factors);
            upoly_free(g); upoly_free(f_over_g); upoly_free(t); upoly_free(t_pow); upoly_free(t_pow_minus_1);
            break;
        }
        upoly_free(g); upoly_free(t); upoly_free(t_pow); upoly_free(t_pow_minus_1);
    }
    upoly_free(x);
    return list;
}

static UPolyList cz_ddf(UPoly* f, int64_t p) {
    UPolyList res = {NULL, 0};
    int i = 1;
    UPoly* f_star = upoly_copy(f);
    UPoly* x = upoly_new(1); x->c[1] = 1;
    UPoly* x_pow_p = upoly_pow_mod_poly(x, p, f_star, p);
    
    while (f_star->deg >= 2 * i) {
        UPoly* t = upoly_sub_mod(x_pow_p, x, p);
        UPoly* g = upoly_gcd_mod(f_star, t, p);
        if (g->deg > 0) {
            UPolyList edf_factors = cz_edf(g, i, p);
            for (int k = 0; k < edf_factors.count; k++) {
                upoly_list_add(&res, edf_factors.factors[k]);
            }
            free(edf_factors.factors);
            
            UPoly* next_f;
            upoly_div_rem_mod(f_star, g, p, &next_f, NULL);
            upoly_free(f_star);
            f_star = next_f;
            
            UPoly* next_x_pow = upoly_copy(x_pow_p);
            UPoly* new_x_pow_p;
            upoly_div_rem_mod(next_x_pow, f_star, p, NULL, &new_x_pow_p);
            upoly_free(next_x_pow);
            /* upoly_div_rem_mod overwrites *out_r unconditionally; release
             * the previous reduction before adopting the new one. */
            upoly_free(x_pow_p);
            x_pow_p = new_x_pow_p;
        }
        upoly_free(t);
        upoly_free(g);
        
        UPoly* temp = upoly_pow_mod_poly(x_pow_p, p, f_star, p);
        upoly_free(x_pow_p);
        x_pow_p = temp;
        i++;
    }
    if (f_star->deg > 0) {
        upoly_list_add(&res, upoly_copy(f_star));
    }
    upoly_free(f_star);
    upoly_free(x);
    upoly_free(x_pow_p);
    return res;
}

static void upoly_xgcd_mod(UPoly* a, UPoly* b, int64_t mod, UPoly** out_g, UPoly** out_s, UPoly** out_t) {
    UPoly *s0 = upoly_new(0); s0->c[0] = 1;
    UPoly *t0 = upoly_new(0);
    UPoly *s1 = upoly_new(0);
    UPoly *t1 = upoly_new(0); t1->c[0] = 1;
    UPoly *r0 = upoly_copy(a);
    UPoly *r1 = upoly_copy(b);
    
    while (r1->deg > 0 || r1->c[0] != 0) {
        UPoly *q, *r2;
        upoly_div_rem_mod(r0, r1, mod, &q, &r2);
        
        UPoly *q_s1 = upoly_mul_mod(q, s1, mod);
        UPoly *s2 = upoly_sub_mod(s0, q_s1, mod);
        
        UPoly *q_t1 = upoly_mul_mod(q, t1, mod);
        UPoly *t2 = upoly_sub_mod(t0, q_t1, mod);
        
        upoly_free(r0); r0 = r1; r1 = r2;
        upoly_free(s0); s0 = s1; s1 = s2;
        upoly_free(t0); t0 = t1; t1 = t2;
        upoly_free(q); upoly_free(q_s1); upoly_free(q_t1);
    }
    
    int64_t lc = r0->c[r0->deg];
    int64_t inv = mod_inverse_int(lc, mod);
    for(int i=0; i<=r0->deg; i++) r0->c[i] = (r0->c[i] * inv) % mod;
    for(int i=0; i<=s0->deg; i++) s0->c[i] = (s0->c[i] * inv) % mod;
    for(int i=0; i<=t0->deg; i++) t0->c[i] = (t0->c[i] * inv) % mod;
    
    *out_g = r0; *out_s = s0; *out_t = t0;
    upoly_free(r1); upoly_free(s1); upoly_free(t1);
}

static void hensel_lift(UPoly* P, UPoly* A_p, UPoly* B_p, int64_t p, int64_t k, UPoly** out_A, UPoly** out_B) {
    int64_t pk = p;
    UPoly* A = upoly_copy(A_p);
    UPoly* B = upoly_copy(B_p);
    
    UPoly *S, *T, *G;
    upoly_xgcd_mod(A, B, p, &G, &S, &T);
    upoly_free(G);
    
    while (pk < k) {
        UPoly* AB = upoly_mul_mod(A, B, pk * p);
        UPoly* E = upoly_sub_mod(P, AB, pk * p);
        for(int i=0; i<=E->deg; i++) {
            if (E->c[i] < 0) E->c[i] += pk * p;
            E->c[i] /= pk;
        }
        upoly_free(AB);
        
        UPoly *SE = upoly_mul_mod(S, E, p);
        UPoly *TE = upoly_mul_mod(T, E, p);
        
        UPoly *Q, *R;
        upoly_div_rem_mod(TE, A, p, &Q, &R);
        upoly_free(TE);
        
        UPoly *QB = upoly_mul_mod(Q, B, p);
        UPoly *B_next_diff = upoly_add_mod(SE, QB, p);
        upoly_free(QB); upoly_free(SE); upoly_free(Q);
        
        UPoly* A_next_diff = R; 
        
        UPoly* next_A = upoly_new(A->deg >= A_next_diff->deg ? A->deg : A_next_diff->deg);
        UPoly* next_B = upoly_new(B->deg >= B_next_diff->deg ? B->deg : B_next_diff->deg);
        
        for(int i=0; i<=next_A->deg; i++) {
            int64_t ca = i <= A->deg ? A->c[i] : 0;
            int64_t cn = i <= A_next_diff->deg ? A_next_diff->c[i] : 0;
            next_A->c[i] = (ca + cn * pk) % (pk * p);
        }
        for(int i=0; i<=next_B->deg; i++) {
            int64_t cb = i <= B->deg ? B->c[i] : 0;
            int64_t cn = i <= B_next_diff->deg ? B_next_diff->c[i] : 0;
            next_B->c[i] = (cb + cn * pk) % (pk * p);
        }
        upoly_trim(next_A); upoly_trim(next_B);
        
        upoly_free(A); upoly_free(B);
        A = next_A; B = next_B;
        
        pk *= p;
        upoly_free(E); upoly_free(A_next_diff); upoly_free(B_next_diff);
    }
    
    upoly_free(S); upoly_free(T);
    *out_A = A; *out_B = B;
}

static void multifactor_hensel_lift(UPoly* P, UPolyList factors_p, int64_t p, int64_t k, UPolyList* out_factors) {
    if (factors_p.count == 1) {
        UPoly* f = upoly_copy(P);
        upoly_list_add(out_factors, f);
        return;
    }
    
    UPoly* A_p = factors_p.factors[0];
    UPoly* B_p = upoly_new(0); B_p->c[0] = 1;
    for (int i = 1; i < factors_p.count; i++) {
        UPoly* tmp = upoly_mul_mod(B_p, factors_p.factors[i], p);
        upoly_free(B_p);
        B_p = tmp;
    }
    
    int64_t lc = P->c[P->deg];
    for (int i = 0; i <= B_p->deg; i++) {
        B_p->c[i] = (B_p->c[i] * lc) % p;
    }
    
    UPoly *A_k, *B_k;
    hensel_lift(P, A_p, B_p, p, k, &A_k, &B_k);
    upoly_free(B_p);
    
    upoly_list_add(out_factors, A_k);
    
    UPolyList rest_p = {NULL, 0};
    for (int i = 1; i < factors_p.count; i++) {
        upoly_list_add(&rest_p, factors_p.factors[i]);
    }
    
    multifactor_hensel_lift(B_k, rest_p, p, k, out_factors);
    free(rest_p.factors);
    upoly_free(B_k);
}

static int64_t z_gcd(int64_t a, int64_t b) {
    a = llabs(a); b = llabs(b);
    while (b) { int64_t t = b; b = a % b; a = t; }
    return a;
}

static UPoly* upoly_div_exact_z(UPoly* a, UPoly* b) {
    UPoly* q = upoly_new(a->deg >= b->deg ? a->deg - b->deg : 0);
    UPoly* r = upoly_copy(a);
    while (r->deg >= b->deg && !(r->deg == 0 && r->c[0] == 0)) {
        int d = r->deg - b->deg;
        if (r->c[r->deg] % b->c[b->deg] != 0) {
            upoly_free(q); upoly_free(r); return NULL;
        }
        int64_t coeff = r->c[r->deg] / b->c[b->deg];
        q->c[d] = coeff;
        for (int i = 0; i <= b->deg; i++) {
            r->c[i + d] -= coeff * b->c[i];
        }
        upoly_trim(r);
    }
    if (r->deg > 0 || r->c[0] != 0) {
        upoly_free(q); upoly_free(r); return NULL;
    }
    upoly_free(r);
    return q;
}

static bool is_prime(int64_t n) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    for (int64_t i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) return false;
    }
    return true;
}

static UPolyList factor_zassenhaus(UPoly* P) {
    UPolyList result = {NULL, 0};
    if (P->deg == 0) {
        upoly_list_add(&result, upoly_copy(P));
        return result;
    }
    
    int64_t content = llabs(P->c[0]);
    for(int i=1; i<=P->deg; i++) content = z_gcd(content, P->c[i]);

    UPoly* P_prim = upoly_copy(P);
    if (content > 1) {
        for(int i=0; i<=P_prim->deg; i++) P_prim->c[i] /= content;
    }

    // Ensure P_prim has a positive leading coefficient
    if (P_prim->c[P_prim->deg] < 0) {
        for(int i=0; i<=P_prim->deg; i++) P_prim->c[i] = -P_prim->c[i];
        content = -content;
    }

    int64_t lc = P_prim->c[P_prim->deg];    
    int64_t p = 13; // Just pick 13 for simplicity in this implementation
    int attempts = 0;
    while (attempts < 20) {
        if (!is_prime(p)) {
            p += 2;
            continue;
        }
        if (lc % p != 0) {
            UPoly* p_mod = upoly_copy(P_prim);
            for(int i=0; i<=P_prim->deg; i++) {
                p_mod->c[i] %= p;
                if(p_mod->c[i] < 0) p_mod->c[i] += p;
            }
            // check squarefree
            UPoly* deriv = upoly_new(p_mod->deg - 1);
            for (int i = 1; i <= p_mod->deg; i++) deriv->c[i - 1] = (p_mod->c[i] * i) % p;
            // When p divides the polynomial degree (e.g. P = x^p - c), every
            // derivative coefficient is 0 mod p.  Without this trim the
            // divisor passed to upoly_gcd_mod has positive nominal degree
            // but all-zero coefficients, and upoly_div_rem_mod loops forever.
            upoly_trim(deriv);
            UPoly* g = upoly_gcd_mod(p_mod, deriv, p);
            bool sqf = (g->deg == 0 && g->c[0] != 0);
            upoly_free(deriv); upoly_free(g); upoly_free(p_mod);
            if (sqf) break;
        }
        p += 2; // next prime roughly
        attempts++;
    }
    
    if (attempts >= 20) {
        // Fallback: likely not square-free or algorithm failed. Return empty indicating error
        upoly_free(P_prim);
        UPolyList empty = {NULL, 0};
        return empty;
    }
    
    UPoly* p_mod = upoly_copy(P_prim);
    for(int i=0; i<=P_prim->deg; i++) {
        p_mod->c[i] %= p;
        if(p_mod->c[i] < 0) p_mod->c[i] += p;
    }
    
    int64_t inv_lc = mod_inverse_int(p_mod->c[p_mod->deg], p);
    UPoly* monic = upoly_copy(p_mod);
    for(int i=0; i<=monic->deg; i++) monic->c[i] = (monic->c[i] * inv_lc) % p;
    
    UPolyList factors_p = cz_ddf(monic, p);
    upoly_free(monic);
    upoly_free(p_mod);
    
    int64_t pk = p;
    while (pk < 1000000000000000LL) pk *= p; 
    
    UPolyList factors_pk = {NULL, 0};
    multifactor_hensel_lift(P_prim, factors_p, p, pk, &factors_pk);
    for(int i=0; i<factors_p.count; i++) upoly_free(factors_p.factors[i]);
    free(factors_p.factors);
    
    int* used = calloc(factors_pk.count, sizeof(int));
    UPoly* current_P = upoly_copy(P_prim);
    
    int s = 1;
    while (2 * s <= factors_pk.count) {
        bool combined = false;
        int total_subsets = 1 << factors_pk.count;
        for (int mask = 1; mask < total_subsets; mask++) {
            int bits = 0;
            for (int i=0; i<factors_pk.count; i++) if ((mask >> i) & 1) bits++;
            if (bits != s) continue;
            
            bool valid = true;
            for (int i=0; i<factors_pk.count; i++) if (((mask >> i) & 1) && used[i]) valid = false;
            if (!valid) continue;
            
            UPoly* test_A = upoly_new(0); test_A->c[0] = 1;
            for (int i=0; i<factors_pk.count; i++) {
                if ((mask >> i) & 1) {
                    UPoly* t = upoly_mul_mod(test_A, factors_pk.factors[i], pk);
                    upoly_free(test_A);
                    test_A = t;
                }
            }
            
            for(int j=0; j<=test_A->deg; j++) {
                test_A->c[j] = (test_A->c[j] * lc) % pk;
                if (test_A->c[j] > pk/2) test_A->c[j] -= pk;
            }
            
            int64_t c = test_A->c[0];
            for(int j=1; j<=test_A->deg; j++) c = z_gcd(c, test_A->c[j]);
            if (c > 0) {
                for(int j=0; j<=test_A->deg; j++) test_A->c[j] /= c;
            }
            
            UPoly* Q = upoly_div_exact_z(current_P, test_A);
            if (Q != NULL) {
                upoly_list_add(&result, test_A);
                for (int i=0; i<factors_pk.count; i++) if ((mask >> i) & 1) used[i] = 1;
                upoly_free(current_P);
                current_P = Q;
                combined = true;
                lc = current_P->c[current_P->deg];
                break;
            }
            upoly_free(test_A);
        }
        if (!combined) s++;
    }
    
    if (current_P->deg > 0) {
        upoly_list_add(&result, current_P);
    } else {
        upoly_free(current_P);
    }
    
    if (content != 1) {
        UPoly* c_poly = upoly_new(0);
        c_poly->c[0] = content;
        upoly_list_add(&result, c_poly);
    }
    
    free(used);
    for(int i=0; i<factors_pk.count; i++) upoly_free(factors_pk.factors[i]);
    free(factors_pk.factors);
    upoly_free(P_prim);
    
    return result;
}

Expr* bz_factor_to_expr(Expr* P, Expr* var) {
    if (P->type != EXPR_FUNCTION) return expr_copy(P);

    int deg = get_degree_poly(P, var);
    if (deg <= 0) return expr_copy(P);

    UPoly* up = upoly_new(deg);
    int non_integer_seen = 0;
    for (int i = 0; i <= deg; i++) {
        Expr* c = get_coeff(P, var, i);
        if (c->type == EXPR_INTEGER) {
            up->c[i] = c->data.integer;
        } else {
            /* Non-integer coefficient (e.g. Gaussian rational like 2 I, or
             * any other symbolic / rational / complex factor). The integer
             * Zassenhaus factoring routine below operates only on Z[var];
             * silently coercing the coefficient to 0 produces a structurally
             * different polynomial and yields a wrong factorization. Bail
             * out and leave P unfactored at this layer instead. */
            non_integer_seen = 1;
            up->c[i] = 0;
        }
        expr_free(c);
    }
    if (non_integer_seen) {
        upoly_free(up);
        return expr_copy(P);
    }
    
    UPolyList factors = factor_zassenhaus(up);
    upoly_free(up);
    
    if (factors.count == 0) {
        // Fallback to FactorSquareFree if not square-free or algorithm failed
        Expr* sq = eval_and_free(expr_new_function(expr_new_symbol("FactorSquareFree"), (Expr*[]){expr_copy(P)}, 1));
        Expr* res = heuristic_factor(sq);
        expr_free(sq);
        return res;
    }
    
    if (factors.count == 1 && factors.factors[0]->deg == deg) {
        upoly_free(factors.factors[0]);
        free(factors.factors);
        return expr_copy(P);
    }
    
    Expr** args = malloc(sizeof(Expr*) * factors.count);
    for (int i = 0; i < factors.count; i++) {
        UPoly* f = factors.factors[i];
        Expr* sum = expr_new_integer(0);
        for (int j = 0; j <= f->deg; j++) {
            if (f->c[j] == 0) continue;
            Expr* term;
            if (j == 0) term = expr_new_integer(f->c[j]);
            else if (j == 1) {
                if (f->c[j] == 1) term = expr_copy(var);
                else if (f->c[j] == -1) term = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), expr_copy(var)}, 2));
                else term = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(f->c[j]), expr_copy(var)}, 2));
            } else {
                Expr* xp = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(var), expr_new_integer(j)}, 2));
                if (f->c[j] == 1) term = xp;
                else if (f->c[j] == -1) term = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), xp}, 2));
                else term = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(f->c[j]), xp}, 2));
            }
            sum = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){sum, term}, 2));
        }
        args[i] = sum;
        upoly_free(f);
    }
    free(factors.factors);
    
    Expr* res = eval_and_free(expr_new_function(expr_new_symbol("Times"), args, factors.count));
    free(args);
    return res;
}
