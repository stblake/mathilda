/* ---- Gaussian content ----------------------------------------------------
 *
 * A polynomial with Gaussian-integer coefficients, at least one of them not
 * real, has its content taken over Z[i], as in the reference:
 * Factor[(1 + I) a + 2 b] is (1 + I) (a + (1 - I) b), and when what is left
 * is purely imaginary the unit I goes into the content too:
 * Factor[4 I a - 4 I b] is (4 I) (a - b), Factor[I x^2 - I] is
 * I (-1 + x) (1 + x).  The content is the gcd's associate in the first
 * quadrant (re > 0, im >= 0); the remaining sign is the ordinary real-content
 * sign of the primitive part.  p = g (p/g) whatever p is, so the rewrite is
 * value-preserving; it declines (NULL) on anything that is not a sum of terms
 * with Gaussian-integer coefficients and no denominators, and when the
 * content is 1 (a real integer content is the classical path's). */
static bool gc_coeff(const Expr* t, mpz_t re, mpz_t im) {
    const Expr* c = t;
    if (t->type == EXPR_FUNCTION && t->data.function.head->type == EXPR_SYMBOL &&
        t->data.function.head->data.symbol.name == SYM_Times && t->data.function.arg_count >= 1) {
        for (size_t i = 1; i < t->data.function.arg_count; i++) {
            const Expr* f = t->data.function.args[i];
            if (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL &&
                f->data.function.head->data.symbol.name == SYM_Power && f->data.function.arg_count == 2) {
                const Expr* ex = f->data.function.args[1];
                if (!(ex->type == EXPR_INTEGER && ex->data.integer > 0)) return false;
            } else if (expr_is_numeric_like(f)) return false;
        }
        c = t->data.function.args[0];
    } else if (t->type == EXPR_FUNCTION && t->data.function.head->type == EXPR_SYMBOL &&
               t->data.function.head->data.symbol.name == SYM_Power) {
        const Expr* ex = t->data.function.args[1];
        if (!(t->data.function.arg_count == 2 && ex->type == EXPR_INTEGER && ex->data.integer > 0)) return false;
    }
    if (expr_is_integer_like(c)) {
        mpz_t v; expr_to_mpz(c, v); mpz_set(re, v); mpz_clear(v); mpz_set_ui(im, 0);
        return true;
    }
    if (c->type == EXPR_FUNCTION && c->data.function.head->type == EXPR_SYMBOL &&
        c->data.function.head->data.symbol.name == SYM_Complex && c->data.function.arg_count == 2 &&
        expr_is_integer_like(c->data.function.args[0]) && expr_is_integer_like(c->data.function.args[1])) {
        mpz_t v;
        expr_to_mpz(c->data.function.args[0], v); mpz_set(re, v); mpz_clear(v);
        expr_to_mpz(c->data.function.args[1], v); mpz_set(im, v); mpz_clear(v);
        return true;
    }
    if (expr_is_numeric_like(c)) return false;   /* a rational, real or radical coefficient */
    mpz_set_ui(re, 1); mpz_set_ui(im, 0);        /* a bare monomial */
    return true;
}

/* q = round(n / d) componentwise, d > 0 */
static void gc_round_div(mpz_t q, const mpz_t n, const mpz_t d) {
    mpz_t t; mpz_init(t);
    mpz_mul_2exp(t, n, 1); mpz_add(t, t, d);     /* floor((2n + d) / 2d) */
    mpz_mul_2exp(q, d, 1);
    mpz_fdiv_q(q, t, q);
    mpz_clear(t);
}

/* (ar + ai i) <- gcd over Z[i] of itself and (br + bi i); b is clobbered */
static void gc_gcd(mpz_t ar, mpz_t ai, mpz_t br, mpz_t bi) {
    mpz_t nr, ni, N, qr, qi, t;
    mpz_inits(nr, ni, N, qr, qi, t, NULL);
    while (mpz_sgn(br) || mpz_sgn(bi)) {
        /* q = round(a conj(b) / |b|^2) */
        mpz_mul(nr, ar, br); mpz_addmul(nr, ai, bi);
        mpz_mul(ni, ai, br); mpz_submul(ni, ar, bi);
        mpz_mul(N, br, br); mpz_addmul(N, bi, bi);
        gc_round_div(qr, nr, N); gc_round_div(qi, ni, N);
        /* r = a - q b; a <- b; b <- r */
        mpz_mul(t, qr, br); mpz_submul(t, qi, bi); mpz_sub(nr, ar, t);
        mpz_mul(t, qr, bi); mpz_addmul(t, qi, br); mpz_sub(ni, ai, t);
        mpz_swap(ar, br); mpz_swap(ai, bi);
        mpz_swap(br, nr); mpz_swap(bi, ni);
    }
    mpz_clears(nr, ni, N, qr, qi, t, NULL);
}

static Expr* gc_number(const mpz_t re, const mpz_t im) {
    Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(re));
    if (mpz_sgn(im) == 0) return r;
    return expr_new_function(expr_new_symbol(SYM_Complex),
        (Expr*[]){ r, expr_bigint_normalize(expr_new_bigint_from_mpz(im)) }, 2);
}

static Expr* factor_gaussian_content(const Expr* arg) {
    /* a sum, or a product that expands to one: Factor[-(16 I) (56 - 30 D +
     * 3 D^2) a b] otherwise went to the classical path expanded, and came
     * back 16 a b (-56 I + (30 I) D - (3 I) D^2) */
    if (arg->type != EXPR_FUNCTION || arg->data.function.head->type != EXPR_SYMBOL ||
        (arg->data.function.head->data.symbol.name != SYM_Plus &&
         arg->data.function.head->data.symbol.name != SYM_Times)) return NULL;
    /* A product with a non-real Gaussian-integer coefficient keeps it as its
     * content and factors the rest, as the reference does: Factor of
     * -(4 I) (a b - c d + ...) leaves the sum's orientation alone (taking the
     * content of the expansion, 4 I, would flip it). */
    if (arg->data.function.head->data.symbol.name == SYM_Times) {
        const Expr* c = arg->data.function.args[0];
        if (c->type == EXPR_FUNCTION && c->data.function.head->type == EXPR_SYMBOL &&
            c->data.function.head->data.symbol.name == SYM_Complex &&
            expr_is_integer_like(c->data.function.args[0]) &&
            expr_is_integer_like(c->data.function.args[1])) {
            size_t n = arg->data.function.arg_count;
            Expr** rest = malloc(sizeof(Expr*) * (n - 1));
            for (size_t i = 1; i < n; i++) rest[i - 1] = expr_copy(arg->data.function.args[i]);
            Expr* r = n - 1 == 1 ? rest[0] : expr_new_function(expr_new_symbol(SYM_Times), rest, n - 1);
            free(rest);
            Expr* fr = eval_and_free(expr_new_function(expr_new_symbol("Factor"), (Expr*[]){ r }, 1));
            return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_copy((Expr*)c), fr }, 2));
        }
    }
    Expr* e = expr_expand((Expr*)arg);
    Expr* out = NULL;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Plus) {
        size_t n = e->data.function.arg_count;
        mpz_t gr, gi, cr, ci;
        mpz_inits(gr, gi, cr, ci, NULL);
        bool ok = true, any_complex = false, all_imag = true;
        for (size_t i = 0; i < n && ok; i++) {
            ok = gc_coeff(e->data.function.args[i], cr, ci);
            if (!ok) break;
            if (mpz_sgn(ci)) any_complex = true;
            gc_gcd(gr, gi, cr, ci);
        }
        if (ok && any_complex && (mpz_sgn(gr) || mpz_sgn(gi))) {
            /* first-quadrant associate: multiply by i until re > 0, im >= 0 */
            for (int k = 0; k < 4 && !(mpz_sgn(gr) > 0 && mpz_sgn(gi) >= 0); k++) {
                mpz_neg(gr, gr); mpz_swap(gr, gi);   /* (r + i m) i = -m + i r */
            }
            /* purely imaginary remainder: every c / g has a zero real part */
            mpz_t N, nr; mpz_inits(N, nr, NULL);
            mpz_mul(N, gr, gr); mpz_addmul(N, gi, gi);
            for (size_t i = 0; i < n && all_imag; i++) {
                gc_coeff(e->data.function.args[i], cr, ci);
                mpz_mul(nr, cr, gr); mpz_addmul(nr, ci, gi);   /* Re(c conj(g)) */
                if (mpz_sgn(nr)) all_imag = false;
            }
            mpz_clears(N, nr, NULL);
            if (all_imag) { mpz_neg(gi, gi); mpz_swap(gr, gi); }   /* g <- g i */
            if (!(mpz_cmp_ui(gr, 1) == 0 && mpz_sgn(gi) == 0)) {
                Expr* g = gc_number(gr, gi);
                Expr* prim = eval_and_free(expr_new_function(expr_new_symbol("Expand"),
                    (Expr*[]){ eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Power),
                                       (Expr*[]){ expr_copy(g), expr_new_integer(-1) }, 2),
                                   expr_copy(e) }, 2)) }, 1));
                /* the sign FactorList gives the primitive part goes into the
                 * content, so -I a + I b is -I (a - b), not I (-a + b).  Only
                 * in at most two variables, where the reference's sign is the
                 * lex one FactorList uses; over more it follows no monomial
                 * order (it keeps (4 I) (-SP[t, z] LC[c, d, k, p] + ...)), and
                 * the input's own orientation is what it keeps there. */
                size_t nv = 0, vcap = 4;
                Expr** vs = malloc(sizeof(Expr*) * vcap);
                collect_variables(prim, &vs, &nv, &vcap);
                for (size_t k = 0; k < nv; k++) expr_free(vs[k]);
                free(vs);
                Expr* fl = nv > 2 ? expr_new_symbol("Null")
                    : eval_and_free(expr_new_function(expr_new_symbol("FactorList"),
                        (Expr*[]){ expr_copy(prim) }, 1));
                if (fl->type == EXPR_FUNCTION && fl->data.function.arg_count >= 1) {
                    const Expr* c0 = fl->data.function.args[0];
                    if (c0->type == EXPR_FUNCTION && c0->data.function.arg_count == 2) {
                        const Expr* c = c0->data.function.args[0];
                        bool neg = (c->type == EXPR_INTEGER && c->data.integer < 0) ||
                                   (c->type == EXPR_BIGINT && mpz_sgn(c->data.bigint) < 0);
                        if (neg) {
                            g = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                    (Expr*[]){ expr_new_integer(-1), g }, 2));
                            prim = eval_and_free(expr_new_function(expr_new_symbol("Expand"),
                                    (Expr*[]){ eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                        (Expr*[]){ expr_new_integer(-1), prim }, 2)) }, 1));
                        }
                    }
                }
                expr_free(fl);
                Expr* fp = eval_and_free(expr_new_function(expr_new_symbol("Factor"),
                    (Expr*[]){ prim }, 1));
                out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ g, fp }, 2));
            }
        }
        mpz_clears(gr, gi, cr, ci, NULL);
    }
    expr_free(e);
    return out;
}

#ifdef USE_FLINT
/* One side (numerator or denominator) of a reduced fraction through FLINT:
 * fmpz_poly for a univariate integer polynomial, fmpq_mpoly otherwise.  NULL
 * when neither applies (not a polynomial over Q in symbol generators). */
static Expr* factor_flint_side(Expr* p) {
    if (p->type == EXPR_INTEGER || p->type == EXPR_BIGINT) return NULL;
    Expr* r = flint_univariate_factor_auto(p);
    if (!r) r = flint_polynomial_factor(p);
    return r;
}
#endif

Expr* builtin_factor(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;

    /* x^(2 + e) - x^e factors as (-1 + x) x^e (1 + x): powers with a symbolic
     * exponent factor over their generators b^t (expsplit.c), as in the
     * reference, and are put back together afterwards. */
    {
        ExpSplitTab st = {0};
        Expr* split = expsplit_abstract(res->data.function.args[0], &st);
        if (split) {
            size_t n = res->data.function.arg_count;
            Expr** fa = malloc(n * sizeof(Expr*));
            fa[0] = split;
            for (size_t i = 1; i < n; i++) fa[i] = expr_copy(res->data.function.args[i]);
            Expr* f = eval_and_free(expr_new_function(expr_new_symbol(SYM_Factor), fa, n));
            free(fa);
            Expr* out = f ? expsplit_restore(f, &st) : NULL;
            expsplit_free(&st);
            if (out) return out;
        }
    }

    /* Modulus -> p asks for the factorisation in F_p[x], and it was being
     * IGNORED -- the option was simply not recognised here, so the polynomial
     * was factored over Q and the answer came back as if the option were absent:
     *
     *     Factor[x^2 + 1, Modulus -> 5]   gave  1 + x^2   (GF(5): (x+2)(x+3))
     *     Factor[x^4 + 1, Modulus -> 2]   gave  1 + x^4   (GF(2): (1+x)^4)
     *
     * A silent wrong answer, and FactorList inherited it verbatim since it
     * delegates here.  FLINT's nmod_poly_factor was already linked and in use;
     * flint_nmod_poly_factor_list exposes it.  When the fast path declines -- a
     * non-word-size or composite modulus, more than one variable, USE_FLINT=0 --
     * the whole call declines (NULL, head unevaluated) rather than answering
     * over Q, which is the one outcome that must not come back. */
    {
        size_t n = res->data.function.arg_count;
        Expr* last = res->data.function.args[n - 1];
        Expr* mod_p = NULL;
        if (n >= 2 && last->type == EXPR_FUNCTION
            && last->data.function.head->type == EXPR_SYMBOL
            && last->data.function.head->data.symbol.name == SYM_Rule
            && last->data.function.arg_count == 2
            && last->data.function.args[0]->type == EXPR_SYMBOL
            && last->data.function.args[0]->data.symbol.name == SYM_Modulus)
            mod_p = last->data.function.args[1];

        if (mod_p) {
            if (n != 2 || mod_p->type != EXPR_INTEGER || mod_p->data.integer < 2)
                return NULL;
            Expr** vars = NULL; size_t nv = 0, cap = 0;
            collect_variables(res->data.function.args[0], &vars, &nv, &cap);
            Expr* x = (nv == 1) ? expr_copy(vars[0]) : NULL;
            for (size_t i = 0; i < nv; i++) expr_free(vars[i]);
            free(vars);
            if (!x) return NULL;                       /* not univariate: decline */

            Expr* fl = flint_nmod_poly_factor_list(res->data.function.args[0], x,
                                                   (unsigned long)mod_p->data.integer);
            expr_free(x);
            if (!fl) return NULL;
            /* {{c,1},{f1,e1},...}  ->  c f1^e1 ... */
            size_t np = fl->data.function.arg_count;
            Expr** terms = malloc(np * sizeof(Expr*));
            if (!terms) { expr_free(fl); return NULL; }
            for (size_t i = 0; i < np; i++) {
                Expr* pr = fl->data.function.args[i];
                Expr* base = expr_copy(pr->data.function.args[0]);
                Expr* e    = pr->data.function.args[1];
                terms[i] = (e->type == EXPR_INTEGER && e->data.integer == 1)
                    ? base
                    : internal_power((Expr*[]){ base, expr_copy(e) }, 2);
            }
            expr_free(fl);
            Expr* prod = internal_times(terms, np);
            free(terms);
            return prod;
        }
    }

#ifdef USE_FLINT
    /* Parametric radical field Q(t..)(sqrt k), k a free symbol: factor over the
     * field via FLINT (sqrt(k) -> S collapse -> multivariate factoring over Q).
     * The classical algebraic-factoring path (qafactor/Trager) is built for
     * constant radicands and blows up / hangs on a symbolic radicand — this is
     * the Apart/Factor blocker in the parametric rational integrator. Only fires
     * when a Sqrt[symbol] is actually present, matching WL's default of factoring
     * over the field the input's radicals generate. */
    if (res->data.function.arg_count == 1) {
        Expr* pf = flint_parametric_sqrt_factor(res->data.function.args[0]);
        if (pf) return pf;
    }
#endif

    /* Phase G5: Extension -> α option.  Walk trailing args looking for
     * Rule[Extension, α] / RuleDelayed[Extension, α].  When present,
     * dispatch to the algebraic-factoring path in qafactor.c (Trager).
     *
     * Phase G9: when α is the symbol `Automatic`, run extension_autodetect
     * on the polynomial argument and use the resulting (tier-1) single
     * generator's surface form as α.
     *
     * The memo is bypassed for this branch: the cache key is the
     * polynomial alone, but the result depends on the extension too. */
    if (res->data.function.arg_count >= 2) {
        Expr* alpha_expr = NULL;
        for (size_t i = 1; i < res->data.function.arg_count; i++) {
            Expr* opt = res->data.function.args[i];
            if (opt->type == EXPR_FUNCTION
                && opt->data.function.head
                && opt->data.function.head->type == EXPR_SYMBOL
                && (opt->data.function.head->data.symbol.name == SYM_Rule
                    || opt->data.function.head->data.symbol.name == SYM_RuleDelayed)
                && opt->data.function.arg_count == 2) {
                Expr* lhs = opt->data.function.args[0];
                Expr* rhs = opt->data.function.args[1];
                if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == SYM_Extension) {
                    alpha_expr = rhs;
                }
            }
        }
        /* Auto-detect when the user passed Extension -> Automatic. */
        Expr* alpha_auto = NULL;
        QATower* auto_tower = NULL;
        if (alpha_expr && alpha_expr->type == EXPR_SYMBOL
            && alpha_expr->data.symbol.name
            && strcmp(alpha_expr->data.symbol.name, "Automatic") == 0) {
            auto_tower = extension_autodetect(res->data.function.args[0]);
            if (auto_tower && auto_tower->n == 1) {
                alpha_auto = expr_copy(auto_tower->alpha_renders[0]);
                alpha_expr = alpha_auto;
            } else {
                /* No detectable generator, or multi-generator tower.
                 * Tier-1: fall through to plain Factor. */
                alpha_expr = NULL;
            }
        }
        /* Skip the qa-factor path when alpha_expr resolves to a no-op
         * sentinel (None or Automatic-with-no-detect). */
        if (alpha_expr && alpha_expr->type == EXPR_SYMBOL
            && alpha_expr->data.symbol.name
            && strcmp(alpha_expr->data.symbol.name, "None") == 0) {
            alpha_expr = NULL;
        }
        if (alpha_expr) {
            Expr* poly = res->data.function.args[0];

            /* Phase G6: when α is a List[α_1, α_2, ...], dispatch
             * to the tower-of-extensions path.  Otherwise (G5) treat
             * it as a single algebraic generator. */
            Expr** alpha_list = NULL;
            int n_alphas = 0;
            bool is_tower = false;
            if (alpha_expr->type == EXPR_FUNCTION
                && alpha_expr->data.function.head
                && alpha_expr->data.function.head->type == EXPR_SYMBOL
                && alpha_expr->data.function.head->data.symbol.name == SYM_List) {
                is_tower = true;
                n_alphas = (int)alpha_expr->data.function.arg_count;
                alpha_list = alpha_expr->data.function.args;  /* borrowed */
            }

            /* Pick the polynomial variable.  Q(α)-factoring works in a
             * single indeterminate; we collect the free symbols of poly
             * and require exactly one (after stripping α's surface
             * form).  Multivariate algebraic factoring is out of MVP
             * scope. */
            Expr** vars = NULL;
            size_t vc = 0, vcap = 8;
            vars = (Expr**)malloc(sizeof(Expr*) * vcap);
            collect_variables(poly, &vars, &vc, &vcap);

            Expr* poly_var = NULL;
            size_t live = 0;
            for (size_t i = 0; i < vc; i++) {
                bool is_alpha = false;
                if (is_tower) {
                    for (int k = 0; k < n_alphas; k++) {
                        if (expr_eq(vars[i], alpha_list[k])) {
                            is_alpha = true;
                            break;
                        }
                    }
                } else if (expr_eq(vars[i], alpha_expr)) {
                    is_alpha = true;
                }
                if (is_alpha) continue;
                poly_var = vars[i];
                live++;
            }
            if (live == 1 && poly_var) {
                Expr* result;
                if (is_tower) {
                    result = qa_factor_with_extension_tower(
                        poly, alpha_list, n_alphas, poly_var);
                } else {
                    result = qa_factor_with_extension(poly, alpha_expr, poly_var);
                }
                for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
                free(vars);
                if (result) {
                    if (alpha_auto) expr_free(alpha_auto);
                    if (auto_tower) qa_tower_free(auto_tower);
                    return result;
                }
                /* Fall through to plain Factor on failure. */
            } else {
                for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
                free(vars);
            }
        }
        /* Cleanup auto-detect state before any fall-through to plain Factor. */
        if (alpha_auto) expr_free(alpha_auto);
        if (auto_tower) qa_tower_free(auto_tower);
    }

    /* Memo lookup: if a Simplify call has installed a Factor memo and
     * we've already factored this exact input, return the cached copy.
     * Note: builtin_factor follows the convention that the evaluator
     * owns `res` (we do not call expr_free(res) anywhere in this
     * function), so the hit path simply returns a fresh copy of the
     * cached value without touching `res`. */
    FactorMemo* memo = factor_memo_top();
    if (memo) {
        const Expr* hit = factor_memo_get(memo, res);
        if (hit) {
            return expr_copy((Expr*)hit);
        }
    }
    Expr* memo_key = memo ? expr_copy(res) : NULL;

    /* Berlekamp–Zassenhaus and the heuristic multivariate paths assume
     * integer (rational) coefficients; route inexact inputs through the
     * standard rationalise / numericalise round-trip.  We don't cache
     * inexact results -- they go through their own evaluation path. */
    if (internal_args_contain_inexact(res)) {
        if (memo_key) expr_free(memo_key);
        return internal_rationalize_then_numericalize(res, builtin_factor);
    }

    /* Thread over logic / comparison heads (Less, Greater, Equal,
     * LessEqual, GreaterEqual, Inequality, And, Or). Mathematica applies
     * Factor elementwise inside these structures, e.g.,
     *   Factor[1 < 1+2x+x^2+1/(1+x) < 2]
     * yields Inequality[1, Less, factored, Less, 2].  We thread by
     * recursively invoking Factor on each comparison argument, then
     * rebuilding the same head with the factored subterms.
     * ATTR_LISTABLE on Factor handles the List case automatically;
     * this hand-threading covers the comparison / logic heads it
     * doesn't reach.  For Inequality specifically, operator-symbol
     * slots (odd indices) are passed through verbatim. */
    {
        Expr* arg0 = res->data.function.args[0];
        if (arg0->type == EXPR_FUNCTION
            && arg0->data.function.head
            && arg0->data.function.head->type == EXPR_SYMBOL) {
            const char* h = arg0->data.function.head->data.symbol.name;
            if (h == SYM_Less ||
                h == SYM_Greater ||
                h == SYM_Equal ||
                h == SYM_Unequal ||
                h == SYM_LessEqual ||
                h == SYM_GreaterEqual ||
                h == SYM_Inequality ||
                h == SYM_And ||
                h == SYM_Or) {
                bool is_ineq = (h == SYM_Inequality);
                size_t n = arg0->data.function.arg_count;
                Expr** new_args = (Expr**)malloc(sizeof(Expr*) * n);
                for (size_t i = 0; i < n; i++) {
                    if (is_ineq && (i & 1u) == 1) {
                        /* Operator-symbol slot: copy through unchanged. */
                        new_args[i] = expr_copy(arg0->data.function.args[i]);
                        continue;
                    }
                    Expr* sub_call = expr_new_function(
                        expr_new_symbol("Factor"),
                        (Expr*[]){expr_copy(arg0->data.function.args[i])}, 1);
                    new_args[i] = evaluate(sub_call);
                    expr_free(sub_call);
                }
                Expr* result = expr_new_function(expr_copy(arg0->data.function.head),
                                                 new_args, n);
                free(new_args);
                Expr* evaluated = eval_and_free(result);
                if (memo_key) {
                    factor_memo_put(memo, memo_key, evaluated);
                    expr_free(memo_key);
                }
                return evaluated;
            }
        }
    }

    Expr* arg = res->data.function.args[0];

    /* Algebraic-generator (radical) pass.  When the input contains a
     * sub-expression u with a fractional rational exponent (u^(p/q),
     * q != 1) — e.g. r^(1/5), Log[r]^(2/3), (x+1)^(1/2) — substitute
     * u -> g^m where m is the lcm of all such q's, factor as a polynomial
     * in g, then back-substitute g -> u^(1/m).  The result is fed back
     * through evaluate() so that any residual fractional bases get
     * processed recursively (each level picks a fresh generator name to
     * avoid colliding with previously-introduced ones). */
    {
        Expr* base = NULL;
        Expr* atom = NULL;
        int64_t m = 1;
        if (poly_find_radical_gen(arg, &base, &atom, &m)) {
            char* gen = poly_make_fresh_gen(arg);
            Expr* substituted = poly_subst_radical_to_gen(arg, base, atom, m, gen);
            Expr* call = expr_new_function(expr_new_symbol("Factor"),
                          (Expr*[]){substituted}, 1);
            Expr* result_in_g = evaluate(call);
            expr_free(call);
            Expr* final = poly_subst_radical_from_gen(result_in_g, base, atom, m, gen);
            expr_free(result_in_g);
            expr_free(base);
            expr_free(atom);
            free(gen);
            if (memo_key) {
                factor_memo_put(memo, memo_key, final);
                expr_free(memo_key);
            }
            return final;
        }
    }

    /* Opaque-kernel pass.  Every non-symbol generator (BesselJ[1/3, x], Log[x],
     * f[x], E^x) is an independent indeterminate to the factoriser, so make it
     * a real one: factor with each kernel replaced by a fresh symbol, then
     * restore.  Without this, the multivariate paths specialise a symbol
     * generator by ReplaceAll (facpoly_z_split.inc's x -> 0, the
     * eval_others_at_alpha probes), which also rewrites the x INSIDE the
     * kernels -- not a specialisation of the polynomial, and for BesselJ[-2/3, x]
     * a 0*ComplexInfinity that leaked Infinity::indet out of Factor and
     * Simplify.  It also lets the FLINT multivariate fast path, which declines
     * on unknown heads, take these inputs.  Runs after the radical pass, whose
     * u^(p/q) generators may themselves be kernels; the recursive Factor call
     * sees symbols only, so this fires at most once per call. */
    {
        PolyKernelMap kmap;
        Expr* kerned = poly_kernels_to_symbols(arg, &kmap);
        if (kerned) {
            size_t n = res->data.function.arg_count;
            Expr** cargs = malloc(sizeof(Expr*) * n);
            cargs[0] = kerned;
            for (size_t i = 1; i < n; i++) cargs[i] = expr_copy(res->data.function.args[i]);
            Expr* in_syms = eval_and_free(expr_new_function(expr_new_symbol("Factor"), cargs, n));
            free(cargs);
            Expr* final = eval_and_free(poly_symbols_to_kernels(in_syms, &kmap));
            expr_free(in_syms);
            poly_kernel_map_free(&kmap);
            if (memo_key) {
                factor_memo_put(memo, memo_key, final);
                expr_free(memo_key);
            }
            return final;
        }
    }

    if (res->data.function.arg_count == 1) {
        Expr* gf = factor_gaussian_content(arg);
        if (gf) {
            if (memo_key) {
                factor_memo_put(memo, memo_key, gf);
                expr_free(memo_key);
            }
            return gf;
        }
    }

#ifdef USE_FLINT
    /* Univariate integer-polynomial fast path via FLINT (fmpz_poly_factor):
     * fast, and bignum-safe. It must run BEFORE the Together / collect_variables
     * machinery below, which misclassifies a bignum coefficient as a variable
     * (is_number omits EXPR_BIGINT) and routes a univariate bignum polynomial
     * into the multivariate path, dropping its constant term and returning it
     * unfactored. flint_univariate_factor_auto self-detects the single variable
     * with a bignum-safe scan and emits Mathematica's primitive-positive-leading
     * factor convention. Returns NULL for anything that is not a univariate
     * polynomial over Z (multivariate, rational/symbolic coefficients, a
     * denominator), which then takes the classical path unchanged. */
    {
        Expr* uf = flint_univariate_factor_auto(arg);
        if (uf) {
            if (memo_key) {
                factor_memo_put(memo, memo_key, uf);
                expr_free(memo_key);
            }
            return uf;
        }
    }

    /* Multivariate integer/rational-polynomial fast path via FLINT
     * (fmpq_mpoly_factor).  The classical heuristic multivariate factorer is
     * exponential in the number of variables and stalls badly on inputs like
     * x^99 - y^99 (>20 s), while FLINT returns in milliseconds.  Only fires for
     * the plain single-argument form (no Extension / GaussianIntegers / Modulus
     * option), so algebraic-field factoring is untouched.  flint_polynomial_factor
     * returns NULL for anything that is not a genuine multivariate polynomial
     * over Q (a denominator, a fractional or symbolic exponent, a non-polynomial
     * head such as Sqrt / Sin / Complex) — those fall through to the classical
     * path below unchanged.  The emitted form is FLINT's primitive, positive-
     * leading-coefficient factor convention with a separated rational content,
     * which matches Mathematica's (and the classical path's) output for
     * factorable polynomials; the sole systematic difference is that a globally
     * negative leading term is pulled out as an explicit -1 (e.g. y^2 - x^2 ->
     * -(x - y)(x + y)), which is Mathematica's convention and more canonical
     * than the classical (x + y)(-x + y).  The evaluator re-sorts the Orderless
     * Times afterwards, so factor ordering is independent of FLINT. */
    if (res->data.function.arg_count == 1) {
        Expr* mf = flint_polynomial_factor(arg);
        if (mf) {
            if (memo_key) {
                factor_memo_put(memo, memo_key, mf);
                expr_free(memo_key);
            }
            return mf;
        }
    }
#endif

    Expr* together = eval_and_free(expr_new_function(expr_new_symbol("Together"), (Expr*[]){expr_copy(arg)}, 1));
    Expr* n_call = expr_new_function(expr_new_symbol("Numerator"), (Expr*[]){expr_copy(together)}, 1);
    Expr* num = evaluate(n_call);
    expr_free(n_call);

    Expr* d_call = expr_new_function(expr_new_symbol("Denominator"), (Expr*[]){expr_copy(together)}, 1);
    Expr* den = evaluate(d_call);
    expr_free(d_call);

    /* Correctness guard: the polynomial factorer below assumes num and den are
     * genuine polynomials (Together, above, is contracted to return a single
     * reduced fraction).  When Together declines to combine -- e.g. it bails on
     * a sum of fractions whose numerators share >=2 dependent trig/hyperbolic
     * kernels over Q(i), to avoid a GCD blow-up -- Numerator returns the whole
     * uncombined sum, which still carries Power[linear, -1] denominators.
     * Factoring that non-polynomial treats 1/(x-I) as an indeterminate and
     * returns a NON-EQUIVALENT result (Factor[(Cosh[x]+Sinh[x])/(-I+x) +
     * (Cosh[x]+Sinh[x])/(I+x)] collapsed to 0, which then propagated a wrong
     * Simplify answer).  Factor must be equivalence-preserving, so decline
     * cleanly: return the input unchanged rather than a false factorisation. */
    if (factor_has_leftover_denominator(num) ||
        factor_has_leftover_denominator(den)) {
        Expr* out = expr_copy(arg);
        expr_free(together);
        expr_free(num);
        expr_free(den);
        if (memo_key) {
            factor_memo_put(memo, memo_key, out);
            expr_free(memo_key);
        }
        return out;
    }

    /* Scope discipline: when called from inside Simplify (where a
     * Factor memo is active), use the SHARED-scope behavior -- both
     * numerator and denominator use the variable list collected from
     * the numerator.  This is conservative: a denominator with extra
     * variables won't have those variables in its factoring scope and
     * may end up unfactored.  But it preserves the historical
     * Simplify behavior, which has the property that some downstream
     * TrigRoundtrip paths only converge when intermediate factors
     * remain in their non-fully-factored form.  See plans/FACTOR_PLAN.md
     * cleanup C4.
     *
     * For direct user Factor calls (no memo active), use the SEPARATE-
     * scope behavior: each of num and den is factored with its own
     * variable list.  This gives the correct fully-factored result
     * for inputs like Factor[(x^3+2x^2)/(x^2-4y^2) - (x+2)/(x^2-4y^2)]
     * where the denominator x^2-4y^2 has variables not in the
     * numerator. */
    bool inside_simplify = (factor_memo_top() != NULL);

    size_t num_vc = 0, num_vcap = 16;
    Expr** num_vars = malloc(sizeof(Expr*) * num_vcap);
    collect_variables(num, &num_vars, &num_vc, &num_vcap);

    Expr** den_vars = NULL;
    size_t den_vc = 0;
    if (!inside_simplify) {
        size_t den_vcap = 16;
        den_vars = malloc(sizeof(Expr*) * den_vcap);
        collect_variables(den, &den_vars, &den_vc, &den_vcap);
    }

    Expr* f_num = NULL;
#ifdef USE_FLINT
    /* The numerator and denominator are genuine polynomials now, so the FLINT
     * multivariate factoriser that a direct Factor[poly] uses applies to each:
     * Factor of a large rational expression in 7 variables spent ~50 s in the
     * classical FactorSquareFree / heuristic path below, where Together plus
     * the two FLINT factorisations take ~1 s (FeynCalc's FactorList2). Plain
     * single-argument form only, as for the fast path above. */
    /* Each side is tried on its own, whatever its variable count: a
     * univariate side used to take the int64 Berlekamp-Zassenhaus path
     * below, which hangs (1/(x^43 - 1)), crashes (1/(x^6001 + x + 1)) or
     * leaves the side unfactored once a coefficient outgrows a machine word. */
    Expr* flint_den = NULL;
    if (res->data.function.arg_count == 1) {
        if (num_vc >= 1) f_num = factor_flint_side(num);
        if (den_vc >= 1 || inside_simplify) flint_den = factor_flint_side(den);
        if (!f_num && num_vc == 0) f_num = expr_copy(num);
    }
    if (f_num) {
        /* fall through to the denominator below with flint_den set */
    } else
#endif
    if (num_vc == 1) {
        f_num = bz_factor_to_expr(num, num_vars[0]);
    } else if (num_vc == 2) {
        /* Fast path for squarefree bivariate inputs: try the
         * bivariate Hensel directly without first running
         * FactorSquareFree.  The Hensel orchestrator skips alphas
         * where the image is not squarefree, so for inputs with
         * repeated factors it will fail and we fall back to the
         * generic FactorSquareFree + heuristic_factor pipeline.
         * For already-squarefree inputs (the common case), this
         * skips ~100+ ms of wasted FactorSquareFree work and the
         * is_likely_irreducible probe. */
        Expr* via_hensel = factor_bivariate_via_hensel(num, num_vars);
        if (via_hensel) {
            f_num = via_hensel;
        } else {
            Expr* sq_num = eval_and_free(expr_new_function(expr_new_symbol("FactorSquareFree"), (Expr*[]){expr_copy(num)}, 1));
            f_num = heuristic_factor(sq_num);
            expr_free(sq_num);
        }
    } else {
        /* For n >= 3, run FactorSquareFree first.  It usually does
         * most of the work for inputs with separable factors (e.g.,
         * (1-x^12)(1+x-y^13)(1-y-z^14) -> Times[squarefree pieces]),
         * after which heuristic_factor's Times-recursion factors
         * each piece independently as a smaller problem. */
        Expr* sq_num = eval_and_free(expr_new_function(expr_new_symbol("FactorSquareFree"), (Expr*[]){expr_copy(num)}, 1));
        f_num = heuristic_factor(sq_num);
        expr_free(sq_num);
    }

    /* Choose den's variable list:
     *  - Simplify-context: reuse num's list (legacy shared scope).
     *  - Direct call: use the den-specific list. */
    Expr** den_eff_vars = inside_simplify ? num_vars : den_vars;
    size_t den_eff_vc   = inside_simplify ? num_vc   : den_vc;

    Expr* f_den = NULL;
#ifdef USE_FLINT
    if (flint_den) f_den = flint_den;
    else
#endif
    if (den_eff_vc == 0) {
        f_den = expr_copy(den);
    } else if (den_eff_vc == 1) {
        f_den = bz_factor_to_expr(den, den_eff_vars[0]);
    } else {
        Expr* sq_den = eval_and_free(expr_new_function(expr_new_symbol("FactorSquareFree"), (Expr*[]){expr_copy(den)}, 1));
        f_den = heuristic_factor(sq_den);
        expr_free(sq_den);
    }

    for(size_t i=0; i<num_vc; i++) expr_free(num_vars[i]);
    free(num_vars);
    if (den_vars) {
        for(size_t i=0; i<den_vc; i++) expr_free(den_vars[i]);
        free(den_vars);
    }

    /* Whether anything had to go over a common denominator; read before `den` is
     * freed, and used by the product guard below. */
    bool no_denominator = (den->type == EXPR_INTEGER && den->data.integer == 1);

    Expr* result;
    if (f_den->type == EXPR_INTEGER && f_den->data.integer == 1) {
        result = f_num;
        expr_free(f_den);
    } else {
        Expr* inv_den = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){f_den, expr_new_integer(-1)}, 2));
        result = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){f_num, inv_den}, 2));
    }

    expr_free(together);
    expr_free(num);
    expr_free(den);

    /* A PRODUCT MUST NOT COME BACK AS A SUM.
     *
     * The classical path normalises through Together, and Together expands a
     * multivariate product: Together[E^x (x^2 - 1)] is -E^x + x^2 E^x. The
     * expansion is then not a polynomial in the single variable
     * collect_variables finds (x) — its coefficients carry E^x — so
     * bz_factor_to_expr has nothing to work with and hands the expansion straight
     * back. Factor[E^x (x^2 - 1)] returned -E^x + x^2 E^x: LESS factored than the
     * input it was given, which no call to Factor should ever be.
     *
     * Factoring is multiplicative, so a product's answer is the product of its
     * factors' answers. That rewrite is only valid when nothing had to be put
     * over a common denominator: Factor[(x^2 - 1)/(x - 1)] is legitimately the
     * sum 1 + x, and threading over that Times would hand back an uncancelled
     * (x - 1)(x + 1)/(x - 1) — hence no_denominator. Under that guard this fires
     * ONLY where the classical path has already failed to factor, so nothing
     * that factors today takes this route.
     *
     * The root defect is Together's expansion, which is deliberately left alone:
     * a great deal of the integrator, DSolve and Simplify is tuned around the
     * shape it returns, so changing it is its own piece of work with its own
     * corpus runs, not a rider on this one. */
    if (no_denominator
        && arg->type == EXPR_FUNCTION && arg->data.function.head
        && arg->data.function.head->type == EXPR_SYMBOL
        && arg->data.function.head->data.symbol.name == SYM_Times
        && arg->data.function.arg_count >= 2
        && result->type == EXPR_FUNCTION && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && result->data.function.head->data.symbol.name == SYM_Plus) {
        size_t nf = arg->data.function.arg_count;
        Expr** parts = (Expr**)malloc(sizeof(Expr*) * nf);
        if (parts) {
            for (size_t i = 0; i < nf; i++) {
                /* Each factor is a strict sub-expression and cannot itself be a
                 * Times (Times is Flat), so this recursion is one level deep. */
                parts[i] = eval_and_free(expr_new_function(expr_new_symbol("Factor"),
                             (Expr*[]){expr_copy(arg->data.function.args[i])}, 1));
            }
            Expr* threaded = eval_and_free(
                expr_new_function(expr_new_symbol("Times"), parts, nf));
            free(parts);
            expr_free(result);
            result = threaded;
        }
    }

    /* Cache the (input, output) pair for the duration of the current
     * Simplify call.  memo_key holds a deep copy of the original input
     * we made before computation; we hand it off into the memo here. */
    if (memo_key) {
        factor_memo_put(memo, memo_key, result);
        expr_free(memo_key);
    }

    return result;
}
