/*
 * facpoly_options.inc -- option parsing for Factor and FactorSquareFree,
 * and factoring over GF(p) (Modulus -> p).
 *
 * Included by facpoly.c.  Every trailing argument must be an option rule
 * (Factor::nonopt otherwise) naming a known option (Factor::optx); the first
 * occurrence of an option wins, and unset options take their value from
 * Options[head], so SetOptions applies.  The parsed settings:
 *
 *   Extension -> None | Automatic | alpha | {alpha1, ...}
 *   GaussianIntegers -> True | False     (True adds I to the extension)
 *   Modulus -> 0 | p                     (p prime; -p means p; a symbolic
 *                                         value means 0)
 *   Trig -> True | False
 *
 * A nonzero modulus excludes an extension (Factor::mext) and Gaussian
 * integers (Factor::facim).
 */

typedef struct {
    const Expr* extension;  /* borrowed; NULL for None */
    bool gaussian;
    Expr* modulus;          /* owned; NULL for 0 (over Q), else a prime > 0 */
    bool trig;
} FactorOpts;

static void factor_opts_free(FactorOpts* o) {
    if (o->modulus) expr_free(o->modulus);
    o->modulus = NULL;
}

static bool fo_is_rule(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol.name == SYM_Rule
            || e->data.function.head->data.symbol.name == SYM_RuleDelayed)
        && e->data.function.arg_count == 2;
}

static bool fo_sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

/* Record one option setting.  Returns false (message emitted) on a value
 * that makes the call invalid. */
static bool fo_set(FactorOpts* o, const char* head, const Expr* key,
                   const Expr* val) {
    if (key->data.symbol.name == SYM_Extension) {
        o->extension = fo_sym_is(val, SYM_None) ? NULL : val;
        return true;
    }
    if (key->data.symbol.name == SYM_GaussianIntegers) {
        o->gaussian = fo_sym_is(val, SYM_True);
        return true;
    }
    if (key->data.symbol.name == SYM_Trig) {
        o->trig = fo_sym_is(val, SYM_True);
        return true;
    }
    /* Modulus */
    if (o->modulus) { expr_free(o->modulus); o->modulus = NULL; }
    if (expr_is_integer_like(val)) {
        mpz_t m;
        mpz_init(m);
        expr_to_mpz(val, m);
        mpz_abs(m, m);
        bool ok = mpz_sgn(m) == 0 || mpz_probab_prime_p(m, 30) > 0;
        if (ok && mpz_sgn(m) != 0)
            o->modulus = expr_bigint_normalize(expr_new_bigint_from_mpz(m));
        mpz_clear(m);
        if (!ok) {
            char* v = expr_to_string((Expr*)val);
            mth_message(head, "modp",
                        "Value of option Modulus -> %s should be a prime number or zero.",
                        v ? v : "?");
            free(v);
        }
        return ok;
    }
    if (val->type == EXPR_SYMBOL) return true;          /* symbolic: over Q */
    char* v = expr_to_string((Expr*)val);
    mth_message(head, "bmod",
                "Modulus %s must be an integer or an irreducible univariate polynomial.",
                v ? v : "?");
    free(v);
    return false;
}

/* Parse the options of `res` (arguments 1.. after the expression) for `head`,
 * which accepts the option names in `allowed` (NULL-terminated).  Returns
 * false, with the message emitted, when the call must stay unevaluated. */
static bool factor_parse_options(Expr* res, const char* head,
                                 const char* const* allowed, FactorOpts* o) {
    memset(o, 0, sizeof *o);
    size_t argc = res->data.function.arg_count;
    bool seen[4] = { false, false, false, false };
    const char* names[4] = { SYM_Extension, SYM_GaussianIntegers, SYM_Modulus, SYM_Trig };

    for (size_t i = 1; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (!fo_is_rule(a)) {
            char* bad = expr_to_string(a);
            char* call = expr_to_string(res);
            mth_message(head, "nonopt",
                        "Options expected (instead of %s) beyond position 1 in %s. "
                        "An option must be a rule or a list of rules.",
                        bad ? bad : "?", call ? call : "?");
            free(bad); free(call);
            factor_opts_free(o);
            return false;
        }
        const Expr* key = a->data.function.args[0];
        int k = -1;
        if (key->type == EXPR_SYMBOL)
            for (const char* const* p = allowed; *p; p++)
                if (key->data.symbol.name == *p)
                    for (int j = 0; j < 4; j++)
                        if (names[j] == *p) k = j;
        if (k < 0) {
            char* kn = expr_to_string((Expr*)key);
            char* call = expr_to_string(res);
            mth_message(head, "optx", "Unknown option %s in %s.",
                        kn ? kn : "?", call ? call : "?");
            free(kn); free(call);
            factor_opts_free(o);
            return false;
        }
        if (seen[k]) continue;                         /* first one wins */
        seen[k] = true;
        if (!fo_set(o, head, key, a->data.function.args[1])) {
            factor_opts_free(o);
            return false;
        }
    }

    /* Unset options default from Options[head] (SetOptions). */
    Expr* defs = symtab_get_options(head);
    if (defs && defs->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < defs->data.function.arg_count; i++) {
            const Expr* r = defs->data.function.args[i];
            if (!fo_is_rule(r) || r->data.function.args[0]->type != EXPR_SYMBOL) continue;
            for (int j = 0; j < 4; j++) {
                if (seen[j] || r->data.function.args[0]->data.symbol.name != names[j]) continue;
                seen[j] = true;
                if (!fo_set(o, head, r->data.function.args[0], r->data.function.args[1])) {
                    factor_opts_free(o);
                    return false;
                }
            }
        }
    }

    if (o->modulus && o->extension) {
        mth_message(head, "mext",
                    "Conflicting option setting. If Modulus is not 0, then "
                    "Extension should be None.");
        factor_opts_free(o);
        return false;
    }
    if (o->modulus && o->gaussian) {
        mth_message(head, "facim",
                    "Factoring non-integer polynomials with respect to a modulus "
                    "is not implemented.");
        factor_opts_free(o);
        return false;
    }
    return true;
}

static bool fo_contains_complex(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_I;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Complex) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (fo_contains_complex(e->data.function.args[i])) return true;
    return false;
}

/* {{c, 1}, {f1, e1}, ...} -> c f1^e1 ... (c dropped when it is 1). */
static Expr* fo_pairs_to_product(Expr* fl) {
    size_t np = fl->data.function.arg_count;
    Expr** terms = malloc(np * sizeof(Expr*));
    size_t nt = 0;
    for (size_t i = 0; i < np; i++) {
        Expr* pr = fl->data.function.args[i];
        Expr* base = pr->data.function.args[0];
        Expr* e = pr->data.function.args[1];
        if (i == 0 && base->type == EXPR_INTEGER && base->data.integer == 1) continue;
        terms[nt++] = (e->type == EXPR_INTEGER && e->data.integer == 1)
            ? expr_copy(base)
            : expr_new_function(expr_new_symbol(SYM_Power),
                                (Expr*[]){ expr_copy(base), expr_copy(e) }, 2);
    }
    Expr* out = nt == 0 ? expr_new_integer(1)
              : nt == 1 ? terms[0]
              : expr_new_function(expr_new_symbol(SYM_Times), terms, nt);
    free(terms);
    return eval_and_free(out);
}

/* One polynomial over GF(p).  *status as flint_modular_factor_list. */
static Expr* fo_modular_poly(const Expr* p, const Expr* m, int squarefree, int* status) {
    Expr* fl = flint_modular_factor_list(p, m, squarefree, status);
    if (!fl) return NULL;
    Expr* out = fo_pairs_to_product(fl);
    expr_free(fl);
    return out;
}

/* The inverse of the integer `c` mod m as an Integer, or NULL when c = 0. */
static Expr* fo_mod_inverse(const Expr* c, const Expr* m) {
    mpz_t a, mm;
    mpz_init(a); mpz_init(mm);
    expr_to_mpz(c, a);
    expr_to_mpz(m, mm);
    Expr* out = NULL;
    if (mpz_invert(a, a, mm)) out = expr_bigint_normalize(expr_new_bigint_from_mpz(a));
    mpz_clear(a); mpz_clear(mm);
    return out;
}

/* Reduce the integer content of a product mod m (the factors themselves are
 * already reduced and monic): 6 (1 + x) over GF(5) is 1 + x. */
static Expr* fo_reduce_content(Expr* e, const Expr* m) {
    Expr* c = NULL;
    if (expr_is_integer_like(e)) c = e;
    else if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
             && e->data.function.head->data.symbol.name == SYM_Times
             && e->data.function.arg_count >= 1
             && expr_is_integer_like(e->data.function.args[0]))
        c = e->data.function.args[0];
    if (!c) return e;
    mpz_t a, mm;
    mpz_init(a); mpz_init(mm);
    expr_to_mpz(c, a);
    expr_to_mpz(m, mm);
    mpz_mod(a, a, mm);
    Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(a));
    mpz_clear(a); mpz_clear(mm);
    if (c == e) { expr_free(e); return r; }
    size_t n = e->data.function.arg_count;
    Expr** args = malloc(n * sizeof(Expr*));
    args[0] = r;
    for (size_t i = 1; i < n; i++) args[i] = expr_copy(e->data.function.args[i]);
    Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), args, n));
    free(args);
    expr_free(e);
    return out;
}

/* Factor (squarefree != 0: square-free factor) `e` over GF(m).  Threads over
 * lists, relations and logical connectives; a rational function factors its
 * numerator and denominator separately.  Returns NULL for an unevaluated call
 * (message emitted where one applies); an expression that is not a polynomial
 * comes back unchanged. */
static Expr* factor_modular(const char* head, const Expr* e, const Expr* m, int squarefree) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_List || h == SYM_Equal || h == SYM_Unequal || h == SYM_Less
            || h == SYM_LessEqual || h == SYM_Greater || h == SYM_GreaterEqual
            || h == SYM_And || h == SYM_Or || h == SYM_Not) {
            size_t n = e->data.function.arg_count;
            Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
            for (size_t i = 0; i < n; i++) {
                a[i] = factor_modular(head, e->data.function.args[i], m, squarefree);
                if (!a[i]) {
                    for (size_t j = 0; j < i; j++) expr_free(a[j]);
                    free(a);
                    return NULL;
                }
            }
            Expr* out = eval_and_free(expr_new_function(expr_copy(e->data.function.head), a, n));
            free(a);
            return out;
        }
    }
    /* A product or a positive integer power factors base by base, so
     * (1 + x)^32768 (1 + x + x^2) is never expanded (over Q its
     * coefficients are enormous); equal factors merge in Times, and the
     * numeric content is reduced mod m at the end. */
    if (!squarefree && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        bool is_pow = h == SYM_Power && e->data.function.arg_count == 2
                   && e->data.function.args[1]->type == EXPR_INTEGER
                   && e->data.function.args[1]->data.integer > 1;
        if (h == SYM_Times || is_pow) {
            size_t n = is_pow ? 1 : e->data.function.arg_count;
            Expr** a = malloc(n * sizeof(Expr*));
            for (size_t i = 0; i < n; i++) {
                a[i] = factor_modular(head, e->data.function.args[i], m, squarefree);
                if (!a[i]) {
                    for (size_t j = 0; j < i; j++) expr_free(a[j]);
                    free(a);
                    return NULL;
                }
            }
            Expr* prod = is_pow
                ? expr_new_function(expr_new_symbol(SYM_Power),
                                    (Expr*[]){ a[0], expr_copy(e->data.function.args[1]) }, 2)
                : expr_new_function(expr_new_symbol(SYM_Times), a, n);
            free(a);
            return fo_reduce_content(eval_and_free(prod), m);
        }
    }
    if (fo_contains_complex(e)) {
        mth_message(head, "facim",
                    "Factoring non-integer polynomials with respect to a modulus "
                    "is not implemented.");
        return NULL;
    }
    if (internal_contains_inexact(e)) return NULL;

    int status = 1;
    Expr* out = NULL;
    if (!factor_has_leftover_denominator(e)) {
        out = fo_modular_poly(e, m, squarefree, &status);
    } else {
        Expr* tog = eval_and_free(expr_new_function(expr_new_symbol(SYM_Together),
                                                    (Expr*[]){ expr_copy((Expr*)e) }, 1));
        Expr* num = eval_and_free(expr_new_function(expr_new_symbol(SYM_Numerator),
                                                    (Expr*[]){ expr_copy(tog) }, 1));
        Expr* den = eval_and_free(expr_new_function(expr_new_symbol(SYM_Denominator),
                                                    (Expr*[]){ expr_copy(tog) }, 1));
        expr_free(tog);
        if (!factor_has_leftover_denominator(num) && !factor_has_leftover_denominator(den)) {
            int sd = 1;
            Expr* fn = fo_modular_poly(num, m, squarefree, &status);
            Expr* fd = fn ? fo_modular_poly(den, m, squarefree, &sd) : NULL;
            if (fn && fd) {
                if (fd->type == EXPR_INTEGER && fd->data.integer == 0) {
                    mth_message("Power", "infy", "Infinite expression 1/0 encountered.");
                    expr_free(fn);
                    status = 3;
                } else if (expr_is_integer_like(fd)) {
                    /* a constant denominator is a unit mod m */
                    Expr* inv = fo_mod_inverse(fd, m);
                    Expr* scaled = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                                     (Expr*[]){ inv, expr_copy(num) }, 2));
                    expr_free(fn);
                    out = fo_modular_poly(scaled, m, squarefree, &status);
                    expr_free(scaled);
                } else {
                    out = eval_and_free(internal_divide((Expr*[]){ fn, expr_copy(fd) }, 2));
                }
                expr_free(fd);
            } else {
                if (fn) expr_free(fn);
                if (sd == 2) status = 2;
            }
        }
        expr_free(num);
        expr_free(den);
    }
    if (out) return out;
    if (status == 2) {
        char* es = expr_to_string((Expr*)e);
        char* ms = expr_to_string((Expr*)m);
        mth_message(head, "nmod", "%s is not valid modulo %s.", es ? es : "?", ms ? ms : "?");
        free(es); free(ms);
        return NULL;
    }
    if (status == 3) return NULL;
    /* Not a polynomial in symbol generators: left as it is. */
    return expr_copy((Expr*)e);
}

bool factor_options_valid(Expr* res, const char* head, int squarefree) {
    const char* const all[] = { SYM_Extension, SYM_GaussianIntegers, SYM_Modulus,
                                SYM_Trig, NULL };
    const char* const sqf[] = { SYM_Extension, SYM_Modulus, SYM_Trig, NULL };
    FactorOpts fo;
    if (!factor_parse_options(res, head, squarefree ? sqf : all, &fo)) return false;
    factor_opts_free(&fo);
    return true;
}
