/*
 * simp_sf_engine.c -- table-driven Simplify pass for special functions.
 *
 * See simp_sf_engine.h for the model (families, affine shift matrices,
 * canonicalisers, partner Casoratians, extra passes) and
 * docs/design/simp_special_functions.md for the identities themselves.
 *
 * Pipeline for one call:
 *   0. prefilter  -- one tree walk; NULL at once when no registered head is
 *                    present (the common case: cost ~ one pointer-hash probe
 *                    per function node).
 *   1. canon      -- bottom-up replacement of atoms by their canonical
 *                    generators (descriptor `canon` hooks).
 *   2. families   -- group atoms whose shift arguments differ by integers
 *                    (head, non-shift args, and per shift argument its
 *                    symbolic part s and rational part mod 1), rewrite each
 *                    family onto the best basis, partner pairs jointly and
 *                    reduced modulo their Casoratian.
 *   3. passes     -- whole-expression candidates, kept on a strict drop.
 *   4. gate       -- the result must score strictly below the input; it is
 *                    then dispatched once more (every family is now on a
 *                    basis, so that recursion does not re-enter step 2's
 *                    rewrite).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "simp.h"
#include "simp_internal.h"
#include "simp_sf_engine.h"

/* ---- module registry ---------------------------------------------------- */

static const SfModule* const SF_MODULES[] = {
    &SF_MODULE_BESSEL,
    &SF_MODULE_GAMMA,
    &SF_MODULE_ZETA,
    &SF_MODULE_ERREXP,
    &SF_MODULE_ORTHOPOLY,
    &SF_MODULE_AIRY,
    &SF_MODULE_HYPER,
    &SF_MODULE_MISC,
};
#define SF_NMODULES (sizeof SF_MODULES / sizeof SF_MODULES[0])

/* Pointer hash keyed by (interned head, arity) -> descriptor, plus a
 * head-only presence set (any arity, and pass-only heads) for the prefilter.
 * Built on first use; interned names are stable for the process lifetime. */
#define SF_HASH_CAP 512u
static const char*   g_desc_key[SF_HASH_CAP];
static size_t        g_desc_arity[SF_HASH_CAP];
static const SfDesc* g_desc_val[SF_HASH_CAP];
static const char*   g_head_key[SF_HASH_CAP];
static bool          g_hash_ready = false;

static size_t ptr_slot(const char* p, size_t salt) {
    uintptr_t h = (uintptr_t)p + (uintptr_t)salt * 0x632BE59BD9B4E019ull;
    h ^= h >> 17; h *= (uintptr_t)0x9E3779B97F4A7C15ull; h ^= h >> 29;
    return (size_t)(h & (SF_HASH_CAP - 1));
}

static void head_put(const char* key) {
    if (!key) return;
    size_t i = ptr_slot(key, 0);
    while (g_head_key[i] && g_head_key[i] != key) i = (i + 1) & (SF_HASH_CAP - 1);
    g_head_key[i] = key;
}

static bool head_present(const char* key) {
    size_t i = ptr_slot(key, 0);
    while (g_head_key[i]) {
        if (g_head_key[i] == key) return true;
        i = (i + 1) & (SF_HASH_CAP - 1);
    }
    return false;
}

static void desc_put(const SfDesc* d) {
    const char* key = *d->head;
    if (!key) return;
    size_t i = ptr_slot(key, d->arity);
    while (g_desc_key[i] && !(g_desc_key[i] == key && g_desc_arity[i] == d->arity))
        i = (i + 1) & (SF_HASH_CAP - 1);
    if (!g_desc_key[i]) { g_desc_key[i] = key; g_desc_arity[i] = d->arity; g_desc_val[i] = d; }
}

static const SfDesc* desc_get(const char* key, size_t arity) {
    size_t i = ptr_slot(key, arity);
    while (g_desc_key[i]) {
        if (g_desc_key[i] == key && g_desc_arity[i] == arity) return g_desc_val[i];
        i = (i + 1) & (SF_HASH_CAP - 1);
    }
    return NULL;
}

static void registry_init(void) {
    if (g_hash_ready) return;
    for (size_t m = 0; m < SF_NMODULES; m++)
        if (SF_MODULES[m]->init) SF_MODULES[m]->init();
    for (size_t m = 0; m < SF_NMODULES; m++) {
        const SfModule* mod = SF_MODULES[m];
        for (size_t i = 0; i < mod->npass_heads; i++) head_put(*mod->pass_heads[i]);
        for (size_t i = 0; i < mod->ndesc; i++) {
            head_put(*mod->descs[i].head);
            desc_put(&mod->descs[i]);
        }
    }
    g_hash_ready = true;
}

/* Look up the descriptor for a call atom (matching arity), else NULL. */
static const SfDesc* desc_of(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return NULL;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return NULL;
    return desc_get(h->data.symbol.name, e->data.function.arg_count);
}

static bool any_registered_head(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && head_present(h->data.symbol.name)) return true;
    if (h->type == EXPR_FUNCTION && any_registered_head(h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (any_registered_head(e->data.function.args[i])) return true;
    return false;
}

/* ---- small builders ------------------------------------------------------ */

Expr* sf_cp(const Expr* e) { return expr_copy((Expr*)e); }
Expr* sf_int(int64_t n) { return expr_new_integer(n); }
Expr* sf_sym(const char* name) { return expr_new_symbol(name); }

Expr* sf_call(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
Expr* sf_mk1(const char* head, Expr* a) {
    return eval_and_free(expr_new_function(expr_new_symbol(head), (Expr*[]){ a }, 1));
}
Expr* sf_mk2(const char* head, Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(expr_new_symbol(head), (Expr*[]){ a, b }, 2));
}
Expr* sf_times(Expr* a, Expr* b) { return sf_mk2(SYM_Times, a, b); }
Expr* sf_plus(Expr* a, Expr* b)  { return sf_mk2(SYM_Plus, a, b); }
Expr* sf_pow(Expr* a, Expr* b)   { return sf_mk2(SYM_Power, a, b); }
Expr* sf_minus(Expr* a, Expr* b) { return sf_plus(a, sf_times(sf_int(-1), b)); }
Expr* sf_div(Expr* a, Expr* b)   { return sf_times(a, sf_pow(b, sf_int(-1))); }
Expr* sf_rat(int64_t n, int64_t d) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Rational),
                                           (Expr*[]){ sf_int(n), sf_int(d) }, 2));
}

bool sf_is_call(const Expr* e, const char* head, size_t n) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == head &&
           e->data.function.arg_count == n;
}

bool sf_contains_head(const Expr* e, const char* head) {
    if (e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == head) return true;
    if (h->type == EXPR_FUNCTION && sf_contains_head(h, head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (sf_contains_head(e->data.function.args[i], head)) return true;
    return false;
}

static bool is_zero(const Expr* e) { return e->type == EXPR_INTEGER && e->data.integer == 0; }
static bool is_one(const Expr* e)  { return e->type == EXPR_INTEGER && e->data.integer == 1; }

static Expr* together(Expr* e) {
    Expr* t = call_unary_copy("Together", e);
    if (!t) return e;
    expr_free(e);
    return t;
}

/* ---- order splitting ----------------------------------------------------- */

static bool read_rational(const Expr* e, int64_t* num, int64_t* den) {
    if (e->type == EXPR_INTEGER) { *num = e->data.integer; *den = 1; return true; }
    if (sf_is_call(e, SYM_Rational, 2) &&
        e->data.function.args[0]->type == EXPR_INTEGER &&
        e->data.function.args[1]->type == EXPR_INTEGER) {
        *num = e->data.function.args[0]->data.integer;
        *den = e->data.function.args[1]->data.integer;
        return *den > 0;
    }
    return false;
}

bool sf_split(const Expr* mu, Expr** s, int64_t* num, int64_t* den) {
    if (read_rational(mu, num, den)) { *s = NULL; return true; }
    if (mu->type == EXPR_BIGINT) return false;
    if (sf_is_call(mu, SYM_Rational, 2)) return false;      /* bignum rational */
    if (mu->type == EXPR_FUNCTION && mu->data.function.head->type == EXPR_SYMBOL &&
        mu->data.function.head->data.symbol.name == SYM_Plus &&
        mu->data.function.arg_count >= 2 &&
        read_rational(mu->data.function.args[0], num, den)) {
        size_t n = mu->data.function.arg_count;
        if (n == 2) {
            *s = sf_cp(mu->data.function.args[1]);
        } else {
            Expr** rest = malloc(sizeof(Expr*) * (n - 1));
            for (size_t i = 1; i < n; i++) rest[i - 1] = sf_cp(mu->data.function.args[i]);
            *s = eval_and_free(sf_call(SYM_Plus, rest, n - 1));
            free(rest);
        }
        return true;
    }
    *s = sf_cp(mu); *num = 0; *den = 1;
    return true;
}

static bool num_negative(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer < 0;
    if (e->type == EXPR_REAL)    return e->data.real < 0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) < 0;
    if (sf_is_call(e, SYM_Rational, 2)) return num_negative(e->data.function.args[0]);
    /* Complex[a, b] reads negative when a < 0, or a = 0 and b < 0 (-I x) */
    if (sf_is_call(e, SYM_Complex, 2)) {
        const Expr* re = e->data.function.args[0];
        bool re_zero = re->type == EXPR_INTEGER && re->data.integer == 0;
        return num_negative(re) || (re_zero && num_negative(e->data.function.args[1]));
    }
    return false;
}

bool sf_reads_negative(const Expr* s) {
    if (!s) return false;
    if (num_negative(s)) return true;
    if (s->type == EXPR_FUNCTION && s->data.function.head->type == EXPR_SYMBOL &&
        s->data.function.arg_count >= 1) {
        const char* h = s->data.function.head->data.symbol.name;
        if (h == SYM_Times) return num_negative(s->data.function.args[0]);
        if (h == SYM_Plus)  return sf_reads_negative(s->data.function.args[0]);
    }
    return false;
}

bool sf_proved_integer(const AssumeCtx* ctx, const Expr* mu) {
    if (mu->type == EXPR_INTEGER || mu->type == EXPR_BIGINT) return true;
    if (!ctx) return false;
    if (prov_int(ctx, mu)) return true;
    Expr* s = NULL; int64_t num, den;
    bool ok = sf_split(mu, &s, &num, &den) && den == 1 && (!s || prov_int(ctx, s));
    if (s) expr_free(s);
    return ok;
}

bool sf_reflection_flips(const Expr* mu, int64_t k) {
    Expr* s = NULL; int64_t num, den;
    if (!sf_split(mu, &s, &num, &den)) return false;
    if (s) { bool neg = sf_reads_negative(s); expr_free(s); return neg; }
    int64_t r = num % den; if (r < 0) r += den;      /* frac(mu) = r/den */
    if (2 * r > den) return true;
    if (2 * r == den || r == 0) return 2 * num < -k * den;
    return false;
}

/* ---- reduction modulo a relation ----------------------------------------- */

/* Replace every occurrence of from[i] in e by to[i] (structural; owned out). */
static Expr* subst_exprs(const Expr* e, Expr* const* from, Expr* const* to, size_t n) {
    for (size_t i = 0; i < n; i++)
        if (expr_eq(e, from[i])) return sf_cp(to[i]);
    if (e->type != EXPR_FUNCTION) return sf_cp(e);
    size_t na = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (na ? na : 1));
    for (size_t i = 0; i < na; i++) a[i] = subst_exprs(e->data.function.args[i], from, to, n);
    Expr* r = expr_new_function(subst_exprs(e->data.function.head, from, to, n), a, na);
    free(a);
    return r;
}

Expr* sf_reduce_modulo(const Expr* e, const Expr* R, Expr* const* vars, size_t nv) {
    /* The generator kernels become fresh symbols for the duration:
     * Together maps into a kernel's arguments (AiryAi[1/3 + x] ->
     * AiryAi[(1 + 3 x)/3]), and the generators must still match `vars`
     * when PolynomialReduce sees them.  The symbols' names keep the
     * variables' order irrelevant (PolynomialReduce takes it explicitly). */
    Expr** sym = malloc(sizeof(Expr*) * (nv ? nv : 1));
    for (size_t i = 0; i < nv; i++) {
        char nm[32];
        snprintf(nm, sizeof nm, "SfEngine$g%zu", i);
        sym[i] = sf_sym(nm);
    }
    Expr* es = eval_and_free(subst_exprs(e, vars, sym, nv));
    Expr* Rs = eval_and_free(subst_exprs(R, vars, sym, nv));
    Expr* t = together(es);
    Expr* num = sf_mk1("Numerator", sf_cp(t));
    Expr* den = sf_mk1("Denominator", t);
    Expr** va = malloc(sizeof(Expr*) * (nv ? nv : 1));
    for (size_t i = 0; i < nv; i++) va[i] = sf_cp(sym[i]);
    Expr* vlist = sf_call(SYM_List, va, nv);
    free(va);
    Expr* rlist = sf_call(SYM_List, (Expr*[]){ Rs }, 1);
    Expr* pr = eval_and_free(sf_call("PolynomialReduce", (Expr*[]){ num, rlist, vlist }, 3));
    Expr* out = NULL;
    if (sf_is_call(pr, SYM_List, 2)) {
        Expr* q = sf_div(sf_cp(pr->data.function.args[1]), den);
        den = NULL;
        out = eval_and_free(subst_exprs(q, sym, vars, nv));
        expr_free(q);
    }
    expr_free(pr);
    if (den) expr_free(den);
    for (size_t i = 0; i < nv; i++) expr_free(sym[i]);
    free(sym);
    return out;
}

/* ---- families ------------------------------------------------------------ */

typedef struct {
    const SfDesc* d;
    const Expr*   proto;                  /* borrowed: first member atom      */
    Expr*         s[SF_MAX_SHIFT];        /* owned symbolic parts or NULL     */
    int64_t       den[SF_MAX_SHIFT], cls[SF_MAX_SHIFT], base_num[SF_MAX_SHIFT];
    int64_t*      off;                    /* n x nshift offsets               */
    const Expr**  atom;                   /* borrowed member atoms            */
    size_t        n, cap;
    bool          done;
} Fam;

typedef struct { Fam* f; size_t n, cap; } FamSet;

static int64_t mod_pos(int64_t a, int64_t m) { int64_t r = a % m; return r < 0 ? r + m : r; }

static bool s_eq(const Expr* a, const Expr* b) {
    return (!a && !b) || (a && b && expr_eq((Expr*)a, (Expr*)b));
}

static bool is_shift_arg(const SfDesc* d, size_t i) {
    for (size_t j = 0; j < d->nshift; j++) if (d->shift_arg[j] == i) return true;
    return false;
}

static void fam_add(Fam* f, const int64_t* num, const Expr* atom) {
    size_t ns = f->d->nshift;
    int64_t k[SF_MAX_SHIFT];
    for (size_t j = 0; j < ns; j++) k[j] = (num[j] - f->base_num[j]) / f->den[j];
    for (size_t i = 0; i < f->n; i++) {
        bool same = true;
        for (size_t j = 0; j < ns && same; j++) same = f->off[i * ns + j] == k[j];
        if (same) return;
    }
    if (f->n == f->cap) {
        f->cap = f->cap ? f->cap * 2 : 4;
        f->off  = realloc(f->off, sizeof(int64_t) * f->cap * (ns ? ns : 1));
        f->atom = realloc((void*)f->atom, sizeof(Expr*) * f->cap);
    }
    for (size_t j = 0; j < ns; j++) f->off[f->n * ns + j] = k[j];
    f->atom[f->n++] = atom;
}

static void collect(const Expr* e, FamSet* fs) {
    if (e->type != EXPR_FUNCTION) return;
    const SfDesc* d = desc_of(e);
    if (d && d->nshift > 0 && d->up) {
        Expr* s[SF_MAX_SHIFT] = { NULL };
        int64_t num[SF_MAX_SHIFT], den[SF_MAX_SHIFT];
        bool ok = true;
        size_t got = 0;
        for (size_t j = 0; j < d->nshift && ok; j++, got++)
            ok = sf_split(e->data.function.args[d->shift_arg[j]], &s[j], &num[j], &den[j]);
        if (ok) {
            Fam* hit = NULL;
            for (size_t i = 0; i < fs->n && !hit; i++) {
                Fam* f = &fs->f[i];
                if (f->d != d) continue;
                bool same = true;
                for (size_t a = 0; a < d->arity && same; a++)
                    if (!is_shift_arg(d, a))
                        same = expr_eq(f->proto->data.function.args[a], e->data.function.args[a]);
                for (size_t j = 0; j < d->nshift && same; j++)
                    same = f->den[j] == den[j] && f->cls[j] == mod_pos(num[j], den[j]) &&
                           s_eq(f->s[j], s[j]);
                if (same) hit = f;
            }
            if (hit) {
                fam_add(hit, num, e);
                for (size_t j = 0; j < d->nshift; j++) if (s[j]) expr_free(s[j]);
            } else {
                if (fs->n == fs->cap) {
                    fs->cap = fs->cap ? fs->cap * 2 : 4;
                    fs->f = realloc(fs->f, sizeof(Fam) * fs->cap);
                }
                Fam* f = &fs->f[fs->n++];
                memset(f, 0, sizeof *f);
                f->d = d; f->proto = e;
                for (size_t j = 0; j < d->nshift; j++) {
                    f->s[j] = s[j]; f->den[j] = den[j];
                    f->cls[j] = mod_pos(num[j], den[j]); f->base_num[j] = num[j];
                }
                fam_add(f, num, e);
            }
        } else {
            for (size_t j = 0; j < got; j++) if (s[j]) expr_free(s[j]);
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect(e->data.function.args[i], fs);
}

static void famset_free(FamSet* fs) {
    for (size_t i = 0; i < fs->n; i++) {
        for (size_t j = 0; j < SF_MAX_SHIFT; j++) if (fs->f[i].s[j]) expr_free(fs->f[i].s[j]);
        free(fs->f[i].off); free((void*)fs->f[i].atom);
    }
    free(fs->f);
}

/* The family's shift argument j at lattice offset k: s + (base + k den)/den. */
static Expr* coord_at(const Fam* f, size_t j, int64_t k) {
    Expr* c = sf_rat(f->base_num[j] + k * f->den[j], f->den[j]);
    return f->s[j] ? sf_plus(sf_cp(f->s[j]), c) : c;
}

/* Arguments of the member at offset vector k (owned array of owned Exprs). */
static Expr** args_at(const Fam* f, const int64_t* k) {
    size_t n = f->d->arity;
    Expr** a = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) a[i] = NULL;
    for (size_t j = 0; j < f->d->nshift; j++) a[f->d->shift_arg[j]] = coord_at(f, j, k[j]);
    for (size_t i = 0; i < n; i++)
        if (!a[i]) a[i] = sf_cp(f->proto->data.function.args[i]);
    return a;
}

static void args_free(Expr** a, size_t n) {
    for (size_t i = 0; i < n; i++) expr_free(a[i]);
    free(a);
}

static Expr* atom_at(const Fam* f, const int64_t* k) {
    Expr** a = args_at(f, k);
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(*f->d->head), a, f->d->arity));
    free(a);
    return r;
}

/* ---- affine matrices ------------------------------------------------------ */

/* An n x n matrix of owned Exprs (n = rank + 1), row-major. */
typedef struct { Expr* m[(SF_MAX_RANK + 1) * (SF_MAX_RANK + 1)]; } Mat;

static void mat_free(Mat* A, size_t n) {
    for (size_t i = 0; i < n * n; i++) if (A->m[i]) { expr_free(A->m[i]); A->m[i] = NULL; }
}

static void mat_identity(Mat* A, size_t n) {
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++) A->m[i * n + j] = sf_int(i == j);
}

/* Matrix entries are Together'd only for rank-2 families, whose Lommel-type
 * coefficients are rational functions that grow without it.  A rank-1
 * entry is a short sum of the descriptor's terms, and Together would only
 * rewrite them (((a + 1/2)^2)^(-3/2) -> 8 (1 + 4 a + 4 a^2)^(-3/2)) so they
 * no longer meet the input's copy. */
static Expr* entry_norm(Expr* e, size_t n) { return n == 3 ? together(e) : e; }

/* C = A B. */
static void mat_mul(Mat* C, const Mat* A, const Mat* B, size_t n) {
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++) {
            Expr* acc = NULL;
            for (size_t k = 0; k < n; k++) {
                const Expr* a = A->m[i * n + k];
                const Expr* b = B->m[k * n + j];
                if (is_zero(a) || is_zero(b)) continue;
                Expr* t = is_one(a) ? sf_cp(b) : is_one(b) ? sf_cp(a)
                                    : sf_times(sf_cp(a), sf_cp(b));
                acc = acc ? sf_plus(acc, t) : t;
            }
            C->m[i * n + j] = acc ? entry_norm(acc, n) : sf_int(0);
        }
}

/* The up matrix at point args `a`, direction d; false on decline. */
static bool mat_up(const SfDesc* d, Expr* const* a, int dir, Mat* M) {
    size_t n = (size_t)d->rank + 1;
    for (size_t i = 0; i < n * n; i++) M->m[i] = NULL;
    if (!d->up(a, dir, M->m)) { mat_free(M, n); return false; }
    for (size_t j = 0; j < n; j++) {
        if (M->m[(n - 1) * n + j]) expr_free(M->m[(n - 1) * n + j]);
        M->m[(n - 1) * n + j] = sf_int(j == n - 1);
    }
    for (size_t i = 0; i < (n - 1) * n; i++) {
        if (!M->m[i]) { mat_free(M, n); return false; }
        M->m[i] = entry_norm(eval_and_free(M->m[i]), n);
    }
    return true;
}

/* Inverse of an affine matrix (last row 0..0 1).  False when singular. */
static bool mat_inv(Mat* R, const Mat* M, size_t n) {
    if (n == 2) {
        const Expr* a = M->m[0];
        const Expr* g = M->m[1];
        if (is_zero(a)) return false;
        R->m[0] = sf_pow(sf_cp(a), sf_int(-1));
        R->m[1] = sf_times(sf_int(-1), sf_div(sf_cp(g), sf_cp(a)));
        R->m[2] = sf_int(0); R->m[3] = sf_int(1);
        return true;
    }
    /* n == 3: A = [[a, b], [c, d]], t = (t0, t1) */
    const Expr *a = M->m[0], *b = M->m[1], *t0 = M->m[2];
    const Expr *c = M->m[3], *d = M->m[4], *t1 = M->m[5];
    Expr* det = together(sf_minus(sf_times(sf_cp(a), sf_cp(d)), sf_times(sf_cp(b), sf_cp(c))));
    if (is_zero(det)) { expr_free(det); return false; }
    Expr* idet = sf_pow(det, sf_int(-1));
    Expr* ia = together(sf_times(sf_cp(d), sf_cp(idet)));
    Expr* ib = together(sf_times(sf_int(-1), sf_times(sf_cp(b), sf_cp(idet))));
    Expr* ic = together(sf_times(sf_int(-1), sf_times(sf_cp(c), sf_cp(idet))));
    Expr* id = together(sf_times(sf_cp(a), sf_cp(idet)));
    expr_free(idet);
    /* t' = -A^-1 t */
    Expr* u0 = together(sf_times(sf_int(-1), sf_plus(sf_times(sf_cp(ia), sf_cp(t0)),
                                                       sf_times(sf_cp(ib), sf_cp(t1)))));
    Expr* u1 = together(sf_times(sf_int(-1), sf_plus(sf_times(sf_cp(ic), sf_cp(t0)),
                                                       sf_times(sf_cp(id), sf_cp(t1)))));
    R->m[0] = ia; R->m[1] = ib; R->m[2] = u0;
    R->m[3] = ic; R->m[4] = id; R->m[5] = u1;
    R->m[6] = sf_int(0); R->m[7] = sf_int(0); R->m[8] = sf_int(1);
    return true;
}

/* ---- lattice walk ---------------------------------------------------------- */

/* Cache of transfer matrices T(k) with v(b + k) = T(k) v(b), for one basis b. */
typedef struct {
    int64_t* key;      /* ncache x nshift */
    Mat*     T;
    size_t   n, cap;
} TCache;

static void tcache_free(TCache* c, size_t n) {
    for (size_t i = 0; i < c->n; i++) mat_free(&c->T[i], n);
    free(c->key); free(c->T);
}

static const Mat* tcache_get(const TCache* c, const int64_t* k, size_t ns) {
    for (size_t i = 0; i < c->n; i++) {
        bool same = true;
        for (size_t j = 0; j < ns && same; j++) same = c->key[i * ns + j] == k[j];
        if (same) return &c->T[i];
    }
    return NULL;
}

/* Takes ownership of *T's entries. */
static void tcache_put(TCache* c, const int64_t* k, size_t ns, Mat* T) {
    if (c->n == c->cap) {
        c->cap = c->cap ? c->cap * 2 : 8;
        c->key = realloc(c->key, sizeof(int64_t) * c->cap * ns);
        c->T   = realloc(c->T, sizeof(Mat) * c->cap);
    }
    for (size_t j = 0; j < ns; j++) c->key[c->n * ns + j] = k[j];
    c->T[c->n++] = *T;
}

/* T(target) relative to basis point b (absolute offsets), walking coordinate
 * 0 first, then 1, ...  Returns a pointer into the cache, or NULL on
 * decline. */
static const Mat* walk(const Fam* f, const int64_t* b, const int64_t* target, TCache* c) {
    const SfDesc* d = f->d;
    size_t ns = d->nshift, n = (size_t)d->rank + 1;
    int64_t rel[SF_MAX_SHIFT] = { 0 };
    if (!tcache_get(c, rel, ns)) { Mat I; mat_identity(&I, n); tcache_put(c, rel, ns, &I); }
    for (size_t j = 0; j < ns; j++) {
        int64_t goal = target[j] - b[j];
        while (rel[j] != goal) {
            int step = goal > rel[j] ? 1 : -1;
            int64_t nxt[SF_MAX_SHIFT];
            memcpy(nxt, rel, sizeof nxt);
            nxt[j] += step;
            if (!tcache_get(c, nxt, ns)) {
                const Mat* Tcur = tcache_get(c, rel, ns);
                Mat S, Tn;
                /* rank-1 homogeneous fast path: one bulk factor for the rest
                 * of this coordinate's run */
                if (d->bulk && d->rank == 1) {
                    int64_t abspt[SF_MAX_SHIFT];
                    for (size_t q = 0; q < ns; q++) abspt[q] = b[q] + rel[q];
                    Expr** a = args_at(f, abspt);
                    Expr* fac = d->bulk(a, (int)j, goal - rel[j]);
                    args_free(a, d->arity);
                    if (fac) {
                        nxt[j] = goal;
                        if (!tcache_get(c, nxt, ns)) {
                            S.m[0] = eval_and_free(fac); S.m[1] = sf_int(0);
                            S.m[2] = sf_int(0); S.m[3] = sf_int(1);
                            Tcur = tcache_get(c, rel, ns);
                            mat_mul(&Tn, &S, Tcur, n);
                            mat_free(&S, n);
                            tcache_put(c, nxt, ns, &Tn);
                        } else {
                            expr_free(fac);
                        }
                        memcpy(rel, nxt, sizeof rel);
                        continue;
                    }
                }
                int64_t at[SF_MAX_SHIFT];
                for (size_t q = 0; q < ns; q++) at[q] = b[q] + (step > 0 ? rel[q] : nxt[q]);
                Expr** a = args_at(f, at);
                Mat M;
                bool ok = mat_up(d, a, (int)j, &M);
                args_free(a, d->arity);
                if (!ok) return NULL;
                if (step > 0) {
                    S = M;
                } else {
                    ok = mat_inv(&S, &M, n);
                    mat_free(&M, n);
                    if (!ok) return NULL;
                }
                Tcur = tcache_get(c, rel, ns);
                mat_mul(&Tn, &S, Tcur, n);
                mat_free(&S, n);
                tcache_put(c, nxt, ns, &Tn);
            }
            memcpy(rel, nxt, sizeof rel);
        }
    }
    return tcache_get(c, rel, ns);
}

/* ---- substitution ----------------------------------------------------------- */

typedef struct { const Expr** from; Expr** to; size_t n; } SubstMap;

static Expr* subst_apply(const Expr* e, const SubstMap* m) {
    for (size_t i = 0; i < m->n; i++)
        if (expr_eq((Expr*)e, (Expr*)m->from[i])) return expr_copy(m->to[i]);
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) args[i] = subst_apply(e->data.function.args[i], m);
    Expr* r = expr_new_function(subst_apply(e->data.function.head, m), args, n);
    free(args);
    return r;
}

/* Rewrite every member of f inside e onto the basis at absolute offset b.
 * Evaluated result, or NULL on decline. */
static Expr* rewrite_family(const Expr* e, const Fam* f, const int64_t* b) {
    const SfDesc* d = f->d;
    size_t ns = d->nshift, n = (size_t)d->rank + 1;
    TCache c = { NULL, NULL, 0, 0 };
    Expr* G[SF_MAX_RANK] = { NULL };
    G[0] = atom_at(f, b);
    if (d->rank == 2) {
        int64_t b1[SF_MAX_SHIFT];
        memcpy(b1, b, sizeof b1);
        b1[0] += 1;
        G[1] = atom_at(f, b1);
    }
    SubstMap m = { f->atom, malloc(sizeof(Expr*) * f->n), 0 };
    bool ok = true;
    for (size_t i = 0; i < f->n && ok; i++) {
        const Mat* T = walk(f, b, &f->off[i * ns], &c);
        if (!T) { ok = false; break; }
        Expr* acc = sf_cp(T->m[n - 1]);              /* inhomogeneous part */
        for (size_t j = 0; j + 1 < n; j++)
            if (!is_zero(T->m[j]))
                acc = sf_plus(acc, sf_times(sf_cp(T->m[j]), sf_cp(G[j])));
        m.to[m.n++] = acc;
    }
    Expr* out = ok ? eval_and_free(subst_apply(e, &m)) : NULL;
    for (size_t i = 0; i < m.n; i++) expr_free(m.to[i]);
    free(m.to);
    for (size_t j = 0; j < SF_MAX_RANK; j++) if (G[j]) expr_free(G[j]);
    tcache_free(&c, n);
    return out;
}

/* Expand the shift arguments of every registered atom in e (owned in,
 * owned out).  Together maps into function arguments, so a coefficient
 * kernel such as BarnesG's Gamma[x + 1/3] comes back as Gamma[1/3 (1 + 3 x)]
 * and would no longer group with its family; Expand restores the s + c
 * spelling sf_split reads.  Non-shift arguments keep the user's form. */
static Expr* normalize_shift_args(Expr* e) {
    if (e->type != EXPR_FUNCTION) return e;
    size_t n = e->data.function.arg_count;
    bool sub = false;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        a[i] = normalize_shift_args(sf_cp(e->data.function.args[i]));
        if (a[i] != e->data.function.args[i]) sub = true;
    }
    const SfDesc* d = desc_of(e);
    if (d) {
        for (size_t j = 0; j < d->nshift; j++) {
            size_t i = d->shift_arg[j];
            if (a[i]->type != EXPR_FUNCTION) continue;
            Expr* x = sf_mk1("Expand", sf_cp(a[i]));
            if (!expr_eq(x, a[i])) { expr_free(a[i]); a[i] = x; sub = true; }
            else expr_free(x);
        }
    }
    if (!sub) {
        for (size_t i = 0; i < n; i++) expr_free(a[i]);
        free(a);
        return e;
    }
    Expr* r = eval_and_free(expr_new_function(sf_cp(e->data.function.head), a, n));
    free(a);
    expr_free(e);
    return r;
}

/* Power normal form.  b^(q + rest) with b a positive rational and q a
 * rational of magnitude >= 1 -> b^t b^(q - t + rest), t = q truncated toward
 * zero, exponent expanded: the spelling
 * the evaluator and Simplify's own passes leave on the input
 * (2^(1 + 2 x) -> 2 2^(2 x), 2^(-5/2 + 8 x) -> 2^(-1/2 + 8 x)/4), so a
 * prefactor the engine builds meets the input's copy of the same power. */
static Expr* normalize_numeric_powers(Expr* e) {
    if (e->type != EXPR_FUNCTION) return e;
    size_t n = e->data.function.arg_count;
    bool sub = false;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        a[i] = normalize_numeric_powers(sf_cp(e->data.function.args[i]));
        if (a[i] != e->data.function.args[i]) sub = true;
    }
    Expr* r;
    if (sub) {
        r = eval_and_free(expr_new_function(sf_cp(e->data.function.head), a, n));
        expr_free(e);
    } else {
        for (size_t i = 0; i < n; i++) expr_free(a[i]);
        r = e;
    }
    free(a);
    if (!sf_is_call(r, SYM_Power, 2)) return r;
    const Expr* b = r->data.function.args[0];
    int64_t bn, bd, qn, qd;
    /* (c w)^e -> c^e w^e for a positive rational c: arg(c w) = arg(w), so it
     * holds on the principal branch ((z/2)^n meets D's 2^-n z^n) */
    if (sf_is_call(b, SYM_Times, b->type == EXPR_FUNCTION ? b->data.function.arg_count : 0) &&
        read_rational(b->data.function.args[0], &bn, &bd) && bn > 0 &&
        r->data.function.args[1]->type != EXPR_INTEGER) {
        size_t nb = b->data.function.arg_count;
        Expr** rest = malloc(sizeof(Expr*) * nb);
        for (size_t i = 1; i < nb; i++) rest[i - 1] = sf_cp(b->data.function.args[i]);
        Expr* w = nb == 2 ? rest[0] : sf_call(SYM_Times, rest, nb - 1);
        free(rest);
        const Expr* ex = r->data.function.args[1];
        Expr* out = sf_times(sf_pow(sf_cp(b->data.function.args[0]), sf_cp(ex)), sf_pow(w, sf_cp(ex)));
        expr_free(r);
        return normalize_numeric_powers(out);
    }
    /* (p/q)^e -> p^e q^-e for a positive rational base and symbolic e
     * ((1/2)^(a + b) meets D's 2^(-a - b)) */
    if (read_rational(b, &bn, &bd) && bn > 0 && bd > 1 && r->data.function.args[1]->type != EXPR_INTEGER &&
        !sf_is_call(r->data.function.args[1], SYM_Rational, 2)) {
        const Expr* ex = r->data.function.args[1];
        Expr* out = sf_times(sf_pow(sf_int(bn), sf_cp(ex)), sf_pow(sf_int(bd), sf_times(sf_int(-1), sf_cp(ex))));
        expr_free(r);
        return normalize_numeric_powers(out);
    }
    if (!read_rational(b, &bn, &bd) || bn <= 0) {
        /* symbolic base: undo Together's c (p + q) spelling of an exponent
         * (w^(7/2 + a) -> w^(1/2 (7 + 2 a))), so equal powers stay equal */
        const Expr* ex0 = r->data.function.args[1];
        if (sf_is_call(ex0, SYM_Times, 2) && read_rational(ex0->data.function.args[0], &qn, &qd) &&
            ex0->data.function.args[1]->type == EXPR_FUNCTION &&
            ex0->data.function.args[1]->data.function.head->type == EXPR_SYMBOL &&
            ex0->data.function.args[1]->data.function.head->data.symbol.name == SYM_Plus) {
            Expr* xe = sf_mk1("Expand", sf_cp(ex0));
            Expr* nr = eval_and_free(sf_call(SYM_Power, (Expr*[]){ sf_cp(b), xe }, 2));
            expr_free(r);
            return nr;
        }
        return r;
    }
    /* one spelling for the exponent: 2^(2 (a + b)) and 2^(2 a + 2 b) */
    if (r->data.function.args[1]->type == EXPR_FUNCTION) {
        Expr* xe = sf_mk1("Expand", sf_cp(r->data.function.args[1]));
        if (!expr_eq(xe, r->data.function.args[1])) {
            Expr* nr = eval_and_free(sf_call(SYM_Power, (Expr*[]){ sf_cp(b), xe }, 2));
            expr_free(r);
            r = nr;
            if (!sf_is_call(r, SYM_Power, 2)) return r;
            b = r->data.function.args[0];
        } else {
            expr_free(xe);
        }
    }
    const Expr* ex = r->data.function.args[1];
    if (!(ex->type == EXPR_FUNCTION && ex->data.function.head->type == EXPR_SYMBOL &&
          ex->data.function.head->data.symbol.name == SYM_Plus &&
          read_rational(ex->data.function.args[0], &qn, &qd)))
        return r;
    int64_t t = qn / qd;                                   /* C99: toward zero */
    if (t == 0) return r;
    Expr* rest = sf_plus(sf_minus(sf_cp(ex), sf_cp(ex->data.function.args[0])), sf_rat(qn - t * qd, qd));
    Expr* out = sf_call(SYM_Times, (Expr*[]){ sf_pow(sf_cp(b), sf_int(t)),
                                              sf_call(SYM_Power, (Expr*[]){ sf_cp(b), rest }, 2) }, 2);
    expr_free(r);
    return eval_and_free(out);
}

static bool sf_contains_symbol_pi(const Expr* e) {
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_Pi;
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (sf_contains_symbol_pi(e->data.function.args[i])) return true;
    return false;
}

static bool has_trig_head(const char* h) {
    return h == SYM_Sin || h == SYM_Cos || h == SYM_Tan || h == SYM_Cot || h == SYM_Sec || h == SYM_Csc;
}

/* Sin/Cos of Pi (k/2 + t) -> +-Sin/Cos[Pi t] (quarter turns; the evaluator
 * leaves Cos[Pi/2 + Pi s] alone, so a descriptor's -Sin[Pi s] would not meet
 * the input's spelling).  Owned in/out. */
static Expr* fold_quarter_turns(Expr* e) {
    if (e->type != EXPR_FUNCTION) return e;
    size_t n = e->data.function.arg_count;
    bool sub = false;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        a[i] = fold_quarter_turns(sf_cp(e->data.function.args[i]));
        if (a[i] != e->data.function.args[i]) sub = true;
    }
    Expr* r;
    if (sub) {
        r = eval_and_free(expr_new_function(sf_cp(e->data.function.head), a, n));
        expr_free(e);
    } else {
        for (size_t i = 0; i < n; i++) expr_free(a[i]);
        r = e;
    }
    free(a);
    /* one spelling for a trig argument that is a multiple of Pi: the
     * natural Pi (a + b) (what the evaluator makes of Pi*(a + b), and what a
     * user writes), never Pi a + Pi b -- the bottom-up pass hands the engine
     * subexpressions, so its output must meet the outer input's spelling */
    if (r->type == EXPR_FUNCTION && r->data.function.arg_count == 1 &&
        r->data.function.head->type == EXPR_SYMBOL && has_trig_head(r->data.function.head->data.symbol.name) &&
        sf_contains_head(r->data.function.args[0], SYM_Plus)) {
        Expr* q = sf_mk1("Expand", sf_div(sf_cp(r->data.function.args[0]), sf_sym(SYM_Pi)));
        if (!sf_contains_symbol_pi(q)) {
            Expr* na = sf_times(sf_sym(SYM_Pi), q);
            if (!expr_eq(na, r->data.function.args[0])) {
                Expr* nr = eval_and_free(expr_new_function(sf_cp(r->data.function.head), (Expr*[]){ na }, 1));
                expr_free(r);
                r = nr;
            } else {
                expr_free(na);
            }
        } else {
            expr_free(q);
        }
    }
    bool is_sin = sf_is_call(r, SYM_Sin, 1), is_cos = sf_is_call(r, SYM_Cos, 1);
    if (!is_sin && !is_cos) return r;
    /* arg / Pi = q + t with q rational */
    Expr* u = sf_mk1("Expand", sf_div(sf_cp(r->data.function.args[0]), sf_sym(SYM_Pi)));
    Expr* t = NULL; int64_t qn, qd;
    bool ok = sf_split(u, &t, &qn, &qd) && t && (2 * qn) % qd == 0 && qn != 0;
    expr_free(u);
    if (!ok) { if (t) expr_free(t); return r; }
    int64_t k = (((2 * qn) / qd) % 4 + 4) % 4;               /* quarter turns */
    Expr* pt = sf_times(sf_sym(SYM_Pi), t);
    /* Sin: 0 Sin, 1 Cos, 2 -Sin, 3 -Cos;  Cos: 0 Cos, 1 -Sin, 2 -Cos, 3 Sin */
    bool to_cos = (k % 2 == 1) ? is_sin : is_cos;
    bool neg = is_sin ? (k == 2 || k == 3) : (k == 1 || k == 2);
    Expr* f = sf_mk1(to_cos ? SYM_Cos : SYM_Sin, pt);
    expr_free(r);
    return neg ? sf_times(sf_int(-1), f) : f;
}

/* Times[-1, Plus[a, b, ...]] -> Plus[-a, -b, ...] everywhere: the evaluator
 * keeps -(a + b) undistributed, so Log[-(-2 + 3 x)] and the engine's
 * Log[2 - 3 x] would be two kernels.  Owned in/out. */
static Expr* distribute_negation(Expr* e) {
    if (e->type != EXPR_FUNCTION) return e;
    size_t n = e->data.function.arg_count;
    bool sub = false;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        a[i] = distribute_negation(sf_cp(e->data.function.args[i]));
        if (a[i] != e->data.function.args[i]) sub = true;
    }
    Expr* r;
    if (sub) {
        r = eval_and_free(expr_new_function(sf_cp(e->data.function.head), a, n));
        expr_free(e);
    } else {
        for (size_t i = 0; i < n; i++) expr_free(a[i]);
        r = e;
    }
    free(a);
    if (sf_is_call(r, SYM_Times, 2) && r->data.function.args[0]->type == EXPR_INTEGER &&
        r->data.function.args[0]->data.integer == -1 &&
        sf_is_call(r->data.function.args[1], SYM_Plus, r->data.function.args[1]->type == EXPR_FUNCTION
                   ? r->data.function.args[1]->data.function.arg_count : 0)) {
        Expr* x = sf_mk1("Expand", sf_cp(r));
        expr_free(r);
        return x;
    }
    return r;
}

/* Odd functions: f[u] with the argument expanded and, when it reads
 * negative, f[u] -> -f[-u].  The evaluator folds Erf[-x] but not Erf[-I x]
 * or FresnelS[-1/3 - x], so the engine's -FresnelS[1/3 + x] and the input's
 * FresnelS[-1/3 - x] would otherwise be two kernels.  Owned in/out. */
static bool is_odd_head(const char* h) {
    return h == SYM_Erf || h == SYM_Erfi || h == SYM_InverseErf || h == SYM_FresnelS ||
           h == SYM_FresnelC || h == SYM_SinIntegral || h == SYM_SinhIntegral;
}

static Expr* normalize_odd_args(Expr* e) {
    if (e->type != EXPR_FUNCTION) return e;
    size_t n = e->data.function.arg_count;
    bool sub = false;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        a[i] = normalize_odd_args(sf_cp(e->data.function.args[i]));
        if (a[i] != e->data.function.args[i]) sub = true;
    }
    Expr* r;
    if (sub) {
        r = eval_and_free(expr_new_function(sf_cp(e->data.function.head), a, n));
        expr_free(e);
    } else {
        for (size_t i = 0; i < n; i++) expr_free(a[i]);
        r = e;
    }
    free(a);
    if (r->type != EXPR_FUNCTION || r->data.function.arg_count != 1 ||
        r->data.function.head->type != EXPR_SYMBOL || !is_odd_head(r->data.function.head->data.symbol.name))
        return r;
    const char* h = r->data.function.head->data.symbol.name;
    const Expr* u = r->data.function.args[0];
    if (u->type != EXPR_FUNCTION) return r;
    Expr* x = sf_mk1("Expand", sf_cp(u));
    Expr* out = NULL;
    if (sf_reads_negative(x)) {
        out = sf_times(sf_int(-1), sf_mk1(h, sf_mk1("Expand", sf_times(sf_int(-1), x))));
    } else if (!expr_eq(x, u)) {
        out = sf_mk1(h, x);
    } else {
        expr_free(x);
    }
    if (!out) return r;
    expr_free(r);
    return out;
}

/* Shift-argument and power normal form of a candidate (owned in/out). */
static Expr* normalize_form(Expr* r) {
    r = normalize_shift_args(r);
    r = normalize_odd_args(r);
    r = distribute_negation(r);
    r = fold_quarter_turns(r);
    Expr* t = transform_prime_rebase(r);
    if (t) { expr_free(r); r = t; }
    t = transform_power_oneify(r);
    if (t) { expr_free(r); r = t; }
    return normalize_numeric_powers(r);
}

/* Lowest-scoring of r, Together[r], Expand[r] and, when `factor` is set,
 * Factor[r] (owned in, owned out).  Together cancels over a common
 * denominator; Expand is needed when coefficients carry symbolic powers
 * (x^n (n/x) Y[n] and n x^(n-1) Y[n] only cancel once distributed); Factor
 * pulls the common kernel out of a collapsed rank-1 sum ((1 + a + n)
 * Gamma[a + n] / Gamma[a]).  Factor is skipped for rank-2 families: their
 * rewrites are large polynomials in two basis kernels, where it costs ~6x
 * the whole pass and never wins. */
static Expr* best_form(Expr* r, const Expr* cf, bool factor, size_t* score) {
    static const char* const forms[] = { "Together", "Expand", "Factor" };
    size_t nforms = factor ? 3 : 2;
    /* Put powers into the normal form simp_dispatch gives its input
     * (PrimeRebase, then PowerOneify: 4^u -> 2^(2u), 2^(1 + 2 x) -> 2 2^(2 x)),
     * so a prefactor the rewrite introduced meets the input's own spelling
     * of the same power and the two cancel under Together. */
    r = normalize_form(r);
    size_t best = score_with_func(r, cf);
    for (size_t i = 0; i < nforms; i++) {
        Expr* t = call_unary_copy(forms[i], r);
        if (!t) continue;
        /* Together/Factor map into arguments (Csc[Pi a + Pi b] ->
         * Csc[Pi (a + b)]): re-normalise, and if that exposes a cancellation
         * take Together once more */
        t = normalize_form(t);
        if (i != 0) {
            Expr* t2 = call_unary_copy("Together", t);
            if (t2) {
                t2 = normalize_form(t2);
                if (score_with_func(t2, cf) < score_with_func(t, cf)) { expr_free(t); t = t2; }
                else expr_free(t2);
            }
        }
        size_t st = score_with_func(t, cf);
        if (st < best) { expr_free(r); r = t; best = st; }
        else expr_free(t);
    }
    *score = best;
    return r;
}

/* Does the family need rewriting on its own (beyond partner reduction)? */
static bool fam_triggers(const Fam* f) {
    size_t ns = f->d->nshift;
    if (f->n < 2) return false;
    if (f->d->rank == 1) return true;
    if (f->n > 2) return true;
    /* two members: already a basis iff they differ by exactly e0 */
    int64_t diff0 = f->off[ns] - f->off[0];
    bool rest_eq = true;
    for (size_t j = 1; j < ns && rest_eq; j++) rest_eq = f->off[ns + j] == f->off[j];
    return !(rest_eq && (diff0 == 1 || diff0 == -1));
}

static bool fam_has_point(const Fam* f, const int64_t* k) {
    size_t ns = f->d->nshift;
    for (size_t i = 0; i < f->n; i++) {
        bool same = true;
        for (size_t j = 0; j < ns && same; j++) same = f->off[i * ns + j] == k[j];
        if (same) return true;
    }
    return false;
}

/* Same family key apart from the head (for partner lookup). */
static bool fam_key_eq_but_head(const Fam* f, const Fam* g) {
    if (f->d->arity != g->d->arity || f->d->nshift != g->d->nshift) return false;
    for (size_t a = 0; a < f->d->arity; a++)
        if (!is_shift_arg(f->d, a) &&
            !expr_eq(f->proto->data.function.args[a], g->proto->data.function.args[a]))
            return false;
    for (size_t j = 0; j < f->d->nshift; j++)
        if (f->den[j] != g->den[j] || f->cls[j] != g->cls[j] || !s_eq(f->s[j], g->s[j]))
            return false;
    return true;
}

/* Offsets of g's members re-expressed relative to f's base_num (so both
 * families share one coordinate system).  `shift[j]` is added to g's
 * offsets. */
static void partner_shift(const Fam* f, const Fam* g, int64_t* shift) {
    for (size_t j = 0; j < f->d->nshift; j++)
        shift[j] = (g->base_num[j] - f->base_num[j]) / f->den[j];
}

/* Candidate bases for f (and optionally partner g): member points, and for
 * rank 2 also member - e0.  Absolute offsets in f's coordinates, deduped. */
static size_t candidate_bases(const Fam* f, const Fam* g, int64_t* out, size_t cap) {
    size_t ns = f->d->nshift, nc = 0;
    int64_t gs[SF_MAX_SHIFT] = { 0 };
    if (g) partner_shift(f, g, gs);
    for (int pass = 0; pass < (g ? 2 : 1); pass++) {
        const Fam* h = pass ? g : f;
        for (size_t i = 0; i < h->n; i++)
            for (int back = 0; back <= (f->d->rank == 2); back++) {
                int64_t k[SF_MAX_SHIFT] = { 0 };
                for (size_t j = 0; j < ns; j++)
                    k[j] = h->off[i * ns + j] + (pass ? gs[j] : 0);
                k[0] -= back;
                bool dup = false;
                for (size_t q = 0; q < nc && !dup; q++) {
                    bool same = true;
                    for (size_t j = 0; j < ns && same; j++) same = out[q * ns + j] == k[j];
                    dup = same;
                }
                if (!dup && nc < cap) { memcpy(&out[nc * ns], k, sizeof(int64_t) * ns); nc++; }
            }
    }
    return nc;
}

/* g re-based into f's coordinate system (a shallow view; members shifted). */
static Fam fam_rebased(const Fam* f, const Fam* g, int64_t* offbuf) {
    Fam v = *g;
    int64_t gs[SF_MAX_SHIFT];
    partner_shift(f, g, gs);
    size_t ns = g->d->nshift;
    for (size_t i = 0; i < g->n; i++)
        for (size_t j = 0; j < ns; j++) offbuf[i * ns + j] = g->off[i * ns + j] + gs[j];
    for (size_t j = 0; j < ns; j++) v.base_num[j] = f->base_num[j];
    v.off = offbuf;
    return v;
}

/* Rewrite f (and partner g) onto basis b and, for a pair, reduce modulo the
 * pair relation. */
static Expr* rewrite_at(const Expr* cur, const Fam* f, const Fam* g, const int64_t* b) {
    Expr* r = rewrite_family(cur, f, b);
    if (!r || !g) return r;
    Expr* r2 = rewrite_family(r, g, b);
    expr_free(r);
    if (!r2) return NULL;
    int64_t b1[SF_MAX_SHIFT];
    memcpy(b1, b, sizeof b1);
    b1[0] += 1;
    Expr* c0 = atom_at(f, b);  Expr* c1 = atom_at(f, b1);
    Expr* d0 = atom_at(g, b);  Expr* d1 = atom_at(g, b1);
    Expr** a = args_at(f, b);
    Expr* R = f->d->pair_relation(a, c0, c1, d0, d1);
    args_free(a, f->d->arity);
    if (R) {
        Expr* vars[4] = { c1, d0, c0, d1 };
        Expr* red = sf_reduce_modulo(r2, R, vars, 4);
        expr_free(R);
        if (red) { expr_free(r2); r2 = red; }
    }
    expr_free(c0); expr_free(c1); expr_free(d0); expr_free(d1);
    return r2;
}

/* ---- canonicalisation ------------------------------------------------------- */

/* Canonical replacements are themselves walked again (a Beta becomes Gamma
 * atoms that may then reflect), to this depth; every canon hook is
 * idempotent on its own output, so the bound only guards a buggy hook. */
#define SF_CANON_DEPTH 4

static Expr* canon_walk_d(const Expr* e, const Expr* root, const AssumeCtx* ctx,
                          bool* changed, int depth);

static Expr* canon_walk(const Expr* e, const Expr* root, const AssumeCtx* ctx, bool* changed) {
    return canon_walk_d(e, root, ctx, changed, 0);
}

static Expr* canon_walk_d(const Expr* e, const Expr* root, const AssumeCtx* ctx,
                          bool* changed, int depth) {
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    bool sub = false;
    for (size_t i = 0; i < n; i++) {
        bool c = false;
        args[i] = canon_walk_d(e->data.function.args[i], root, ctx, &c, depth);
        sub |= c;
    }
    Expr* r = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    if (sub) { r = eval_and_free(r); *changed = true; }
    const SfDesc* d = desc_of(r);
    if (d && d->canon) {
        Expr* c = d->canon(r, root, ctx);
        if (c) {
            expr_free(r);
            r = eval_and_free(c);
            *changed = true;
            if (depth < SF_CANON_DEPTH) {
                bool again = false;
                Expr* r2 = canon_walk_d(r, root, ctx, &again, depth + 1);
                expr_free(r);
                r = again ? eval_and_free(r2) : r2;
            }
        }
    }
    return r;
}

/* Number of registered special-function atoms in e (any arity). */
static size_t sf_atom_count(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return 0;
    const Expr* h = e->data.function.head;
    size_t c = (h->type == EXPR_SYMBOL && head_present(h->data.symbol.name)) ? 1 : 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) c += sf_atom_count(e->data.function.args[i]);
    return c;
}

/* Does e hold a registered special-function atom with head h? */
static bool sf_has_head(const Expr* e, const char* h) {
    if (e->type != EXPR_FUNCTION) return false;
    const Expr* eh = e->data.function.head;
    if (eh->type == EXPR_SYMBOL && eh->data.symbol.name == h) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (sf_has_head(e->data.function.args[i], h)) return true;
    return false;
}

/* Does e contain the symbol named s? */
static bool sf_has_symbol(const Expr* e, const char* s) {
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == s;
    if (e->type != EXPR_FUNCTION) return false;
    if (sf_has_symbol(e->data.function.head, s)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (sf_has_symbol(e->data.function.args[i], s)) return true;
    return false;
}

/* Does some special-function head of o vanish from c entirely? */
static bool sf_head_eliminated(const Expr* c, const Expr* o) {
    if (o->type != EXPR_FUNCTION) return false;
    const Expr* h = o->data.function.head;
    if (h->type == EXPR_SYMBOL && head_present(h->data.symbol.name) &&
        !sf_has_head(c, h->data.symbol.name))
        return true;
    for (size_t i = 0; i < o->data.function.arg_count; i++)
        if (sf_head_eliminated(c, o->data.function.args[i])) return true;
    return false;
}

/* Is `c` (score cs) an improvement over `o` (score s0)?  A strict score drop
 * always is.  Eliminating special-function kernels is a simplification the
 * leaf score cannot see -- Gamma[1/3] Gamma[2/3] and 2 Pi/Sqrt[3] tie,
 * (n - k)/(1 + k) outscores Binomial[n, k + 1]/Binomial[n, k],
 * 1/x + 1/(1 + x) + 1/(2 + x) outscores PolyGamma[0, x + 3] - PolyGamma[0, x]
 * -- so a strict drop in the kernel count is accepted while the score stays
 * within double (plus 4) when every kernel is eliminated, and within half
 * (plus 4) otherwise (BesselJ[1/3, x] + BesselJ[7/3, x] -> 8 BesselJ[4/3, x]/(3 x)).
 * A drop that keeps every family and brings in Pi must not raise the score:
 * that is a kernel traded for constants, not removed.  Legendre's duplication
 * takes Gamma[-2 + 2 e]/Gamma[-2 + e] (score 17) to
 * 2^(-3 + 2 e) (-2 + e) Gamma[-1/2 + e]/Sqrt[Pi] (26), which the reference's
 * Simplify does not do, and FeynCalc compares prefactors in the Gamma-ratio
 * form (FCFeynmanParametrize ID11); three StruveH becoming one Gamma over
 * Sqrt[Pi] still counts, the whole Struve family going.  Both measures are
 * bounded below and one strictly decreases, so the rule cannot cycle. */
static bool sf_improves(const Expr* c, size_t cs, const Expr* o, size_t s0) {
    if (cs < s0) return true;
    size_t cn = sf_atom_count(c), n0 = sf_atom_count(o);
    if (cn == 0 && n0 > 0) return cs <= 2 * s0 + 4;   /* every kernel eliminated */
    if (cn >= n0) return false;
    if (cs <= s0) return true;
    if (cs > s0 + s0 / 2 + 4) return false;
    return sf_head_eliminated(c, o) || !(sf_has_symbol(c, SYM_Pi) && !sf_has_symbol(o, SYM_Pi));
}

/* Is f's partner family (rank-2 pair relation) present in fs? */
static bool has_partner(const FamSet* fs, const Fam* f) {
    if (!f->d->partner || !f->d->pair_relation) return false;
    for (size_t i = 0; i < fs->n; i++)
        if (fs->f[i].d->head == f->d->partner && fam_key_eq_but_head(f, &fs->f[i])) return true;
    return false;
}

/* Family-stage repetitions (see simp_special_functions). */
#define SF_FAMILY_ROUNDS 4

/* Candidate bases tried per family (each a full-expression rewrite). */
#define SF_MAX_BASES 64

/* The family stage.  Every family is rewritten onto the basis that scores
 * best for the whole expression.  Unconditionally (greedy = false) this
 * puts all families on bases at once -- the canonical form, needed when only
 * the combination cancels (Pascal's rule spans three Gamma families, each
 * rewrite alone raises the score).  Greedily (greedy = true) a family's
 * rewrite is adopted only when it improves the current form, which keeps a
 * collapsing family's win when another family's rewrite would cost more
 * than it saves.  Takes ownership of `cur`; returns the new form. */
static Expr* run_families(Expr* cur, size_t* cur_score, FamSet* fs, const Expr* cf,
                          bool greedy, bool* changed) {
    for (size_t fi = 0; fi < fs->n; fi++) fs->f[fi].done = false;
    for (size_t fi = 0; fi < fs->n; fi++) {
        Fam* f = &fs->f[fi];
        if (f->done) continue;
        f->done = true;
        /* partner (rank 2 only) */
        Fam* g = NULL;
        if (f->d->partner && f->d->pair_relation && f->d->rank == 2) {
            for (size_t gi = 0; gi < fs->n && !g; gi++)
                if (!fs->f[gi].done && fs->f[gi].d->head == f->d->partner &&
                    fam_key_eq_but_head(f, &fs->f[gi]))
                    g = &fs->f[gi];
        }
        if (!g) {
            /* the partner may own the Casoratian instead */
            for (size_t gi = 0; gi < fs->n && !g; gi++) {
                Fam* h = &fs->f[gi];
                if (!h->done && h->d->partner == f->d->head && h->d->pair_relation &&
                    h->d->rank == 2 && fam_key_eq_but_head(h, f)) {
                    g = f; f = h; f->done = true;
                }
            }
        }
        if (g) g->done = true;
        if (!g && !fam_triggers(f)) continue;

        int64_t* gbuf = NULL;
        Fam gv;
        if (g) {
            gbuf = malloc(sizeof(int64_t) * (g->n ? g->n : 1) * g->d->nshift);
            gv = fam_rebased(f, g, gbuf);
        }
        /* candidate_bases shifts the raw partner into f's coordinates itself;
         * the rebased view gv is what the rewrite consumes */
        int64_t bases[SF_MAX_BASES * SF_MAX_SHIFT];
        size_t nb = candidate_bases(f, g, bases, SF_MAX_BASES);
        Expr* fam_best = NULL;
        size_t fam_score = 0;
        size_t ns = f->d->nshift;
        for (size_t bi = 0; bi < nb; bi++) {
            const int64_t* b = &bases[bi * ns];
            if (f->d->rank == 2 && !g) {
                int64_t b1[SF_MAX_SHIFT];
                memcpy(b1, b, sizeof(int64_t) * ns);
                b1[0] += 1;
                if (!fam_has_point(f, b) && !fam_has_point(f, b1)) continue;
            }
            Expr* r = rewrite_at(cur, f, g ? &gv : NULL, b);
            if (!r) continue;
            size_t s;
            r = best_form(r, cf, f->d->rank == 1, &s);
            if (!fam_best || s < fam_score) {
                if (fam_best) expr_free(fam_best);
                fam_best = r; fam_score = s;
            } else {
                expr_free(r);
            }
        }
        free(gbuf);
        if (fam_best && (!greedy || fam_score <= *cur_score ||
                         sf_improves(fam_best, fam_score, cur, *cur_score))) {
            expr_free(cur);
            cur = fam_best; *cur_score = fam_score;
            *changed = true;
        } else if (fam_best) {
            expr_free(fam_best);
        }
    }
    return cur;
}

/* The family stage repeated on its own output while a family still
 * triggers: a rewrite can introduce atoms of another family (BarnesG[z+1] ->
 * Gamma[z] BarnesG[z]).  Family rewrites only ever emit heads lower in the
 * hierarchy (BarnesG -> Gamma -> none), so this terminates; the round bound
 * only guards a descriptor that breaks that rule.  Family atoms are borrowed
 * from a per-round `root` copy, matched structurally in the rewrites. */
static Expr* families_rounds(Expr* cur, size_t* cur_score, const Expr* cf,
                             bool greedy, bool* changed) {
    for (int round = 0; round < SF_FAMILY_ROUNDS; round++) {
        Expr* root = expr_copy(cur);
        FamSet fs = { NULL, 0, 0 };
        collect(root, &fs);
        bool any = false;
        for (size_t i = 0; i < fs.n && !any; i++)
            any = fam_triggers(&fs.f[i]) || has_partner(&fs, &fs.f[i]);
        bool ch = false;
        if (any) cur = run_families(cur, cur_score, &fs, cf, greedy, &ch);
        famset_free(&fs);
        bool same = expr_eq(root, cur);
        expr_free(root);
        if (!ch || same) break;
        *changed = true;
    }
    return cur;
}

/* ---- entry ------------------------------------------------------------------ */

static bool module_pass_relevant(const SfModule* mod, const Expr* e) {
    for (size_t i = 0; i < mod->npass_heads; i++)
        if (sf_contains_head(e, *mod->pass_heads[i])) return true;
    return false;
}

static bool has_trig(const Expr* e) {
    static const char** heads[] = { &SYM_Sin, &SYM_Cos, &SYM_Tan, &SYM_Cot, &SYM_Sec, &SYM_Csc };
    for (size_t i = 0; i < sizeof heads / sizeof heads[0]; i++)
        if (sf_contains_head(e, *heads[i])) return true;
    return false;
}

/* MATHILDA_NO_SF_ENGINE=1 disables the pass (read once): an A/B switch for
 * isolating the engine's effect on a Simplify result or its cost. */
static int sf_engine_disabled(void) {
    static int cached = -1;
    if (cached < 0) {
        const char* v = getenv("MATHILDA_NO_SF_ENGINE");
        cached = (v && *v && *v != '0') ? 1 : 0;
    }
    return cached;
}

static Expr* simp_special_functions_impl(const Expr* input, const AssumeCtx* ctx,
                                         const Expr* complexity_func);

/* Entry: the kill switch and the prefilter (one tree walk), then the
 * pipeline. */
Expr* simp_special_functions(const Expr* input, const AssumeCtx* ctx,
                             const Expr* complexity_func) {
    if (sf_engine_disabled()) return NULL;
    registry_init();
    if (!any_registered_head(input)) return NULL;
    return simp_special_functions_impl(input, ctx, complexity_func);
}

static Expr* simp_special_functions_impl(const Expr* input, const AssumeCtx* ctx,
                                         const Expr* complexity_func) {
    size_t in_score = score_with_func(input, complexity_func);

    /* 1. canonical generators */
    bool canon_changed = false;
    /* Iterated to a fixed point against the CURRENT form: one rewrite can
     * expose another's partner (the zeta functional equation introduces
     * Gamma[s] beside the input's Gamma[1 - s]). */
    Expr* start = normalize_odd_args(normalize_shift_args(expr_copy((Expr*)input)));  /* 1/2 (-4 + n) -> -2 + n/2 */
    Expr* cur = canon_walk(start, start, ctx, &canon_changed);
    if (start != input) canon_changed = true;
    expr_free(start);
    if (canon_changed) {
        cur = normalize_form(eval_and_free(cur));
        for (int it = 1; it < SF_CANON_DEPTH; it++) {
            bool again = false;
            Expr* nxt = canon_walk(cur, cur, ctx, &again);
            if (!again) { expr_free(nxt); break; }
            nxt = normalize_form(eval_and_free(nxt));
            bool same = expr_eq(nxt, cur);
            expr_free(cur);
            cur = nxt;
            if (same) break;
        }
    }
    size_t cur_score;
    if (canon_changed) cur = best_form(cur, complexity_func, false, &cur_score);  /* I - I E - I (1 - E) */
    else cur_score = score_with_func(cur, complexity_func);
    bool changed = canon_changed;

    /* 2. families: the canonical all-families form first; the greedy form
     * only when that does not pass the final gate (see run_families). */
    {
        size_t s_all = cur_score;
        bool ch_all = false;
        Expr* all = families_rounds(expr_copy(cur), &s_all, complexity_func, false, &ch_all);
        if (ch_all && sf_improves(all, s_all, input, in_score)) {
            expr_free(cur);
            cur = all; cur_score = s_all; changed = true;
        } else {
            if (ch_all && simp_debug_enabled()) {
                char tag[96];
                snprintf(tag, sizeof tag, "SpecialFunctions(all-families rejected: score %zu vs %zu)",
                         s_all, in_score);
                simp_debug_log(tag, input, all, 0.0);
            }
            expr_free(all);
            bool ch_g = false;
            cur = families_rounds(cur, &cur_score, complexity_func, true, &ch_g);
            changed |= ch_g;
        }
    }

    /* 3. whole-expression passes */
    for (size_t m = 0; m < SF_NMODULES; m++) {
        const SfModule* mod = SF_MODULES[m];
        if (!mod->npass || !module_pass_relevant(mod, cur)) continue;
        for (size_t p = 0; p < mod->npass; p++) {
            Expr* c = mod->passes[p](cur, input, ctx);
            if (!c) continue;
            size_t s;
            c = best_form(eval_and_free(c), complexity_func, true, &s);
            /* ties are accepted: a pass restores vocabulary (Gamma[2 + n] ->
             * Factorial[1 + n] when the input used Factorial) at equal score */
            if (s <= cur_score || sf_improves(c, s, cur, cur_score)) {
                expr_free(cur);
                cur = c; cur_score = s;
                changed = true;
            } else {
                expr_free(c);
            }
        }
    }

    /* private heads back to public ones before the result is scored */
    for (size_t m = 0; m < SF_NMODULES; m++)
        if (SF_MODULES[m]->finalize) {
            Expr* f = SF_MODULES[m]->finalize(cur);
            if (f != cur) { cur = f; cur_score = score_with_func(cur, complexity_func); }
        }

    /* 4. gate */
    if (!changed || !sf_improves(cur, cur_score, input, in_score)) {
        if (changed && simp_debug_enabled()) {
            char tag[96];
            snprintf(tag, sizeof tag, "SpecialFunctions(rejected: score %zu vs %zu, kernels %zu vs %zu)",
                     cur_score, in_score, sf_atom_count(cur), sf_atom_count(input));
            simp_debug_log(tag, input, cur, 0.0);
        }
        expr_free(cur);
        return NULL;
    }
    if (simp_debug_enabled()) simp_debug_log("SpecialFunctions", input, cur, 0.0);
    /* A reflection leaves a trig residue (Sec[Pi nu] Sin[Pi/2 + Pi nu]) that
     * only Simplify's bottom-up trig pass closes; anything else needs just
     * one more dispatch (every family is on a basis, so neither re-enters
     * the family rewrite). */
    Expr* fin;
    if (has_trig(cur)) {
        SimpMemo memo;
        simp_memo_init(&memo);
        fin = simp_bottomup(cur, ctx, complexity_func, &memo, 0);
        simp_memo_free(&memo);
    } else {
        fin = simp_dispatch(cur, ctx, complexity_func);
    }
    /* the trig pass may spell an argument Pi a + Pi b next to the input's
     * Pi (a + b); one more normal form lets Csc[u] Sin[u] cancel */
    if (fin) {
        Expr* nf = normalize_form(expr_copy(fin));
        if (score_with_func(nf, complexity_func) <= score_with_func(fin, complexity_func)) {
            expr_free(fin);
            fin = nf;
        } else {
            expr_free(nf);
        }
    }
    if (fin && (score_with_func(fin, complexity_func) <= cur_score ||
                sf_improves(fin, score_with_func(fin, complexity_func), cur, cur_score))) {
        expr_free(cur);
        return fin;
    }
    if (fin) expr_free(fin);
    return cur;
}
