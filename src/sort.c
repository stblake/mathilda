#include "sort.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "assoc.h"
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* ====================================================================
 * Canonical comparison of expressions (Mathematica-compatible `Order`).
 *
 * `expr_compare` lives here -- rather than in expr.c -- so the canonical
 * order policy is co-located with the user-facing Sort/OrderedQ builtins
 * that surface it.  The Orderless-attribute argument-sorting pipeline in
 * eval.c calls into this same comparator via expr.h.
 *
 * Algorithm, applied uniformly to every pair:
 *   1. Numeric atoms (Integer, Real, Rational, BigInt, MPFR) come first
 *      and compare by value.  Numeric < non-numeric.
 *   2. Strings come next, compared lexicographically.
 *   3. Everything else (Symbols and Functions) is compared by polynomial
 *      degree vector:
 *        - collect every Symbol name appearing in either expression,
 *        - sort those names in REVERSE alphabetical order so the
 *          lex-last variable is the most significant dimension,
 *        - for each name v, compute the polynomial degree of each side
 *          in v as a double (see expr_poly_degree),
 *        - lex-compare the degree vectors ascending.
 *      This makes `Plus[a x^2, b x]` sort as `b x + a x^2` and gives
 *      the multivariate grevlex-with-reverse-alpha display order on
 *      `Expand[(1 + x y^2 + y x^3)^3]`.
 *   4. Equal degree vectors fall back to a structural tiebreak: bare
 *      Symbol before compound, then head/arity/args recursively (each
 *      arg re-enters this same comparator).
 * ==================================================================== */

/* True for numeric atom kinds with an extractable double value. */
static bool is_atomic_numeric(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (is_rational((Expr*)e, NULL, NULL)) return true;
    return false;
}

static double get_numeric_value(const Expr* e) {
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    if (e->type == EXPR_BIGINT) return mpz_get_d(e->data.bigint);
    if (e->type == EXPR_REAL) return e->data.real;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
#endif
    int64_t n, d;
    if (is_rational((Expr*)e, &n, &d)) return (double)n / d;
    return 0;
}

/* Case-insensitive then case-sensitive lex (lowercase wins on case ties);
 * preserves the behaviour the previous expr.c implementation pinned for
 * string-vs-string ordering. */
static int string_compare_canonical(const char* s1, const char* s2) {
    const char* p1 = s1;
    const char* p2 = s2;
    while (*p1 && *p2) {
        char l1 = tolower((unsigned char)*p1);
        char l2 = tolower((unsigned char)*p2);
        if (l1 != l2) return (l1 < l2) ? -1 : 1;
        p1++; p2++;
    }
    if (*p1) return 1;
    if (*p2) return -1;

    p1 = s1; p2 = s2;
    while (*p1 && *p2) {
        if (*p1 != *p2) {
            if (islower((unsigned char)*p1) && !islower((unsigned char)*p2)) return -1;
            if (!islower((unsigned char)*p1) && islower((unsigned char)*p2)) return 1;
            return (*p1 < *p2) ? -1 : 1;
        }
        p1++; p2++;
    }
    return 0;
}

/* Return true iff `e` contains a Symbol with the interned name `sym`
 * anywhere in its tree (heads and args included). */
static bool expr_contains_symbol(const Expr* e, const char* sym) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol == sym;
    if (e->type != EXPR_FUNCTION) return false;
    if (expr_contains_symbol(e->data.function.head, sym)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (expr_contains_symbol(e->data.function.args[i], sym)) return true;
    }
    return false;
}

/* Polynomial degree of `e` in the variable `v`, as a double:
 *   0.0    -- e does not contain v.
 *   k > 0  -- e is v (k=1), Power[v, k] for positive integer or rational
 *             k, or a Times product whose direct factors sum to k.
 *   +inf   -- e contains v in a non-polynomial position: inside a
 *             compound that isn't bare-v / Power[v, positive rational]
 *             (e.g. Sin[v], Plus[..., v, ...], Power[v, negative],
 *             Power[non-v-containing-v, _]).
 *
 * Only DIRECT factors of a Times contribute to the sum; a non-poly
 * direct factor that itself contains v (e.g. Sin[v]) saturates the
 * result at +inf.  This matches Mathematica:
 *   - `Sin[x] + x`      displays as `x + Sin[x]`     (poly-in-x before non-poly).
 *   - `Sqrt[x] + x^2`   displays as `Sqrt[x] + x^2`  (deg 1/2 < 2).
 *   - `1/x + 1`         displays as `1 + 1/x`        (negative-exp x -> +inf).
 *   - `(3 + 6 a + 3 a^2) x` sorts by deg_x = 1 regardless of the Plus
 *      coefficient's internal a-degree. */
static double expr_poly_degree(const Expr* e, const char* v) {
    if (!e) return 0.0;
    if (e->type == EXPR_SYMBOL) return (e->data.symbol == v) ? 1.0 : 0.0;
    if (e->type != EXPR_FUNCTION) return 0.0;
    if (e->data.function.head->type != EXPR_SYMBOL) {
        return expr_contains_symbol(e, v) ? INFINITY : 0.0;
    }
    const char* h = e->data.function.head->data.symbol;
    if (h == SYM_Power && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol == v) {
            if (exp->type == EXPR_INTEGER) {
                int64_t k = exp->data.integer;
                if (k < 0) return INFINITY;
                return (double)k;
            }
            int64_t num, den;
            if (is_rational((Expr*)exp, &num, &den) && den > 0) {
                if (num < 0) return INFINITY;
                return (double)num / (double)den;
            }
            return INFINITY;
        }
        return (expr_contains_symbol(base, v) || expr_contains_symbol(exp, v))
               ? INFINITY : 0.0;
    }
    if (h == SYM_Times) {
        double sum = 0.0;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            double d = expr_poly_degree(e->data.function.args[i], v);
            if (isinf(d)) return INFINITY;
            sum += d;
        }
        return sum;
    }
    return expr_contains_symbol(e, v) ? INFINITY : 0.0;
}

/* Growable buffer of distinct interned-symbol-name pointers.  We need
 * collection to be ORDER-INDEPENDENT across a vs b: a silent truncation
 * after `a` would leave `b`'s unique symbols unrepresented, making the
 * resulting comparator asymmetric (compare(a,b) != -compare(b,a)) and
 * causing qsort on Orderless heads with hundreds of unknowns to
 * oscillate forever between two non-converging orderings. */
typedef struct {
    const char* stack_buf[64];  /* common-case avoid-malloc path */
    const char** items;          /* points at stack_buf or heap     */
    size_t count;
    size_t cap;
    bool ok;
} SymSet;

static void symset_init(SymSet* s) {
    s->items = s->stack_buf;
    s->count = 0;
    s->cap = sizeof(s->stack_buf) / sizeof(s->stack_buf[0]);
    s->ok = true;
}

static void symset_free(SymSet* s) {
    if (s->items && s->items != s->stack_buf) free((void*)s->items);
    s->items = NULL;
}

static bool symset_grow(SymSet* s) {
    size_t new_cap = s->cap * 2;
    const char** new_buf = (const char**)malloc(sizeof(const char*) * new_cap);
    if (!new_buf) { s->ok = false; return false; }
    memcpy(new_buf, s->items, sizeof(const char*) * s->count);
    if (s->items != s->stack_buf) free((void*)s->items);
    s->items = new_buf;
    s->cap = new_cap;
    return true;
}

static void symset_add(SymSet* s, const char* sym) {
    if (!s->ok) return;
    for (size_t i = 0; i < s->count; i++) {
        if (s->items[i] == sym) return;
    }
    if (s->count >= s->cap && !symset_grow(s)) return;
    s->items[s->count++] = sym;
}

/* Append every distinct Symbol-name pointer appearing in `e` into `s`,
 * DFS-visiting heads and args.  Grows the buffer as needed; never
 * truncates. */
static void collect_symbols_in(const Expr* e, SymSet* s) {
    if (!e || !s->ok) return;
    if (e->type == EXPR_SYMBOL) {
        symset_add(s, e->data.symbol);
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    collect_symbols_in(e->data.function.head, s);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        collect_symbols_in(e->data.function.args[i], s);
    }
}

/* qsort comparator that sorts char* pointers in REVERSE lex order, so
 * the lex-LAST variable name (e.g. `y` in {a, b, x, y}) ends up first
 * in the sorted array and drives the most significant dimension of the
 * degree-vector compare. */
static int symbol_reverse_cmp(const void* pa, const void* pb) {
    const char* a = *(const char* const*)pa;
    const char* b = *(const char* const*)pb;
    return strcmp(b, a);
}

int expr_compare(const Expr* a, const Expr* b) {
    if (a == b) return 0;
    if (!a) return -1;
    if (!b) return 1;

    /* 1. Numeric atoms. */
    bool a_atomic = is_atomic_numeric(a);
    bool b_atomic = is_atomic_numeric(b);
    if (a_atomic && b_atomic) {
        if (a->type == EXPR_INTEGER && b->type == EXPR_INTEGER) {
            if (a->data.integer < b->data.integer) return -1;
            if (a->data.integer > b->data.integer) return 1;
            return 0;
        }
        if (expr_is_integer_like(a) && expr_is_integer_like(b)) {
            mpz_t ma2, mb2;
            expr_to_mpz(a, ma2);
            expr_to_mpz(b, mb2);
            int cmp = mpz_cmp(ma2, mb2);
            mpz_clear(ma2);
            mpz_clear(mb2);
            return cmp;
        }
        double va = get_numeric_value(a);
        double vb = get_numeric_value(b);
        if (va < vb) return -1;
        if (va > vb) return 1;
        if (a->type != b->type) return (int)a->type - (int)b->type;
        return 0;
    }
    if (a_atomic) return -1;
    if (b_atomic) return 1;

    /* 2. Strings. */
    if (a->type == EXPR_STRING && b->type == EXPR_STRING) {
        return string_compare_canonical(a->data.string, b->data.string);
    }
    if (a->type == EXPR_STRING) return -1;
    if (b->type == EXPR_STRING) return 1;

    /* 3. Polynomial degree vector comparison. */
    SymSet vs;
    symset_init(&vs);
    collect_symbols_in(a, &vs);
    collect_symbols_in(b, &vs);

    qsort((void*)vs.items, vs.count, sizeof(const char*), symbol_reverse_cmp);

    int deg_cmp = 0;
    for (size_t i = 0; i < vs.count; i++) {
        double da = expr_poly_degree(a, vs.items[i]);
        double db = expr_poly_degree(b, vs.items[i]);
        if (da < db) { deg_cmp = -1; break; }
        if (da > db) { deg_cmp = 1; break; }
    }
    symset_free(&vs);
    if (deg_cmp != 0) return deg_cmp;

    /* 4. Structural tiebreak. */
    if (a->type == EXPR_SYMBOL && b->type == EXPR_SYMBOL) {
        if (a->data.symbol == b->data.symbol) return 0;
        return strcmp(a->data.symbol, b->data.symbol);
    }
    if (a->type == EXPR_SYMBOL) return -1;
    if (b->type == EXPR_SYMBOL) return 1;

    if (a->type == EXPR_FUNCTION && b->type == EXPR_FUNCTION) {
        int head_cmp = expr_compare(a->data.function.head, b->data.function.head);
        if (head_cmp != 0) return head_cmp;
        size_t ac = a->data.function.arg_count;
        size_t bc = b->data.function.arg_count;
        if (ac != bc) return (ac < bc) ? -1 : 1;
        for (size_t i = 0; i < ac; i++) {
            int c = expr_compare(a->data.function.args[i], b->data.function.args[i]);
            if (c != 0) return c;
        }
        return 0;
    }

    return 0;
}

/* ==================================================================== */

static Expr* current_sort_p = NULL;

static int custom_compare(const void* a, const void* b) {
    Expr* ea = *(Expr**)a;
    Expr* eb = *(Expr**)b;
    
    if (current_sort_p == NULL) {
        return expr_compare(ea, eb);
    }
    
    // Call p[ea, eb]
    Expr* args[2] = { expr_copy(ea), expr_copy(eb) };
    Expr* call = expr_new_function(expr_copy(current_sort_p), args, 2);
    Expr* result = evaluate(call);
    expr_free(call);
    
    int cmp = -1; // Default to "in order" (e1 < e2)
    
    if (result->type == EXPR_INTEGER) {
        // 1: e1 before e2 (e1 < e2)
        // 0: same
        // -1: e1 after e2 (e1 > e2)
        if (result->data.integer == 1) cmp = -1;
        else if (result->data.integer == 0) cmp = 0;
        else if (result->data.integer == -1) cmp = 1;
    } else if (result->type == EXPR_SYMBOL) {
        if (result->data.symbol == SYM_True) cmp = -1;
        else if (result->data.symbol == SYM_False) cmp = 1;
    }
    
    expr_free(result);
    return cmp;
}

Expr* builtin_sort(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) {
        return NULL;
    }
    
    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION) {
        return expr_copy(list);
    }

    /* Sort[assoc] orders the entries by their values (keys follow). */
    if (res->data.function.arg_count == 1 && is_association(list))
        return assoc_sort_by_value(list);

    Expr* p = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    size_t count = list->data.function.arg_count;
    if (count == 0) {
        return expr_copy(list);
    }
    
    Expr** sorted_args = malloc(sizeof(Expr*) * count);
    for (size_t i = 0; i < count; i++) {
        sorted_args[i] = expr_copy(list->data.function.args[i]);
    }
    
    // Set global context for comparison
    Expr* old_p = current_sort_p;
    current_sort_p = p;
    
    qsort(sorted_args, count, sizeof(Expr*), custom_compare);
    
    current_sort_p = old_p;
    
    Expr* result = expr_new_function(expr_copy(list->data.function.head), sorted_args, count);
    free(sorted_args);

    return result;
}

/* ------------------- SortBy ------------------- */

/* A payload paired with the sort key f[payload], for a key-once qsort
 * (avoids re-evaluating f during every comparison). */
typedef struct { Expr* payload; Expr* key; } SortByPair;

static int sortby_pair_cmp(const void* a, const void* b) {
    return expr_compare(((const SortByPair*)a)->key, ((const SortByPair*)b)->key);
}

/* SortBy[list, f]  — sort list elements by canonical order of f[element].
 * SortBy[assoc, f] — sort association entries by f[value] (keys follow).
 * SortBy[f]        — operator form: SortBy[f][expr] == SortBy[expr, f].
 * The key f[...] is evaluated once per element. */
Expr* builtin_sort_by(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 1) {
        /* Operator form: SortBy[f] -> Function[SortBy[#, f]]. */
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner = expr_new_function(expr_new_symbol(SYM_SortBy), inner_args, 2);
        Expr* func_args[1] = { inner };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }
    if (argc != 2) return NULL;

    Expr* coll = res->data.function.args[0];
    Expr* f    = res->data.function.args[1];
    if (coll->type != EXPR_FUNCTION) return expr_copy(coll);

    bool assoc = is_association(coll);
    size_t n = coll->data.function.arg_count;
    if (n == 0) return expr_copy(coll);

    SortByPair* pairs = malloc(sizeof(SortByPair) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* elem = coll->data.function.args[i];
        /* For an association, sort by f applied to the value; else by f[elem]. */
        Expr* subject = elem;
        if (assoc && elem->type == EXPR_FUNCTION && elem->data.function.arg_count == 2)
            subject = elem->data.function.args[1];
        Expr* call_args[1] = { expr_copy(subject) };
        Expr* call = expr_new_function(expr_copy(f), call_args, 1);
        pairs[i].payload = expr_copy(elem);
        pairs[i].key = evaluate(call);
        expr_free(call);
    }

    qsort(pairs, n, sizeof(SortByPair), sortby_pair_cmp);

    Expr** out = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) { out[i] = pairs[i].payload; expr_free(pairs[i].key); }
    free(pairs);

    Expr* result = expr_new_function(expr_copy(coll->data.function.head), out, n);
    free(out);
    return result;
}


/* ------------------- MaximalBy / MinimalBy ------------------- */

/* mode 0 = MaximalBy, mode 1 = MinimalBy.
 * MaximalBy[list, f]  -> the element(s) e maximising f[e] (all ties, in order).
 * MaximalBy[assoc, f] -> the entries whose value maximises f[value] (an assoc).
 * MaximalBy[f]        -> operator form. */
static Expr* maximal_minimal_by(Expr* res, int mode) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 1) {
        const char* head = (mode == 0) ? SYM_MaximalBy : SYM_MinimalBy;
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner = expr_new_function(expr_new_symbol(head), inner_args, 2);
        Expr* func_args[1] = { inner };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }
    if (argc != 2) return NULL;

    Expr* coll = res->data.function.args[0];
    Expr* f    = res->data.function.args[1];
    if (coll->type != EXPR_FUNCTION) return NULL;

    bool assoc = is_association(coll);
    size_t n = coll->data.function.arg_count;
    if (n == 0) return expr_copy(coll);  /* empty in, empty out */

    Expr** keys = malloc(sizeof(Expr*) * n);   /* f-value per element (owned) */
    for (size_t i = 0; i < n; i++) {
        Expr* elem = coll->data.function.args[i];
        Expr* subject = elem;
        if (assoc && elem->type == EXPR_FUNCTION && elem->data.function.arg_count == 2)
            subject = elem->data.function.args[1];
        Expr* call_args[1] = { expr_copy(subject) };
        Expr* call = expr_new_function(expr_copy(f), call_args, 1);
        keys[i] = evaluate(call);
        expr_free(call);
    }

    size_t best = 0;
    for (size_t i = 1; i < n; i++) {
        int c = expr_compare(keys[i], keys[best]);
        if ((mode == 0 && c > 0) || (mode == 1 && c < 0)) best = i;
    }

    Expr** out = malloc(sizeof(Expr*) * n);
    size_t nout = 0;
    for (size_t i = 0; i < n; i++)
        if (expr_compare(keys[i], keys[best]) == 0)
            out[nout++] = expr_copy(coll->data.function.args[i]);

    for (size_t i = 0; i < n; i++) expr_free(keys[i]);
    free(keys);

    Expr* result = expr_new_function(expr_copy(coll->data.function.head), out, nout);
    free(out);
    return result;
}

Expr* builtin_maximal_by(Expr* res) { return maximal_minimal_by(res, 0); }
Expr* builtin_minimal_by(Expr* res) { return maximal_minimal_by(res, 1); }

Expr* builtin_orderedq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) {
        return NULL;
    }
    
    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION) {
        return expr_new_symbol(SYM_True);
    }
    
    Expr* p = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    size_t count = list->data.function.arg_count;
    
    if (count <= 1) {
        return expr_new_symbol(SYM_True);
    }
    
    Expr* old_p = current_sort_p;
    current_sort_p = p;
    
    bool ordered = true;
    for (size_t i = 0; i < count - 1; i++) {
        Expr* ea = list->data.function.args[i];
        Expr* eb = list->data.function.args[i+1];
        
        int cmp = 0;
        if (current_sort_p == NULL) {
            cmp = expr_compare(ea, eb);
        } else {
            Expr* args[2] = { expr_copy(ea), expr_copy(eb) };
            Expr* call = expr_new_function(expr_copy(current_sort_p), args, 2);
            Expr* result = evaluate(call);
            expr_free(call);
            
            cmp = -1; // Default to "in order" (ea <= eb)
            if (result->type == EXPR_INTEGER) {
                if (result->data.integer == 1) cmp = -1;
                else if (result->data.integer == 0) cmp = 0;
                else if (result->data.integer == -1) cmp = 1;
            } else if (result->type == EXPR_SYMBOL) {
                if (result->data.symbol == SYM_True) cmp = -1;
                else if (result->data.symbol == SYM_False) cmp = 1;
            }
            expr_free(result);
        }
        
        if (cmp > 0) {
            ordered = false;
            break;
        }
    }
    
    current_sort_p = old_p;
    
    return expr_new_symbol(ordered ? "True" : "False");
}
