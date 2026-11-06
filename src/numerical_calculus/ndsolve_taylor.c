/* Mathilda — Taylor-mode automatic differentiation of an ODE right-hand side
 * (MPFR).  See ndsolve_taylor.h for the contract.
 *
 * Tape.  Each RHS expression is flattened into nodes evaluated in creation
 * order (children precede parents).  A node holds the Taylor coefficients
 * v[0..p] of its value along the solution.  Kinds and their recurrences, for
 * the k-th coefficient (a, b the operand series, c a constant):
 *
 *   CONST   c, 0, 0, ...                     T      t0, 1, 0, ...
 *   Y_i     the state's own coefficients     ADD    sum of operands
 *   SCALE   c a_k                            MUL    sum_{j=0..k} a_j b_{k-j}
 *   POW     u = a^c:  u_k = 1/(k a_0) sum_{j<k} (c (k-j) - j) a_{k-j} u_j
 *   EXP     e_k = 1/k sum_{j=1..k} j a_j e_{k-j}
 *   LOG     l_k = (a_k - 1/k sum_{j=1..k-1} j l_j a_{k-j}) / a_0
 *   SIN/COS s_k =  1/k sum j a_j c_{k-j},  c_k = -1/k sum j a_j s_{k-j}
 *   SINH/COSH  the same with both signs +
 *   ATAN    u = ArcTan[a], w = 1 + a^2:
 *           u_k = (k a_k - sum_{j=1..k-1} j u_j w_{k-j}) / (k w_0)
 *   ABS     u = Abs[a] = s a, s = Sign[a] on the side of the expansion point
 *           being integrated into (Sign[a_0], or, when a_0 = 0, the sign of the
 *           first nonzero a_k (-1)^k along a backward step).  The series is
 *           exact up to the next zero of a -- the kink, which the driver
 *           locates (nd_taylor_kink_*) and lands on, so no step straddles it.
 *
 * At a kink the driver marks the argument as exactly zero.  A fractional
 * power c > 0 of an argument that is then exactly zero, |a|^c ~ |a_1 s|^c, has
 * no Taylor series at all; its coefficients are set to those of the regular
 * part (zero) and nd_taylor_singular() reports it, so the driver crosses that
 * point with a Runge-Kutta step instead of the series.
 *
 * The state recurrence is then Y_{k+1} = F_k / (k + 1).  A positive integer
 * power is expanded into multiplications (so y^3 stays defined where y = 0,
 * which the a^c recurrence, dividing by a_0, is not). */
#include "ndsolve_taylor.h"

#ifdef USE_MPFR
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../sym_names.h"
#include "../numeric.h"
#include "../common.h"

#define RND MPFR_RNDN

enum {
    TN_CONST, TN_T, TN_Y, TN_ADD, TN_SCALE, TN_MUL, TN_POW,
    TN_EXP, TN_LOG, TN_SIN, TN_COS, TN_SINH, TN_COSH, TN_ATAN, TN_ABS
};

typedef struct {
    int    kind;
    int    a, b;        /* operand node indices (b: MUL rhs, ATAN's 1 + a^2) */
    int*   args;        /* ADD operands */
    int    nargs;
    int    yi;          /* TN_Y: state component */
    int    pair;        /* SIN/SINH: index of the companion COS/COSH node */
    mpfr_t c;           /* CONST value, SCALE factor, POW exponent */
    double cd;          /* c rounded to double (machine evaluation) */
    bool   has_c;
    int    sgn;         /* TN_ABS: sign in use (+-1), 0 while undetermined */
    bool   force_zero;  /* TN_ABS: the argument is exactly zero (a kink) here */
} TNode;

struct NdTaylorTape {
    long     bits;
    size_t   d;
    TNode*   nodes;
    int      n, cap;
    int*     out;       /* out[i]: node holding f_i */
    int*     ynode;     /* ynode[i]: the TN_Y node of component i */
    /* coefficient storage, n * (pcap + 1); dv: the machine twin */
    double*  dv;
    int      dpcap;
    mpfr_t*  v;
    int      pcap;
    int      vn;        /* number of nodes v was sized for */
    /* kinks: the TN_ABS nodes, the direction of integration, and whether the
     * last nd_taylor_coeffs met a fractional power of an exact zero */
    int*     absn; int nabs, capabs;
    int      dir;
    bool     any_forced;
    bool     singular;
    /* build context */
    Expr* const* ysym;
    const char*  tvar;
    /* build-time common-subexpression memo: a subexpression that occurs more
     * than once (Sin[p], rho[p]^2, ...) is one tape node, and a^n reuses the
     * powers of a already built (a^4 = (a^2)^2 shares a^2 with a^3). */
    const Expr** memo_e; uint64_t* memo_h; int* memo_n; int nmemo, capmemo;
    int* pw_base; int64_t* pw_exp; int* pw_node; int npw, cappw;
};

#define V(T, node, k) ((T)->v[(size_t)(node) * (size_t)((T)->pcap + 1) + (size_t)(k)])

static int tn_new(NdTaylorTape* T, int kind) {
    if (T->n == T->cap) {
        int nc = T->cap ? 2 * T->cap : 32;
        TNode* nn = realloc(T->nodes, sizeof(TNode) * (size_t)nc);
        if (!nn) return -1;
        T->nodes = nn; T->cap = nc;
    }
    TNode* nd = &T->nodes[T->n];
    memset(nd, 0, sizeof(*nd));
    nd->kind = kind; nd->a = nd->b = nd->pair = -1;
    return T->n++;
}

static int tn_const(NdTaylorTape* T, const mpfr_t val) {
    int i = tn_new(T, TN_CONST);
    if (i < 0) return -1;
    mpfr_init2(T->nodes[i].c, T->bits); T->nodes[i].has_c = true;
    mpfr_set(T->nodes[i].c, val, RND);
    return i;
}

static int tn_const_si(NdTaylorTape* T, long v) {
    mpfr_t x; mpfr_init2(x, T->bits); mpfr_set_si(x, v, RND);
    int i = tn_const(T, x);
    mpfr_clear(x);
    return i;
}

static int tn_unary(NdTaylorTape* T, int kind, int a) {
    if (a < 0) return -1;
    int i = tn_new(T, kind);
    if (i >= 0) T->nodes[i].a = a;
    return i;
}

static int tn_mul(NdTaylorTape* T, int a, int b) {
    if (a < 0 || b < 0) return -1;
    int i = tn_new(T, TN_MUL);
    if (i >= 0) { T->nodes[i].a = a; T->nodes[i].b = b; }
    return i;
}

static int tn_scale(NdTaylorTape* T, const mpfr_t c, int a) {
    if (a < 0) return -1;
    int i = tn_new(T, TN_SCALE);
    if (i < 0) return -1;
    T->nodes[i].a = a;
    mpfr_init2(T->nodes[i].c, T->bits); T->nodes[i].has_c = true;
    mpfr_set(T->nodes[i].c, c, RND);
    return i;
}

static int tn_pow(NdTaylorTape* T, int a, const mpfr_t alpha) {
    int i = tn_scale(T, alpha, a);
    if (i >= 0) T->nodes[i].kind = TN_POW;
    return i;
}

/* SIN/COS (or SINH/COSH) pair over operand a; returns the first (SIN) index. */
static int tn_pair(NdTaylorTape* T, int kind_first, int a) {
    int s = tn_unary(T, kind_first, a);
    if (s < 0) return -1;
    int c = tn_unary(T, kind_first == TN_SIN ? TN_COS : TN_COSH, a);
    if (c < 0) return -1;
    T->nodes[s].pair = c;
    return s;
}

/* True when e mentions t or a state symbol. */
static bool tn_depends(const NdTaylorTape* T, const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol.name == T->tvar) return true;
        for (size_t i = 0; i < T->d; i++)
            if (T->ysym[i]->data.symbol.name == e->data.symbol.name) return true;
        return false;
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (tn_depends(T, e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (tn_depends(T, e->data.function.args[i])) return true;
    return false;
}

/* Numeric value of a t/Y-free subexpression at the tape precision. */
static bool tn_eval_const(const NdTaylorTape* T, const Expr* e, mpfr_t out) {
    NumericSpec spec; spec.mode = NUMERIC_MODE_MPFR; spec.bits = T->bits;
    Expr* num = numericalize(e, spec);
    if (!num) return false;
    mpfr_t im; mpfr_init2(im, T->bits);
    bool inexact = false;
    bool ok = get_approx_mpfr(num, out, im, &inexact) && mpfr_number_p(out)
              && mpfr_zero_p(im);
    mpfr_clear(im);
    expr_free(num);
    return ok;
}

static int tn_build(NdTaylorTape* T, const Expr* e);

static int tn_pow_lookup(const NdTaylorTape* T, int a, int64_t n) {
    for (int i = 0; i < T->npw; i++)
        if (T->pw_base[i] == a && T->pw_exp[i] == n) return T->pw_node[i];
    return -1;
}

static void tn_pow_store(NdTaylorTape* T, int a, int64_t n, int node) {
    if (node < 0) return;
    if (T->npw == T->cappw) {
        int nc = T->cappw ? 2 * T->cappw : 16;
        int* b = realloc(T->pw_base, sizeof(int) * (size_t)nc);
        if (!b) return;
        T->pw_base = b;
        int64_t* e = realloc(T->pw_exp, sizeof(int64_t) * (size_t)nc);
        if (!e) return;
        T->pw_exp = e;
        int* nd = realloc(T->pw_node, sizeof(int) * (size_t)nc);
        if (!nd) return;
        T->pw_node = nd;
        T->cappw = nc;
    }
    T->pw_base[T->npw] = a; T->pw_exp[T->npw] = n; T->pw_node[T->npw] = node; T->npw++;
}

/* a^n for an integer n != 0: shared binary powering (n > 0) or a reciprocal. */
static int tn_int_power(NdTaylorTape* T, int a, int64_t n) {
    if (a < 0) return -1;
    if (n == 1) return a;
    int memo = tn_pow_lookup(T, a, n);
    if (memo >= 0) return memo;
    if (n < 0) {
        if (n == INT64_MIN) return -1;
        int q = tn_int_power(T, a, -n);
        mpfr_t m1; mpfr_init2(m1, T->bits); mpfr_set_si(m1, -1, RND);
        int r = tn_pow(T, q, m1);
        mpfr_clear(m1);
        return r;
    }
    if (n > 64) {   /* large power: the a^c recurrence (needs a_0 != 0) */
        mpfr_t al; mpfr_init2(al, T->bits); mpfr_set_si(al, (long)n, RND);
        int r = tn_pow(T, a, al);
        mpfr_clear(al);
        return r;
    }
    int r;
    if (n % 2 == 0) { int q = tn_int_power(T, a, n / 2); r = tn_mul(T, q, q); }
    else            r = tn_mul(T, tn_int_power(T, a, n - 1), a);
    tn_pow_store(T, a, n, r);
    return r;
}

static int tn_build_node(NdTaylorTape* T, const Expr* e);

/* Abs[a]: a kink node, registered so a step ends on a zero of a. */
static int tn_abs(NdTaylorTape* T, int a) {
    int r = tn_unary(T, TN_ABS, a);
    if (r < 0) return -1;
    if (T->nabs == T->capabs) {
        int nc = T->capabs ? 2 * T->capabs : 4;
        int* na = realloc(T->absn, sizeof(int) * (size_t)nc);
        if (!na) return -1;
        T->absn = na; T->capabs = nc;
    }
    T->absn[T->nabs++] = r;
    return r;
}

/* a + b + s |a - b| (s = +1 for Max, -1 for Min), halved. */
static int tn_maxmin(NdTaylorTape* T, int a, int b, bool is_max) {
    if (a < 0 || b < 0) return -1;
    mpfr_t c; mpfr_init2(c, T->bits);
    mpfr_set_si(c, -1, RND);
    int nb = tn_scale(T, c, b);
    int dif = tn_new(T, TN_ADD);
    int* da = malloc(sizeof(int) * 2);
    int r = -1;
    if (nb >= 0 && dif >= 0 && da) {
        da[0] = a; da[1] = nb;
        T->nodes[dif].args = da; T->nodes[dif].nargs = 2; da = NULL;
        int ab = tn_abs(T, dif);
        if (ab >= 0 && !is_max) ab = tn_scale(T, c, ab);
        int sum = ab >= 0 ? tn_new(T, TN_ADD) : -1;
        int* sa = sum >= 0 ? malloc(sizeof(int) * 3) : NULL;
        if (sa) {
            sa[0] = a; sa[1] = b; sa[2] = ab;
            T->nodes[sum].args = sa; T->nodes[sum].nargs = 3;
            mpfr_set_d(c, 0.5, RND);
            r = tn_scale(T, c, sum);
        }
    }
    free(da);
    mpfr_clear(c);
    return r;
}

static int tn_build(NdTaylorTape* T, const Expr* e) {
    uint64_t h = expr_hash(e);
    for (int i = 0; i < T->nmemo; i++)
        if (T->memo_h[i] == h && expr_eq(T->memo_e[i], e)) return T->memo_n[i];
    int r = tn_build_node(T, e);
    if (r < 0) return r;
    if (T->nmemo == T->capmemo) {
        int nc = T->capmemo ? 2 * T->capmemo : 32;
        const Expr** me = realloc(T->memo_e, sizeof(Expr*) * (size_t)nc);
        if (!me) return r;
        T->memo_e = me;
        uint64_t* mh = realloc(T->memo_h, sizeof(uint64_t) * (size_t)nc);
        if (!mh) return r;
        T->memo_h = mh;
        int* mn = realloc(T->memo_n, sizeof(int) * (size_t)nc);
        if (!mn) return r;
        T->memo_n = mn;
        T->capmemo = nc;
    }
    T->memo_e[T->nmemo] = e; T->memo_h[T->nmemo] = h; T->memo_n[T->nmemo] = r; T->nmemo++;
    return r;
}

static int tn_build_node(NdTaylorTape* T, const Expr* e) {
    if (!tn_depends(T, e)) {
        mpfr_t c; mpfr_init2(c, T->bits);
        int r = tn_eval_const(T, e, c) ? tn_const(T, c) : -1;
        mpfr_clear(c);
        return r;
    }
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol.name == T->tvar) return tn_new(T, TN_T);
        for (size_t i = 0; i < T->d; i++)
            if (T->ysym[i]->data.symbol.name == e->data.symbol.name) return T->ynode[i];
        return -1;
    }
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return -1;
    const char* h = e->data.function.head->data.symbol.name;
    size_t n = e->data.function.arg_count;
    Expr* const* A = e->data.function.args;

    if (h == SYM_Plus && n >= 1) {
        int* args = malloc(sizeof(int) * n);
        if (!args) return -1;
        for (size_t i = 0; i < n; i++) {
            args[i] = tn_build(T, A[i]);
            if (args[i] < 0) { free(args); return -1; }
        }
        int r = tn_new(T, TN_ADD);
        if (r < 0) { free(args); return -1; }
        T->nodes[r].args = args; T->nodes[r].nargs = (int)n;
        return r;
    }
    if (h == SYM_Times && n >= 1) {
        /* Fold the constant factors into one scale, multiply the rest. */
        mpfr_t c; mpfr_init2(c, T->bits); mpfr_set_ui(c, 1, RND);
        bool scaled = false;
        int prod = -1;
        for (size_t i = 0; i < n; i++) {
            if (!tn_depends(T, A[i])) {
                mpfr_t f; mpfr_init2(f, T->bits);
                bool ok = tn_eval_const(T, A[i], f);
                if (ok) mpfr_mul(c, c, f, RND);
                mpfr_clear(f);
                if (!ok) { mpfr_clear(c); return -1; }
                scaled = true;
                continue;
            }
            int x = tn_build(T, A[i]);
            if (x < 0) { mpfr_clear(c); return -1; }
            prod = (prod < 0) ? x : tn_mul(T, prod, x);
            if (prod < 0) { mpfr_clear(c); return -1; }
        }
        int r = (scaled && prod >= 0) ? tn_scale(T, c, prod) : prod;
        mpfr_clear(c);
        return r;
    }
    if (h == SYM_Power && n == 2) {
        const Expr* b = A[0];
        const Expr* x = A[1];
        if (!tn_depends(T, x)) {
            if (x->type == EXPR_INTEGER) {
                if (x->data.integer == 0) return tn_const_si(T, 1);
                return tn_int_power(T, tn_build(T, b), x->data.integer);
            }
            mpfr_t al; mpfr_init2(al, T->bits);
            int r = -1;
            if (tn_eval_const(T, x, al)) {
                /* an integer-valued exponent (3., 6/2) is a product, defined
                 * where the base vanishes */
                if (mpfr_integer_p(al) && mpfr_fits_slong_p(al, RND) && !mpfr_zero_p(al))
                    r = tn_int_power(T, tn_build(T, b), (int64_t)mpfr_get_si(al, RND));
                else r = tn_pow(T, tn_build(T, b), al);
            }
            mpfr_clear(al);
            return r;
        }
        /* variable exponent: b^x = Exp[x Log[b]] (b = E: Exp[x]) */
        int xe = tn_build(T, x);
        if (b->type == EXPR_SYMBOL && b->data.symbol.name == SYM_E)
            return tn_unary(T, TN_EXP, xe);
        if (!tn_depends(T, b)) {
            mpfr_t c; mpfr_init2(c, T->bits);
            int r = -1;
            if (tn_eval_const(T, b, c) && mpfr_sgn(c) > 0) {
                mpfr_log(c, c, RND);
                r = tn_unary(T, TN_EXP, tn_scale(T, c, xe));
            }
            mpfr_clear(c);
            return r;
        }
        return tn_unary(T, TN_EXP, tn_mul(T, xe, tn_unary(T, TN_LOG, tn_build(T, b))));
    }
    /* Max / Min fold pairwise through the kink of Abs:
     *   Max[a, b] = (a + b + Abs[a - b])/2,  Min[a, b] = (a + b - Abs[a - b])/2,
     * so a step ends where a - b changes sign and each branch is the smooth
     * series of the larger (smaller) argument -- Max[y, 0]^(3/2) keeps the
     * high-order series method instead of falling to Runge-Kutta, whose
     * MPFR steps (41728 of them where the series takes 469) dominated a
     * Thomas-Fermi shooting run approaching its movable singularity. */
    if ((h == SYM_Max || h == SYM_Min) && n >= 2) {
        int acc = tn_build(T, A[0]);
        for (size_t i = 1; i < n && acc >= 0; i++)
            acc = tn_maxmin(T, acc, tn_build(T, A[i]), h == SYM_Max);
        return acc;
    }
    if (n != 1) return -1;
    if (h == SYM_Exp) return tn_unary(T, TN_EXP, tn_build(T, A[0]));
    if (h == SYM_Abs) return tn_abs(T, tn_build(T, A[0]));
    if (h == SYM_Log) return tn_unary(T, TN_LOG, tn_build(T, A[0]));
    if (h == SYM_Sqrt) {
        mpfr_t al; mpfr_init2(al, T->bits); mpfr_set_d(al, 0.5, RND);
        int r = tn_pow(T, tn_build(T, A[0]), al);
        mpfr_clear(al);
        return r;
    }
    if (h == SYM_Sin || h == SYM_Cos || h == SYM_Tan) {
        int s = tn_pair(T, TN_SIN, tn_build(T, A[0]));
        if (s < 0) return -1;
        if (h == SYM_Sin) return s;
        int c = T->nodes[s].pair;
        if (h == SYM_Cos) return c;
        return tn_mul(T, s, tn_int_power(T, c, -1));
    }
    if (h == SYM_Sinh || h == SYM_Cosh || h == SYM_Tanh) {
        int s = tn_pair(T, TN_SINH, tn_build(T, A[0]));
        if (s < 0) return -1;
        if (h == SYM_Sinh) return s;
        int c = T->nodes[s].pair;
        if (h == SYM_Cosh) return c;
        return tn_mul(T, s, tn_int_power(T, c, -1));
    }
    if (h == SYM_ArcTan) {
        int a = tn_build(T, A[0]);
        if (a < 0) return -1;
        int sq = tn_mul(T, a, a);
        int one = tn_const_si(T, 1);
        if (sq < 0 || one < 0) return -1;
        int w = tn_new(T, TN_ADD);
        if (w < 0) return -1;
        T->nodes[w].args = malloc(sizeof(int) * 2);
        if (!T->nodes[w].args) return -1;
        T->nodes[w].args[0] = one; T->nodes[w].args[1] = sq; T->nodes[w].nargs = 2;
        int r = tn_unary(T, TN_ATAN, a);
        if (r >= 0) T->nodes[r].b = w;
        return r;
    }
    return -1;
}

void nd_taylor_free(NdTaylorTape* T) {
    if (!T) return;
    for (int i = 0; i < T->n; i++) {
        if (T->nodes[i].has_c) mpfr_clear(T->nodes[i].c);
        free(T->nodes[i].args);
    }
    free(T->nodes);
    if (T->v) {
        size_t nv = (size_t)T->vn * (size_t)(T->pcap + 1);
        for (size_t i = 0; i < nv; i++) mpfr_clear(T->v[i]);
        free(T->v);
    }
    free(T->dv);
    free(T->memo_e); free(T->memo_h); free(T->memo_n);
    free(T->pw_base); free(T->pw_exp); free(T->pw_node);
    free(T->out); free(T->ynode); free(T->absn);
    free(T);
}

NdTaylorTape* nd_taylor_build(Expr* const* f, size_t d, Expr* const* ysym,
                              const char* tvar, long bits) {
    if (d == 0) return NULL;
    NdTaylorTape* T = calloc(1, sizeof(NdTaylorTape));
    if (!T) return NULL;
    T->bits = bits; T->d = d; T->ysym = ysym; T->tvar = tvar; T->dir = 1;
    T->out = malloc(sizeof(int) * d);
    T->ynode = malloc(sizeof(int) * d);
    if (!T->out || !T->ynode) { nd_taylor_free(T); return NULL; }
    for (size_t i = 0; i < d; i++) {
        if (!ysym[i] || ysym[i]->type != EXPR_SYMBOL) { nd_taylor_free(T); return NULL; }
        int y = tn_new(T, TN_Y);
        if (y < 0) { nd_taylor_free(T); return NULL; }
        T->nodes[y].yi = (int)i;
        T->ynode[i] = y;
    }
    for (size_t i = 0; i < d; i++) {
        T->out[i] = tn_build(T, f[i]);
        if (T->out[i] < 0) { nd_taylor_free(T); return NULL; }
    }
    T->ysym = NULL;   /* build-only context */
    free(T->memo_e); free(T->memo_h); free(T->memo_n);
    T->memo_e = NULL; T->memo_h = NULL; T->memo_n = NULL; T->nmemo = T->capmemo = 0;
    free(T->pw_base); free(T->pw_exp); free(T->pw_node);
    T->pw_base = NULL; T->pw_exp = NULL; T->pw_node = NULL; T->npw = T->cappw = 0;
    for (int i = 0; i < T->n; i++)
        if (T->nodes[i].has_c) T->nodes[i].cd = mpfr_get_d(T->nodes[i].c, RND);
    return T;
}

static bool tn_reserve(NdTaylorTape* T, int p) {
    if (T->v && p <= T->pcap && T->vn == T->n) return true;
    if (T->v) {
        size_t nv = (size_t)T->vn * (size_t)(T->pcap + 1);
        for (size_t i = 0; i < nv; i++) mpfr_clear(T->v[i]);
        free(T->v);
        T->v = NULL;
    }
    int pc = p < 8 ? 8 : p;
    size_t nv = (size_t)T->n * (size_t)(pc + 1);
    T->v = malloc(sizeof(mpfr_t) * nv);
    if (!T->v) return false;
    for (size_t i = 0; i < nv; i++) mpfr_init2(T->v[i], T->bits);
    T->pcap = pc; T->vn = T->n;
    return true;
}

bool nd_taylor_coeffs(NdTaylorTape* T, const mpfr_t t0, mpfr_t* Y0, int p,
                      mpfr_t* C) {
    if (!T || p < 1 || !tn_reserve(T, p)) return false;
    size_t d = T->d;
    mpfr_t acc, tmp, s2;
    mpfr_inits2(T->bits, acc, tmp, s2, (mpfr_ptr)0);

    /* Constant and t nodes are fixed series. */
    for (int n = 0; n < T->n; n++) {
        TNode* nd = &T->nodes[n];
        if (nd->kind == TN_CONST) {
            mpfr_set(V(T, n, 0), nd->c, RND);
            for (int k = 1; k <= p; k++) mpfr_set_zero(V(T, n, k), 1);
        } else if (nd->kind == TN_T) {
            mpfr_set(V(T, n, 0), t0, RND);
            mpfr_set_ui(V(T, n, 1), 1, RND);
            for (int k = 2; k <= p; k++) mpfr_set_zero(V(T, n, k), 1);
        }
    }
    for (size_t i = 0; i < d; i++) mpfr_set(C[i], Y0[i], RND);
    T->singular = false;
    T->any_forced = false;
    for (int j = 0; j < T->nabs; j++) if (T->nodes[T->absn[j]].force_zero) T->any_forced = true;

    bool ok = true;
    for (int k = 0; k < p && ok; k++) {
        for (size_t i = 0; i < d; i++) mpfr_set(V(T, T->ynode[i], k), C[(size_t)k * d + i], RND);
        for (int n = 0; n < T->n && ok; n++) {
            TNode* nd = &T->nodes[n];
            mpfr_ptr out = V(T, n, k);
            switch (nd->kind) {
            case TN_CONST: case TN_T: case TN_Y: case TN_COS: case TN_COSH:
                break;   /* fixed, supplied, or computed with its pair */
            case TN_ADD:
                mpfr_set(out, V(T, nd->args[0], k), RND);
                for (int j = 1; j < nd->nargs; j++) mpfr_add(out, out, V(T, nd->args[j], k), RND);
                break;
            case TN_SCALE:
                mpfr_mul(out, nd->c, V(T, nd->a, k), RND);
                break;
            case TN_MUL:
                mpfr_set_zero(acc, 1);
                for (int j = 0; j <= k; j++) {
                    mpfr_mul(tmp, V(T, nd->a, j), V(T, nd->b, k - j), RND);
                    mpfr_add(acc, acc, tmp, RND);
                }
                mpfr_set(out, acc, RND);
                break;
            case TN_POW: {
                mpfr_ptr a0 = V(T, nd->a, 0);
                if (mpfr_zero_p(a0) && T->any_forced && mpfr_sgn(nd->c) > 0
                    && !mpfr_integer_p(nd->c)) {
                    /* |a_1 s|^c at a kink: no series; the regular part is 0 */
                    mpfr_set_zero(out, 1);
                    T->singular = true;
                    break;
                }
                if (k == 0) { mpfr_pow(out, a0, nd->c, RND); break; }
                if (mpfr_zero_p(a0) && mpfr_sgn(nd->c) > 0) {
                    /* A base that vanishes identically (to this order) has
                     * the zero series as its positive power: Max[y, 0]^c on
                     * the branch y < 0, where the base is y + Abs[y] = 0
                     * exactly. */
                    bool allz = true;
                    for (int j = 1; j <= k && allz; j++)
                        if (!mpfr_zero_p(V(T, nd->a, j))) allz = false;
                    if (allz) { mpfr_set_zero(out, 1); break; }
                }
                if (mpfr_zero_p(a0)) { ok = false; break; }
                mpfr_set_zero(acc, 1);
                for (int j = 0; j < k; j++) {
                    /* (c (k - j) - j) a_{k-j} u_j */
                    mpfr_mul_si(tmp, nd->c, k - j, RND);
                    mpfr_sub_si(tmp, tmp, j, RND);
                    mpfr_mul(tmp, tmp, V(T, nd->a, k - j), RND);
                    mpfr_mul(tmp, tmp, V(T, n, j), RND);
                    mpfr_add(acc, acc, tmp, RND);
                }
                mpfr_mul_si(tmp, a0, k, RND);
                mpfr_div(out, acc, tmp, RND);
                break;
            }
            case TN_EXP:
                if (k == 0) { mpfr_exp(out, V(T, nd->a, 0), RND); break; }
                mpfr_set_zero(acc, 1);
                for (int j = 1; j <= k; j++) {
                    mpfr_mul(tmp, V(T, nd->a, j), V(T, n, k - j), RND);
                    mpfr_mul_si(tmp, tmp, j, RND);
                    mpfr_add(acc, acc, tmp, RND);
                }
                mpfr_div_si(out, acc, k, RND);
                break;
            case TN_LOG: {
                mpfr_ptr a0 = V(T, nd->a, 0);
                if (k == 0) {
                    if (mpfr_sgn(a0) <= 0) { ok = false; break; }
                    mpfr_log(out, a0, RND);
                    break;
                }
                mpfr_set_zero(acc, 1);
                for (int j = 1; j < k; j++) {
                    mpfr_mul(tmp, V(T, n, j), V(T, nd->a, k - j), RND);
                    mpfr_mul_si(tmp, tmp, j, RND);
                    mpfr_add(acc, acc, tmp, RND);
                }
                mpfr_div_si(acc, acc, k, RND);
                mpfr_sub(acc, V(T, nd->a, k), acc, RND);
                mpfr_div(out, acc, a0, RND);
                break;
            }
            case TN_SIN: case TN_SINH: {
                int c = nd->pair;
                bool hyp = (nd->kind == TN_SINH);
                if (k == 0) {
                    if (hyp) mpfr_sinh_cosh(out, V(T, c, 0), V(T, nd->a, 0), RND);
                    else     mpfr_sin_cos(out, V(T, c, 0), V(T, nd->a, 0), RND);
                    break;
                }
                mpfr_set_zero(acc, 1); mpfr_set_zero(s2, 1);
                for (int j = 1; j <= k; j++) {
                    mpfr_mul_si(tmp, V(T, nd->a, j), j, RND);
                    mpfr_ptr cc = V(T, c, k - j), ss = V(T, n, k - j);
                    mpfr_fma(acc, tmp, cc, acc, RND);
                    mpfr_fma(s2, tmp, ss, s2, RND);
                }
                mpfr_div_si(out, acc, k, RND);
                mpfr_div_si(V(T, c, k), s2, hyp ? k : -k, RND);
                break;
            }
            case TN_ABS: {
                mpfr_ptr ak = V(T, nd->a, k);
                if (k == 0) {
                    nd->sgn = (nd->force_zero || mpfr_zero_p(ak)) ? 0 : mpfr_sgn(ak);
                    if (nd->sgn == 0) mpfr_set_zero(out, 1);
                    else mpfr_abs(out, ak, RND);
                    break;
                }
                if (nd->sgn == 0 && !mpfr_zero_p(ak)) {
                    int sg = mpfr_sgn(ak);
                    nd->sgn = (T->dir < 0 && (k & 1)) ? -sg : sg;
                }
                if (nd->sgn == 0) mpfr_set_zero(out, 1);
                else if (nd->sgn > 0) mpfr_set(out, ak, RND);
                else mpfr_neg(out, ak, RND);
                break;
            }
            case TN_ATAN: {
                if (k == 0) { mpfr_atan(out, V(T, nd->a, 0), RND); break; }
                mpfr_mul_si(acc, V(T, nd->a, k), k, RND);
                for (int j = 1; j < k; j++) {
                    mpfr_mul(tmp, V(T, n, j), V(T, nd->b, k - j), RND);
                    mpfr_mul_si(tmp, tmp, j, RND);
                    mpfr_sub(acc, acc, tmp, RND);
                }
                mpfr_mul_si(tmp, V(T, nd->b, 0), k, RND);
                mpfr_div(out, acc, tmp, RND);
                break;
            }
            default:
                ok = false;
            }
            if (ok && !mpfr_number_p(out)) ok = false;
        }
        if (!ok) break;
        for (size_t i = 0; i < d; i++) {
            mpfr_ptr dst = C[(size_t)(k + 1) * d + i];
            mpfr_div_si(dst, V(T, T->out[i], k), k + 1, RND);
            if (!mpfr_number_p(dst)) ok = false;
        }
    }
    mpfr_clears(acc, tmp, s2, (mpfr_ptr)0);
    for (int j = 0; j < T->nabs; j++) T->nodes[T->absn[j]].force_zero = false;
    return ok;
}

void nd_taylor_set_direction(NdTaylorTape* T, int dir) { if (T) T->dir = dir < 0 ? -1 : 1; }
int  nd_taylor_kink_count(const NdTaylorTape* T) { return T ? T->nabs : 0; }
int  nd_taylor_kink_sign(const NdTaylorTape* T, int j) { return T->nodes[T->absn[j]].sgn; }
mpfr_ptr nd_taylor_kink_coeff(NdTaylorTape* T, int j, int k) {
    return V(T, T->nodes[T->absn[j]].a, k);
}
void nd_taylor_force_kink_zero(NdTaylorTape* T, int j) { T->nodes[T->absn[j]].force_zero = true; }
bool nd_taylor_singular(const NdTaylorTape* T) { return T && T->singular; }

/* ------------------------------------------------------------------ *
 *  Machine-precision evaluation of the same tape                      *
 * ------------------------------------------------------------------ */
#define DV(T, node, k) ((T)->dv[(size_t)(node) * (size_t)((T)->dpcap + 1) + (size_t)(k)])

bool nd_taylor_coeffs_d(NdTaylorTape* T, double t0, const double* Y0, int p,
                        double* C) {
    if (!T || p < 1) return false;
    if (!T->dv || p > T->dpcap) {
        free(T->dv);
        T->dpcap = p < 4 ? 4 : p;
        T->dv = malloc(sizeof(double) * (size_t)T->n * (size_t)(T->dpcap + 1));
        if (!T->dv) { T->dpcap = 0; return false; }
    }
    size_t d = T->d;
    for (int n = 0; n < T->n; n++) {
        TNode* nd = &T->nodes[n];
        if (nd->kind == TN_CONST) {
            DV(T, n, 0) = nd->cd;
            for (int k = 1; k <= p; k++) DV(T, n, k) = 0.0;
        } else if (nd->kind == TN_T) {
            DV(T, n, 0) = t0; DV(T, n, 1) = 1.0;
            for (int k = 2; k <= p; k++) DV(T, n, k) = 0.0;
        }
    }
    for (size_t i = 0; i < d; i++) C[i] = Y0[i];
    for (int k = 0; k < p; k++) {
        for (size_t i = 0; i < d; i++) DV(T, T->ynode[i], k) = C[(size_t)k * d + i];
        for (int n = 0; n < T->n; n++) {
            TNode* nd = &T->nodes[n];
            double r = 0.0;
            switch (nd->kind) {
            case TN_CONST: case TN_T: case TN_Y: case TN_COS: case TN_COSH:
                continue;
            case TN_ADD:
                for (int j = 0; j < nd->nargs; j++) r += DV(T, nd->args[j], k);
                break;
            case TN_SCALE: r = nd->cd * DV(T, nd->a, k); break;
            case TN_MUL:
                for (int j = 0; j <= k; j++) r += DV(T, nd->a, j) * DV(T, nd->b, k - j);
                break;
            case TN_POW: {
                double a0 = DV(T, nd->a, 0);
                if (k == 0) { r = pow(a0, nd->cd); break; }
                if (a0 == 0.0 && nd->cd > 0.0) {
                    /* identically vanishing base: zero series (as above) */
                    bool allz = true;
                    for (int j = 1; j <= k && allz; j++)
                        if (DV(T, nd->a, j) != 0.0) allz = false;
                    if (allz) break;   /* r stays 0 */
                }
                if (a0 == 0.0) return false;
                for (int j = 0; j < k; j++)
                    r += (nd->cd * (k - j) - j) * DV(T, nd->a, k - j) * DV(T, n, j);
                r /= k * a0;
                break;
            }
            case TN_EXP:
                if (k == 0) { r = exp(DV(T, nd->a, 0)); break; }
                for (int j = 1; j <= k; j++) r += j * DV(T, nd->a, j) * DV(T, n, k - j);
                r /= k;
                break;
            case TN_LOG: {
                double a0 = DV(T, nd->a, 0);
                if (k == 0) { if (!(a0 > 0.0)) return false; r = log(a0); break; }
                for (int j = 1; j < k; j++) r += j * DV(T, n, j) * DV(T, nd->a, k - j);
                r = (DV(T, nd->a, k) - r / k) / a0;
                break;
            }
            case TN_SIN: case TN_SINH: {
                int c = nd->pair;
                bool hyp = (nd->kind == TN_SINH);
                if (k == 0) {
                    double a0 = DV(T, nd->a, 0);
                    r = hyp ? sinh(a0) : sin(a0);
                    DV(T, c, 0) = hyp ? cosh(a0) : cos(a0);
                    break;
                }
                double s2 = 0.0;
                for (int j = 1; j <= k; j++) {
                    double ja = j * DV(T, nd->a, j);
                    r  += ja * DV(T, c, k - j);
                    s2 += ja * DV(T, n, k - j);
                }
                r /= k;
                DV(T, c, k) = hyp ? s2 / k : -s2 / k;
                break;
            }
            case TN_ABS: {
                double ak = DV(T, nd->a, k);
                if (k == 0) nd->sgn = (ak > 0.0) - (ak < 0.0);
                else if (nd->sgn == 0 && ak != 0.0) nd->sgn = (ak > 0.0) ? 1 : -1;
                r = nd->sgn * ak;
                break;
            }
            case TN_ATAN: {
                if (k == 0) { r = atan(DV(T, nd->a, 0)); break; }
                r = k * DV(T, nd->a, k);
                for (int j = 1; j < k; j++) r -= j * DV(T, n, j) * DV(T, nd->b, k - j);
                r /= k * DV(T, nd->b, 0);
                break;
            }
            default:
                return false;
            }
            if (!isfinite(r)) return false;
            DV(T, n, k) = r;
        }
        for (size_t i = 0; i < d; i++) {
            double v = DV(T, T->out[i], k) / (k + 1);
            if (!isfinite(v)) return false;
            C[(size_t)(k + 1) * d + i] = v;
        }
    }
    return true;
}

#endif /* USE_MPFR */
