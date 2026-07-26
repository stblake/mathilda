/* Mathilda — numeric compiler for the NDSolve reduced RHS (see header). */
#include "ndsolve_compile.h"
#include "../expr.h"
#include "../arithmetic.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.7182818284590452354
#endif

/* ------------------------------------------------------------------ *
 *  Bytecode                                                           *
 * ------------------------------------------------------------------ */
enum {
    OP_CONST, OP_VAR, OP_TVAR,
    OP_ADD, OP_SUB, OP_MUL, OP_DIV, OP_NEG, OP_INV,
    OP_POW, OP_POWI, OP_SQRT, OP_EXP, OP_LOG,
    OP_SIN, OP_COS, OP_TAN, OP_COT, OP_SEC, OP_CSC,
    OP_SINH, OP_COSH, OP_TANH, OP_ASIN, OP_ACOS, OP_ATAN,
    OP_ABS, OP_SIGN, OP_ERF, OP_ERFC, OP_FLOOR, OP_CEIL, OP_ROUND,
    OP_MAX, OP_MIN, OP_ATAN2
};

typedef struct { uint8_t op; int32_t i; double d; } NdInstr;

/* Interned-name -> reduced-state-index map (open addressing on the interned
 * pointer).  Built from P->ysym so the compiler recognizes the state symbols of
 * ANY front-end — MoL PDEs (NDSolve`w<k>), ODEs (NDSolve`y<k>), etc. — instead
 * of a hardcoded naming convention. */
typedef struct { const char** key; int* val; size_t cap; } NameMap;

static bool namemap_init(NameMap* m, const NdProblem* P) {
    size_t cap = 8;
    while (cap < P->d * 2) cap <<= 1;
    m->cap = cap;
    m->key = calloc(cap, sizeof(const char*));
    m->val = malloc(cap * sizeof(int));
    if (!m->key || !m->val) { free(m->key); free(m->val); m->key = NULL; return false; }
    for (size_t k = 0; k < P->d; k++) {
        if (!P->ysym[k] || P->ysym[k]->type != EXPR_SYMBOL) continue;
        const char* nm = P->ysym[k]->data.symbol.name;
        size_t h = ((uintptr_t)nm >> 4) & (cap - 1);
        while (m->key[h]) h = (h + 1) & (cap - 1);
        m->key[h] = nm; m->val[h] = (int)k;
    }
    return true;
}
static int namemap_get(const NameMap* m, const char* nm) {
    size_t h = ((uintptr_t)nm >> 4) & (m->cap - 1);
    while (m->key[h]) {
        if (m->key[h] == nm) return m->val[h];
        h = (h + 1) & (m->cap - 1);
    }
    return -1;
}
static void namemap_free(NameMap* m) { free(m->key); free(m->val); }

typedef struct {
    NdInstr* code;
    size_t   n;          /* instruction count                 */
    int      stackmax;   /* max value-stack depth             */
    int*     dep;        /* sorted state indices read (owned) */
    int      ndep;
} NdProg;

struct NdCompiled {
    size_t   d;
    NdProg*  prog;       /* [d] */
    double*  stack;      /* [global stackmax] scratch         */
    /* CPR coloring for the finite-difference Jacobian */
    int*     color;      /* [d] column -> color (>=0)         */
    int      ncolor;
    /* Jacobian scratch (allocated lazily) */
    double  *Yp, *Ym, *fp, *fm, *hj;
};

/* ------------------------------------------------------------------ *
 *  Emitter                                                            *
 * ------------------------------------------------------------------ */
typedef struct {
    NdInstr* code;
    size_t   n, cap;
    int      depth, maxdepth;
    bool     ok;
    /* per-component variable dependency set */
    bool*    seen;       /* [d] */
    size_t   d;
    const char* tvar;
    const NameMap* map;  /* state-symbol name -> index */
} Emitter;

static void em_ins(Emitter* E, uint8_t op, int32_t iv, double dv, int delta) {
    if (!E->ok) return;
    if (E->n == E->cap) {
        size_t nc = E->cap ? E->cap * 2 : 32;
        NdInstr* nb = realloc(E->code, nc * sizeof(NdInstr));
        if (!nb) { E->ok = false; return; }
        E->code = nb; E->cap = nc;
    }
    E->code[E->n].op = op; E->code[E->n].i = iv; E->code[E->n].d = dv;
    E->n++;
    E->depth += delta;
    if (E->depth > E->maxdepth) E->maxdepth = E->depth;
}
#define PUSH(E, op, iv, dv) em_ins(E, op, iv, dv, +1)
#define UNARY(E, op)        em_ins(E, op, 0, 0.0, 0)
#define BINARY(E, op)       em_ins(E, op, 0, 0.0, -1)

static bool named_const(const char* nm, double* out) {
    if (strcmp(nm, "Pi") == 0)         { *out = M_PI; return true; }
    if (strcmp(nm, "E") == 0)          { *out = M_E;  return true; }
    if (strcmp(nm, "EulerGamma") == 0) { *out = 0.57721566490153286061; return true; }
    if (strcmp(nm, "Degree") == 0)     { *out = M_PI / 180.0; return true; }
    if (strcmp(nm, "GoldenRatio") == 0){ *out = 1.61803398874989484820; return true; }
    if (strcmp(nm, "Catalan") == 0)    { *out = 0.91596559417721901505; return true; }
    return false;
}

/* Map a unary-function head name to an opcode; returns -1 if not unary. */
static int unary_op(const char* h) {
    if (strcmp(h, "Sqrt") == 0)    return OP_SQRT;
    if (strcmp(h, "Exp") == 0)     return OP_EXP;
    if (strcmp(h, "Sin") == 0)     return OP_SIN;
    if (strcmp(h, "Cos") == 0)     return OP_COS;
    if (strcmp(h, "Tan") == 0)     return OP_TAN;
    if (strcmp(h, "Cot") == 0)     return OP_COT;
    if (strcmp(h, "Sec") == 0)     return OP_SEC;
    if (strcmp(h, "Csc") == 0)     return OP_CSC;
    if (strcmp(h, "Sinh") == 0)    return OP_SINH;
    if (strcmp(h, "Cosh") == 0)    return OP_COSH;
    if (strcmp(h, "Tanh") == 0)    return OP_TANH;
    if (strcmp(h, "ArcSin") == 0)  return OP_ASIN;
    if (strcmp(h, "ArcCos") == 0)  return OP_ACOS;
    if (strcmp(h, "Abs") == 0)     return OP_ABS;
    if (strcmp(h, "Sign") == 0)    return OP_SIGN;
    if (strcmp(h, "Erf") == 0)     return OP_ERF;
    if (strcmp(h, "Erfc") == 0)    return OP_ERFC;
    if (strcmp(h, "Floor") == 0)   return OP_FLOOR;
    if (strcmp(h, "Ceiling") == 0) return OP_CEIL;
    if (strcmp(h, "Round") == 0)   return OP_ROUND;
    return -1;
}

static bool emit(Emitter* E, const Expr* e);

/* Emit an n-ary reduction (Plus/Times/Max/Min): fold left with `binop`. */
static bool emit_nary(Emitter* E, const Expr* e, uint8_t binop, double ident) {
    size_t n = e->data.function.arg_count;
    if (n == 0) { PUSH(E, OP_CONST, 0, ident); return E->ok; }
    if (!emit(E, e->data.function.args[0])) return false;
    for (size_t i = 1; i < n; i++) {
        if (!emit(E, e->data.function.args[i])) return false;
        BINARY(E, binop);
    }
    return E->ok;
}

static bool emit(Emitter* E, const Expr* e) {
    if (!e || !E->ok) { E->ok = false; return false; }

    /* numeric leaf (Integer/Real/BigInt/Rational/MPFR) */
    double v;
    if (nd_to_double(e, &v)) { PUSH(E, OP_CONST, 0, v); return E->ok; }

    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        int k = namemap_get(E->map, nm);
        if (k >= 0 && (size_t)k < E->d) {
            E->seen[k] = true;
            PUSH(E, OP_VAR, k, 0.0);
            return E->ok;
        }
        if (E->tvar && strcmp(nm, E->tvar) == 0) { PUSH(E, OP_TVAR, 0, 0.0); return E->ok; }
        double c;
        if (named_const(nm, &c)) { PUSH(E, OP_CONST, 0, c); return E->ok; }
        E->ok = false; return false;                 /* free/unknown symbol */
    }

    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) {
        E->ok = false; return false;
    }
    const char* h = e->data.function.head->data.symbol.name;
    Expr** A = e->data.function.args;
    size_t na = e->data.function.arg_count;

    if (strcmp(h, "Plus") == 0)  return emit_nary(E, e, OP_ADD, 0.0);
    if (strcmp(h, "Times") == 0) return emit_nary(E, e, OP_MUL, 1.0);
    if (strcmp(h, "Max") == 0)   return emit_nary(E, e, OP_MAX, -HUGE_VAL);
    if (strcmp(h, "Min") == 0)   return emit_nary(E, e, OP_MIN,  HUGE_VAL);

    if (strcmp(h, "Power") == 0 && na == 2) {
        const Expr* base = A[0]; const Expr* ex = A[1];
        if (ex->type == EXPR_INTEGER) {
            long n = ex->data.integer;
            if (!emit(E, base)) return false;
            if (n == -1) { UNARY(E, OP_INV); }
            else         { em_ins(E, OP_POWI, (int32_t)n, 0.0, 0); }
            return E->ok;
        }
        /* Rational ±1/2 -> sqrt fast paths */
        int64_t rn, rd;
        if (is_rational(ex, &rn, &rd) && (rd == 2) && (rn == 1 || rn == -1)) {
            if (!emit(E, base)) return false;
            UNARY(E, OP_SQRT);
            if (rn == -1) UNARY(E, OP_INV);
            return E->ok;
        }
        if (!emit(E, base)) return false;
        if (!emit(E, ex)) return false;
        BINARY(E, OP_POW);
        return E->ok;
    }

    if (strcmp(h, "Subtract") == 0 && na == 2) {
        if (!emit(E, A[0]) || !emit(E, A[1])) return false;
        BINARY(E, OP_SUB); return E->ok;
    }
    if (strcmp(h, "Divide") == 0 && na == 2) {
        if (!emit(E, A[0]) || !emit(E, A[1])) return false;
        BINARY(E, OP_DIV); return E->ok;
    }
    if (strcmp(h, "Minus") == 0 && na == 1) {
        if (!emit(E, A[0])) return false;
        UNARY(E, OP_NEG); return E->ok;
    }
    if (strcmp(h, "Log") == 0) {
        if (na == 1) {
            if (!emit(E, A[0])) return false;
            UNARY(E, OP_LOG);
            return E->ok;
        }
        if (na == 2) {                     /* Log[b, x] = log(x)/log(b) */
            if (!emit(E, A[1])) return false;
            UNARY(E, OP_LOG);
            if (!emit(E, A[0])) return false;
            UNARY(E, OP_LOG);
            BINARY(E, OP_DIV);
            return E->ok;
        }
        E->ok = false; return false;
    }
    if (strcmp(h, "ArcTan") == 0) {
        if (na == 1) {
            if (!emit(E, A[0])) return false;
            UNARY(E, OP_ATAN);
            return E->ok;
        }
        if (na == 2) {                     /* ArcTan[x, y] = atan2(y, x) */
            if (!emit(E, A[0]) || !emit(E, A[1])) return false;
            BINARY(E, OP_ATAN2); return E->ok;
        }
        E->ok = false; return false;
    }

    if (na == 1) {
        int u = unary_op(h);
        if (u >= 0) {
            if (!emit(E, A[0])) return false;
            UNARY(E, (uint8_t)u);
            return E->ok;
        }
    }

    E->ok = false; return false;               /* unsupported head -> bail */
}

/* ------------------------------------------------------------------ *
 *  Compile                                                            *
 * ------------------------------------------------------------------ */
NdCompiled* nd_compile_rhs(const NdProblem* P) {
    size_t d = P->d;
    if (d == 0) return NULL;
    NdCompiled* C = calloc(1, sizeof(*C));
    if (!C) return NULL;
    C->d = d;
    C->prog = calloc(d, sizeof(NdProg));
    bool* seen = malloc(sizeof(bool) * d);
    NameMap map;
    if (!C->prog || !seen || !P->ysym || !namemap_init(&map, P)) {
        free(seen); nd_compiled_free(C); return NULL;
    }

    int gstackmax = 1;
    bool ok = true;
    for (size_t i = 0; i < d && ok; i++) {
        Emitter E; memset(&E, 0, sizeof(E));
        E.ok = true; E.seen = seen; E.d = d; E.tvar = P->tvar; E.map = &map;
        memset(seen, 0, sizeof(bool) * d);
        if (!emit(&E, P->f[i]) || !E.ok) { free(E.code); ok = false; break; }
        C->prog[i].code = E.code;
        C->prog[i].n = E.n;
        C->prog[i].stackmax = E.maxdepth;
        if (E.maxdepth > gstackmax) gstackmax = E.maxdepth;
        int nd = 0;
        for (size_t j = 0; j < d; j++) if (seen[j]) nd++;
        C->prog[i].dep = malloc(sizeof(int) * (nd ? nd : 1));
        C->prog[i].ndep = nd;
        int p = 0;
        for (size_t j = 0; j < d; j++) if (seen[j]) C->prog[i].dep[p++] = (int)j;
    }
    free(seen);
    namemap_free(&map);
    if (!ok) { nd_compiled_free(C); return NULL; }

    C->stack = malloc(sizeof(double) * (size_t)gstackmax);
    if (!C->stack) { nd_compiled_free(C); return NULL; }

    /* CPR column coloring: two columns share a color only if no row reads both.
     * Greedy over columns, using each row's dependency list as the constraint. */
    C->color = malloc(sizeof(int) * d);
    if (!C->color) { nd_compiled_free(C); return NULL; }
    for (size_t j = 0; j < d; j++) C->color[j] = -1;
    /* row lists indexed by column: which rows read column j (transpose of dep) */
    /* For each column j, forbid colors already used by any column that co-occurs
     * with j in some row.  Build co-occurrence lazily per column. */
    int maxcolor = -1;
    bool* used = calloc(d, sizeof(bool));   /* color-used scratch, size <= d */
    /* For quick "columns sharing a row with j", precompute per-column the rows
     * that touch it, then walk those rows' deps. */
    int* col_rows_off = calloc(d + 1, sizeof(int));
    if (!used || !col_rows_off) { free(used); free(col_rows_off); nd_compiled_free(C); return NULL; }
    for (size_t i = 0; i < d; i++)
        for (int t = 0; t < C->prog[i].ndep; t++)
            col_rows_off[C->prog[i].dep[t] + 1]++;
    for (size_t j = 0; j < d; j++) col_rows_off[j + 1] += col_rows_off[j];
    int* col_rows = malloc(sizeof(int) * (size_t)(col_rows_off[d] ? col_rows_off[d] : 1));
    int* fillpos = malloc(sizeof(int) * d);
    if (!col_rows || !fillpos) { free(used); free(col_rows_off); free(col_rows); free(fillpos); nd_compiled_free(C); return NULL; }
    for (size_t j = 0; j < d; j++) fillpos[j] = col_rows_off[j];
    for (size_t i = 0; i < d; i++)
        for (int t = 0; t < C->prog[i].ndep; t++) {
            int j = C->prog[i].dep[t];
            col_rows[fillpos[j]++] = (int)i;
        }
    for (size_t j = 0; j < d; j++) {
        /* mark colors used by columns co-occurring with j */
        for (int rp = col_rows_off[j]; rp < col_rows_off[j + 1]; rp++) {
            int i = col_rows[rp];
            for (int t = 0; t < C->prog[i].ndep; t++) {
                int j2 = C->prog[i].dep[t];
                if (j2 != (int)j && C->color[j2] >= 0) used[C->color[j2]] = true;
            }
        }
        int c = 0; while (c <= maxcolor && used[c]) c++;
        C->color[j] = c;
        if (c > maxcolor) maxcolor = c;
        /* clear the scratch we touched */
        for (int rp = col_rows_off[j]; rp < col_rows_off[j + 1]; rp++) {
            int i = col_rows[rp];
            for (int t = 0; t < C->prog[i].ndep; t++) {
                int j2 = C->prog[i].dep[t];
                if (j2 != (int)j && C->color[j2] >= 0) used[C->color[j2]] = false;
            }
        }
    }
    C->ncolor = maxcolor + 1;
    free(used); free(col_rows_off); free(col_rows); free(fillpos);
    if (C->ncolor < 1) C->ncolor = 1;
    return C;
}

/* ------------------------------------------------------------------ *
 *  VM                                                                 *
 * ------------------------------------------------------------------ */
static double ipow(double b, long n) {
    if (n < 0) { b = 1.0 / b; n = -n; }
    double r = 1.0;
    while (n) { if (n & 1) r *= b; b *= b; n >>= 1; }
    return r;
}

static double run_prog(const NdProg* pr, const double* Y, double t, double* S) {
    const NdInstr* c = pr->code;
    size_t n = pr->n, sp = 0;
    for (size_t k = 0; k < n; k++) {
        double x;
        switch (c[k].op) {
            case OP_CONST: S[sp++] = c[k].d; break;
            case OP_VAR:   S[sp++] = Y[c[k].i]; break;
            case OP_TVAR:  S[sp++] = t; break;
            case OP_ADD:   S[sp-2] += S[sp-1]; sp--; break;
            case OP_SUB:   S[sp-2] -= S[sp-1]; sp--; break;
            case OP_MUL:   S[sp-2] *= S[sp-1]; sp--; break;
            case OP_DIV:   S[sp-2] /= S[sp-1]; sp--; break;
            case OP_NEG:   S[sp-1] = -S[sp-1]; break;
            case OP_INV:   S[sp-1] = 1.0 / S[sp-1]; break;
            case OP_POW:   S[sp-2] = pow(S[sp-2], S[sp-1]); sp--; break;
            case OP_POWI:  S[sp-1] = ipow(S[sp-1], c[k].i); break;
            case OP_SQRT:  S[sp-1] = sqrt(S[sp-1]); break;
            case OP_EXP:   S[sp-1] = exp(S[sp-1]); break;
            case OP_LOG:   S[sp-1] = log(S[sp-1]); break;
            case OP_SIN:   S[sp-1] = sin(S[sp-1]); break;
            case OP_COS:   S[sp-1] = cos(S[sp-1]); break;
            case OP_TAN:   S[sp-1] = tan(S[sp-1]); break;
            case OP_COT:   S[sp-1] = 1.0 / tan(S[sp-1]); break;
            case OP_SEC:   S[sp-1] = 1.0 / cos(S[sp-1]); break;
            case OP_CSC:   S[sp-1] = 1.0 / sin(S[sp-1]); break;
            case OP_SINH:  S[sp-1] = sinh(S[sp-1]); break;
            case OP_COSH:  S[sp-1] = cosh(S[sp-1]); break;
            case OP_TANH:  S[sp-1] = tanh(S[sp-1]); break;
            case OP_ASIN:  S[sp-1] = asin(S[sp-1]); break;
            case OP_ACOS:  S[sp-1] = acos(S[sp-1]); break;
            case OP_ATAN:  S[sp-1] = atan(S[sp-1]); break;
            case OP_ABS:   S[sp-1] = fabs(S[sp-1]); break;
            case OP_SIGN:  x = S[sp-1]; S[sp-1] = (x > 0) - (x < 0); break;
            case OP_ERF:   S[sp-1] = erf(S[sp-1]); break;
            case OP_ERFC:  S[sp-1] = erfc(S[sp-1]); break;
            case OP_FLOOR: S[sp-1] = floor(S[sp-1]); break;
            case OP_CEIL:  S[sp-1] = ceil(S[sp-1]); break;
            case OP_ROUND: S[sp-1] = round(S[sp-1]); break;
            case OP_MAX:   S[sp-2] = S[sp-2] >= S[sp-1] ? S[sp-2] : S[sp-1]; sp--; break;
            case OP_MIN:   S[sp-2] = S[sp-2] <= S[sp-1] ? S[sp-2] : S[sp-1]; sp--; break;
            case OP_ATAN2: S[sp-2] = atan2(S[sp-1], S[sp-2]); sp--; break;
            default: return NAN;
        }
    }
    return sp ? S[0] : NAN;
}

bool nd_compiled_eval(NdCompiled* C, double t, const double* Y, double* out) {
    for (size_t i = 0; i < C->d; i++) {
        out[i] = run_prog(&C->prog[i], Y, t, C->stack);
        if (!isfinite(out[i])) return false;
    }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Colored finite-difference Jacobian                                 *
 * ------------------------------------------------------------------ */
static bool jac_alloc(NdCompiled* C) {
    if (C->Yp) return true;
    size_t d = C->d;
    C->Yp = malloc(sizeof(double) * d);
    C->Ym = malloc(sizeof(double) * d);
    C->fp = malloc(sizeof(double) * d);
    C->fm = malloc(sizeof(double) * d);
    C->hj = malloc(sizeof(double) * d);
    if (!C->Yp || !C->Ym || !C->fp || !C->fm || !C->hj) return false;
    return true;
}

bool nd_compiled_jacobian(NdCompiled* C, double t, const double* Y, double* Jout) {
    size_t d = C->d;
    if (!jac_alloc(C)) return false;
    memset(Jout, 0, sizeof(double) * d * d);
    for (int col = 0; col < C->ncolor; col++) {
        /* perturb every column of this color simultaneously */
        memcpy(C->Yp, Y, sizeof(double) * d);
        memcpy(C->Ym, Y, sizeof(double) * d);
        for (size_t j = 0; j < d; j++) {
            if (C->color[j] != col) continue;
            double hj = (fabs(Y[j]) + 1.0) * 1.0e-7;
            C->hj[j] = hj;
            C->Yp[j] += hj;
            C->Ym[j] -= hj;
        }
        if (!nd_compiled_eval(C, t, C->Yp, C->fp)) return false;
        if (!nd_compiled_eval(C, t, C->Ym, C->fm)) return false;
        /* extract: within a color, each row depends on at most one such column */
        for (size_t i = 0; i < d; i++) {
            const NdProg* pr = &C->prog[i];
            for (int tt = 0; tt < pr->ndep; tt++) {
                int j = pr->dep[tt];
                if (C->color[j] != col) continue;
                Jout[i * d + (size_t)j] = (C->fp[i] - C->fm[i]) / (2.0 * C->hj[j]);
            }
        }
    }
    return true;
}

int nd_compiled_ncolor(const NdCompiled* C) { return C ? C->ncolor : 0; }

/* ------------------------------------------------------------------ *
 *  Free                                                               *
 * ------------------------------------------------------------------ */
void nd_compiled_free(NdCompiled* C) {
    if (!C) return;
    if (C->prog) {
        for (size_t i = 0; i < C->d; i++) { free(C->prog[i].code); free(C->prog[i].dep); }
        free(C->prog);
    }
    free(C->stack); free(C->color);
    free(C->Yp); free(C->Ym); free(C->fp); free(C->fm); free(C->hj);
    free(C);
}
