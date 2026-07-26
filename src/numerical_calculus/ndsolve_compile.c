/* Mathilda — NDSolve reduced-RHS compiler (see header).
 *
 * This is now a thin client of the general Compile[] engine (src/compile): each
 * reduced component f_i(t, Y) is compiled to a CompiledProgram over the state
 * symbols and the time variable, and evaluated via the shared typed VM — no
 * bespoke bytecode here.  What remains NDSolve-specific is the Curtis–Powell–Reid
 * colored finite-difference Jacobian, built on the per-component variable
 * dependencies the engine exposes (compiled_arg_deps).
 *
 * The public API (nd_compile_rhs / nd_compiled_eval / nd_compiled_jacobian /
 * nd_compiled_ncolor / nd_compiled_free) is unchanged, so ndsolve_common.c and
 * ndsolve_mol.c are untouched. */
#include "ndsolve_compile.h"
#include "../compile/compile.h"
#include "../expr.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct NdCompiled {
    size_t d;
    CompiledProgram** prog;   /* [d] one program per reduced component      */
    double*  args;            /* [d+1] scratch: [Y_0..Y_{d-1}, t]           */
    int**    dep;             /* [d] dep[i] = sorted state indices row i reads */
    int*     ndep;            /* [d]                                        */
    /* CPR coloring for the finite-difference Jacobian */
    int*     color;           /* [d] column -> color (>=0)                  */
    int      ncolor;
    /* Jacobian scratch (allocated lazily) */
    double  *Yp, *Ym, *fp, *fm, *hj;
};

/* ------------------------------------------------------------------ *
 *  Compile: one program per component, then CPR coloring              *
 * ------------------------------------------------------------------ */
NdCompiled* nd_compile_rhs(const NdProblem* P) {
    size_t d = P->d;
    if (d == 0 || !P->ysym || !P->f || !P->tvar) return NULL;

    /* argument set shared by every component: the d reduced-state symbols, then
     * the time variable — all real (complex PDEs are realified before compile). */
    const char** names = malloc(sizeof(char*) * (d + 1));
    CompileType* types = malloc(sizeof(CompileType) * (d + 1));
    if (!names || !types) { free(names); free(types); return NULL; }
    for (size_t k = 0; k < d; k++) {
        if (!P->ysym[k] || P->ysym[k]->type != EXPR_SYMBOL) { free(names); free(types); return NULL; }
        names[k] = P->ysym[k]->data.symbol.name;
        types[k] = CT_REAL;
    }
    names[d] = P->tvar; types[d] = CT_REAL;

    NdCompiled* C = calloc(1, sizeof(*C));
    if (!C) { free(names); free(types); return NULL; }
    C->d = d;
    C->prog = calloc(d, sizeof(CompiledProgram*));
    C->dep = calloc(d, sizeof(int*));
    C->ndep = calloc(d, sizeof(int));
    C->args = malloc(sizeof(double) * (d + 1));
    if (!C->prog || !C->dep || !C->ndep || !C->args) { free(names); free(types); nd_compiled_free(C); return NULL; }

    int* depbuf = malloc(sizeof(int) * (d + 1));
    bool ok = depbuf != NULL;
    for (size_t i = 0; i < d && ok; i++) {
        C->prog[i] = compile_expr(P->f[i], names, types, d + 1);
        /* compiled_eval_real needs an all-real signature; a non-real result
         * (an oddity for a real PDE) is bailed to the symbolic sampler. */
        if (!C->prog[i] || compiled_result_type(C->prog[i]) != CT_REAL) { ok = false; break; }
        /* per-row state dependencies (drop the time column d) for the Jacobian */
        size_t nraw = compiled_arg_deps(C->prog[i], depbuf, d + 1);
        int nd = 0;
        for (size_t r = 0; r < nraw; r++) if (depbuf[r] < (int)d) nd++;
        C->dep[i] = malloc(sizeof(int) * (nd ? nd : 1));
        if (!C->dep[i]) { ok = false; break; }
        int p = 0;
        for (size_t r = 0; r < nraw; r++) if (depbuf[r] < (int)d) C->dep[i][p++] = depbuf[r];
        C->ndep[i] = nd;
    }
    free(depbuf); free(names); free(types);
    if (!ok) { nd_compiled_free(C); return NULL; }

    /* CPR column coloring: two columns share a color only if no row reads both. */
    C->color = malloc(sizeof(int) * d);
    bool* used = calloc(d, sizeof(bool));
    int* col_rows_off = calloc(d + 1, sizeof(int));
    if (!C->color || !used || !col_rows_off) { free(used); free(col_rows_off); nd_compiled_free(C); return NULL; }
    for (size_t j = 0; j < d; j++) C->color[j] = -1;
    for (size_t i = 0; i < d; i++)
        for (int t = 0; t < C->ndep[i]; t++) col_rows_off[C->dep[i][t] + 1]++;
    for (size_t j = 0; j < d; j++) col_rows_off[j + 1] += col_rows_off[j];
    int* col_rows = malloc(sizeof(int) * (size_t)(col_rows_off[d] ? col_rows_off[d] : 1));
    int* fillpos = malloc(sizeof(int) * d);
    if (!col_rows || !fillpos) { free(used); free(col_rows_off); free(col_rows); free(fillpos); nd_compiled_free(C); return NULL; }
    for (size_t j = 0; j < d; j++) fillpos[j] = col_rows_off[j];
    for (size_t i = 0; i < d; i++)
        for (int t = 0; t < C->ndep[i]; t++) { int j = C->dep[i][t]; col_rows[fillpos[j]++] = (int)i; }
    int maxcolor = -1;
    for (size_t j = 0; j < d; j++) {
        for (int rp = col_rows_off[j]; rp < col_rows_off[j + 1]; rp++) {
            int i = col_rows[rp];
            for (int t = 0; t < C->ndep[i]; t++) { int j2 = C->dep[i][t]; if (j2 != (int)j && C->color[j2] >= 0) used[C->color[j2]] = true; }
        }
        int c = 0; while (c <= maxcolor && used[c]) c++;
        C->color[j] = c; if (c > maxcolor) maxcolor = c;
        for (int rp = col_rows_off[j]; rp < col_rows_off[j + 1]; rp++) {
            int i = col_rows[rp];
            for (int t = 0; t < C->ndep[i]; t++) { int j2 = C->dep[i][t]; if (j2 != (int)j && C->color[j2] >= 0) used[C->color[j2]] = false; }
        }
    }
    C->ncolor = maxcolor + 1;
    free(used); free(col_rows_off); free(col_rows); free(fillpos);
    if (C->ncolor < 1) C->ncolor = 1;
    return C;
}

/* ------------------------------------------------------------------ *
 *  Evaluate + colored FD Jacobian                                     *
 * ------------------------------------------------------------------ */
bool nd_compiled_eval(NdCompiled* C, double t, const double* Y, double* out) {
    memcpy(C->args, Y, sizeof(double) * C->d);
    C->args[C->d] = t;
    /* batch: load the shared state+t once, run every component on one frame */
    return compiled_eval_real_batch((const CompiledProgram* const*)C->prog, C->d,
                                    C->args, C->d + 1, out);
}

static bool jac_alloc(NdCompiled* C) {
    if (C->Yp) return true;
    size_t d = C->d;
    C->Yp = malloc(sizeof(double) * d); C->Ym = malloc(sizeof(double) * d);
    C->fp = malloc(sizeof(double) * d); C->fm = malloc(sizeof(double) * d);
    C->hj = malloc(sizeof(double) * d);
    return C->Yp && C->Ym && C->fp && C->fm && C->hj;
}

bool nd_compiled_jacobian(NdCompiled* C, double t, const double* Y, double* Jout) {
    size_t d = C->d;
    if (!jac_alloc(C)) return false;
    memset(Jout, 0, sizeof(double) * d * d);
    for (int col = 0; col < C->ncolor; col++) {
        memcpy(C->Yp, Y, sizeof(double) * d);
        memcpy(C->Ym, Y, sizeof(double) * d);
        for (size_t j = 0; j < d; j++) {
            if (C->color[j] != col) continue;
            double hj = (fabs(Y[j]) + 1.0) * 1.0e-7;
            C->hj[j] = hj; C->Yp[j] += hj; C->Ym[j] -= hj;
        }
        if (!nd_compiled_eval(C, t, C->Yp, C->fp)) return false;
        if (!nd_compiled_eval(C, t, C->Ym, C->fm)) return false;
        for (size_t i = 0; i < d; i++)
            for (int tt = 0; tt < C->ndep[i]; tt++) {
                int j = C->dep[i][tt];
                if (C->color[j] != col) continue;
                Jout[i * d + (size_t)j] = (C->fp[i] - C->fm[i]) / (2.0 * C->hj[j]);
            }
    }
    return true;
}

int nd_compiled_ncolor(const NdCompiled* C) { return C ? C->ncolor : 0; }

void nd_compiled_free(NdCompiled* C) {
    if (!C) return;
    if (C->prog) { for (size_t i = 0; i < C->d; i++) compiled_free(C->prog[i]); free(C->prog); }
    if (C->dep)  { for (size_t i = 0; i < C->d; i++) free(C->dep[i]); free(C->dep); }
    free(C->ndep); free(C->args); free(C->color);
    free(C->Yp); free(C->Ym); free(C->fp); free(C->fm); free(C->hj);
    free(C);
}
