/*
 * hrep_internal.h -- shared internals of the formula-representation engine
 * (hrep_core.c, hrep_certs.c). Not a public interface.
 */
#ifndef MATHILDA_HREP_INTERNAL_H
#define MATHILDA_HREP_INTERNAL_H

#include "expr.h"
#include <stdbool.h>
#include <stdint.h>

/* Largest ODE order handled (d = max(p, q+1)) and the most derivative-formula
 * factors one solve may accumulate. */
#define HR_MAXD    8
#define HR_MAXPOST 64

/* A parameter as symbolic part + rational part: value = sym + n/d, sym NULL
 * when the parameter is rational. Two parameters are in one residue class
 * when their sym parts agree and their rational parts differ by an integer. */
typedef struct { Expr* sym; int64_t n, d; } HrPar;

bool hr_par_split(const Expr* e, HrPar* out);
void hr_par_free(HrPar* p);

/* A lookup certificate: F(a0; b0; w) = f0, with w the symbol Hrep`w.
 * A family certificate has parameters linear in Hrep`al and stands for a
 * whole lattice of origins; `cond` (an expression in Hrep`al, or NULL) is a
 * validity condition, and an instance where it evaluates to False is skipped.
 * Entries are parsed on use. */
typedef struct {
    int p, q;
    const char* a[4];
    const char* b[4];
    const char* f0;
    const char* cond;
    int family;   /* number of family parameters: 0 static, 1 (Hrep`al), 2 (+ Hrep`be) */
} HrCert;

const HrCert* hr_certificates(int* n);

/* e /. sym -> val, evaluated; borrows. */
Expr* hr_subst(const Expr* e, const char* sym, const Expr* val);

/* Present an answer: consumes r, returns the tidied (owned) expression. */
Expr* hr_tidy(Expr* r);

#endif /* MATHILDA_HREP_INTERNAL_H */
