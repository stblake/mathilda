/*
 * nc_phase.h -- linear-phase splitting shared by NSum and NIntegrate (see
 * nc_phase.c).
 */
#ifndef MATHILDA_NC_PHASE_H
#define MATHILDA_NC_PHASE_H

#include <stdbool.h>
#include "expr.h"
#include "arithmetic.h"   /* is_complex */

/* e free of the symbol var (interned name)? */
bool nc_free_of(const Expr* e, const char* var);
/* A var-dependent phase anywhere in e: b^ex (b free of var, ex not), a
 * circular function of var, or BesselJ/BesselY of a var-dependent argument. */
bool nc_has_phase(const Expr* e, const char* var);
/* Split one evaluated product term into F(var) e^(lam var), F phase-free:
 * returns F and sets *lam (both new), or NULL. */
Expr* nc_strip_phase(const Expr* term, const char* var, Expr** lam);
/* Copy of e with BesselJ/BesselY of an argument c var + d (c > 0, Re > 0 at
 * var = n_check) written through the scaled-Hankel placeholders; *bad is set
 * (NULL returned) for a Bessel factor of another form. */
Expr* nc_bessel_to_hankel(const Expr* e, const char* var, double n_check, bool* bad);
/* Copy of e with the placeholders replaced by their BesselK form. */
Expr* nc_hankel_expand(const Expr* e);
/* Does e contain a scaled-Hankel placeholder? */
bool nc_has_hankel(const Expr* e);

#endif /* MATHILDA_NC_PHASE_H */
