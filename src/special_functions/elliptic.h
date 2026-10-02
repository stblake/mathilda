#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include "expr.h"
#include <stdbool.h>

/* The Legendre elliptic integrals, in the PARAMETER convention m = k^2 that the
 * Wolfram Language uses -- NOT the modulus k. Getting this wrong is silent, so
 * every entry point below names its integral explicitly:
 *
 *   EllipticK[m]          = Int_0^(Pi/2) dt / Sqrt(1 - m Sin[t]^2)
 *   EllipticF[phi, m]     = Int_0^phi    dt / Sqrt(1 - m Sin[t]^2)
 *   EllipticE[m]          = Int_0^(Pi/2) Sqrt(1 - m Sin[t]^2) dt
 *   EllipticE[phi, m]     = Int_0^phi    Sqrt(1 - m Sin[t]^2) dt
 *   EllipticPi[n, m]      = Int_0^(Pi/2) dt / ((1 - n Sin[t]^2) Sqrt(1 - m Sin[t]^2))
 *   EllipticPi[n, phi, m] = Int_0^phi    dt / ((1 - n Sin[t]^2) Sqrt(1 - m Sin[t]^2))
 *
 * EllipticE and EllipticPi are arity-overloaded exactly as in Wolfram: one
 * argument is complete for E, two arguments are complete for Pi. */
Expr* builtin_ellipticf(Expr* res);
Expr* builtin_elliptice(Expr* res);
Expr* builtin_ellipticpi(Expr* res);
Expr* builtin_elliptick(Expr* res);
void  elliptic_init(void);

/* Machine-precision real kernels in the shared ND kernel ABI: Carlson symmetric
 * forms in `double`, for the packed/NDArray element-wise paths. Each returns
 * false to DECLINE — outside the real principal domain, at a pole, or wherever
 * the answer is not a real machine number (EllipticPi with n > 1 past the pole
 * is genuinely complex, so it declines here and the FLINT path takes it). */
bool elliptic_machine_k(double m, double* out);
bool elliptic_machine_e_complete(double m, double* out);
bool elliptic_machine_f(double phi, double m, double* out);
bool elliptic_machine_e_inc(double phi, double m, double* out);
/* The third kind, via Carlson R_J. Both DECLINE where p = 1 - n (complete) or
 * p = 1 - n Sin[phi]^2 (incomplete) is <= 0: there the value is genuinely
 * complex (Pi[3/2 | 1/2] = -0.4567 - 2.7207 I), not a real principal value, so
 * a double out parameter cannot carry it and Arb must place the branch. */
bool elliptic_machine_pi(double n, double m, double* out);
bool elliptic_machine_pi_inc(double n, double phi, double m, double* out);

#endif /* ELLIPTIC_H */
