/* fit.h — linear least-squares regression (Fit) and DesignMatrix.
 *
 * Registers the `Fit` and `DesignMatrix` builtins.  See fit.c for the
 * full description of the supported call forms, data shapes, options
 * (WorkingPrecision, FitRegularization, NormFunction), and solvers.
 */
#ifndef MATHILDA_FIT_H
#define MATHILDA_FIT_H

void fit_init(void);

#endif /* MATHILDA_FIT_H */
