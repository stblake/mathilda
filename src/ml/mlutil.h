/* mlutil.h -- shared plumbing for the src/ml builtins.
 *
 * Promoted here when the THIRD copy appeared, which was the stated trigger: pca.c had
 * ml_list_of_reals / ml_list_matrix / ml_read_data as statics and predict.c had grown
 * its own ml_reals_list. Two copies is a coincidence; three is a shared need.
 *
 * Everything here is about the boundary between an Expr and a machine buffer, which is
 * the one thing every ml builtin does and none of them should each solve.
 */
#ifndef ML_UTIL_H
#define ML_UTIL_H

#include <stddef.h>
#include <stdbool.h>
#include "expr.h"

/* Build a plain List of machine reals, and a plain List of such Lists.
 *
 * NOT na_build_vector / na_build_matrix, deliberately: those return a VISIBLE NDArray,
 * whose head is NDArray, so `result === {{...}}` compares False against the literal a
 * user would write -- while Inverse, Dot and LinearSolve all compare True. Building
 * Lists keeps the ml surface consistent with the rest of the system and lets the
 * evaluator's own packing gate decide whether the result is held as a buffer. */
Expr* ml_list_of_reals(const double* v, size_t n);
Expr* ml_list_matrix(const double* x, size_t n, size_t dim);

/* Read an Expr as an n x dim machine matrix.
 *
 * A flat list is n observations of ONE variable, not one observation of n -- rows are
 * examples and columns are features throughout src/ml. `was_vector` reports which
 * shape arrived, for the callers that must treat them differently (PrincipalComponents
 * declines a single variable; Standardize does not). */
bool ml_read_data(Expr* e, size_t* n, size_t* dim, double** buf, bool* was_vector);

/* Squared Euclidean distance between two dim-vectors.
 *
 * A local four-liner rather than a shared metric layer, and that is a deliberate
 * scoping call. find_clusters.c has a real machine-precision metric layer
 * (fc_dist_pos / fc_dist_to_point) that honours a DistanceFunction, but it is private
 * and takes that file's FcData. Extracting it means introducing a metric enum in
 * src/ml and rewriting a 2600-line file whose 22 one-dimensional answers are pinned --
 * a genuine risk to a finished family, bought for a function this size.
 *
 * The trigger for doing that extraction properly is a src/ml consumer that needs
 * METRIC SELECTION, not one that needs Euclidean: the moment Predict or Classify gains
 * a DistanceFunction option, the private layer has its second real consumer and the
 * refactor pays for itself. Until then this is the honest size of the need. */
double ml_sqdist(const double* a, const double* b, size_t dim);

/* Lower Cholesky factor of a dim x dim row-major symmetric matrix into `l`. Returns
 * false if the matrix is not positive definite -- which is information, not merely an
 * error: it says the covariance is singular, i.e. fewer points than dimensions or
 * collinear ones.
 *
 * Written in-house rather than routed through mat_lapack_dpotrf because `dim` is a
 * feature count -- single digits in practice -- and the callers run this once per
 * component per iteration, where the cost is elsewhere entirely. */
bool ml_chol(const double* a, double* l, size_t dim);

/* Squared Mahalanobis distance (x - mu)' S^-1 (x - mu), via forward substitution
 * against the Cholesky factor `l`. `y` is scratch of length dim. */
double ml_mahalanobis(const double* l, const double* mu, const double* x,
                      size_t dim, double* y);

/* A real feature value: a number (na_read_scalar) or an exact numeric constant
 * (1/2, Sqrt[2], Pi) at its machine value. False for anything else -- g[9] is
 * not the number 9. */
bool ml_read_real(const Expr* e, double* out);

/* Batch application of a dim-feature model, model[{x1, x2, ...}]: when args[0]
 * is a List of examples -- each a List of dim features, or (dim == 1) at least
 * two scalars -- returns the List of apply(head, {example, args[1..]}) answers
 * (NULL if any example declines). Sets *batch to whether args[0] was a batch;
 * when it is not, returns NULL and the caller reads a single example. */
Expr* ml_apply_batch(Expr* head, Expr** args, size_t argc, size_t dim,
                     Expr* (*apply)(Expr*, Expr**, size_t), bool* batch);

#endif /* ML_UTIL_H */
