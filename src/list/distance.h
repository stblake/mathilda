#ifndef DISTANCE_H
#define DISTANCE_H

#include "expr.h"

Expr* builtin_euclidean_distance(Expr* res);
Expr* builtin_squared_euclidean_distance(Expr* res);
Expr* builtin_manhattan_distance(Expr* res);
Expr* builtin_cosine_distance(Expr* res);
Expr* builtin_edit_distance(Expr* res);
Expr* builtin_hamming_distance(Expr* res);

/* Exact pairwise distances, shared with FindClusters. NULL when the two points
 * have incompatible shapes. */
Expr* distance_squared_euclidean(Expr* u, Expr* v);
Expr* distance_manhattan(Expr* u, Expr* v);
Expr* distance_cosine(Expr* u, Expr* v);
Expr* distance_edit(Expr* u, Expr* v);
Expr* distance_hamming(Expr* u, Expr* v);

#endif /* DISTANCE_H */
