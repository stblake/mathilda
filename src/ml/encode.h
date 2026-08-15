/* encode.h -- the categorical label vocabulary.
 *
 * The piece every earlier family avoided. Clustering, dimensionality reduction,
 * prediction and distributions all took numeric features and produced numeric answers;
 * Classify is where a CLASS arrives, and a class is a string, a symbol, or anything else
 * a user cares to name. So there has to be one place that turns arbitrary expressions
 * into the dense indices an algorithm can count with, and turns them back again.
 *
 * Designed once and put in its own module for the same reason the model representation
 * was: a later ContingencyTable, a categorical FEATURE encoder, and any second classifier
 * all need exactly this, and three copies of it would be three chances to order the
 * classes differently.
 *
 * ORDER IS FIRST APPEARANCE in the training data, not sorted. Sorting would need a total
 * order on arbitrary expressions -- which expr_compare does provide -- but first
 * appearance is the more useful contract: it is stable under adding classes later, it
 * matches how the rest of this codebase numbers clusters (see fc_assign_from_uf), and it
 * makes the vocabulary readable next to the data that produced it. What matters is only
 * that it is DETERMINISTIC, and it is.
 */
#ifndef ML_ENCODE_H
#define ML_ENCODE_H

#include <stddef.h>
#include <stdbool.h>
#include "expr.h"

/* A label vocabulary: the distinct classes, in first-appearance order.
 *
 * `label[i]` is OWNED (a copy), because the training expression it came from may be freed
 * while the fitted model outlives it. */
typedef struct {
    size_t  count;
    Expr**  label;
} MlLabels;

/* Collect the distinct values of `vals` (n entries) into `out`, in first-appearance
 * order. Comparison is expr_eq, so two labels are the same class exactly when they are
 * structurally identical -- "a" and "a" are one class, "a" and a are two, which is the
 * same distinction the pattern matcher makes everywhere else.
 *
 * `index[i]` (n entries, may be NULL) receives the class index of vals[i]. */
bool ml_labels_build(Expr** vals, size_t n, MlLabels* out, size_t* index);

/* Index of `v` in the vocabulary, or SIZE_MAX if absent -- which is a real answer for a
 * classifier asked about a class it never saw, not an error. */
size_t ml_labels_index(const MlLabels* l, const Expr* v);

void ml_labels_free(MlLabels* l);

/* The vocabulary as a plain List, for storing in a model payload. */
Expr* ml_labels_to_list(const MlLabels* l);

/* Read a vocabulary back out of such a List. Borrows nothing: the labels are copied. */
bool ml_labels_from_list(Expr* list, MlLabels* out);

#endif /* ML_ENCODE_H */
