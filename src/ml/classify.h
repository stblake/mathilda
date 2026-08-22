/* classify.h -- Classify. See classify.c for the design notes. */
#ifndef ML_CLASSIFY_H
#define ML_CLASSIFY_H
#include "expr.h"
Expr* ml_classifier_apply(Expr* head, Expr** args, size_t argc);
void ml_classify_init(void);
#endif /* ML_CLASSIFY_H */
