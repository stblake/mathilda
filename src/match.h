#ifndef MATCH_H
#define MATCH_H

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

// Environment to store pattern bindings (e.g., x -> 4)
typedef struct MatchEnv {
    char** symbols;
    Expr** values;
    size_t count;
    size_t capacity;
    bool (*callback)(struct MatchEnv*, void*);
    void* callback_data;
    Expr* top_pattern;   /* the pattern of the innermost match() in progress
                            (borrowed); used to bind leftover names on success */
} MatchEnv;

// Create a new empty match environment
MatchEnv* env_new(void);

// Free the match environment (does not free the Expr* values, they are typically shared or cloned by the user)
void env_free(MatchEnv* env);

// Add a binding to the environment
void env_set(MatchEnv* env, const char* symbol, Expr* value);

// Get a binding from the environment, returns NULL if not found
Expr* env_get(MatchEnv* env, const char* symbol);

// Match an expression against a pattern. Returns true if match succeeds,
// and populates the env with bindings.
bool match(Expr* expr, Expr* pattern, MatchEnv* env);

/* Publish the head of the definition whose rule is about to be matched, so a
 * Condition guard can resolve OptionValue[...] against the rule's
 * OptionsPattern; returns the previous head for the caller to restore. */
const char* match_set_rule_head(const char* head);
/* The current rule head (saved/restored across a TimeConstrained unwind). */
const char* match_get_rule_head(void);

// Replace bound variables in expr with their values from env.
// Returns a new expression (caller must free).
Expr* replace_bindings(Expr* expr, MatchEnv* env);

// Flat-head "leftover" rule application, Phase 2. Call this ONLY after an
// ordinary match(expr, pattern, env) has already failed (keep that match inline
// in the caller so the common path pays nothing). If `expr`'s head is Flat
// (Plus/Times/...) and it has MORE operands than the pattern's flat-call arity,
// a Flat-head pattern matches a SUBSET and the surplus is re-wrapped in the head
// (Mathematica semantics); returns the replaced result (env left populated for
// post-processing) or NULL. `env` must be a fresh per-rule environment.
Expr* match_flat_leftover(Expr* expr, Expr* pattern, Expr* replacement, MatchEnv* env);

Expr* builtin_matchq(Expr* res);

#endif // MATCH_H
