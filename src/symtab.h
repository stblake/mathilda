#ifndef SYMTAB_H
#define SYMTAB_H

#include "expr.h"
#include "eval.h"

// A simple rule in a rewrite system (e.g. pattern -> replacement)
typedef struct Rule {
    Expr* pattern;
    Expr* replacement;
    /* §3.5 DownValue dispatch index. Cheap pre-computed key used by
     * apply_down_values to skip rules whose top-level shape cannot match
     * the input before invoking the matcher.
     *
     *   dispatch_arity  : number of args the pattern's top-level call
     *                     accepts. -1 means "variable arity" because the
     *                     pattern contains a top-level BlankSequence,
     *                     BlankNullSequence, Optional, Repeated,
     *                     RepeatedNull or OptionsPattern. A wildcard arity
     *                     also forces first_arg_head_canon to NULL.
     *   first_arg_head_canon : canonical (interned) head of the first
     *                     argument's pattern -- e.g. "Plus" for
     *                     f[Plus[x_,y_], _], "Integer" for f[5, _], or
     *                     "Integer" for f[Pattern[n, Blank[Integer]], _].
     *                     NULL means the first arg accepts any head, so
     *                     the rule passes the filter unconditionally.
     *
     * Both keys may also be set to wildcard for OwnValues (where dispatch
     * does not apply); the filter then degrades to a no-op. */
    int32_t dispatch_arity;
    const char* first_arg_head_canon;
    /* §3.6 specificity score. Larger = more specific. Rules in the list
     * are sorted by descending specificity, with insertion order
     * preserved as the tie-breaker for equal scores. */
    int32_t specificity;
    struct Rule* next;
} Rule;

// Typedef for built-in C functions evaluating expressions
typedef Expr* (*BuiltinFunc)(Expr* res);

// A symbol definition containing rules associated with the symbol
typedef struct {
    char* symbol_name;
    Rule* own_values;   // x = 4
    Rule* down_values;  // f[x_] = x + 1
    BuiltinFunc builtin_func; // C function evaluating this symbol
    uint32_t attributes;
    char* docstring;
} SymbolDef;

// Initialize the symbol table
void symtab_init(void);

// Clear and free the entire symbol table
void symtab_clear(void);

// Clear own_values and down_values for a specific symbol
void symtab_clear_symbol(const char* symbol_name);

// Get definition for a symbol, or create it if it doesn't exist
SymbolDef* symtab_get_def(const char* symbol_name);

// Look up a symbol without creating one. Returns NULL if not present.
// Used by the context resolver to test whether a qualified name exists.
SymbolDef* symtab_lookup(const char* symbol_name);

// Register a built-in C function for a symbol
void symtab_add_builtin(const char* symbol_name, BuiltinFunc func);

// Set/get docstring for a symbol
void symtab_set_docstring(const char* symbol_name, const char* doc);
const char* symtab_get_docstring(const char* symbol_name);

// Add an OwnValue rule (e.g., x = 4)
void symtab_add_own_value(const char* symbol_name, Expr* pattern, Expr* replacement);

// Add a DownValue rule (e.g., f[x_] = x + 1)
void symtab_add_down_value(const char* symbol_name, Expr* pattern, Expr* replacement);

// Get rules
Rule* symtab_get_own_values(const char* symbol_name);
Rule* symtab_get_down_values(const char* symbol_name);

// Apply DownValues to an expression f[...]
// Returns new evaluated expression if a rule applied, else NULL
Expr* apply_down_values(Expr* expr);

// Apply OwnValues to a symbol x
// Returns new evaluated expression if a rule applied, else NULL
Expr* apply_own_values(Expr* expr);

#endif // SYMTAB_H
