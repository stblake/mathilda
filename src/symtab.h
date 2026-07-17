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

// A symbol definition containing rules associated with the symbol.
//
// Phase 2 (EVAL_SYMTAB_IMPROVEMENTS): the string interner and the symbol table
// are unified -- this node IS the intern entry. A single hash table keyed on the
// name backs both interning (canonical `const char*`) and definition lookup, so
// resolving a symbol to its definition no longer hashes twice across two tables.
typedef struct SymbolDef {
    char* symbol_name;  // owned: the canonical interned name (freed by intern_clear)
    Rule* own_values;   // x = 4
    Rule* down_values;  // f[x_] = x + 1
    BuiltinFunc builtin_func; // C function evaluating this symbol
    uint32_t attributes;
    /* Phase 1 (EVAL_SYMTAB_IMPROVEMENTS): cached immutable "base" attributes --
     * the floor from attr.c's builtin_attrs[] table. get_attributes() returns
     * (base_attributes | attributes) as an O(1) field read + OR instead of the
     * old ~140-entry strcmp scan on every call. The static table never changes,
     * so once folded the base is permanently valid; base_seeded guards the
     * one-time fold (lazily on first get_attributes, eagerly in attr_init for
     * table symbols). `attributes` keeps its exact meaning (dynamic bits, also
     * carrying the base bits that attr_init copies in for match.c's direct
     * reads), so no mutation site needs to change. */
    uint32_t base_attributes;
    uint8_t  base_seeded;
    /* Phase 2 unified-table bookkeeping.
     *   is_system : migrated from the old InternEntry; set by
     *               intern_mark_all_system() for every kernel-interned name.
     *   has_def   : 1 once the name has been materialized as a real definition
     *               (symtab_get_def called). Distinguishes a defined symbol
     *               (enumerated by Names[], visible to symtab_lookup) from an
     *               interned-only node whose name was merely canonicalized.
     *               symtab_remove_symbol / symtab_clear reset this to 0 (the
     *               node and its name survive so live Exprs' pointers stay
     *               valid); intern_clear frees the node outright. */
    uint8_t  is_system;
    uint8_t  has_def;
    char* docstring;
    /* Default option settings for this symbol, the DefaultValues-equivalent
     * backing Options[f] / SetOptions[f]. A List[Rule[name,val], ...] owned
     * by the SymbolDef, or NULL when the symbol has no registered options.
     * Survives Clear[f] (only rules are cleared), freed on symbol removal. */
    Expr* default_options;
    /* Optional machine-precision element-wise kernels for NDArray arguments.
     * Opaque (const NDUnaryKernel* / const NDBinaryKernel* from ndarray.h) to
     * keep symtab.h free of an ndarray.h include cycle; the evaluator casts.
     * NULL (the calloc default) means "no NDArray fast path for this head".
     * Point at static-const descriptors, so the SymbolDef never owns them. */
    const void* ndarray_unary_kernel;
    const void* ndarray_binary_kernel;
    struct SymbolDef* next;   // hash-bucket chain (replaces the old SymEntry)
} SymbolDef;

// Initialize the symbol table
void symtab_init(void);

// Clear and free the entire symbol table
void symtab_clear(void);

// Clear own_values and down_values for a specific symbol
void symtab_clear_symbol(const char* symbol_name);

// Remove a symbol entirely: free its definition (values, attributes,
// docstring) and unlink it from the table. The interned name itself is
// owned by sym_intern and is NOT freed, so a later reference recreates a
// fresh, empty definition. No-op if the symbol is not present.
void symtab_remove_symbol(const char* symbol_name);

// Remove the single OwnValue (own_value=true) or DownValue (own_value=false)
// rule on `symbol_name` whose left-hand side is identical to `lhs` up to
// renaming of bound pattern variables (alpha-equivalence). Backs Unset /
// `lhs =.`. Returns true iff a matching rule was found and removed.
bool symtab_remove_matching_rule(const char* symbol_name, const Expr* lhs,
                                 bool own_value);

// Get definition for a symbol, or create it if it doesn't exist
SymbolDef* symtab_get_def(const char* symbol_name);

// Look up a symbol without creating one. Returns NULL if not present.
// Used by the context resolver to test whether a qualified name exists.
SymbolDef* symtab_lookup(const char* symbol_name);

// Visit every symbol currently in the table, in unspecified order. `name` is
// the interned (borrowed) symbol name -- do not free it. `user` is passed
// through untouched. Used by Names[] to enumerate the symbol namespace.
void symtab_for_each(void (*visit)(const char* name, SymbolDef* def,
                                   void* user), void* user);

// Register a built-in C function for a symbol
void symtab_add_builtin(const char* symbol_name, BuiltinFunc func);

// Set/get docstring for a symbol
void symtab_set_docstring(const char* symbol_name, const char* doc);
const char* symtab_get_docstring(const char* symbol_name);

// Register an NDArray element-wise machine kernel for a function head. `kernel`
// points at a static-const NDUnaryKernel / NDBinaryKernel (see ndarray.h); the
// SymbolDef stores it opaquely and never owns it. The evaluator maps it over an
// NDArray argument at C speed instead of falling back to List threading.
void symtab_set_ndarray_unary_kernel(const char* symbol_name, const void* kernel);
void symtab_set_ndarray_binary_kernel(const char* symbol_name, const void* kernel);

// Set/get the default option settings (a List[Rule[name,val], ...]) for a
// symbol -- the store behind Options[f] and SetOptions[f]. symtab_set_options
// takes ownership of `list` (freeing any previously stored options); pass NULL
// to clear. symtab_get_options returns a borrowed pointer (do not free),
// or NULL when the symbol has no registered options.
void symtab_set_options(const char* symbol_name, Expr* list);
Expr* symtab_get_options(const char* symbol_name);

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

// Phase 3a: def-threaded variant. `def` must be the definition of `expr`'s
// (symbol) head; lets the evaluator avoid re-resolving the head. NULL def or
// non-function expr returns NULL.
Expr* apply_down_values_def(SymbolDef* def, Expr* expr);

// Apply OwnValues to a symbol x
// Returns new evaluated expression if a rule applied, else NULL
Expr* apply_own_values(Expr* expr);

#endif // SYMTAB_H
