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
    /* §3.6 specificity score: a cheap key that rules out most candidates in
     * add_rule's alpha-equivalent-redefinition search.  It does not order
     * the list -- pattern generality does (see pattern_subsumes). */
    int32_t specificity;
    /* Hash of the pattern up to renaming of its pattern variables, so an
     * alpha-equivalent redefinition is found without normalising every
     * stored rule (see add_rule). Computed lazily; 0 = not yet computed. */
    uint64_t alpha_hash;
    /* 1 when the rule was made by a delayed assignment (:=, ^:=, /: :=), 0
     * for an immediate one (=, ^=, /: =). Presentation only: Definition
     * renders it; matching and rewriting never read it. */
    uint8_t delayed;
    /* 1 when the replacement is a Module / With / Block whose body ends (in
     * tail position) in `val /; cond`: the rule is CONDITIONAL -- the guard
     * is tested inside the scope when the rule is tried, and a failing guard
     * falls through to the next rule.  Such a rule is not replaced by a later
     * definition with the same LHS. */
    uint8_t scoped_cond;
    /* 1 when the replacement mentions Return: applying the rule then
     * evaluates the result and yields v for a top-level Return[v], as the
     * reference returns from a definition (f[x_] := If[x > 0, Return[1], 2]). */
    uint8_t has_return;
    /* Ordering facts about the pattern, fixed at insertion: 1 when it holds
     * no pattern object at all (f[1]; such rules lead the list), and 1 when
     * it names a pattern variable twice (f[x_, x_]; never taken as more
     * general than another pattern). */
    uint8_t literal;
    uint8_t nonlinear;
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
    Rule* sub_values;   // f[x_][y_] := ... -- keyed on the innermost head symbol
#if UP_VALUES
    /* UpValues: rules keyed on THIS symbol but matched against an ENCLOSING
     * call in which this symbol appears at level one. Installed by UpSet (^=),
     * UpSetDelayed (^:=), TagSet / TagSetDelayed (g /: lhs = / := rhs). Tried
     * before the enclosing head's DownValues/builtin -- see apply_up_values. */
    Rule* up_values;
#endif
    /* FormatValues: display rules keyed on this symbol, installed by
     * Format[lhs] := rhs and tag /: MakeBoxes[lhs, form] := rhs (or Format).
     * Never consulted by evaluation; the OutputForm printer applies the ones
     * that are not MakeBoxes rules. */
    Rule* format_values;
    /* NValues: N[x] = v, N[f[...], p_] := ..., x /: N[x] = v -- the values N
     * gives the symbol x or calls of f. An N[lhs] rule is stored as
     * N[lhs, {MachinePrecision, MachinePrecision}], as the reference shows
     * it. Consulted only by N (numeric.c n_value). */
    Rule* n_values;
    BuiltinFunc builtin_func; // C function evaluating this symbol
    uint32_t attributes;
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
    uint8_t  is_system;      /* 1: System`, 2: Mathilda` (Mathilda-only), 0: neither */
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
    const void* ndarray_nary_kernel;   /* const NDNaryKernel*, arity >= 3 */
    /* Packed-list opt-in. A PACKED LIST is an ordinary List stored as a dense
     * buffer (EXPR_NDARRAY, present_as == NDA_HEAD_LIST -- see pack.h). It is
     * a single node with no args[] array, so a builtin that walks
     * arg->data.function.args past one either crashes or, worse, takes the
     * "not a function" branch and returns a confidently wrong answer:
     * Count[packed, _Real] would say 0.
     *
     * There are ~7100 sites in the tree that touch data.function.args
     * directly, so this cannot be retrofitted behind an accessor. Instead the
     * evaluator MATERIALISES every packed argument before calling any head
     * that has not set this flag. That makes the ~400 files that know nothing
     * about packing correct by construction, at the cost of one unpack -- and
     * it means opting a head in is a claim that it has actually been checked,
     * not a default.
     *
     * Implied for any head with an ndarray kernel above; set explicitly with
     * symtab_set_packed_aware for the rest. */
    uint8_t  packed_aware;
    /* Narrower still: this head is packed-aware AND exact on an NDT_INT64
     * buffer. The generic ndt_get/ndt_set pair routes through `double` and is
     * exact only to 2^53, and most of the ND layer was written when only
     * Compile[] could create an integer buffer, so it uses that pair freely.
     * Packed lists infer int64 from an all-Integer list, which puts integer
     * buffers in front of that code for the first time.
     *
     * Until each site has an exact int64 path, the gate materialises an int64
     * packed argument for every head NOT marked here -- so `Total[{1,2,3}]`
     * stays the Integer 6 rather than 6., and `Sin[{1,2,3}]` stays symbolic
     * rather than becoming {0,0,0}. Both of those were live before this bit
     * existed. Marking a head here is a claim it reads the buffer through
     * ndt_get_i / memcpy only. */
    uint8_t  packed_int64_ok;
    /* This head's buffer path handles operands of DIFFERENT RANK, by pairing
     * the lower-rank one against the leading axes of the higher (Listable
     * threading across levels: {{1,2,3},{4,5,6}} + {10,20}). Only Plus and
     * Times, via ndarray_elementwise's broadcast pre-pass.
     *
     * Consulted only by the gate in evaluate_step, and only to decide whether a
     * plain numeric List argument may be PACKED UP to meet a buffer instead of
     * the buffer being materialised down to meet it. Same-shape operands need no
     * flag -- every ND kernel handles those -- but a rank mismatch is dispatched
     * differently by different heads (ndarray_elementwise_power still requires
     * same_shape), and packing there would route the call somewhere its kernel
     * declines. So this is the narrow claim: "a rank mismatch is still mine". */
    uint8_t  packed_broadcast_ok;
    /* How this head treats a SparseArray (EXPR_SPARSE) argument: 0 leaves it an
     * atom, otherwise an SPG_* mode from sparse/sparse.h that the evaluator's
     * sparse transparency gate applies (sparse/sparse_gate.c). Opt-in, because
     * the language itself treats a SparseArray as an atom for pattern and
     * container heads. */
    uint8_t  sparse_mode;
    /* A native SparseArray kernel for this head (an SpNativeFn from
     * sparse/sparse.h), tried before the sparse_mode gate whenever an argument
     * is a SparseArray. It borrows the call and returns the answer, or NULL to
     * decline -- and then the gate (if any) answers. NULL: no kernel. */
    const void* sparse_native;
    /* Hash index of the leading pattern-free DownValues (f[1] = ..., the
     * memo entries of f[n_] := f[n] = ...), private to symtab.c: lookup and
     * redefinition of such a rule are O(1) instead of a scan of every stored
     * rule.  Valid only while it records the current down_values head;
     * dropped at every removal. NULL: not built. */
    struct LitIndex* lit_index;
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

// SubValues: rules whose left-hand side has a compound head, f[x_][y_] (any
// depth), keyed on the innermost head symbol f.  Tried for a call whose head
// is compound, after the head has been evaluated (see apply_sub_values).
/* A delayed replacement that is a Module/With/Block ending in val /; cond (a
 * scoped guard): does `rhs` have one, and -- for a bound replacement -- its
 * evaluated value, or NULL when the guard did not give True. Shared by the
 * DownValue dispatcher and the rule-applying builtins (ReplaceAll, Cases). */
bool  symtab_rhs_has_scoped_condition(Expr* rhs);
Expr* symtab_try_scoped_condition(Expr* bound);
void symtab_add_sub_value(const char* symbol_name, Expr* pattern, Expr* replacement,
                          bool delayed);
Rule* symtab_get_sub_values(const char* symbol_name);
Expr* apply_sub_values(SymbolDef* def, Expr* expr);
// Unset of a SubValue (alpha-equivalent LHS); true iff one was removed.
bool symtab_remove_matching_sub_value(const char* symbol_name, const Expr* lhs);
// The innermost symbol of a compound head g[x][y] -> "g", or NULL.
const char* symtab_innermost_head(const Expr* e);

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
void symtab_set_ndarray_nary_kernel(const char* symbol_name, const void* kernel);

/* Declare that `symbol_name` handles a PACKED LIST argument correctly -- see
 * SymbolDef.packed_aware. Implied by the kernel setters above. Setting this on
 * a head that has not actually been checked is a silent-wrong-answer bug, so
 * the list of heads that use it lives in one place (pack_mark_aware_heads) and
 * is reviewed as a unit rather than sprinkled across modules. */
void symtab_set_packed_aware(const char* symbol_name);
/* Clear packed_aware. Needed because the ndarray kernel setters set it
 * implicitly, and a handful of kernels change an element's HEAD (Floor,
 * Ceiling, Round, IntegerPart, Sign, Im, Clip) or collapse a Listable answer to
 * a scalar (Precision, Accuracy) -- correct for a visible NDArray[...], wrong
 * for a packed List, which must answer exactly as the plain List does. */
void symtab_clear_packed_aware(const char* symbol_name);

/* Additionally declare that `symbol_name` is EXACT on an int64 buffer -- see
 * SymbolDef.packed_int64_ok. Implies packed_aware. */
void symtab_set_packed_int64_ok(const char* symbol_name);
/* SymbolDef.packed_broadcast_ok. Implies packed_aware. */
void symtab_set_packed_broadcast_ok(const char* symbol_name);
/* SymbolDef.sparse_mode (an SPG_* value from sparse/sparse.h). */
void symtab_set_sparse_mode(const char* symbol_name, int mode);
/* SymbolDef.sparse_native (an SpNativeFn from sparse/sparse.h). */
void symtab_set_sparse_native(const char* symbol_name, const void* fn);

// Set/get the default option settings (a List[Rule[name,val], ...]) for a
// symbol -- the store behind Options[f] and SetOptions[f]. symtab_set_options
// takes ownership of `list` (freeing any previously stored options); pass NULL
// to clear. symtab_get_options returns a borrowed pointer (do not free),
// or NULL when the symbol has no registered options.
void symtab_set_options(const char* symbol_name, Expr* list);
Expr* symtab_get_options(const char* symbol_name);

// Add an OwnValue rule (e.g., x = 4)
void symtab_add_own_value(const char* symbol_name, Expr* pattern, Expr* replacement);
/* Count of symtab_init calls so far; a load-once memo keys itself on it. */
unsigned symtab_generation(void);

// Add a DownValue rule (e.g., f[x_] = x + 1)
void symtab_add_down_value(const char* symbol_name, Expr* pattern, Expr* replacement);

// The same, recording whether the rule came from a delayed assignment
// (SetDelayed / UpSetDelayed / TagSetDelayed) or an immediate one. The flag only
// affects how Definition renders the rule (`=` vs `:=`, `^=` vs `^:=`). The
// plain forms above default to immediate for OwnValues, delayed otherwise.
void symtab_add_own_value_as(const char* symbol_name, Expr* pattern, Expr* replacement, bool delayed);
void symtab_add_down_value_as(const char* symbol_name, Expr* pattern, Expr* replacement, bool delayed);

// Get rules
Rule* symtab_get_own_values(const char* symbol_name);
Rule* symtab_get_down_values(const char* symbol_name);
// Unlink the DownValue whose LHS is exactly `lhs` (canonicalised); caller owns
// the returned rule (or NULL). Free a detached rule with symtab_rule_free.
Rule* symtab_detach_down_value_exact(const char* symbol_name, Expr* lhs);
void  symtab_rule_free(Rule* r);

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

#if UP_VALUES
/* ============================================================
 * UpValues (trial, gated on UP_VALUES).
 *
 * An UpValue associated with symbol `s` fires when `s` appears at LEVEL ONE of
 * an enclosing call g[..., s, ...] / g[..., s[...], ...] -- i.e. `s` is a direct
 * argument, or the head of a direct argument. Matching UpValues are tried BEFORE
 * the enclosing head g's DownValues and builtin.
 * ============================================================ */

/* Running total of UpValue rules across all symbols. The evaluator's per-call
 * hook early-outs on `== 0`, so a tree with no UpValues pays a single load+branch
 * before DownValue dispatch. Maintained by the add/remove/set helpers below. */
extern size_t symtab_up_value_count;

// Add an UpValue rule on `symbol_name` (the tag), with `pattern` the enclosing
// call to match and `replacement` its rewrite. Bumps the eval clock + rule epoch.
void symtab_add_up_value(const char* symbol_name, Expr* pattern, Expr* replacement);
void symtab_add_up_value_as(const char* symbol_name, Expr* pattern, Expr* replacement, bool delayed);

// Borrowed UpValue list for a symbol (may be NULL).
Rule* symtab_get_up_values(const char* symbol_name);

// Apply UpValues during evaluation of the function call `expr`. Borrows `expr`;
// returns a freshly-owned rewrite if some level-one symbol's UpValue matched,
// else NULL (same ownership contract as apply_down_values_def). First match
// across the level-one candidates, left-to-right, wins.
Expr* apply_up_values(Expr* expr);

// Remove the single UpValue rule on `symbol_name` whose LHS is alpha-equivalent
// to `lhs`. Backs TagUnset (f /: lhs =.). Returns true iff one was removed.
bool symtab_remove_matching_up_value(const char* symbol_name, const Expr* lhs);

/* Replace a symbol's entire OwnValue / DownValue / UpValue list from a List of
 * rules (Rule/RuleDelayed[lhs, rhs], an outer HoldPattern on lhs stripped).
 * Backs OwnValues[f]=list / DownValues[f]=list / UpValues[f]=list. The previous
 * list is freed. Bumps the eval clock + rule epoch. */
void symtab_set_own_values(const char* symbol_name, const Expr* list);
void symtab_set_down_values(const char* symbol_name, const Expr* list);
void symtab_set_up_values(const char* symbol_name, const Expr* list);
/* FormatValues[f] = list / NValues[f] = list (NULL or {} clears); share the
 * backend above. */
void symtab_set_format_values(const char* symbol_name, const Expr* list);
void symtab_set_n_values(const char* symbol_name, const Expr* list);
#endif // UP_VALUES

/* FormatValues. symtab_format_value_count is the total number of rules across
 * every symbol, so the printer's format hook early-outs on 0. Add copies its
 * arguments (a same-LHS rule is replaced). */
extern size_t symtab_format_value_count;
void  symtab_add_format_value(const char* symbol_name, Expr* pattern, Expr* replacement, bool delayed);
Rule* symtab_get_format_values(const char* symbol_name);

/* NValues: as FormatValues (see SymbolDef.n_values). */
extern size_t symtab_n_value_count;
void  symtab_add_n_value(const char* symbol_name, Expr* pattern, Expr* replacement, bool delayed);
Rule* symtab_get_n_values(const char* symbol_name);

#endif // SYMTAB_H
