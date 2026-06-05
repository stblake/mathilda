#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "symtab.h"
#include "match.h"
#include "print.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "eval.h"
#include "arithmetic.h"
#include <string.h>
#include <stdlib.h>

#define SYMTAB_SIZE 65535

typedef struct SymEntry {
    SymbolDef* def;
    struct SymEntry* next;
} SymEntry;

static SymEntry* symtab[SYMTAB_SIZE] = {0};

static unsigned int hash(const char* s) {
    unsigned int h = 5381;
    while (*s) h = ((h << 5) + h) + *s++;
    return h % SYMTAB_SIZE;
}

void symtab_init(void) {
    memset(symtab, 0, sizeof(symtab));
    /* Populate the SYM_* canonical-pointer cache so that any code path
     * that compares e->data.symbol against SYM_* (match.c, eval.c, the
     * arithmetic heads, etc.) works as soon as the symtab is live -- not
     * only after core_init() has finished registering builtins. Tests
     * that initialise just the symtab depend on this. */
    sym_names_init();
}

/* Both lookup helpers below first intern the input name, then key on
 * the canonical pointer. This collapses the per-bucket strcmp into a
 * pointer compare and is correct because the bucket entries' symbol_name
 * fields are themselves interned at creation time. */
SymbolDef* symtab_get_def(const char* symbol_name) {
    const char* canon = intern_symbol(symbol_name);
    unsigned int idx = hash(canon);
    SymEntry* entry = symtab[idx];
    while (entry) {
        if (entry->def->symbol_name == canon) {
            return entry->def;
        }
        entry = entry->next;
    }

    // Create new symbol definition
    SymbolDef* def = malloc(sizeof(SymbolDef));
    def->symbol_name = (char*)canon;  /* interned; symtab does not own the storage */
    def->own_values = NULL;
    def->down_values = NULL;
    def->builtin_func = NULL;
    def->attributes = 0;
    def->docstring = NULL;

    SymEntry* new_entry = malloc(sizeof(SymEntry));
    new_entry->def = def;
    new_entry->next = symtab[idx];
    symtab[idx] = new_entry;

    return def;
}

SymbolDef* symtab_lookup(const char* symbol_name) {
    const char* canon = intern_symbol(symbol_name);
    unsigned int idx = hash(canon);
    SymEntry* entry = symtab[idx];
    while (entry) {
        if (entry->def->symbol_name == canon) {
            return entry->def;
        }
        entry = entry->next;
    }
    return NULL;
}

void symtab_add_builtin(const char* symbol_name, BuiltinFunc func) {
    SymbolDef* def = symtab_get_def(symbol_name);
    def->builtin_func = func;
}

void symtab_set_docstring(const char* symbol_name, const char* doc) {
    SymbolDef* def = symtab_get_def(symbol_name);
    if (def->docstring) free(def->docstring);
    def->docstring = doc ? strdup(doc) : NULL;
}

const char* symtab_get_docstring(const char* symbol_name) {
    SymbolDef* def = symtab_get_def(symbol_name);
    return def->docstring;
}

/* ============================================================
 * §3.5/§3.6: dispatch-key extraction and specificity scoring.
 *
 * These run once per rule at insertion time. The output is a cheap
 * pre-computed key (dispatch_arity, first_arg_head_canon) and an
 * integer specificity score that lets apply_down_values short-circuit
 * unmatchable rules and the rule list stay sorted "most specific first."
 * ============================================================ */

/* Strip transparent pattern wrappers down to the inner pattern.
 *
 * - HoldPattern / Verbatim are transparent for matching, so the
 *   dispatch key should look through them.
 * - Condition[p, c] and PatternTest[p, t] are guards that gate the
 *   match without altering its STRUCTURAL shape — the dispatch
 *   filter wants the arity and first-arg-head of the inner pattern
 *   `p`, not of the Condition/PatternTest wrapper itself.  The
 *   matcher (match.c:278-313 for Condition, 394-417 for PatternTest)
 *   evaluates the guard post-match; we just need dispatch not to
 *   pre-filter the rule away on a confused arity/head key.
 *
 * Stripping Condition here was the missing piece that prevented
 * every `f[...] /; cond := rhs` DownValue from firing (pre-fix the
 * rule's dispatch_arity was the Condition's arg count and its
 * first_arg_head_canon was the inner pattern's HEAD symbol — the
 * filter then routinely disagreed with an input's actual first-arg
 * head and skipped the rule entirely).  See tests/test_condition_
 * downvalue.c for the pinning suite. */
static const Expr* strip_pattern_wrappers(const Expr* p) {
    while (p && p->type == EXPR_FUNCTION && p->data.function.head &&
           p->data.function.head->type == EXPR_SYMBOL &&
           p->data.function.arg_count >= 1) {
        const char* h = p->data.function.head->data.symbol;
        if (h == SYM_HoldPattern || h == SYM_Verbatim) {
            p = p->data.function.args[0];
            continue;
        }
        if ((h == SYM_Condition || h == SYM_PatternTest) &&
            p->data.function.arg_count == 2) {
            p = p->data.function.args[0];
            continue;
        }
        break;
    }
    return p;
}

/* True when `p`, considered as a top-level argument slot in a DownValue
 * pattern, can absorb a variable number of input arguments. Patterns
 * that satisfy this force dispatch_arity = -1. */
static bool slot_is_variable_length(const Expr* p) {
    p = strip_pattern_wrappers(p);
    if (!p) return false;
    if (p->type != EXPR_FUNCTION || !p->data.function.head ||
        p->data.function.head->type != EXPR_SYMBOL) {
        return false;
    }
    const char* h = p->data.function.head->data.symbol;
    if (h == SYM_BlankSequence) return true;
    if (h == SYM_BlankNullSequence) return true;
    if (h == SYM_Optional) return true;
    if (h == SYM_Repeated) return true;
    if (h == SYM_RepeatedNull) return true;
    if (h == SYM_OptionsPattern) return true;
    /* Pattern[x, q] is variable iff q is. */
    if (h == SYM_Pattern && p->data.function.arg_count == 2) {
        return slot_is_variable_length(p->data.function.args[1]);
    }
    return false;
}

/* True when `p` is an Optional pattern (`x_.`, `Optional[p]`,
 * `Optional[p, def]`), possibly under a `Pattern[name, ...]` wrapper.
 * An Optional argument in a function-call pattern can be omitted, which
 * collapses the call to its base form (`Power[b_, n_.]` matches a bare
 * `b` as `b^1`; `Times[a_., x_]` matches `x`), so the concrete value
 * need not share the pattern's head. */
static bool pattern_arg_is_optional(const Expr* p) {
    if (!p || p->type != EXPR_FUNCTION || !p->data.function.head ||
        p->data.function.head->type != EXPR_SYMBOL) {
        return false;
    }
    const char* h = p->data.function.head->data.symbol;
    if (h == SYM_Optional) return true;
    if (h == SYM_Pattern && p->data.function.arg_count == 2) {
        return pattern_arg_is_optional(p->data.function.args[1]);
    }
    return false;
}

/* Extract the canonical head symbol that any concrete value matching
 * this pattern slot must have. Returns NULL when the slot is wildcard
 * with respect to head (Blank[], _, anything sequence-like, or anything
 * we are not willing to specialize on). */
static const char* pattern_arg_head_canon(const Expr* p) {
    p = strip_pattern_wrappers(p);
    if (!p) return NULL;

    /* Atoms in patterns match themselves -- their head equals the head
     * of the corresponding atomic input value. */
    if (p->type == EXPR_INTEGER || p->type == EXPR_BIGINT) {
        return intern_symbol("Integer");
    }
    if (p->type == EXPR_REAL) return intern_symbol("Real");
    if (p->type == EXPR_STRING) return intern_symbol("String");
    if (p->type == EXPR_SYMBOL) {
        /* Bare symbol used as a literal (e.g. f[True, x_]). Head is
         * "Symbol" -- the filter will only reject inputs whose first arg
         * is not a symbol (e.g. an integer or a function call). */
        return intern_symbol("Symbol");
    }
    if (p->type != EXPR_FUNCTION) return NULL;

    Expr* head = p->data.function.head;
    if (!head || head->type != EXPR_SYMBOL) return NULL;
    const char* h = head->data.symbol;

    /* Blank[] -> wildcard; Blank[h] -> h. */
    if (h == SYM_Blank) {
        if (p->data.function.arg_count == 0) return NULL;
        Expr* head_arg = p->data.function.args[0];
        if (head_arg && head_arg->type == EXPR_SYMBOL) {
            return intern_symbol(head_arg->data.symbol);
        }
        return NULL;
    }
    /* Pattern[x, q] -> use q's head. */
    if (h == SYM_Pattern && p->data.function.arg_count == 2) {
        return pattern_arg_head_canon(p->data.function.args[1]);
    }
    /* Sequence-like or optional patterns are wildcards w.r.t. head. */
    if (h == SYM_BlankSequence || h == SYM_BlankNullSequence ||
        h == SYM_Optional || h == SYM_Repeated ||
        h == SYM_RepeatedNull || h == SYM_OptionsPattern) {
        /* For typed BlankSequence[h] etc., picking up h would still be
         * unsafe because the matcher binds whole sequences. Stay safe. */
        return NULL;
    }
    /* Condition[p, c] / PatternTest[p, t] are transparent for head. */
    if ((h == SYM_Condition || h == SYM_PatternTest) &&
        p->data.function.arg_count >= 1) {
        return pattern_arg_head_canon(p->data.function.args[0]);
    }
    /* Alternatives[a, b, ...] -- only safe to specialize when every
     * branch agrees on a head. */
    if (h == SYM_Alternatives) {
        const char* acc = NULL;
        for (size_t i = 0; i < p->data.function.arg_count; i++) {
            const char* hi = pattern_arg_head_canon(p->data.function.args[i]);
            if (!hi) return NULL;          /* one branch is wildcard */
            if (acc && acc != hi) return NULL; /* branches disagree */
            acc = hi;
        }
        return acc;
    }
    /* Plain function call in a pattern -- head must match, UNLESS one of
     * its arguments is Optional. An omitted Optional collapses the call to
     * its base (e.g. Power[b_, n_.] matches a bare `b` as b^1, Times[a_., x_]
     * matches `x`), so a concrete value need not share this head. Returning
     * the head here would make the §3.5 dispatch filter skip such rules for
     * any input that is not literally headed by `h`. */
    for (size_t i = 0; i < p->data.function.arg_count; i++) {
        if (pattern_arg_is_optional(p->data.function.args[i])) return NULL;
    }
    return intern_symbol(h);
}

/* §3.5 dispatch_arity: count of top-level arg slots, or -1 when any
 * slot is variable-length. Pattern is the WHOLE LHS, e.g. f[a, b]; we
 * inspect the top-level function's args. */
static int32_t pattern_dispatch_arity(const Expr* pattern) {
    pattern = strip_pattern_wrappers(pattern);
    if (!pattern || pattern->type != EXPR_FUNCTION) return -1;
    size_t n = pattern->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        if (slot_is_variable_length(pattern->data.function.args[i])) return -1;
    }
    return (int32_t)n;
}

/* §3.5 first_arg_head_canon: NULL when arity is variable (a leading
 * sequence pattern would absorb the first arg) OR the first slot is
 * wildcard. */
static const char* pattern_first_arg_head(const Expr* pattern) {
    pattern = strip_pattern_wrappers(pattern);
    if (!pattern || pattern->type != EXPR_FUNCTION) return NULL;
    if (pattern->data.function.arg_count == 0) return NULL;
    /* If any slot at or before the first is variable-length, the "first
     * arg" of the input might end up bound to the sequence. Be safe. */
    if (slot_is_variable_length(pattern->data.function.args[0])) return NULL;
    return pattern_arg_head_canon(pattern->data.function.args[0]);
}

/* §3.6 specificity score. Computed recursively over the pattern tree.
 *
 *   literal atoms                  : +100 (most discriminating)
 *   typed Blank[h]                 :  +20
 *   plain Blank[]                  :  -10
 *   BlankSequence[]                : -100
 *   BlankNullSequence[]            : -200
 *   Optional[...]                  :  -50
 *   Repeated[p]                    : score(p) - 50
 *   RepeatedNull[p]                : score(p) - 100
 *   Condition / PatternTest        : score(inner) + 10 (extra constraint)
 *   HoldPattern / Verbatim         : score(inner) (transparent)
 *   Alternatives[a, b, ...]        : min(score(a_i)) (weakest branch wins)
 *   plain function call f[a, b...] :  +50 + sum(score(arg_i))
 *   Pattern[name, q]               : score(q) (named binding adds nothing)
 *
 * Scores are added across siblings, so a 3-arg literal call beats a
 * 2-arg literal call. Insertion order is used as the tie-break, so two
 * structurally-similar pattern rules keep the user's chosen ordering. */
static int32_t pattern_specificity(const Expr* p) {
    if (!p) return 0;
    if (p->type == EXPR_INTEGER || p->type == EXPR_BIGINT ||
        p->type == EXPR_REAL || p->type == EXPR_STRING) {
        return 100;
    }
    if (p->type == EXPR_SYMBOL) {
        return 100;
    }
    if (p->type != EXPR_FUNCTION) return 0;

    Expr* head = p->data.function.head;
    if (head && head->type == EXPR_SYMBOL) {
        const char* h = head->data.symbol;
        size_t ac = p->data.function.arg_count;

        if (h == SYM_HoldPattern || h == SYM_Verbatim) {
            return ac >= 1 ? pattern_specificity(p->data.function.args[0]) : 0;
        }
        if (h == SYM_Pattern) {
            /* Pattern[name, q] -- the name is irrelevant for matching power. */
            return ac >= 2 ? pattern_specificity(p->data.function.args[1]) : 0;
        }
        if (h == SYM_Blank) {
            if (ac == 0) return -10;
            return 20; /* typed blank */
        }
        if (h == SYM_BlankSequence) return -100;
        if (h == SYM_BlankNullSequence) return -200;
        if (h == SYM_Optional) return -50;
        if (h == SYM_Repeated) {
            return (ac >= 1 ? pattern_specificity(p->data.function.args[0]) : 0) - 50;
        }
        if (h == SYM_RepeatedNull) {
            return (ac >= 1 ? pattern_specificity(p->data.function.args[0]) : 0) - 100;
        }
        if (h == SYM_Condition || h == SYM_PatternTest) {
            return (ac >= 1 ? pattern_specificity(p->data.function.args[0]) : 0) + 10;
        }
        if (h == SYM_Alternatives) {
            if (ac == 0) return 0;
            int32_t best = pattern_specificity(p->data.function.args[0]);
            for (size_t i = 1; i < ac; i++) {
                int32_t s = pattern_specificity(p->data.function.args[i]);
                if (s < best) best = s; /* weakest branch wins */
            }
            return best;
        }
        if (h == SYM_OptionsPattern) return -150;
    }

    /* Plain function call: literal head + sum of children. */
    int32_t score = 50;
    if (head) score += pattern_specificity(head);
    for (size_t i = 0; i < p->data.function.arg_count; i++) {
        score += pattern_specificity(p->data.function.args[i]);
    }
    return score;
}

/* ============================================================
 * §3.4 pattern canonicalization
 *
 * Stored DownValue / OwnValue patterns must share canonical form with
 * the runtime expressions they will eventually be matched against;
 * otherwise the §3.5 dispatch filter -- and in some cases the matcher
 * itself -- can silently reject a rule that the abstract pattern would
 * accept (i.e. one whose MatchQ form returns True against the input).
 *
 * The asymmetry exists because `Set` and `SetDelayed` deliberately hold
 * the LHS: every runtime input flowing into `apply_down_values` has
 * already been canonicalized by `evaluate_step` -- Flat heads flattened,
 * Orderless heads sorted, `Power[Times[a,...,z], n_Integer]` distributed
 * over the factors (src/power.c integer-exponent fast path), nested
 * `Power[Power[base, e1], e2]` collapsed -- but the stored pattern was
 * never put through that pass. The pre-2026-05 dispatch filter then
 * compared the parsed-but-uncanonicalized first-arg head of the rule
 * (e.g. `Power`) against the evaluated first-arg head of the input
 * (e.g. `Times`), found them unequal, and skipped the rule before the
 * matcher ran.
 *
 * The concrete bug that motivated this pass was the CRC integral table:
 * rules like
 *
 *     IntegrateTable[1/(x_ Sqrt[a_ + b_. x_ + c_. x_^2]), x_] /; ...
 *
 * parse as `Power[Times[x_, Sqrt[...]], -1]` at the first-arg slot, but
 * every runtime call arrives as `Times[Power[x,-1], Power[Sqrt[...],-1]]`
 * because Power-over-Times-integer distribution runs during normal
 * evaluation. The dispatch filter cached `Power` as the rule's first-arg
 * head, the input presented `Times`, and ~80 CRC formulas with similar
 * `1/(...)` LHSs went silently unmatched (see tasks/crc_corpus_2026-05-15.md
 * "Out-of-scope findings" -- which misattributed the failure to a 3-term
 * Plus matcher gap).
 *
 * What this pass structurally changes in the pattern matcher:
 *
 *   - Stored patterns now share canonical form with runtime inputs.
 *     `Power[Times[t1,...,tn], k]` and `Times[Power[t1,k],...,Power[tn,k]]`
 *     are no longer two distinguishable LHSs when `k` is an integer
 *     literal; the second form is the only one that survives storage.
 *   - Likewise, `Power[Power[b, e1], e2]` collapses to `Power[b, e1*e2]`
 *     when both exponents are numeric literals, matching the evaluator's
 *     nested-Power normalization.
 *   - Flat-headed pattern children (Plus, Times, ...) are flattened, and
 *     Orderless-headed pattern children are sorted with
 *     `eval_compare_expr_ptrs` -- the same comparator that orders input
 *     arguments -- so the §3.5 dispatch filter sees consistent shapes on
 *     both sides.
 *   - Pattern-marker wrappers (Pattern, Blank, BlankSequence,
 *     BlankNullSequence, Optional, Condition, PatternTest, HoldPattern,
 *     Verbatim, Alternatives, Repeated, RepeatedNull, OptionsPattern)
 *     are NOT structurally rewritten; the canonicalizer recurses into the
 *     inner pattern slot of each but never folds the wrapper itself. This
 *     preserves match-time semantics in full (a `/; cond` guard still
 *     evaluates with bindings in scope, an Optional default still binds
 *     when the slot is absent, HoldPattern still forces literal matching,
 *     etc.).
 *
 * The pass is idempotent: applying it twice to a pattern yields the same
 * tree as applying it once, so re-storing an already-canonicalized rule
 * is safe.
 *
 * Callers transfer ownership of `p` to `pattern_canonicalize`; the
 * returned tree may share or replace `p`. ============================================================ */
static Expr* pattern_canonicalize(Expr* p);

/* True for an exponent we are willing to fold across nested Powers.
 * Restricted to exact integer / rational literals so the collapse is
 * always sound -- floating-point or symbolic exponents are left alone. */
static bool is_foldable_numeric_exponent(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    int64_t n, d;
    return is_rational(e, &n, &d);
}

/* Returns Times[copies-of-a-and-b] reduced through `evaluate` so the
 * result is in canonical integer/rational form. Used to compose nested
 * Power exponents; both inputs must satisfy `is_foldable_numeric_exponent`. */
static Expr* compose_numeric_exponents(const Expr* a, const Expr* b) {
    Expr* args[2] = { expr_copy((Expr*)a), expr_copy((Expr*)b) };
    Expr* prod = expr_new_function(expr_new_symbol("Times"), args, 2);
    Expr* res = evaluate(prod);
    /* evaluate() makes its own copy of the input — release ours. */
    expr_free(prod);
    return res;
}

static Expr* pattern_canonicalize(Expr* p) {
    if (!p || p->type != EXPR_FUNCTION) return p;

    Expr* head = p->data.function.head;
    if (head && head->type == EXPR_SYMBOL) {
        const char* h = head->data.symbol;

        /* HoldPattern / Verbatim explicitly opt out of canonicalization;
         * their whole purpose is to match a literal subtree without
         * pattern-semantic interpretation. */
        if (h == SYM_HoldPattern || h == SYM_Verbatim) return p;

        /* Blank-family and OptionsPattern carry only type/head metadata
         * in their children -- no nested pattern that wants canonicalization. */
        if (h == SYM_Blank || h == SYM_BlankSequence ||
            h == SYM_BlankNullSequence || h == SYM_OptionsPattern) {
            return p;
        }

        /* Pattern[name, q]: name is held literally; canonicalize q only. */
        if (h == SYM_Pattern) {
            if (p->data.function.arg_count >= 2) {
                p = expr_unshare(p);
                p->data.function.args[1] =
                    pattern_canonicalize(p->data.function.args[1]);
            }
            return p;
        }

        /* Optional[p] / Optional[p, default]: canonicalize the inner
         * pattern. The default is a value, not a pattern; leave it alone
         * to avoid firing evaluator side effects before bindings exist. */
        if (h == SYM_Optional) {
            if (p->data.function.arg_count >= 1) {
                p = expr_unshare(p);
                p->data.function.args[0] =
                    pattern_canonicalize(p->data.function.args[0]);
            }
            return p;
        }

        /* Condition[p, guard] / PatternTest[p, test]: canonicalize the
         * inner pattern. The guard / test runs at match time with the
         * bindings substituted -- canonicalizing it here would force
         * evaluation against still-unbound variables. */
        if (h == SYM_Condition || h == SYM_PatternTest) {
            if (p->data.function.arg_count >= 1) {
                p = expr_unshare(p);
                p->data.function.args[0] =
                    pattern_canonicalize(p->data.function.args[0]);
            }
            return p;
        }

        /* Repeated[p] / RepeatedNull[p]: canonicalize the inner pattern. */
        if (h == SYM_Repeated || h == SYM_RepeatedNull) {
            if (p->data.function.arg_count >= 1) {
                p = expr_unshare(p);
                p->data.function.args[0] =
                    pattern_canonicalize(p->data.function.args[0]);
            }
            return p;
        }

        /* Alternatives[a, b, ...]: each branch is itself a pattern. */
        if (h == SYM_Alternatives) {
            p = expr_unshare(p);
            for (size_t i = 0; i < p->data.function.arg_count; i++) {
                p->data.function.args[i] =
                    pattern_canonicalize(p->data.function.args[i]);
            }
            return p;
        }
    }

    /* Plain function call: canonicalize all children, then apply the
     * structural rewrites that mirror evaluate_step. */
    p = expr_unshare(p);
    if (p->data.function.head) {
        p->data.function.head = pattern_canonicalize(p->data.function.head);
    }
    for (size_t i = 0; i < p->data.function.arg_count; i++) {
        p->data.function.args[i] =
            pattern_canonicalize(p->data.function.args[i]);
    }

    /* Look up structural attributes by head symbol. */
    head = p->data.function.head;
    uint32_t attrs = 0;
    if (head && head->type == EXPR_SYMBOL) {
        SymbolDef* def = symtab_lookup(head->data.symbol);
        if (def) attrs = def->attributes;
    }

    /* Flat: flatten nested same-head args (Plus[Plus[a,b],c] -> Plus[a,b,c]). */
    if (head && head->type == EXPR_SYMBOL && (attrs & ATTR_FLAT)) {
        eval_flatten_args(p, head->data.symbol);
    }

    /* Orderless: sort args using the same comparator the evaluator uses
     * for input expressions, so the §3.5 first-arg-head dispatch key and
     * the matcher both see identical orderings. */
    if (head && head->type == EXPR_SYMBOL && (attrs & ATTR_ORDERLESS)) {
        eval_sort_args(p);
    }

    /* Power[Times[t1,...,tn], k_Integer] -> Times[Power[t1,k],...,Power[tn,k]].
     * Mirrors src/power.c's integer-exponent fast path that fires
     * unconditionally during normal evaluation. */
    if (head && head->type == EXPR_SYMBOL && head->data.symbol == SYM_Power
        && p->data.function.arg_count == 2) {
        Expr* base = p->data.function.args[0];
        Expr* exp  = p->data.function.args[1];
        if (exp && (exp->type == EXPR_INTEGER || exp->type == EXPR_BIGINT)
            && base && base->type == EXPR_FUNCTION
            && base->data.function.head
            && base->data.function.head->type == EXPR_SYMBOL
            && base->data.function.head->data.symbol == SYM_Times) {
            size_t bc = base->data.function.arg_count;
            Expr** new_args = malloc(sizeof(Expr*) * bc);
            for (size_t i = 0; i < bc; i++) {
                Expr* pow_args[2] = {
                    expr_copy(base->data.function.args[i]),
                    expr_copy(exp)
                };
                new_args[i] = expr_new_function(
                    expr_new_symbol("Power"), pow_args, 2);
            }
            Expr* times = expr_new_function(
                expr_new_symbol("Times"), new_args, bc);
            free(new_args);
            expr_free(p);
            /* Recurse on the new Times so its children (each a Power
             * we just synthesized) themselves get folded if applicable
             * (e.g. Power[Power[plus, 1/2], -1] -> Power[plus, -1/2]). */
            return pattern_canonicalize(times);
        }
    }

    /* Power[Power[base, e1], e2] -> Power[base, e1*e2] when both
     * exponents are numeric literals. Matches the evaluator's nested-
     * Power collapse for sound bases. */
    if (head && head->type == EXPR_SYMBOL && head->data.symbol == SYM_Power
        && p->data.function.arg_count == 2) {
        Expr* base = p->data.function.args[0];
        Expr* outer_exp = p->data.function.args[1];
        if (base && base->type == EXPR_FUNCTION
            && base->data.function.head
            && base->data.function.head->type == EXPR_SYMBOL
            && base->data.function.head->data.symbol == SYM_Power
            && base->data.function.arg_count == 2) {
            Expr* inner_base = base->data.function.args[0];
            Expr* inner_exp  = base->data.function.args[1];
            if (is_foldable_numeric_exponent(inner_exp)
                && is_foldable_numeric_exponent(outer_exp)) {
                Expr* new_exp = compose_numeric_exponents(inner_exp, outer_exp);
                Expr* pow_args[2] = { expr_copy(inner_base), new_exp };
                Expr* folded = expr_new_function(
                    expr_new_symbol("Power"), pow_args, 2);
                expr_free(p);
                /* The collapsed Power may now expose a further Power-of-
                 * Times opportunity if `inner_base` itself is a Times;
                 * recurse to let the previous branch fire. */
                return pattern_canonicalize(folded);
            }
        }
    }

    return p;
}

// Helper to add a rule to a list
static void add_rule(Rule** list, Expr* pattern, Expr* replacement) {
    /* Any change to an OwnValue/DownValue list could change the meaning
     * of any expression that touches this symbol, so bump the global
     * eval clock. This invalidates every cached evaluation in one shot
     * (conservative variant of EVAL_IMPROVEMENTS_PLAN §3.3). */
    eval_clock_bump();

    /* §3.4 canonicalize the LHS before deriving any dispatch keys or
     * doing duplicate-rule detection. Held LHSs would otherwise sit in
     * parser-shape while every runtime input arrives in evaluated-shape;
     * see the long comment above pattern_canonicalize for the full
     * rationale. We retain ownership of `canon`. */
    Expr* canon = pattern_canonicalize(expr_copy(pattern));

    // Check if rule with same pattern already exists
    for (Rule* curr = *list; curr; curr = curr->next) {
        if (expr_eq(curr->pattern, canon)) {
            // Replace replacement in place; sort key is unchanged.
            expr_free(curr->replacement);
            curr->replacement = expr_copy(replacement);
            expr_free(canon);
            return;
        }
    }

    Rule* new_rule = malloc(sizeof(Rule));
    new_rule->pattern = canon;
    new_rule->replacement = expr_copy(replacement);
    new_rule->dispatch_arity = pattern_dispatch_arity(canon);
    new_rule->first_arg_head_canon =
        (new_rule->dispatch_arity == -1) ? NULL : pattern_first_arg_head(canon);
    new_rule->specificity = pattern_specificity(canon);
    new_rule->next = NULL;

    /* §3.6 stable insertion: walk until we find a rule with strictly
     * lower specificity, then insert before it. Equal-specificity rules
     * keep insertion order (earlier additions stay earlier). */
    Rule* prev = NULL;
    Rule* curr = *list;
    while (curr && curr->specificity >= new_rule->specificity) {
        prev = curr;
        curr = curr->next;
    }
    new_rule->next = curr;
    if (prev) {
        prev->next = new_rule;
    } else {
        *list = new_rule;
    }
}
void symtab_add_own_value(const char* symbol_name, Expr* pattern, Expr* replacement) {
    SymbolDef* def = symtab_get_def(symbol_name);
    add_rule(&def->own_values, pattern, replacement);
}

void symtab_add_down_value(const char* symbol_name, Expr* pattern, Expr* replacement) {
    SymbolDef* def = symtab_get_def(symbol_name);
    add_rule(&def->down_values, pattern, replacement);
}

void symtab_clear_symbol(const char* symbol_name) {
    /* Removing rules invalidates cached evaluations that might have
     * relied on them; bump the eval clock unconditionally. */
    eval_clock_bump();

    SymbolDef* def = symtab_get_def(symbol_name);
    Rule* curr = def->own_values;
    while (curr) {
        Rule* next = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        curr = next;
    }
    def->own_values = NULL;

    curr = def->down_values;
    while (curr) {
        Rule* next = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        curr = next;
    }
    def->down_values = NULL;
}

void symtab_clear(void) {
    for (int i = 0; i < SYMTAB_SIZE; i++) {
        SymEntry* entry = symtab[i];
        while (entry) {
            SymEntry* next = entry->next;
            
            // Inline symtab_clear_symbol logic without lookup
            Rule* curr = entry->def->own_values;
            while (curr) {
                Rule* r_next = curr->next;
                expr_free(curr->pattern);
                expr_free(curr->replacement);
                free(curr);
                curr = r_next;
            }
            entry->def->own_values = NULL;

            curr = entry->def->down_values;
            while (curr) {
                Rule* r_next = curr->next;
                expr_free(curr->pattern);
                expr_free(curr->replacement);
                free(curr);
                curr = r_next;
            }
            entry->def->down_values = NULL;

            if (entry->def->docstring) free(entry->def->docstring);
            /* symbol_name is interned and owned by sym_intern; do not free. */
            free(entry->def);
            free(entry);
            entry = next;
        }
        symtab[i] = NULL;
    }
}


Rule* symtab_get_own_values(const char* symbol_name) {
    const char* canon = intern_symbol(symbol_name);
    unsigned int idx = hash(canon);
    SymEntry* entry = symtab[idx];
    while (entry) {
        if (entry->def->symbol_name == canon) {
            return entry->def->own_values;
        }
        entry = entry->next;
    }
    return NULL;
}

Rule* symtab_get_down_values(const char* symbol_name) {
    const char* canon = intern_symbol(symbol_name);
    unsigned int idx = hash(canon);
    SymEntry* entry = symtab[idx];
    while (entry) {
        if (entry->def->symbol_name == canon) {
            return entry->def->down_values;
        }
        entry = entry->next;
    }
    return NULL;
}

/* §3.5 input-side dispatch key. Returns the canonical head of an
 * input-call's first argument, or NULL when we can't safely
 * specialize on it. The classification mirrors pattern_arg_head_canon
 * so a head set on a rule's pattern will line up with the head we
 * compute here for an input expression. */
static const char* input_arg_head_canon(const Expr* arg) {
    if (!arg) return NULL;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_BIGINT) {
        return intern_symbol("Integer");
    }
    if (arg->type == EXPR_REAL) return intern_symbol("Real");
    if (arg->type == EXPR_STRING) return intern_symbol("String");
    if (arg->type == EXPR_SYMBOL) return intern_symbol("Symbol");
    if (arg->type == EXPR_FUNCTION && arg->data.function.head &&
        arg->data.function.head->type == EXPR_SYMBOL) {
        return intern_symbol(arg->data.function.head->data.symbol);
    }
    return NULL;
}

Expr* apply_down_values(Expr* expr) {
    if (!expr || expr->type != EXPR_FUNCTION) return NULL;

    Expr* head = expr->data.function.head;
    if (head->type != EXPR_SYMBOL) return NULL;

    SymbolDef* def = symtab_get_def(head->data.symbol);
    if (!def->down_values) return NULL;

    /* §3.5 dispatch filter: pre-compute the input's shape once, then
     * compare against each rule's cached key before invoking the
     * (much more expensive) matcher. A NULL key on either side means
     * "wildcard" -- we always scan in that case. */
    const int32_t input_arity = (int32_t)expr->data.function.arg_count;
    const char* input_first_head =
        (input_arity > 0)
            ? input_arg_head_canon(expr->data.function.args[0])
            : NULL;

    for (Rule* rule = def->down_values; rule; rule = rule->next) {
        /* Arity filter: only skip when both sides committed to a
         * specific arity that disagrees. A rule with -1 arity (variable)
         * always passes; an input never has variable arity. */
        if (rule->dispatch_arity != -1 && rule->dispatch_arity != input_arity) {
            continue;
        }
        /* First-arg head filter: only skip when the rule expects a
         * specific head and we determined the input's first-arg head and
         * the two disagree. Pointer compare is sound because both are
         * interned. */
        if (rule->first_arg_head_canon && input_first_head &&
            rule->first_arg_head_canon != input_first_head) {
            continue;
        }

        MatchEnv* env = env_new();
        if (match(expr, rule->pattern, env)) {
            Expr* result = replace_bindings(rule->replacement, env);
            env_free(env);
            return result;
        }
        env_free(env);
    }
    return NULL;
}

Expr* apply_own_values(Expr* expr) {
    if (!expr || expr->type != EXPR_SYMBOL) return NULL;
    
    SymbolDef* def = symtab_get_def(expr->data.symbol);
    Rule* rule = def->own_values;
    while (rule) {
        MatchEnv* env = env_new();
        if (match(expr, rule->pattern, env)) {
            Expr* result = replace_bindings(rule->replacement, env);
            env_free(env);
            return result;
        }
        env_free(env);
        rule = rule->next;
    }
    return NULL;
}
