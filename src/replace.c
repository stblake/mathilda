#include "replace.h"
#include "symtab.h"
#include "eval.h"
#include "match.h"
#include "sym_names.h"
#include "assoc.h"
#include "assoc_struct.h"   /* associations are atomic: rules reach values only */
#include "common.h"         /* head_is */
#include "part.h"
#include "print.h"          /* expr_to_string */
#include "message.h"        /* mth_message */
#include <stdlib.h>
#include <string.h>

typedef struct {
    Expr* pos_spec;
    Expr* replacement;
    bool delayed;
} RawRule;

static void add_raw_rule(RawRule** rules, size_t* cap, size_t* count, Expr* path, Expr* replacement, bool delayed) {
    if (*count >= *cap) {
        *cap = (*cap == 0) ? 4 : (*cap * 2);
        *rules = realloc(*rules, sizeof(RawRule) * (*cap));
    }
    (*rules)[*count].pos_spec = expr_copy(path);
    (*rules)[*count].replacement = expr_copy(replacement);
    (*rules)[*count].delayed = delayed;
    (*count)++;
}

static void parse_single_rule(Expr* rule_expr, RawRule** rules, size_t* cap, size_t* count) {
    if (rule_expr->type != EXPR_FUNCTION) return;
    bool delayed = false;
    if (rule_expr->data.function.head->data.symbol.name == SYM_Rule) delayed = false;
    else if (rule_expr->data.function.head->data.symbol.name == SYM_RuleDelayed) delayed = true;
    else return;

    if (rule_expr->data.function.arg_count != 2) return;
    Expr* lhs = rule_expr->data.function.args[0];
    Expr* rhs = rule_expr->data.function.args[1];

    if (lhs->type == EXPR_FUNCTION && lhs->data.function.head->data.symbol.name == SYM_List) {
        if (lhs->data.function.arg_count > 0 && 
            lhs->data.function.args[0]->type == EXPR_FUNCTION &&
            lhs->data.function.args[0]->data.function.head->data.symbol.name == SYM_List) {
            for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
                add_raw_rule(rules, cap, count, lhs->data.function.args[i], rhs, delayed);
            }
        } else {
            add_raw_rule(rules, cap, count, lhs, rhs, delayed);
        }
    } else {
        Expr* lhs_copy = expr_copy(lhs);
        Expr* wrapped = expr_new_function(expr_new_symbol(SYM_List), &lhs_copy, 1);
        add_raw_rule(rules, cap, count, wrapped, rhs, delayed);
        expr_free(wrapped);
    }
}

typedef struct {
    int64_t pos_index;
    int64_t neg_index;
} PathElement;

/*
 * normalize_pos_spec:
 * Convert negative integer indices in a rule's position spec to their positive
 * equivalents using the known path lengths. Returns a new expression that must
 * be freed, or NULL if no normalization was needed.
 */
static Expr* normalize_pos_spec(Expr* pos_spec, PathElement* path, size_t depth) {
    if (pos_spec->type != EXPR_FUNCTION || depth == 0) return NULL;
    if (pos_spec->data.function.arg_count != depth) return NULL;

    bool needs_norm = false;
    for (size_t d = 0; d < depth; d++) {
        Expr* elem = pos_spec->data.function.args[d];
        if (elem->type == EXPR_INTEGER && elem->data.integer < 0) {
            needs_norm = true;
            break;
        }
    }
    if (!needs_norm) return NULL;

    Expr** new_args = malloc(sizeof(Expr*) * depth);
    for (size_t d = 0; d < depth; d++) {
        Expr* elem = pos_spec->data.function.args[d];
        if (elem->type == EXPR_INTEGER && elem->data.integer < 0) {
            /* Recover length from pos_index and neg_index:
             * neg_index = pos_index - length - 1, so length = pos_index - neg_index - 1 */
            int64_t length = path[d].pos_index - path[d].neg_index - 1;
            int64_t positive = elem->data.integer + length + 1;
            new_args[d] = expr_new_integer(positive);
        } else {
            new_args[d] = expr_copy(elem);
        }
    }
    Expr* result = expr_new_function(expr_copy(pos_spec->data.function.head), new_args, depth);
    free(new_args);
    return result;
}

static bool check_match(RawRule* rules, size_t num_rules, PathElement* path, size_t depth,
                        Expr** out_replacement, bool* out_delayed, MatchEnv** out_env) {
    for (int64_t r = (int64_t)num_rules - 1; r >= 0; r--) {
        /* Build the current position as a List of positive indices */
        Expr** alias_args = NULL;
        if (depth > 0) alias_args = malloc(sizeof(Expr*) * depth);
        for (size_t d = 0; d < depth; d++) {
            alias_args[d] = expr_new_integer(path[d].pos_index);
        }
        Expr* alias_expr = expr_new_function(expr_new_symbol(SYM_List), alias_args, depth);

        /* Normalize any negative indices in the rule's pos_spec to positive */
        Expr* norm_spec = normalize_pos_spec(rules[r].pos_spec, path, depth);
        Expr* spec_to_match = norm_spec ? norm_spec : rules[r].pos_spec;

        MatchEnv* env = env_new();
        if (match(alias_expr, spec_to_match, env)) {
            *out_replacement = rules[r].replacement;
            *out_delayed = rules[r].delayed;
            *out_env = env;
            expr_free(alias_expr);
            if (alias_args) free(alias_args);
            if (norm_spec) expr_free(norm_spec);
            return true;
        }
        env_free(env);
        expr_free(alias_expr);
        if (alias_args) free(alias_args);
        if (norm_spec) expr_free(norm_spec);
    }
    return false;
}

/* A rule's replacement is inserted with its outer Unevaluated layers removed,
 * for Rule and RuleDelayed alike (a Rule keeps an
 * Unevaluated RHS as an argument): x /. d -> Unevaluated[Sequence[]] deletes
 * d (the bare Sequence[] then splices, even into Power[d, 2] -> Power[2] ->
 * 2), and Hold[d] /. d :> Unevaluated[1 + 1] is Hold[1 + 1]. Consumes `e`. */
static Expr* strip_unevaluated(Expr* e) {
    while (e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Unevaluated) {
        Expr* inner;
        if (e->data.function.arg_count == 1) {
            inner = expr_copy(e->data.function.args[0]);
        } else {
            /* Binding already spliced a Sequence into the wrapper:
             * Unevaluated[Sequence[]] arrives as Unevaluated[], and
             * Unevaluated[Sequence[a, b]] as Unevaluated[a, b]. */
            size_t n = e->data.function.arg_count;
            Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
            for (size_t i = 0; i < n; i++) a[i] = expr_copy(e->data.function.args[i]);
            inner = expr_new_function(expr_new_symbol(SYM_Sequence), a, n);
            free(a);
        }
        expr_free(e);
        e = inner;
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
            e->data.function.head->data.symbol.name == SYM_Sequence) break;
    }
    return e;
}

static Expr* replace_part_rec(Expr* expr, RawRule* rules, size_t num_rules, PathElement* current_path, size_t depth) {
    Expr* replacement_expr = NULL;
    bool delayed = false;
    MatchEnv* env = NULL;

    if (check_match(rules, num_rules, current_path, depth, &replacement_expr, &delayed, &env)) {
        Expr* res = replace_bindings(replacement_expr, env);
        if (delayed) {
            Expr* eval_res = evaluate(res);
            expr_free(res);
            res = eval_res;
        }
        res = strip_unevaluated(res);
        env_free(env);
        return res;
    }

    if (expr->type == EXPR_FUNCTION) {
        size_t len = expr->data.function.arg_count;
        Expr** new_args = NULL;
        if (len > 0) new_args = malloc(sizeof(Expr*) * len);

        PathElement* next_path = malloc(sizeof(PathElement) * (depth + 1));
        if (depth > 0) memcpy(next_path, current_path, sizeof(PathElement) * depth);

        /* Check head (index 0) -- only replace if a rule explicitly targets index 0
         * with a literal integer. Default Heads -> False means patterns like _ or
         * Except[...] should not match the head position. */
        next_path[depth].pos_index = 0;
        next_path[depth].neg_index = 0;
        Expr* new_head = NULL;
        bool head_replaced = false;
        for (size_t r = num_rules; r > 0; r--) {
            Expr* spec = rules[r - 1].pos_spec;
            /* Only consider this rule for the head if it has a literal 0 at this depth */
            if (spec->type == EXPR_FUNCTION && spec->data.function.arg_count == depth + 1) {
                Expr* idx = spec->data.function.args[depth];
                if (idx->type == EXPR_INTEGER && idx->data.integer == 0) {
                    new_head = replace_part_rec(expr->data.function.head, rules, num_rules, next_path, depth + 1);
                    head_replaced = true;
                    break;
                }
            }
        }
        if (!head_replaced) {
            new_head = expr_copy(expr->data.function.head);
        }

        for (size_t i = 0; i < len; i++) {
            next_path[depth].pos_index = i + 1;
            next_path[depth].neg_index = (int64_t)(i + 1) - (int64_t)len - 1;
            new_args[i] = replace_part_rec(expr->data.function.args[i], rules, num_rules, next_path, depth + 1);
        }
        
        free(next_path);
        
        Expr* result = expr_new_function(new_head, new_args, len);
        if (new_args) free(new_args);
        return result;
    }
    
    return expr_copy(expr);
}

/* ---------------------------------------------------------------------------
 * Legacy positional forms (used by old packages, e.g. DifferentialForms'
 * Indices):
 *
 *   ReplacePart[expr, new, pos]         parts at pos become new
 *   ReplacePart[expr, new, pos, npos]   parts at pos become new[[npos]]
 *
 * pos is an integer n, a single position {i, j, ...}, or a list of positions
 * {{...}, {...}} (so {} names NO position and leaves expr alone, while {{}}
 * names expr itself). Unlike the rule form, a position that does not exist is
 * an error: ReplacePart::partw, and the call stays unevaluated. The forms are
 * translated into one `position -> value` RawRule per position and handed to
 * the rule engine (replace_part_rec), so the two surfaces cannot drift.
 * ------------------------------------------------------------------------- */

static bool rp_is_int_list(const Expr* e) {
    if (!head_is(e, SYM_List)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (e->data.function.args[i]->type != EXPR_INTEGER) return false;
    return true;
}

/* Normalise a legacy position spec into a List of positions (each a List of
 * integers). Returns an owned List, or NULL when the spec is malformed. */
static Expr* rp_positions(Expr* spec) {
    if (spec->type == EXPR_INTEGER) {
        Expr* one = expr_copy(spec);
        Expr* pos = expr_new_function(expr_new_symbol(SYM_List), &one, 1);
        return expr_new_function(expr_new_symbol(SYM_List), &pos, 1);
    }
    if (!head_is(spec, SYM_List)) return NULL;
    size_t n = spec->data.function.arg_count;
    if (n > 0 && rp_is_int_list(spec)) {          /* one position {i, j, ...} */
        Expr* pos = expr_copy(spec);
        return expr_new_function(expr_new_symbol(SYM_List), &pos, 1);
    }
    for (size_t i = 0; i < n; i++)                /* a list of positions */
        if (!rp_is_int_list(spec->data.function.args[i])) return NULL;
    return expr_copy(spec);
}

/* The part of `e` at integer position `pos` (borrowed), or NULL if it does not
 * exist. Index 0 is the head; negative indices count from the end. */
static Expr* rp_part_at(Expr* e, const Expr* pos) {
    for (size_t d = 0; d < pos->data.function.arg_count; d++) {
        if (e->type != EXPR_FUNCTION) return NULL;
        int64_t k = pos->data.function.args[d]->data.integer;
        int64_t len = (int64_t)e->data.function.arg_count;
        if (k == 0) { e = e->data.function.head; continue; }
        if (k < 0) k += len + 1;
        if (k < 1 || k > len) return NULL;
        e = e->data.function.args[k - 1];
    }
    return e;
}

static void rp_partw(Expr* pos, Expr* e) {
    char* ps = expr_to_string(pos);
    char* es = expr_to_string(e);
    mth_message("ReplacePart", "partw", "Part %s of %s does not exist.", ps, es);
    free(ps);
    free(es);
}

static Expr* replace_part_legacy(Expr* res) {
    size_t argc = res->data.function.arg_count;
    Expr* expr = res->data.function.args[0];
    Expr* new_e = res->data.function.args[1];
    Expr* poss = rp_positions(res->data.function.args[2]);
    Expr* nposs = NULL;
    if (!poss || (argc == 4 && !(nposs = rp_positions(res->data.function.args[3])))) {
        Expr* bad = poss ? res->data.function.args[3] : res->data.function.args[2];
        char* bs = expr_to_string(bad);
        char* rs = expr_to_string(res);
        mth_message("ReplacePart", "psl",
                    "Position specification %s in %s is not a machine-sized "
                    "integer or a list of machine-sized integers.", bs, rs);
        free(bs);
        free(rs);
        if (poss) expr_free(poss);
        return NULL;
    }
    size_t np = poss->data.function.arg_count;
    /* npos pairs with pos position by position; a single npos serves all. */
    if (nposs && nposs->data.function.arg_count != np &&
        nposs->data.function.arg_count != 1) {
        expr_free(poss);
        expr_free(nposs);
        return NULL;
    }

    RawRule* rules = malloc(sizeof(RawRule) * (np ? np : 1));
    size_t nr = 0;
    bool ok = true;
    for (size_t i = 0; i < np && ok; i++) {
        Expr* pos = poss->data.function.args[i];
        if (!rp_part_at(expr, pos)) { rp_partw(pos, expr); ok = false; break; }
        Expr* val = new_e;
        if (nposs) {
            Expr* npos = nposs->data.function.args[nposs->data.function.arg_count == 1 ? 0 : i];
            val = rp_part_at(new_e, npos);
            if (!val) { rp_partw(npos, new_e); ok = false; break; }
        }
        rules[nr].pos_spec = expr_copy(pos);
        rules[nr].replacement = expr_copy(val);
        rules[nr].delayed = false;
        nr++;
    }

    Expr* result = NULL;
    if (ok) result = nr ? replace_part_rec(expr, rules, nr, NULL, 0) : expr_copy(expr);
    for (size_t i = 0; i < nr; i++) {
        expr_free(rules[i].pos_spec);
        expr_free(rules[i].replacement);
    }
    free(rules);
    expr_free(poss);
    if (nposs) expr_free(nposs);
    return result;
}

Expr* builtin_replace_part(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count == 3 || res->data.function.arg_count == 4)
        return replace_part_legacy(res);
    if (res->data.function.arg_count == 1) {
        // Operator form ReplacePart[rules]
        return NULL; // Stays unevaluated
    } else if (res->data.function.arg_count == 2) {
        Expr* expr = res->data.function.args[0];
        Expr* rules_expr = res->data.function.args[1];

        size_t cap = 4;
        size_t num_rules = 0;
        RawRule* rules = malloc(sizeof(RawRule) * cap);

        if (rules_expr->type == EXPR_FUNCTION && rules_expr->data.function.head->data.symbol.name == SYM_List) {
            for (size_t i = 0; i < rules_expr->data.function.arg_count; i++) {
                parse_single_rule(rules_expr->data.function.args[i], &rules, &cap, &num_rules);
            }
        } else {
            parse_single_rule(rules_expr, &rules, &cap, &num_rules);
        }

        if (num_rules == 0) {
            free(rules);
            return expr_copy(expr);
        }

        /* Association: replace values by key position (single-level Key[k] /
         * "str" / positional i, or wrapped in a one-element list). Replace-only,
         * matching Wolfram (absent keys are left unchanged, not added). */
        if (is_association(expr)) {
            /* Rebuild over a private array. expr_copy is a refcount bump, so
             * `result` would BE the caller's association and rewriting an
             * entry's value in place used to corrupt it:
             *   w = <|"x" -> 1|>; ReplacePart[w, "x" -> 9]  mutated w. */
            size_t len = expr->data.function.arg_count;
            Expr** new_args = NULL;
            if (len > 0) {
                new_args = malloc(sizeof(Expr*) * len);
                for (size_t i = 0; i < len; i++)
                    new_args[i] = expr_copy(expr->data.function.args[i]);
            }
            for (size_t r = 0; r < num_rules; r++) {
                Expr* s = rules[r].pos_spec;
                if (s->type == EXPR_FUNCTION && s->data.function.head->type == EXPR_SYMBOL &&
                    s->data.function.head->data.symbol.name == SYM_List && s->data.function.arg_count == 1)
                    s = s->data.function.args[0];   /* unwrap {X} */
                Expr* key = NULL; bool positional = false; int64_t p = 0;
                if (s->type == EXPR_FUNCTION && s->data.function.head->type == EXPR_SYMBOL &&
                    s->data.function.head->data.symbol.name == SYM_Key && s->data.function.arg_count == 1) {
                    key = s->data.function.args[0];
                } else if (s->type == EXPR_INTEGER) {
                    positional = true; p = s->data.integer;
                } else if (!(s->type == EXPR_FUNCTION && s->data.function.head->type == EXPR_SYMBOL &&
                             s->data.function.head->data.symbol.name == SYM_List)) {
                    key = s;   /* a literal key */
                } else {
                    continue;  /* multi-element/nested spec: unsupported here */
                }
                int64_t t = -1;
                if (positional) {
                    if (p < 0) p = (int64_t)len + p + 1;
                    if (p >= 1 && p <= (int64_t)len) t = p - 1;
                } else {
                    for (size_t i = 0; i < len; i++) {
                        Expr* rule = new_args[i];
                        if (is_rule2(rule) &&
                            expr_eq(rule->data.function.args[0], key)) { t = (int64_t)i; break; }
                    }
                }
                if (t < 0) continue;   /* replace-only: leave absent keys */
                Expr* rule = new_args[t];
                if (!is_rule2(rule)) continue;
                Expr* val = rules[r].delayed ? evaluate(rules[r].replacement)
                                             : expr_copy(rules[r].replacement);
                new_args[t] = assoc_entry_with_value(rule, val);   /* adopts val */
                expr_free(rule);
            }
            for (size_t i = 0; i < num_rules; i++) {
                expr_free(rules[i].pos_spec);
                expr_free(rules[i].replacement);
            }
            free(rules);
            Expr* result = expr_new_function(expr_copy(expr->data.function.head),
                                             new_args, len);
            free(new_args);
            return result;
        }

        Expr* result = replace_part_rec(expr, rules, num_rules, NULL, 0);

        for (size_t i = 0; i < num_rules; i++) {
            expr_free(rules[i].pos_spec);
            expr_free(rules[i].replacement);
        }
        free(rules);

        return result;
    }
    return NULL;
}

/* One rule of a parsed rule set. `pattern == NULL` marks an ASSOCIATION used
 * as a rule set (ReplaceAll[{"a", "b"}, <|"a" -> 1|>] is {1, "b"}), and then
 * `replacement` is that association: its keys match LITERALLY, never as
 * patterns (x /. <|x_ -> 1|> is x in Mathematica), so a whole association is
 * one entry answered by an O(1) keyed lookup instead of n pattern matches. */
typedef struct {
    Expr* pattern;
    Expr* replacement;
    bool delayed;
    Expr* owned;     /* a hoisted pattern lhs /; t (owned), else NULL */
} ReplaceRule;

/* Append the rule(s) `r` to a growable rule array: a two-argument Rule or
 * RuleDelayed, or a well-formed association (one keyed entry, see above).
 * Anything else is ignored. Borrows from `r`. */
static void rules_push(ReplaceRule** rules, size_t* cap, size_t* count, Expr* r) {
    bool dict = assoc_is_wellformed(r);
    if (!dict && !is_rule2(r)) return;
    if (*count == *cap) {
        *cap = *cap ? *cap * 2 : 4;
        *rules = realloc(*rules, sizeof(ReplaceRule) * (*cap));
    }
    ReplaceRule* slot = &(*rules)[(*count)++];
    slot->owned = NULL;
    if (dict) {
        slot->pattern = NULL;
        slot->replacement = r;
        slot->delayed = false;   /* a RuleDelayed value is returned held */
    } else {
        slot->pattern = r->data.function.args[0];
        slot->replacement = r->data.function.args[1];
        slot->delayed = r->data.function.head->data.symbol.name == SYM_RuleDelayed;
        /* lhs :> body /; test is the rule lhs /; test :> body: the guard is
         * part of the match, so a failing guard backtracks into the pattern's
         * other assignments (which pair of a Plus) and then lets the next rule
         * try. A Rule (->) keeps a Condition value. */
        if (slot->delayed) {
            Expr* body = slot->replacement;
            Expr* pat = eval_hoist_rhs_conditions(slot->pattern, slot->replacement, &body);
            if (pat != slot->pattern) {
                slot->owned = pat;
                slot->pattern = pat;
                slot->replacement = body;
            }
        }
    }
}

/* Release a rule array from rules_collect, with any hoisted patterns. */
static void rules_free(ReplaceRule* rules, size_t count) {
    for (size_t i = 0; i < count; i++)
        if (rules[i].owned) expr_free(rules[i].owned);
    free(rules);
}

/* A delayed rule's bound replacement, honouring a scoped guard
 * (lhs :> Module[{t = ...}, val /; cond]): its value, or NULL when the guard
 * fails and the rule does not apply. Consumes `bound`. */
static Expr* rule_finish(const ReplaceRule* r, Expr* bound) {
    if (!bound) return NULL;
    if (!r->delayed || !symtab_rhs_has_scoped_condition(r->replacement))
        return strip_unevaluated(bound);
    Expr* v = symtab_try_scoped_condition(bound);
    expr_free(bound);
    return v;
}

/* Collect the rules of `rules_expr` (a rule, an association, or a List of
 * them) into a fresh array the caller frees; *out_count receives the count. */
static ReplaceRule* rules_collect(Expr* rules_expr, size_t* out_count) {
    size_t cap = 4, count = 0;
    ReplaceRule* rules = malloc(sizeof(ReplaceRule) * cap);
    if (head_is(rules_expr, SYM_List)) {
        for (size_t i = 0; i < rules_expr->data.function.arg_count; i++)
            rules_push(&rules, &cap, &count, rules_expr->data.function.args[i]);
    } else {
        rules_push(&rules, &cap, &count, rules_expr);
    }
    *out_count = count;
    return rules;
}

/* Try one rule at `e`: the owned, bound (but unevaluated) replacement, or NULL
 * when the rule does not apply. */
static Expr* rule_try(const ReplaceRule* r, Expr* e) {
    if (!r->pattern) {
        Expr* v = assoc_lookup_value(r->replacement, e);
        return v ? expr_copy(v) : NULL;
    }
    /* Flat-aware: a Flat-head (Plus/Times/...) rule rewrites a SUBSET of a
     * longer operand sequence, keeping the rest (e.g. 1+x+y /. x+y->z gives
     * 1+z). The ordinary match stays inline (so the common path is unchanged);
     * only a failed match falls through to Flat-head leftover matching, which
     * adds genuinely-new matches without altering existing ones. */
    MatchEnv* env = env_new();
    Expr* out = match(e, r->pattern, env)
        ? replace_bindings(r->replacement, env)
        : match_flat_leftover(e, r->pattern, r->replacement, env);
    env_free(env);
    return rule_finish(r, out);
}

/* A Rational or Complex number is an atom to every replacement traversal: its
 * parts and head are never rewritten, as in the reference -- {1/2} /. 2 -> 3
 * stays {1/2} (it was {1/3}), {1/2} /. 1 -> 5 stays {1/2}, and
 * {I} /. Complex -> h stays {I}.  A rule for the number as a whole
 * (Rational[a_, b_] -> g[a, b]) still matches it. */
static bool rp_number_atom(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rational ||
            e->data.function.head->data.symbol.name == SYM_Complex);
}

static int64_t get_expr_depth_replace(Expr* e, bool heads) {
    if (e->type != EXPR_FUNCTION) return 1;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Rational || h == SYM_Complex) return 1;
    }
    bool assoc = assoc_is_wellformed(e);     /* parts = values */
    int64_t max_d = assoc ? 1 : 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int64_t d = get_expr_depth_replace(struct_part(e, i, assoc), heads);
        if (d > max_d) max_d = d;
    }
    if (heads) {
        int64_t d_head = get_expr_depth_replace(e->data.function.head, heads);
        if (d_head > max_d) max_d = d_head;
    }
    return max_d + 1;
}

static Expr* do_replace_at_level(Expr* e, int64_t current_level, int64_t min_l, int64_t max_l, bool heads, ReplaceRule* rules, size_t num_rules) {
    if (max_l >= 0 && current_level > max_l) return expr_copy(e);

    int64_t d = get_expr_depth_replace(e, heads);

    bool match_level = true;
    if (min_l >= 0) {
        if (current_level < min_l) match_level = false;
    } else {
        if (d > -min_l) match_level = false;
    }

    if (max_l >= 0) {
        if (current_level > max_l) match_level = false;
    } else {
        if (d < -max_l) match_level = false;
    }

    Expr* new_e = NULL;
    if (e->type == EXPR_FUNCTION && !rp_number_atom(e)) {
        Expr* new_head = e->data.function.head;
        if (heads) {
            new_head = do_replace_at_level(e->data.function.head, current_level + 1, min_l, max_l, heads, rules, num_rules);
        } else {
            new_head = expr_copy(new_head);
        }

        /* An atomic association's level-1 parts are its values; keys are
         * never offered to the rules (assoc_struct.h). */
        bool assoc = assoc_is_wellformed(e);
        size_t count = e->data.function.arg_count;
        Expr** new_args = NULL;
        if (count > 0) new_args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            new_args[i] = do_replace_at_level(struct_part(e, i, assoc), current_level + 1, min_l, max_l, heads, rules, num_rules);
        }

        new_e = struct_rebuild(e, assoc, new_head, new_args, count, true);
        if (new_args) free(new_args);
    } else {
        new_e = expr_copy(e);
    }

    if (match_level) {
        for (size_t i = 0; i < num_rules; i++) {
            Expr* repl = rule_try(&rules[i], new_e);
            if (repl) {
                expr_free(new_e);
                return repl;
            }
        }
    }

    return new_e;
}

static Expr* dispatch_unwrap(Expr* rules_expr, bool* indexed);

static bool is_rule(Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    return (h == SYM_Rule || h == SYM_RuleDelayed);
}

static Expr* apply_replace_nested(Expr* expr, Expr* rules_expr, int64_t min_l, int64_t max_l, bool heads) {
    if (rules_expr->type == EXPR_FUNCTION && rules_expr->data.function.head->data.symbol.name == SYM_List) {
        if (rules_expr->data.function.arg_count > 0 && !is_rule(rules_expr->data.function.args[0])) {
            size_t count = rules_expr->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = apply_replace_nested(expr, rules_expr->data.function.args[i], min_l, max_l, heads);
            }
            Expr* res = expr_new_function(expr_new_symbol(SYM_List), args, count);
            free(args);
            return res;
        }
    }
    
    size_t num_rules = 0;
    ReplaceRule* rules = rules_collect(rules_expr, &num_rules);

    Expr* result;
    if (num_rules > 0) {
        result = do_replace_at_level(expr, 0, min_l, max_l, heads, rules, num_rules);
    } else {
        result = expr_copy(expr);
    }
    rules_free(rules, num_rules);
    return result;
}

Expr* builtin_replace(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* rules_expr = res->data.function.args[1];
    
    int64_t min_l = 0, max_l = 0;
    bool heads = false;
    
    if (res->data.function.arg_count >= 3) {
        Expr* ls = res->data.function.args[2];
        if (ls->type == EXPR_INTEGER) {
            min_l = 1; max_l = ls->data.integer;
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_All) {
            min_l = 0; max_l = 1000000;
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_Infinity) {
            min_l = 1; max_l = 1000000;
        } else if (ls->type == EXPR_FUNCTION && ls->data.function.head->data.symbol.name == SYM_List) {
            if (ls->data.function.arg_count == 1 && ls->data.function.args[0]->type == EXPR_INTEGER) {
                min_l = max_l = ls->data.function.args[0]->data.integer;
            } else if (ls->data.function.arg_count == 2) {
                if (ls->data.function.args[0]->type == EXPR_INTEGER) min_l = ls->data.function.args[0]->data.integer;
                if (ls->data.function.args[1]->type == EXPR_INTEGER) max_l = ls->data.function.args[1]->data.integer;
                else if (ls->data.function.args[1]->type == EXPR_SYMBOL && ls->data.function.args[1]->data.symbol.name == SYM_Infinity) max_l = 1000000;
            }
        }
    }
    
    for (size_t i = 2; i < res->data.function.arg_count; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol.name == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol.name == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_True) heads = true;
                else if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_False) heads = false;
            }
        }
    }
    
    return apply_replace_nested(expr, dispatch_unwrap(rules_expr, NULL), min_l, max_l, heads);
}

static Expr* do_replace_all(Expr* e, ReplaceRule* rules, size_t num_rules) {
    for (size_t i = 0; i < num_rules; i++) {
        Expr* repl = rule_try(&rules[i], e);
        if (repl) return repl;
    }

    if (e->type == EXPR_FUNCTION && !rp_number_atom(e)) {
        Expr* new_head = do_replace_all(e->data.function.head, rules, num_rules);
        /* Inside an atomic association the rules reach the values only: the
         * keys and the k -> v entries are never matched, so
         * <|a -> 1, b -> 2|> /. b -> a keeps both keys (assoc_struct.h). The
         * head is still tried, as in Mathematica: <|a -> 1|> /. Association
         * -> g is g[a -> 1]. */
        bool assoc = assoc_is_wellformed(e);
        size_t count = e->data.function.arg_count;
        Expr** new_args = NULL;
        if (count > 0) new_args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            new_args[i] = do_replace_all(struct_part(e, i, assoc), rules, num_rules);
        }
        Expr* new_e = struct_rebuild(e, assoc, new_head, new_args, count, true);
        if (new_args) free(new_args);
        return new_e;
    }

    return expr_copy(e);
}

/* ---- Dispatch: a rule set hashed by left-hand side --------------------------
 *
 * ReplaceAll normally tries every rule at every node, O(nodes x rules).  When
 * the rules come wrapped in Dispatch[...] -- the user's explicit request for a
 * dispatch table -- each rule is keyed by what its LHS can possibly match:
 *   ATOM  a literal atom (x, 1, "s"): matches only that atom;
 *   HEAD  a literal symbol head f[...]: matches only an f[...] node.
 * HoldPattern / Condition / PatternTest / Pattern wrappers are looked through
 * (the guard still runs) and Verbatim[x] keys on x.  Everything else -- a
 * pattern-construct LHS (x_, _h, a | b, ...), a compound head, an association
 * entry, or an LHS with an Optional / x_. argument (x^n_. matches a bare x,
 * which has no Power head) -- is WILD and tried at every node.  At a node the
 * candidates are its bucket merged with the wild rules in ORIGINAL rule order,
 * so the first rule that matches is exactly the one the linear scan finds. */

enum { RK_WILD = 0, RK_ATOM = 1, RK_HEAD = 2 };

typedef struct {
    Expr* key;          /* borrowed: an atom, or a head symbol */
    int kind;
    size_t* idx;        /* rule indices, ascending */
    size_t n, cap;
} RuleBucket;

typedef struct {
    RuleBucket* tab;    /* open addressing; key == NULL marks an empty slot */
    size_t cap;         /* power of two */
    size_t* wild;       /* ascending */
    size_t nwild;
} RuleIndex;

static bool rk_pattern_head(const char* h) {
    return h == SYM_Blank || h == SYM_BlankSequence || h == SYM_BlankNullSequence ||
           h == SYM_Pattern || h == SYM_Alternatives || h == SYM_Except ||
           h == SYM_PatternTest || h == SYM_Condition || h == SYM_Optional ||
           h == SYM_Repeated || h == SYM_RepeatedNull || SYM_IS_HOLDPATTERN(h) ||
           h == SYM_Verbatim || h == SYM_Longest || h == SYM_Shortest ||
           h == SYM_OptionsPattern || h == SYM_PatternSequence ||
           h == SYM_OrderlessPatternSequence || h == SYM_KeyValuePattern;
}

/* Strip the wrappers that do not change what a pattern can match. */
static Expr* rk_unwrap(Expr* p) {
    for (;;) {
        if (p->type != EXPR_FUNCTION || p->data.function.head->type != EXPR_SYMBOL) return p;
        const char* h = p->data.function.head->data.symbol.name;
        size_t n = p->data.function.arg_count;
        if (SYM_IS_HOLDPATTERN(h) && n == 1)                     p = p->data.function.args[0];
        else if ((h == SYM_Condition || h == SYM_PatternTest) && n == 2) p = p->data.function.args[0];
        else if (h == SYM_Pattern && n == 2)                    p = p->data.function.args[1];
        else return p;
    }
}

/* The key of a rule's LHS (see above); RK_WILD when it has none. */
static int rk_rule_key(const ReplaceRule* r, Expr** key) {
    if (!r->pattern) return RK_WILD;
    Expr* p = rk_unwrap(r->pattern);
    bool verbatim = false;
    if (p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL &&
        p->data.function.head->data.symbol.name == SYM_Verbatim && p->data.function.arg_count == 1) {
        p = p->data.function.args[0];
        verbatim = true;
    }
    if (p->type != EXPR_FUNCTION) { *key = p; return RK_ATOM; }
    Expr* h = p->data.function.head;
    if (h->type != EXPR_SYMBOL) return RK_WILD;
    if (!verbatim) {
        if (rk_pattern_head(h->data.symbol.name)) return RK_WILD;
        for (size_t i = 0; i < p->data.function.arg_count; i++) {
            Expr* a = rk_unwrap(p->data.function.args[i]);
            if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL &&
                a->data.function.head->data.symbol.name == SYM_Optional)
                return RK_WILD;
        }
    }
    *key = h;
    return RK_HEAD;
}

static size_t rk_slot(const RuleIndex* ix, Expr* key, int kind) {
    size_t mask = ix->cap - 1;
    size_t h = (size_t)(expr_hash(key) * 2u + (uint64_t)kind) & mask;
    while (ix->tab[h].key &&
           !(ix->tab[h].kind == kind && expr_eq(ix->tab[h].key, key)))
        h = (h + 1) & mask;
    return h;
}

static void rule_index_build(RuleIndex* ix, const ReplaceRule* rules, size_t n) {
    ix->cap = 16;
    while (ix->cap < 2 * n) ix->cap *= 2;
    ix->tab = calloc(ix->cap, sizeof(RuleBucket));
    ix->wild = malloc((n ? n : 1) * sizeof(size_t));
    ix->nwild = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* key = NULL;
        int kind = rk_rule_key(&rules[i], &key);
        if (kind == RK_WILD) { ix->wild[ix->nwild++] = i; continue; }
        RuleBucket* b = &ix->tab[rk_slot(ix, key, kind)];
        if (!b->key) { b->key = key; b->kind = kind; }
        if (b->n == b->cap) {
            b->cap = b->cap ? 2 * b->cap : 2;
            b->idx = realloc(b->idx, b->cap * sizeof(size_t));
        }
        b->idx[b->n++] = i;
    }
}

static void rule_index_free(RuleIndex* ix) {
    for (size_t i = 0; i < ix->cap; i++) free(ix->tab[i].idx);
    free(ix->tab);
    free(ix->wild);
}

/* do_replace_all over an indexed rule set: only the node's bucket and the wild
 * rules are tried, in original order. */
static Expr* do_replace_all_idx(Expr* e, ReplaceRule* rules, const RuleIndex* ix) {
    const RuleBucket* b = NULL;
    if (e->type != EXPR_FUNCTION) {
        const RuleBucket* s = &ix->tab[rk_slot(ix, e, RK_ATOM)];
        if (s->key) b = s;
    } else if (e->data.function.head->type == EXPR_SYMBOL) {
        const RuleBucket* s = &ix->tab[rk_slot(ix, e->data.function.head, RK_HEAD)];
        if (s->key) b = s;
    }
    size_t i = 0, j = 0, nb = b ? b->n : 0;
    while (i < nb || j < ix->nwild) {
        size_t r = (j >= ix->nwild || (i < nb && b->idx[i] < ix->wild[j])) ? b->idx[i++]
                                                                           : ix->wild[j++];
        Expr* repl = rule_try(&rules[r], e);
        if (repl) return repl;
    }

    if (e->type == EXPR_FUNCTION && !rp_number_atom(e)) {
        Expr* new_head = do_replace_all_idx(e->data.function.head, rules, ix);
        bool assoc = assoc_is_wellformed(e);
        size_t count = e->data.function.arg_count;
        Expr** new_args = count ? malloc(sizeof(Expr*) * count) : NULL;
        for (size_t k = 0; k < count; k++)
            new_args[k] = do_replace_all_idx(struct_part(e, k, assoc), rules, ix);
        Expr* new_e = struct_rebuild(e, assoc, new_head, new_args, count, true);
        free(new_args);
        return new_e;
    }
    return expr_copy(e);
}

/* Dispatch[rules] -> rules (borrowed), noting that an index was asked for. */
static Expr* dispatch_unwrap(Expr* rules_expr, bool* indexed) {
    if (rules_expr->type == EXPR_FUNCTION && rules_expr->data.function.head->type == EXPR_SYMBOL &&
        rules_expr->data.function.head->data.symbol.name == SYM_Dispatch &&
        rules_expr->data.function.arg_count == 1) {
        if (indexed) *indexed = true;
        return rules_expr->data.function.args[0];
    }
    return rules_expr;
}

static Expr* apply_replace_all_nested(Expr* expr, Expr* rules_expr, bool indexed) {
    if (rules_expr->type == EXPR_FUNCTION && rules_expr->data.function.head->data.symbol.name == SYM_List) {
        if (rules_expr->data.function.arg_count > 0 && !is_rule(rules_expr->data.function.args[0])) {
            size_t count = rules_expr->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = apply_replace_all_nested(expr, rules_expr->data.function.args[i], indexed);
            }
            Expr* res = expr_new_function(expr_new_symbol(SYM_List), args, count);
            free(args);
            return res;
        }
    }

    size_t num_rules = 0;
    ReplaceRule* rules = rules_collect(rules_expr, &num_rules);

    Expr* result;
    if (num_rules == 0) {
        result = expr_copy(expr);
    } else if (indexed) {
        RuleIndex ix;
        rule_index_build(&ix, rules, num_rules);
        result = do_replace_all_idx(expr, rules, &ix);
        rule_index_free(&ix);
    } else {
        result = do_replace_all(expr, rules, num_rules);
    }
    rules_free(rules, num_rules);
    return result;
}

/* Dispatch[rules] -- a dispatch table over a rule set: ReplaceAll and
 * ReplaceRepeated index it by left-hand side (see above), and Replace,
 * ReplaceList and ReplaceAt accept it as the rules.  Normalisation only:
 * Dispatch[rule] is Dispatch[{rule}], Dispatch[Dispatch[r]] is Dispatch[r];
 * a list of rules (or of rule lists) stays as it is; anything else stays
 * unevaluated.  Normal[Dispatch[r]] gives r. */
Expr* builtin_dispatch(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* a = res->data.function.args[0];
    if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL &&
        a->data.function.head->data.symbol.name == SYM_Dispatch && a->data.function.arg_count == 1)
        return expr_copy(a);
    if (is_rule(a) && a->data.function.arg_count == 2) {
        Expr* l = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ expr_copy(a) }, 1);
        return expr_new_function(expr_new_symbol(SYM_Dispatch), &l, 1);
    }
    return NULL;
}

Expr* builtin_replace_all(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* expr = res->data.function.args[0];
    bool indexed = false;
    Expr* rules_expr = dispatch_unwrap(res->data.function.args[1], &indexed);

    return apply_replace_all_nested(expr, rules_expr, indexed);
}

/* ReplaceRepeated to a fixed point, at most `max_scans` scans (-1: no limit).
 * Stopping at the limit while the expression is still changing emits
 * ReplaceRepeated::rrlim naming the ORIGINAL expression, and returns the
 * expression reached. */
static Expr* apply_replace_repeated_nested(Expr* expr, Expr* rules_expr, int64_t max_scans,
                                           bool indexed) {
    if (rules_expr->type == EXPR_FUNCTION && rules_expr->data.function.head->data.symbol.name == SYM_List) {
        if (rules_expr->data.function.arg_count > 0 && !is_rule(rules_expr->data.function.args[0])) {
            size_t count = rules_expr->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = apply_replace_repeated_nested(expr, rules_expr->data.function.args[i], max_scans,
                                                        indexed);
            }
            Expr* res = expr_new_function(expr_new_symbol(SYM_List), args, count);
            free(args);
            return res;
        }
    }

    Expr* current = expr_copy(expr);
    for (int64_t scans = 0; max_scans < 0 || scans < max_scans; scans++) {
        Expr* next_raw = apply_replace_all_nested(current, rules_expr, indexed);
        Expr* next = eval_and_free(next_raw);
        if (expr_eq(current, next)) {
            expr_free(next);
            return current;
        }
        expr_free(current);
        current = next;
    }
    char* shown = expr_to_string_outputform(expr);
    mth_message("ReplaceRepeated", "rrlim", "Exiting after %s scanned %lld times.",
                shown ? shown : "", (long long)max_scans);
    free(shown);
    return current;
}

/* Decode a MaxIterations value: a non-negative integer, or Infinity (-1).
 * Returns false (with ReplaceRepeated::iopnf) for anything else. */
static bool rr_max_iterations(const Expr* v, int64_t* out) {
    if (v->type == EXPR_INTEGER && v->data.integer >= 0) { *out = v->data.integer; return true; }
    if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Infinity) { *out = -1; return true; }
    char* shown = expr_to_string_outputform((Expr*)v);
    mth_message("ReplaceRepeated", "iopnf",
                "Value of option MaxIterations -> %s should be a non-negative integer or Infinity.",
                shown ? shown : "");
    free(shown);
    return false;
}

/* ReplaceRepeated[expr, rules] / ReplaceRepeated[expr, rules, MaxIterations -> n]:
 * apply rules until the expression stops changing, scanning at most n times
 * (default Options[ReplaceRepeated], MaxIterations -> 65536). */
Expr* builtin_replace_repeated(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2 && argc != 3) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* rules_expr = res->data.function.args[1];

    int64_t max_scans = 65536;
    Expr* defs = symtab_get_options("ReplaceRepeated");        /* borrowed */
    for (size_t i = 0; defs && defs->type == EXPR_FUNCTION && i < defs->data.function.arg_count; i++) {
        Expr* o = defs->data.function.args[i];
        if (is_rule(o) && o->data.function.arg_count == 2 && o->data.function.args[0]->type == EXPR_SYMBOL &&
            o->data.function.args[0]->data.symbol.name == SYM_MaxIterations &&
            !rr_max_iterations(o->data.function.args[1], &max_scans))
            return NULL;
    }
    if (argc == 3) {
        Expr* o = res->data.function.args[2];
        if (!is_rule(o) || o->data.function.arg_count != 2 ||
            o->data.function.args[0]->type != EXPR_SYMBOL) return NULL;
        if (o->data.function.args[0]->data.symbol.name != SYM_MaxIterations) {
            char* shown = expr_to_string_outputform(res);
            mth_message("ReplaceRepeated", "optx", "Unknown option %s in %s.",
                        o->data.function.args[0]->data.symbol.name, shown ? shown : "");
            free(shown);
            return NULL;
        }
        if (!rr_max_iterations(o->data.function.args[1], &max_scans)) return NULL;
    }
    bool indexed = false;
    rules_expr = dispatch_unwrap(rules_expr, &indexed);
    return apply_replace_repeated_nested(expr, rules_expr, max_scans, indexed);
}


typedef struct {
    const ReplaceRule* rule;
    Expr* replacement;
    bool delayed;
    Expr** results;
    size_t count;
    size_t cap;
    int64_t limit;
} ReplaceListState;

static bool replacelist_callback(struct MatchEnv* env, void* user_data) {
    ReplaceListState* state = (ReplaceListState*)user_data;
    
    if (state->limit >= 0 && (int64_t)state->count >= state->limit) return true;
    
    Expr* repl = replace_bindings(state->replacement, env);
    if (state->delayed && symtab_rhs_has_scoped_condition(state->replacement)) {
        repl = rule_finish(state->rule, repl);   /* already evaluated */
        if (!repl) return false;                  /* guard failed: next match */
    } else if (state->delayed) {
        repl = eval_and_free(repl);
    }
    
    if (state->count == state->cap) {
        state->cap *= 2;
        state->results = realloc(state->results, sizeof(Expr*) * state->cap);
    }
    state->results[state->count++] = repl;
    
    return false;
}

Expr* builtin_replacelist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;

    Expr* expr = res->data.function.args[0];
    Expr* rules_expr = dispatch_unwrap(res->data.function.args[1], NULL);

    int64_t limit = -1;
    if (argc == 3) {
        Expr* l_expr = res->data.function.args[2];
        if (l_expr->type == EXPR_INTEGER) {
            limit = l_expr->data.integer;
        }
    }
    
    size_t num_rules = 0;
    ReplaceRule* rules = rules_collect(rules_expr, &num_rules);

    ReplaceListState state;
    state.results = malloc(sizeof(Expr*) * 16);
    state.count = 0;
    state.cap = 16;
    state.limit = limit;

    for (size_t i = 0; i < num_rules; i++) {
        if (state.limit >= 0 && (int64_t)state.count >= state.limit) break;
        if (!rules[i].pattern) {
            /* An association rule set: its literal keys match at most once. */
            Expr* repl = rule_try(&rules[i], expr);
            if (repl) {
                if (state.count == state.cap) {
                    state.cap *= 2;
                    state.results = realloc(state.results, sizeof(Expr*) * state.cap);
                }
                state.results[state.count++] = repl;
            }
            continue;
        }
        state.rule = &rules[i];
        state.replacement = rules[i].replacement;
        state.delayed = rules[i].delayed;
        
        MatchEnv* env = env_new();
        env->callback = replacelist_callback;
        env->callback_data = &state;
        
        match(expr, rules[i].pattern, env);
        env_free(env);
        
        if (state.limit >= 0 && (int64_t)state.count >= state.limit) {
            break;
        }
    }
    
    rules_free(rules, num_rules);
    
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), state.results, state.count);
    free(state.results);
    return list;
}

/* ------------------- ReplaceAt ------------------- */

/*
 * Parse rules_expr (a Rule, RuleDelayed, or List of such) into a freshly
 * allocated array of ReplaceRule. The pattern/replacement pointers borrow
 * from rules_expr; the array itself must be freed by the caller. Sets
 * *out_count to the number of valid rules collected.
 */
static ReplaceRule* parse_replace_rules(Expr* rules_expr, size_t* out_count) {
    return rules_collect(rules_expr, out_count);
}

/*
 * Try each rule in order against target. The first rule that matches wins.
 * Returns a freshly-owned expression: either the replacement (with bindings
 * substituted, evaluated for delayed rules) or expr_copy(target) when no
 * rule matches.
 */
static Expr* replaceat_apply_rules(ReplaceRule* rules, size_t num_rules, Expr* target) {
    for (size_t i = 0; i < num_rules; i++) {
        Expr* repl = rule_try(&rules[i], target);
        if (repl) return rules[i].delayed ? eval_and_free(repl) : repl;
    }
    return expr_copy(target);
}

/* Leaf action for ReplaceAt: try the rules at the addressed part. `ctx` is the
 * ReplaceAtCtx holding the parsed rules. */
typedef struct { ReplaceRule* rules; size_t count; } ReplaceAtCtx;

static Expr* replaceat_leaf(void* ctx, Expr* leaf) {
    ReplaceAtCtx* c = (ReplaceAtCtx*)ctx;
    return replaceat_apply_rules(c->rules, c->count, leaf);
}

/*
 * ReplaceAt[expr, rules, pos]   apply rules at one position, or at several.
 *
 * Position resolution -- integers (negatives from the end, 0 for the head),
 * All, Span, and Key[k]/literal keys on associations -- lives in the shared
 * walker expr_apply_at_path (src/part.c), which MapAt uses too. A position
 * that does not exist makes the walker return NULL and ReplaceAt stay
 * unevaluated, matching Part[{a, b, c}, 5].
 */
Expr* builtin_replace_at(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;

    Expr* expr = res->data.function.args[0];
    Expr* rules_expr = dispatch_unwrap(res->data.function.args[1], NULL);
    Expr* pos = res->data.function.args[2];

    size_t num_rules = 0;
    ReplaceRule* rules = parse_replace_rules(rules_expr, &num_rules);

    if (num_rules == 0) {
        rules_free(rules, num_rules);
        return expr_copy(expr);
    }

    ReplaceAtCtx ctx = { rules, num_rules };
    Expr* result = expr_apply_at_positions(expr, pos, replaceat_leaf, &ctx);

    rules_free(rules, num_rules);
    return result;
}

void replace_init(void) {
    symtab_add_builtin("ReplacePart", builtin_replace_part);
    symtab_add_builtin("Replace", builtin_replace);
    symtab_get_def("Replace")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ReplaceAll", builtin_replace_all);
    symtab_get_def("ReplaceAll")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Dispatch", builtin_dispatch);
    symtab_get_def("Dispatch")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Dispatch",
        "Dispatch[{lhs1 -> rhs1, lhs2 -> rhs2, ...}]\n"
        "\tgenerates an optimized dispatch table of rules: ReplaceAll (/.) and\n"
        "\tReplaceRepeated (//.) index the rules by left-hand side, so each subexpression\n"
        "\ttries only the rules that could match it.  Replace, ReplaceList and ReplaceAt\n"
        "\taccept it too, and Normal[Dispatch[rules]] gives the rules.");
    symtab_add_builtin("ReplaceRepeated", builtin_replace_repeated);
    symtab_get_def("ReplaceRepeated")->attributes |= ATTR_PROTECTED;
    {
        Expr* ra[2] = { expr_new_symbol("MaxIterations"), expr_new_integer(65536) };
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
        symtab_set_options("ReplaceRepeated", expr_new_function(expr_new_symbol(SYM_List), &rule, 1));
    }
    symtab_add_builtin("ReplaceList", builtin_replacelist);
    symtab_get_def("ReplaceList")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ReplaceAt", builtin_replace_at);
    symtab_get_def("ReplaceAt")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ReplaceAt",
        "ReplaceAt[expr, rules, n]\n"
        "\ttransforms expr by replacing the n-th element using rules.\n"
        "ReplaceAt[expr, rules, {i, j, ...}]\n"
        "\treplaces the part of expr at position {i, j, ...}.\n"
        "ReplaceAt[expr, rules, {{i1, j1, ...}, {i2, j2, ...}, ...}]\n"
        "\treplaces parts at several positions.\n"
        "\n"
        "Rules may be a single Rule/RuleDelayed or a list of them; rules are tried\n"
        "in order and the first match wins. Negative indices count from the end;\n"
        "0 targets the head. All and Span specifications are supported. On an\n"
        "association a position is a key, Key[k], or a positional index over the\n"
        "entries, and the rules are tried against the value. Repeated positions\n"
        "cause rules to be applied repeatedly to that part. ReplaceAt[expr, rules,\n"
        "{}] is an empty list of positions and replaces nothing, while {{}} is the\n"
        "position of expr itself. A position that does not exist leaves ReplaceAt\n"
        "unevaluated.");
}
