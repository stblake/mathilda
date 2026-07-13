#include "purefunc.h"
#include "symtab.h"
#include "eval.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "attr.h"
#include "sym_names.h"

static uint32_t parse_pure_attr(Expr* attr_expr) {
    if (attr_expr->type == EXPR_SYMBOL) {
        const char* name = attr_expr->data.symbol.name;
        if (name == SYM_Listable) return ATTR_LISTABLE;
        if (name == SYM_Flat) return ATTR_FLAT;
        if (name == SYM_Orderless) return ATTR_ORDERLESS;
        if (name == SYM_NumericFunction) return ATTR_NUMERICFUNCTION;
        if (name == SYM_OneIdentity) return ATTR_ONEIDENTITY;
        if (name == SYM_HoldFirst) return ATTR_HOLDFIRST;
        if (name == SYM_HoldRest) return ATTR_HOLDREST;
        if (name == SYM_HoldAll) return ATTR_HOLDALL;
        if (name == SYM_HoldAllComplete) return ATTR_HOLDALLCOMPLETE;
        if (name == SYM_SequenceHold) return ATTR_SEQUENCEHOLD;
        if (name == SYM_NHoldRest) return ATTR_NHOLDREST;
        if (name == SYM_Protected) return ATTR_PROTECTED;
    }
    return ATTR_NONE;
}

uint32_t pure_function_attributes(Expr* func) {
    if (func->type == EXPR_FUNCTION && func->data.function.head->type == EXPR_SYMBOL &&
        func->data.function.head->data.symbol.name == SYM_Function) {
        /* Mathematica default: Function[params, body] has no Hold attributes,
         * so its arguments are evaluated normally before being substituted
         * into body. The 3-arg form Function[params, body, attrs] can opt in
         * to HoldAll / HoldFirst / HoldRest / HoldAllComplete and a handful
         * of structural attributes. */
        uint32_t attrs = ATTR_NONE;
        if (func->data.function.arg_count >= 3) {
            Expr* attr_spec = func->data.function.args[2];
            if (attr_spec->type == EXPR_SYMBOL) {
                attrs |= parse_pure_attr(attr_spec);
            } else if (attr_spec->type == EXPR_FUNCTION &&
                       attr_spec->data.function.head->type == EXPR_SYMBOL &&
                       attr_spec->data.function.head->data.symbol.name == SYM_List) {
                for (size_t i = 0; i < attr_spec->data.function.arg_count; i++) {
                    attrs |= parse_pure_attr(attr_spec->data.function.args[i]);
                }
            }
        }
        return attrs;
    }
    return ATTR_NONE;
}

void purefunc_init(void) {
    symtab_add_builtin("Slot", builtin_slot);
    symtab_add_builtin("SlotSequence", builtin_slotsequence);
    symtab_get_def("Slot")->attributes |= ATTR_PROTECTED;
    symtab_get_def("SlotSequence")->attributes |= ATTR_PROTECTED;
}

Expr* builtin_slot(Expr* res) {
    (void)res;
    // Slot stays unevaluated. 
    // It's only replaced during apply_pure_function.
    return NULL;
}

Expr* builtin_slotsequence(Expr* res) {
    (void)res;
    // SlotSequence stays unevaluated.
    return NULL;
}

static Expr* substitute_slots(Expr* e, Expr** args, size_t arg_count) {
    if (!e) return NULL;
    
    if (e->type == EXPR_FUNCTION) {
        // If it's another Function, don't recurse into it. 
        // We only substitute slots for the current level.
        if (e->data.function.head->type == EXPR_SYMBOL && 
            e->data.function.head->data.symbol.name == SYM_Function) {
            return expr_copy(e);
        }
        
        // Check if this function is Slot[n] or SlotSequence[n]
        if (e->data.function.head->type == EXPR_SYMBOL) {
            const char* head = e->data.function.head->data.symbol.name;
            
            if (head == SYM_Slot && e->data.function.arg_count == 1) {
                Expr* arg = e->data.function.args[0];
                if (arg->type == EXPR_INTEGER) {
                    int64_t idx = arg->data.integer;
                    if (idx >= 1 && (size_t)idx <= arg_count) {
                        return expr_copy(args[idx - 1]);
                    }
                }
            }
            
            if (head == SYM_SlotSequence && e->data.function.arg_count == 1) {
                Expr* arg = e->data.function.args[0];
                if (arg->type == EXPR_INTEGER) {
                    int64_t idx = arg->data.integer;
                    if (idx >= 1 && (size_t)idx <= arg_count) {
                        size_t seq_count = arg_count - (size_t)idx + 1;
                        Expr** seq_args = malloc(sizeof(Expr*) * seq_count);
                        for (size_t i = 0; i < seq_count; i++) {
                            seq_args[i] = expr_copy(args[(size_t)idx - 1 + i]);
                        }
                        Expr* seq = expr_new_function(expr_new_symbol(SYM_Sequence), seq_args, seq_count);
                        free(seq_args);
                        return seq;
                    }
                }
            }
        }
        
        // Standard recursion for other functions
        Expr** new_args = malloc(sizeof(Expr*) * e->data.function.arg_count);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            new_args[i] = substitute_slots(e->data.function.args[i], args, arg_count);
        }
        Expr* new_head = substitute_slots(e->data.function.head, args, arg_count);
        Expr* res = expr_new_function(new_head, new_args, e->data.function.arg_count);
        free(new_args);
        return res;
    }
    
    return expr_copy(e);
}

/*
 * substitute_names:
 * Walks `e` and replaces free occurrences of any parameter symbol in
 * `names` with the corresponding expression in `vals`. The substitution is
 * lexical (no evaluation) and does not recurse into nested Function
 * expressions -- those shadow the names in their own right, so the outer
 * Function's parameters must not leak into them.
 * The returned expression is freshly allocated; the caller owns it and
 * `vals` is not consumed (copies are made as needed).
 */
static Expr* substitute_names(Expr* e, char** names, Expr** vals, size_t count) {
    if (!e) return NULL;

    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < count; i++) {
            if (names[i] && strcmp(e->data.symbol.name, names[i]) == 0) {
                return expr_copy(vals[i]);
            }
        }
        return expr_copy(e);
    }

    if (e->type == EXPR_FUNCTION) {
        /* Do not descend into a nested Function -- its body scopes new
         * parameters and the outer parameters may be shadowed. We copy
         * the whole nested Function as-is. */
        if (e->data.function.head->type == EXPR_SYMBOL &&
            e->data.function.head->data.symbol.name == SYM_Function) {
            return expr_copy(e);
        }

        Expr** new_args = malloc(sizeof(Expr*) * e->data.function.arg_count);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            new_args[i] = substitute_names(e->data.function.args[i], names, vals, count);
        }
        Expr* new_head = substitute_names(e->data.function.head, names, vals, count);
        Expr* res = expr_new_function(new_head, new_args, e->data.function.arg_count);
        free(new_args);
        return res;
    }

    return expr_copy(e);
}

/*
 * apply_pure_function:
 * Applies a Function[...] expression to the given arguments.
 *
 * Four input shapes are supported:
 *   Function[body]                    -- body uses # / ## slots
 *   Function[Null, body, attrs]       -- slot form with attributes
 *   Function[x, body]                 -- single named parameter
 *   Function[{x1, ...}, body]         -- list of named parameters
 *   Function[{x1, ...}, body, attrs]  -- named parameters with attributes
 *
 * Parameter binding is performed by lexical substitution (not by mutating
 * the global symbol table) so that references wrapped in Unevaluated still
 * see the substituted expression. Nested Function expressions are treated
 * as opaque to avoid capture.
 */
/*
 * trap_return:
 * Strip a Return[v] / Return[] marker that targets *this* Function
 * boundary. Return[v, h] is consumed only when h == Function; otherwise
 * the marker is handed upward unchanged so an enclosing Module / Block
 * / loop with the matching head can claim it.
 *
 * Takes ownership of `result` and returns either the original pointer
 * (when the marker is absent or propagating) or a freshly-owned
 * unwrapped value (when the marker is consumed). Frees `result` on the
 * consume path.
 */
static Expr* trap_return(Expr* result) {
    Expr* rv = NULL;
    EvalReturnAction ra = eval_classify_return(result, SYM_Function, &rv);
    if (ra == EVAL_RETURN_CONSUME) {
        expr_free(result);
        return rv;
    }
    return result;
}

Expr* apply_pure_function(Expr* head, Expr** args, size_t arg_count) {
    if (head->type != EXPR_FUNCTION) return NULL;

    size_t head_argc = head->data.function.arg_count;

    /* Function[body] -- slot form with 1 argument */
    if (head_argc == 1) {
        Expr* body = head->data.function.args[0];
        Expr* substituted = substitute_slots(body, args, arg_count);
        Expr* result = evaluate(substituted);
        expr_free(substituted);
        return trap_return(result);
    }

    if (head_argc < 2) return NULL;

    Expr* params = head->data.function.args[0];
    Expr* body = head->data.function.args[1];

    /* Function[Null, body, ...] -- slot form with attributes.
     * (Null indicates "use slots" rather than named parameters.) */
    if (params->type == EXPR_SYMBOL && params->data.symbol.name == SYM_Null) {
        Expr* substituted = substitute_slots(body, args, arg_count);
        Expr* result = evaluate(substituted);
        expr_free(substituted);
        return trap_return(result);
    }

    /* Function[x, body] -- single named parameter */
    if (params->type == EXPR_SYMBOL) {
        char* name = params->data.symbol.name;
        Expr* val = (arg_count >= 1) ? args[0] : params;
        Expr* substituted = substitute_names(body, &name, &val, 1);
        Expr* result = evaluate(substituted);
        expr_free(substituted);
        return trap_return(result);
    }

    /* Function[{x1, x2, ...}, body] -- list of named parameters */
    if (params->type == EXPR_FUNCTION &&
        params->data.function.head->type == EXPR_SYMBOL &&
        params->data.function.head->data.symbol.name == SYM_List) {
        size_t var_count = params->data.function.arg_count;
        if (var_count == 0) {
            Expr* result = evaluate(body);
            return trap_return(result);
        }

        char** names = malloc(sizeof(char*) * var_count);
        Expr** vals = malloc(sizeof(Expr*) * var_count);
        for (size_t i = 0; i < var_count; i++) {
            Expr* var = params->data.function.args[i];
            if (var->type == EXPR_SYMBOL) {
                names[i] = var->data.symbol.name;
                /* Missing arguments leave the parameter symbolic, which is
                 * the same behaviour as unbound slots. */
                vals[i] = (i < arg_count) ? args[i] : var;
            } else {
                names[i] = NULL;
                vals[i] = var;
            }
        }

        Expr* substituted = substitute_names(body, names, vals, var_count);
        Expr* result = evaluate(substituted);
        expr_free(substituted);
        free(names);
        free(vals);
        return trap_return(result);
    }

    return NULL;
}
