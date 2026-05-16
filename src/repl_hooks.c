/*
 * REPL session hooks: $PreRead, $Pre, $Post, $PrePrint, $Epilog.
 *
 * Each hook is a global symbol whose OwnValue, if assigned, is applied
 * at a specific phase of the REPL loop. The implementation is split
 * into small helpers so the read/eval/print phases of repl.c stay
 * legible and so each hook can be exercised directly from unit tests
 * without having to drive the readline loop.
 *
 * Detection rule: a hook is "set" iff symtab_get_own_values(sym)
 * returns a non-empty rule list. The init function below registers the
 * symbols (so docstrings work) but installs no default OwnValue, so
 * out-of-the-box every hook is a no-op.
 *
 * All evaluation is funnelled through the standard evaluator: the
 * helper builds, e.g., $Pre[expr] and calls evaluate(). This means
 * pattern matching, attribute handling and recursion guards behave
 * exactly as they would for any user-issued call.
 */

#include "repl_hooks.h"

#include "eval.h"
#include "expr.h"
#include "symtab.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local strdup so the hook module does not depend on the POSIX
 * implementation; CLAUDE.md forbids unguarded strdup elsewhere but
 * src/expr.c already relies on it via the platform's <string.h>. We
 * keep our own small copy here purely so this file's intent is
 * self-contained. */
static char* hooks_strdup(const char* s) {
    if (!s) return NULL;
    size_t n = strlen(s);
    char* out = malloc(n + 1);
    if (!out) return NULL;
    memcpy(out, s, n + 1);
    return out;
}

/* True iff the named hook symbol has at least one OwnValue assigned. */
static bool hook_is_set(const char* hook_name) {
    return symtab_get_own_values(hook_name) != NULL;
}

/*
 * Build hook_name[arg] and evaluate it.
 *
 * Takes ownership of `arg`. The arg is wrapped into a freshly built
 * function call whose head is a fresh symbol node naming the hook.
 * After evaluation the wrapper call (and transitively `arg`) is
 * released; the freshly evaluated result is returned to the caller,
 * who must expr_free() it.
 */
static Expr* hook_call_eval(const char* hook_name, Expr* arg) {
    Expr* head = expr_new_symbol(hook_name);
    Expr* args[1];
    args[0] = arg; /* ownership transfers to the call node */
    Expr* call = expr_new_function(head, args, 1);
    Expr* result = evaluate(call);
    expr_free(call);
    return result;
}

void repl_hooks_init(void) {
    /* Touch each symbol so it exists in the symbol table even when no
     * OwnValue is ever assigned; this lets ?$PreRead and friends
     * surface the docstring without first requiring an assignment. */
    symtab_get_def("$PreRead");
    symtab_get_def("$Pre");
    symtab_get_def("$Post");
    symtab_get_def("$PrePrint");
    symtab_get_def("$Epilog");

    symtab_set_docstring("$PreRead",
        "$PreRead\n"
        "\tis a global variable whose value, if set, is applied to the\n"
        "\ttext of every input expression before it is fed to Mathilda.\n"
        "\n"
        "The value of $PreRead must be a function that takes one string\n"
        "argument and returns a string. If unset, the input is parsed\n"
        "without modification. If the hook returns a non-string value\n"
        "the original input is used and a $PreRead::strret diagnostic\n"
        "is emitted.");

    symtab_set_docstring("$Pre",
        "$Pre\n"
        "\tis a global variable whose value, if set, is applied to every\n"
        "\tinput expression after parsing and before standard evaluation.\n"
        "\n"
        "Unless $Pre is assigned to a head with HoldAll, the wrapped\n"
        "expression is evaluated before $Pre sees it -- which makes the\n"
        "effect indistinguishable from $Post.");

    symtab_set_docstring("$Post",
        "$Post\n"
        "\tis a global variable whose value, if set, is applied to every\n"
        "\toutput expression after evaluation.");

    symtab_set_docstring("$PrePrint",
        "$PrePrint\n"
        "\tis a global variable whose value, if set, is applied to every\n"
        "\texpression just before it is printed. Out[n] is assigned the\n"
        "\tunmodified result, but the printed form reflects the value\n"
        "\treturned by $PrePrint.");

    symtab_set_docstring("$Epilog",
        "$Epilog\n"
        "\tis a symbol whose value, if any, is evaluated once when the\n"
        "\tMathilda session terminates (via Quit[] or EOF).");
}

char* repl_apply_pre_read(const char* input) {
    if (!input) return NULL;
    if (!hook_is_set("$PreRead")) {
        return hooks_strdup(input);
    }

    Expr* str_arg = expr_new_string(input);
    Expr* result = hook_call_eval("$PreRead", str_arg);

    char* out;
    if (result && result->type == EXPR_STRING) {
        out = hooks_strdup(result->data.string);
    } else {
        fprintf(stderr,
                "$PreRead::strret: $PreRead returned a non-string value; "
                "using original input.\n");
        out = hooks_strdup(input);
    }
    if (result) expr_free(result);
    return out;
}

Expr* repl_apply_pre(Expr* expr) {
    if (!expr) return NULL;
    if (!hook_is_set("$Pre")) return expr;
    return hook_call_eval("$Pre", expr);
}

Expr* repl_apply_post(Expr* expr) {
    if (!expr) return NULL;
    if (!hook_is_set("$Post")) return expr;
    return hook_call_eval("$Post", expr);
}

Expr* repl_apply_pre_print(Expr* expr) {
    if (!expr) return NULL;
    if (!hook_is_set("$PrePrint")) return expr;
    return hook_call_eval("$PrePrint", expr);
}

void repl_apply_epilog(void) {
    if (!hook_is_set("$Epilog")) return;
    Expr* sym = expr_new_symbol("$Epilog");
    Expr* result = evaluate(sym);
    expr_free(sym);
    if (result) expr_free(result);
}
