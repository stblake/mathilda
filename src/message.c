/*
 * message.c -- global message-suppression depth and the Quiet / Check /
 * Message builtins.  See message.h for rationale.
 */
#include "message.h"
#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "attr.h"
#include "sym_names.h"
#include <stdio.h>

/* Suppression nesting depth.  Zero => messages print; > 0 => suppressed. */
static int g_msg_suppress_depth = 0;

void mth_msg_suppress_push(void) { g_msg_suppress_depth++; }

void mth_msg_suppress_pop(void) {
    if (g_msg_suppress_depth > 0) g_msg_suppress_depth--;
}

int mth_msg_suppressed(void) { return g_msg_suppress_depth > 0; }

/* Save / restore the suppression depth as one value.  TimeConstrained's timeout
 * siglongjmp can unwind out of a Quiet[] region (its mth_msg_suppress_pop then
 * never runs, leaving messages silenced for the rest of the session); the guard
 * captures the depth before the body and restores it after the jump. */
int  mth_msg_suppress_depth_save(void)     { return g_msg_suppress_depth; }
void mth_msg_suppress_depth_load(int d)    { g_msg_suppress_depth = d; }

/* Depth for the `Solve::ifun` advisory scope (see message.h).  Zero => the
 * message may print; > 0 => a complete-solution caller (Reduce) is on the
 * stack, so the "use Reduce" advice is moot and suppressed. */
static int g_ifun_suppress_depth = 0;

void mth_msg_ifun_suppress_push(void) { g_ifun_suppress_depth++; }

void mth_msg_ifun_suppress_pop(void) {
    if (g_ifun_suppress_depth > 0) g_ifun_suppress_depth--;
}

int mth_msg_ifun_suppressed(void) { return g_ifun_suppress_depth > 0; }

/* Message-fired counter (see message.h). */
static unsigned long g_msg_fired = 0;
void          mth_msg_note_fired(void)  { g_msg_fired++; }
unsigned long mth_msg_fired_count(void) { return g_msg_fired; }

/* ------------------------------------------------- Quiet / Check / Message */

/* Quiet[expr] / Quiet[expr, spec]: evaluate expr with diagnostics suppressed
 * and return its value.  HoldAll -- the argument arrives unevaluated.  The
 * `spec` of the two-argument form (a message name or list) is IGNORED: all
 * messages are suppressed, a harmless superset of the requested set. */
Expr* builtin_quiet(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    mth_msg_suppress_push();
    Expr* val = evaluate(res->data.function.args[0]);
    mth_msg_suppress_pop();
    return val;   /* an in-flight Throw sentinel propagates unchanged */
}

/* Check[expr, failexpr] / Check[expr, failexpr, spec]: evaluate expr; if a
 * message fired during that evaluation, return the evaluated failexpr,
 * otherwise expr's value.  HoldAll.  A Throw propagates (it is not a
 * message).  The `spec` argument is ignored: any message counts. */
Expr* builtin_check(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    unsigned long before = mth_msg_fired_count();
    Expr* val = evaluate(res->data.function.args[0]);
    if (val && eval_is_inflight_throw(val)) return val;
    if (mth_msg_fired_count() != before) {
        if (val) expr_free(val);
        return evaluate(res->data.function.args[1]);
    }
    return val;
}

/* Message[sym::tag, e1, e2, ...]: note that a diagnostic fired (so an enclosing
 * Check sees it) and, unless suppressed, print its template.  HoldFirst -- the
 * message name is held; when a template string is defined for it, evaluating it
 * yields that string.  Argument substitution into the template is not
 * performed: the port uses Message only in error branches immediately followed
 * by Throw, so the printed text is diagnostic, not load-bearing.  Returns Null. */
Expr* builtin_message(Expr* res) {
    mth_msg_note_fired();
    if (!mth_msg_suppressed() && res->type == EXPR_FUNCTION &&
        res->data.function.arg_count >= 1) {
        Expr* mn = evaluate(res->data.function.args[0]);
        if (mn) {
            if (mn->type == EXPR_STRING)
                fprintf(stderr, "%s\n", mn->data.string);
            expr_free(mn);
        }
    }
    return expr_new_symbol(SYM_Null);
}

void message_init(void) {
    symtab_add_builtin("Quiet", builtin_quiet);
    symtab_get_def("Quiet")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Quiet",
        "Quiet[expr] evaluates expr with messages suppressed and returns its value. "
        "Quiet[expr, spec] suppresses all messages regardless of spec.");

    symtab_add_builtin("Check", builtin_check);
    symtab_get_def("Check")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Check",
        "Check[expr, failexpr] returns failexpr if a message is generated while "
        "evaluating expr, otherwise the value of expr.");

    symtab_add_builtin("Message", builtin_message);
    symtab_get_def("Message")->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_set_docstring("Message",
        "Message[sym::tag, args] prints the named message unless messages are suppressed.");
}
