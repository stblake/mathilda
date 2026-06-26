#include "expr.h"
#include "parse.h"
#include "print.h"
#include "part.h"
#include "eval.h"
#include "symtab.h"
#include "repl_hooks.h"
#include "sym_names.h"
#include "show.h"
#include "render3d.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <readline/readline.h>
#include <readline/history.h>

#define MAX_INPUT_LEN 10240

/* True if `s` consists only of whitespace and (* ... *) comments
 * (nested comments allowed). Used to distinguish a no-op line from a
 * genuine parse failure, so the REPL doesn't shout "Parse error" at a
 * stray comment. */
static int is_blank_or_comment_only(const char* s) {
    while (*s) {
        if (isspace((unsigned char)*s)) {
            s++;
        } else if (s[0] == '(' && s[1] == '*') {
            int depth = 1;
            s += 2;
            while (*s && depth > 0) {
                if (s[0] == '(' && s[1] == '*') { depth++; s += 2; }
                else if (s[0] == '*' && s[1] == ')') { depth--; s += 2; }
                else { s++; }
            }
            if (depth > 0) return 0;  /* unterminated comment is a real error */
        } else {
            return 0;
        }
    }
    return 1;
}

void process_input(const char* input, int line_number) {
    if (strlen(input) == 0) return;
    if (is_blank_or_comment_only(input)) return;

    // Update $Line
    Expr* line_sym = expr_new_symbol(SYM_DollarLine);
    Expr* line_val = expr_new_integer(line_number);
    symtab_add_own_value("$Line", line_sym, line_val);
    expr_free(line_sym);
    expr_free(line_val);

    /* $PreRead: text-level hook applied to the raw input string
     * before parsing. Pass-through (a strdup of `input`) when unset. */
    char* cooked_input = repl_apply_pre_read(input);
    if (!cooked_input) {
        printf("Parse error\n\n");
        return;
    }

    // Parse the (possibly hook-transformed) input
    Expr* parsed = parse_expression(cooked_input);
    free(cooked_input);
    if (!parsed) {
        printf("Parse error\n\n");
        return;
    }

    /* Store In[line_number] = parsed BEFORE running $Pre/$Post so that
     * In[n] reflects what the user typed, not whatever a hook did. */
    Expr* in_sym = expr_new_symbol(SYM_In);
    Expr* in_arg = expr_new_integer(line_number);
    Expr* in_args[] = {in_arg};
    Expr* in_pattern = expr_new_function(in_sym, in_args, 1);
    symtab_add_down_value("In", in_pattern, parsed);
    expr_free(in_pattern);

    /* $Pre: applied to the parsed expression. Consumes our reference
     * to `parsed`; we treat the result as the new input to evaluate. */
    Expr* pre_input = repl_apply_pre(expr_copy(parsed));
    Expr* evaluated = evaluate(pre_input);
    expr_free(pre_input);

    /* $Post: applied to the evaluator's result. Consumes our reference
     * to `evaluated`. */
    evaluated = repl_apply_post(evaluated);

    /* A NULL result means the evaluation produced nothing displayable
     * (e.g. a hook absorbed the value). Skip Out[n] storage and the
     * "Out[n]= " banner rather than crashing in expr_copy. */
    if (!evaluated) {
        expr_free(parsed);
        return;
    }

    /* Store Out[line_number] = evaluated (post-$Post, pre-$PrePrint:
     * Mathematica's documented ordering). */
    Expr* out_sym = expr_new_symbol(SYM_Out);
    Expr* out_arg = expr_new_integer(line_number);
    Expr* out_args[] = {out_arg};
    Expr* out_pattern = expr_new_function(out_sym, out_args, 1);
    symtab_add_down_value("Out", out_pattern, evaluated);
    expr_free(out_pattern);

    /* $PrePrint: applied only for display. Out[n] keeps the
     * pre-$PrePrint value above; here we render a possibly modified
     * copy. */
    Expr* to_print = repl_apply_pre_print(expr_copy(evaluated));

    printf("Out[%d]= ", line_number);
    expr_print(to_print);
    printf("\n"); // extra blank line

    /* Mathematica's front end auto-displays a top-level Graphics[...] (or
     * Graphics3D[...], from Plot3D) result. This REPL is the sole "front
     * end", so it owns rendering: Show[]/Plot[]/Plot3D[] merely return such
     * an object and we render it here. Routing every display through one
     * path means `g // Graphics`, Show[...], Plot[...] and Plot3D[...] all
     * render identically, and a trailing `;` (which yields Null) correctly
     * suppresses the window. graphics_show/graphics3d_show borrow the expr
     * (no ownership transfer); on a non-graphics build their stubs print a
     * one-line "install raylib" hint instead. */
    if (to_print && to_print->type == EXPR_FUNCTION
        && to_print->data.function.head->type == EXPR_SYMBOL
        && to_print->data.function.arg_count >= 1) {
        if (to_print->data.function.head->data.symbol == SYM_Graphics) graphics_show(to_print);
        else if (to_print->data.function.head->data.symbol == SYM_Graphics3D) graphics3d_show(to_print);
    }

    expr_free(to_print);
    expr_free(parsed);
    expr_free(evaluated);
}

void repl_loop() {
    printf("\nMathilda - A tiny, LLM-generated, Mathematica-like computer algebra system.\n\n");
    printf("This program is free, open source software and comes with ABSOLUTELY NO WARRANTY.\n\n");
    printf("End a line with '\\' to enter a multiline expression. Press Return to evaluate.\n");
    printf("Exit by evaluating Quit[] or CONTROL-C.\n\n");
    
    int line_number = 1;
    char prompt[64];
    char full_input[MAX_INPUT_LEN] = {0};
    int in_multiline = 0;
    
    while (true) {
        int submit_now = 0;
        if (!in_multiline) {
            snprintf(prompt, sizeof(prompt), "In[%d]:= ", line_number);
        } else {
            prompt[0] = '\0';
        }
        
        char* line = readline(prompt);
        if (!line) {
            printf("\n");
            /* EOF: run $Epilog before tearing down. */
            repl_apply_epilog();
            break;
        }
        
        size_t line_len = strlen(line);
        
        // Remove trailing whitespace to properly check for backslash
        size_t check_len = line_len;
        while (check_len > 0 && (line[check_len - 1] == ' ' || line[check_len - 1] == '\t')) {
            check_len--;
        }
        
        int has_backslash = 0;
        if (check_len > 0 && line[check_len - 1] == '\\') {
            has_backslash = 1;
            // Remove the backslash from the line
            line[check_len - 1] = '\0';
            line_len = strlen(line);
        }
        
        // Check buffer limits
        if (strlen(full_input) + line_len + 2 >= MAX_INPUT_LEN) {
            printf("Input too long!\n");
            full_input[0] = '\0';
            in_multiline = 0;
            free(line);
            continue;
        }
        
        if (in_multiline) {
            strcat(full_input, "\n");
        }
        strcat(full_input, line);
        
        if (has_backslash) {
            in_multiline = 1;
        } else {
            submit_now = 1;
        }
        
        if (submit_now) {
            if (strlen(full_input) == 0) {
                free(line);
                continue;
            }
            
            add_history(full_input);
            
            if (strcmp(full_input, "Quit[]") == 0) {
                /* User-requested shutdown: run $Epilog first. */
                repl_apply_epilog();
                free(line);
                break;
            }
            
            process_input(full_input, line_number);
            
            full_input[0] = '\0';
            in_multiline = 0;
            line_number++;
        }
        
        free(line);
    }
    
    printf("\n");
}

#include "core.h"

int main() {
    symtab_init();
    core_init();
    
    // Load internal init.m silently if it exists
    FILE* check_fp = fopen("./src/internal/init.m", "r");
    if (check_fp) {
        fclose(check_fp);
        Expr* init_call = parse_expression("Get[\"./src/internal/init.m\"]");
        if (init_call) {
            Expr* res = evaluate(init_call);
            expr_free(init_call);
            if (res) expr_free(res);
        }
    }
    
    repl_loop();
    return 0;
}
