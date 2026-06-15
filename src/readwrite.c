/*
 * readwrite.c - File I/O builtins (Get, Put).
 *
 * Get reads Mathilda source from a file and evaluates each expression,
 * returning the last value (used by the REPL bootstrap to load
 * the internal .m initialization files).
 *
 * Put writes one or more expressions to a file in InputForm so the
 * output can be read back with Get.  The parser also recognises the
 * infix shorthand `expr >> "file"` and lowers it to `Put[expr, "file"]`.
 */

#include "readwrite.h"
#include "sym_names.h"
#include "loadmodule.h"
#include "expr.h"
#include "symtab.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Expr* builtin_get(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;

    Expr* file_arg = res->data.function.args[0];
    if (file_arg->type != EXPR_STRING) return NULL;

    const char* filename = file_arg->data.string;
    int opened = 0;
    Expr* last_eval = mathilda_run_file(filename, &opened);
    if (!opened) {
        printf("Get::noopen: Cannot open %s.\n", filename);
        return expr_new_symbol("$Failed");
    }
    /* mathilda_run_file returns Null for an empty file; never NULL when
     * opened succeeded. */
    return last_eval ? last_eval : expr_new_symbol(SYM_Null);
}

/* Shared implementation for Put / PutAppend.  `mode` selects fopen
 * mode: "w" truncates ("Put"), "a" appends ("PutAppend").  `who` names
 * the calling builtin for diagnostics. */
static Expr* put_common(Expr* res, const char* mode, const char* who) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return NULL;

    Expr* last = res->data.function.args[argc - 1];
    if (last->type != EXPR_STRING) return NULL;

    const char* filename = last->data.string;
    FILE* fp = fopen(filename, mode);
    if (!fp) {
        printf("%s::noopen: Cannot open %s.\n", who, filename);
        return expr_new_symbol("$Failed");
    }

    /* When only the filename is supplied (Put["file"]) the loop runs
     * zero times — leaving an empty file in "w" mode and leaving the
     * file untouched (besides creating it if absent) in "a" mode. */
    for (size_t i = 0; i + 1 < argc; i++) {
        char* s = expr_to_string(res->data.function.args[i]);
        if (s) {
            fputs(s, fp);
            free(s);
        }
        fputc('\n', fp);
    }

    fclose(fp);
    return expr_new_symbol(SYM_Null);
}

/* Put[expr_1, ..., expr_n, "file"] truncates "file" and writes each
 * expr_i as standard-form text followed by a newline.  The single-arg
 * form Put["file"] creates an empty file. */
Expr* builtin_put(Expr* res) {
    return put_common(res, "w", "Put");
}

/* PutAppend[expr_1, ..., expr_n, "file"] appends each expr_i to the
 * end of "file", creating the file if it does not exist.  Single-arg
 * form PutAppend["file"] creates an empty file if it doesn't exist
 * but does not truncate an existing file. */
Expr* builtin_putappend(Expr* res) {
    return put_common(res, "a", "PutAppend");
}

void readwrite_init(void) {
    /* Docstrings for Get / Put / PutAppend live in info_init() — keep
     * this module focused on builtin registration and attributes. */
    symtab_add_builtin("Get", builtin_get);
    symtab_get_def("Get")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("Put", builtin_put);
    symtab_get_def("Put")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("PutAppend", builtin_putappend);
    symtab_get_def("PutAppend")->attributes |= ATTR_PROTECTED;
}
