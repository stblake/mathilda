/*
 * loadmodule.c -- runtime loading of Mathilda (.m) source modules.
 *
 * See loadmodule.h for the contract. The file-reading loop here is the same
 * one Get[] uses (builtin_get in readwrite.c is now a thin wrapper over
 * mathilda_run_file); mathilda_load_module adds working-directory-independent
 * path resolution and load-once bookkeeping on top.
 */

#include "loadmodule.h"
#include "sym_names.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "parse.h"
#include "eval.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Expr* mathilda_run_file(const char* path, int* opened) {
    FILE* fp = fopen(path, "r");
    if (!fp) {
        if (opened) *opened = 0;
        return NULL;
    }
    if (opened) *opened = 1;

    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    if (fsize < 0) {
        fclose(fp);
        return expr_new_symbol(SYM_DollarFailed);
    }

    char* buffer = malloc((size_t)fsize + 1);
    if (!buffer) {
        fclose(fp);
        return expr_new_symbol(SYM_DollarFailed);
    }

    size_t read_len = fread(buffer, 1, (size_t)fsize, fp);
    buffer[read_len] = '\0';
    fclose(fp);

    Expr* last_eval = expr_new_symbol(SYM_Null);
    const char* ptr = buffer;
    while (*ptr != '\0') {
        Expr* parsed = parse_next_expression(&ptr);
        if (!parsed) break; /* EOF or unparseable: stop. */
        Expr* evaluated = evaluate(parsed);
        if (evaluated) {
            expr_free(last_eval);
            last_eval = evaluated;
        }
        expr_free(parsed);
    }

    free(buffer);
    return last_eval;
}

/* Modules loaded so far, keyed by the relpath requested. Bounded; a kernel
 * has a small fixed set of internal modules. */
#define LM_MAX_LOADED 256
#define LM_RELPATH_MAX 256
static char lm_loaded[LM_MAX_LOADED][LM_RELPATH_MAX];
static int  lm_loaded_count = 0;

static int lm_already_loaded(const char* relpath) {
    for (int i = 0; i < lm_loaded_count; i++) {
        if (strcmp(lm_loaded[i], relpath) == 0) return 1;
    }
    return 0;
}

static void lm_mark_loaded(const char* relpath) {
    if (lm_loaded_count >= LM_MAX_LOADED) return; /* silently cap */
    size_t n = strlen(relpath);
    if (n >= LM_RELPATH_MAX) return;
    memcpy(lm_loaded[lm_loaded_count], relpath, n + 1);
    lm_loaded_count++;
}

/* Try one candidate path; on success (file opened) load it and return 1. */
static int lm_try_path(const char* full) {
    int opened = 0;
    Expr* last = mathilda_run_file(full, &opened);
    if (last) expr_free(last);
    return opened;
}

int mathilda_load_module(const char* relpath) {
    if (!relpath || !*relpath) return 0;
    if (lm_already_loaded(relpath)) return 1;

    char full[1024];

    /* 1. $MATHILDA_HOME/<relpath> */
    const char* home = getenv("MATHILDA_HOME");
    if (home && *home) {
        snprintf(full, sizeof(full), "%s/%s", home, relpath);
        if (lm_try_path(full)) { lm_mark_loaded(relpath); return 1; }
    }

    /* 2..5. Relative ladder rooted at the source tree's internal directory. */
    static const char* const bases[] = {
        "src/internal/",
        "../src/internal/",
        "../../src/internal/",
        "../../../src/internal/",
        NULL
    };
    for (int i = 0; bases[i]; i++) {
        snprintf(full, sizeof(full), "%s%s", bases[i], relpath);
        if (lm_try_path(full)) { lm_mark_loaded(relpath); return 1; }
    }

    fprintf(stderr,
        "LoadModule::nofile: cannot locate src/internal/%s on disk.\n",
        relpath);
    return 0;
}

/* LoadModule["relpath"] -> True if the module was found and loaded (or had
 * already been loaded), False otherwise. */
Expr* builtin_loadmodule(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING) return NULL;
    int ok = mathilda_load_module(arg->data.string);
    return expr_new_symbol(ok ? "True" : "False");
}

void loadmodule_init(void) {
    symtab_add_builtin("LoadModule", builtin_loadmodule);
    symtab_get_def("LoadModule")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LoadModule",
        "LoadModule[\"relpath\"]\n"
        "\tloads the internal Mathilda source module at relpath (relative to\n"
        "\tsrc/internal), resolving the location independently of the current\n"
        "\tworking directory. Each module is loaded at most once. Returns True\n"
        "\tif the module was located and loaded, False otherwise.");
}
