/*
 * parallel.c -- ParallelMap, ParallelTable, ParallelEvaluate, Kernels,
 * LaunchKernels, DistributeDefinitions: the parallel-computation interface of
 * a kernel with no subkernels.
 *
 * Mathilda runs a single kernel ($KernelCount is 0, $KernelID is 0).  The
 * faithful meaning of the parallel functions there is the serial one:
 *   ParallelMap[f, expr, ...]   is Map[f, expr, ...]
 *   ParallelTable[body, it...]  is Table[body, it...]
 *   ParallelEvaluate[expr]      evaluates expr on every subkernel -- of which
 *                               there are none -- so it gives {}
 *   Kernels[]                   is {}
 *   DistributeDefinitions[...]  distributes to no kernel and gives {}
 *   LaunchKernels[...]          launches nothing: a message and {}
 * Trailing options of the parallel heads (DistributedContexts, Method ->
 * "CoarsestGrained" | "ItemsPerEvaluation" -> ..., ...) only steer how work is
 * split between subkernels, so they are dropped.  Packages that guard their
 * parallel paths on $KernelCount / $KernelID (FeynCalc's $ParallelizeFeynCalc)
 * therefore run unchanged.
 */
#include "parallel.h"
#include "expr.h"
#include "symtab.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "attr.h"
#include "message.h"
#include <stdlib.h>

/* Is e a Rule / RuleDelayed (an option)? */
static int is_option(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rule ||
            e->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
           e->data.function.arg_count == 2;
}

/* head[args of res without its trailing options] -- copies, so the evaluator
 * is free to release res.  The first `min` arguments are positional whatever
 * they look like -- ParallelMap[f, a -> b] maps over the rule, it has no
 * options -- so only arguments after them can be options.  NULL when fewer
 * than `min` arguments are given. */
static Expr* rebuild_without_options(Expr* res, const char* head, size_t min) {
    size_t n = res->data.function.arg_count;
    while (n > min && is_option(res->data.function.args[n - 1])) n--;
    if (n < min) return NULL;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    if (!a) return NULL;
    for (size_t i = 0; i < n; i++) a[i] = expr_copy(res->data.function.args[i]);
    Expr* out = expr_new_function(expr_new_symbol(head), a, n);
    free(a);
    return out;
}

static Expr* empty_list(void) {
    return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
}

/* ParallelMap[f, expr] / ParallelMap[f, expr, levelspec] -> Map[...]. */
static Expr* builtin_parallelmap(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    return rebuild_without_options(res, "Map", 2);
}

/* ParallelTable[body, iterators...] -> Table[...] (HoldAll, as Table). */
static Expr* builtin_paralleltable(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    return rebuild_without_options(res, "Table", 1);
}

/* ParallelEvaluate[expr] / ParallelEvaluate[expr, kernels]: no subkernels. */
static Expr* builtin_parallelevaluate(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    return empty_list();
}

/* Kernels[] -- the running subkernels: none. */
static Expr* builtin_kernels(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 0) return NULL;
    return empty_list();
}

/* DistributeDefinitions[s1, ...] -- nothing to distribute to. */
static Expr* builtin_distributedefinitions(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    return empty_list();
}

/* LaunchKernels[...] -- subkernels are not available. */
static Expr* builtin_launchkernels(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    mth_message("LaunchKernels", "nopar",
                "Parallel subkernels are not available; evaluation remains serial.");
    return empty_list();
}

static void reg(const char* name, Expr* (*fn)(Expr*), unsigned attrs, const char* doc) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= attrs | ATTR_PROTECTED;
    symtab_set_docstring(name, doc);
}

void parallel_init(void) {
    reg("ParallelMap", builtin_parallelmap, 0,
        "ParallelMap[f, expr]\n"
        "\tapplies f to each element of expr.  Mathilda has no parallel subkernels, so this is\n"
        "\tMap[f, expr] (a levelspec is passed on; options such as DistributedContexts are ignored).");
    reg("ParallelTable", builtin_paralleltable, ATTR_HOLDALL,
        "ParallelTable[expr, {i, imax}, ...]\n"
        "\tbuilds a table of values of expr.  Without parallel subkernels this is Table[...];\n"
        "\toptions such as Method are ignored.");
    reg("ParallelEvaluate", builtin_parallelevaluate, ATTR_HOLDFIRST,
        "ParallelEvaluate[expr]\n"
        "\tevaluates expr on every parallel subkernel and gives the list of results.  Mathilda\n"
        "\thas no subkernels ($KernelCount is 0), so expr is not evaluated and the result is {}.");
    reg("Kernels", builtin_kernels, 0,
        "Kernels[]\n\tgives the list of running parallel subkernels: always {} in Mathilda.");
    reg("DistributeDefinitions", builtin_distributedefinitions, ATTR_HOLDALL,
        "DistributeDefinitions[s1, s2, ...]\n"
        "\tdistributes the definitions of the si to all parallel subkernels.  There are none in\n"
        "\tMathilda, so nothing is distributed and the result is {}.");
    reg("LaunchKernels", builtin_launchkernels, 0,
        "LaunchKernels[n]\n"
        "\twould launch n parallel subkernels.  Mathilda has none: it emits LaunchKernels::nopar\n"
        "\tand gives {}, and evaluation remains serial.");
    /* Option names of the parallel heads: System` symbols, inert. */
    intern_symbol("ParallelKernels");
    intern_symbol("DistributedContexts");
}
