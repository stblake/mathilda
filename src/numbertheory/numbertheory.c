/* numbertheory.c -- registration hub for the number-theory subsystem.
 *
 * The individual builtins live one-per-file under src/numbertheory/
 * (mirroring src/linalg/); cross-cutting private helpers are in
 * numbertheory_internal.h.  This file owns only numbertheory_init(), which
 * registers the builtins, sets their attributes, and installs docstrings.
 * Called from core_init() alongside the other subsystem initialisers. */

#include "numbertheory.h"
#include "symtab.h"
#include "attr.h"

void numbertheory_init(void) {
    symtab_add_builtin("GCD", builtin_gcd);
    symtab_add_builtin("LCM", builtin_lcm);
    symtab_add_builtin("ExtendedGCD", builtin_extendedgcd);
    symtab_add_builtin("PowerMod", builtin_powermod);
    symtab_add_builtin("PrimitiveRoot", builtin_primitiveroot);
    symtab_add_builtin("PrimitiveRootList", builtin_primitiverootlist);
    symtab_add_builtin("MultiplicativeOrder", builtin_multiplicativeorder);
    symtab_add_builtin("Factorial", builtin_factorial);
    symtab_add_builtin("Factorial2", builtin_factorial2);
    symtab_add_builtin("FactorialPower", builtin_factorialpower);
    symtab_add_builtin("Binomial", builtin_binomial);
    symtab_add_builtin("Divisible", builtin_divisible);
    symtab_add_builtin("CoprimeQ", builtin_coprimeq);
    symtab_add_builtin("Divisors", builtin_divisors);
    symtab_add_builtin("DivisorSigma", builtin_divisorsigma);
    symtab_add_builtin("MoebiusMu", builtin_moebiusmu);
    symtab_add_builtin("Prime", builtin_prime);
    symtab_add_builtin("PrimePi", builtin_primepi);

    symtab_get_def("GCD")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY);
    symtab_get_def("LCM")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY);
    symtab_get_def("ExtendedGCD")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    symtab_get_def("PowerMod")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_get_def("PrimitiveRoot")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_get_def("PrimitiveRootList")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_get_def("MultiplicativeOrder")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Factorial")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Factorial2")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("FactorialPower")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_set_docstring("FactorialPower",
        "FactorialPower[n, k]\n\tThe falling factorial n (n - 1) (n - 2) ... (n - k + 1).\n\tFor non-negative integer k, expands to a product of k linear factors.\n\tEquivalent to n! / (n - k)! when both n and k are non-negative integers.");
    symtab_get_def("Binomial")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Divisible")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    symtab_get_def("CoprimeQ")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE | ATTR_ORDERLESS);
    symtab_get_def("Divisors")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    symtab_get_def("DivisorSigma")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE | ATTR_NHOLDALL);
    symtab_get_def("MoebiusMu")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    symtab_get_def("Prime")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    symtab_get_def("PrimePi")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
}
