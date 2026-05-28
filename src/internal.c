#include "internal.h"
#include "eval.h"
#include "symtab.h"
#include <string.h>
#include <stdlib.h>
extern Expr* builtin_table(Expr* res);
extern Expr* builtin_range(Expr* res);
extern Expr* builtin_array(Expr* res);
extern Expr* builtin_take(Expr* res);
extern Expr* builtin_drop(Expr* res);
extern Expr* builtin_flatten(Expr* res);
extern Expr* builtin_partition(Expr* res);
extern Expr* builtin_rotateleft(Expr* res);
extern Expr* builtin_rotateright(Expr* res);
extern Expr* builtin_reverse(Expr* res);
extern Expr* builtin_transpose(Expr* res);
extern Expr* builtin_tally(Expr* res);
extern Expr* builtin_union(Expr* res);
extern Expr* builtin_deleteduplicates(Expr* res);
extern Expr* builtin_split(Expr* res);
extern Expr* builtin_total(Expr* res);
extern Expr* builtin_commonest(Expr* res);
extern Expr* builtin_min(Expr* res);
extern Expr* builtin_max(Expr* res);
extern Expr* builtin_listq(Expr* res);
extern Expr* builtin_vectorq(Expr* res);
extern Expr* builtin_matrixq(Expr* res);
extern Expr* builtin_and(Expr* res);
extern Expr* builtin_or(Expr* res);
extern Expr* builtin_not(Expr* res);
extern Expr* builtin_polynomialq(Expr* res);
extern Expr* builtin_variables(Expr* res);
extern Expr* builtin_coefficient(Expr* res);
extern Expr* builtin_coefficientlist(Expr* res);
extern Expr* builtin_polynomialgcd(Expr* res);
extern Expr* builtin_polynomiallcm(Expr* res);
extern Expr* builtin_polynomialquotient(Expr* res);
extern Expr* builtin_polynomialremainder(Expr* res);
extern Expr* builtin_polynomialmod(Expr* res);
extern Expr* builtin_collect(Expr* res);
extern Expr* builtin_decompose(Expr* res);
extern Expr* builtin_hornerform(Expr* res);
extern Expr* builtin_resultant(Expr* res);
extern Expr* builtin_polynomialextendedgcd(Expr* res);
extern Expr* builtin_discriminant(Expr* res);
extern Expr* builtin_dot(Expr* res);
extern Expr* builtin_det(Expr* res);
extern Expr* builtin_cross(Expr* res);
extern Expr* builtin_norm(Expr* res);
extern Expr* builtin_tr(Expr* res);
extern Expr* builtin_rowreduce(Expr* res);
extern Expr* builtin_identitymatrix(Expr* res);
extern Expr* builtin_diagonalmatrix(Expr* res);
extern Expr* builtin_factorsquarefree(Expr* res);
extern Expr* builtin_factor(Expr* res);
extern Expr* builtin_sameq(Expr* res);
extern Expr* builtin_unsameq(Expr* res);
extern Expr* builtin_equal(Expr* res);
extern Expr* builtin_unequal(Expr* res);
extern Expr* builtin_less(Expr* res);
extern Expr* builtin_greater(Expr* res);
extern Expr* builtin_lessequal(Expr* res);
extern Expr* builtin_greaterequal(Expr* res);
extern Expr* builtin_sinh(Expr* res);
extern Expr* builtin_cosh(Expr* res);
extern Expr* builtin_tanh(Expr* res);
extern Expr* builtin_coth(Expr* res);
extern Expr* builtin_sech(Expr* res);
extern Expr* builtin_csch(Expr* res);
extern Expr* builtin_arcsinh(Expr* res);
extern Expr* builtin_arccosh(Expr* res);
extern Expr* builtin_arctanh(Expr* res);
extern Expr* builtin_arccoth(Expr* res);
extern Expr* builtin_arcsech(Expr* res);
extern Expr* builtin_arccsch(Expr* res);
extern Expr* builtin_primeq(Expr* res);
extern Expr* builtin_primepi(Expr* res);
extern Expr* builtin_memberq(Expr* res);
extern Expr* builtin_factorinteger(Expr* res);
extern Expr* builtin_nextprime(Expr* res);
extern Expr* builtin_eulerphi(Expr* res);
extern Expr* builtin_re(Expr* res);
extern Expr* builtin_im(Expr* res);
extern Expr* builtin_reim(Expr* res);
extern Expr* builtin_abs(Expr* res);
extern Expr* builtin_conjugate(Expr* res);
extern Expr* builtin_arg(Expr* res);
extern Expr* builtin_mean(Expr* res);
extern Expr* builtin_rootmeansquare(Expr* res);
extern Expr* builtin_variance(Expr* res);
extern Expr* builtin_standard_deviation(Expr* res);
extern Expr* builtin_slot(Expr* res);
extern Expr* builtin_slotsequence(Expr* res);
extern Expr* builtin_numerator(Expr* res);
extern Expr* builtin_denominator(Expr* res);
extern Expr* builtin_cancel(Expr* res);
extern Expr* builtin_together(Expr* res);
extern Expr* builtin_apart(Expr* res);
extern Expr* builtin_atomq(Expr* res);
extern Expr* builtin_length(Expr* res);
extern Expr* builtin_dimensions(Expr* res);
extern Expr* builtin_clear(Expr* res);
extern Expr* builtin_part(Expr* res);
extern Expr* builtin_head(Expr* res);
extern Expr* builtin_first(Expr* res);
extern Expr* builtin_last(Expr* res);
extern Expr* builtin_most(Expr* res);
extern Expr* builtin_rest(Expr* res);
extern Expr* builtin_insert(Expr* res);
extern Expr* builtin_delete(Expr* res);
extern Expr* builtin_append(Expr* res);
extern Expr* builtin_prepend(Expr* res);
extern Expr* builtin_append_to(Expr* res);
extern Expr* builtin_prepend_to(Expr* res);
extern Expr* builtin_attributes(Expr* res);
extern Expr* builtin_set_attributes(Expr* res);
extern Expr* builtin_own_values(Expr* res);
extern Expr* builtin_down_values(Expr* res);
extern Expr* builtin_out(Expr* res);
extern Expr* builtin_plus(Expr* res);
extern Expr* builtin_times(Expr* res);
extern Expr* builtin_divide(Expr* res);
extern Expr* builtin_subtract(Expr* res);
extern Expr* builtin_complex(Expr* res);
extern Expr* builtin_rational(Expr* res);
extern Expr* builtin_power(Expr* res);
extern Expr* builtin_sqrt(Expr* res);
extern Expr* builtin_apply(Expr* res);
extern Expr* builtin_map(Expr* res);
extern Expr* builtin_map_all(Expr* res);
extern Expr* builtin_through(Expr* res);
extern Expr* builtin_select(Expr* res);
extern Expr* builtin_freeq(Expr* res);
extern Expr* builtin_sort(Expr* res);
extern Expr* builtin_orderedq(Expr* res);
extern Expr* builtin_level(Expr* res);
extern Expr* builtin_depth(Expr* res);
extern Expr* builtin_leafcount(Expr* res);
extern Expr* builtin_bytecount(Expr* res);
extern Expr* builtin_matchq(Expr* res);
extern Expr* builtin_compoundexpression(Expr* res);
extern Expr* builtin_numberq(Expr* res);
extern Expr* builtin_integerq(Expr* res);
extern Expr* builtin_evenq(Expr* res);
extern Expr* builtin_oddq(Expr* res);
extern Expr* builtin_mod(Expr* res);
extern Expr* builtin_quotient(Expr* res);
extern Expr* builtin_quotientremainder(Expr* res);
extern Expr* builtin_gcd(Expr* res);
extern Expr* builtin_lcm(Expr* res);
extern Expr* builtin_powermod(Expr* res);
extern Expr* builtin_factorial(Expr* res);
extern Expr* builtin_binomial(Expr* res);
extern Expr* builtin_information(Expr* res);
extern Expr* builtin_floor(Expr* res);
extern Expr* builtin_ceiling(Expr* res);
extern Expr* builtin_round(Expr* res);
extern Expr* builtin_integerpart(Expr* res);
extern Expr* builtin_fractionalpart(Expr* res);
extern Expr* builtin_replace_part(Expr* res);
extern Expr* builtin_replace(Expr* res);
extern Expr* builtin_replace_all(Expr* res);
extern Expr* builtin_replace_repeated(Expr* res);
extern Expr* builtin_replacelist(Expr* res);
extern Expr* builtin_module(Expr* res);
extern Expr* builtin_block(Expr* res);
extern Expr* builtin_with(Expr* res);
extern Expr* builtin_expand(Expr* res);
extern Expr* builtin_sin(Expr* res);
extern Expr* builtin_cos(Expr* res);
extern Expr* builtin_tan(Expr* res);
extern Expr* builtin_cot(Expr* res);
extern Expr* builtin_sec(Expr* res);
extern Expr* builtin_csc(Expr* res);
extern Expr* builtin_arcsin(Expr* res);
extern Expr* builtin_arccos(Expr* res);
extern Expr* builtin_arctan(Expr* res);
extern Expr* builtin_arccot(Expr* res);
extern Expr* builtin_arcsec(Expr* res);
extern Expr* builtin_arccsc(Expr* res);
extern Expr* builtin_log(Expr* res);
extern Expr* builtin_exp(Expr* res);
extern Expr* builtin_timing(Expr* res);
extern Expr* builtin_repeated_timing(Expr* res);
extern Expr* builtin_possible_zero_q(Expr* res);

Expr* internal_call_impl(const char* name, Expr* (*builtin_func)(Expr*), Expr** args, size_t count) {
    Expr* res = expr_new_function(expr_new_symbol(name), args, count);
    SymbolDef* def = symtab_get_def(name);
    uint32_t attrs = def ? def->attributes : 0;
    if (attrs & ATTR_FLAT) eval_flatten_args(res, name);
    if (attrs & ATTR_ORDERLESS) eval_sort_args(res);
    Expr* evaluated = builtin_func(res);
    if (evaluated) {
        expr_free(res);
        return evaluated;
    }
    return res;
}
Expr* internal_table(Expr** args, size_t count) { return internal_call_impl("Table", builtin_table, args, count); }
Expr* internal_range(Expr** args, size_t count) { return internal_call_impl("Range", builtin_range, args, count); }
Expr* internal_array(Expr** args, size_t count) { return internal_call_impl("Array", builtin_array, args, count); }
Expr* internal_take(Expr** args, size_t count) { return internal_call_impl("Take", builtin_take, args, count); }
Expr* internal_drop(Expr** args, size_t count) { return internal_call_impl("Drop", builtin_drop, args, count); }
Expr* internal_flatten(Expr** args, size_t count) { return internal_call_impl("Flatten", builtin_flatten, args, count); }
Expr* internal_partition(Expr** args, size_t count) { return internal_call_impl("Partition", builtin_partition, args, count); }
Expr* internal_rotateleft(Expr** args, size_t count) { return internal_call_impl("RotateLeft", builtin_rotateleft, args, count); }
Expr* internal_rotateright(Expr** args, size_t count) { return internal_call_impl("RotateRight", builtin_rotateright, args, count); }
Expr* internal_reverse(Expr** args, size_t count) { return internal_call_impl("Reverse", builtin_reverse, args, count); }
Expr* internal_transpose(Expr** args, size_t count) { return internal_call_impl("Transpose", builtin_transpose, args, count); }
Expr* internal_tally(Expr** args, size_t count) { return internal_call_impl("Tally", builtin_tally, args, count); }
Expr* internal_union(Expr** args, size_t count) { return internal_call_impl("Union", builtin_union, args, count); }
Expr* internal_deleteduplicates(Expr** args, size_t count) { return internal_call_impl("DeleteDuplicates", builtin_deleteduplicates, args, count); }
Expr* internal_split(Expr** args, size_t count) { return internal_call_impl("Split", builtin_split, args, count); }
Expr* internal_total(Expr** args, size_t count) { return internal_call_impl("Total", builtin_total, args, count); }
Expr* internal_commonest(Expr** args, size_t count) { return internal_call_impl("Commonest", builtin_commonest, args, count); }
Expr* internal_min(Expr** args, size_t count) { return internal_call_impl("Min", builtin_min, args, count); }
Expr* internal_max(Expr** args, size_t count) { return internal_call_impl("Max", builtin_max, args, count); }
Expr* internal_listq(Expr** args, size_t count) { return internal_call_impl("ListQ", builtin_listq, args, count); }
Expr* internal_vectorq(Expr** args, size_t count) { return internal_call_impl("VectorQ", builtin_vectorq, args, count); }
Expr* internal_matrixq(Expr** args, size_t count) { return internal_call_impl("MatrixQ", builtin_matrixq, args, count); }
Expr* internal_and(Expr** args, size_t count) { return internal_call_impl("And", builtin_and, args, count); }
Expr* internal_or(Expr** args, size_t count) { return internal_call_impl("Or", builtin_or, args, count); }
Expr* internal_not(Expr** args, size_t count) { return internal_call_impl("Not", builtin_not, args, count); }
Expr* internal_polynomialq(Expr** args, size_t count) { return internal_call_impl("PolynomialQ", builtin_polynomialq, args, count); }
Expr* internal_variables(Expr** args, size_t count) { return internal_call_impl("Variables", builtin_variables, args, count); }
Expr* internal_coefficient(Expr** args, size_t count) { return internal_call_impl("Coefficient", builtin_coefficient, args, count); }
Expr* internal_coefficientlist(Expr** args, size_t count) { return internal_call_impl("CoefficientList", builtin_coefficientlist, args, count); }
Expr* internal_polynomialgcd(Expr** args, size_t count) { return internal_call_impl("PolynomialGCD", builtin_polynomialgcd, args, count); }
Expr* internal_polynomiallcm(Expr** args, size_t count) { return internal_call_impl("PolynomialLCM", builtin_polynomiallcm, args, count); }
Expr* internal_polynomialquotient(Expr** args, size_t count) { return internal_call_impl("PolynomialQuotient", builtin_polynomialquotient, args, count); }
Expr* internal_polynomialremainder(Expr** args, size_t count) { return internal_call_impl("PolynomialRemainder", builtin_polynomialremainder, args, count); }
Expr* internal_polynomialmod(Expr** args, size_t count) { return internal_call_impl("PolynomialMod", builtin_polynomialmod, args, count); }
Expr* internal_collect(Expr** args, size_t count) { return internal_call_impl("Collect", builtin_collect, args, count); }
Expr* internal_decompose(Expr** args, size_t count) { return internal_call_impl("Decompose", builtin_decompose, args, count); }
Expr* internal_hornerform(Expr** args, size_t count) { return internal_call_impl("HornerForm", builtin_hornerform, args, count); }
Expr* internal_resultant(Expr** args, size_t count) { return internal_call_impl("Resultant", builtin_resultant, args, count); }
Expr* internal_polynomialextendedgcd(Expr** args, size_t count) { return internal_call_impl("PolynomialExtendedGCD", builtin_polynomialextendedgcd, args, count); }
Expr* internal_discriminant(Expr** args, size_t count) { return internal_call_impl("Discriminant", builtin_discriminant, args, count); }
Expr* internal_dot(Expr** args, size_t count) { return internal_call_impl("Dot", builtin_dot, args, count); }
Expr* internal_det(Expr** args, size_t count) { return internal_call_impl("Det", builtin_det, args, count); }
Expr* internal_cross(Expr** args, size_t count) { return internal_call_impl("Cross", builtin_cross, args, count); }
Expr* internal_norm(Expr** args, size_t count) { return internal_call_impl("Norm", builtin_norm, args, count); }
Expr* internal_tr(Expr** args, size_t count) { return internal_call_impl("Tr", builtin_tr, args, count); }
Expr* internal_rowreduce(Expr** args, size_t count) { return internal_call_impl("RowReduce", builtin_rowreduce, args, count); }
Expr* internal_identitymatrix(Expr** args, size_t count) { return internal_call_impl("IdentityMatrix", builtin_identitymatrix, args, count); }
Expr* internal_diagonalmatrix(Expr** args, size_t count) { return internal_call_impl("DiagonalMatrix", builtin_diagonalmatrix, args, count); }
Expr* internal_factorsquarefree(Expr** args, size_t count) { return internal_call_impl("FactorSquareFree", builtin_factorsquarefree, args, count); }
Expr* internal_factor(Expr** args, size_t count) { return internal_call_impl("Factor", builtin_factor, args, count); }
Expr* internal_sameq(Expr** args, size_t count) { return internal_call_impl("SameQ", builtin_sameq, args, count); }
Expr* internal_unsameq(Expr** args, size_t count) { return internal_call_impl("UnsameQ", builtin_unsameq, args, count); }
Expr* internal_equal(Expr** args, size_t count) { return internal_call_impl("Equal", builtin_equal, args, count); }
Expr* internal_unequal(Expr** args, size_t count) { return internal_call_impl("Unequal", builtin_unequal, args, count); }
Expr* internal_less(Expr** args, size_t count) { return internal_call_impl("Less", builtin_less, args, count); }
Expr* internal_greater(Expr** args, size_t count) { return internal_call_impl("Greater", builtin_greater, args, count); }
Expr* internal_lessequal(Expr** args, size_t count) { return internal_call_impl("LessEqual", builtin_lessequal, args, count); }
Expr* internal_greaterequal(Expr** args, size_t count) { return internal_call_impl("GreaterEqual", builtin_greaterequal, args, count); }
Expr* internal_sinh(Expr** args, size_t count) { return internal_call_impl("Sinh", builtin_sinh, args, count); }
Expr* internal_cosh(Expr** args, size_t count) { return internal_call_impl("Cosh", builtin_cosh, args, count); }
Expr* internal_tanh(Expr** args, size_t count) { return internal_call_impl("Tanh", builtin_tanh, args, count); }
Expr* internal_coth(Expr** args, size_t count) { return internal_call_impl("Coth", builtin_coth, args, count); }
Expr* internal_sech(Expr** args, size_t count) { return internal_call_impl("Sech", builtin_sech, args, count); }
Expr* internal_csch(Expr** args, size_t count) { return internal_call_impl("Csch", builtin_csch, args, count); }
Expr* internal_arcsinh(Expr** args, size_t count) { return internal_call_impl("ArcSinh", builtin_arcsinh, args, count); }
Expr* internal_arccosh(Expr** args, size_t count) { return internal_call_impl("ArcCosh", builtin_arccosh, args, count); }
Expr* internal_arctanh(Expr** args, size_t count) { return internal_call_impl("ArcTanh", builtin_arctanh, args, count); }
Expr* internal_arccoth(Expr** args, size_t count) { return internal_call_impl("ArcCoth", builtin_arccoth, args, count); }
Expr* internal_arcsech(Expr** args, size_t count) { return internal_call_impl("ArcSech", builtin_arcsech, args, count); }
Expr* internal_arccsch(Expr** args, size_t count) { return internal_call_impl("ArcCsch", builtin_arccsch, args, count); }
Expr* internal_primeq(Expr** args, size_t count) { return internal_call_impl("PrimeQ", builtin_primeq, args, count); }
Expr* internal_primepi(Expr** args, size_t count) { return internal_call_impl("PrimePi", builtin_primepi, args, count); }
Expr* internal_memberq(Expr** args, size_t count) { return internal_call_impl("MemberQ", builtin_memberq, args, count); }
Expr* internal_factorinteger(Expr** args, size_t count) { return internal_call_impl("FactorInteger", builtin_factorinteger, args, count); }
Expr* internal_nextprime(Expr** args, size_t count) { return internal_call_impl("NextPrime", builtin_nextprime, args, count); }
Expr* internal_eulerphi(Expr** args, size_t count) { return internal_call_impl("EulerPhi", builtin_eulerphi, args, count); }
Expr* internal_re(Expr** args, size_t count) { return internal_call_impl("Re", builtin_re, args, count); }
Expr* internal_im(Expr** args, size_t count) { return internal_call_impl("Im", builtin_im, args, count); }
Expr* internal_reim(Expr** args, size_t count) { return internal_call_impl("ReIm", builtin_reim, args, count); }
Expr* internal_abs(Expr** args, size_t count) { return internal_call_impl("Abs", builtin_abs, args, count); }
Expr* internal_conjugate(Expr** args, size_t count) { return internal_call_impl("Conjugate", builtin_conjugate, args, count); }
Expr* internal_arg(Expr** args, size_t count) { return internal_call_impl("Arg", builtin_arg, args, count); }
Expr* internal_mean(Expr** args, size_t count) { return internal_call_impl("Mean", builtin_mean, args, count); }
Expr* internal_rootmeansquare(Expr** args, size_t count) { return internal_call_impl("RootMeanSquare", builtin_rootmeansquare, args, count); }
Expr* internal_variance(Expr** args, size_t count) { return internal_call_impl("Variance", builtin_variance, args, count); }
Expr* internal_standard_deviation(Expr** args, size_t count) { return internal_call_impl("StandardDeviation", builtin_standard_deviation, args, count); }
Expr* internal_slot(Expr** args, size_t count) { return internal_call_impl("Slot", builtin_slot, args, count); }
Expr* internal_slotsequence(Expr** args, size_t count) { return internal_call_impl("SlotSequence", builtin_slotsequence, args, count); }
Expr* internal_numerator(Expr** args, size_t count) { return internal_call_impl("Numerator", builtin_numerator, args, count); }
Expr* internal_denominator(Expr** args, size_t count) { return internal_call_impl("Denominator", builtin_denominator, args, count); }
Expr* internal_cancel(Expr** args, size_t count) { return internal_call_impl("Cancel", builtin_cancel, args, count); }
Expr* internal_together(Expr** args, size_t count) { return internal_call_impl("Together", builtin_together, args, count); }
Expr* internal_apart(Expr** args, size_t count) { return internal_call_impl("Apart", builtin_apart, args, count); }
Expr* internal_atomq(Expr** args, size_t count) { return internal_call_impl("AtomQ", builtin_atomq, args, count); }
Expr* internal_length(Expr** args, size_t count) { return internal_call_impl("Length", builtin_length, args, count); }
Expr* internal_dimensions(Expr** args, size_t count) { return internal_call_impl("Dimensions", builtin_dimensions, args, count); }
Expr* internal_clear(Expr** args, size_t count) { return internal_call_impl("Clear", builtin_clear, args, count); }
Expr* internal_part(Expr** args, size_t count) { return internal_call_impl("Part", builtin_part, args, count); }
Expr* internal_head(Expr** args, size_t count) { return internal_call_impl("Head", builtin_head, args, count); }
Expr* internal_first(Expr** args, size_t count) { return internal_call_impl("First", builtin_first, args, count); }
Expr* internal_last(Expr** args, size_t count) { return internal_call_impl("Last", builtin_last, args, count); }
Expr* internal_most(Expr** args, size_t count) { return internal_call_impl("Most", builtin_most, args, count); }
Expr* internal_rest(Expr** args, size_t count) { return internal_call_impl("Rest", builtin_rest, args, count); }
Expr* internal_insert(Expr** args, size_t count) { return internal_call_impl("Insert", builtin_insert, args, count); }
Expr* internal_delete(Expr** args, size_t count) { return internal_call_impl("Delete", builtin_delete, args, count); }
Expr* internal_append(Expr** args, size_t count) { return internal_call_impl("Append", builtin_append, args, count); }
Expr* internal_prepend(Expr** args, size_t count) { return internal_call_impl("Prepend", builtin_prepend, args, count); }
Expr* internal_append_to(Expr** args, size_t count) { return internal_call_impl("AppendTo", builtin_append_to, args, count); }
Expr* internal_prepend_to(Expr** args, size_t count) { return internal_call_impl("PrependTo", builtin_prepend_to, args, count); }
Expr* internal_attributes(Expr** args, size_t count) { return internal_call_impl("Attributes", builtin_attributes, args, count); }
Expr* internal_set_attributes(Expr** args, size_t count) { return internal_call_impl("SetAttributes", builtin_set_attributes, args, count); }
Expr* internal_own_values(Expr** args, size_t count) { return internal_call_impl("OwnValues", builtin_own_values, args, count); }
Expr* internal_down_values(Expr** args, size_t count) { return internal_call_impl("DownValues", builtin_down_values, args, count); }
Expr* internal_out(Expr** args, size_t count) { return internal_call_impl("Out", builtin_out, args, count); }
Expr* internal_plus(Expr** args, size_t count) { return internal_call_impl("Plus", builtin_plus, args, count); }
Expr* internal_times(Expr** args, size_t count) { return internal_call_impl("Times", builtin_times, args, count); }
Expr* internal_divide(Expr** args, size_t count) { return internal_call_impl("Divide", builtin_divide, args, count); }
Expr* internal_subtract(Expr** args, size_t count) { return internal_call_impl("Subtract", builtin_subtract, args, count); }
Expr* internal_complex(Expr** args, size_t count) { return internal_call_impl("Complex", builtin_complex, args, count); }
Expr* internal_rational(Expr** args, size_t count) { return internal_call_impl("Rational", builtin_rational, args, count); }
Expr* internal_power(Expr** args, size_t count) { return internal_call_impl("Power", builtin_power, args, count); }
Expr* internal_sqrt(Expr** args, size_t count) { return internal_call_impl("Sqrt", builtin_sqrt, args, count); }
Expr* internal_apply(Expr** args, size_t count) { return internal_call_impl("Apply", builtin_apply, args, count); }
Expr* internal_map(Expr** args, size_t count) { return internal_call_impl("Map", builtin_map, args, count); }
Expr* internal_map_all(Expr** args, size_t count) { return internal_call_impl("MapAll", builtin_map_all, args, count); }
Expr* internal_through(Expr** args, size_t count) { return internal_call_impl("Through", builtin_through, args, count); }
Expr* internal_select(Expr** args, size_t count) { return internal_call_impl("Select", builtin_select, args, count); }
Expr* internal_freeq(Expr** args, size_t count) { return internal_call_impl("FreeQ", builtin_freeq, args, count); }
Expr* internal_sort(Expr** args, size_t count) { return internal_call_impl("Sort", builtin_sort, args, count); }
Expr* internal_orderedq(Expr** args, size_t count) { return internal_call_impl("OrderedQ", builtin_orderedq, args, count); }
Expr* internal_level(Expr** args, size_t count) { return internal_call_impl("Level", builtin_level, args, count); }
Expr* internal_depth(Expr** args, size_t count) { return internal_call_impl("Depth", builtin_depth, args, count); }
Expr* internal_leafcount(Expr** args, size_t count) { return internal_call_impl("LeafCount", builtin_leafcount, args, count); }
Expr* internal_bytecount(Expr** args, size_t count) { return internal_call_impl("ByteCount", builtin_bytecount, args, count); }
Expr* internal_matchq(Expr** args, size_t count) { return internal_call_impl("MatchQ", builtin_matchq, args, count); }
Expr* internal_compoundexpression(Expr** args, size_t count) { return internal_call_impl("CompoundExpression", builtin_compoundexpression, args, count); }
Expr* internal_numberq(Expr** args, size_t count) { return internal_call_impl("NumberQ", builtin_numberq, args, count); }
Expr* internal_integerq(Expr** args, size_t count) { return internal_call_impl("IntegerQ", builtin_integerq, args, count); }
Expr* internal_evenq(Expr** args, size_t count) { return internal_call_impl("EvenQ", builtin_evenq, args, count); }
Expr* internal_oddq(Expr** args, size_t count) { return internal_call_impl("OddQ", builtin_oddq, args, count); }
Expr* internal_mod(Expr** args, size_t count) { return internal_call_impl("Mod", builtin_mod, args, count); }
Expr* internal_quotient(Expr** args, size_t count) { return internal_call_impl("Quotient", builtin_quotient, args, count); }
Expr* internal_quotientremainder(Expr** args, size_t count) { return internal_call_impl("QuotientRemainder", builtin_quotientremainder, args, count); }
Expr* internal_gcd(Expr** args, size_t count) { return internal_call_impl("GCD", builtin_gcd, args, count); }
Expr* internal_lcm(Expr** args, size_t count) { return internal_call_impl("LCM", builtin_lcm, args, count); }
Expr* internal_powermod(Expr** args, size_t count) { return internal_call_impl("PowerMod", builtin_powermod, args, count); }
Expr* internal_factorial(Expr** args, size_t count) { return internal_call_impl("Factorial", builtin_factorial, args, count); }
Expr* internal_binomial(Expr** args, size_t count) { return internal_call_impl("Binomial", builtin_binomial, args, count); }
Expr* internal_information(Expr** args, size_t count) { return internal_call_impl("Information", builtin_information, args, count); }
Expr* internal_floor(Expr** args, size_t count) { return internal_call_impl("Floor", builtin_floor, args, count); }
Expr* internal_ceiling(Expr** args, size_t count) { return internal_call_impl("Ceiling", builtin_ceiling, args, count); }
Expr* internal_round(Expr** args, size_t count) { return internal_call_impl("Round", builtin_round, args, count); }
Expr* internal_integerpart(Expr** args, size_t count) { return internal_call_impl("IntegerPart", builtin_integerpart, args, count); }
Expr* internal_fractionalpart(Expr** args, size_t count) { return internal_call_impl("FractionalPart", builtin_fractionalpart, args, count); }
Expr* internal_replace_part(Expr** args, size_t count) { return internal_call_impl("ReplacePart", builtin_replace_part, args, count); }
Expr* internal_replace(Expr** args, size_t count) { return internal_call_impl("Replace", builtin_replace, args, count); }
Expr* internal_replace_all(Expr** args, size_t count) { return internal_call_impl("ReplaceAll", builtin_replace_all, args, count); }
Expr* internal_replace_repeated(Expr** args, size_t count) { return internal_call_impl("ReplaceRepeated", builtin_replace_repeated, args, count); }
Expr* internal_replacelist(Expr** args, size_t count) { return internal_call_impl("ReplaceList", builtin_replacelist, args, count); }
Expr* internal_module(Expr** args, size_t count) { return internal_call_impl("Module", builtin_module, args, count); }
Expr* internal_block(Expr** args, size_t count) { return internal_call_impl("Block", builtin_block, args, count); }
Expr* internal_with(Expr** args, size_t count) { return internal_call_impl("With", builtin_with, args, count); }
Expr* internal_expand(Expr** args, size_t count) { return internal_call_impl("Expand", builtin_expand, args, count); }
Expr* internal_sin(Expr** args, size_t count) { return internal_call_impl("Sin", builtin_sin, args, count); }
Expr* internal_cos(Expr** args, size_t count) { return internal_call_impl("Cos", builtin_cos, args, count); }
Expr* internal_tan(Expr** args, size_t count) { return internal_call_impl("Tan", builtin_tan, args, count); }
Expr* internal_cot(Expr** args, size_t count) { return internal_call_impl("Cot", builtin_cot, args, count); }
Expr* internal_sec(Expr** args, size_t count) { return internal_call_impl("Sec", builtin_sec, args, count); }
Expr* internal_csc(Expr** args, size_t count) { return internal_call_impl("Csc", builtin_csc, args, count); }
Expr* internal_arcsin(Expr** args, size_t count) { return internal_call_impl("ArcSin", builtin_arcsin, args, count); }
Expr* internal_arccos(Expr** args, size_t count) { return internal_call_impl("ArcCos", builtin_arccos, args, count); }
Expr* internal_arctan(Expr** args, size_t count) { return internal_call_impl("ArcTan", builtin_arctan, args, count); }
Expr* internal_arccot(Expr** args, size_t count) { return internal_call_impl("ArcCot", builtin_arccot, args, count); }
Expr* internal_arcsec(Expr** args, size_t count) { return internal_call_impl("ArcSec", builtin_arcsec, args, count); }
Expr* internal_arccsc(Expr** args, size_t count) { return internal_call_impl("ArcCsc", builtin_arccsc, args, count); }
Expr* internal_log(Expr** args, size_t count) { return internal_call_impl("Log", builtin_log, args, count); }
Expr* internal_exp(Expr** args, size_t count) { return internal_call_impl("Exp", builtin_exp, args, count); }
Expr* internal_timing(Expr** args, size_t count) { return internal_call_impl("Timing", builtin_timing, args, count); }
Expr* internal_repeated_timing(Expr** args, size_t count) { return internal_call_impl("RepeatedTiming", builtin_repeated_timing, args, count); }
Expr* internal_possible_zero_q(Expr** args, size_t count) { return internal_call_impl("PossibleZeroQ", builtin_possible_zero_q, args, count); }
