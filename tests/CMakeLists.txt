cmake_minimum_required(VERSION 3.10)
project(Mathilda_tests)

option(USE_ECM "Use ECM library for factorisation" ON)
option(USE_MPFR "Use MPFR for arbitrary-precision reals (N[expr, prec])" ON)

if(USE_MPFR)
    add_compile_definitions(USE_MPFR)
    link_libraries(mpfr)
endif()

include_directories(/usr/local/include)
if(USE_ECM)
    include_directories(../src/external/ecm)
endif()

if(EXISTS "/opt/homebrew/include")
    include_directories(/opt/homebrew/include)
endif()
link_directories(/usr/local/lib)
if(EXISTS "/opt/homebrew/lib")
    link_directories(/opt/homebrew/lib)
endif()

if(USE_ECM)
    add_custom_command(
        OUTPUT ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a
        COMMAND ./build_ecm.sh
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/..
        COMMENT "Building ECM library"
    )
    add_custom_target(build_ecm DEPENDS ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    link_libraries(gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    add_compile_definitions(NO_ECM)
    link_libraries(gmp m)
endif()

set(COMMON_SRC
    ../src/expr.c
    ../src/parse.c
    ../src/eval.c
    ../src/internal.c
    ../src/attr.c
    ../src/core.c
    ../src/comparisons.c
    ../src/symtab.c
    ../src/match.c
    ../src/part.c
    ../src/print.c
    ../src/plus.c
    ../src/times.c
    ../src/arithmetic.c
    ../src/facint.c
    ../src/power.c
    ../src/funcprog.c
    ../src/boolean.c
    ../src/list.c
    ../src/patterns.c
    ../src/patterns.c
    ../src/replace.c
    ${CMAKE_SOURCE_DIR}/../src/trigsimp.c
    ${CMAKE_SOURCE_DIR}/../src/trigrat.c
    ${CMAKE_SOURCE_DIR}/../src/simp.c
    ../src/complex.c
    ../src/trig.c
    ../src/trig_canon.c
    ../src/hyperbolic.c
    ../src/logexp.c
    ../src/piecewise.c
    ../src/purefunc.c
    ../src/modular.c
    ../src/sort.c
    ../src/stats.c
    ../src/info.c
    ../src/expand.c
    ../src/poly.c
    ../src/rat.c
    ../src/facpoly.c
    ../src/parfrac.c
    ../src/datetime.c
    ../src/linalg.c
    ../src/matsol.c
    ../src/matinv.c
    ../src/readwrite.c
    ../src/files.c
    ../src/iter.c
    ../src/cond.c
    ../src/default_helper.c
    ../src/random.c
    ../src/picostrings.c
    ../src/series.c
    ../src/deriv.c
    ../src/limit.c
    ../src/context.c
    ../src/numeric.c
    ../src/precision.c
    ../src/rationalize.c
    ../src/common.c
    ../src/zupoly.c
    ../src/bpoly.c
    ../src/mpoly.c
    ../src/mvfactor.c
    ../src/mvfactor3.c
    ../src/qa.c
    ../src/qaupoly.c
    ../src/qafactor.c
    ../src/sym_intern.c
    ../src/sym_names.c
    ../src/repl_hooks.c
    ../src/options.c
    ../src/integrate.c
    ../src/intrat.c
    ../src/intrischnorman.c
    ../src/intsimp.c
    ../src/root.c
    ../src/solve.c
    ../src/solvepoly.c
    ../src/solvelinsys.c
)

add_executable(expr_tests test_expr.c test_part.c ${COMMON_SRC})
target_include_directories(expr_tests PRIVATE ../src)

add_executable(parse_tests test_parse.c ${COMMON_SRC})
target_include_directories(parse_tests PRIVATE ../src)

add_executable(eval_tests test_eval.c ${COMMON_SRC})
target_include_directories(eval_tests PRIVATE ../src)

add_executable(match_tests test_match.c ${COMMON_SRC})
target_include_directories(match_tests PRIVATE ../src)

add_executable(match_extensive_tests test_match_extensive.c ${COMMON_SRC})
target_include_directories(match_extensive_tests PRIVATE ../src)

add_executable(symtab_tests test_symtab.c ${COMMON_SRC})
target_include_directories(symtab_tests PRIVATE ../src)

add_executable(regression_tests test_regression.c ${COMMON_SRC})
target_include_directories(regression_tests PRIVATE ../src)

add_executable(comparisons_tests test_comparisons.c ${COMMON_SRC})
target_include_directories(comparisons_tests PRIVATE ../src)

add_executable(list_tests test_list.c ${COMMON_SRC})
target_include_directories(list_tests PRIVATE ../src)

add_executable(iter_tests test_iter.c ${COMMON_SRC})
target_link_libraries(iter_tests m)
target_include_directories(iter_tests PRIVATE ../src)

add_executable(cond_tests test_cond.c ${COMMON_SRC})
target_link_libraries(cond_tests m)
target_include_directories(cond_tests PRIVATE ../src)

add_executable(patterns_tests test_patterns.c ${COMMON_SRC})
target_link_libraries(patterns_tests m)
target_include_directories(patterns_tests PRIVATE ../src)

add_executable(replace_tests test_replace.c ${COMMON_SRC})
target_include_directories(replace_tests PRIVATE ../src)

add_executable(core_tests test_core.c ${COMMON_SRC})
target_include_directories(core_tests PRIVATE ../src)

add_executable(chop_tests test_chop.c ${COMMON_SRC})
target_include_directories(chop_tests PRIVATE ../src)

add_executable(clip_tests test_clip.c ${COMMON_SRC})
target_include_directories(clip_tests PRIVATE ../src)

add_executable(symbol_tests test_symbol.c ${COMMON_SRC})
target_include_directories(symbol_tests PRIVATE ../src)

add_executable(tostring_toexpression_tests test_tostring_toexpression.c ${COMMON_SRC})
target_include_directories(tostring_toexpression_tests PRIVATE ../src)

add_executable(trig_tests test_trig.c ${COMMON_SRC})
target_include_directories(trig_tests PRIVATE ../src)

add_executable(trig_canon_tests test_trig_canon.c ${COMMON_SRC})
target_include_directories(trig_canon_tests PRIVATE ../src)

add_executable(hyperbolic_tests test_hyperbolic.c ${COMMON_SRC})
target_include_directories(hyperbolic_tests PRIVATE ../src)

add_executable(logexp_tests test_logexp.c ${COMMON_SRC})
target_include_directories(logexp_tests PRIVATE ../src)

add_executable(piecewise_tests test_piecewise.c ${COMMON_SRC})
target_include_directories(piecewise_tests PRIVATE ../src)

add_executable(purefunc_tests test_purefunc.c ${COMMON_SRC})
target_include_directories(purefunc_tests PRIVATE ../src)

add_executable(modular_tests test_modular.c ${COMMON_SRC})
target_include_directories(modular_tests PRIVATE ../src)

add_executable(funcprog_through_tests test_funcprog_through.c ${COMMON_SRC})
target_include_directories(funcprog_through_tests PRIVATE ../src)

add_executable(mapat_tests test_mapat.c ${COMMON_SRC})
target_include_directories(mapat_tests PRIVATE ../src)

add_executable(replaceat_tests test_replaceat.c ${COMMON_SRC})
target_include_directories(replaceat_tests PRIVATE ../src)

add_executable(nest_tests test_nest.c ${COMMON_SRC})
target_include_directories(nest_tests PRIVATE ../src)

add_executable(nestlist_tests test_nestlist.c ${COMMON_SRC})
target_include_directories(nestlist_tests PRIVATE ../src)

add_executable(nestwhile_tests test_nestwhile.c ${COMMON_SRC})
target_include_directories(nestwhile_tests PRIVATE ../src)

add_executable(nestwhilelist_tests test_nestwhilelist.c ${COMMON_SRC})
target_include_directories(nestwhilelist_tests PRIVATE ../src)

add_executable(fixedpointlist_tests test_fixedpointlist.c ${COMMON_SRC})
target_include_directories(fixedpointlist_tests PRIVATE ../src)

add_executable(fixedpoint_tests test_fixedpoint.c ${COMMON_SRC})
target_include_directories(fixedpoint_tests PRIVATE ../src)

add_executable(fold_tests test_fold.c ${COMMON_SRC})
target_include_directories(fold_tests PRIVATE ../src)

add_executable(foldlist_tests test_foldlist.c ${COMMON_SRC})
target_include_directories(foldlist_tests PRIVATE ../src)

add_executable(split_tests test_split.c ${COMMON_SRC})
target_include_directories(split_tests PRIVATE ../src)

add_executable(sort_tests test_sort.c ${COMMON_SRC})
target_include_directories(sort_tests PRIVATE ../src)

add_executable(backtrack_tests test_backtrack.c ${COMMON_SRC})
target_include_directories(backtrack_tests PRIVATE ../src)

add_executable(poly_tests test_poly.c ${COMMON_SRC})
target_include_directories(poly_tests PRIVATE ../src)

add_executable(expand_tests test_expand.c ${COMMON_SRC})
target_include_directories(expand_tests PRIVATE ../src)

add_executable(expandfrac_tests test_expandfrac.c ${COMMON_SRC})
target_include_directories(expandfrac_tests PRIVATE ../src)

add_executable(stats_tests test_stats.c ${COMMON_SRC})
target_include_directories(stats_tests PRIVATE ../src)

add_executable(datetime_tests test_datetime.c ${COMMON_SRC})
target_include_directories(datetime_tests PRIVATE ../src)

add_executable(facpoly_tests test_facpoly.c ${COMMON_SRC})
target_include_directories(facpoly_tests PRIVATE ../src)

add_executable(solve_tests test_solve.c ${COMMON_SRC})
target_include_directories(solve_tests PRIVATE ../src)

add_executable(factor_phase0_tests test_factor_phase0.c ${COMMON_SRC})
target_include_directories(factor_phase0_tests PRIVATE ../src)

add_executable(factor_baseline_tests test_factor_baseline.c ${COMMON_SRC})
target_include_directories(factor_baseline_tests PRIVATE ../src)

add_executable(factor_recombine_tests test_factor_recombine.c ${COMMON_SRC})
target_include_directories(factor_recombine_tests PRIVATE ../src)

add_executable(zupoly_tests test_zupoly.c ${COMMON_SRC})
target_include_directories(zupoly_tests PRIVATE ../src /usr/local/include)
target_link_directories(zupoly_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(zupoly_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(zupoly_tests gmp m)
endif()

# Phase G1: Q(α) algebraic-number type (Trager's algebraic factoring).
add_executable(qa_tests test_qa.c ${COMMON_SRC})
target_include_directories(qa_tests PRIVATE ../src /usr/local/include)
target_link_directories(qa_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(qa_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(qa_tests gmp m)
endif()

# Phase G2: Q(α)[x] univariate polynomial type.
add_executable(qaupoly_tests test_qaupoly.c ${COMMON_SRC})
target_include_directories(qaupoly_tests PRIVATE ../src /usr/local/include)
target_link_directories(qaupoly_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(qaupoly_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(qaupoly_tests gmp m)
endif()

# Phase G3: Norm via resultant (qaupoly_norm).
add_executable(qafactor_tests test_qafactor.c ${COMMON_SRC})
target_include_directories(qafactor_tests PRIVATE ../src /usr/local/include)
target_link_directories(qafactor_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(qafactor_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(qafactor_tests build_ecm)
else()
    target_link_libraries(qafactor_tests gmp m)
endif()

# Phase G9: extension_autodetect — walk an expression and build its
# algebraic-extension compositum automatically.
add_executable(extension_autodetect_tests test_extension_autodetect.c ${COMMON_SRC})
target_include_directories(extension_autodetect_tests PRIVATE ../src /usr/local/include)
target_link_directories(extension_autodetect_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(extension_autodetect_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(extension_autodetect_tests build_ecm)
else()
    target_link_libraries(extension_autodetect_tests gmp m)
endif()

# Phase G9: builtin integration tests for Extension -> Automatic on
# Cancel, Together, Apart, Factor, PolynomialGCD/LCM/Quotient/Remainder.
add_executable(extension_auto_builtins_tests test_extension_auto_builtins.c ${COMMON_SRC})
target_include_directories(extension_auto_builtins_tests PRIVATE ../src /usr/local/include)
target_link_directories(extension_auto_builtins_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(extension_auto_builtins_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(extension_auto_builtins_tests build_ecm)
else()
    target_link_libraries(extension_auto_builtins_tests gmp m)
endif()

# Phase G9 / Phase C: Simplify with cube-root and higher-order
# algebraic generators (q > 2 inside simp_algebraic).
add_executable(simp_algebraic_cuberoot_tests test_simp_algebraic_cuberoot.c ${COMMON_SRC})
target_include_directories(simp_algebraic_cuberoot_tests PRIVATE ../src /usr/local/include)
target_link_directories(simp_algebraic_cuberoot_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(simp_algebraic_cuberoot_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(simp_algebraic_cuberoot_tests build_ecm)
else()
    target_link_libraries(simp_algebraic_cuberoot_tests gmp m)
endif()

add_executable(bpoly_tests test_bpoly.c ${COMMON_SRC})
target_include_directories(bpoly_tests PRIVATE ../src /usr/local/include)
target_link_directories(bpoly_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(bpoly_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(bpoly_tests gmp m)
endif()

add_executable(mpoly_tests test_mpoly.c ${COMMON_SRC})
target_include_directories(mpoly_tests PRIVATE ../src /usr/local/include)
target_link_directories(mpoly_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(mpoly_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(mpoly_tests gmp m)
endif()

add_executable(mvfactor3_tests test_mvfactor3.c ${COMMON_SRC})
target_include_directories(mvfactor3_tests PRIVATE ../src /usr/local/include)
target_link_directories(mvfactor3_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(mvfactor3_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(mvfactor3_tests gmp m)
endif()

add_executable(mvfactor_tests test_mvfactor.c ${COMMON_SRC})
target_include_directories(mvfactor_tests PRIVATE ../src /usr/local/include)
target_link_directories(mvfactor_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(mvfactor_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(mvfactor_tests gmp m)
endif()

add_executable(factor_terms_tests test_factor_terms.c ${COMMON_SRC})
target_include_directories(factor_terms_tests PRIVATE ../src)

add_executable(rat_tests test_rat.c ${COMMON_SRC})
target_include_directories(rat_tests PRIVATE ../src)

add_executable(linalg_tests test_linalg.c ${COMMON_SRC})
target_include_directories(linalg_tests PRIVATE ../src)

add_executable(eigen_tests test_eigen.c ${COMMON_SRC})
target_include_directories(eigen_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(eigen_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(eigen_tests build_ecm)
else()
    target_link_libraries(eigen_tests gmp m)
endif()
add_test(NAME eigen_tests COMMAND eigen_tests)

add_executable(linearsolve_tests test_linearsolve.c ${COMMON_SRC})
target_include_directories(linearsolve_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(linearsolve_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(linearsolve_tests build_ecm)
else()
    target_link_libraries(linearsolve_tests gmp m)
endif()
add_test(NAME linearsolve_tests COMMAND linearsolve_tests)

add_executable(matsol_methods_tests test_matsol_methods.c ${COMMON_SRC})
target_include_directories(matsol_methods_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(matsol_methods_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(matsol_methods_tests build_ecm)
else()
    target_link_libraries(matsol_methods_tests gmp m)
endif()
add_test(NAME matsol_methods_tests COMMAND matsol_methods_tests)

add_executable(matinv_tests test_matinv.c ${COMMON_SRC})
target_include_directories(matinv_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(matinv_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(matinv_tests build_ecm)
else()
    target_link_libraries(matinv_tests gmp m)
endif()
add_test(NAME matinv_tests COMMAND matinv_tests)

add_executable(matinv_methods_tests test_matinv_methods.c ${COMMON_SRC})
target_include_directories(matinv_methods_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(matinv_methods_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(matinv_methods_tests build_ecm)
else()
    target_link_libraries(matinv_methods_tests gmp m)
endif()
add_test(NAME matinv_methods_tests COMMAND matinv_methods_tests)

add_executable(parfrac_tests test_parfrac.c ${COMMON_SRC})
target_include_directories(parfrac_tests PRIVATE ../src)

add_executable(polymod_tests test_polymod.c ${COMMON_SRC})
target_include_directories(polymod_tests PRIVATE ../src)

add_executable(bigint_tests test_bigint.c ${COMMON_SRC})
target_include_directories(bigint_tests PRIVATE ../src /usr/local/include)
target_link_directories(bigint_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(bigint_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(bigint_tests gmp m)
endif()

add_executable(expr_sharing_tests test_expr_sharing.c ${COMMON_SRC})
target_include_directories(expr_sharing_tests PRIVATE ../src /usr/local/include)
target_link_directories(expr_sharing_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(expr_sharing_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(expr_sharing_tests gmp m)
endif()

add_executable(eval_timestamps_tests test_eval_timestamps.c ${COMMON_SRC})
target_include_directories(eval_timestamps_tests PRIVATE ../src /usr/local/include)
target_link_directories(eval_timestamps_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(eval_timestamps_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(eval_timestamps_tests gmp m)
endif()

add_executable(eval_eager_exit_tests test_eval_eager_exit.c ${COMMON_SRC})
target_include_directories(eval_eager_exit_tests PRIVATE ../src /usr/local/include)
target_link_directories(eval_eager_exit_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(eval_eager_exit_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(eval_eager_exit_tests gmp m)
endif()

add_executable(rule_dispatch_tests test_rule_dispatch.c ${COMMON_SRC})
target_include_directories(rule_dispatch_tests PRIVATE ../src /usr/local/include)
target_link_directories(rule_dispatch_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(rule_dispatch_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(rule_dispatch_tests gmp m)
endif()

add_executable(integrals_tests test_integrals.c ${COMMON_SRC})
target_include_directories(integrals_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrals_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(integrals_tests gmp m)
endif()

add_executable(integrate_dispatch_tests test_integrate_dispatch.c ${COMMON_SRC})
target_include_directories(integrate_dispatch_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_dispatch_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(integrate_dispatch_tests gmp m)
endif()

add_executable(condition_downvalue_tests test_condition_downvalue.c ${COMMON_SRC})
target_include_directories(condition_downvalue_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(condition_downvalue_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(condition_downvalue_tests gmp m)
endif()

add_executable(print_tests test_print.c ${COMMON_SRC})
target_include_directories(print_tests PRIVATE ../src)

add_executable(context_tests test_context.c ${COMMON_SRC})
target_include_directories(context_tests PRIVATE ../src)

add_executable(list_set_tests test_list_set.c ${COMMON_SRC})
target_include_directories(list_set_tests PRIVATE ../src)

add_executable(random_tests test_random.c ${COMMON_SRC})
target_include_directories(random_tests PRIVATE ../src)

add_executable(evaluate_tests test_evaluate.c ${COMMON_SRC})
target_include_directories(evaluate_tests PRIVATE ../src)

add_executable(releasehold_tests test_releasehold.c ${COMMON_SRC})
target_include_directories(releasehold_tests PRIVATE ../src)

add_executable(unevaluated_tests test_unevaluated.c ${COMMON_SRC})
target_include_directories(unevaluated_tests PRIVATE ../src)

add_executable(strings_tests test_strings.c ${COMMON_SRC})
target_include_directories(strings_tests PRIVATE ../src)

add_executable(series_tests test_series.c ${COMMON_SRC})
target_include_directories(series_tests PRIVATE ../src)

add_executable(series_twoterm_tests test_series_twoterm.c ${COMMON_SRC})
target_include_directories(series_twoterm_tests PRIVATE ../src)

add_executable(deriv_tests test_deriv.c ${COMMON_SRC})
target_include_directories(deriv_tests PRIVATE ../src)

add_executable(deriv_symbolic_order_tests test_deriv_symbolic_order.c ${COMMON_SRC})
target_include_directories(deriv_symbolic_order_tests PRIVATE ../src)

add_executable(limit_tests test_limit.c ${COMMON_SRC})
target_include_directories(limit_tests PRIVATE ../src)

add_executable(limit_assumptions_tests test_limit_assumptions.c ${COMMON_SRC})
target_include_directories(limit_assumptions_tests PRIVATE ../src)

add_executable(collect_corpus_tests test_collect_corpus.c ${COMMON_SRC})
target_include_directories(collect_corpus_tests PRIVATE ../src)

add_executable(power_corpus_tests test_power_corpus.c ${COMMON_SRC})
target_include_directories(power_corpus_tests PRIVATE ../src)

add_executable(radical_fusion_tests test_radical_fusion.c ${COMMON_SRC})
target_include_directories(radical_fusion_tests PRIVATE ../src)

add_executable(radical_polyops_tests test_radical_polyops.c ${COMMON_SRC})
target_include_directories(radical_polyops_tests PRIVATE ../src)

add_executable(radical_canonical_tests test_radical_canonical.c ${COMMON_SRC})
target_include_directories(radical_canonical_tests PRIVATE ../src)

add_executable(infinity_tests test_infinity.c ${COMMON_SRC})
target_include_directories(infinity_tests PRIVATE ../src)

add_executable(increment_tests test_increment.c ${COMMON_SRC})
target_include_directories(increment_tests PRIVATE ../src)

add_executable(numeric_tests test_numeric.c ${COMMON_SRC})
target_include_directories(numeric_tests PRIVATE ../src)

add_executable(simplify_tests test_simplify.c ${COMMON_SRC})
target_include_directories(simplify_tests PRIVATE ../src)

add_executable(trigrat_tests test_trigrat.c ${COMMON_SRC})
target_include_directories(trigrat_tests PRIVATE ../src)

add_executable(factorial_simplify_tests test_factorial_simplify.c ${COMMON_SRC})
target_include_directories(factorial_simplify_tests PRIVATE ../src)

add_executable(assuming_tests test_assuming.c ${COMMON_SRC})
target_include_directories(assuming_tests PRIVATE ../src)

add_executable(element_tests test_element.c ${COMMON_SRC})
target_include_directories(element_tests PRIVATE ../src)

add_executable(logexp_simplify_tests test_logexp_simplify.c ${COMMON_SRC})
target_include_directories(logexp_simplify_tests PRIVATE ../src)

add_executable(invtrig_simplify_tests test_invtrig_simplify.c ${COMMON_SRC})
target_include_directories(invtrig_simplify_tests PRIVATE ../src)

add_executable(radical_simplify_tests test_radical_simplify.c ${COMMON_SRC})
target_include_directories(radical_simplify_tests PRIVATE ../src)

add_executable(simp_denest_inv_sqrt_tests test_simp_denest_inv_sqrt.c ${COMMON_SRC})
target_include_directories(simp_denest_inv_sqrt_tests PRIVATE ../src)

add_executable(simp_denest_phase2_tests test_simp_denest_phase2.c ${COMMON_SRC})
target_include_directories(simp_denest_phase2_tests PRIVATE ../src)

add_executable(rationalize_tests test_rationalize.c ${COMMON_SRC})
target_include_directories(rationalize_tests PRIVATE ../src)

add_executable(inexact_dispatch_tests test_inexact_dispatch.c ${COMMON_SRC})
target_include_directories(inexact_dispatch_tests PRIVATE ../src)

if(USE_ECM)
    add_dependencies(expr_tests build_ecm)
    add_dependencies(parse_tests build_ecm)
    add_dependencies(eval_tests build_ecm)
    add_dependencies(match_tests build_ecm)
    add_dependencies(match_extensive_tests build_ecm)
    add_dependencies(symtab_tests build_ecm)
    add_dependencies(regression_tests build_ecm)
    add_dependencies(comparisons_tests build_ecm)
    add_dependencies(list_tests build_ecm)
    add_dependencies(iter_tests build_ecm)
    add_dependencies(cond_tests build_ecm)
    add_dependencies(patterns_tests build_ecm)
    add_dependencies(replace_tests build_ecm)
    add_dependencies(core_tests build_ecm)
    add_dependencies(trig_tests build_ecm)
    add_dependencies(trig_canon_tests build_ecm)
    add_dependencies(hyperbolic_tests build_ecm)
    add_dependencies(logexp_tests build_ecm)
    add_dependencies(piecewise_tests build_ecm)
    add_dependencies(purefunc_tests build_ecm)
    add_dependencies(modular_tests build_ecm)
    add_dependencies(funcprog_through_tests build_ecm)
    add_dependencies(mapat_tests build_ecm)
    add_dependencies(replaceat_tests build_ecm)
    add_dependencies(nest_tests build_ecm)
    add_dependencies(nestlist_tests build_ecm)
    add_dependencies(nestwhile_tests build_ecm)
    add_dependencies(nestwhilelist_tests build_ecm)
    add_dependencies(fixedpointlist_tests build_ecm)
    add_dependencies(fixedpoint_tests build_ecm)
    add_dependencies(fold_tests build_ecm)
    add_dependencies(foldlist_tests build_ecm)
    add_dependencies(split_tests build_ecm)
    add_dependencies(sort_tests build_ecm)
    add_dependencies(backtrack_tests build_ecm)
    add_dependencies(poly_tests build_ecm)
    add_dependencies(expand_tests build_ecm)
    add_dependencies(expandfrac_tests build_ecm)
    add_dependencies(stats_tests build_ecm)
    add_dependencies(datetime_tests build_ecm)
    add_dependencies(facpoly_tests build_ecm)
    add_dependencies(factor_terms_tests build_ecm)
    add_dependencies(rat_tests build_ecm)
    add_dependencies(linalg_tests build_ecm)
    add_dependencies(parfrac_tests build_ecm)
    add_dependencies(polymod_tests build_ecm)
    add_dependencies(bigint_tests build_ecm)
    add_dependencies(expr_sharing_tests build_ecm)
    add_dependencies(eval_timestamps_tests build_ecm)
    add_dependencies(eval_eager_exit_tests build_ecm)
    add_dependencies(rule_dispatch_tests build_ecm)
    add_dependencies(integrals_tests build_ecm)
    add_dependencies(print_tests build_ecm)
    add_dependencies(random_tests build_ecm)
    add_dependencies(evaluate_tests build_ecm)
    add_dependencies(releasehold_tests build_ecm)
    add_dependencies(unevaluated_tests build_ecm)
    add_dependencies(strings_tests build_ecm)
    add_dependencies(series_tests build_ecm)
    add_dependencies(series_twoterm_tests build_ecm)
    add_dependencies(deriv_tests build_ecm)
    add_dependencies(limit_tests build_ecm)
    add_dependencies(radical_fusion_tests build_ecm)
    add_dependencies(radical_polyops_tests build_ecm)
    add_dependencies(radical_canonical_tests build_ecm)
    add_dependencies(infinity_tests build_ecm)
    add_dependencies(increment_tests build_ecm)
    add_dependencies(numeric_tests build_ecm)
    add_dependencies(simplify_tests build_ecm)
    add_dependencies(factorial_simplify_tests build_ecm)
    add_dependencies(assuming_tests build_ecm)
    add_dependencies(element_tests build_ecm)
    add_dependencies(logexp_simplify_tests build_ecm)
    add_dependencies(invtrig_simplify_tests build_ecm)
    add_dependencies(radical_simplify_tests build_ecm)
    add_dependencies(rationalize_tests build_ecm)
    add_dependencies(inexact_dispatch_tests build_ecm)
endif()

add_executable(simp_tests test_simp.c ${COMMON_SRC})
target_include_directories(simp_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(simp_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(simp_tests gmp m)
endif()
add_test(NAME simp_tests COMMAND simp_tests)

add_executable(distribute_tests test_distribute.c ${COMMON_SRC})
target_include_directories(distribute_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(distribute_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
else()
    target_link_libraries(distribute_tests gmp m)
endif()
add_test(NAME distribute_tests COMMAND distribute_tests)

add_executable(trigexpand_tests test_trigexpand.c ${COMMON_SRC})
target_include_directories(trigexpand_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(trigexpand_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(trigexpand_tests build_ecm)
else()
    target_link_libraries(trigexpand_tests gmp m)
endif()
add_test(NAME trigexpand_tests COMMAND trigexpand_tests)

add_executable(trigfactor_tests test_trigfactor.c ${COMMON_SRC})
target_include_directories(trigfactor_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(trigfactor_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(trigfactor_tests build_ecm)
else()
    target_link_libraries(trigfactor_tests gmp m)
endif()
add_test(NAME trigfactor_tests COMMAND trigfactor_tests)

add_executable(trigreduce_tests test_trigreduce.c ${COMMON_SRC})
target_include_directories(trigreduce_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(trigreduce_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(trigreduce_tests build_ecm)
else()
    target_link_libraries(trigreduce_tests gmp m)
endif()
add_test(NAME trigreduce_tests COMMAND trigreduce_tests)

add_executable(recursion_limit_tests test_recursion_limit.c ${COMMON_SRC})
target_include_directories(recursion_limit_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(recursion_limit_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(recursion_limit_tests build_ecm)
else()
    target_link_libraries(recursion_limit_tests gmp m)
endif()
add_test(NAME recursion_limit_tests COMMAND recursion_limit_tests)

add_executable(return_tests test_return.c ${COMMON_SRC})
target_include_directories(return_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(return_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(return_tests build_ecm)
else()
    target_link_libraries(return_tests gmp m)
endif()
add_test(NAME return_tests COMMAND return_tests)

add_executable(repl_hooks_tests test_repl_hooks.c ${COMMON_SRC})
target_include_directories(repl_hooks_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(repl_hooks_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(repl_hooks_tests build_ecm)
else()
    target_link_libraries(repl_hooks_tests gmp m)
endif()
add_test(NAME repl_hooks_tests COMMAND repl_hooks_tests)

add_executable(extension_options_tests test_extension_options.c ${COMMON_SRC})
target_include_directories(extension_options_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(extension_options_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(extension_options_tests build_ecm)
else()
    target_link_libraries(extension_options_tests gmp m)
endif()
add_test(NAME extension_options_tests COMMAND extension_options_tests)

add_executable(intrat_tests test_intrat.c ${COMMON_SRC})
target_include_directories(intrat_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(intrat_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(intrat_tests build_ecm)
else()
    target_link_libraries(intrat_tests gmp m)
endif()
add_test(NAME intrat_tests COMMAND intrat_tests)

add_executable(intrischnorman_tests test_intrischnorman.c ${COMMON_SRC})
target_include_directories(intrischnorman_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(intrischnorman_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(intrischnorman_tests build_ecm)
else()
    target_link_libraries(intrischnorman_tests gmp m)
endif()
add_test(NAME intrischnorman_tests COMMAND intrischnorman_tests)

# Wall-clock benchmark harness for Integrate`RischNorman.  Not a test
# (no correctness asserts); built alongside the test suite so it is
# easy to rebuild and run for profiling.
add_executable(bench_intrischnorman bench_intrischnorman.c ${COMMON_SRC})
target_include_directories(bench_intrischnorman PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(bench_intrischnorman gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(bench_intrischnorman build_ecm)
else()
    target_link_libraries(bench_intrischnorman gmp m)
endif()

# Comprehensive corpus runner for IntegrateRationalTests.m and
# IntegrateRationalInlineCases.m.  Loads the .m file at runtime and
# runs every {integrand, var, ...} entry through the rational
# integrator with a differential check; uses fork-per-case to
# isolate hangs and crashes from later cases.
add_executable(intrat_corpus_tests test_intrat_corpus.c ${COMMON_SRC})
target_include_directories(intrat_corpus_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(intrat_corpus_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(intrat_corpus_tests build_ecm)
else()
    target_link_libraries(intrat_corpus_tests gmp m)
endif()
add_test(NAME intrat_corpus_tests
         COMMAND intrat_corpus_tests ../../IntegrateRationalTests.m)
add_test(NAME intrat_inline_corpus_tests
         COMMAND intrat_corpus_tests ../../IntegrateRationalInlineCases.m)

# Comprehensive corpus runner for the CRC integral table.  Loads
# CRCIntegralsCorpus.m (a flat List of {integrand, var} pairs, ~3 per
# rule in src/internal/CRCMathTablesIntegrals.m) and validates each via
# Integrate[..., Method -> "CRCTable"] + differential check.
add_executable(crc_corpus_tests test_crc_corpus.c ${COMMON_SRC})
target_include_directories(crc_corpus_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(crc_corpus_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(crc_corpus_tests build_ecm)
else()
    target_link_libraries(crc_corpus_tests gmp m)
endif()
add_test(NAME crc_corpus_tests
         COMMAND crc_corpus_tests ../../CRCIntegralsCorpus.m)

add_executable(readwrite_tests test_readwrite.c ${COMMON_SRC})
target_include_directories(readwrite_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(readwrite_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(readwrite_tests build_ecm)
else()
    target_link_libraries(readwrite_tests gmp m)
endif()
add_test(NAME readwrite_tests COMMAND readwrite_tests)

add_executable(files_tests test_files.c ${COMMON_SRC})
target_include_directories(files_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(files_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(files_tests build_ecm)
else()
    target_link_libraries(files_tests gmp m)
endif()
add_test(NAME files_tests COMMAND files_tests)

add_executable(time_constrained_tests test_time_constrained.c ${COMMON_SRC})
target_include_directories(time_constrained_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(time_constrained_tests gmp m ${CMAKE_SOURCE_DIR}/../src/external/ecm/.libs/libecm.a)
    add_dependencies(time_constrained_tests build_ecm)
else()
    target_link_libraries(time_constrained_tests gmp m)
endif()
add_test(NAME time_constrained_tests COMMAND time_constrained_tests)
