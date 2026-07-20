cmake_minimum_required(VERSION 3.10)
project(Mathilda_tests)

option(USE_ECM "Use ECM library for factorisation" ON)
option(USE_MPFR "Use MPFR for arbitrary-precision reals (N[expr, prec])" ON)
option(USE_LAPACK "Use BLAS/LAPACK for fast machine-precision linear algebra" ON)
option(USE_GRAPHICS "Use Raylib for the 2D graphics engine (Graphics[]/Show[]/Plot[] rendering)" ON)
option(USE_FLINT "Use FLINT (>= 3.0) for algebraic-extension polynomial arithmetic" ON)
option(USE_REGEX "Use PCRE2 for RegularExpression[] and regex string functions" ON)

if(USE_MPFR)
    add_compile_definitions(USE_MPFR)
    link_libraries(mpfr)
endif()

# BLAS/LAPACK autodetection — mirror the four-tier ladder used by the
# top-level makefile.  Auto-degrades to USE_LAPACK=OFF when not found
# so the test build never fails purely on a missing system library.
if(USE_LAPACK)
    if(APPLE)
        find_library(MATHILDA_ACCELERATE Accelerate)
        if(MATHILDA_ACCELERATE)
            add_compile_definitions(USE_LAPACK MATHILDA_USE_ACCELERATE)
            link_libraries(${MATHILDA_ACCELERATE})
        else()
            message(WARNING "Accelerate framework not found; building tests with USE_LAPACK=OFF")
            set(USE_LAPACK OFF)
        endif()
    else()
        find_package(LAPACK)
        find_path(MATHILDA_LAPACKE_INCLUDE_DIR lapacke.h)
        find_library(MATHILDA_LAPACKE_LIB lapacke)
        if(LAPACK_FOUND AND MATHILDA_LAPACKE_INCLUDE_DIR AND MATHILDA_LAPACKE_LIB)
            add_compile_definitions(USE_LAPACK)
            include_directories(${MATHILDA_LAPACKE_INCLUDE_DIR})
            link_libraries(${MATHILDA_LAPACKE_LIB} ${LAPACK_LIBRARIES})
        else()
            message(WARNING "LAPACK/LAPACKE not detected; building tests with USE_LAPACK=OFF")
            message(WARNING "  Ubuntu/Debian:  sudo apt install liblapacke-dev libopenblas-dev")
            message(WARNING "  Fedora/RHEL:    sudo dnf install lapack-devel lapacke-devel openblas-devel")
            message(WARNING "  Arch:           sudo pacman -S openblas lapacke")
            set(USE_LAPACK OFF)
        endif()
    endif()
endif()

# Raylib autodetection for the graphics engine.  Mirrors the makefile's
# pkg-config probe; degrades to USE_GRAPHICS=OFF (sampling.c and the
# symbolic Graphics[]/Plot[] construction path are still fully testable —
# only render.c/hershey_font.c, which need a live display, are excluded).
if(USE_GRAPHICS)
    # On macOS, Homebrew installs its .pc files outside the search path used by
    # a MacPorts/system pkg-config, so the probe below fails on a stock shell
    # even when raylib is installed. Prepend the Homebrew prefixes (Apple-
    # Silicon + Intel) so graphics is detected by default, matching the makefile.
    if(APPLE)
        set(ENV{PKG_CONFIG_PATH}
            "/opt/homebrew/lib/pkgconfig:/usr/local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(RAYLIB QUIET raylib)
    endif()
    if(RAYLIB_FOUND)
        add_compile_definitions(USE_GRAPHICS)
        include_directories(${RAYLIB_INCLUDE_DIRS})
        # On Linux, libraylib is commonly a static archive (Ubuntu's
        # libraylib-dev, or a from-source `make install`) whose transitive
        # windowing deps — X11, GL, pthread, dl, rt, ... — live in the .pc
        # file's `Libs.private:` line. pkg_check_modules exposes those only via
        # its _STATIC_ variables; linking the plain _LIBRARIES set drops them,
        # so every test executable fails with undefined `XInternAtom`, `XSync`,
        # `XFree`, etc. (issue #18). The _STATIC_ set is a superset that equals
        # the plain set on macOS (Libs.private is empty there — frameworks are
        # recorded in the dylib), so prefer it unconditionally, keeping
        # pkg-config's link order (raylib before its deps) that static archives
        # require. Fall back to the plain set if _STATIC_ ever comes back empty.
        if(RAYLIB_STATIC_LIBRARIES)
            link_directories(${RAYLIB_STATIC_LIBRARY_DIRS})
            link_libraries(${RAYLIB_STATIC_LIBRARIES})
        else()
            link_directories(${RAYLIB_LIBRARY_DIRS})
            link_libraries(${RAYLIB_LIBRARIES})
        endif()
    else()
        message(WARNING "Raylib not detected; building tests with USE_GRAPHICS=OFF")
        message(WARNING "  macOS (Homebrew): brew install raylib")
        message(WARNING "  Ubuntu/Debian:    sudo apt install libraylib-dev")
        set(USE_GRAPHICS OFF)
    endif()
endif()

# FLINT (>= 3.0) autodetection for the algebraic-extension engine. Mirrors the
# makefile's pkg-config probe with the same >= 3.0 version floor (ANTIC) and
# graceful degrade: when absent, flint_bridge.c compiles to stubs and the
# flint_bridge_tests binary skips itself at runtime.
if(USE_FLINT)
    if(APPLE)
        set(ENV{PKG_CONFIG_PATH}
            "/opt/homebrew/lib/pkgconfig:/usr/local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(FLINT QUIET flint>=3.0)
    endif()
    if(FLINT_FOUND)
        add_compile_definitions(USE_FLINT)
        include_directories(${FLINT_INCLUDE_DIRS})
        link_directories(${FLINT_LIBRARY_DIRS})
        link_libraries(${FLINT_LIBRARIES})
    else()
        message(WARNING "FLINT >= 3.0 not detected; building tests with USE_FLINT=OFF")
        message(WARNING "  macOS (Homebrew): brew install flint")
        message(WARNING "  Ubuntu/Debian:    sudo apt install libflint-dev   (needs >= 3.0)")
        set(USE_FLINT OFF)
    endif()
endif()

# PCRE2 autodetection — mirror the makefile's libpcre2-8 pkg-config probe with
# graceful degrade: when absent, the regex engine wrapper compiles to stubs and
# the regex string builtins warn and stay unevaluated at runtime.
if(USE_REGEX)
    if(APPLE)
        set(ENV{PKG_CONFIG_PATH}
            "/opt/homebrew/lib/pkgconfig:/usr/local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(PCRE2 QUIET libpcre2-8)
    endif()
    if(PCRE2_FOUND)
        add_compile_definitions(USE_REGEX PCRE2_CODE_UNIT_WIDTH=8)
        include_directories(${PCRE2_INCLUDE_DIRS})
        link_directories(${PCRE2_LIBRARY_DIRS})
        link_libraries(${PCRE2_LIBRARIES})
    else()
        message(WARNING "PCRE2 not detected; building tests with USE_REGEX=OFF")
        message(WARNING "  macOS (Homebrew): brew install pcre2")
        message(WARNING "  Ubuntu/Debian:    sudo apt install libpcre2-dev")
        set(USE_REGEX OFF)
    endif()
endif()

# FFTW autodetection — mirror the makefile's fftw3 pkg-config probe. When
# absent, the Fourier machine path compiles to a naive O(n^2) DFT fallback.
option(USE_FFTW "Use FFTW for the machine-precision Fourier transform" ON)
if(USE_FFTW)
    if(APPLE)
        set(ENV{PKG_CONFIG_PATH}
            "/opt/homebrew/lib/pkgconfig:/usr/local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(FFTW QUIET fftw3)
    endif()
    if(FFTW_FOUND)
        add_compile_definitions(USE_FFTW)
        include_directories(${FFTW_INCLUDE_DIRS})
        link_directories(${FFTW_LIBRARY_DIRS})
        link_libraries(${FFTW_LIBRARIES})
    else()
        message(WARNING "FFTW not detected; building tests with USE_FFTW=OFF")
        message(WARNING "  macOS (Homebrew): brew install fftw")
        message(WARNING "  Ubuntu/Debian:    sudo apt install libfftw3-dev")
        set(USE_FFTW OFF)
    endif()
endif()

include_directories(/usr/local/include)
include_directories(../src/list)
include_directories(../src/linalg)
include_directories(../src/numbertheory)
include_directories(../src/poly)
include_directories(../src/simp)
include_directories(../src/calculus)
include_directories(../src/sum)
include_directories(../src/special_functions)
include_directories(../src/numerical_calculus)
include_directories(../src/numerical_roots)
include_directories(../src/graphics)
include_directories(../src/graph)
include_directories(../src/strings)
include_directories(../src/strings/regex)
if(EXISTS "/opt/homebrew/include")
    include_directories(/opt/homebrew/include)
endif()
link_directories(/usr/local/lib)
if(EXISTS "/opt/homebrew/lib")
    link_directories(/opt/homebrew/lib)
endif()

# GMP-ECM is a system library (brew install gmp-ecm / apt install libecm-dev),
# no longer vendored. Detect the public header + shared lib and link -lecm.
if(USE_ECM)
    find_path(ECM_INCLUDE_DIR ecm.h PATHS /usr/local/include /opt/homebrew/include)
    find_library(ECM_LIBRARY NAMES ecm PATHS /usr/local/lib /opt/homebrew/lib)
    if(ECM_INCLUDE_DIR AND ECM_LIBRARY)
        include_directories(${ECM_INCLUDE_DIR})
        set(ECM_LIB ${ECM_LIBRARY})
    else()
        message(WARNING "GMP-ECM not found; building tests with NO_ECM. Install it: brew install gmp-ecm / sudo apt install libecm-dev")
        set(USE_ECM OFF)
    endif()
endif()

# A no-op target so the many `add_dependencies(<test> build_ecm)` calls below
# keep resolving whether or not ECM is present (nothing to build any more).
add_custom_target(build_ecm)

if(USE_ECM)
    link_libraries(gmp m ${ECM_LIB})
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
    ../src/int.c
    ../src/real.c
    ../src/comparisons.c
    ../src/symtab.c
    ../src/match.c
    ../src/part.c
    ../src/vectors.c
    ../src/assoc.c
    ../src/ndarray.c
    ../src/ndkernels.c
    ../src/print.c
    ../src/plus.c
    ../src/times.c
    ../src/arithmetic.c
    ../src/numbertheory/numbertheory.c
    ../src/numbertheory/nt_util.c
    ../src/numbertheory/nt_gaussian.c
    ../src/numbertheory/gcd.c
    ../src/numbertheory/lcm.c
    ../src/numbertheory/extendedgcd.c
    ../src/numbertheory/powermod.c
    ../src/numbertheory/factorial.c
    ../src/numbertheory/factorialpower.c
    ../src/numbertheory/factorial2.c
    ../src/numbertheory/binomial.c
    ../src/numbertheory/primitiveroot.c
    ../src/numbertheory/multiplicativeorder.c
    ../src/numbertheory/divisible.c
    ../src/numbertheory/coprimeq.c
    ../src/numbertheory/divisors.c
    ../src/numbertheory/divisorsigma.c
    ../src/numbertheory/moebiusmu.c
    ../src/numbertheory/liouvillelambda.c
    ../src/numbertheory/primeomega.c
    ../src/numbertheory/primenu.c
    ../src/numbertheory/jacobisymbol.c
    ../src/numbertheory/prime.c
    ../src/numbertheory/primecount.c
    ../src/partitions.c
    ../src/special_functions/gamma.c
    ../src/special_functions/beta.c
    ../src/special_functions/erf.c
    ../src/special_functions/erfc.c
    ../src/special_functions/erfi.c
    ../src/special_functions/expintegralei.c
    ../src/special_functions/logintegral.c
    ../src/special_functions/sinintegral.c
    ../src/special_functions/cosintegral.c
    ../src/special_functions/sinhintegral.c
    ../src/special_functions/coshintegral.c
    ../src/special_functions/fresnel.c
    ../src/special_functions/sinc.c
    ../src/special_functions/inverf.c
    ../src/special_functions/inverfc.c
    ../src/special_functions/loggamma.c
    ../src/special_functions/polygamma.c
    ../src/special_functions/harmonicnumber.c
    ../src/special_functions/pochhammer.c
    ../src/special_functions/eulergamma.c
    ../src/special_functions/zeta.c
    ../src/special_functions/hurwitzzeta.c
    ../src/special_functions/stieltjesgamma.c
    ../src/special_functions/bernoullib.c
    ../src/special_functions/eulere.c
    ../src/special_functions/polylog.c
    ../src/special_functions/lerchphi.c
    ../src/special_functions/productlog.c
    ../src/special_functions/legendre.c
    ../src/special_functions/airyai.c
    ../src/special_functions/airybi.c
    ../src/special_functions/bessel.c
    ../src/facint.c
    ../src/special_functions/fibonacci.c
    ../src/special_functions/lucas.c
    ../src/special_functions/hyperfactorial.c
    ../src/special_functions/barnesg.c
    ../src/special_functions/qpochhammer.c
    ../src/special_functions/besseljzero.c
    ../src/power.c
    ../src/funcprog.c
    ../src/boolean.c
    ../src/names.c
    ../src/list/table.c
    ../src/list/range.c
    ../src/list/array.c
    ../src/list/constant_array.c
    ../src/list/array_flatten.c
    ../src/list/take_drop.c
    ../src/list/flatten.c
    ../src/list/partition.c
    ../src/list/rotate.c
    ../src/list/reverse.c
    ../src/list/transpose.c
    ../src/list/setops.c
    ../src/list/split.c
    ../src/list/rescale.c
    ../src/list/pad.c
    ../src/list/total.c
    ../src/list/accumulate.c
    ../src/list/differences.c
    ../src/list/ratios.c
    ../src/list/listpredicates.c
    ../src/list/matrixq.c
    ../src/list/minmax.c
    ../src/list/join.c
    ../src/list/list_common.c
    ../src/list/list_init.c
    ../src/patterns.c
    ../src/patterns.c
    ../src/replace.c
    ${CMAKE_SOURCE_DIR}/../src/simp/trigsimp.c
    ${CMAKE_SOURCE_DIR}/../src/simp/trigrat.c
    ${CMAKE_SOURCE_DIR}/../src/simp/radrat.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_log.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_complexity.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_assume.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_util.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_assume_rewrite.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_trig_roundtrip.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_denest.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_cuberoot.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_rationalize.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_algebraic.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_canon.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_trig_pi.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_power.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_tan_add.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_logexp_abs.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_search.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_factorial.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_bottomup.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_builtins.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_trigexp_zero.c
    ${CMAKE_SOURCE_DIR}/../src/simp/simp_transform.c
    ../src/complex.c
    ../src/complex_expand.c
    ../src/fourier.c
    ../src/convolutions.c
    ../src/trig.c
    ../src/trig_canon.c
    ../src/hyperbolic.c
    ../src/logexp.c
    ../src/special_functions/hypergeopfq.c
    ../src/piecewise.c
    ../src/interp.c
    ../src/interp_mpfr.c
    ../src/purefunc.c
    ../src/modular.c
    ../src/sort.c
    ../src/stats.c
    ../src/fit.c
    ../src/info.c
    ../src/expand.c
    ../src/expand_power.c
    ../src/poly/poly.c
    ../src/poly/minpoly.c
    ../src/poly/subresultants.c
    ../src/poly/subresultantpoly.c
    ../src/rat.c
    ../src/rootreduce.c
    ../src/poly/facpoly.c
    ../src/parfrac.c
    ../src/contfrac.c
    ../src/datetime.c
    ../src/linalg/linalg.c
    ../src/linalg/util.c
    ../src/linalg/dot.c
    ../src/linalg/det.c
    ../src/linalg/flint_mat_bridge.c
    ../src/linalg/cross.c
    ../src/linalg/norm.c
    ../src/linalg/normalize.c
    ../src/linalg/tr.c
    ../src/linalg/construct.c
    ../src/linalg/hankelmat.c
    ../src/linalg/toeplitzmat.c
    ../src/linalg/vandermondemat.c
    ../src/linalg/hilbertmat.c
    ../src/linalg/matpow.c
    ../src/linalg/linsolve.c
    ../src/linalg/inv.c
    ../src/linalg/lstsq.c
    ../src/linalg/nullspace.c
    ../src/linalg/matrank.c
    ../src/linalg/latticereduce.c
    ../src/linalg/qrdecomp.c
    ../src/linalg/qrdecomp_machine.c
    ../src/linalg/qrdecomp_mpfr.c
    ../src/linalg/ludecomp.c
    ../src/linalg/ludecomp_machine.c
    ../src/linalg/ludecomp_mpfr.c
    ../src/linalg/svdecomp.c
    ../src/linalg/svdecomp_machine.c
    ../src/linalg/svdecomp_mpfr.c
    ../src/linalg/lapack.c
    ../src/linalg/numarray.c
    ../src/linalg/ndlinalg.c
    ../src/linalg/blas_bridge.c
    ../src/linalg/lapack_bridge.c
    ../src/linalg/posdef_q.c
    ../src/linalg/negdef_q.c
    ../src/linalg/eigen.c
    ../src/linalg/eigen_common.c
    ../src/linalg/eigen_direct.c
    ../src/linalg/eigen_arnoldi.c
    ../src/linalg/eigen_banded.c
    ../src/linalg/eigen_feast.c
    ../src/readwrite.c
    ../src/loadmodule.c
    ../src/files.c
    ../src/iter.c
    ../src/cond.c
    ../src/default_helper.c
    ../src/random.c
    ../src/strings/stringlength.c
    ../src/strings/characters.c
    ../src/strings/stringjoin.c
    ../src/strings/stringpart.c
    ../src/strings/stringtake.c
    ../src/strings/stringdrop.c
    ../src/strings/stringreverse.c
    ../src/strings/stringinsert.c
    ../src/strings/stringreplacepart.c
    ../src/strings/stringsplit.c
    ../src/strings/stringtrim.c
    ../src/strings/stringextract.c
    ../src/strings/stringpartition.c
    ../src/strings/stringrepeat.c
    ../src/strings/strings_init.c
    ../src/strings/regex/regex_engine.c
    ../src/strings/regex/regex_common.c
    ../src/strings/regex/string_pattern.c
    ../src/strings/regex/regularexpression.c
    ../src/strings/regex/stringmatchq.c
    ../src/strings/regex/stringcases.c
    ../src/strings/regex/stringreplace.c
    ../src/strings/regex/regex_init.c
    ../src/calculus/series.c
    ../src/calculus/deriv.c
    ../src/calculus/limit.c
    ../src/calculus/residue.c
    ../src/context.c
    ../src/numeric.c
    ../src/precision.c
    ../src/rationalize.c
    ../src/common.c
    ../src/poly/zupoly.c
    ../src/poly/flint_bridge.c
    ../src/poly/flint_qqbar.c
    ../src/poly/bpoly.c
    ../src/poly/mpoly.c
    ../src/poly/mvfactor.c
    ../src/poly/mvfactor3.c
    ../src/poly/qa.c
    ../src/poly/qaupoly.c
    ../src/poly/qafactor.c
    ../src/poly/rootofunity.c
    ../src/poly/squarefreeq.c
    ../src/poly/exponent.c
    ../src/poly/facpoly_list.c
    ../src/poly/irrpolyq.c
    ../src/poly/groebner.c
    ../src/poly/groebnerbasis.c
    ../src/poly/groebnerwalk.c
    ../src/poly/eliminate.c
    ../src/sym_intern.c
    ../src/sym_names.c
    ../src/version.c
    ../src/repl_hooks.c
    ../src/options.c
    ../src/options_builtin.c
    ../src/calculus/integrate.c
    ../src/calculus/integrate_interp.c
    ../src/calculus/integrate_unknown.c
    ../src/calculus/integrate_derivdivides.c
    ../src/calculus/integrate_linrad.c
    ../src/calculus/integrate_quadrad.c
    ../src/calculus/integrate_linratiorad.c
    ../src/calculus/integrate_chebychev.c
    ../src/calculus/integrate_goursat.c
    ../src/calculus/integrate_fresnel.c
    ../src/calculus/integrate_jeffrey.c
    ../src/calculus/integrate_newton_leibniz.c
    ../src/calculus/integrate_symmetry.c
    ../src/calculus/integrate_beta.c
    ../src/calculus/integrate_line.c
    ../src/calculus/integrate_residue.c
    ../src/calculus/integrate_diffunderint.c
    ../src/calculus/integrate_ramanujan.c
    ../src/calculus/intrat.c
    ../src/calculus/intrischnorman.c
    ../src/calculus/integrate_risch_transcendental.c
    ../src/calculus/risch_util.c
    ../src/calculus/integrate_risch_rde.c
    ../src/calculus/risch_tower.c
    ../src/calculus/risch_field_integrate.c
    ../src/calculus/risch_singleext.c
    ../src/calculus/risch_trig_frontend.c
    ../src/calculus/risch_special.c
    ../src/calculus/cherry_driver.c
    ../src/calculus/cherry_ei.c
    ../src/calculus/cherry_li.c
    ../src/calculus/cherry_dilog.c
    ../src/calculus/cherry_sigma_decomp.c
    ../src/calculus/risch_field.c
    ../src/calculus/risch_canonical.c
    ../src/calculus/risch_structure.c
    ../src/calculus/risch_hypertangent.c
    ../src/calculus/risch_coupled.c
    ../src/calculus/risch_hermite.c
    ../src/calculus/intsimp.c
    ../src/sum/sum.c
    ../src/sum/sum_util.c
    ../src/sum/sum_polynomial.c
    ../src/sum/sum_geometric.c
    ../src/sum/sum_gosper.c
    ../src/sum/sum_hypergeometric.c
    ../src/sum/sum_rational.c
    ../src/sum/sum_euler.c
    ../src/sum/sum_euler_nonlinear.c
    ../src/sum/sum_dirichlet.c
    ../src/sum/sum_zetaseries.c
    ../src/sum/sum_alternating.c
    ../src/sum/sum_trigonometric.c
    ../src/sum/sum_logzeta.c
    ../src/sum/sum_logrational.c
    ../src/product/product.c
    ../src/product/product_util.c
    ../src/product/product_telescoping.c
    ../src/product/product_rational.c
    ../src/product/product_geometric.c
    ../src/product/product_qproduct.c
    ../src/product/product_special.c
    ../src/product/product_infinite.c
    ../src/product/product_cantor.c
    ../src/product/product_viete.c
    ../src/product/product_eulerprime.c
    ../src/product/product_rational_infinite.c
    ../src/product/product_besselzero.c
    ../src/product/product_logsum.c
    ../src/root.c
    ../src/root_numeric.c
    ../src/radicals.c
    ../src/poly/poly_eval_mpfr.c
    ../src/solve.c
    ../src/solvealways.c
    ../src/numerical_roots/findroot.c
    ../src/findmin.c
    ../src/poly/solvepoly.c
    ../src/solvelinsys.c
    ../src/solvenlsys.c
    ../src/solverad.c
    ../src/solveinv.c
    ../src/solvetrig.c
    ../src/zero_test.c
    ../src/numeric_complex.c
    ../src/flint_num_bridge.c
    ../src/numerical_calculus/quadrature.c
    ../src/numerical_calculus/seqaccel.c
    ../src/numerical_calculus/dequad.c
    ../src/numerical_calculus/nresidue.c
    ../src/numerical_calculus/nderiv.c
    ../src/numerical_calculus/nseries.c
    ../src/numerical_calculus/nlimit.c
    ../src/numerical_calculus/nsum.c
    ../src/numerical_calculus/nprod.c
    ../src/numerical_calculus/gkadapt.c
    ../src/numerical_calculus/denint.c
    ../src/numerical_calculus/oscint.c
    ../src/numerical_calculus/oscde.c
    ../src/numerical_calculus/mcint.c
    ../src/numerical_calculus/cubature.c
    ../src/numerical_calculus/ncrule.c
    ../src/numerical_calculus/levincoll.c
    ../src/numerical_calculus/nint.c
    ../src/numerical_roots/nroots.c
    ../src/numerical_roots/nroots_aberth.c
    ../src/numerical_roots/nroots_jt.c
    ../src/numerical_roots/nsolve.c
    ../src/numerical_roots/nsolve_system.c
    ../src/graphics/graphics_init.c
    ../src/graphics/graphics_json.c
    ../src/graphics/sampling.c
    ../src/graphics/plot_common.c
    ../src/graphics/plot.c
    ../src/graphics/plot3d.c
    ../src/graphics/parametricplot.c
    ../src/graphics/parametricplot3d.c
    ../src/graphics/contourplot.c
    ../src/graphics/densityplot.c
    ../src/graphics/complexplot.c
    ../src/graphics/streamplot.c
    ../src/graphics/vectorplot.c
    ../src/graphics/polarplot.c
    ../src/graphics/barchart.c
    ../src/graphics/animate.c
    ../src/graphics/listplot.c
    ../src/graphics/show.c
    ../src/graph/graph.c
    ../src/graph/graph_util.c
    ../src/graph/construct.c
    ../src/graph/graphq.c
    ../src/graph/vertexlist.c
    ../src/graph/edgelist.c
    ../src/graph/counts.c
    ../src/graph/adjlist.c
    ../src/graph/degree.c
    ../src/graph/directedq.c
    ../src/graph/adjmat.c
    ../src/graph/incmat.c
    ../src/graph/adjgraph.c
    ../src/graph/generators.c
    ../src/graph/shortestpath.c
    ../src/graph/components.c
    ../src/graph/spanningtree.c
    ../src/graph/connectivity.c
    ../src/graph/graphplot.c
)

# render.c / render3d.c / hershey_font.c need a live Raylib + display; only
# compiled into COMMON_SRC when USE_GRAPHICS actually resolved to ON above.
if(USE_GRAPHICS)
    list(APPEND COMMON_SRC
        ../src/graphics/render.c
        ../src/graphics/render3d.c
        ../src/graphics/hershey_font.c
    )
endif()

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

add_executable(graph_tests test_graph.c ${COMMON_SRC})
target_link_libraries(graph_tests m)
target_include_directories(graph_tests PRIVATE ../src ../src/graph)
add_test(NAME graph_tests COMMAND graph_tests)

add_executable(zero_test_tests test_zero_test.c ${COMMON_SRC})
target_link_libraries(zero_test_tests m)
target_include_directories(zero_test_tests PRIVATE ../src)
add_test(NAME zero_test_tests COMMAND zero_test_tests)

add_executable(trigexp_zero_tests test_trigexp_zero.c ${COMMON_SRC})
target_link_libraries(trigexp_zero_tests m)
target_include_directories(trigexp_zero_tests PRIVATE ../src)
add_test(NAME trigexp_zero_tests COMMAND trigexp_zero_tests)

add_executable(rootreduce_tests test_rootreduce.c ${COMMON_SRC})
target_include_directories(rootreduce_tests PRIVATE ../src)
add_test(NAME rootreduce_tests COMMAND rootreduce_tests)

add_executable(flint_bridge_tests test_flint_bridge.c ${COMMON_SRC})
target_include_directories(flint_bridge_tests PRIVATE ../src)
add_test(NAME flint_bridge_tests COMMAND flint_bridge_tests)

add_executable(list_tests test_list.c ${COMMON_SRC})
target_include_directories(list_tests PRIVATE ../src)

add_executable(intersection_tests test_intersection.c ${COMMON_SRC})
target_include_directories(intersection_tests PRIVATE ../src)
add_test(NAME intersection_tests COMMAND intersection_tests)

add_executable(fourier_tests test_fourier.c ${COMMON_SRC})
target_include_directories(fourier_tests PRIVATE ../src)
add_test(NAME fourier_tests COMMAND fourier_tests)

add_executable(convolutions_tests test_convolutions.c ${COMMON_SRC})
target_include_directories(convolutions_tests PRIVATE ../src)
add_test(NAME convolutions_tests COMMAND convolutions_tests)

add_executable(correlations_tests test_correlations.c ${COMMON_SRC})
target_include_directories(correlations_tests PRIVATE ../src)
add_test(NAME correlations_tests COMMAND correlations_tests)

add_executable(complement_tests test_complement.c ${COMMON_SRC})
target_include_directories(complement_tests PRIVATE ../src)
add_test(NAME complement_tests COMMAND complement_tests)

add_executable(constant_array_tests test_constant_array.c ${COMMON_SRC})
target_include_directories(constant_array_tests PRIVATE ../src)
add_test(NAME constant_array_tests COMMAND constant_array_tests)

add_executable(unit_vector_tests test_unit_vector.c ${COMMON_SRC})
target_include_directories(unit_vector_tests PRIVATE ../src)
add_test(NAME unit_vector_tests COMMAND unit_vector_tests)

add_executable(array_flatten_tests test_array_flatten.c ${COMMON_SRC})
target_include_directories(array_flatten_tests PRIVATE ../src)
add_test(NAME array_flatten_tests COMMAND array_flatten_tests)

add_executable(rescale_tests test_rescale.c ${COMMON_SRC})
target_include_directories(rescale_tests PRIVATE ../src)
add_test(NAME rescale_tests COMMAND rescale_tests)

add_executable(padright_tests test_padright.c ${COMMON_SRC})
target_include_directories(padright_tests PRIVATE ../src)
add_test(NAME padright_tests COMMAND padright_tests)

add_executable(padleft_tests test_padleft.c ${COMMON_SRC})
target_include_directories(padleft_tests PRIVATE ../src)
add_test(NAME padleft_tests COMMAND padleft_tests)

add_executable(divisible_tests test_divisible.c ${COMMON_SRC})
target_include_directories(divisible_tests PRIVATE ../src)
add_test(NAME divisible_tests COMMAND divisible_tests)

add_executable(coprimeq_tests test_coprimeq.c ${COMMON_SRC})
target_include_directories(coprimeq_tests PRIVATE ../src)
add_test(NAME coprimeq_tests COMMAND coprimeq_tests)

add_executable(exponent_tests test_exponent.c ${COMMON_SRC})
target_include_directories(exponent_tests PRIVATE ../src)
add_test(NAME exponent_tests COMMAND exponent_tests)

add_executable(factorlist_tests test_factorlist.c ${COMMON_SRC})
target_include_directories(factorlist_tests PRIVATE ../src)
add_test(NAME factorlist_tests COMMAND factorlist_tests)

add_executable(divisors_tests test_divisors.c ${COMMON_SRC})
target_include_directories(divisors_tests PRIVATE ../src)
add_test(NAME divisors_tests COMMAND divisors_tests)

add_executable(divisorsigma_tests test_divisorsigma.c ${COMMON_SRC})
target_include_directories(divisorsigma_tests PRIVATE ../src)
add_test(NAME divisorsigma_tests COMMAND divisorsigma_tests)

add_executable(moebiusmu_tests test_moebiusmu.c ${COMMON_SRC})
target_include_directories(moebiusmu_tests PRIVATE ../src)
add_test(NAME moebiusmu_tests COMMAND moebiusmu_tests)

add_executable(liouvillelambda_tests test_liouvillelambda.c ${COMMON_SRC})
target_include_directories(liouvillelambda_tests PRIVATE ../src)
add_test(NAME liouvillelambda_tests COMMAND liouvillelambda_tests)

add_executable(primeomega_tests test_primeomega.c ${COMMON_SRC})
target_include_directories(primeomega_tests PRIVATE ../src)
add_test(NAME primeomega_tests COMMAND primeomega_tests)

add_executable(primenu_tests test_primenu.c ${COMMON_SRC})
target_include_directories(primenu_tests PRIVATE ../src)
add_test(NAME primenu_tests COMMAND primenu_tests)

add_executable(prime_tests test_prime.c ${COMMON_SRC})
target_include_directories(prime_tests PRIVATE ../src)
add_test(NAME prime_tests COMMAND prime_tests)

add_executable(iter_tests test_iter.c ${COMMON_SRC})
target_link_libraries(iter_tests m)
target_include_directories(iter_tests PRIVATE ../src)

add_executable(sum_tests test_sum.c ${COMMON_SRC})
target_link_libraries(sum_tests m)
target_include_directories(sum_tests PRIVATE ../src)

add_executable(sum_rational_tests test_sum_rational.c ${COMMON_SRC})
target_link_libraries(sum_rational_tests m)
target_include_directories(sum_rational_tests PRIVATE ../src)

add_executable(sum_product_families_tests test_sum_product_families.c ${COMMON_SRC})
target_link_libraries(sum_product_families_tests m)
target_include_directories(sum_product_families_tests PRIVATE ../src)

add_executable(sum_euler_tests test_sum_euler.c ${COMMON_SRC})
target_link_libraries(sum_euler_tests m)
target_include_directories(sum_euler_tests PRIVATE ../src)
add_test(NAME sum_euler_tests COMMAND sum_euler_tests)

add_executable(sum_alternating_tests test_sum_alternating.c ${COMMON_SRC})
target_link_libraries(sum_alternating_tests m)
target_include_directories(sum_alternating_tests PRIVATE ../src)
add_test(NAME sum_alternating_tests COMMAND sum_alternating_tests)

add_executable(sum_trigonometric_tests test_sum_trigonometric.c ${COMMON_SRC})
target_link_libraries(sum_trigonometric_tests m)
target_include_directories(sum_trigonometric_tests PRIVATE ../src)
add_test(NAME sum_trigonometric_tests COMMAND sum_trigonometric_tests)

add_executable(product_tests test_product.c ${COMMON_SRC})
target_link_libraries(product_tests m)
target_include_directories(product_tests PRIVATE ../src ../src/product)

add_executable(product_stress_tests test_product_stress.c ${COMMON_SRC})
target_link_libraries(product_stress_tests m)
target_include_directories(product_stress_tests PRIVATE ../src ../src/product)

add_executable(product_special_tests test_product_special.c ${COMMON_SRC})
target_link_libraries(product_special_tests m)
target_include_directories(product_special_tests PRIVATE ../src ../src/special_functions)

add_executable(product_infinite_families_tests test_product_infinite_families.c ${COMMON_SRC})
target_link_libraries(product_infinite_families_tests m)
target_include_directories(product_infinite_families_tests PRIVATE ../src ../src/product)

add_executable(sum_log_tests test_sum_log.c ${COMMON_SRC})
target_link_libraries(sum_log_tests m)
target_include_directories(sum_log_tests PRIVATE ../src ../src/sum)

add_executable(gamma_transforms_tests test_gamma_transforms.c ${COMMON_SRC})
target_link_libraries(gamma_transforms_tests m)
target_include_directories(gamma_transforms_tests PRIVATE ../src)

add_executable(hypergeopfq_tests test_hypergeopfq.c ${COMMON_SRC})
target_link_libraries(hypergeopfq_tests m)
target_include_directories(hypergeopfq_tests PRIVATE ../src)

add_executable(cond_tests test_cond.c ${COMMON_SRC})
target_link_libraries(cond_tests m)
target_include_directories(cond_tests PRIVATE ../src)

add_executable(catch_throw_tests test_catch_throw.c ${COMMON_SRC})
target_link_libraries(catch_throw_tests m)
target_include_directories(catch_throw_tests PRIVATE ../src)

add_executable(goto_label_tests test_goto_label.c ${COMMON_SRC})
target_link_libraries(goto_label_tests m)
target_include_directories(goto_label_tests PRIVATE ../src)
add_test(NAME goto_label_tests COMMAND goto_label_tests)

add_executable(boolean_tests test_boolean.c ${COMMON_SRC})
target_link_libraries(boolean_tests m)
target_include_directories(boolean_tests PRIVATE ../src)

add_executable(names_tests test_names.c ${COMMON_SRC})
target_link_libraries(names_tests m)
target_include_directories(names_tests PRIVATE ../src)

add_executable(association_tests test_association.c ${COMMON_SRC})
target_link_libraries(association_tests m)
target_include_directories(association_tests PRIVATE ../src)
add_test(NAME association_tests COMMAND association_tests)

# Association performance regression gate: fails if a hash-backed op scales
# worse than O(n) (doubling-ratio check; machine-independent).
add_executable(bench_assoc bench_assoc.c ${COMMON_SRC})
target_link_libraries(bench_assoc m)
target_include_directories(bench_assoc PRIVATE ../src)
add_test(NAME bench_assoc COMMAND bench_assoc)

# Evaluator baseline benchmark: Phase 0 of EVAL_SYMTAB_IMPROVEMENTS.md. Measures
# per-node dispatch cost so the symbol-object rework can be validated against a
# committed baseline. Machine-normalized "much slower" gate like bench_assoc.
add_executable(bench_eval bench_eval.c ${COMMON_SRC})
target_link_libraries(bench_eval m)
target_include_directories(bench_eval PRIVATE ../src)
add_test(NAME bench_eval COMMAND bench_eval)

add_executable(bench_ndarray_linalg bench_ndarray_linalg.c ${COMMON_SRC})
target_link_libraries(bench_ndarray_linalg m)
target_include_directories(bench_ndarray_linalg PRIVATE ../src)
add_test(NAME bench_ndarray_linalg COMMAND bench_ndarray_linalg)

add_executable(ndarray_tests test_ndarray.c ${COMMON_SRC})
target_link_libraries(ndarray_tests m)
target_include_directories(ndarray_tests PRIVATE ../src)
add_test(NAME ndarray_tests COMMAND ndarray_tests)

add_executable(ndarray_functions_tests test_ndarray_functions.c ${COMMON_SRC})
target_link_libraries(ndarray_functions_tests m)
target_include_directories(ndarray_functions_tests PRIVATE ../src)
add_test(NAME ndarray_functions_tests COMMAND ndarray_functions_tests)

add_executable(ndarray_linalg_tests test_ndarray_linalg.c ${COMMON_SRC})
target_link_libraries(ndarray_linalg_tests m)
target_include_directories(ndarray_linalg_tests PRIVATE ../src)
add_test(NAME ndarray_linalg_tests COMMAND ndarray_linalg_tests)

add_executable(patterns_tests test_patterns.c ${COMMON_SRC})
target_link_libraries(patterns_tests m)
target_include_directories(patterns_tests PRIVATE ../src)

add_executable(replace_tests test_replace.c ${COMMON_SRC})
target_include_directories(replace_tests PRIVATE ../src)

add_executable(sequence_tests test_sequence.c ${COMMON_SRC})
target_include_directories(sequence_tests PRIVATE ../src)
add_test(NAME sequence_tests COMMAND sequence_tests)

add_executable(options_tests test_options.c ${COMMON_SRC})
target_include_directories(options_tests PRIVATE ../src)

add_executable(core_tests test_core.c ${COMMON_SRC})
target_include_directories(core_tests PRIVATE ../src)

add_executable(chop_tests test_chop.c ${COMMON_SRC})
target_include_directories(chop_tests PRIVATE ../src)

add_executable(nresidue_tests test_nresidue.c ${COMMON_SRC})
target_include_directories(nresidue_tests PRIVATE ../src ../src/numerical_calculus)

add_executable(nderiv_tests test_nderiv.c ${COMMON_SRC})
target_include_directories(nderiv_tests PRIVATE ../src ../src/numerical_calculus)

add_executable(nseries_tests test_nseries.c ${COMMON_SRC})
target_include_directories(nseries_tests PRIVATE ../src ../src/numerical_calculus)

add_executable(nlimit_tests test_nlimit.c ${COMMON_SRC})
target_include_directories(nlimit_tests PRIVATE ../src ../src/numerical_calculus)

add_executable(nsum_tests test_nsum.c ${COMMON_SRC})
target_include_directories(nsum_tests PRIVATE ../src ../src/numerical_calculus)

add_executable(nprod_tests test_nprod.c ${COMMON_SRC})
target_include_directories(nprod_tests PRIVATE ../src ../src/numerical_calculus)

add_executable(nint_tests test_nint.c ${COMMON_SRC})
target_include_directories(nint_tests PRIVATE ../src ../src/numerical_calculus)

add_executable(levincoll_tests test_levincoll.c ${COMMON_SRC})
target_include_directories(levincoll_tests PRIVATE ../src ../src/numerical_calculus)

add_executable(nroots_tests test_nroots.c ${COMMON_SRC})
target_include_directories(nroots_tests PRIVATE ../src ../src/numerical_roots)

add_executable(nsolve_tests test_nsolve.c ${COMMON_SRC})
target_include_directories(nsolve_tests PRIVATE ../src ../src/numerical_roots)

add_executable(nsolve_stress_tests test_nsolve_stress.c ${COMMON_SRC})
target_include_directories(nsolve_stress_tests PRIVATE ../src ../src/numerical_roots)

add_executable(clearall_remove_protect_tests test_clearall_remove_protect.c ${COMMON_SRC})
target_include_directories(clearall_remove_protect_tests PRIVATE ../src)

add_executable(unset_tests test_unset.c ${COMMON_SRC})
target_include_directories(unset_tests PRIVATE ../src)

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

add_executable(primitive_root_tests test_primitive_root.c ${COMMON_SRC})
target_include_directories(primitive_root_tests PRIVATE ../src)

add_executable(multiplicative_order_tests test_multiplicative_order.c ${COMMON_SRC})
target_include_directories(multiplicative_order_tests PRIVATE ../src)

add_executable(findroot_tests test_findroot.c ${COMMON_SRC})
target_include_directories(findroot_tests PRIVATE ../src ../src/numerical_roots)

add_executable(root_numeric_tests test_root_numeric.c ${COMMON_SRC})
target_include_directories(root_numeric_tests PRIVATE ../src)

add_executable(radicals_tests test_radicals.c ${COMMON_SRC})
target_include_directories(radicals_tests PRIVATE ../src)

add_executable(groebner_tests test_groebner.c ${COMMON_SRC})
target_include_directories(groebner_tests PRIVATE ../src)

add_executable(eliminate_tests test_eliminate.c ${COMMON_SRC})
target_include_directories(eliminate_tests PRIVATE ../src)

add_executable(findmin_tests test_findmin.c ${COMMON_SRC})
target_include_directories(findmin_tests PRIVATE ../src)

add_executable(fit_tests test_fit.c ${COMMON_SRC})
target_include_directories(fit_tests PRIVATE ../src)
add_test(NAME fit_tests COMMAND fit_tests)

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

add_executable(differences_tests test_differences.c ${COMMON_SRC})
target_include_directories(differences_tests PRIVATE ../src)

add_executable(ratios_tests test_ratios.c ${COMMON_SRC})
target_include_directories(ratios_tests PRIVATE ../src)

add_executable(sort_tests test_sort.c ${COMMON_SRC})
target_include_directories(sort_tests PRIVATE ../src)

add_executable(backtrack_tests test_backtrack.c ${COMMON_SRC})
target_include_directories(backtrack_tests PRIVATE ../src)

add_executable(poly_tests test_poly.c ${COMMON_SRC})
target_include_directories(poly_tests PRIVATE ../src)

add_executable(risch_field_tests test_risch_field.c ${COMMON_SRC})
target_include_directories(risch_field_tests PRIVATE ../src ../src/calculus)
add_test(NAME risch_field_tests COMMAND risch_field_tests)

add_executable(risch_canonical_tests test_risch_canonical.c ${COMMON_SRC})
target_include_directories(risch_canonical_tests PRIVATE ../src ../src/calculus)
add_test(NAME risch_canonical_tests COMMAND risch_canonical_tests)

add_executable(risch_canonical_wiring_tests test_risch_canonical_wiring.c ${COMMON_SRC})
target_include_directories(risch_canonical_wiring_tests PRIVATE ../src ../src/calculus)
add_test(NAME risch_canonical_wiring_tests COMMAND risch_canonical_wiring_tests)

add_executable(risch_residue_split_tests test_risch_residue_split.c ${COMMON_SRC})
target_include_directories(risch_residue_split_tests PRIVATE ../src ../src/calculus)
add_test(NAME risch_residue_split_tests COMMAND risch_residue_split_tests)

add_executable(risch_elementaryq_tests test_risch_elementaryq.c ${COMMON_SRC})
target_include_directories(risch_elementaryq_tests PRIVATE ../src ../src/calculus)
add_test(NAME risch_elementaryq_tests COMMAND risch_elementaryq_tests)

add_executable(risch_rde_tower_tests test_risch_rde_tower.c ${COMMON_SRC})
target_include_directories(risch_rde_tower_tests PRIVATE ../src ../src/calculus)
add_test(NAME risch_rde_tower_tests COMMAND risch_rde_tower_tests)

add_executable(risch_structure_tests test_risch_structure.c ${COMMON_SRC})
target_include_directories(risch_structure_tests PRIVATE ../src ../src/calculus)
add_test(NAME risch_structure_tests COMMAND risch_structure_tests)

add_executable(risch_structure_real_tests test_risch_structure_real.c ${COMMON_SRC})
target_include_directories(risch_structure_real_tests PRIVATE ../src ../src/calculus)
add_test(NAME risch_structure_real_tests COMMAND risch_structure_real_tests)

add_executable(risch_logderiv_radical_tests test_risch_logderiv_radical.c ${COMMON_SRC})
target_include_directories(risch_logderiv_radical_tests PRIVATE ../src ../src/calculus)
add_test(NAME risch_logderiv_radical_tests COMMAND risch_logderiv_radical_tests)

add_executable(risch_hypertangent_tests test_risch_hypertangent.c ${COMMON_SRC})
target_include_directories(risch_hypertangent_tests PRIVATE ../src ../src/calculus)
add_test(NAME risch_hypertangent_tests COMMAND risch_hypertangent_tests)

add_executable(risch_coupled_tests test_risch_coupled.c ${COMMON_SRC})
target_include_directories(risch_coupled_tests PRIVATE ../src ../src/calculus)
add_test(NAME risch_coupled_tests COMMAND risch_coupled_tests)

add_executable(risch_hermite_tests test_risch_hermite.c ${COMMON_SRC})
target_include_directories(risch_hermite_tests PRIVATE ../src ../src/calculus)
add_test(NAME risch_hermite_tests COMMAND risch_hermite_tests)

add_executable(subresultants_tests test_subresultants.c ${COMMON_SRC})
target_include_directories(subresultants_tests PRIVATE ../src)

add_executable(subresultantpolynomials_tests test_subresultantpolynomials.c ${COMMON_SRC})
target_include_directories(subresultantpolynomials_tests PRIVATE ../src)

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

add_executable(squarefreeq_tests test_squarefreeq.c ${COMMON_SRC})
target_include_directories(squarefreeq_tests PRIVATE ../src)

add_executable(polynomialsqrt_tests test_polynomialsqrt.c ${COMMON_SRC})
target_include_directories(polynomialsqrt_tests PRIVATE ../src)
add_test(NAME polynomialsqrt_tests COMMAND polynomialsqrt_tests)

add_executable(rt_resultant_tests test_rt_resultant.c ${COMMON_SRC})
target_include_directories(rt_resultant_tests PRIVATE ../src ../src/calculus)
add_test(NAME rt_resultant_tests COMMAND rt_resultant_tests)

add_executable(cherry_ei_tests test_cherry_ei.c ${COMMON_SRC})
target_include_directories(cherry_ei_tests PRIVATE ../src ../src/calculus)
add_test(NAME cherry_ei_tests COMMAND cherry_ei_tests)

add_executable(cherry_stress_tests test_cherry_stress.c ${COMMON_SRC})
target_include_directories(cherry_stress_tests PRIVATE ../src ../src/calculus)
add_test(NAME cherry_stress_tests COMMAND cherry_stress_tests)

add_executable(cherry_li_tests test_cherry_li.c ${COMMON_SRC})
target_include_directories(cherry_li_tests PRIVATE ../src ../src/calculus)
add_test(NAME cherry_li_tests COMMAND cherry_li_tests)

add_executable(cherry_sigma_tests test_cherry_sigma.c ${COMMON_SRC})
target_include_directories(cherry_sigma_tests PRIVATE ../src)
add_test(NAME cherry_sigma_tests COMMAND cherry_sigma_tests)

add_executable(cherry_dilog_tests test_cherry_dilog.c ${COMMON_SRC})
target_include_directories(cherry_dilog_tests PRIVATE ../src ../src/calculus)
add_test(NAME cherry_dilog_tests COMMAND cherry_dilog_tests)

add_executable(integrate_fresnel_tests test_integrate_fresnel.c ${COMMON_SRC})
target_include_directories(integrate_fresnel_tests PRIVATE ../src ../src/calculus)
add_test(NAME integrate_fresnel_tests COMMAND integrate_fresnel_tests)

add_executable(irreduciblepolynomialq_tests test_irreduciblepolynomialq.c ${COMMON_SRC})
target_include_directories(irreduciblepolynomialq_tests PRIVATE ../src)

add_executable(solve_tests test_solve.c ${COMMON_SRC})
target_include_directories(solve_tests PRIVATE ../src)

add_executable(solvealways_tests test_solvealways.c ${COMMON_SRC})
target_include_directories(solvealways_tests PRIVATE ../src)

add_executable(solveradicals_tests test_solveradicals.c ${COMMON_SRC})
target_include_directories(solveradicals_tests PRIVATE ../src)

add_executable(solvenlsys_tests test_solvenlsys.c ${COMMON_SRC})
target_include_directories(solvenlsys_tests PRIVATE ../src)

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
    target_link_libraries(zupoly_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(zupoly_tests gmp m)
endif()

# Phase G1: Q(α) algebraic-number type (Trager's algebraic factoring).
add_executable(qa_tests test_qa.c ${COMMON_SRC})
target_include_directories(qa_tests PRIVATE ../src /usr/local/include)
target_link_directories(qa_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(qa_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(qa_tests gmp m)
endif()

# Phase G2: Q(α)[x] univariate polynomial type.
add_executable(qaupoly_tests test_qaupoly.c ${COMMON_SRC})
target_include_directories(qaupoly_tests PRIVATE ../src /usr/local/include)
target_link_directories(qaupoly_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(qaupoly_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(qaupoly_tests gmp m)
endif()

# Phase G3: Norm via resultant (qaupoly_norm).
add_executable(qafactor_tests test_qafactor.c ${COMMON_SRC})
target_include_directories(qafactor_tests PRIVATE ../src /usr/local/include)
target_link_directories(qafactor_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(qafactor_tests gmp m ${ECM_LIB})
    add_dependencies(qafactor_tests build_ecm)
else()
    target_link_libraries(qafactor_tests gmp m)
endif()

# Phase 1 (SIMPLIFY_IMPROVEMENT_PLAN): cyclotomic-extension primitives.
add_executable(rootofunity_tests test_rootofunity.c ${COMMON_SRC})
target_include_directories(rootofunity_tests PRIVATE ../src /usr/local/include)
target_link_directories(rootofunity_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(rootofunity_tests gmp m ${ECM_LIB})
    add_dependencies(rootofunity_tests build_ecm)
else()
    target_link_libraries(rootofunity_tests gmp m)
endif()

# Phase G9: extension_autodetect — walk an expression and build its
# algebraic-extension compositum automatically.
add_executable(extension_autodetect_tests test_extension_autodetect.c ${COMMON_SRC})
target_include_directories(extension_autodetect_tests PRIVATE ../src /usr/local/include)
target_link_directories(extension_autodetect_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(extension_autodetect_tests gmp m ${ECM_LIB})
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
    target_link_libraries(extension_auto_builtins_tests gmp m ${ECM_LIB})
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
    target_link_libraries(simp_algebraic_cuberoot_tests gmp m ${ECM_LIB})
    add_dependencies(simp_algebraic_cuberoot_tests build_ecm)
else()
    target_link_libraries(simp_algebraic_cuberoot_tests gmp m)
endif()

add_executable(bpoly_tests test_bpoly.c ${COMMON_SRC})
target_include_directories(bpoly_tests PRIVATE ../src /usr/local/include)
target_link_directories(bpoly_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(bpoly_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(bpoly_tests gmp m)
endif()

add_executable(mpoly_tests test_mpoly.c ${COMMON_SRC})
target_include_directories(mpoly_tests PRIVATE ../src /usr/local/include)
target_link_directories(mpoly_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(mpoly_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(mpoly_tests gmp m)
endif()

add_executable(mvfactor3_tests test_mvfactor3.c ${COMMON_SRC})
target_include_directories(mvfactor3_tests PRIVATE ../src /usr/local/include)
target_link_directories(mvfactor3_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(mvfactor3_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(mvfactor3_tests gmp m)
endif()

add_executable(mvfactor_tests test_mvfactor.c ${COMMON_SRC})
target_include_directories(mvfactor_tests PRIVATE ../src /usr/local/include)
target_link_directories(mvfactor_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(mvfactor_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(mvfactor_tests gmp m)
endif()

add_executable(factor_terms_tests test_factor_terms.c ${COMMON_SRC})
target_include_directories(factor_terms_tests PRIVATE ../src)

add_executable(rat_tests test_rat.c ${COMMON_SRC})
target_include_directories(rat_tests PRIVATE ../src)

add_executable(linalg_tests test_linalg.c ${COMMON_SRC})
target_include_directories(linalg_tests PRIVATE ../src)

add_executable(hilbertmatrix_tests test_hilbertmatrix.c ${COMMON_SRC})
target_include_directories(hilbertmatrix_tests PRIVATE ../src)
if(USE_ECM)
    add_dependencies(hilbertmatrix_tests build_ecm)
endif()
add_test(NAME hilbertmatrix_tests COMMAND hilbertmatrix_tests)

add_executable(hankelmatrix_tests test_hankelmatrix.c ${COMMON_SRC})
target_include_directories(hankelmatrix_tests PRIVATE ../src)
if(USE_ECM)
    add_dependencies(hankelmatrix_tests build_ecm)
endif()
add_test(NAME hankelmatrix_tests COMMAND hankelmatrix_tests)

add_executable(toeplitzmatrix_tests test_toeplitzmatrix.c ${COMMON_SRC})
target_include_directories(toeplitzmatrix_tests PRIVATE ../src)
if(USE_ECM)
    add_dependencies(toeplitzmatrix_tests build_ecm)
endif()
add_test(NAME toeplitzmatrix_tests COMMAND toeplitzmatrix_tests)

add_executable(vandermondematrix_tests test_vandermondematrix.c ${COMMON_SRC})
target_include_directories(vandermondematrix_tests PRIVATE ../src)
if(USE_ECM)
    add_dependencies(vandermondematrix_tests build_ecm)
endif()
add_test(NAME vandermondematrix_tests COMMAND vandermondematrix_tests)

add_executable(nullspace_tests test_nullspace.c ${COMMON_SRC})
target_include_directories(nullspace_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(nullspace_tests gmp m ${ECM_LIB})
    add_dependencies(nullspace_tests build_ecm)
else()
    target_link_libraries(nullspace_tests gmp m)
endif()
add_test(NAME nullspace_tests COMMAND nullspace_tests)

add_executable(matrank_tests test_matrank.c ${COMMON_SRC})
target_include_directories(matrank_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(matrank_tests gmp m ${ECM_LIB})
    add_dependencies(matrank_tests build_ecm)
else()
    target_link_libraries(matrank_tests gmp m)
endif()
add_test(NAME matrank_tests COMMAND matrank_tests)

add_executable(latticereduce_tests test_latticereduce.c ${COMMON_SRC})
target_include_directories(latticereduce_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(latticereduce_tests gmp m ${ECM_LIB})
    add_dependencies(latticereduce_tests build_ecm)
else()
    target_link_libraries(latticereduce_tests gmp m)
endif()
add_test(NAME latticereduce_tests COMMAND latticereduce_tests)

add_executable(minimalpolynomial_tests test_minimalpolynomial.c ${COMMON_SRC})
target_include_directories(minimalpolynomial_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(minimalpolynomial_tests gmp m ${ECM_LIB})
    add_dependencies(minimalpolynomial_tests build_ecm)
else()
    target_link_libraries(minimalpolynomial_tests gmp m)
endif()
add_test(NAME minimalpolynomial_tests COMMAND minimalpolynomial_tests)

add_executable(interp_tests test_interp.c ${COMMON_SRC})
target_include_directories(interp_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(interp_tests gmp m ${ECM_LIB})
    add_dependencies(interp_tests build_ecm)
else()
    target_link_libraries(interp_tests gmp m)
endif()
add_test(NAME interp_tests COMMAND interp_tests)

add_executable(findintegernullvector_tests test_findintegernullvector.c ${COMMON_SRC})
target_include_directories(findintegernullvector_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(findintegernullvector_tests gmp m ${ECM_LIB})
    add_dependencies(findintegernullvector_tests build_ecm)
else()
    target_link_libraries(findintegernullvector_tests gmp m)
endif()
add_test(NAME findintegernullvector_tests COMMAND findintegernullvector_tests)

add_executable(qrdecomposition_tests test_qrdecomposition.c ${COMMON_SRC})
target_include_directories(qrdecomposition_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(qrdecomposition_tests gmp m ${ECM_LIB})
    add_dependencies(qrdecomposition_tests build_ecm)
else()
    target_link_libraries(qrdecomposition_tests gmp m)
endif()
add_test(NAME qrdecomposition_tests COMMAND qrdecomposition_tests)

# Phase-3 LAPACK machine-precision QR kernel coverage.  The shape /
# aspect-ratio / pivoting / fuzz matrix from tasks/qr_perf_plan.md §5.2.
# Same source set as qrdecomposition_tests; routed through qr_dispatch
# which picks the LAPACK kernel when USE_LAPACK is on.
add_executable(qrdecomposition_machine_tests test_qrdecomposition_machine.c ${COMMON_SRC})
target_include_directories(qrdecomposition_machine_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(qrdecomposition_machine_tests gmp m ${ECM_LIB})
    add_dependencies(qrdecomposition_machine_tests build_ecm)
else()
    target_link_libraries(qrdecomposition_machine_tests gmp m)
endif()
add_test(NAME qrdecomposition_machine_tests COMMAND qrdecomposition_machine_tests)

# Phase-4 MPFR Householder kernel tests.  Cover the > 53-bit precision
# axis from tasks/qr_perf_plan.md §5.3.  Skip-cleanly under USE_MPFR=OFF.
add_executable(qrdecomposition_mpfr_tests test_qrdecomposition_mpfr.c ${COMMON_SRC})
target_include_directories(qrdecomposition_mpfr_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(qrdecomposition_mpfr_tests gmp m ${ECM_LIB})
    add_dependencies(qrdecomposition_mpfr_tests build_ecm)
else()
    target_link_libraries(qrdecomposition_mpfr_tests gmp m)
endif()
add_test(NAME qrdecomposition_mpfr_tests COMMAND qrdecomposition_mpfr_tests)

# LUDecomposition unit tests -- symbolic core + the public {lu, p, c}
# contract.  Covers integer / rational / complex / symbolic / singular
# inputs through the symbolic Doolittle kernel.
add_executable(ludecomposition_tests test_ludecomposition.c ${COMMON_SRC})
target_include_directories(ludecomposition_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(ludecomposition_tests gmp m ${ECM_LIB})
    add_dependencies(ludecomposition_tests build_ecm)
else()
    target_link_libraries(ludecomposition_tests gmp m)
endif()
add_test(NAME ludecomposition_tests COMMAND ludecomposition_tests)

# LUDecomposition LAPACK machine-precision kernel coverage.
add_executable(ludecomposition_machine_tests test_ludecomposition_machine.c ${COMMON_SRC})
target_include_directories(ludecomposition_machine_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(ludecomposition_machine_tests gmp m ${ECM_LIB})
    add_dependencies(ludecomposition_machine_tests build_ecm)
else()
    target_link_libraries(ludecomposition_machine_tests gmp m)
endif()
add_test(NAME ludecomposition_machine_tests COMMAND ludecomposition_machine_tests)

# LUDecomposition MPFR kernel tests.  Skip-cleanly under USE_MPFR=OFF.
add_executable(ludecomposition_mpfr_tests test_ludecomposition_mpfr.c ${COMMON_SRC})
target_include_directories(ludecomposition_mpfr_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(ludecomposition_mpfr_tests gmp m ${ECM_LIB})
    add_dependencies(ludecomposition_mpfr_tests build_ecm)
else()
    target_link_libraries(ludecomposition_mpfr_tests gmp m)
endif()
add_test(NAME ludecomposition_mpfr_tests COMMAND ludecomposition_mpfr_tests)

# SingularValueDecomposition symbolic-path tests.
add_executable(singularvaluedecomposition_tests test_singularvaluedecomposition.c ${COMMON_SRC})
target_include_directories(singularvaluedecomposition_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(singularvaluedecomposition_tests gmp m ${ECM_LIB})
    add_dependencies(singularvaluedecomposition_tests build_ecm)
else()
    target_link_libraries(singularvaluedecomposition_tests gmp m)
endif()
add_test(NAME singularvaluedecomposition_tests COMMAND singularvaluedecomposition_tests)

# SingularValueDecomposition LAPACK (machine-precision) tests.
add_executable(singularvaluedecomposition_machine_tests test_singularvaluedecomposition_machine.c ${COMMON_SRC})
target_include_directories(singularvaluedecomposition_machine_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(singularvaluedecomposition_machine_tests gmp m ${ECM_LIB})
    add_dependencies(singularvaluedecomposition_machine_tests build_ecm)
else()
    target_link_libraries(singularvaluedecomposition_machine_tests gmp m)
endif()
add_test(NAME singularvaluedecomposition_machine_tests COMMAND singularvaluedecomposition_machine_tests)

# SingularValueDecomposition MPFR Jacobi tests.  Skip-cleanly under USE_MPFR=OFF.
add_executable(singularvaluedecomposition_mpfr_tests test_singularvaluedecomposition_mpfr.c ${COMMON_SRC})
target_include_directories(singularvaluedecomposition_mpfr_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(singularvaluedecomposition_mpfr_tests gmp m ${ECM_LIB})
    add_dependencies(singularvaluedecomposition_mpfr_tests build_ecm)
else()
    target_link_libraries(singularvaluedecomposition_mpfr_tests gmp m)
endif()
add_test(NAME singularvaluedecomposition_mpfr_tests COMMAND singularvaluedecomposition_mpfr_tests)

# Phase-5 QR microbench.  NOT wired into add_test -- this is a spot-check
# tool that prints per-shape x precision wall time, not a CI gate.  Build
# it with `make bench_qr` then run `./bench_qr` after a perf-relevant
# change.
add_executable(bench_qr EXCLUDE_FROM_ALL bench_qr.c ${COMMON_SRC})
target_include_directories(bench_qr PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(bench_qr gmp m ${ECM_LIB})
    add_dependencies(bench_qr build_ecm)
else()
    target_link_libraries(bench_qr gmp m)
endif()

# LAPACK linkage smoke test (Phase 1 of QR perf plan).  Confirms that
# the configured BLAS/LAPACK distribution is callable from C and that
# mathilda_lapack_probe() returns the expected value under both
# USE_LAPACK=ON and OFF.
add_executable(lapack_tests test_lapack.c ${COMMON_SRC})
target_include_directories(lapack_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(lapack_tests gmp m ${ECM_LIB})
    add_dependencies(lapack_tests build_ecm)
else()
    target_link_libraries(lapack_tests gmp m)
endif()
add_test(NAME lapack_tests COMMAND lapack_tests)

# BLAS` / LAPACK` context builtins (blas_bridge.c / lapack_bridge.c).
add_executable(blas_tests test_blas.c ${COMMON_SRC})
target_include_directories(blas_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(blas_tests gmp m ${ECM_LIB})
    add_dependencies(blas_tests build_ecm)
else()
    target_link_libraries(blas_tests gmp m)
endif()
add_test(NAME blas_tests COMMAND blas_tests)

add_executable(lapack_builtin_tests test_lapack_builtin.c ${COMMON_SRC})
target_include_directories(lapack_builtin_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(lapack_builtin_tests gmp m ${ECM_LIB})
    add_dependencies(lapack_builtin_tests build_ecm)
else()
    target_link_libraries(lapack_builtin_tests gmp m)
endif()
add_test(NAME lapack_builtin_tests COMMAND lapack_builtin_tests)

add_executable(conjugate_transpose_tests test_conjugate_transpose.c ${COMMON_SRC})
target_include_directories(conjugate_transpose_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(conjugate_transpose_tests gmp m ${ECM_LIB})
    add_dependencies(conjugate_transpose_tests build_ecm)
else()
    target_link_libraries(conjugate_transpose_tests gmp m)
endif()
add_test(NAME conjugate_transpose_tests COMMAND conjugate_transpose_tests)

add_executable(hermitian_matrix_q_tests test_hermitian_matrix_q.c ${COMMON_SRC})
target_include_directories(hermitian_matrix_q_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(hermitian_matrix_q_tests gmp m ${ECM_LIB})
    add_dependencies(hermitian_matrix_q_tests build_ecm)
else()
    target_link_libraries(hermitian_matrix_q_tests gmp m)
endif()
add_test(NAME hermitian_matrix_q_tests COMMAND hermitian_matrix_q_tests)

add_executable(symmetric_matrix_q_tests test_symmetric_matrix_q.c ${COMMON_SRC})
target_include_directories(symmetric_matrix_q_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(symmetric_matrix_q_tests gmp m ${ECM_LIB})
    add_dependencies(symmetric_matrix_q_tests build_ecm)
else()
    target_link_libraries(symmetric_matrix_q_tests gmp m)
endif()
add_test(NAME symmetric_matrix_q_tests COMMAND symmetric_matrix_q_tests)

add_executable(square_matrix_q_tests test_square_matrix_q.c ${COMMON_SRC})
target_include_directories(square_matrix_q_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(square_matrix_q_tests gmp m ${ECM_LIB})
    add_dependencies(square_matrix_q_tests build_ecm)
else()
    target_link_libraries(square_matrix_q_tests gmp m)
endif()
add_test(NAME square_matrix_q_tests COMMAND square_matrix_q_tests)

add_executable(diagonal_matrix_q_tests test_diagonal_matrix_q.c ${COMMON_SRC})
target_include_directories(diagonal_matrix_q_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(diagonal_matrix_q_tests gmp m ${ECM_LIB})
    add_dependencies(diagonal_matrix_q_tests build_ecm)
else()
    target_link_libraries(diagonal_matrix_q_tests gmp m)
endif()
add_test(NAME diagonal_matrix_q_tests COMMAND diagonal_matrix_q_tests)

add_executable(upper_triangular_matrix_q_tests test_upper_triangular_matrix_q.c ${COMMON_SRC})
target_include_directories(upper_triangular_matrix_q_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(upper_triangular_matrix_q_tests gmp m ${ECM_LIB})
    add_dependencies(upper_triangular_matrix_q_tests build_ecm)
else()
    target_link_libraries(upper_triangular_matrix_q_tests gmp m)
endif()
add_test(NAME upper_triangular_matrix_q_tests COMMAND upper_triangular_matrix_q_tests)

add_executable(positive_definite_matrix_q_tests test_positive_definite_matrix_q.c ${COMMON_SRC})
target_include_directories(positive_definite_matrix_q_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(positive_definite_matrix_q_tests gmp m ${ECM_LIB})
    add_dependencies(positive_definite_matrix_q_tests build_ecm)
else()
    target_link_libraries(positive_definite_matrix_q_tests gmp m)
endif()
add_test(NAME positive_definite_matrix_q_tests COMMAND positive_definite_matrix_q_tests)

add_executable(negative_definite_matrix_q_tests test_negative_definite_matrix_q.c ${COMMON_SRC})
target_include_directories(negative_definite_matrix_q_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(negative_definite_matrix_q_tests gmp m ${ECM_LIB})
    add_dependencies(negative_definite_matrix_q_tests build_ecm)
else()
    target_link_libraries(negative_definite_matrix_q_tests gmp m)
endif()
add_test(NAME negative_definite_matrix_q_tests COMMAND negative_definite_matrix_q_tests)

add_executable(normalize_tests test_normalize.c ${COMMON_SRC})
target_include_directories(normalize_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(normalize_tests gmp m ${ECM_LIB})
    add_dependencies(normalize_tests build_ecm)
else()
    target_link_libraries(normalize_tests gmp m)
endif()
add_test(NAME normalize_tests COMMAND normalize_tests)

add_executable(integer_digits_tests test_integer_digits.c ${COMMON_SRC})
target_include_directories(integer_digits_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integer_digits_tests gmp m ${ECM_LIB})
    add_dependencies(integer_digits_tests build_ecm)
else()
    target_link_libraries(integer_digits_tests gmp m)
endif()
add_test(NAME integer_digits_tests COMMAND integer_digits_tests)

add_executable(real_digits_tests test_real_digits.c ${COMMON_SRC})
target_include_directories(real_digits_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(real_digits_tests gmp m ${ECM_LIB})
    add_dependencies(real_digits_tests build_ecm)
else()
    target_link_libraries(real_digits_tests gmp m)
endif()
add_test(NAME real_digits_tests COMMAND real_digits_tests)

add_executable(mantissa_exponent_tests test_mantissa_exponent.c ${COMMON_SRC})
target_include_directories(mantissa_exponent_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(mantissa_exponent_tests gmp m ${ECM_LIB})
    add_dependencies(mantissa_exponent_tests build_ecm)
else()
    target_link_libraries(mantissa_exponent_tests gmp m)
endif()
add_test(NAME mantissa_exponent_tests COMMAND mantissa_exponent_tests)

add_executable(real_exponent_tests test_real_exponent.c ${COMMON_SRC})
target_include_directories(real_exponent_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(real_exponent_tests gmp m ${ECM_LIB})
    add_dependencies(real_exponent_tests build_ecm)
else()
    target_link_libraries(real_exponent_tests gmp m)
endif()
add_test(NAME real_exponent_tests COMMAND real_exponent_tests)

add_executable(parser_precision_tests test_parser_precision.c ${COMMON_SRC})
target_include_directories(parser_precision_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(parser_precision_tests gmp m ${ECM_LIB})
    add_dependencies(parser_precision_tests build_ecm)
else()
    target_link_libraries(parser_precision_tests gmp m)
endif()
add_test(NAME parser_precision_tests COMMAND parser_precision_tests)

add_executable(integer_length_tests test_integer_length.c ${COMMON_SRC})
target_include_directories(integer_length_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integer_length_tests gmp m ${ECM_LIB})
    add_dependencies(integer_length_tests build_ecm)
else()
    target_link_libraries(integer_length_tests gmp m)
endif()
add_test(NAME integer_length_tests COMMAND integer_length_tests)

add_executable(digit_count_tests test_digit_count.c ${COMMON_SRC})
target_include_directories(digit_count_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(digit_count_tests gmp m ${ECM_LIB})
    add_dependencies(digit_count_tests build_ecm)
else()
    target_link_libraries(digit_count_tests gmp m)
endif()
add_test(NAME digit_count_tests COMMAND digit_count_tests)

add_executable(digit_sum_tests test_digit_sum.c ${COMMON_SRC})
target_include_directories(digit_sum_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(digit_sum_tests gmp m ${ECM_LIB})
    add_dependencies(digit_sum_tests build_ecm)
else()
    target_link_libraries(digit_sum_tests gmp m)
endif()
add_test(NAME digit_sum_tests COMMAND digit_sum_tests)

add_executable(from_digits_tests test_from_digits.c ${COMMON_SRC})
target_include_directories(from_digits_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(from_digits_tests gmp m ${ECM_LIB})
    add_dependencies(from_digits_tests build_ecm)
else()
    target_link_libraries(from_digits_tests gmp m)
endif()
add_test(NAME from_digits_tests COMMAND from_digits_tests)

add_executable(integer_exponent_tests test_integer_exponent.c ${COMMON_SRC})
target_include_directories(integer_exponent_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integer_exponent_tests gmp m ${ECM_LIB})
    add_dependencies(integer_exponent_tests build_ecm)
else()
    target_link_libraries(integer_exponent_tests gmp m)
endif()
add_test(NAME integer_exponent_tests COMMAND integer_exponent_tests)

add_executable(extended_gcd_tests test_extended_gcd.c ${COMMON_SRC})
target_include_directories(extended_gcd_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(extended_gcd_tests gmp m ${ECM_LIB})
    add_dependencies(extended_gcd_tests build_ecm)
else()
    target_link_libraries(extended_gcd_tests gmp m)
endif()
add_test(NAME extended_gcd_tests COMMAND extended_gcd_tests)

add_executable(jacobisymbol_tests test_jacobisymbol.c ${COMMON_SRC})
target_include_directories(jacobisymbol_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(jacobisymbol_tests gmp m ${ECM_LIB})
    add_dependencies(jacobisymbol_tests build_ecm)
else()
    target_link_libraries(jacobisymbol_tests gmp m)
endif()
add_test(NAME jacobisymbol_tests COMMAND jacobisymbol_tests)

add_executable(integer_partitions_tests test_integer_partitions.c ${COMMON_SRC})
target_include_directories(integer_partitions_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integer_partitions_tests gmp m ${ECM_LIB})
    add_dependencies(integer_partitions_tests build_ecm)
else()
    target_link_libraries(integer_partitions_tests gmp m)
endif()
add_test(NAME integer_partitions_tests COMMAND integer_partitions_tests)

add_executable(partitionsp_tests test_partitionsp.c ${COMMON_SRC})
target_include_directories(partitionsp_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(partitionsp_tests gmp m ${ECM_LIB})
    add_dependencies(partitionsp_tests build_ecm)
else()
    target_link_libraries(partitionsp_tests gmp m)
endif()
add_test(NAME partitionsp_tests COMMAND partitionsp_tests)

add_executable(partitionsq_tests test_partitionsq.c ${COMMON_SRC})
target_include_directories(partitionsq_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(partitionsq_tests gmp m ${ECM_LIB})
    add_dependencies(partitionsq_tests build_ecm)
else()
    target_link_libraries(partitionsq_tests gmp m)
endif()
add_test(NAME partitionsq_tests COMMAND partitionsq_tests)

add_executable(integer_string_tests test_integer_string.c ${COMMON_SRC})
target_include_directories(integer_string_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integer_string_tests gmp m ${ECM_LIB})
    add_dependencies(integer_string_tests build_ecm)
else()
    target_link_libraries(integer_string_tests gmp m)
endif()
add_test(NAME integer_string_tests COMMAND integer_string_tests)

add_executable(machine_number_q_tests test_machine_number_q.c ${COMMON_SRC})
target_include_directories(machine_number_q_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(machine_number_q_tests gmp m ${ECM_LIB})
    add_dependencies(machine_number_q_tests build_ecm)
else()
    target_link_libraries(machine_number_q_tests gmp m)
endif()
add_test(NAME machine_number_q_tests COMMAND machine_number_q_tests)

add_executable(eigen_tests test_eigen.c ${COMMON_SRC})
target_include_directories(eigen_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(eigen_tests gmp m ${ECM_LIB})
    add_dependencies(eigen_tests build_ecm)
else()
    target_link_libraries(eigen_tests gmp m)
endif()
add_test(NAME eigen_tests COMMAND eigen_tests)

add_executable(mateigen_direct_tests test_mateigen_direct.c eigen_corpus.c ${COMMON_SRC})
target_include_directories(mateigen_direct_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(mateigen_direct_tests gmp m ${ECM_LIB})
    add_dependencies(mateigen_direct_tests build_ecm)
else()
    target_link_libraries(mateigen_direct_tests gmp m)
endif()
add_test(NAME mateigen_direct_tests COMMAND mateigen_direct_tests)

add_executable(mateigen_arnoldi_tests test_mateigen_arnoldi.c eigen_corpus.c ${COMMON_SRC})
target_include_directories(mateigen_arnoldi_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(mateigen_arnoldi_tests gmp m ${ECM_LIB})
    add_dependencies(mateigen_arnoldi_tests build_ecm)
else()
    target_link_libraries(mateigen_arnoldi_tests gmp m)
endif()
add_test(NAME mateigen_arnoldi_tests COMMAND mateigen_arnoldi_tests)

add_executable(mateigen_banded_tests test_mateigen_banded.c eigen_corpus.c ${COMMON_SRC})
target_include_directories(mateigen_banded_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(mateigen_banded_tests gmp m ${ECM_LIB})
    add_dependencies(mateigen_banded_tests build_ecm)
else()
    target_link_libraries(mateigen_banded_tests gmp m)
endif()
add_test(NAME mateigen_banded_tests COMMAND mateigen_banded_tests)

add_executable(mateigen_feast_tests test_mateigen_feast.c eigen_corpus.c ${COMMON_SRC})
target_include_directories(mateigen_feast_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(mateigen_feast_tests gmp m ${ECM_LIB})
    add_dependencies(mateigen_feast_tests build_ecm)
else()
    target_link_libraries(mateigen_feast_tests gmp m)
endif()
add_test(NAME mateigen_feast_tests COMMAND mateigen_feast_tests)

add_executable(linearsolve_tests test_linearsolve.c ${COMMON_SRC})
target_include_directories(linearsolve_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(linearsolve_tests gmp m ${ECM_LIB})
    add_dependencies(linearsolve_tests build_ecm)
else()
    target_link_libraries(linearsolve_tests gmp m)
endif()
add_test(NAME linearsolve_tests COMMAND linearsolve_tests)

add_executable(matsol_methods_tests test_matsol_methods.c ${COMMON_SRC})
target_include_directories(matsol_methods_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(matsol_methods_tests gmp m ${ECM_LIB})
    add_dependencies(matsol_methods_tests build_ecm)
else()
    target_link_libraries(matsol_methods_tests gmp m)
endif()
add_test(NAME matsol_methods_tests COMMAND matsol_methods_tests)

add_executable(matinv_tests test_matinv.c ${COMMON_SRC})
target_include_directories(matinv_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(matinv_tests gmp m ${ECM_LIB})
    add_dependencies(matinv_tests build_ecm)
else()
    target_link_libraries(matinv_tests gmp m)
endif()
add_test(NAME matinv_tests COMMAND matinv_tests)

add_executable(matinv_methods_tests test_matinv_methods.c ${COMMON_SRC})
target_include_directories(matinv_methods_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(matinv_methods_tests gmp m ${ECM_LIB})
    add_dependencies(matinv_methods_tests build_ecm)
else()
    target_link_libraries(matinv_methods_tests gmp m)
endif()
add_test(NAME matinv_methods_tests COMMAND matinv_methods_tests)

add_executable(matlstsq_tests test_matlstsq.c ${COMMON_SRC})
target_include_directories(matlstsq_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(matlstsq_tests gmp m ${ECM_LIB})
    add_dependencies(matlstsq_tests build_ecm)
else()
    target_link_libraries(matlstsq_tests gmp m)
endif()
add_test(NAME matlstsq_tests COMMAND matlstsq_tests)

add_executable(parfrac_tests test_parfrac.c ${COMMON_SRC})
target_include_directories(parfrac_tests PRIVATE ../src)

add_executable(contfrac_tests test_contfrac.c ${COMMON_SRC})
target_include_directories(contfrac_tests PRIVATE ../src)

add_executable(polymod_tests test_polymod.c ${COMMON_SRC})
target_include_directories(polymod_tests PRIVATE ../src)

add_executable(bigint_tests test_bigint.c ${COMMON_SRC})
target_include_directories(bigint_tests PRIVATE ../src /usr/local/include)
target_link_directories(bigint_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(bigint_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(bigint_tests gmp m)
endif()

add_executable(expr_sharing_tests test_expr_sharing.c ${COMMON_SRC})
target_include_directories(expr_sharing_tests PRIVATE ../src /usr/local/include)
target_link_directories(expr_sharing_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(expr_sharing_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(expr_sharing_tests gmp m)
endif()

add_executable(eval_timestamps_tests test_eval_timestamps.c ${COMMON_SRC})
target_include_directories(eval_timestamps_tests PRIVATE ../src /usr/local/include)
target_link_directories(eval_timestamps_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(eval_timestamps_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(eval_timestamps_tests gmp m)
endif()

add_executable(eval_eager_exit_tests test_eval_eager_exit.c ${COMMON_SRC})
target_include_directories(eval_eager_exit_tests PRIVATE ../src /usr/local/include)
target_link_directories(eval_eager_exit_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(eval_eager_exit_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(eval_eager_exit_tests gmp m)
endif()

add_executable(rule_dispatch_tests test_rule_dispatch.c ${COMMON_SRC})
target_include_directories(rule_dispatch_tests PRIVATE ../src /usr/local/include)
target_link_directories(rule_dispatch_tests PRIVATE /usr/local/lib)
if(USE_ECM)
    target_link_libraries(rule_dispatch_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(rule_dispatch_tests gmp m)
endif()

add_executable(integrals_tests test_integrals.c ${COMMON_SRC})
target_include_directories(integrals_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrals_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrals_tests gmp m)
endif()

add_executable(integrate_dispatch_tests test_integrate_dispatch.c ${COMMON_SRC})
target_include_directories(integrate_dispatch_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_dispatch_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_dispatch_tests gmp m)
endif()

add_executable(integrate_derivdivides_tests test_integrate_derivdivides.c ${COMMON_SRC})
target_include_directories(integrate_derivdivides_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_derivdivides_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_derivdivides_tests gmp m)
endif()
add_test(NAME integrate_derivdivides_tests COMMAND integrate_derivdivides_tests)

add_executable(integrate_risch_transcendental_tests test_integrate_risch_transcendental.c ${COMMON_SRC})
target_include_directories(integrate_risch_transcendental_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_risch_transcendental_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_risch_transcendental_tests gmp m)
endif()
add_test(NAME integrate_risch_transcendental_tests COMMAND integrate_risch_transcendental_tests)

add_executable(integrate_linrad_tests test_integrate_linrad.c ${COMMON_SRC})
target_include_directories(integrate_linrad_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_linrad_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_linrad_tests gmp m)
endif()
add_test(NAME integrate_linrad_tests COMMAND integrate_linrad_tests)

add_executable(integrate_quadrad_tests test_integrate_quadrad.c ${COMMON_SRC})
target_include_directories(integrate_quadrad_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_quadrad_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_quadrad_tests gmp m)
endif()
add_test(NAME integrate_quadrad_tests COMMAND integrate_quadrad_tests)

add_executable(integrate_linratiorad_tests test_integrate_linratiorad.c ${COMMON_SRC})
target_include_directories(integrate_linratiorad_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_linratiorad_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_linratiorad_tests gmp m)
endif()
add_test(NAME integrate_linratiorad_tests COMMAND integrate_linratiorad_tests)

add_executable(integrate_chebychev_tests test_integrate_chebychev.c ${COMMON_SRC})
target_include_directories(integrate_chebychev_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_chebychev_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_chebychev_tests gmp m)
endif()
add_test(NAME integrate_chebychev_tests COMMAND integrate_chebychev_tests)

add_executable(integrate_goursat_tests test_integrate_goursat.c ${COMMON_SRC})
target_include_directories(integrate_goursat_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_goursat_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_goursat_tests gmp m)
endif()
add_test(NAME integrate_goursat_tests COMMAND integrate_goursat_tests)

add_executable(integrate_newton_leibniz_tests test_integrate_newton_leibniz.c ${COMMON_SRC})
target_include_directories(integrate_newton_leibniz_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_newton_leibniz_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_newton_leibniz_tests gmp m)
endif()
add_test(NAME integrate_newton_leibniz_tests COMMAND integrate_newton_leibniz_tests)

add_executable(integrate_line_tests test_integrate_line.c ${COMMON_SRC})
target_include_directories(integrate_line_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_line_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_line_tests gmp m)
endif()
add_test(NAME integrate_line_tests COMMAND integrate_line_tests)

add_executable(integrate_residue_tests test_integrate_residue.c ${COMMON_SRC})
target_include_directories(integrate_residue_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_residue_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_residue_tests gmp m)
endif()
add_test(NAME integrate_residue_tests COMMAND integrate_residue_tests)

add_executable(integrate_symmetry_tests test_integrate_symmetry.c ${COMMON_SRC})
target_include_directories(integrate_symmetry_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_symmetry_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_symmetry_tests gmp m)
endif()
add_test(NAME integrate_symmetry_tests COMMAND integrate_symmetry_tests)

add_executable(integrate_principalvalue_tests test_integrate_principalvalue.c ${COMMON_SRC})
target_include_directories(integrate_principalvalue_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_principalvalue_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_principalvalue_tests gmp m)
endif()
add_test(NAME integrate_principalvalue_tests COMMAND integrate_principalvalue_tests)

add_executable(integrate_osclog_tests test_integrate_osclog.c ${COMMON_SRC})
target_include_directories(integrate_osclog_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_osclog_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_osclog_tests gmp m)
endif()
add_test(NAME integrate_osclog_tests COMMAND integrate_osclog_tests)

add_executable(integrate_beta_tests test_integrate_beta.c ${COMMON_SRC})
target_include_directories(integrate_beta_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_beta_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_beta_tests gmp m)
endif()
add_test(NAME integrate_beta_tests COMMAND integrate_beta_tests)

add_executable(integrate_diffunderint_tests test_integrate_diffunderint.c ${COMMON_SRC})
target_include_directories(integrate_diffunderint_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_diffunderint_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_diffunderint_tests gmp m)
endif()
add_test(NAME integrate_diffunderint_tests COMMAND integrate_diffunderint_tests)

add_executable(integrate_ramanujan_tests test_integrate_ramanujan.c ${COMMON_SRC})
target_include_directories(integrate_ramanujan_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_ramanujan_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_ramanujan_tests gmp m)
endif()
add_test(NAME integrate_ramanujan_tests COMMAND integrate_ramanujan_tests)

add_executable(residue_tests test_residue.c ${COMMON_SRC})
target_include_directories(residue_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(residue_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(residue_tests gmp m)
endif()
add_test(NAME residue_tests COMMAND residue_tests)

add_executable(integrate_jeffrey_tests test_integrate_jeffrey.c ${COMMON_SRC})
target_include_directories(integrate_jeffrey_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_jeffrey_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(integrate_jeffrey_tests gmp m)
endif()
add_test(NAME integrate_jeffrey_tests COMMAND integrate_jeffrey_tests)

add_executable(condition_downvalue_tests test_condition_downvalue.c ${COMMON_SRC})
target_include_directories(condition_downvalue_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(condition_downvalue_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(condition_downvalue_tests gmp m)
endif()

add_executable(print_tests test_print.c ${COMMON_SRC})
target_include_directories(print_tests PRIVATE ../src)

add_executable(graphics_sampling_tests test_graphics_sampling.c ${COMMON_SRC})
target_link_libraries(graphics_sampling_tests m)
target_include_directories(graphics_sampling_tests PRIVATE ../src ../src/graphics)
add_test(NAME graphics_sampling_tests COMMAND graphics_sampling_tests)

add_executable(graphics_tests test_graphics.c ${COMMON_SRC})
target_include_directories(graphics_tests PRIVATE ../src ../src/graphics)
add_test(NAME graphics_tests COMMAND graphics_tests)

add_executable(plot3d_tests test_plot3d.c ${COMMON_SRC})
target_include_directories(plot3d_tests PRIVATE ../src ../src/graphics)
add_test(NAME plot3d_tests COMMAND plot3d_tests)

add_executable(parametricplot_tests test_parametricplot.c ${COMMON_SRC})
target_link_libraries(parametricplot_tests m)
target_include_directories(parametricplot_tests PRIVATE ../src ../src/graphics)
add_test(NAME parametricplot_tests COMMAND parametricplot_tests)

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

add_executable(regex_tests test_regex.c ${COMMON_SRC})
target_include_directories(regex_tests PRIVATE ../src ../src/strings ../src/strings/regex)

add_executable(stringfns_tests test_stringfns.c ${COMMON_SRC})
target_include_directories(stringfns_tests PRIVATE ../src ../src/strings ../src/strings/regex)

add_executable(stringpartition_tests test_stringpartition.c ${COMMON_SRC})
target_include_directories(stringpartition_tests PRIVATE ../src ../src/strings ../src/strings/regex)
if(USE_ECM)
    target_link_libraries(stringpartition_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(stringpartition_tests gmp m)
endif()
add_test(NAME stringpartition_tests COMMAND stringpartition_tests)

add_executable(stringrepeat_tests test_stringrepeat.c ${COMMON_SRC})
target_include_directories(stringrepeat_tests PRIVATE ../src ../src/strings ../src/strings/regex)
if(USE_ECM)
    target_link_libraries(stringrepeat_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(stringrepeat_tests gmp m)
endif()
add_test(NAME stringrepeat_tests COMMAND stringrepeat_tests)

add_executable(series_tests test_series.c ${COMMON_SRC})
target_include_directories(series_tests PRIVATE ../src)

add_executable(series_twoterm_tests test_series_twoterm.c ${COMMON_SRC})
target_include_directories(series_twoterm_tests PRIVATE ../src)

add_executable(deriv_tests test_deriv.c ${COMMON_SRC})
target_include_directories(deriv_tests PRIVATE ../src)

add_executable(deriv_symbolic_order_tests test_deriv_symbolic_order.c ${COMMON_SRC})
target_include_directories(deriv_symbolic_order_tests PRIVATE ../src)

add_executable(deriv_array_tests test_deriv_array.c ${COMMON_SRC})
target_include_directories(deriv_array_tests PRIVATE ../src)

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

add_executable(gamma_tests test_gamma.c ${COMMON_SRC})
target_include_directories(gamma_tests PRIVATE ../src)

add_executable(beta_tests test_beta.c ${COMMON_SRC})
target_include_directories(beta_tests PRIVATE ../src)

add_executable(erf_tests test_erf.c ${COMMON_SRC})
target_include_directories(erf_tests PRIVATE ../src)

add_executable(erfc_tests test_erfc.c ${COMMON_SRC})
target_include_directories(erfc_tests PRIVATE ../src)

add_executable(erfi_tests test_erfi.c ${COMMON_SRC})
target_include_directories(erfi_tests PRIVATE ../src)

add_executable(expintegralei_tests test_expintegralei.c ${COMMON_SRC})
target_include_directories(expintegralei_tests PRIVATE ../src)

add_executable(logintegral_tests test_logintegral.c ${COMMON_SRC})
target_include_directories(logintegral_tests PRIVATE ../src)

add_executable(sinintegral_tests test_sinintegral.c ${COMMON_SRC})
target_include_directories(sinintegral_tests PRIVATE ../src)

add_executable(cosintegral_tests test_cosintegral.c ${COMMON_SRC})
target_include_directories(cosintegral_tests PRIVATE ../src)

add_executable(sinhintegral_tests test_sinhintegral.c ${COMMON_SRC})
target_include_directories(sinhintegral_tests PRIVATE ../src)

add_executable(coshintegral_tests test_coshintegral.c ${COMMON_SRC})
target_include_directories(coshintegral_tests PRIVATE ../src)

add_executable(fresnelc_tests test_fresnelc.c ${COMMON_SRC})
target_include_directories(fresnelc_tests PRIVATE ../src)

add_executable(fresnels_tests test_fresnels.c ${COMMON_SRC})
target_include_directories(fresnels_tests PRIVATE ../src)

add_executable(sinc_tests test_sinc.c ${COMMON_SRC})
target_include_directories(sinc_tests PRIVATE ../src)

add_executable(positive_tests test_positive.c ${COMMON_SRC})
target_include_directories(positive_tests PRIVATE ../src)

add_executable(negative_tests test_negative.c ${COMMON_SRC})
target_include_directories(negative_tests PRIVATE ../src)

add_executable(nonnegative_tests test_nonnegative.c ${COMMON_SRC})
target_include_directories(nonnegative_tests PRIVATE ../src)

add_executable(nonpositive_tests test_nonpositive.c ${COMMON_SRC})
target_include_directories(nonpositive_tests PRIVATE ../src)

add_executable(inverf_tests test_inverf.c ${COMMON_SRC})
target_include_directories(inverf_tests PRIVATE ../src)

add_executable(inverfc_tests test_inverfc.c ${COMMON_SRC})
target_include_directories(inverfc_tests PRIVATE ../src)

add_executable(polygamma_tests test_polygamma.c ${COMMON_SRC})
target_include_directories(polygamma_tests PRIVATE ../src)

add_executable(harmonicnumber_tests test_harmonicnumber.c ${COMMON_SRC})
target_include_directories(harmonicnumber_tests PRIVATE ../src)

add_executable(pochhammer_tests test_pochhammer.c ${COMMON_SRC})
target_include_directories(pochhammer_tests PRIVATE ../src)

add_executable(zeta_tests test_zeta.c ${COMMON_SRC})
target_include_directories(zeta_tests PRIVATE ../src)

add_executable(hurwitzzeta_tests test_hurwitzzeta.c ${COMMON_SRC})
target_include_directories(hurwitzzeta_tests PRIVATE ../src)

add_executable(stieltjesgamma_tests test_stieltjesgamma.c ${COMMON_SRC})
target_include_directories(stieltjesgamma_tests PRIVATE ../src)

add_executable(bernoullib_tests test_bernoullib.c ${COMMON_SRC})
target_include_directories(bernoullib_tests PRIVATE ../src)

add_executable(eulere_tests test_eulere.c ${COMMON_SRC})
target_include_directories(eulere_tests PRIVATE ../src)

add_executable(polylog_tests test_polylog.c ${COMMON_SRC})
target_include_directories(polylog_tests PRIVATE ../src)

add_executable(lerchphi_tests test_lerchphi.c ${COMMON_SRC})
target_include_directories(lerchphi_tests PRIVATE ../src)

add_executable(productlog_tests test_productlog.c ${COMMON_SRC})
target_include_directories(productlog_tests PRIVATE ../src)

add_executable(legendre_tests test_legendre.c ${COMMON_SRC})
target_include_directories(legendre_tests PRIVATE ../src)

add_executable(airyai_tests test_airyai.c ${COMMON_SRC})
target_include_directories(airyai_tests PRIVATE ../src)

add_executable(airybi_tests test_airybi.c ${COMMON_SRC})
target_include_directories(airybi_tests PRIVATE ../src)

add_executable(besselj_tests test_besselj.c ${COMMON_SRC})
target_include_directories(besselj_tests PRIVATE ../src)

add_executable(besselk_tests test_besselk.c ${COMMON_SRC})
target_include_directories(besselk_tests PRIVATE ../src)

add_executable(besseli_tests test_besseli.c ${COMMON_SRC})
target_include_directories(besseli_tests PRIVATE ../src)

add_executable(bessely_tests test_bessely.c ${COMMON_SRC})
target_include_directories(bessely_tests PRIVATE ../src)

add_executable(loggamma_tests test_loggamma.c ${COMMON_SRC})
target_include_directories(loggamma_tests PRIVATE ../src)

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

add_executable(simp_log_tests test_simp_log.c ${COMMON_SRC})
target_include_directories(simp_log_tests PRIVATE ../src)

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

add_executable(numeric_complex_tests test_numeric_complex.c ${COMMON_SRC})
target_include_directories(numeric_complex_tests PRIVATE ../src)

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
    add_dependencies(catch_throw_tests build_ecm)
    add_dependencies(boolean_tests build_ecm)
    add_dependencies(names_tests build_ecm)
    add_dependencies(association_tests build_ecm)
    add_dependencies(bench_assoc build_ecm)
    add_dependencies(bench_ndarray_linalg build_ecm)
    add_dependencies(ndarray_tests build_ecm)
    add_dependencies(ndarray_linalg_tests build_ecm)
    add_dependencies(patterns_tests build_ecm)
    add_dependencies(replace_tests build_ecm)
    add_dependencies(sequence_tests build_ecm)
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
    add_dependencies(subresultants_tests build_ecm)
    add_dependencies(subresultantpolynomials_tests build_ecm)
    add_dependencies(expand_tests build_ecm)
    add_dependencies(expandfrac_tests build_ecm)
    add_dependencies(stats_tests build_ecm)
    add_dependencies(fit_tests build_ecm)
    add_dependencies(datetime_tests build_ecm)
    add_dependencies(facpoly_tests build_ecm)
    add_dependencies(factor_terms_tests build_ecm)
    add_dependencies(rat_tests build_ecm)
    add_dependencies(linalg_tests build_ecm)
    add_dependencies(parfrac_tests build_ecm)
    add_dependencies(contfrac_tests build_ecm)
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
    add_dependencies(regex_tests build_ecm)
    add_dependencies(stringfns_tests build_ecm)
    add_dependencies(stringpartition_tests build_ecm)
    add_dependencies(stringrepeat_tests build_ecm)
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
    add_dependencies(simp_log_tests build_ecm)
    add_dependencies(invtrig_simplify_tests build_ecm)
    add_dependencies(radical_simplify_tests build_ecm)
    add_dependencies(rationalize_tests build_ecm)
    add_dependencies(inexact_dispatch_tests build_ecm)
    add_dependencies(numeric_complex_tests build_ecm)
endif()

add_executable(simp_tests test_simp.c ${COMMON_SRC})
target_include_directories(simp_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(simp_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(simp_tests gmp m)
endif()
add_test(NAME simp_tests COMMAND simp_tests)

add_executable(fibonacci_tests test_fibonacci.c ${COMMON_SRC})
target_include_directories(fibonacci_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(fibonacci_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(fibonacci_tests gmp m)
endif()
add_test(NAME fibonacci_tests COMMAND fibonacci_tests)

add_executable(lucas_tests test_lucas.c ${COMMON_SRC})
target_include_directories(lucas_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(lucas_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(lucas_tests gmp m)
endif()
add_test(NAME lucas_tests COMMAND lucas_tests)

add_executable(distribute_tests test_distribute.c ${COMMON_SRC})
target_include_directories(distribute_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(distribute_tests gmp m ${ECM_LIB})
else()
    target_link_libraries(distribute_tests gmp m)
endif()
add_test(NAME distribute_tests COMMAND distribute_tests)

add_executable(thread_tests test_thread.c ${COMMON_SRC})
target_include_directories(thread_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(thread_tests gmp m ${ECM_LIB})
    add_dependencies(thread_tests build_ecm)
else()
    target_link_libraries(thread_tests gmp m)
endif()
add_test(NAME thread_tests COMMAND thread_tests)

add_executable(mapthread_tests test_mapthread.c ${COMMON_SRC})
target_include_directories(mapthread_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(mapthread_tests gmp m ${ECM_LIB})
    add_dependencies(mapthread_tests build_ecm)
else()
    target_link_libraries(mapthread_tests gmp m)
endif()
add_test(NAME mapthread_tests COMMAND mapthread_tests)

add_executable(trigexpand_tests test_trigexpand.c ${COMMON_SRC})
target_include_directories(trigexpand_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(trigexpand_tests gmp m ${ECM_LIB})
    add_dependencies(trigexpand_tests build_ecm)
else()
    target_link_libraries(trigexpand_tests gmp m)
endif()
add_test(NAME trigexpand_tests COMMAND trigexpand_tests)

add_executable(complexexpand_tests test_complexexpand.c ${COMMON_SRC})
target_include_directories(complexexpand_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(complexexpand_tests gmp m ${ECM_LIB})
    add_dependencies(complexexpand_tests build_ecm)
else()
    target_link_libraries(complexexpand_tests gmp m)
endif()
add_test(NAME complexexpand_tests COMMAND complexexpand_tests)

add_executable(trigfactor_tests test_trigfactor.c ${COMMON_SRC})
target_include_directories(trigfactor_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(trigfactor_tests gmp m ${ECM_LIB})
    add_dependencies(trigfactor_tests build_ecm)
else()
    target_link_libraries(trigfactor_tests gmp m)
endif()
add_test(NAME trigfactor_tests COMMAND trigfactor_tests)

add_executable(trigreduce_tests test_trigreduce.c ${COMMON_SRC})
target_include_directories(trigreduce_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(trigreduce_tests gmp m ${ECM_LIB})
    add_dependencies(trigreduce_tests build_ecm)
else()
    target_link_libraries(trigreduce_tests gmp m)
endif()
add_test(NAME trigreduce_tests COMMAND trigreduce_tests)

add_executable(recursion_limit_tests test_recursion_limit.c ${COMMON_SRC})
target_include_directories(recursion_limit_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(recursion_limit_tests gmp m ${ECM_LIB})
    add_dependencies(recursion_limit_tests build_ecm)
else()
    target_link_libraries(recursion_limit_tests gmp m)
endif()
add_test(NAME recursion_limit_tests COMMAND recursion_limit_tests)

add_executable(return_tests test_return.c ${COMMON_SRC})
target_include_directories(return_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(return_tests gmp m ${ECM_LIB})
    add_dependencies(return_tests build_ecm)
else()
    target_link_libraries(return_tests gmp m)
endif()
add_test(NAME return_tests COMMAND return_tests)

add_executable(repl_hooks_tests test_repl_hooks.c ${COMMON_SRC})
target_include_directories(repl_hooks_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(repl_hooks_tests gmp m ${ECM_LIB})
    add_dependencies(repl_hooks_tests build_ecm)
else()
    target_link_libraries(repl_hooks_tests gmp m)
endif()
add_test(NAME repl_hooks_tests COMMAND repl_hooks_tests)

add_executable(extension_options_tests test_extension_options.c ${COMMON_SRC})
target_include_directories(extension_options_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(extension_options_tests gmp m ${ECM_LIB})
    add_dependencies(extension_options_tests build_ecm)
else()
    target_link_libraries(extension_options_tests gmp m)
endif()
add_test(NAME extension_options_tests COMMAND extension_options_tests)

add_executable(intrat_tests test_intrat.c ${COMMON_SRC})
target_include_directories(intrat_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(intrat_tests gmp m ${ECM_LIB})
    add_dependencies(intrat_tests build_ecm)
else()
    target_link_libraries(intrat_tests gmp m)
endif()
add_test(NAME intrat_tests COMMAND intrat_tests)

add_executable(integrate_unknown_tests test_integrate_unknown.c ${COMMON_SRC})
target_include_directories(integrate_unknown_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(integrate_unknown_tests gmp m ${ECM_LIB})
    add_dependencies(integrate_unknown_tests build_ecm)
else()
    target_link_libraries(integrate_unknown_tests gmp m)
endif()
add_test(NAME integrate_unknown_tests COMMAND integrate_unknown_tests)

add_executable(intrischnorman_tests test_intrischnorman.c ${COMMON_SRC})
target_include_directories(intrischnorman_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(intrischnorman_tests gmp m ${ECM_LIB})
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
    target_link_libraries(bench_intrischnorman gmp m ${ECM_LIB})
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
    target_link_libraries(intrat_corpus_tests gmp m ${ECM_LIB})
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
    target_link_libraries(crc_corpus_tests gmp m ${ECM_LIB})
    add_dependencies(crc_corpus_tests build_ecm)
else()
    target_link_libraries(crc_corpus_tests gmp m)
endif()
add_test(NAME crc_corpus_tests
         COMMAND crc_corpus_tests ../../CRCIntegralsCorpus.m)

add_executable(readwrite_tests test_readwrite.c ${COMMON_SRC})
target_include_directories(readwrite_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(readwrite_tests gmp m ${ECM_LIB})
    add_dependencies(readwrite_tests build_ecm)
else()
    target_link_libraries(readwrite_tests gmp m)
endif()
add_test(NAME readwrite_tests COMMAND readwrite_tests)

add_executable(files_tests test_files.c ${COMMON_SRC})
target_include_directories(files_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(files_tests gmp m ${ECM_LIB})
    add_dependencies(files_tests build_ecm)
else()
    target_link_libraries(files_tests gmp m)
endif()
add_test(NAME files_tests COMMAND files_tests)

add_executable(time_constrained_tests test_time_constrained.c ${COMMON_SRC})
target_include_directories(time_constrained_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(time_constrained_tests gmp m ${ECM_LIB})
    add_dependencies(time_constrained_tests build_ecm)
else()
    target_link_libraries(time_constrained_tests gmp m)
endif()
add_test(NAME time_constrained_tests COMMAND time_constrained_tests)

add_executable(powerexpand_tests test_powerexpand.c ${COMMON_SRC})
target_include_directories(powerexpand_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(powerexpand_tests gmp m ${ECM_LIB})
    add_dependencies(powerexpand_tests build_ecm)
else()
    target_link_libraries(powerexpand_tests gmp m)
endif()
add_test(NAME powerexpand_tests COMMAND powerexpand_tests)

add_executable(fullsimplify_tests test_fullsimplify.c ${COMMON_SRC})
target_include_directories(fullsimplify_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(fullsimplify_tests gmp m ${ECM_LIB})
    add_dependencies(fullsimplify_tests build_ecm)
else()
    target_link_libraries(fullsimplify_tests gmp m)
endif()
add_test(NAME fullsimplify_tests COMMAND fullsimplify_tests)

# Per-rule corpus: loads tests/fullsimplify_corpus.m and checks
# FullSimplify[input] == expected for each {input, expected} pair.
add_executable(fullsimplify_corpus_tests test_fullsimplify_corpus.c ${COMMON_SRC})
target_include_directories(fullsimplify_corpus_tests PRIVATE ../src)
if(USE_ECM)
    target_link_libraries(fullsimplify_corpus_tests gmp m ${ECM_LIB})
    add_dependencies(fullsimplify_corpus_tests build_ecm)
else()
    target_link_libraries(fullsimplify_corpus_tests gmp m)
endif()
add_test(NAME fullsimplify_corpus_tests
         COMMAND fullsimplify_corpus_tests ../fullsimplify_corpus.m)
