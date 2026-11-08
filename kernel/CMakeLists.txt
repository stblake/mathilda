#############################################################################
# Mathilda xeus Jupyter kernel (xmathilda)
#
# Builds a native C++17 Jupyter kernel that embeds Mathilda's evaluator in
# process via its C ABI (src/ffi/mathilda_ffi.h). No Python is needed at
# runtime — only xeus and the Mathilda static library.
#
# Quick start (with the xeus toolchain from conda-forge on PATH / CMAKE_PREFIX_PATH):
#   cmake -S kernel -B kernel/build -DCMAKE_BUILD_TYPE=Release
#   cmake --build kernel/build -j
#   cmake --install kernel/build --prefix <env>     # installs the kernelspec
#############################################################################
cmake_minimum_required(VERSION 3.16)
project(xmathilda VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# --- The Mathilda source tree (one level up by default) --------------------
get_filename_component(MATHILDA_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
set(MATHILDA_ROOT "${MATHILDA_ROOT}" CACHE PATH "Path to the Mathilda source tree")

# Extra `make` flags used when (re)building libmathilda.a, e.g. to drop
# optional dependencies for a leaner, headless kernel:
#   -DMATHILDA_MAKE_FLAGS="USE_GRAPHICS=0;USE_FFTW=0;USE_ECM=0"
set(MATHILDA_MAKE_FLAGS "" CACHE STRING "Extra flags passed to `make libmathilda.a`")

# --- xeus toolchain --------------------------------------------------------
# Both xeus 5.x and 6.x are supported; src/mathilda_interpreter.hpp picks the
# xinterpreter API branch from XEUS_VERSION_MAJOR (see MATHILDA_XEUS_6 there).
# The floor is a real one: xeus 4's execute_request_impl took an xrequest_context
# rather than a send_reply_callback, which is more than a signature away. Say so
# at configure time — a version error here is far easier to act on than a wall of
# "marked 'override', but does not override" (issue #85).
find_package(xeus 5.0 REQUIRED)
find_package(xeus-zmq 3.0 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(Threads REQUIRED)

if(xeus_VERSION VERSION_GREATER_EQUAL 6.0)
    message(STATUS "xmathilda: xeus ${xeus_VERSION} — xinterpreter API 6 "
                   "(shutdown_request_impl(bool), interrupt_request_impl)")
else()
    message(STATUS "xmathilda: xeus ${xeus_VERSION} — xinterpreter API 5 "
                   "(void shutdown_request_impl(), no interrupt hook)")
endif()

# --- Build (or locate) libmathilda.a ---------------------------------------
# The archive is the whole evaluator minus main() (repl.o), plus the FFI entry
# points. Built by the Mathilda makefile with its own GCC toolchain — which is
# independent of this CMake project's C++ compiler (they meet only at the C
# ABI, which is portable).
#
# The target ALWAYS runs `make`: the makefile, not CMake, knows libmathilda.a's
# real prerequisites (every .o and, via -MMD, every header), and is a no-op when
# the archive is current. An add_custom_command(OUTPUT …) with no DEPENDS would
# build the archive once and then never again — the kernel silently embedded a
# weeks-old evaluator while the REPL moved on. xmathilda links the archive by
# full path, so CMake relinks it whenever make rewrites the file.
set(MATHILDA_LIB "${MATHILDA_ROOT}/libmathilda.a")
add_custom_target(mathilda_lib
    COMMAND make -C "${MATHILDA_ROOT}" libmathilda.a ${MATHILDA_MAKE_FLAGS}
    BYPRODUCTS "${MATHILDA_LIB}"
    WORKING_DIRECTORY "${MATHILDA_ROOT}"
    COMMENT "make libmathilda.a (incremental)"
    VERBATIM)

# --- Transitive C dependencies of libmathilda.a ----------------------------
# The makefile, not CMake, knows them: it probes for SuiteSparse, ARPACK,
# LAPACKE and its BLAS provider (FlexiBLAS on Fedora/RHEL), FLINT, PCRE2, ...
# with the same USE_* flags the archive is built with. So ask it — `make
# print-libs` prints the archive's exact link line, the words the installed
# mathilda.pc records as Libs.private. Guessing the set here with find_library
# missed every library the makefile autodetected beyond the core few, and the
# link failed with undefined references to umfpack_*, cholmod_*, dsaupd_c,
# LAPACKE_* (issue #131, Fedora).
#
# Set MATHILDA_EXTRA_LIBS to pin an exact list instead (a ;-list of -l/-L
# flags, library names or full paths). Empty (the default) means "query the
# makefile on every configure", so changing MATHILDA_MAKE_FLAGS re-derives it.
set(MATHILDA_EXTRA_LIBS "" CACHE STRING
    "Transitive libs for libmathilda.a (empty: ask `make print-libs`)")

if(MATHILDA_EXTRA_LIBS)
    set(_mth_libs ${MATHILDA_EXTRA_LIBS})
    message(STATUS "xmathilda: MATHILDA_EXTRA_LIBS pinned: ${_mth_libs}")
else()
    execute_process(
        COMMAND make -s --no-print-directory -C "${MATHILDA_ROOT}"
                print-libs ${MATHILDA_MAKE_FLAGS}
        OUTPUT_VARIABLE _mth_line
        ERROR_QUIET
        RESULT_VARIABLE _mth_rc
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT _mth_rc EQUAL 0 OR _mth_line STREQUAL "")
        message(FATAL_ERROR
            "xmathilda: `make -C ${MATHILDA_ROOT} print-libs` failed (exit "
            "${_mth_rc}); cannot determine libmathilda.a's link libraries. "
            "Pass them explicitly with -DMATHILDA_EXTRA_LIBS=\"gmp;mpfr;...\".")
    endif()
    # Split the shell words, then rejoin each two-word `-framework X` into one
    # item: CMake would otherwise read a bare `X` as a library name (-lX).
    separate_arguments(_mth_words UNIX_COMMAND "${_mth_line}")
    set(_mth_libs "")
    set(_mth_pending "")
    foreach(_w IN LISTS _mth_words)
        if(_mth_pending)
            list(APPEND _mth_libs "-framework ${_w}")
            set(_mth_pending "")
        elseif(_w STREQUAL "-framework")
            set(_mth_pending 1)
        else()
            list(APPEND _mth_libs "${_w}")
        endif()
    endforeach()
    message(STATUS "xmathilda: libmathilda.a links: ${_mth_line}")
endif()

# --- The kernel executable -------------------------------------------------
add_executable(xmathilda
    src/main.cpp
    src/mathilda_interpreter.cpp)
add_dependencies(xmathilda mathilda_lib)

target_include_directories(xmathilda PRIVATE
    "${MATHILDA_ROOT}/src"
    "${MATHILDA_ROOT}/src/ffi")

target_link_libraries(xmathilda PRIVATE
    xeus xeus-zmq nlohmann_json::nlohmann_json
    "${MATHILDA_LIB}"
    ${_mth_libs}
    Threads::Threads
    m)

# libxeus and libxeus-zmq are shared and their install_names are @rpath-relative,
# so the INSTALLED binary needs an rpath naming the toolchain's lib directory.
# CMake strips the build-tree rpath on install and adds nothing by default, so
# without this the kernelspec launches a binary that dies at dyld
# ("Library not loaded: @rpath/libxeus-zmq.N.dylib ... no LC_RPATH's found") and
# Jupyter reports only "Kernel died before replying to kernel_info". It happened
# to work when the build ran inside an activated conda env, which injects
# -Wl,-rpath via LDFLAGS — i.e. it worked by accident, and not for a packager.
# Distro builds that manage rpaths themselves can override with
# -DCMAKE_SKIP_INSTALL_RPATH=ON.
set_target_properties(xmathilda PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

# --- Kernelspec installation ----------------------------------------------
# kernel.json tells Jupyter how to launch the kernel; MATHILDA_HOME points the
# loader at the internal/ module tree (init.m etc.) shipped in the source.
#
# The absolute path in "argv" must be where the binary actually LANDS, so
# kernel.json has to be generated at install time, not configure time: the
# documented flow is `cmake --install … --prefix <env>`, and --prefix overrides
# CMAKE_INSTALL_PREFIX only for that install step. Baking it in at configure time
# silently wrote /usr/local/bin/xmathilda (the default prefix) into the
# kernelspec, so Jupyter launched the wrong binary or none at all. install(CODE)
# runs with the effective prefix; install rules fire in declaration order, so the
# regenerated file is the one install(FILES) picks up.
# $MATHILDA_HOME names the internal/ directory ITSELF (the one holding init.m),
# not its parent: pointing it at src/ loaded init.m only by luck of the CWD and
# missed the integral tables (issue #96). The loader now accepts either.
set(XEUS_MATHILDA_HOME "${MATHILDA_ROOT}/src/internal" CACHE PATH "The internal/ module dir (holds init.m)")

install(TARGETS xmathilda RUNTIME DESTINATION bin)
install(CODE "
    set(XKERNEL_BIN \"\${CMAKE_INSTALL_PREFIX}/bin/xmathilda\")
    set(XEUS_MATHILDA_HOME \"${XEUS_MATHILDA_HOME}\")
    configure_file(\"${CMAKE_CURRENT_SOURCE_DIR}/kernel.json.in\"
                   \"${CMAKE_CURRENT_BINARY_DIR}/kernel.json\" @ONLY)
    message(STATUS \"xmathilda kernelspec argv[0]: \${XKERNEL_BIN}\")
")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/kernel.json"
        DESTINATION share/jupyter/kernels/xmathilda)
