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
find_package(xeus REQUIRED)
find_package(xeus-zmq REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(Threads REQUIRED)

# --- Build (or locate) libmathilda.a ---------------------------------------
# The archive is the whole evaluator minus main() (repl.o), plus the FFI entry
# points. Built by the Mathilda makefile with its own GCC toolchain — which is
# independent of this CMake project's C++ compiler (they meet only at the C
# ABI, which is portable).
set(MATHILDA_LIB "${MATHILDA_ROOT}/libmathilda.a")
add_custom_command(
    OUTPUT "${MATHILDA_LIB}"
    COMMAND ${CMAKE_COMMAND} -E echo "Building libmathilda.a via the Mathilda makefile..."
    COMMAND make -C "${MATHILDA_ROOT}" libmathilda.a ${MATHILDA_MAKE_FLAGS}
    WORKING_DIRECTORY "${MATHILDA_ROOT}"
    COMMENT "make libmathilda.a"
    VERBATIM)
add_custom_target(mathilda_lib DEPENDS "${MATHILDA_LIB}")

# --- Transitive C dependencies of libmathilda.a ----------------------------
# Which of these are actually referenced depends on the USE_* flags the archive
# was built with. find_library is tolerant: a missing optional lib is simply
# not linked. Override MATHILDA_EXTRA_LIBS to pin an exact set.
set(_mth_libs "")
foreach(_name gmp mpfr flint pcre2-8 fftw3 ecm readline raylib)
    find_library(_lib_${_name} NAMES ${_name})
    if(_lib_${_name})
        list(APPEND _mth_libs ${_lib_${_name}})
    endif()
endforeach()
if(APPLE)
    # LAPACK/BLAS via Accelerate; CoreGraphics/CoreFoundation for the raylib
    # offscreen path. On Linux use -DMATHILDA_EXTRA_LIBS="lapacke;lapack;blas;gfortran".
    list(APPEND _mth_libs "-framework Accelerate"
                          "-framework CoreGraphics"
                          "-framework CoreFoundation")
endif()
set(MATHILDA_EXTRA_LIBS "${_mth_libs}" CACHE STRING "Transitive libs for libmathilda.a")

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
    ${MATHILDA_EXTRA_LIBS}
    Threads::Threads
    m)

# --- Kernelspec installation ----------------------------------------------
# kernel.json tells Jupyter how to launch the kernel; MATHILDA_HOME points the
# loader at the internal/ module tree (init.m etc.) shipped in the source.
set(XEUS_MATHILDA_HOME "${MATHILDA_ROOT}/src" CACHE PATH "Dir containing internal/")
set(XKERNEL_BIN "${CMAKE_INSTALL_PREFIX}/bin/xmathilda")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/kernel.json.in"
               "${CMAKE_CURRENT_BINARY_DIR}/kernel.json" @ONLY)

install(TARGETS xmathilda RUNTIME DESTINATION bin)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/kernel.json"
        DESTINATION share/jupyter/kernels/xmathilda)
