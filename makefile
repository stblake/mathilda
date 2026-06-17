USE_ECM ?= 1
USE_MPFR ?= 1
USE_LAPACK ?= 1

CC = gcc
CFLAGS = -O3 -std=c99 -Wall -Wextra -g -I./src -I./src/linalg -I./src/poly -I./src/simp -I./src/calculus -I./src/sum -I./src/special_functions -I./src/numerical_calculus -I./src/numerical_roots -I/usr/local/include
LDFLAGS = -lreadline -L/usr/local/lib -lgmp -lm

ifeq ($(USE_ECM), 1)
CFLAGS += -I./src/external/ecm
LDFLAGS := src/external/ecm/.libs/libecm.a $(LDFLAGS)
ECM_TARGET = src/external/ecm/.libs/libecm.a
else
CFLAGS += -DNO_ECM
ECM_TARGET =
endif

# Arbitrary-precision reals (MPFR) — enables N[expr, prec], Precision/
# Accuracy/SetPrecision/SetAccuracy, and `3.98`50` precision literals.
# Disable with `USE_MPFR=0` to build without the MPFR dependency; the
# machine-precision path continues to work with a runtime warning on the
# unsupported operations.
ifeq ($(USE_MPFR), 1)
CFLAGS  += -DUSE_MPFR
LDFLAGS += -lmpfr
endif

ifneq ($(wildcard /opt/homebrew/include),)
CFLAGS += -I/opt/homebrew/include
endif

ifneq ($(wildcard /opt/homebrew/lib),)
LDFLAGS += -L/opt/homebrew/lib
endif

# BLAS/LAPACK for fast machine-precision linear-algebra kernels
# (machine-precision QRDecomposition; later: Inverse, LinearSolve, Det,
# Eigenvalues, LeastSquares, SVD).  Four-tier autodetection:
#   1. Darwin              -> Apple Accelerate framework (zero install).
#   2. pkg-config lapacke  -> use pkg-config flags (OpenBLAS / MKL / etc).
#   3. /usr/include/lapacke.h or /usr/local/include/lapacke.h
#                          -> link -llapacke -llapack -lblas.
#   4. nothing found       -> warn and override USE_LAPACK := 0.
# When USE_LAPACK is off, machine-precision linalg falls back to the
# MPFR / symbolic path with a one-time runtime warning.  This matches
# the existing USE_MPFR=0 / USE_ECM=0 graceful-degrade policy so that
# `git clone && make` always succeeds, no matter the host environment.
ifeq ($(USE_LAPACK), 1)
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S),Darwin)
    CFLAGS  += -DUSE_LAPACK -DMATHILDA_USE_ACCELERATE
    LDFLAGS += -framework Accelerate
  else ifneq ($(shell pkg-config --exists lapacke 2>/dev/null && echo y),)
    CFLAGS  += -DUSE_LAPACK $(shell pkg-config --cflags lapacke)
    LDFLAGS += $(shell pkg-config --libs lapacke)
  else ifneq ($(wildcard /usr/include/lapacke.h)$(wildcard /usr/local/include/lapacke.h),)
    CFLAGS  += -DUSE_LAPACK
    LDFLAGS += -llapacke -llapack -lblas
  else
    $(warning LAPACK/LAPACKE not detected; building with USE_LAPACK=0)
    $(warning   Ubuntu/Debian:  sudo apt install liblapacke-dev libopenblas-dev)
    $(warning   Fedora/RHEL:    sudo dnf install lapack-devel lapacke-devel openblas-devel)
    $(warning   Arch:           sudo pacman -S openblas lapacke)
    override USE_LAPACK := 0
  endif
endif

SRC_DIR = src
SRC = $(wildcard $(SRC_DIR)/*.c) $(wildcard $(SRC_DIR)/linalg/*.c) $(wildcard $(SRC_DIR)/poly/*.c) $(wildcard $(SRC_DIR)/simp/*.c) $(wildcard $(SRC_DIR)/calculus/*.c) $(wildcard $(SRC_DIR)/sum/*.c) $(wildcard $(SRC_DIR)/special_functions/*.c) $(wildcard $(SRC_DIR)/numerical_calculus/*.c) $(wildcard $(SRC_DIR)/numerical_roots/*.c)
OBJ = $(SRC:.c=.o)
TARGET = Mathilda

TEST_BINARIES = eval_tests expr_tests parse_tests test_ld test_ops test_pattern list_tests stats_tests expand_tests
TEST_DIR = tests
CMAKE_TEST_BINARIES = comparisons_tests eval_tests expr_tests match_tests match_extensive_tests parse_tests regression_tests symtab_tests list_tests trig_tests hyperbolic_tests logexp_tests piecewise_tests purefunc_tests stats_tests expand_tests

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

src/external/ecm/.libs/libecm.a:
	./build_ecm.sh

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(ECM_TARGET)
	$(CC) $(CFLAGS) -c $< -o $@

$(SRC_DIR)/boolean.o: $(SRC_DIR)/boolean.c
	$(CC) $(CFLAGS) -c $< -o $@

$(SRC_DIR)/list.o: $(SRC_DIR)/list.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
	rm -rf *.dSYM
	rm -f $(TEST_BINARIES)
	rm -f *~
	rm -f $(SRC_DIR)/*~
	rm -f $(TEST_DIR)/*~
	rm -f *.o
	if [ -f $(TEST_DIR)/Makefile ]; then $(MAKE) -C $(TEST_DIR) clean; fi
	rm -f $(addprefix $(TEST_DIR)/, $(CMAKE_TEST_BINARIES))

# Regenerate the documentation website's per-builtin pages from the docstrings,
# attributes, and spec examples. Requires the built ./Mathilda binary (examples
# are verified against it). The generated Markdown is committed; CI only builds
# the MkDocs site from it.
docs: $(TARGET)
	python3 site/generate.py

# Build the static site locally (needs `pip install -r site/requirements.txt`).
docs-build:
	mkdocs build --strict -f site/mkdocs.yml

# Serve the site with live reload at http://127.0.0.1:8000
docs-serve:
	mkdocs serve -f site/mkdocs.yml

.PHONY: all clean docs docs-build docs-serve
