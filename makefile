USE_ECM ?= 1
USE_MPFR ?= 1

CC = gcc
CFLAGS = -O3 -std=c99 -Wall -Wextra -g -I./src -I./src/linalg -I./src/poly -I/usr/local/include
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

SRC_DIR = src
SRC = $(wildcard $(SRC_DIR)/*.c) $(wildcard $(SRC_DIR)/linalg/*.c) $(wildcard $(SRC_DIR)/poly/*.c)
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

.PHONY: all clean
