COMPILER ?= gfortran
PREC ?= 8
CONFIG ?= release

SRC_DIR := src
QRUPDATE_DIR := $(SRC_DIR)/qrupdate
BUILD_DIR := build/$(CONFIG)-wp$(PREC)
LIB := $(BUILD_DIR)/libqrlinalg.a

ifneq ($(PREC),$(filter $(PREC),8 10 16))
  $(error PREC must be 8, 10, or 16)
endif

ifeq ($(COMPILER),gfortran)
  FC := gfortran
  MODULE_FLAGS := -I$(BUILD_DIR) -J$(BUILD_DIR)
  ORIG_MODULE_FLAGS = -I$(DIFF_ORIG_DIR) -J$(DIFF_ORIG_DIR)
  FIXED_FLAGS := -ffixed-line-length-none
  PREPROCESS_FLAGS := -cpp
  WP_FLAGS := $(PREPROCESS_FLAGS) -DQRLINALG_WP=$(PREC)
  RELEASE_FLAGS := -O3 -march=native
  DEBUG_FLAGS := -O0 -g -fcheck=all -fbacktrace
  OWN_WARNING_FLAGS := -Wall -Wextra -Wno-unused-dummy-argument
else ifeq ($(COMPILER),ifort)
  FC := ifort
  MODULE_FLAGS := -I$(BUILD_DIR) -module $(BUILD_DIR)
  ORIG_MODULE_FLAGS = -I$(DIFF_ORIG_DIR) -module $(DIFF_ORIG_DIR)
  FIXED_FLAGS := -extend-source
  PREPROCESS_FLAGS := -fpp
  WP_FLAGS := $(PREPROCESS_FLAGS) -DQRLINALG_WP=$(PREC)
  RELEASE_FLAGS := -O3 -ip -fp-model precise
  DEBUG_FLAGS := -O0 -g -check all -traceback
  OWN_WARNING_FLAGS := -warn all
else ifeq ($(COMPILER),ifx)
  FC := ifx
  MODULE_FLAGS := -I$(BUILD_DIR) -module $(BUILD_DIR)
  ORIG_MODULE_FLAGS = -I$(DIFF_ORIG_DIR) -module $(DIFF_ORIG_DIR)
  FIXED_FLAGS := -extend-source
  PREPROCESS_FLAGS := -fpp
  WP_FLAGS := $(PREPROCESS_FLAGS) -DQRLINALG_WP=$(PREC)
  RELEASE_FLAGS := -O3 -fp-model precise
  DEBUG_FLAGS := -O0 -g -check all -traceback
  OWN_WARNING_FLAGS := -warn all
else ifeq ($(COMPILER),nvfortran)
  FC := nvfortran
  MODULE_FLAGS := -I$(BUILD_DIR) -module $(BUILD_DIR)
  ORIG_MODULE_FLAGS = -I$(DIFF_ORIG_DIR) -module $(DIFF_ORIG_DIR)
  FIXED_FLAGS := -Mextend
  PREPROCESS_FLAGS := -Mpreprocess
  WP_FLAGS := $(PREPROCESS_FLAGS) -DQRLINALG_WP=$(PREC)
  RELEASE_FLAGS := -O3 -tp=native
  DEBUG_FLAGS := -O0 -g -Mbounds -traceback
  OWN_WARNING_FLAGS := -Minform=warn
else
  $(error unsupported COMPILER=$(COMPILER))
endif

ifeq ($(CONFIG),release)
  CONFIG_FLAGS := $(RELEASE_FLAGS)
else ifeq ($(CONFIG),debug)
  CONFIG_FLAGS := $(DEBUG_FLAGS)
else
  $(error CONFIG must be release or debug)
endif

FFLAGS := $(MODULE_FLAGS) $(CONFIG_FLAGS)
ARFLAGS = rcs
QRLINALG_TEST_FLAGS ?=
TEST_NAMES := test_initialization test_factorization test_dgemm test_replacement \
	test_append test_delete test_inverse_iteration \
	test_inverse_iteration_failures
TEST_EXES := $(addprefix $(BUILD_DIR)/,$(TEST_NAMES))
TEST_SUPPORT_OBJECT := $(BUILD_DIR)/test_support.o

PYTHON ?= python3
DATA_DIR ?= data
DIFF_MANIFEST ?= $(DATA_DIR)/cases.csv
DIFF_DIR := build/differential-$(CONFIG)-wp$(PREC)
DIFF_NEW_DIR := $(DIFF_DIR)/new
DIFF_ORIG_DIR := $(DIFF_DIR)/orig
DIFF_NEW_EXE := $(DIFF_NEW_DIR)/qrlinalg_driver
DIFF_ORIG_EXE := $(DIFF_ORIG_DIR)/orig_driver
DIFF_OUTPUT_DIR := $(DIFF_DIR)/results
DIFF_SOURCE_DIR := test/differential
ORIG_SOURCE_DIR := orig/claude
ORIG_MPI_OBJECT := $(DIFF_ORIG_DIR)/mpi_serial.o
ORIG_WP_OBJECT := $(DIFF_ORIG_DIR)/wp_def.o
ORIG_GLOBVARS_OBJECT := $(DIFF_ORIG_DIR)/globvars.o
ORIG_BLAS_OBJECT := $(DIFF_ORIG_DIR)/blas.o
ORIG_LINALG_OBJECT := $(DIFF_ORIG_DIR)/linalg.o
ORIG_DRIVER_OBJECT := $(DIFF_ORIG_DIR)/orig_driver.o
ORIG_OBJECTS := $(ORIG_MPI_OBJECT) $(ORIG_WP_OBJECT) $(ORIG_GLOBVARS_OBJECT) \
	$(ORIG_BLAS_OBJECT) $(ORIG_LINALG_OBJECT) $(ORIG_DRIVER_OBJECT)

OBJECTS := \
	$(BUILD_DIR)/wp_def.o \
	$(BUILD_DIR)/qrupdate_blas.o \
	$(BUILD_DIR)/qrupdate_lapack.o \
	$(BUILD_DIR)/qrupdate_linalg.o \
	$(BUILD_DIR)/qrupdate_error.o \
	$(BUILD_DIR)/qrupdate_real.o \
	$(BUILD_DIR)/qrupdate_complex.o \
	$(BUILD_DIR)/qrupdate.o \
	$(BUILD_DIR)/qrlinalg.o

.PHONY: all release debug build check check-one test test-one differential-build \
	differential-test compare-orig clean

all: release

release:
	$(MAKE) CONFIG=release build

debug:
	$(MAKE) CONFIG=debug build

build: $(LIB)

check: test

check-one:
	$(MAKE) CONFIG=debug PREC=$(PREC) BUILD_DIR=build/test-wp$(PREC) \
		QRLINALG_TEST_FLAGS=-DQRLINALG_TESTING test-one

test:
	$(MAKE) check-one PREC=8
	$(MAKE) check-one PREC=10
	$(MAKE) check-one PREC=16

test-one: $(TEST_EXES)
	@set -e; for test_exe in $(TEST_EXES); do $$test_exe; done

differential-build: $(DIFF_NEW_EXE) $(DIFF_ORIG_EXE)

differential-test compare-orig: differential-build
	$(PYTHON) $(DIFF_SOURCE_DIR)/compare_outputs.py \
		--manifest "$(DIFF_MANIFEST)" \
		--new-exe "$(DIFF_NEW_EXE)" \
		--orig-exe "$(DIFF_ORIG_EXE)" \
		--output-dir "$(DIFF_OUTPUT_DIR)"

$(BUILD_DIR):
	mkdir -p $@

$(DIFF_NEW_DIR) $(DIFF_ORIG_DIR):
	mkdir -p $@

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

$(BUILD_DIR)/wp_def.o: $(SRC_DIR)/wp_def_$(PREC).f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(WP_FLAGS) -c $< -o $@

$(BUILD_DIR)/qrupdate_blas.o: $(QRUPDATE_DIR)/BLAS.f $(BUILD_DIR)/wp_def.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(FIXED_FLAGS) $(WP_FLAGS) -c $< -o $@

$(BUILD_DIR)/qrupdate_lapack.o: $(QRUPDATE_DIR)/LAPACK.f $(BUILD_DIR)/wp_def.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(FIXED_FLAGS) -c $< -o $@

$(BUILD_DIR)/qrupdate_linalg.o: $(QRUPDATE_DIR)/qrupdate_linalg.f90 $(BUILD_DIR)/wp_def.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/qrupdate_error.o: $(QRUPDATE_DIR)/qrupdate_error.f90 $(BUILD_DIR)/qrupdate_linalg.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/qrupdate_real.o: $(QRUPDATE_DIR)/qrupdate_real.f90 $(BUILD_DIR)/qrupdate_error.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/qrupdate_complex.o: $(QRUPDATE_DIR)/qrupdate_complex.f90 $(BUILD_DIR)/qrupdate_error.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/qrupdate.o: $(QRUPDATE_DIR)/qrupdate.f90 $(BUILD_DIR)/qrupdate_real.o $(BUILD_DIR)/qrupdate_complex.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/qrlinalg.o: $(SRC_DIR)/qrlinalg.f90 $(BUILD_DIR)/wp_def.o \
		$(BUILD_DIR)/qrupdate.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(PREPROCESS_FLAGS) $(QRLINALG_TEST_FLAGS) \
		$(OWN_WARNING_FLAGS) -c $< -o $@

$(TEST_SUPPORT_OBJECT): test/test_support.f90 $(LIB)
	$(FC) $(FFLAGS) $(OWN_WARNING_FLAGS) -c $< -o $@

$(BUILD_DIR)/test_%: test/test_%.f90 $(TEST_SUPPORT_OBJECT) $(LIB)
	$(FC) $(FFLAGS) $(OWN_WARNING_FLAGS) $< $(TEST_SUPPORT_OBJECT) $(LIB) -o $@

$(DIFF_NEW_EXE): $(DIFF_SOURCE_DIR)/qrlinalg_driver.f90 $(LIB) | $(DIFF_NEW_DIR)
	$(FC) $(FFLAGS) $(OWN_WARNING_FLAGS) $< $(LIB) -o $@

$(ORIG_MPI_OBJECT): $(DIFF_SOURCE_DIR)/mpi_serial.f90 | $(DIFF_ORIG_DIR)
	$(FC) $(ORIG_MODULE_FLAGS) $(CONFIG_FLAGS) -c $< -o $@

$(ORIG_WP_OBJECT): $(ORIG_SOURCE_DIR)/wp_def_$(PREC).f90 $(ORIG_MPI_OBJECT) | $(DIFF_ORIG_DIR)
	$(FC) $(ORIG_MODULE_FLAGS) $(CONFIG_FLAGS) -c $< -o $@

$(ORIG_GLOBVARS_OBJECT): $(ORIG_SOURCE_DIR)/globvars.f90 $(ORIG_WP_OBJECT) | $(DIFF_ORIG_DIR)
	$(FC) $(ORIG_MODULE_FLAGS) $(CONFIG_FLAGS) -c $< -o $@

$(ORIG_BLAS_OBJECT): $(QRUPDATE_DIR)/BLAS.f $(ORIG_WP_OBJECT) | $(DIFF_ORIG_DIR)
	$(FC) $(ORIG_MODULE_FLAGS) $(CONFIG_FLAGS) $(FIXED_FLAGS) $(WP_FLAGS) \
		-c $< -o $@

$(ORIG_LINALG_OBJECT): $(ORIG_SOURCE_DIR)/linalg.f90 $(ORIG_GLOBVARS_OBJECT) | $(DIFF_ORIG_DIR)
	$(FC) $(ORIG_MODULE_FLAGS) $(CONFIG_FLAGS) -c $< -o $@

$(ORIG_DRIVER_OBJECT): $(DIFF_SOURCE_DIR)/orig_driver.f90 $(ORIG_LINALG_OBJECT) | $(DIFF_ORIG_DIR)
	$(FC) $(ORIG_MODULE_FLAGS) $(CONFIG_FLAGS) -c $< -o $@

$(DIFF_ORIG_EXE): $(ORIG_OBJECTS)
	$(FC) $(ORIG_MODULE_FLAGS) $(CONFIG_FLAGS) $^ -o $@

clean:
	rm -rf build
