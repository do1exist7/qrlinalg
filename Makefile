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
  FIXED_FLAGS := -ffixed-line-length-none
  PREPROCESS_FLAGS := -cpp
  WP_FLAGS := $(PREPROCESS_FLAGS) -DQRLINALG_WP=$(PREC)
  RELEASE_FLAGS := -O3 -march=native
  DEBUG_FLAGS := -O0 -g -fcheck=all -fbacktrace
  OWN_WARNING_FLAGS := -Wall -Wextra -Wno-unused-dummy-argument
else ifeq ($(COMPILER),ifort)
  FC := ifort
  MODULE_FLAGS := -I$(BUILD_DIR) -module $(BUILD_DIR)
  FIXED_FLAGS := -extend-source
  PREPROCESS_FLAGS := -fpp
  WP_FLAGS := $(PREPROCESS_FLAGS) -DQRLINALG_WP=$(PREC)
  RELEASE_FLAGS := -O3 -ip -fp-model precise
  DEBUG_FLAGS := -O0 -g -check all -traceback
  OWN_WARNING_FLAGS := -warn all
else ifeq ($(COMPILER),ifx)
  FC := ifx
  MODULE_FLAGS := -I$(BUILD_DIR) -module $(BUILD_DIR)
  FIXED_FLAGS := -extend-source
  PREPROCESS_FLAGS := -fpp
  WP_FLAGS := $(PREPROCESS_FLAGS) -DQRLINALG_WP=$(PREC)
  RELEASE_FLAGS := -O3 -fp-model precise
  DEBUG_FLAGS := -O0 -g -check all -traceback
  OWN_WARNING_FLAGS := -warn all
else ifeq ($(COMPILER),nvfortran)
  FC := nvfortran
  MODULE_FLAGS := -I$(BUILD_DIR) -module $(BUILD_DIR)
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
TEST_EXE := $(BUILD_DIR)/test_qrlinalg

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

.PHONY: all release debug build check test test-one clean

all: release

release:
	$(MAKE) CONFIG=release build

debug:
	$(MAKE) CONFIG=debug build

build: $(LIB)

check: test

test:
	$(MAKE) CONFIG=debug PREC=8 BUILD_DIR=build/test-wp8 QRLINALG_TEST_FLAGS=-DQRLINALG_TESTING test-one
	$(MAKE) CONFIG=debug PREC=10 BUILD_DIR=build/test-wp10 QRLINALG_TEST_FLAGS=-DQRLINALG_TESTING test-one
	$(MAKE) CONFIG=debug PREC=16 BUILD_DIR=build/test-wp16 QRLINALG_TEST_FLAGS=-DQRLINALG_TESTING test-one

test-one: $(TEST_EXE)
	$(TEST_EXE)

$(BUILD_DIR):
	mkdir -p $@

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

$(BUILD_DIR)/wp_def.o: $(SRC_DIR)/wp_def_$(PREC).f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(WP_FLAGS) -c $< -o $@

$(BUILD_DIR)/qrupdate_blas.o: $(QRUPDATE_DIR)/BLAS.f $(BUILD_DIR)/wp_def.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(FIXED_FLAGS) -c $< -o $@

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

$(BUILD_DIR)/qrlinalg.o: $(SRC_DIR)/qrlinalg.f90 $(BUILD_DIR)/wp_def.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(PREPROCESS_FLAGS) $(QRLINALG_TEST_FLAGS) \
		$(OWN_WARNING_FLAGS) -c $< -o $@

$(TEST_EXE): test/test_qrlinalg.f90 $(LIB)
	$(FC) $(FFLAGS) $(OWN_WARNING_FLAGS) $< $(LIB) -o $@

clean:
	rm -rf build
