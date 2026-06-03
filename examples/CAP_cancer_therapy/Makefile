THIS_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
ABM4BIO_ROOT := $(realpath $(THIS_DIR)/../..)
ABM_BINARY := $(ABM4BIO_ROOT)/build/ABM4bio
BDM_ENV_SCRIPT := $(ABM4BIO_ROOT)/libs/biodynamo-v1.05.143/bin/thisbdm.sh
EXCEL_DATA := $(ABM4BIO_ROOT)/libs/experimental_data/Data Modeling CCA+PDAC 0 to 72h post Gorjet - WITH MEAN 4 experiments.xlsx
# Use miniconda Python explicitly (has scikit-optimize installed)
PYTHON := /home/aiwsif/miniconda/bin/python3

# Cap BLAS/OpenMP threads for ALL python invocations so numpy import
# never tries to fork 64 threads and exhausts RLIMIT_NPROC.
export OPENBLAS_NUM_THREADS := 1
export MKL_NUM_THREADS := 1
export NUMEXPR_NUM_THREADS := 1
export VECLIB_MAXIMUM_THREADS := 1
export BLIS_NUM_THREADS := 1
export OMP_NUM_THREADS := 1
SCRIPTS := scripts

INPUT_CONTROL ?= input_control.csv
RESULTS_CONTROL ?= results_control
SEED ?= 1234
ABM_INITIAL_CELLS ?= 100
GRID_RESOLUTION ?= 10
VIZ_INTERVAL_H ?= 1
ABM_REPLICATES ?= 1
ABM_REPLICATE_SEED_STEP ?= 1
OMP_THREADS ?= 2
ABM_DA_MAXITER ?= 50
ABM_NM_MAXITER ?= 200
ABM_METHODS ?= bayesian powell
ABM_FIT_PARAMS ?= cancer_cell/can_divide/probability cancer_cell/can_apoptose/probability cancer_cell/can_grow/diameter_rate cancer_cell/can_grow/probability
ABM_OVERGROWTH_FACTOR ?= 4.0
ABM_REPLICATE_TIMEOUT_S ?= 900
ABM_MAX_EVALS ?= 200
ABM_BO_CALLS ?= 80
ABM_BO_INIT ?= 20
ABM_BAND_STD_SCALE ?= 2.0
ABM_OUT_DIR ?= calibration_outputs/EGI1_control_0to72h
CELL_LINE ?= EGI1

.PHONY: help build show-config prepare_input run_control check_trajectory calibrate_control calibrate_control_clean clean clean_calibration clean_all

help:
	@echo ""
	@echo "=== CAP in-vitro: untreated control calibration (0–72 h) ==="
	@echo ""
	@echo "  make build              Rebuild ABM4bio"
	@echo "  make prepare_input      Generate $(INPUT_CONTROL) (coarse grid, CAP off)"
	@echo "  make run_control        Single 0–72 h untreated simulation"
	@echo "  make check_trajectory   Print viable counts vs t0 from stats.csv"
	@echo "  make calibrate_control  Fit untreated control vs Excel (clean wrapper)"
	@echo ""
	@echo "Options (examples):"
	@echo "  make prepare_input GRID_RESOLUTION=10 VIZ_INTERVAL_H=1 ABM_INITIAL_CELLS=100"
	@echo "  make run_control OMP_THREADS=8 SEED=1234"
	@echo "  make calibrate_control ABM_BO_CALLS=100 ABM_BO_INIT=30 OMP_THREADS=2"
	@echo "  make calibrate_control ABM_METHODS='bayesian powell' ABM_BO_CALLS=80 ABM_BO_INIT=20"
	@echo ""
	@echo "Calibration pipeline (default: bayesian -> powell):"
	@echo "  ABM_METHODS='$(ABM_METHODS)'"
	@echo "    bayesian  - GP surrogate BO (global, noise-aware); use for exploration"
	@echo "    powell    - Powell refinement from BO best (local); use for polishing"
	@echo "    dual_annealing - optional global search (add if BO budget is large)"
	@echo "    NOTE: differential_evolution needs popsize*n_params*maxiter evals - too"
	@echo "          expensive for small budgets; excluded from default pipeline."
	@echo "  ABM_FIT_PARAMS='$(ABM_FIT_PARAMS)'"
	@echo "    Default 4-parameter set: divide_p, apoptose_p, grow_rate, grow_p"
	@echo "  OMP_THREADS=$(OMP_THREADS)  - Threads per ABM run (keep low to avoid fork exhaustion)"
	@echo "  ABM_OVERGROWTH_FACTOR=$(ABM_OVERGROWTH_FACTOR)  - Kill at N0 x max_viability x factor cells"
	@echo "  ABM_REPLICATE_TIMEOUT_S=$(ABM_REPLICATE_TIMEOUT_S)  - Per-replicate wall-time kill"
	@echo "  ABM_BO_CALLS=$(ABM_BO_CALLS)  - Total BO function evaluations"
	@echo "  ABM_BO_INIT=$(ABM_BO_INIT)     - Random BO initialization points"
	@echo ""
	@echo "  make clean              Remove $(RESULTS_CONTROL)"
	@echo "  make clean_calibration  Remove calibration_outputs/"
	@echo "  make show-config        Print active variables"
	@echo ""

show-config:
	@echo "ABM_BINARY=$(ABM_BINARY)"
	@echo "EXCEL_DATA=$(EXCEL_DATA)"
	@echo "INPUT_CONTROL=$(INPUT_CONTROL)"
	@echo "RESULTS_CONTROL=$(RESULTS_CONTROL)"
	@echo "OMP_THREADS=$(OMP_THREADS)"
	@echo "GRID_RESOLUTION=$(GRID_RESOLUTION)"
	@echo "ABM_INITIAL_CELLS=$(ABM_INITIAL_CELLS)"
	@echo "ABM_FIT_PARAMS=$(ABM_FIT_PARAMS)"
	@echo "ABM_METHODS=$(ABM_METHODS)"
	@echo "ABM_BO_CALLS=$(ABM_BO_CALLS)"
	@echo "ABM_BO_INIT=$(ABM_BO_INIT)"
	@echo "ABM_OVERGROWTH_FACTOR=$(ABM_OVERGROWTH_FACTOR)"
	@echo "ABM_REPLICATE_TIMEOUT_S=$(ABM_REPLICATE_TIMEOUT_S)"
	@echo "ABM_OUT_DIR=$(ABM_OUT_DIR)"

build:
	cd $(ABM4BIO_ROOT) && make fresh -j10

prepare_input:
	cd $(THIS_DIR) && OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 $(PYTHON) $(SCRIPTS)/prepare_control_input.py \
		--base-input input_control_base.csv \
		--output $(INPUT_CONTROL) \
		--initial-cells $(ABM_INITIAL_CELLS) \
		--grid-resolution $(GRID_RESOLUTION) \
		--viz-interval-h $(VIZ_INTERVAL_H)

define REMOVE_RESULTS_DIR
	@dir="$(1)"; \
	if [ -z "$$dir" ] || [ "$$dir" = "/" ]; then echo "refusing unsafe path" >&2; exit 1; fi; \
	if [ ! -e "$$dir" ]; then exit 0; fi; \
	echo "Removing $$dir ..."; \
	chmod -R u+w "$$dir" 2>/dev/null || true; \
	find "$$dir" -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true; \
	rm -rf "$$dir" 2>/dev/null || true; \
	for attempt in 1 2 3 4 5; do \
		[ ! -e "$$dir" ] && exit 0; sleep 1; chmod -R u+w "$$dir" 2>/dev/null || true; rm -rf "$$dir" 2>/dev/null || true; \
	done; \
	[ ! -e "$$dir" ] || (echo "ERROR: could not remove $$dir" >&2; exit 1)
endef

run_control: prepare_input
	$(call REMOVE_RESULTS_DIR,$(RESULTS_CONTROL))
	bash -c 'export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 OMP_NUM_THREADS=$(OMP_THREADS); \
		export BDM_THISBDM_SILENT=true; source "$(BDM_ENV_SCRIPT)"; \
		"$(ABM_BINARY)" $(INPUT_CONTROL) $(SEED)'

check_trajectory:
	cd $(THIS_DIR) && $(PYTHON) $(SCRIPTS)/check_trajectory.py $(RESULTS_CONTROL)/stats.csv

calibrate_control: prepare_input
	$(MAKE) --no-print-directory calibrate_control_clean

calibrate_control_clean:
	@if [ "$(SEED)" -le 0 ]; then echo "SEED must be > 0"; exit 2; fi
	@echo "=== Checking Python environment ==="
	@echo "Using Python: $(PYTHON)"
	@$(PYTHON) -c "import sys; print('Python:', sys.executable)" && $(PYTHON) -c "import skopt; print('scikit-optimize:', skopt.__version__)"
	cd $(THIS_DIR) && PYTHON_BIN=$(PYTHON) BDM_ENV_SCRIPT=$(BDM_ENV_SCRIPT) $(SCRIPTS)/run_control_calibration_clean.sh \
		--excel "$(EXCEL_DATA)" \
		--cell-line $(CELL_LINE) \
		--template-input $(INPUT_CONTROL) \
		--abm-binary "$(ABM_BINARY)" \
		--bdm-env-script "$(BDM_ENV_SCRIPT)" \
		--omp-threads $(OMP_THREADS) \
		--seed $(SEED) \
		--replicates $(ABM_REPLICATES) \
		--replicate-seed-step $(ABM_REPLICATE_SEED_STEP) \
		--fit-params $(ABM_FIT_PARAMS) \
		--methods $(ABM_METHODS) \
		--da-maxiter $(ABM_DA_MAXITER) \
		--nm-maxiter $(ABM_NM_MAXITER) \
		--overgrowth-factor $(ABM_OVERGROWTH_FACTOR) \
		--max-evaluations $(ABM_MAX_EVALS) \
		--population-band-std-scale $(ABM_BAND_STD_SCALE) \
		--replicate-timeout-s $(ABM_REPLICATE_TIMEOUT_S) \
		--bo-calls $(ABM_BO_CALLS) \
		--bo-init $(ABM_BO_INIT) \
		--out-dir $(ABM_OUT_DIR)

clean:
	$(call REMOVE_RESULTS_DIR,$(RESULTS_CONTROL))

clean_calibration:
	rm -rf calibration_outputs

clean_all: clean clean_calibration
