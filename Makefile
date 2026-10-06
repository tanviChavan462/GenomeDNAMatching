# Makefile for Parallel k-mer-Based AMR Screening Project
CC = gcc
CFLAGS = -O3 -Wall
OMP_FLAGS = -fopenmp

SRC_DIR = src
BIN_DIR = bin
DATA_DIR = data
RESULTS_DIR = results

TARGET_SEQ = $(BIN_DIR)/kmer_sequential.exe
TARGET_OMP = $(BIN_DIR)/kmer_openmp.exe

all: $(BIN_DIR) $(RESULTS_DIR) $(TARGET_SEQ) $(TARGET_OMP)

$(BIN_DIR):
	@mkdir $(BIN_DIR) 2>nul || exit 0

$(RESULTS_DIR):
	@mkdir $(RESULTS_DIR) 2>nul || exit 0

$(TARGET_SEQ): $(SRC_DIR)/kmer_sequential.c $(SRC_DIR)/amr_common.h
	$(CC) $(CFLAGS) $< -o $@

$(TARGET_OMP): $(SRC_DIR)/kmer_openmp.c $(SRC_DIR)/amr_common.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) $< -o $@

# Run sequential screening
run-seq: $(TARGET_SEQ)
	$(TARGET_SEQ) -b 5

# Run OpenMP screening with 8 threads
run-omp: $(TARGET_OMP)
	$(TARGET_OMP) -threads 8 -b 5

# Generate scaled datasets
datasets:
	python scripts/generate_datasets.py

# Run full benchmark experiments
benchmark: all datasets
	python scripts/benchmark.py

# Generate performance graphs
plots:
	python scripts/plot_results.py

# Clean build artifacts
clean:
	@del /Q $(BIN_DIR)\*.exe 2>nul || exit 0

.PHONY: all run-seq run-omp datasets benchmark plots clean
