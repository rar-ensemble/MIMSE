CXX		  := g++
CXX_FLAGS := -Wall -Wextra -O3 -Wno-deprecated-copy -DHAVE_INLINE -std=c++17#-fopenmp -D_GLIBCXX_PARALLEL #-frename-registers -march=native # -fno-signed-zeros -fno-trapping-math -funroll-loops #-fprofile-dir="./gcda" -fprofile-use#-fprofile-generate #-fprofile-use #-Ofast #-std=c++17
PY		  := python3

BIN		:= bin
SRC		:= src
TEST_SRC:= test
INCLUDE	:= include
EIGEN_INCLUDE := include/eigen
SPECTRA_INCLUDE := include/spectra/include
FUNC_INCLUDE := include/func
LIB		:= lib
OUTPUT	:= output
GRAPHS	:= output/graphs
GCDA	:= gcda
JOBS	:= Q_job_creation/pysrc

LIBRARIES	:= -lgsl -lgslcblas
EXECUTABLE	:= main
TEST_EXECUTABLE := test
JOBS_CREATE	:= create_jobs_Q.py


main: $(BIN)/$(EXECUTABLE)

test: $(BIN)/$(TEST_EXECUTABLE)

build: clean_main clean_output main
	./$(BIN)/$(EXECUTABLE)

run: clean_output main
	./$(BIN)/$(EXECUTABLE)

check: clean_test test
	./$(BIN)/$(TEST_EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -I$(EIGEN_INCLUDE) -I$(SPECTRA_INCLUDE) -I$(FUNC_INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(TEST_EXECUTABLE): $(TEST_SRC)/*.cpp
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -I$(EIGEN_INCLUDE) -I$(SPECTRA_INCLUDE) -I$(FUNC_INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

clean_main:
	-rm -f $(BIN)/$(EXECUTABLE)

clean_test:
	-rm -f $(BIN)/$(TEST_EXECUTABLE)

clean:
	-rm -f $(BIN)/* $(OUTPUT)/* $(GRAPHS)/*

clean_output:
	-rm -f $(OUTPUT)/* $(GRAPHS)/*

clean_gcda:
	-rm -f $(GCDA)/*

job: main
	$(PY) $(JOBS)/$(JOBS_CREATE)


