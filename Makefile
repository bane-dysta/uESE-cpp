CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS += -Iinclude

SOURCES = src/main.cpp src/io.cpp src/model.cpp src/output.cpp src/parameters.cpp src/surface.cpp
HEADERS = include/uese/io.hpp include/uese/model.hpp include/uese/parameters.hpp include/uese/surface.hpp src/solvents_data.inc

.PHONY: all test clean
all: uese

uese: $(SOURCES) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SOURCES) -o $@

test: uese
	python3 tests/run_tests.py ./uese

clean:
	rm -f uese uese.exe
