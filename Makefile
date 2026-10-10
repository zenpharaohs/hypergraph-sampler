# Builds the C core as a shared library (used by the Python binding) and a static library.
#   make            -> build/libhypersample.{so|dylib}, build/libhypersample.a
#   make test       -> C smoke test + Python tests
#   make clean
CC      ?= cc
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra -Wno-unused-function
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SHLIB := build/libhypersample.dylib
  SHFLAGS := -dynamiclib -install_name @rpath/libhypersample.dylib
else
  SHLIB := build/libhypersample.so
  SHFLAGS := -shared
endif
SRC := $(wildcard src/*.c)
OBJ := $(patsubst src/%.c,build/obj/%.o,$(SRC))

all: $(SHLIB) build/libhypersample.a

build/obj/%.o: src/%.c include/hypersample.h src/hs_internal.h
	@mkdir -p build/obj
	$(CC) $(CFLAGS) -fPIC -Iinclude -c $< -o $@

$(SHLIB): $(OBJ)
	$(CC) $(SHFLAGS) -o $@ $(OBJ) -lm

build/libhypersample.a: $(OBJ)
	ar rcs $@ $(OBJ)

build/smoke: tests/smoke.c build/libhypersample.a
	$(CC) $(CFLAGS) -Iinclude tests/smoke.c build/libhypersample.a -lm -o $@

test: all build/smoke
	./build/smoke
	PYTHONPATH=python python3 -m pytest -q tests

clean:
	rm -rf build

.PHONY: all test clean
