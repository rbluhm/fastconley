# Builds the fastconley Stata plugin from the engine header shared with the
# R package (../../src/conley_core.h). Targets: linux, macos, windows.
# Needs the Armadillo headers (header-only use; no BLAS/LAPACK is linked):
#   make linux ARMA_INC=/path/to/armadillo/include
# Outputs go to ../src/fastconley_<platform>.plugin, which fastconley.pkg
# installs as fastconley.plugin on the matching platform.
CORE_DIR ?= ../../src
ARMA_INC ?= /usr/include
OUT_DIR  ?= ../src
BUILD_ID ?= $(shell git rev-parse --short HEAD 2>/dev/null || echo local)

COMMON = -std=c++11 -O3 -DNDEBUG -DARMA_64BIT_WORD -DARMA_DONT_USE_WRAPPER \
         -DARMA_DONT_USE_BLAS -DARMA_DONT_USE_LAPACK -DARMA_DONT_USE_SUPERLU \
         -DFASTCONLEY_BUILD_ID=\"$(BUILD_ID)\" -I$(CORE_DIR) -I$(ARMA_INC)

.PHONY: linux macos windows clean

linux:
	gcc -O2 -fPIC -DSYSTEM=OPUNIX -c stplugin.c -o stplugin_linux.o
	g++ $(COMMON) -DSYSTEM=OPUNIX -shared -fPIC -pthread -static-libstdc++ -static-libgcc \
	    stata_plugin.cpp stplugin_linux.o -o $(OUT_DIR)/fastconley_linux64.plugin

# Universal Mach-O (Apple Silicon + Intel).
macos:
	clang -O2 -DSYSTEM=APPLEMAC -arch arm64 -arch x86_64 -mmacosx-version-min=11.0 \
	    -c stplugin.c -o stplugin_macos.o
	clang++ $(COMMON) -DSYSTEM=APPLEMAC -bundle -arch arm64 -arch x86_64 -mmacosx-version-min=11.0 \
	    stata_plugin.cpp stplugin_macos.o -o $(OUT_DIR)/fastconley_macosx.plugin

# mingw-w64, fully static (winpthreads, libstdc++, libgcc).
windows:
	$(WIN_CC) -O2 -DSYSTEM=STWIN32 -c stplugin.c -o stplugin_win.o
	$(WIN_CXX) $(COMMON) -DSYSTEM=STWIN32 -shared -static -static-libgcc -static-libstdc++ \
	    stata_plugin.cpp stplugin_win.o -o $(OUT_DIR)/fastconley_win64.plugin

WIN_CC  ?= x86_64-w64-mingw32-gcc
WIN_CXX ?= x86_64-w64-mingw32-g++

# Load smoke test: dlopen the freshly built plugin and run its "check"
# subcommand through a minimal SPI table (no Stata needed). Usage:
#   make smoke PLUGIN=../src/fastconley_linux64.plugin
PLUGIN ?= $(OUT_DIR)/fastconley_linux64.plugin
smoke:
	$(CC) -O2 -DSYSTEM=$(SMOKE_SYSTEM) -I. load_check.c -o load_check$(EXE) $(SMOKE_LIBS)
	./load_check$(EXE) $(PLUGIN)

SMOKE_SYSTEM ?= OPUNIX
SMOKE_LIBS ?= -ldl
EXE ?=

# Standalone engine benchmark (no Stata): see bench_engine.cpp. On MSYS2 use
# the same command; BENCH_CXX defaults to g++.
BENCH_CXX ?= g++
bench:
	$(BENCH_CXX) $(COMMON) -pthread bench_engine.cpp -o bench_engine$(EXE)
	./bench_engine$(EXE) 1 4

clean:
	rm -f stplugin_*.o load_check load_check.exe bench_engine bench_engine.exe
