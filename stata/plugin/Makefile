# Builds the fastconley Stata plugin from the engine header shared with the
# R package (../../src/conley_core.h). Armadillo is used header-only.
CORE_DIR ?= ../../src
ARMA_INC ?= /usr/include
OUT_DIR  ?= ../src
BUILD_ID ?= $(shell git rev-parse --short HEAD 2>/dev/null || echo local)

COMMON = -std=c++14 -O3 -fno-math-errno -DNDEBUG -DARMA_64BIT_WORD -DARMA_DONT_USE_WRAPPER \
         -DARMA_DONT_USE_BLAS -DARMA_DONT_USE_LAPACK -DARMA_DONT_USE_SUPERLU \
         -DFASTCONLEY_BUILD_ID=\"$(BUILD_ID)\" -I$(CORE_DIR) -I$(ARMA_INC)

LINUX_CC    ?= gcc
LINUX_CXX   ?= g++
LINUX_STRIP ?= strip
WIN_CC      ?= x86_64-w64-mingw32-gcc
WIN_CXX     ?= x86_64-w64-mingw32-g++
WIN_STRIP   ?= x86_64-w64-mingw32-strip
MAC_CC      ?= clang
MAC_CXX     ?= clang++
MAC_STRIP   ?= strip
LIPO        ?= lipo
CODESIGN    ?= codesign

LINUX_PLUGIN = $(OUT_DIR)/fastconley_linux64.plugin
WIN_PLUGIN   = $(OUT_DIR)/fastconley_win64.plugin
MAC_PLUGIN   = $(OUT_DIR)/fastconley_macosx.plugin

.PHONY: linux macos windows smoke golden_check golden check bench clean

# Keep the static C++ runtime while hiding its symbols. The version script is
# also the authoritative ELF export allow-list.
linux:
	mkdir -p $(OUT_DIR)
	$(LINUX_CC) -O2 -fPIC -DSYSTEM=OPUNIX -c stplugin.c -o stplugin_linux.o
	$(LINUX_CXX) $(COMMON) -DSYSTEM=OPUNIX -shared -fPIC -pthread \
	    -static-libstdc++ -static-libgcc \
	    -Wl,--version-script=exports.map,--exclude-libs,ALL \
	    stata_plugin.cpp stplugin_linux.o -o $(LINUX_PLUGIN).debug
	cp $(LINUX_PLUGIN).debug $(LINUX_PLUGIN)
	$(LINUX_STRIP) --strip-unneeded $(LINUX_PLUGIN)

# Build the two Mach-O slices separately so the deployment target is explicit
# on each, combine the unstripped slices, strip the release copy, then ad-hoc
# sign the final universal binary.
macos:
	mkdir -p $(OUT_DIR)
	$(MAC_CC) -O2 -DSYSTEM=APPLEMAC -arch arm64 -mmacosx-version-min=11.0 \
	    -c stplugin.c -o stplugin_macos_arm64.o
	$(MAC_CXX) $(COMMON) -DSYSTEM=APPLEMAC -bundle -arch arm64 \
	    -mmacosx-version-min=11.0 -Wl,-exported_symbols_list,exports.txt \
	    stata_plugin.cpp stplugin_macos_arm64.o -o $(MAC_PLUGIN).arm64.debug
	$(MAC_CC) -O2 -DSYSTEM=APPLEMAC -arch x86_64 -mmacosx-version-min=11.0 \
	    -c stplugin.c -o stplugin_macos_x86_64.o
	$(MAC_CXX) $(COMMON) -DSYSTEM=APPLEMAC -bundle -arch x86_64 \
	    -mmacosx-version-min=11.0 -Wl,-exported_symbols_list,exports.txt \
	    stata_plugin.cpp stplugin_macos_x86_64.o -o $(MAC_PLUGIN).x86_64.debug
	$(LIPO) -create $(MAC_PLUGIN).arm64.debug $(MAC_PLUGIN).x86_64.debug \
	    -output $(MAC_PLUGIN).debug
	rm -f $(MAC_PLUGIN).arm64.debug $(MAC_PLUGIN).x86_64.debug
	cp $(MAC_PLUGIN).debug $(MAC_PLUGIN)
	$(MAC_STRIP) -x $(MAC_PLUGIN)
	$(CODESIGN) -s - --force $(MAC_PLUGIN)

# mingw-w64: a fully static runtime, no auto-exported implementation symbols,
# and a stripped release copy. __declspec(dllexport) keeps the two SPI entries.
windows:
	mkdir -p $(OUT_DIR)
	$(WIN_CC) -O2 -DSYSTEM=STWIN32 -c stplugin.c -o stplugin_win.o
	$(WIN_CXX) $(COMMON) -DSYSTEM=STWIN32 -shared -static -static-libgcc \
	    -static-libstdc++ -Wl,--exclude-all-symbols \
	    stata_plugin.cpp stplugin_win.o -o $(WIN_PLUGIN).debug
	cp $(WIN_PLUGIN).debug $(WIN_PLUGIN)
	$(WIN_STRIP) --strip-unneeded $(WIN_PLUGIN)

# dlopen/LoadLibrary the already-stripped plugin and exercise the full SPI
# protocol against deterministic data. GOLDEN_FILE is relative to this dir.
PLUGIN ?= $(LINUX_PLUGIN)
GOLDEN_FILE ?= golden.txt
SMOKE_SYSTEM ?= OPUNIX
SMOKE_LIBS ?= -ldl
EXE ?=
smoke:
	$(CC) -O2 -DSYSTEM=$(SMOKE_SYSTEM) -I. load_check.c -o load_check$(EXE) $(SMOKE_LIBS)
	./load_check$(EXE) $(PLUGIN) $(GOLDEN_FILE)

# Header-only deterministic numerical oracle. `make golden` is intentionally
# the single regeneration command after a shared-engine change.
GOLDEN_CXX ?= g++
golden_check:
	$(GOLDEN_CXX) $(COMMON) -pthread golden_check.cpp -o golden_check$(EXE)

golden: golden_check
	./golden_check$(EXE) --write $(GOLDEN_FILE)

check: golden_check
	./golden_check$(EXE) $(GOLDEN_FILE)
	$(MAKE) smoke PLUGIN=$(PLUGIN) GOLDEN_FILE=$(GOLDEN_FILE) \
	    SMOKE_SYSTEM=$(SMOKE_SYSTEM) SMOKE_LIBS="$(SMOKE_LIBS)" EXE=$(EXE)

# Standalone engine benchmark (no Stata).
BENCH_CXX ?= g++
bench:
	$(BENCH_CXX) $(COMMON) -pthread bench_engine.cpp -o bench_engine$(EXE)
	./bench_engine$(EXE) 1 4

clean:
	rm -f stplugin_*.o load_check load_check.exe golden_check golden_check.exe \
	    bench_engine bench_engine.exe
