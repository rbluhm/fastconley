/* load_check.c -- loads a fastconley Stata plugin the way Stata does
 * (dlopen / LoadLibrary), hands it a minimal ST_plugin table, runs the
 * "check" subcommand, and prints the engine version it reports. Used by the
 * CI workflow to prove each binary loads and links on its platform without
 * Stata. Build: see the Makefile target `smoke`. */
#include <stdio.h>
#include <string.h>
#include "stplugin.h"

#if SYSTEM == STWIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

static char g_version[256] = "";
static char g_build[256] = "";

static ST_retcode my_macresave(char *name, char *value) {
  if (strcmp(name, "FASTCONLEY_ENGINE_VERSION") == 0) strncpy(g_version, value, 255);
  if (strcmp(name, "FASTCONLEY_ENGINE_BUILD") == 0) strncpy(g_build, value, 255);
  return 0;
}
static ST_retcode my_spouterr(char *msg) { fputs(msg, stderr); return 0; }
static ST_retcode my_spoutsml(char *msg) { fputs(msg, stdout); return 0; }

int main(int argc, char **argv) {
  ST_plugin table;
  ST_retcode (*init)(ST_plugin *);
  ST_retcode (*call)(int, char **);
  char *args[1] = { (char *)"check" };
  ST_retcode rc;
  if (argc < 2) { fprintf(stderr, "usage: load_check <plugin>\n"); return 2; }
  memset(&table, 0, sizeof(table));
  table.macresave = my_macresave;
  table.spouterr = my_spouterr;
  table.spoutsml = my_spoutsml;
#if SYSTEM == STWIN32
  HMODULE h = LoadLibraryA(argv[1]);
  if (!h) { fprintf(stderr, "LoadLibrary failed: %lu\n", GetLastError()); return 1; }
  init = (ST_retcode (*)(ST_plugin *))GetProcAddress(h, "pginit");
  call = (ST_retcode (*)(int, char **))GetProcAddress(h, "stata_call");
#else
  void *h = dlopen(argv[1], RTLD_NOW);
  if (!h) { fprintf(stderr, "dlopen failed: %s\n", dlerror()); return 1; }
  init = (ST_retcode (*)(ST_plugin *))dlsym(h, "pginit");
  call = (ST_retcode (*)(int, char **))dlsym(h, "stata_call");
#endif
  if (!init || !call) { fprintf(stderr, "pginit/stata_call not exported\n"); return 1; }
  rc = init(&table);
  printf("pginit -> SPI %d.%d\n", (int)(rc & 0xffff), (int)(rc >> 16));
  rc = call(1, args);
  printf("check -> rc %d, engine version %s, build %s\n", (int)rc, g_version, g_build);
  return (rc == 0 && g_version[0]) ? 0 : 1;
}
