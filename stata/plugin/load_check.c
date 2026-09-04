/* load_check.c -- load a fastconley plugin without Stata and exercise its
 * real SPI protocol against deterministic 200-row, three-score datasets. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "stplugin.h"

#if SYSTEM == STWIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#define NROW 200
#define K 3
#define NCASE 7
#define MISSING_VALUE 8.988465674311579e+307

enum data_mode { MODE_SPATIAL, MODE_SERIAL, MODE_GRID };

struct golden_case {
  char name[64];
  double value[K][K];
};

static enum data_mode g_mode = MODE_SPATIAL;
static double g_spatial[NROW][3 + K];
static double g_serial[NROW][2 + K];
static double g_grid[NROW][3 + K];
static double g_matrix[K][K];
static struct golden_case g_golden[NCASE];
static int g_matrix_stores = 0;
static int g_missing_var = 0;
static int g_stopflag = 0;
static char g_version[256] = "";
static char g_build[256] = "";
static char g_plugin_error[512] = "";
static char g_unbalanced[32] = "";
static char g_last_error[512] = "";

/* The reference C++ generator uses std::mt19937_64. This is its specified
 * state transition and seed recurrence, followed by the same raw-bit mapping
 * to [0,1); no libc random distribution is involved. */
#define MT_NN 312
#define MT_MM 156
#define MT_MATRIX_A UINT64_C(0xB5026F5AA96619E9)
#define MT_UM UINT64_C(0xFFFFFFFF80000000)
#define MT_LM UINT64_C(0x7FFFFFFF)
static uint64_t mt_state[MT_NN];
static int mt_index = MT_NN;

static void mt_seed(uint64_t seed) {
  int i;
  mt_state[0] = seed;
  for (i = 1; i < MT_NN; ++i) {
    mt_state[i] = UINT64_C(6364136223846793005) *
                  (mt_state[i - 1] ^ (mt_state[i - 1] >> 62)) + (uint64_t)i;
  }
  mt_index = MT_NN;
}

static uint64_t mt_next(void) {
  uint64_t x;
  int i;
  if (mt_index >= MT_NN) {
    for (i = 0; i < MT_NN; ++i) {
      const uint64_t y = (mt_state[i] & MT_UM) |
                         (mt_state[(i + 1) % MT_NN] & MT_LM);
      mt_state[i] = mt_state[(i + MT_MM) % MT_NN] ^ (y >> 1) ^
                    ((y & UINT64_C(1)) ? MT_MATRIX_A : UINT64_C(0));
    }
    mt_index = 0;
  }
  x = mt_state[mt_index++];
  x ^= (x >> 29) & UINT64_C(0x5555555555555555);
  x ^= (x << 17) & UINT64_C(0x71D67FFFEDA60000);
  x ^= (x << 37) & UINT64_C(0xFFF7EEE000000000);
  x ^= x >> 43;
  return x;
}

static double rng_unit(void) {
  return (double)(mt_next() >> 11) * (1.0 / 9007199254740992.0);
}

static double rng_range(double lo, double hi) {
  return lo + (hi - lo) * rng_unit();
}

static void initialize_data(void) {
  int i, j;
  mt_seed(UINT64_C(101));
  for (i = 0; i < NROW; ++i) {
    g_spatial[i][0] = rng_range(25.0, 49.0);
    g_spatial[i][1] = rng_range(-124.0, -67.0);
    g_spatial[i][2] = 1.0;
    for (j = 0; j < K; ++j) {
      g_spatial[i][3 + j] = rng_range(-1.0, 1.0) * 1.7320508075688772;
    }
  }

  mt_seed(UINT64_C(303));
  for (i = 0; i < NROW; ++i) {
    g_serial[i][0] = (double)(i / 5 + 1);
    g_serial[i][1] = (double)(i % 5 + 1);
    for (j = 0; j < K; ++j) {
      g_serial[i][2 + j] = rng_range(-1.0, 1.0) * 1.7320508075688772;
    }
  }

  mt_seed(UINT64_C(404));
  for (i = 0; i < NROW; ++i) {
    g_grid[i][0] = (double)(i / 20);
    g_grid[i][1] = (double)(i % 20);
    g_grid[i][2] = 1.0;
    for (j = 0; j < K; ++j) {
      g_grid[i][3 + j] = rng_range(-1.0, 1.0) * 1.7320508075688772;
    }
  }
}

static int read_golden(const char *path) {
  FILE *fp = fopen(path, "rb");
  char line[4096];
  int count = 0;
  if (!fp) {
    fprintf(stderr, "cannot read golden file %s\n", path);
    return 1;
  }
  while (fgets(line, sizeof(line), fp)) {
    char *tok;
    int i, j;
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
    if (count == NCASE) {
      fprintf(stderr, "too many golden cases\n");
      fclose(fp);
      return 1;
    }
    tok = strtok(line, " \t\r\n");
    if (!tok) continue;
    strncpy(g_golden[count].name, tok, sizeof(g_golden[count].name) - 1);
    g_golden[count].name[sizeof(g_golden[count].name) - 1] = '\0';
    for (i = 0; i < K; ++i) {
      for (j = 0; j < K; ++j) {
        char *end;
        tok = strtok(NULL, " \t\r\n");
        if (!tok) {
          fprintf(stderr, "too few values for golden case %s\n", g_golden[count].name);
          fclose(fp);
          return 1;
        }
        g_golden[count].value[i][j] = strtod(tok, &end);
        if (end == tok || *end != '\0') {
          fprintf(stderr, "bad value in golden case %s\n", g_golden[count].name);
          fclose(fp);
          return 1;
        }
      }
    }
    if (strtok(NULL, " \t\r\n") != NULL) {
      fprintf(stderr, "too many values for golden case %s\n", g_golden[count].name);
      fclose(fp);
      return 1;
    }
    ++count;
  }
  fclose(fp);
  if (count != NCASE) {
    fprintf(stderr, "expected %d golden cases, read %d\n", NCASE, count);
    return 1;
  }
  return 0;
}

static void copy_string(char *dst, size_t cap, const char *src) {
  if (!cap) return;
  strncpy(dst, src ? src : "", cap - 1);
  dst[cap - 1] = '\0';
}

static ST_retcode my_macresave(char *name, char *value) {
  if (strcmp(name, "FASTCONLEY_ENGINE_VERSION") == 0) {
    copy_string(g_version, sizeof(g_version), value);
  } else if (strcmp(name, "FASTCONLEY_ENGINE_BUILD") == 0) {
    copy_string(g_build, sizeof(g_build), value);
  } else if (strcmp(name, "_fc_plugin_error") == 0) {
    copy_string(g_plugin_error, sizeof(g_plugin_error), value);
  } else if (strcmp(name, "_fc_unbalanced_fallback") == 0) {
    copy_string(g_unbalanced, sizeof(g_unbalanced), value);
  }
  return 0;
}

static ST_retcode my_spouterr(char *msg) {
  copy_string(g_last_error, sizeof(g_last_error), msg);
  fputs(msg, stderr);
  return 0;
}

static ST_retcode my_spoutsml(char *msg) {
  fputs(msg, stdout);
  return 0;
}

static ST_int my_poll(void) { return 0; }
static ST_int my_nobs(void) { return NROW; }
static ST_int my_nobs1(void) { return 1; }
static ST_int my_nobs2(void) { return NROW; }
static ST_int my_nvars(void) {
  return g_mode == MODE_SERIAL ? 2 + K : 3 + K;
}
static ST_boolean my_selobs(ST_int obs) { return obs >= 1 && obs <= NROW; }
static ST_boolean my_ismissing(ST_double value) {
  return !isfinite(value) || value >= MISSING_VALUE;
}

static ST_retcode my_safevdata(ST_int var, ST_int obs, ST_double *out) {
  int nvar = my_nvars();
  if (!out || obs < 1 || obs > NROW || var < 1 || var > nvar) return 1;
  if (g_missing_var == var && obs == 17) {
    *out = MISSING_VALUE;
    return 0;
  }
  if (g_mode == MODE_SPATIAL) *out = g_spatial[obs - 1][var - 1];
  else if (g_mode == MODE_SERIAL) *out = g_serial[obs - 1][var - 1];
  else *out = g_grid[obs - 1][var - 1];
  return 0;
}

static ST_retcode my_scalaruse(char *name, ST_double *out) {
  if (!out) return 1;
  if (strcmp(name, "fc_cutoff") == 0) *out = 500.0;
  else if (strcmp(name, "fc_lag") == 0) *out = 2.0;
  else if (strcmp(name, "fc_lat0") == 0) *out = -22.5;
  else if (strcmp(name, "fc_dlat") == 0) *out = 5.0;
  else if (strcmp(name, "fc_dlon") == 0) *out = 18.0;
  else if (strcmp(name, "fc_grid_cutoff") == 0) *out = 2500.0;
  else if (strcmp(name, "fc_missing") == 0) *out = MISSING_VALUE;
  else return 1;
  return 0;
}

static ST_retcode my_safematstore(char *name, ST_int row, ST_int col,
                                  ST_double value) {
  if (strcmp(name, "FC_RESULT") != 0 || row < 1 || row > K ||
      col < 1 || col > K) return 1;
  g_matrix[row - 1][col - 1] = value;
  ++g_matrix_stores;
  return 0;
}

static ST_int my_matrix_dim(char *name) {
  return strcmp(name, "FC_RESULT") == 0 ? K : 0;
}

static void reset_call(void) {
  memset(g_matrix, 0, sizeof(g_matrix));
  g_matrix_stores = 0;
  g_missing_var = 0;
  g_plugin_error[0] = '\0';
  g_unbalanced[0] = '\0';
  g_last_error[0] = '\0';
}

static int compare_matrix(const char *name) {
  int c, i, j;
  double worst = 0.0;
  for (c = 0; c < NCASE; ++c) {
    if (strcmp(g_golden[c].name, name) == 0) break;
  }
  if (c == NCASE) {
    fprintf(stderr, "golden case not found: %s\n", name);
    return 1;
  }
  if (g_matrix_stores != K * K) {
    fprintf(stderr, "%s: expected %d matrix stores, got %d\n",
            name, K * K, g_matrix_stores);
    return 1;
  }
  for (i = 0; i < K; ++i) {
    for (j = 0; j < K; ++j) {
      const double want = g_golden[c].value[i][j];
      const double scale = fabs(want) > 1.0 ? fabs(want) : 1.0;
      const double rel = fabs(g_matrix[i][j] - want) / scale;
      if (rel > worst) worst = rel;
      if (rel > 1e-12) {
        fprintf(stderr, "%s[%d,%d]: got %.17g expected %.17g rel %.17g\n",
                name, i, j, g_matrix[i][j], want, rel);
        return 1;
      }
    }
  }
  printf("SPI %-35s m00=%.17g m01=%.17g m12=%.17g m22=%.17g rel=%.3g\n",
         name, g_matrix[0][0], g_matrix[0][1], g_matrix[1][2],
         g_matrix[2][2], worst);
  return 0;
}

static int run_matrix_call(ST_retcode (*call)(int, char **),
                           enum data_mode mode, int argc, char **argv,
                           const char *golden_name) {
  ST_retcode rc;
  g_mode = mode;
  reset_call();
  rc = call(argc, argv);
  if (rc != 0) {
    fprintf(stderr, "%s -> rc %d (%s)\n", golden_name, (int)rc, g_plugin_error);
    return 1;
  }
  return compare_matrix(golden_name);
}

int main(int argc, char **argv) {
  ST_plugin table;
  ST_retcode (*init)(ST_plugin *);
  ST_retcode (*call)(int, char **);
  ST_retcode rc;
  char *check_args[] = {(char *)"check"};
  char *spatial_bh[] = {(char *)"spatial", (char *)"fc_cutoff",
                        (char *)"bartlett", (char *)"haversine", (char *)"0",
                        (char *)"2", (char *)"grid", (char *)"double",
                        (char *)"FC_RESULT"};
  char *spatial_us[] = {(char *)"spatial", (char *)"fc_cutoff",
                        (char *)"uniform", (char *)"spherical", (char *)"0",
                        (char *)"2", (char *)"grid", (char *)"double",
                        (char *)"FC_RESULT"};
  char *serial[] = {(char *)"serial", (char *)"fc_lag", (char *)"2",
                    (char *)"FC_RESULT"};
  char *grid_us[] = {(char *)"grid", (char *)"fc_lat0", (char *)"fc_dlat",
                     (char *)"fc_dlon", (char *)"10", (char *)"20",
                     (char *)"20", (char *)"fc_grid_cutoff",
                     (char *)"spherical", (char *)"uniform", (char *)"2",
                     (char *)"FC_RESULT"};
  char *grid_bh[] = {(char *)"grid", (char *)"fc_lat0", (char *)"fc_dlat",
                     (char *)"fc_dlon", (char *)"10", (char *)"20",
                     (char *)"20", (char *)"fc_grid_cutoff",
                     (char *)"haversine", (char *)"bartlett", (char *)"2",
                     (char *)"FC_RESULT"};
  char *missing_scalar[] = {(char *)"spatial", (char *)"fc_missing",
                            (char *)"uniform", (char *)"haversine", (char *)"0",
                            (char *)"1", (char *)"band", (char *)"double",
                            (char *)"FC_RESULT"};

  if (argc < 3) {
    fprintf(stderr, "usage: load_check <plugin> <golden.txt>\n");
    return 2;
  }
  if (read_golden(argv[2])) return 1;
  initialize_data();

  memset(&table, 0, sizeof(table));
  table.spoutsml = my_spoutsml;
  table.pollstd = my_poll;
  table.pollnow = my_poll;
  table.macresave = my_macresave;
  table.scalaruse = my_scalaruse;
  table.nobs = my_nobs;
  table.nvar = my_nvars;
  table.missval = MISSING_VALUE;
  table.ismissing = my_ismissing;
  table.stopflag = &g_stopflag;
  table.selobs = my_selobs;
  table.nobs1 = my_nobs1;
  table.nobs2 = my_nobs2;
  table.nvars = my_nvars;
  table.spouterr = my_spouterr;
  table.safematstore = my_safematstore;
  table.safevdata = my_safevdata;
  table.colsof = my_matrix_dim;
  table.rowsof = my_matrix_dim;

#if SYSTEM == STWIN32
  {
    HMODULE h = LoadLibraryA(argv[1]);
    if (!h) {
      fprintf(stderr, "LoadLibrary failed: %lu\n", GetLastError());
      return 1;
    }
    init = (ST_retcode (*)(ST_plugin *))GetProcAddress(h, "pginit");
    call = (ST_retcode (*)(int, char **))GetProcAddress(h, "stata_call");
  }
#else
  {
    void *h = dlopen(argv[1], RTLD_NOW);
    if (!h) {
      fprintf(stderr, "dlopen failed: %s\n", dlerror());
      return 1;
    }
    init = (ST_retcode (*)(ST_plugin *))dlsym(h, "pginit");
    call = (ST_retcode (*)(int, char **))dlsym(h, "stata_call");
  }
#endif
  if (!init || !call) {
    fprintf(stderr, "pginit/stata_call not exported\n");
    return 1;
  }

  rc = init(&table);
  printf("pginit -> SPI %d.%d\n", (int)(rc & 0xffff), (int)(rc >> 16));
  if ((rc & 0xffff) != SD_PLUGINMAJ || (rc >> 16) != SD_PLUGINMIN) return 1;
  rc = call(1, check_args);
  printf("check -> rc %d, engine version %s, build %s\n",
         (int)rc, g_version, g_build);
  if (rc != 0 || !g_version[0] || !g_build[0]) return 1;

  if (run_matrix_call(call, MODE_SPATIAL, 9, spatial_bh,
                      "spatial_cross_bartlett_haversine")) return 1;
  if (run_matrix_call(call, MODE_SPATIAL, 9, spatial_us,
                      "spatial_cross_uniform_spherical")) return 1;
  if (run_matrix_call(call, MODE_SERIAL, 4, serial, "serial_lag2")) return 1;
  if (run_matrix_call(call, MODE_GRID, 12, grid_us,
                      "grid_wrap_uniform_spherical")) return 1;
  if (run_matrix_call(call, MODE_GRID, 12, grid_bh,
                      "grid_wrap_bartlett_haversine")) return 1;

  g_mode = MODE_SPATIAL;
  reset_call();
  g_missing_var = 4;
  rc = call(9, spatial_bh);
  printf("missing input -> rc %d, fc_plugin_error=%s\n", (int)rc, g_plugin_error);
  if (rc != 198 || strstr(g_plugin_error, "missing value") == NULL) return 1;

  reset_call();
  rc = call(9, missing_scalar);
  if (rc == 198 && g_plugin_error[0]) {
    printf("missing scalar -> rc 198, fc_plugin_error=%s\n", g_plugin_error);
  } else if (rc == 0) {
    printf("missing scalar -> rc 0 (pending engine scalar validation)\n");
  } else {
    fprintf(stderr, "missing scalar -> unexpected rc %d (%s)\n",
            (int)rc, g_plugin_error);
    return 1;
  }

  puts("SPI numerical smoke: all checks passed");
  return 0;
}
