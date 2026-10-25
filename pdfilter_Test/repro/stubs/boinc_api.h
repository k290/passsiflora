// Minimal stand-in for BOINC's boinc_api.h, so that the application can be built without the BOINC libraries.
#pragma once
#include <cstdlib>
#include <cstring>
#include "app_ipc.h"
static inline int  boinc_time_to_checkpoint() { return 0; }
static inline void boinc_checkpoint_completed() {}
static inline void boinc_fraction_done(double) {}
static inline void boinc_begin_critical_section() {}
static inline void boinc_end_critical_section() {}
static inline void boinc_finish(int rc) { std::exit(rc); }
static inline void boinc_get_init_data(APP_INIT_DATA &a) {
  memset(&a, 0, sizeof a);
  const char *t = std::getenv("NF_STUB_GPU_TYPE");
  strcpy(a.gpu_type, t ? t : "NVIDIA");
  a.gpu_device_num = 0;
}
