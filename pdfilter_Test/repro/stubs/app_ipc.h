// Minimal stand-in for BOINC's app_ipc.h (only what boinc_api.h above needs).
#pragma once
struct APP_INIT_DATA { char gpu_type[256]; int gpu_device_num; };
