// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 Alperen Yavuz
// nf_port_cpu.h -- maps the portable kernel (nf_kernel.h) onto plain C++ on the host.  The GPU self-check compares the device
// with the kernel's own filter compiled for the CPU (nf_gpu_common.cpp), and pdfilter_Test/test_gpu_filter.cpp tests it
// there.  Only the filter functions and the self-test kernel are called on the host (NF_GID() is then a plain variable).
#ifndef NF_PORT_CPU_H
#define NF_PORT_CPU_H

#include <stdint.h>
#include <math.h>

typedef int64_t i64;
typedef uint64_t u64;
typedef unsigned int nfuint;

#define NF_INLINE static inline
#define NF_G
#define NF_C const
#define NF_KERNEL static
#define NF_I64_MIN INT64_MIN
static u64 nf_cpu_gid;
#define NF_GID() nf_cpu_gid
#define NF_ATOMIC_INC(p) ((*(p))++)

#endif
