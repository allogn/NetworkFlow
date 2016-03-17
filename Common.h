//
// Created by alvis on 16.03.16.
//

#ifndef NETWORKFLOW_COMMON_H
#define NETWORKFLOW_COMMON_H

#include <limits.h>
#include <assert.h>
#include "nheap.h"

#if defined(LONG)
typedef long intT;
typedef unsigned long uintT;
#define INT_T_MAX LONG_MAX
#define UINT_T_MAX ULONG_MAX
#else
typedef int intT;
typedef unsigned int uintT;
#define INT_T_MAX INT_MAX
#define UINT_T_MAX UINT_MAX
#endif

//edges store 32-bit quantities unless EDGELONG is defined
#if defined(EDGELONG)
typedef long intE;
typedef unsigned long uintE;
#define INT_E_MAX LONG_MAX
#define UINT_E_MAX ULONG_MAX
#else
typedef int intE;
typedef unsigned int uintE;
#define INT_E_MAX INT_MAX
#define UINT_E_MAX UINT_MAX
#endif

typedef fHeap<intT> mmHeap;

#endif //NETWORKFLOW_COMMON_H
