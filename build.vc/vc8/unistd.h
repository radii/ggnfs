
#ifndef UNISTD_H
#define UNISTD_H

#include "getopt.h"

#define _USE_MATH_DEFINES
#define strtoull   _strtoui64
#define popen      _popen
#define pclose     _pclose
#define strcasecmp _strcmpi

#ifndef inline
#define inline __inline
#endif

#define asm(x)

typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

#endif