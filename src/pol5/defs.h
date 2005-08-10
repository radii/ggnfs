#ifndef __DEFS_H__
#define __DEFS_H__
#define uchar unsigned char
#define uint  unsigned int
#define ulong unsigned long
#define ushort unsigned short
#define ull unsigned long long

#define __STDC_FORMAT_MACROS
#if defined (_MSC_VER)
#include <basetsd.h>

#define int8_t  INT8
#define uint8_t UINT8

#define int16_t  INT16
#define uint16_t UINT16

#define int32_t INT32
#define uint32_t UINT32

#define int64_t INT64
#define uint64_t UINT64

#define PRId8 "I8d"
#define PRIi8 "I8i"
#define PRIo8 "I8o"
#define PRIu8 "I8u"
#define PRIx8 "I8x"
#define PRIX8 "I8X"

#define PRId16 "I16d"
#define PRIi16 "I16i"
#define PRIo16 "I16o"
#define PRIu16 "I16u"
#define PRIx16 "I16x"
#define PRIX16 "I16X"

#define PRId32 "I32d"
#define PRIi32 "I32i"
#define PRIo32 "I32o"
#define PRIu32 "I32u"
#define PRIx32 "I32x"
#define PRIX32 "I32X"

#define PRId64 "I64d"
#define PRIi64 "I64i"
#define PRIo64 "I64o"
#define PRIu64 "I64u"
#define PRIx64 "I64x"
#define PRIX64 "I64X"

#define SCNd32 "I32d"
#define SCNi32 "I32i"
#define SCNo32 "I32o"
#define SCNu32 "I32u"
#define SCNx32 "I32x"
#define SCNX32 "I32X"

#define SCNd64 "I64d"
#define SCNi64 "I64i"
#define SCNo64 "I64o"
#define SCNu64 "I64u"
#define SCNx64 "I64x"
#define SCNX64 "I64X"

#else
#include <inttypes.h>
#endif
#endif
