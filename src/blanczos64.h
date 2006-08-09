#ifndef __BLANCZOS64_H__
#define __BLANCZOS64_H__

#define _ASMx86_32 /* Please disable it for 64 bit mode binary.  */
#define _MMX_REGS
#define _HW_BSFL   /* Work very slow on Opteron CPU.  */

#if defined (_MMX_REGS)
  #include <xmmintrin.h>
  #define m64 __m64
  #define m_clear() _mm_setzero_si64 ()
  #define m_empty() _m_empty ()
  static __inline \
  m64 mem64xor (m64 v, const void* p)
  {
    return _m_pxor (v, *(m64*)p);
  }
  static __inline \
  m64 var64xor (m64 a, m64 b)
  {
    return _m_pxor (a, b);
  }

#else // _MMX_REGS
  #define m64 u64
  #define m_clear() (0ULL)
  #define m_empty()

  static __inline \
  m64 mem64xor (m64 v, const void* p)
  {
    return (v ^ *(m64*)p);
  }
  static __inline \
  m64 var64xor (m64 a, m64 b)
  {
    return (a ^ b);
  }
#endif // _MMX_REGS

#if defined (_HW_BSFL)
  #if defined (_MSC_VER)
    #include <intrin.h>
    #pragma intrinsic(_BitScanForward)
    #define bsfl _BitScanForward
  #elif defined (__GNUC__)
    static __inline \
    void bsfl (u32 *i, u32 x)
    {
      u32 r;
      asm ("bsfl %1, %0" : "=r" (r): "r" (x));
      *i = r;
    }
  #endif
#else
  /* Assemler insruction 'bsfl' work very slow on Opteron CPU,
     I realy not know why... :-(
     Next implementation of 'De Bruijn sequences', as alternative.  */
  static const u32 db = 0x4653ADFUL; /* == 00000100011001010011101011011111.  */
  static const u32 s  = 32-5;
  static const u32 dbt[32] = {
    0x00,0x01,0x02,0x06,0x03,0x0B,0x07,0x10,
    0x04,0x0E,0x0C,0x15,0x08,0x17,0x11,0x1A,
    0x1F,0x05,0x0A,0x0F,0x0D,0x14,0x16,0x19,
    0x1E,0x09,0x13,0x18,0x1D,0x12,0x1C,0x1B
  };
  static __inline \
  void bsfl (u32 *i, u32 x)
  {
    x &= -(s32)x;  /* isolate lowest bit.  */
    x *= db;       /* multiplication by a power of two is a shift.  */
    x >>= s;       /* use log_2(BITS_PER_LONG) highest bits.  */

    *i = dbt[x];   /* lookup.  */
  }
#endif

#endif // __BLANCZOS64_H__
