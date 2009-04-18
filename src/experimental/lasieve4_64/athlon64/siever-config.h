/*1:*/
#line 12 "siever-config.w"

#ifndef __SIEVER_CONFIG_H__
#define __SIEVER_CONFIG_H__
#include <stddef.h> 
#include <limits.h> 
#define L1_BITS 15
#define ULONG_RI
typedef unsigned u32_t;
typedef int i32_t;
typedef short int i16_t;
typedef unsigned short u16_t;
typedef unsigned long u64_t;
typedef long int i64_t;
#define U32_MAX 0xffffffff
#define I32_MAX INT_MAX

#define HAVE_ASM_GETBC
void asm_getbc(u32_t,u32_t,u32_t,u32_t*,u32_t*,u32_t*,u32_t*);
#define ASM_SCHEDSIEVE
#define ASM_SCHEDTDSIEVE2
u16_t**tdsieve_sched2buf(u16_t*,u16_t*,unsigned char*,u16_t**,u16_t**);

#define ASM_MPZ_TD
#define PREINVERT
#if 1
void MMX_TdAllocate(int,size_t,size_t);
u16_t*MMX_TdInit(int,u16_t*,u16_t*,u32_t*,int);
void MMX_TdUpdate(int,int);
#define MMX_TD
#define MMX_REGW 8
#endif
#define ASM_LINESIEVER
#define ASM_LINESIEVER3
#define ASM_LINESIEVER2
#define ASM_LINESIEVER1
#define ASM_TDSLINIE1
#define ASM_TDSLINIE2
#define ASM_TDSLINIE3
#define ASM_TDSLINIE
#define ASM_SEARCH0
#define HAVE_MM_LASIEVE_SETUP
#define HAVE_MM_LASIEVE_SETUP2
#define  MAX_FB_PER_P 2

u32_t*asm_lasieve_mm_setup0(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*asm_lasieve_mm_setup1(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*asm_lasieve_mm_setup2(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*asm_lasieve_mm_setup3(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*asm_lasieve_mm_setup20(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*asm_lasieve_mm_setup21(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*asm_lasieve_mm_setup22(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*asm_lasieve_mm_setup23(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);

/*:1*//*2:*/
#line 66 "siever-config.w"

#define FB_RAS 3

/*:2*//*4:*/
#line 79 "siever-config.w"

#define N_PRIMEBOUNDS 9
#endif
/*:4*/
