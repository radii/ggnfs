# temporary
HOST = x86
CC = gcc
CFLAGS = -O3 -funroll-loops -finline-functions -ftracer \
	 -fomit-frame-pointer -W -Wall -Wconversion -march=prescott -DGMP_BUG

INC = -Iinclude -I../.. -I../../../include
LIBS= -lm -lgmp

ifeq ($(HOST),x86)
	# if not using GCC 3.4, you may remove the -O0
	CFLAGS += -O0 -DHAVE_FLOAT64 -DHAVE_ASM_INTEL
endif

#---------------------------------------------
COMMON_HDRS = \
       ../../if.h \
       include/profile.h

COMMON_OBJS = \
       ../../if.o \
       common/profile.o

#---------------------------------------------
STAGE1_HDRS = \
       include/poly_stage1.h \
       stage1/stage1_impl.h \
       stage1/stage1_inline.h

STAGE1_OBJS = \
       stage1/check.o \
       stage1/knapsack.o \
       stage1/knapsack_approx.o \
       stage1/stage1.o \
       stage1/struct.o

ifeq ($(HOST),x86)
	STAGE1_OBJS += stage1/asm_hash5.o
endif

#---------------------------------------------
STAGE2_HDRS = \
	include/poly_stage2.h \
	stage2/stage2_impl.h \
	stage2/stage2_inline.h

STAGE2_OBJS = \
	stage2/assess.o \
	stage2/check.o \
	stage2/dickman.o \
	stage2/primes.o \
	stage2/root_sieve.o \
	stage2/roots.o \
	stage2/stage2.o \
	stage2/translate.o

ifeq ($(HOST),x86)
	STAGE2_OBJS += stage2/asm_rs.o
endif

#---------------------------------------------


all: pol51m0b pol51opt

clean:
	rm -f pol51m0b pol51m0b.exe pol51opt pol51opt.exe \
		pol51m0b.o libstage1.a pol51opt.o libstage2.a \
		$(COMMON_OBJS) $(STAGE1_OBJS) $(STAGE2_OBJS) core

%.o: %.c $(COMMON_HDRS) $(STAGE1_HDRS) $(STAGE2_HDRS)
	$(CC) $(INC) -c $(CFLAGS) -o $@ $<

%.o: %.s
	$(CC) $(INC) -c $(CFLAGS) -o $@ $<

pol51m0b: pol51m0b.o $(COMMON_OBJS) libstage1.a
	$(CC) $(CFLAGS) $(INC) -o $@ pol51m0b.o \
		$(COMMON_OBJS) libstage1.a $(LIBS)

pol51opt: pol51opt.o $(COMMON_OBJS) libstage2.a
	$(CC) $(CFLAGS) $(INC) -o $@ pol51opt.o \
		$(COMMON_OBJS) libstage2.a $(LIBS)

libstage1.a: $(STAGE1_OBJS)
	ar r libstage1.a $(STAGE1_OBJS)

libstage2.a: $(STAGE2_OBJS)
	ar r libstage2.a $(STAGE2_OBJS)
