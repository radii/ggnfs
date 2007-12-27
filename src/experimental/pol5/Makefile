# temporary
HOST = x86
CC = gcc
CFLAGS = -O3 -funroll-loops -finline-functions -ftracer \
	 -fomit-frame-pointer -W -Wall -Wconversion -march=prescott -DGMP_BUG

INC = -I../.. -I../../../include
LIBS= -lm -lgmp

STAGE1_HDRS = \
       ../../if.h \
       stage1/poly_stage1.h \
       stage1/stage1_impl.h \
       stage1/stage1_inline.h \
       stage1/profile.h

STAGE1_OBJS = \
       ../../if.o \
       stage1/check.o \
       stage1/knapsack.o \
       stage1/knapsack_approx.o \
       stage1/profile.o \
       stage1/stage1.o \
       stage1/struct.o

ifeq ($(HOST),x86)
	# if not using GCC 3.4, you may remove the -O0
	CFLAGS += -O0 -DHAVE_FLOAT64 -DHAVE_ASM_INTEL
	STAGE1_OBJS += stage1/asm_hash5.o
endif

all: pol51m0b

clean:
	rm -f pol51m0b pol51m0b.exe pol51m0b.o libstage1.a stage1/*.o core

%.o: %.c $(STAGE1_HDRS)
	$(CC) $(INC) -c $(CFLAGS) -o $@ $<

%.o: %.s
	$(CC) $(INC) -c $(CFLAGS) -o $@ $<

pol51m0b: pol51m0b.o libstage1.a
	$(CC) $(CFLAGS) $(INC) -o $@ pol51m0b.o libstage1.a $(LIBS)

libstage1.a: $(STAGE1_OBJS)
	ar r libstage1.a $(STAGE1_OBJS)
