BINDIR=../../bin
INC=-I.. -I../../include $(LOCALINC)
LIBFLAGS=$(LOCALLIB)
LIBS=-lm -lgmp
CFLAGS2=

# passed in from higher-level directory now
#CFLAGS = -Wall -O2 -I/usr/local/include
#CFLAGS = -g -Wall -Wno-implicit
#LIBFLAGS = -I. -I/usr/local/include -L/usr/local/lib

POL5_SOURCEFILES = fnmatch.c pol51m0b.c pol51m0n.c pol51opt.c \
  if.c ../if.h assess.c primes.c roots.c zeit.c dickman.tab \
  asm_hash5.asm asm_hash5n.asm asm_rs.asm \
  asm_hash5.s asm_hash5n.s asm_rs.s \
  asm_hash5_alpha.s asm_hash5n_alpha.s \
  Example Readme parameters rsa512.data

BINS=$(BINDIR)/pol51m0b $(BINDIR)/pol51m0n $(BINDIR)/pol51opt

OBJS=../if.o zeit.o fnmatch.o primes.o assess.o roots.o
OBJS2=
OBJS3=
OBJS4=

ifeq ($(HOST),x86)
	CFLAGS2+=-DHAVE_FLOAT64 -DHAVE_ASM_INTEL
	OBJS2=asm_hash5.o
	OBJS3=asm_hash5n.o
	OBJS4=asm_rs.o
endif

ifeq ($(HOST),alpha)
	CFLAGS2+=-DHAVE_ASM_ALPHA
	OBJS2=asm_hash5_alpha.o
	OBJS3=asm_hash5n_alpha.o
endif

all: $(BINS)

%.o: %.c ../if.h
	$(CC) $(INC) -c $(CFLAGS) $(CFLAGS2) -o $@ $<

%.o: %.s
	$(CC) $(INC) -c $(CFLAGS) $(CFLAGS2) $^

%.s: %.asm
	gasp -c '#' -o $@ $^

$(BINDIR)/pol51m0b: pol51m0b.o $(OBJS) $(OBJS2)
	$(CC) $(INC) -o $@ $(CFLAGS) $(CFLAGS2) $^ $(LIBFLAGS) $(LIBS)

$(BINDIR)/pol51m0n: pol51m0n.o $(OBJS) $(OBJS3)
	$(CC) $(INC) $(CFLAGS) $(CFLAGS2) -o $@ $^ $(LIBFLAGS) $(LIBS)

$(BINDIR)/pol51opt: pol51opt.o $(OBJS) $(OBJS4)
	$(CC) $(INC) $(CFLAGS) $(CFLAGS2) -o $@ $^ $(LIBFLAGS) $(LIBS)

clean : ;
	rm -f $(BINS) *.o core

