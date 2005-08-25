BINDIR=../../bin
INC=-I.. -I../../include $(LOCALINC)
LIBFLAGS=$(LOCALLIB)

# passed in from higher-level directory now
#CFLAGS = -Wall -O2 -I/usr/local/include
#CFLAGS = -g -Wall -Wno-implicit
#LIBFLAGS = -I. -I/usr/local/include -L/usr/local/lib

LIBS       = -lm -lgmp
LIBSPARI       = -lm -lgmp -lpari

POL5_SOURCEFILES = fnmatch.c pol51m0b.c pol51m0n.c pol51opt.c \
  if.c ../if.h assess.c primes.c roots.c zeit.c dickman.tab \
  asm_hash5.asm asm_hash5n.asm asm_rs.asm \
  asm_hash5.s asm_hash5n.s asm_rs.s \
  asm_hash5_alpha.s asm_hash5n_alpha.s \
  Example Readme parameters rsa512.data

BINS=$(BINDIR)/pol51m0b $(BINDIR)/pol51m0n $(BINDIR)/pol51opt
BINS_INTEL=$(BINDIR)/pol51m0b-intel $(BINDIR)/pol51m0n-intel $(BINDIR)/pol51opt-intel

all: $(BINS)
intel: $(BINS_INTEL)

%.o: %.c ../if.h
	$(CC) $(INC) -c $(CFLAGS) -o $@ $<

%.o: %.s
	$(CC) $(INC) -c $(CFLAGS) $^

%.s: %.asm
	gasp -c '#' -o $@ $^

$(BINDIR)/pol51m0b: pol51m0b.o ../if.o zeit.o fnmatch.o primes.o assess.o roots.o
	$(CC) $(INC) -o $@ $(CFLAGS) $^ $(LIBFLAGS) $(LIBS)

pol51m0b-intel.o: pol51m0b.c ../if.h
	$(CC) $(INC) -c $(CFLAGS) -DHAVE_FLOAT64 -DHAVE_ASM_INTEL -o $@ $<

zeit-intel.o: zeit.c ../if.h
	$(CC) $(INC) -c $(CFLAGS) -DHAVE_ASM_INTEL -o $@ $<

zeit-alpha.o: zeit.c ../if.h
	$(CC) $(INC) -c $(CFLAGS) -DHAVE_ASM_ALPHA -o $@ $<

pol51m0b-alpha.o: pol51m0b.c ../if.h
	$(CC) $(INC) -c $(CFLAGS) -DHAVE_ASM_ALPHA -o $@ $<

$(BINDIR)/pol51m0b-intel: pol51m0b-intel.o ../if.o zeit-intel.o asm_hash5.o fnmatch.o primes.o assess.o roots.o
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(LIBFLAGS) $(LIBS)
	cp -f $(BINDIR)/pol51m0b-intel $(BINDIR)/pol51m0b

$(BINDIR)/pol51m0b-alpha: pol51m0b-alpha.o ../if.o zeit-alpha.o asm_hash5_alpha.o primes.o assess.o roots.o
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(LIBFLAGS) $(LIBS)

$(BINDIR)/pol51m0n: pol51m0n.o ../if.o fnmatch.o zeit.o assess.o primes.o roots.o primes.o assess.o roots.o
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(LIBFLAGS) $(LIBS)

pol51m0n-intel.o: pol51m0n.c ../if.h
	$(CC) $(INC) -c $(CFLAGS) -DHAVE_FLOAT64 -DHAVE_ASM_INTEL -o $@ $<

pol51m0n-alpha.o: pol51m0n.c ../if.h
	$(CC) $(INC) -c $(CFLAGS) -DHAVE_ASM_ALPHA -o $@ $<

$(BINDIR)/pol51m0n-intel: pol51m0n-intel.o ../if.o zeit-intel.o fnmatch.o asm_hash5n.o primes.o assess.o roots.o
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(LIBFLAGS) $(LIBS)
	cp -f $(BINDIR)/pol51m0n-intel $(BINDIR)/pol51m0n

$(BINDIR)/pol51m0n-alpha: pol51m0n-alpha.o ../if.o zeit-alpha.o asm_hash5n_alpha.o primes.o assess.o roots.o
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(LIBFLAGS) $(LIBS)



$(BINDIR)/pol51opt: pol51opt.o ../if.o fnmatch.o primes.o assess.o roots.o zeit.o
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(LIBFLAGS) $(LIBS)

pol51opt-intel.o: pol51opt.c ../if.h
	$(CC) $(INC) -c $(CFLAGS) -DHAVE_ASM_INTEL -o $@ $<

$(BINDIR)/pol51opt-intel: pol51opt-intel.o ../if.o primes.o assess.o roots.o zeit-intel.o fnmatch.o asm_rs.o
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(LIBFLAGS) $(LIBS)
	cp -f $(BINDIR)/pol51opt-intel $(BINDIR)/pol51opt

clean : ;
	rm -f $(BINS) $(BINS_INTEL) *.o core

