
CC = gcc 
# passed in from higher-level directory now
#CFLAGS = -Wall -O2 -I/usr/local/include
#LIBFLAGS = -I. -I/usr/local/include -L/usr/local/lib
LIBS       = -lm -lgmp
LIBSPARI       = -lm -lgmp -lpari
#CFLAGS = -g -Wall -Wno-implicit

POL5_SOURCEFILES = fnmatch.c pol51m0b.c pol51m0n.c pol51opt.c \
  if.c if.h assess.c primes.c roots.c zeit.c dickman.tab \
  asm_hash5.asm asm_hash5n.asm asm_rs.asm \
  asm_hash5.s asm_hash5n.s asm_rs.s \
  asm_hash5_alpha.s asm_hash5n_alpha.s \
  Example Readme parameters rsa512.data

BINS=pol51m0b pol51m0n pol51opt
BINS_INTEL=pol51m0b-intel pol51m0n-intel pol51opt-intel

all: $(BINS)
intel: $(BINS_INTEL)

%.o: %.c if.h
	$(CC) -c $(CFLAGS) -o $@ $<

%.o: %.s
	$(CC) -c $(CFLAGS) $^

%.s: %.asm
	gasp -c '#' -o $@ $^

pol51m0b: pol51m0b.o if.o zeit.o fnmatch.o
	$(CC) -o pol51m0b $(CFLAGS) $^ $(LIBFLAGS) $(LIBS)

pol51m0b-intel.o: pol51m0b.c if.h
	$(CC) -c $(CFLAGS) -DHAVE_FLOAT64 -DHAVE_ASM_INTEL -o $@ $<

zeit-intel.o: zeit.c if.h
	$(CC) -c $(CFLAGS) -DHAVE_ASM_INTEL -o $@ $<

zeit-alpha.o: zeit.c if.h
	$(CC) -c $(CFLAGS) -DHAVE_ASM_ALPHA -o $@ $<

pol51m0b-alpha.o: pol51m0b.c if.h
	$(CC) -c $(CFLAGS) -DHAVE_ASM_ALPHA -o $@ $<

pol51m0b-intel: pol51m0b-intel.o if.o zeit-intel.o asm_hash5.o fnmatch.o
	$(CC) -o pol51m0b-intel $(CFLAGS) $^ $(LIBFLAGS) $(LIBS)

pol51m0b-alpha: pol51m0b-alpha.o if.o zeit-alpha.o asm_hash5_alpha.o
	$(CC) -o pol51m0b-alpha $(CFLAGS) $^ $(LIBFLAGS) $(LIBS)

pol51m0n: pol51m0n.o if.o fnmatch.o zeit.o
	$(CC) -o pol51m0n $(CFLAGS) $^ $(LIBFLAGS) $(LIBS)

pol51m0n-intel.o: pol51m0n.c if.h
	$(CC) -c $(CFLAGS) -DHAVE_FLOAT64 -DHAVE_ASM_INTEL -o $@ $<

pol51m0n-alpha.o: pol51m0n.c if.h
	$(CC) -c $(CFLAGS) -DHAVE_ASM_ALPHA -o $@ $<

pol51m0n-intel: pol51m0n-intel.o if.o zeit-intel.o fnmatch.o asm_hash5n.o
	$(CC) -o pol51m0n-intel $(CFLAGS) $^ $(LIBFLAGS) $(LIBS)

pol51m0n-alpha: pol51m0n-alpha.o if.o zeit-alpha.o asm_hash5n_alpha.o
	$(CC) -o pol51m0n-alpha $(CFLAGS) $^ $(LIBFLAGS) $(LIBS)



pol51opt: pol51opt.o if.o fnmatch.o primes.o assess.o roots.o zeit.o
	$(CC) -o pol51opt $(CFLAGS) $^ $(LIBFLAGS) $(LIBS)

pol51opt-intel.o: pol51opt.c if.h
	$(CC) -c $(CFLAGS) -DHAVE_ASM_INTEL -o $@ $<

pol51opt-intel: pol51opt-intel.o if.o primes.o assess.o roots.o zeit-intel.o fnmatch.o asm_rs.o
	$(CC) -o pol51opt-intel $(CFLAGS) $^ $(LIBFLAGS) $(LIBS)



dist:
	tar cO $(POL5_SOURCEFILES) \
	Makefile | gzip --best --stdout > pol5.tar.gz

clean : ;
ifeq ($(OS),Windows_NT)
	rm -f *.exe *.o core
else
	rm -f $(BINS) $(BINS_INTEL) *.o core
endif

