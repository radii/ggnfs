ALLOPT=-march=$(ARCH) -ftracer -pipe
INC=-I. -I.. -I../include $(LOCALINC)
LIBS=-lgmp -lm
BINDIR=../bin
LIBFLAGS=$(LOCALLIB)

DLMALLOC=0
MATBUILD_TPIE=0
TPIE_DIR=f:/tmp/nfs/tpie

#Hope, we're running Mingw32 ;)
ifeq ($(OSTYPE),msys)
    DLMALLOC=1
endif

ifeq ($(HOST),generic)
    ALLOPT+= -DGGNFS_HOST_GENERIC
endif

# -DL2_CACHE_SIZE=512 works well on some Pentium 4 processors.
# It brings the opposite effect on AMD chips and some Pentium 4 processors.
ifeq ($(ARCH),pentium4)
    ALLOPT+= -DL2_CACHE_SIZE=512
endif

DEBUGOPT=-W -Wall -g -pg

CFLAGS=-O3 -ffast-math -funroll-loops -finline-functions \
       -fomit-frame-pointer -W -Wall $(ALLOPT)
#CFLAGS=$(DEBUGOPT) $(ALLOPT)


export ARCH HOST ALLOPT CFLAGS DEBUGOPT

OBJS=getprimes.o fbmisc.o squfof.o rels.o $(LANCZOS).o poly.o mpz_poly.o \
     mpz_mat.o smintfact.o misc.o ecm4c.o nfmisc.o montgomery_sqrt.o \
     matstuff.o dickman.o fbgen.o llist.o if.o lasieve4/mpz-ull.o

BINS=$(BINDIR)/sieve $(BINDIR)/procrels $(BINDIR)/sqrt $(BINDIR)/polyselect \
     $(BINDIR)/makefb $(BINDIR)/matsolve $(BINDIR)/matbuild

LSBINS=latsiever polsel

ifeq ($(HOST),x86)
  LANCZOS=blanczos64
  OBJS+=modinv1002.o mulmod32.o
else
  LANCZOS=blanczos64-no-mmx
  OBJS+=modinv32.o
endif

ifeq ($(DLMALLOC),1)
  ALLOPT+= -DMALLOC_REPORTING
  OBJS+= dlmalloc.o
endif

ifeq ($(MATBUILD_TPIE),1)
  BINS+= $(BINDIR)/matbuild-tpie
  INC+= -I$(TPIE_DIR)/include
  LIBFLAGS+= -L$(TPIE_DIR)/lib
  LIBS+=-ltpie
endif

.PHONY: all tests bins latsiever polsel strip clean lasieve-clean \
        polsel-clean squeaky

all : $(OBJS) $(BINS) $(LSBINS) strip
# all : $(OBJS) $(BINS) $(LSBINS)
tests : $(TESTS)
bins : $(BINS)

.c.o :
	$(CC) $(INC) $(CFLAGS) -o $@ -c $*.c

.s.o :
	$(AS) -o $@ $*.s

$(BINDIR)/sieve : sieve.c clsieve.c $(OBJS) makefb.o
	$(CC) $(INC) $(CFLAGS) $(LIBFLAGS) -o $@ $^ $(LIBS)

$(BINDIR)/makefb : makefb.c $(OBJS)
	$(CC) $(INC) $(CFLAGS) $(LIBFLAGS) -D_MAKEFB_STANDALONE -o $@ \
	makefb.c $(OBJS) $(LIBS)

$(BINDIR)/procrels : procrels.c $(OBJS)
	$(CC) $(INC) $(CFLAGS) $(LIBFLAGS) -o $@ procrels.c $(OBJS) $(LIBS)

$(BINDIR)/matbuild : matbuild.c $(OBJS)
	$(CC) $(INC) $(CFLAGS) $(LIBFLAGS) -o $@ matbuild.c combparts.c \
	  $(OBJS) $(LIBS)

$(BINDIR)/matbuild-tpie : matbuild-tpie.c $(OBJS)
	$(CPP) $(INC) $(CFLAGS) $(LIBFLAGS) -o $@ matbuild-tpie.c \
	  experimental/combparts_tpie.cpp experimental/llist_tpie.cpp \
	  $(OBJS) $(LIBS)

$(BINDIR)/matsolve : matsolve.c $(OBJS)
	$(CC) $(INC) $(CFLAGS) $(LIBFLAGS) -o $@ matsolve.c $(OBJS) $(LIBS)

$(BINDIR)/sqrt : sqrt.c $(OBJS)
	$(CC) $(INC) $(CFLAGS) $(LIBFLAGS) -o $@ sqrt.c $(OBJS) $(LIBS)

$(BINDIR)/polyselect : polyselect.c $(OBJS)
	$(CC) $(INC) $(CFLAGS) $(LIBFLAGS) -o $@ polyselect.c $(OBJS) $(LIBS)

latsiever :
	$(MAKE) -C lasieve4

polsel :
	$(MAKE) -C pol5

strip :
ifeq ($(OS),Windows_NT)
	strip $(BINDIR)/*.exe
else
	strip $(BINS)
endif

clean : lasieve-clean polsel-clean
	-rm -f *.o $(BINDIR)/*.exe core $(BINS)

lasieve-clean :
	$(MAKE) -C lasieve4 clean

polsel-clean :
	$(MAKE) -C pol5 clean

squeaky : clean

