CC=gcc

ALLOPT=-march=$(ARCH) -mtune=$(ARCH)

ifeq ($(ARCH),k8)
    ALLOPT+= -D__ppc__
# -D__ppc__ should be removed later
# TODO: carefully inspect all ``#if defined(__ppc__)'' and
# add ``#if defined(__x86_64__)'' where required
endif

DEBUGOPT=-W -Wall -g -pg

CFLAGS=-O3 -fomit-frame-pointer -W -Wall $(ALLOPT)
#CFLAGS=$(DEBUGOPT) $(ALLOPT)

export CC ALLOPT CFLAGS DEBUGOPT

INC=-I. -I.. -I../include
LIBS=-lgmp -lm
BINDIR=../bin
LANCZOS=blanczos64-no-mmx

OBJS=getprimes.o modinv32.o \
     fbmisc.o squfof.o rels.o \
     $(LANCZOS).o poly.o mpz_poly.o \
     mpz_mat.o smintfact.o misc.o ecm4c.o nfmisc.o \
     montgomery_sqrt.o matstuff.o dickman.o makefb.o fbgen.o\
     llist.o combparts.o if.o
BINS=$(BINDIR)/sieve $(BINDIR)/procrels $(BINDIR)/sqrt $(BINDIR)/polyselect \
     $(BINDIR)/makefb $(BINDIR)/matsolve $(BINDIR)/matbuild
LSBINS=latsiever polsel
TESTBIN=./tsrc
TESTCODE=./tsrc

.PHONY: all tests bins latsiever polsel clean lasieve-clean polsel-clean \
        squeaky

all : $(OBJS) $(BINS) $(LSBINS)
tests : $(TESTS)
bins : $(BINS)

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $*.c

smintfact.o : smintfact.c
	$(CC) $(CFLAGS) $(INC) -c smintfact.c -o $@

fbmisc.o : fbmisc.c
	$(CC) $(CFLAGS) $(INC) -c fbmisc.c -o $@

rels.o : rels.c
	$(CC) $(CFLAGS) $(INC) -c rels.c -o $@

latsieve.o : latsieve.c
	$(CC) $(CFLAGS) $(INC) -c latsieve.c -o $@

$(BINDIR)/sieve : sieve.c clsieve.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) -o $@ sieve.c clsieve.c $(OBJS) $(LIBS)

$(BINDIR)/makefb : makefb.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) -D_MAKEFB_STANDALONE -o $@ makefb.c getprimes.o \
             fbmisc.o poly.o mpz_poly.o smintfact.o modinv32.c if.o \
             misc.o nfmisc.o mpz_mat.o ecm4c.o fbgen.o rels.o $(LIBS)

$(BINDIR)/procrels : procrels.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) -o $@ procrels.c $(OBJS) $(LIBS)

$(BINDIR)/matsolve : matsolve.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) -o $@ matsolve.c $(OBJS) $(LIBS)

$(BINDIR)/matbuild : matbuild.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) -o $@ matbuild.c $(OBJS) $(LIBS)

$(BINDIR)/sqrt : sqrt.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) -o $@ sqrt.c $(OBJS) $(LIBS)

$(BINDIR)/polyselect : polyselect.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) -o $@ polyselect.c $(OBJS) $(LIBS)

latsiever :
	$(MAKE) -C lasieve4

polsel :
	$(MAKE) -C pol5

clean : lasieve-clean polsel-clean
	-rm -f *.o core $(BINS)

lasieve-clean :
	$(MAKE) -C lasieve4 clean

polsel-clean :
	$(MAKE) -C pol5 clean

squeaky : clean

