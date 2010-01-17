# Copyright (C) 2001,2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#
# 6/13/04: Hacked up for use in GGNFS by Chris Monico.

INC=-I.. -I../../include $(LOCALINC)
LIBFLAGS=-L. $(LOCALLIB)
CFLAGS+=-Os

BINDIR=../../bin
FACT=prho.o
.SUFFIXES:

.SECONDARY: *.c *.o

.PHONY: all clean

SRCFILES=fbgen.c gnfs-lasieve4e.c input-poly.c mpqs.c mpz-ull.c \
         real-poly-aux.c redu2.c gmp-aux.c if.c lasieve-prepn.c \
	 primgen32.c recurrence6.c lasieve.h asm/siever-config.h

OBJS=../if.o input-poly.o redu2.o recurrence6.o ../fbgen.o \
     real-poly-aux.o primgen32.o lasieve-prepn.o mpqs.o

LIBS=-lgmp-aux -lgmp -lm

ASMDIRS=piii ppc32 itanium generic mips

list_asm_files = \
     $(foreach file, $(shell make -s -C $(asm_dir) bup),$(asm_dir)/$(file))

ASMFILES=$(foreach asm_dir,$(ASMDIRS),$(list_asm_files))

all : $(BINDIR)/gnfs-lasieve4I12e $(BINDIR)/gnfs-lasieve4I13e \
      $(BINDIR)/gnfs-lasieve4I14e $(BINDIR)/gnfs-lasieve4I15e \
      $(BINDIR)/gnfs-lasieve4I16e $(BINDIR)/gnfs-lasieve4I11e

asm/%:
	$(MAKE) -C asm $*

%.o: %.c lasieve.h asm/siever-config.h
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

%.o: %.s
	$(CC) $(INC) -c $^

libgmp-aux.a: gmp-aux.o mpz-ull.o
	$(AR) rcs $@ $^

gnfs-lasieve4eI%.o: gnfs-lasieve4e.c lasieve.h asm/siever-config.h
	$(CC) $(CFLAGS) $(INC) -c -DI_bits=$* -o $@ $<

$(BINDIR)/gnfs-lasieve4I%e: gnfs-lasieve4eI%.o $(OBJS) libgmp-aux.a \
                             asm/liblasieve.a asm/liblasieveI%.a $(FACT)
	$(CC) $(CFLAGS) $(INC) $(LIBFLAGS) -o $@ $^ $(LIBS)

asm/lib%.a:
	$(MAKE) -C asm

clean:
ifeq ($(OS),Windows_NT)
	-rm -f *.o *.a $(BINDIR)/gnfs-lasieve4I1?e.exe
else
	-rm -f *.o *.a $(BINDIR)/gnfs-lasieve4I1?e
endif
	(test -d asm && $(MAKE) -C asm clean) || exit 0
