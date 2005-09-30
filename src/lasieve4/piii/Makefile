# Copyright (C) 2001,2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#
# 6/13/04: Hacked up for use in GGNFS by Chris Monico.

INC=-I.. -I../.. -I../../../include $(LOCALINC)

SRCFILES=modinv1002.s gcd32.c ri-aux.s psp.c montgomery_mul.c basemath.c \
         asm_mul64.S asm_mul96.s asm_mul128.s sieve-from-sched.S lasched.c \
         medsched.c Makefile siever-config.c siever-config.h

.SECONDARY: *.s *.S

.PHONY:	all bup clean

all: liblasieveI12.a liblasieveI13.a liblasieveI14.a liblasieve.a

%.o: %.c siever-config.h ../lasieve.h
	$(CC) -c $(CFLAGS) $(INC) $< -o $@

%.o: %.s
	$(CC) $(CFLAGS) $(INC) -c $^

%.o: %.S
	$(CC) $(CFLAGS) $(INC) -c $^

liblasieve.a: modinv1002.o gcd32.o ri-aux.o psp.o montgomery_mul.o \
              basemath.o sieve-from-sched.o siever-config.o \
              asm_mul64.o asm_mul96.o asm_mul128.o
	$(AR) rcs $@ $^

laschedI%.o: lasched.c siever-config.h ../lasieve.h
	$(CC) $(CFLAGS) $(INC) -DI_bits=$* -c -o $@ $<

medschedI%.o: medsched.c siever-config.h ../lasieve.h
	$(CC) $(CFLAGS) $(INC) -DI_bits=$* -c -o $@ $<

liblasieveI%.a: laschedI%.o medschedI%.o
	$(AR) rcs $@ $^

bup:
	echo $(SRCFILES)

clean:
	-rm -f *.o *.a

