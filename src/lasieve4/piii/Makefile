# Copyright (C) 2001,2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#
# 6/13/04: Hacked up for use in GGNFS by Chris Monico.

INC=-I.. -I../.. -I../../../include

SRCFILES=modinv1002.s gcd32.c ri-aux.s psp.c montgomery_mul.c basemath.c \
         asm_mul64.S asm_mul96.s asm_mul128.s sieve-from-sched.S lasched.c \
         medsched.c Makefile siever-config.c lasieve-asm.h

.SECONDARY: *.s

.PHONY:	all bup clean

all: liblasieveI12.a liblasieveI13.a liblasieveI14.a liblasieve.a

%.o: %.c
	$(CC) -c $(CFLAGS) $(INC) $< -o $@

%.o: %.s
	$(CC) $(CFLAGS) $(INC) -c $^

montgomery_mul.o: lasieve-asm.h ../lasieve.h

liblasieve.a: modinv1002.o gcd32.o ri-aux.o psp.o montgomery_mul.o \
              basemath.o sieve-from-sched.o siever-config.o \
              asm_mul64.o asm_mul96.o asm_mul128.o
	ar rcs $@ $^

laschedI12.o: lasched.c lasieve-asm.h ../lasieve.h
	$(CC) $(CFLAGS) $(INC) -DI_bits=12 -c -o $@ $<

laschedI13.o: lasched.c lasieve-asm.h ../lasieve.h
	$(CC) $(CFLAGS) $(INC) -DI_bits=13 -c -o $@ $<

laschedI14.o: lasched.c lasieve-asm.h ../lasieve.h
	$(CC) $(CFLAGS) $(INC) -DI_bits=14 -c -o $@ $<

medschedI12.o: medsched.c lasieve-asm.h ../lasieve.h
	$(CC) $(CFLAGS) $(INC) -DI_bits=12 -c -o $@ $<

medschedI13.o: medsched.c lasieve-asm.h ../lasieve.h
	$(CC) $(CFLAGS) $(INC) -DI_bits=13 -c -o $@ $<

medschedI14.o: medsched.c lasieve-asm.h ../lasieve.h
	$(CC) $(CFLAGS) $(INC) -DI_bits=14 -c -o $@ $<

liblasieveI12.a: laschedI12.o medschedI12.o
	ar rcs $@ $^

liblasieveI13.a: laschedI13.o medschedI13.o
	ar rcs $@ $^

liblasieveI14.a: laschedI14.o medschedI14.o
	ar rcs $@ $^

bup:
	echo $(SRCFILES)

clean:
	-rm -f *.o *.a

