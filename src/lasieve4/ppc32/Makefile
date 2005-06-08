# Copyright (C) 2001,2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the 
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

CFLAGS= -O3

CC=gcc

SRCFILES=gcd32.c 32bit.h modinv32.c mpqs.c psp.c Makefile \
	siever-config32.h lasched.c medsched.c \

all: liblasieve.a liblasieveI12.a liblasieveI13.a liblasieveI14.a

liblasieve.a: gcd32.o psp.o mpqs.o modinv32.o
	ar rcs $@ $^

laschedI%.o: lasched.c siever-config.h
	$(CC) $(CFLAGS) -DI_bits=$* -c -o $@ $<

medschedI%.o: medsched.c siever-config.h
	$(CC) $(CFLAGS) -DI_bits=$* -c -o $@ $<

liblasieveI%.a:	laschedI%.o medschedI%.o
	ar rcs $@ $^

bup:
	echo $(SRCFILES)
