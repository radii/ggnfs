# Copyright (C) 2001,2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the 
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

CFLAGS= -O3

CC=gcc
CTANGLE=ctangle

SRCFILES=gcd32.c 32bit.h modinv32.c psp.c lasched.w medsched.w siever-config.w

%.c: %.w %.ch
	$(CTANGLE) $*.w $*.ch

%.c: %.w
	$(CTANGLE) $*.w

%.h: %.w %.ch
	$(CTANGLE) $*.w $*.ch

%.h: %.w
	$(CTANGLE) $*.w

%.tex: %.w %.ch
	cweave $*.w $*.ch

%.tex: %.w
	cweave $*.w

%.dvi: %.tex
	tex $<

%.o: %.s
	cc -c $^

.SECONDARY:  *.s

.PHONY: bup

if.h:	../if.h
	cp -f $^ $@

gmp-aux.h: ../gmp-aux.h
	cp -f $^ $@

liblasieve.a: gcd32.o psp.o modinv32.o
	ar rcs $@ $^

laschedI%.o: lasched.c siever-config.h if.h
	$(CC) $(CFLAGS) -DI_bits=$* -g -c -o $@ $<

medschedI%.o: medsched.c siever-config.h if.h
	$(CC) $(CFLAGS) -DI_bits=$* -g -c -o $@ $<

liblasieveI%.a:	laschedI%.o medschedI%.o # slinieI%.o
	ar rcs $@ $^

bup:
	echo $(SRCFILES)
