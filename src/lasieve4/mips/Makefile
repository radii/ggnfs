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

SRCFILES=gcd32.c 32bit.h psp.c siever-config.w \
	ula2.asm medsched2.asm lasched3.asm pt64.asm lasieve-prepn.asm \
	sieve-from-sched4p.asm tdsieve-from-sched4q.asm ri-aux.asm \
	lasched.h medsched.h slinie.asm ls-defs.asm gcd32.c modinv1002.asm

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

../gmp-aux.h: ../gmp-aux.w
	$(CTANGLE) ../gmp-aux.w

../if.h: ../if.w
	$(CTANGLE) ../if.w

gcd32.o psp.o: siever-config.h ../if.h 

liblasieve.a: ri-aux.o sieve-from-sched4p.o tdsieve-from-sched4q.o \
	modinv1002.o psp.o pt64.o lasieve-prepn.o gcd32.o
	ar rcs $@ $^

ls-defsI%.asm:
	echo i_bits .ASSIGNA $* - 1 > $@

laschedI%.s: ls-defsI%.asm ls-defs.asm ula2.asm lasched3.asm
	gasp -c '#' -o $@ $^

medschedI%.s: ls-defsI%.asm ls-defs.asm ula2.asm medsched2.asm
	gasp -c '#' -o $@ $^

slinieI%.s: ls-defsI%.asm ls-defs.asm slinie.asm
	gasp -c '#' -o $@ $^


liblasieveI%.a:	laschedI%.o medschedI%.o slinieI%.o
	ar rcs $@ $^

bup:
	echo $(SRCFILES)
