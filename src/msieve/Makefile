# ----------------------------------------------------------------------
# Copyright 2007, Jason Papadopoulos
# 
# This file is part of GGNFS.
# 
# GGNFS is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# GGNFS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GGNFS; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# 
# $Id: Makefile,v 1.1 2007-12-16 03:54:46 jasonp_sf Exp $
# ----------------------------------------------------------------------

CC = gcc
WARN_FLAGS = -Wall -W -Wconversion
OPT_FLAGS = -O3 -fomit-frame-pointer -DNDEBUG -D_FILE_OFFSET_BITS=64
MACHINE_FLAGS = -DCPU_GENERIC

ifeq ($(ARCH),pentium2)
	MACHINE_FLAGS = -march=pentium2 -DCPU_PENTIUM2
endif
ifeq ($(ARCH),pentium3)
	MACHINE_FLAGS = -march=pentium3 -DCPU_PENTIUM3
endif
ifeq ($(ARCH),pentium4)
	MACHINE_FLAGS = -march=pentium4 -DCPU_PENTIUM4
endif
ifeq ($(ARCH),prescott)
	MACHINE_FLAGS = -march=prescott -DCPU_CORE
endif
ifeq ($(ARCH),pentium-m)
	MACHINE_FLAGS = -march=pentium-m -DCPU_PENTIUM_M
endif
ifeq ($(ARCH),athlon)
	MACHINE_FLAGS = -march=athlon -DCPU_ATHLON
endif
ifeq ($(ARCH),athlon-xp)
	MACHINE_FLAGS = -march=athlon-xp -DCPU_ATHLON_XP
endif
ifeq ($(ARCH),x86_64)
	MACHINE_FLAGS = -march=k8 -DCPU_OPTERON
endif
ifeq ($(ARCH),nocona)
	MACHINE_FLAGS = -march=nocona -DCPU_PENTIUM4
endif
ifeq ($(ARCH),970)
	MACHINE_FLAGS = -mcpu=970 -mtune=970 -m64 -pipe -DCPU_GENERIC
endif
ifeq ($(ARCH),7450)
	MACHINE_FLAGS = -march=7450 -pipe -DCPU_GENERIC
endif

LIB_THREAD = -lpthread
ifeq ($(OSTYPE),win32)
	LIB_THREAD = 
endif
ifeq ($(OSTYPE),msys)
	LIB_THREAD = 
endif

CFLAGS = $(OPT_FLAGS) $(MACHINE_FLAGS) $(WARN_FLAGS) -Iinclude

#---------------------------------- Generic file lists -------------------

MISC_HDR = \
	include/ap.h \
	include/common.h \
	include/fastmult.h \
	include/msieve.h \
	include/mp.h \
	include/mp_int.h \
	include/util.h

MISC_SRCS = \
	misc/ap.c \
	misc/driver.c \
	misc/expr_eval.c \
	misc/fastmult.c \
	misc/mp.c \
	misc/prime_delta.c \
	misc/prime_sieve.c \
	misc/squfof.c \
	misc/tinyqs.c \
	misc/util.c

MISC_OBJS = $(MISC_SRCS:.c=.o)

#---------------------------------- NFS file lists -------------------------

NFS_HDR = \
	gnfs/filter/filter.h \
	gnfs/filter/merge_util.h \
	gnfs/lanczos/lanczos.h \
	gnfs/sieve/sieve.h \
	gnfs/sqrt/sqrt.h \
	gnfs/gnfs.h

NFS_SRCS = \
	gnfs/filter/clique.c \
	gnfs/filter/duplicate.c \
	gnfs/filter/filter.c \
	gnfs/filter/merge.c \
	gnfs/filter/merge_post.c \
	gnfs/filter/merge_pre.c \
	gnfs/filter/merge_util.c \
	gnfs/filter/singleton.c \
	gnfs/lanczos/gf2.c \
	gnfs/lanczos/lanczos.c \
	gnfs/lanczos/lanczos_matmul0.c \
	gnfs/lanczos/lanczos_matmul1.c \
	gnfs/lanczos/lanczos_matmul2.c \
	gnfs/lanczos/lanczos_pre.c \
	gnfs/sieve/batch_factor.c \
	gnfs/sieve/sieve_line.c \
	gnfs/sieve/sieve_util.c \
	gnfs/sqrt/sqrt.c \
	gnfs/sqrt/sqrt_a.c \
	gnfs/fb.c \
	gnfs/ffpoly.c \
	gnfs/gnfs.c \
	gnfs/relation.c

NFS_OBJS = $(NFS_SRCS:.c=.no)

#---------------------------------- make targets -------------------------

all: $(MISC_OBJS) $(NFS_OBJS)
	ar r libmsieve.a $(MISC_OBJS) $(NFS_OBJS)
	ranlib libmsieve.a
	$(CC) $(CFLAGS) demo.c -o msieve -I../../include \
			libmsieve.a -lm $(LIB_THREAD)

clean:
	rm -f msieve msieve.exe libmsieve.a $(MISC_OBJS) $(NFS_OBJS)

#----------------------------------------- build rules ----------------------

# misc file build rules

%.o: %.c $(MISC_HDR)
	$(CC) $(CFLAGS) -c -o $@ $<

# NFS build rules

%.no: %.c $(MISC_HDR) $(NFS_HDR)
	$(CC) $(CFLAGS) -Ignfs -c -o $@ $<
