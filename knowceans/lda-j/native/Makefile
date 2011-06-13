###############################################################################
# (C) Copyright 2005, Gregor Heinrich (gregor :: arbylon : net) (This file is
# part of the lda-j (org.knowceans.lda.*) experimental software package.)
#
# lda-j is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# lda-j is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA 02111-1307 USA
###############################################################################

.SUFFIXES: .c .u
CC= gcc
#CFLAGS= -O3
CFLAGS= -g -mno-cygwin -I ${JAVA_HOME}/include -I ${JAVA_HOME}/include/win32 -shared -Wl,--add-stdcall-alias 
LDFLAGS= -lm -mno-cygwin

LOBJECTS= cokus.o

LSOURCE= cokus.c

all: $(LOBJECTS)
	$(CC) $(CFLAGS) $(LOBJECTS) -o ../mtcokus.dll $(LDFLAGS)

clean:
	-rm -f *.o
