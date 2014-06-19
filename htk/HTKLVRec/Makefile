# HTKLVRec/Makefile.  Generated from Makefile.in by configure.

############################################################
# Makefile for HTKLVRec
# 
# N.B. you must make HTKLiblv.a (in the HTKLib directory first)
#
# "make" or "make all": to build all all binaries
# "make install": install all binaries in ${exec_prefix}/bin
# "make clean": remove object files
# "make cleanup": remove object files and binaries
############################################################

SHELL   = 	/bin/sh
inc     = 	../HTKLib
HTKLIB  = 	$(inc)/HTKLib.a 
srcdir = .
top_srcdir = ..

prefix = /usr/local
exec_prefix = ${prefix}
bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var
libdir = ${exec_prefix}/lib
infodir = ${prefix}/share/info
mandir = ${prefix}/share/man
includedir = ${prefix}/include
oldincludedir = /usr/include

CC      = 	g++ -std=c++0x -g
CFLAGS  := 	-DNO_LAT_LM -m32 -D_SVID_SOURCE -DOSS_AUDIO -D'ARCH="x86_64"' -w -I$(inc) -I./ -I../others/
LDFLAGS = 	-L/usr/X11R6/lib  -lm
INSTALL = 	/usr/bin/install -c
HTKLIB = $(inc)/HTKLiblv.a
HEADER = HLVLM.h  HLVModel.h  HLVNet.h  HLVRec.h config.h misc.h

debug: CFLAGS+=-DDEBUG_TRACE
debug: all
all: HDecode.mod


# binaries
HDecode.mod: HDecode.mod.o HLVNet.mod.o HLVRec.mod.o HLVLM.mod.o HLVModel.mod.o misc.mod.o $(HTKLIB)
	$(CC) $(CFLAGS)	-o HDecode.mod $^ $(LDFLAGS)

# .mod version object files
HLVModel.mod.o: HLVModel.c $(HEADER)
	$(CC) -c -DHDECODE_MOD $(CFLAGS) $<
	mv HLVModel.o $@

HLVNet.mod.o: HLVNet.cpp $(HEADER)
	$(CC) -c -DHDECODE_MOD $(CFLAGS) $<
	mv HLVNet.o $@

HLVLM.mod.o: HLVLM.c $(HEADER)
	$(CC) -c -DHDECODE_MOD $(CFLAGS) $<
	mv HLVLM.o $@

HLVRec.mod.o: HLVRec.c HLVRec-GC.c	HLVRec-LM.c	 HLVRec-misc.c	HLVRec-outP.c \
		  HLVRec-propagate.c  HLVRec-traceback.c $(HEADER)
	$(CC) -c -DHDECODE_MOD $(CFLAGS) HLVRec.c 
	mv HLVRec.o $@

HDecode.mod.o: HDecode.cpp $(HEADER)
	$(CC) -c -DHDECODE_MOD $(CFLAGS) $<
	mv HDecode.o $@

misc.mod.o: misc.cpp $(HEADER)
	$(CC) -c -DHDECODE_MOD $(CFLAGS) $<
	mv misc.o $@


# housekeeping rules
strip: HDecode
	-strip HDecode HDecode.mod

clean:
	-rm -f *.o 

cleanup:
	-rm -f *.o HDecode HDecode.mod

distclean:
	-rm -f *.o HDecode HDecode.mod Makefile

install: mkinstalldir
	$(INSTALL) -m 755 HDecode $(bindir)
	$(INSTALL) -m 755 HDecode.mod $(bindir)

mkinstalldir:
	-mkdir -p $(bindir)

.PHONY: std mod all strip clean cleanup distclean install mkinstalldir


