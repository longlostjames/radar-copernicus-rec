# Makefile for Radar Serial Message Package (RSM)
# Chris Walden 20101014

# This is the base path for the package
ROOTPATH=.

# This is where the bits and pieces are kept
SRCDIR = $(ROOTPATH)/src
LIBDIR = $(ROOTPATH)/lib
BINDIR = $(ROOTPATH)/bin
INCDIR = $(ROOTPATH)/include

CC     = gcc
CFLAGS = -Wall -O3 -ffast-math -I$(INCDIR)
LIBS   = -lm -lpthread -lrt

# Top level rule
all : $(LIBDIR)/librsm.a
test: $(BINDIR)/RSM_PositionTest

# The main library
$(LIBDIR)/librsm.a : $(BINDIR)/RSM_SerialMessage.o $(BINDIR)/RSM_SerialPLC.o \
		     $(BINDIR)/RSM_NetworkMessage.o
	ar r $@ $(BINDIR)/RSM_SerialMessage.o $(BINDIR)/RSM_SerialPLC.o \
		$(BINDIR)/RSM_NetworkMessage.o

$(BINDIR)/RSM_SerialMessage.o : $(SRCDIR)/RSM_SerialMessage.c $(INCDIR)/RSM.h
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSM_SerialMessage.c

$(BINDIR)/RSM_SerialPLC.o : $(SRCDIR)/RSM_SerialPLC.c $(INCDIR)/RSM_SerialPLC.h
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSM_SerialPLC.c

$(BINDIR)/RSM_NetworkMessage.o : $(SRCDIR)/RSM_NetworkMessage.c $(INCDIR)/RSM.h
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSM_NetworkMessage.c

clean :
	$(RM) $(BINDIR)/*.[doa]
	$(RM) $(LIBDIR)/*.[doa]

$(BINDIR)/RSM_PositionTest:  $(SRCDIR)/RSM_PositionTest.c $(LIBDIR)/librsm.a
	$(CC) $(CFLAGS) -o $@ $(SRCDIR)/RSM_PositionTest.c \
		-L$(LIBDIR) -lrsm $(LIBS)
