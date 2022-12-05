# Makefile for Radar GPIB ENET Package (REN)
# Owain Davies 04/02/2004

# This is the base path for the package
ROOTPATH=.

# This is where the bits and pieces are kept
SRCDIR = $(ROOTPATH)/src
LIBDIR = $(ROOTPATH)/lib
BINDIR = $(ROOTPATH)/bin
INCDIR = $(ROOTPATH)/include

CC     = gcc
CFLAGS = -Wall -O3 -ffast-math -I$(INCDIR)
LIBS   = -lm

# The master header file
INC = $(INCDIR)/REN.h

# Top level rule
all : $(LIBDIR)/libren.a

# The main library
$(LIBDIR)/libren.a : $(BINDIR)/REN_GpibEnet.o
	ar r $@ $(BINDIR)/REN_GpibEnet.o

$(BINDIR)/REN_GpibEnet.o : $(SRCDIR)/REN_GpibEnet.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/REN_GpibEnet.c

clean :
	$(RM) $(BINDIR)/REN_GpibEnet.o
	$(RM) $(LIBDIR)/libren.a
