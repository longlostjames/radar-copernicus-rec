# Makefile for Radar NetCDF Package (RNC)
# Owain Davies 04/02/2004

# This is the base path for the package
ROOTPATH=.

# This is where the bits and pieces are kept
SRCDIR = $(ROOTPATH)/src
LIBDIR = $(ROOTPATH)/lib
BINDIR = $(ROOTPATH)/bin
INCDIR = $(ROOTPATH)/include

CC     = gcc
CFLAGS = -Wall -O3 -ffast-math -I$(INCDIR) -I../RSP/include
LIBS   = -lm

# The master header file
INC = $(INCDIR)/RLC.h

# Top level rule
all : $(LIBDIR)/librlc.a

# The main library
$(LIBDIR)/librlc.a : $(BINDIR)/RLC_ProgramLecroy.o
	ar r $@ $(BINDIR)/RLC_ProgramLecroy.o

$(BINDIR)/RLC_ProgramLecroy.o : $(SRCDIR)/RLC_ProgramLecroy.c $(INC)
	$(CC) $(CFLAGS) -c $(SRCDIR)/RLC_ProgramLecroy.c -o $(BINDIR)/RLC_ProgramLecroy.o

clean :
	$(RM) $(BINDIR)/RLC_ProgramLecroy.o
	$(RM) $(LIBDIR)/librlc.a
