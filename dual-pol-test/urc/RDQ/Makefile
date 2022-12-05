# Makefile for Radar Data Acquisition Package (RDQ)
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
INC = $(INCDIR)/RDQ.h

# Top level rule
all : $(LIBDIR)/librdq12.a

# The main library
$(LIBDIR)/librdq12.a : $(BINDIR)/RDQ_DataAcquisition.o
	ar r $@ $(BINDIR)/RDQ_DataAcquisition.o

$(BINDIR)/RDQ_DataAcquisition.o : $(SRCDIR)/RDQ_DataAcquisition.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RDQ_DataAcquisition.c

clean :
	$(RM) $(BINDIR)/*.[doa]
	$(RM) $(LIBDIR)/*.[doa]
