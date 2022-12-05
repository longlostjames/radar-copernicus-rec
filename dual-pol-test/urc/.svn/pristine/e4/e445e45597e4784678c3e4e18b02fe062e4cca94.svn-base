# Makefile for Radar Data Acquisition Package (RDQ)
# Owain Davies 04/02/2004

# This is the base path for the package
ROOTPATH=.

# This is where the bits and pieces are kept
SRCDIR = $(ROOTPATH)/src
LIBDIR = $(ROOTPATH)/lib
BINDIR = $(ROOTPATH)/bin
INCDIR = $(ROOTPATH)/include

DSAKDIR = /root/d2k-dask_175
CONDIR  = $(DSAKDIR)/samples/conio

PROG    = RAT_AdlinkDataAcquisition

CC     = gcc
CFLAGS = -Wall -O3 -ffast-math -I$(INCDIR) -I$(DSAKDIR)/include
LIBS   = -L$(DSAKDIR)/lib -lpci_dask2k -lm -lpthread -lrt

# The master header file
INC = $(INCDIR)/RAT.h

# Top level rule
all : $(LIBDIR)/librat.a
diag: $(BINDIR)/RAT_diag

$(BINDIR)/RAT_diag : $(BINDIR)/RAT_diag.o $(BINDIR)/conio.o $(LIBDIR)/librat.a
	$(CC) $(CFLAGS) -o $@ $(BINDIR)/RAT_diag.o $(BINDIR)/conio.o -L $(LIBDIR) -lrat $(LIBS)

$(BINDIR)/RAT_diag.o : $(SRCDIR)/RAT_diag.c
	$(CC) $(CFLAGS) -o $(BINDIR)/RAT_diag.o -c $(SRCDIR)/RAT_diag.c

# The main library
$(LIBDIR)/librat.a : $(BINDIR)/$(PROG).o
	ar r $@ $(BINDIR)/$(PROG).o

$(BINDIR)/$(PROG).o : $(SRCDIR)/$(PROG).c
	$(CC) $(CFLAGS) -c $(SRCDIR)/$(PROG).c -o $(BINDIR)/$(PROG).o

$(BINDIR)/conio.o : $(CONDIR)/conio.c
	$(CC) $(CFLAGS) -o $(BINDIR)/conio.o -c $(CONDIR)/conio.c

clean :
	$(RM) $(BINDIR)/*.[doa]
	$(RM) $(LIBDIR)/*.[doa]
	$(RM) $(BINDIR)/RAT_diag
