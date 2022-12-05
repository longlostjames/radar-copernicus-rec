# Makefile for Radar Signal Processing Package (RSP)
# Ed Pavelin 05/11/03

# This is the base path for the package
ROOTPATH=.

# This is where the bits and pieces are kept
SRCDIR = $(ROOTPATH)/src
LIBDIR = $(ROOTPATH)/lib
BINDIR = $(ROOTPATH)/bin
INCDIR = $(ROOTPATH)/include

CC     = gcc
CFLAGS = -Wall -O3 -ffast-math -I$(INCDIR)
LIBS   = -lfftw3 -lm

# The master header file
INC = $(INCDIR)/RSP.h

# Top level rule
all : $(LIBDIR)/librsp.a

# The main library
$(LIBDIR)/librsp.a : $(BINDIR)/RSP_CalcSpecMom.o $(BINDIR)/RSP_FindPeaks.o \
	$(BINDIR)/RSP_CalcPSD.o $(BINDIR)/RSP_Initialise.o \
	$(BINDIR)/RSP_Correlate.o $(BINDIR)/RSP_ClutterInterp.o \
	$(BINDIR)/RSP_FreeMemory.o $(BINDIR)/RSP_CalcPhase.o \
	$(BINDIR)/RSP_Observables.o $(BINDIR)/RSP_DisplayParams.o
	ar r $@ $(BINDIR)/RSP_CalcSpecMom.o \
		$(BINDIR)/RSP_FindPeaks.o $(BINDIR)/RSP_CalcPSD.o \
		$(BINDIR)/RSP_Initialise.o $(BINDIR)/RSP_Correlate.o \
		$(BINDIR)/RSP_ClutterInterp.o $(BINDIR)/RSP_FreeMemory.o \
		$(BINDIR)/RSP_CalcPhase.o $(BINDIR)/RSP_Observables.o \
		$(BINDIR)/RSP_DisplayParams.o

$(BINDIR)/RSP_DisplayParams.o : $(SRCDIR)/RSP_DisplayParams.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSP_DisplayParams.c

$(BINDIR)/RSP_Observables.o : $(SRCDIR)/RSP_Observables.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSP_Observables.c

$(BINDIR)/RSP_CalcPhase.o : $(SRCDIR)/RSP_CalcPhase.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSP_CalcPhase.c

$(BINDIR)/RSP_Initialise.o : $(SRCDIR)/RSP_Initialise.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSP_Initialise.c

$(BINDIR)/RSP_FreeMemory.o : $(SRCDIR)/RSP_FreeMemory.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSP_FreeMemory.c

$(BINDIR)/RSP_CalcSpecMom.o : $(SRCDIR)/RSP_CalcSpecMom.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSP_CalcSpecMom.c

$(BINDIR)/RSP_FindPeaks.o : $(SRCDIR)/RSP_FindPeaks.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSP_FindPeaks.c

$(BINDIR)/RSP_CalcPSD.o : $(SRCDIR)/RSP_CalcPSD.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSP_CalcPSD.c

$(BINDIR)/RSP_Correlate.o : $(SRCDIR)/RSP_Correlate.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSP_Correlate.c

$(BINDIR)/RSP_ClutterInterp.o : $(SRCDIR)/RSP_ClutterInterp.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RSP_ClutterInterp.c

clean :
	$(RM) $(BINDIR)/*.[doa]
	$(RM) $(LIBDIR)/*.[doa]
