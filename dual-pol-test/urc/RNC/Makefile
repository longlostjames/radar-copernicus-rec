# Makefile for Radar NetCDF Package (RNC)
# Owain Davies 04/02/2004

# This is the base path for the package
ROOTPATH=.

# This is where the bits and pieces are kept
SRCDIR = $(ROOTPATH)/src
LIBDIR = $(ROOTPATH)/lib
BINDIR = $(ROOTPATH)/bin
INCDIR = $(ROOTPATH)/include

CC=gcc
CFLAGS = -Wall -O3 -ffast-math -I$(INCDIR) -I$(ROOTPATH)/../RSP/include \
				-I$(ROOTPATH)/../include
LIBS= -lm

# The master header file
INC = $(INCDIR)/RNC.h

# Top level rule
all : $(LIBDIR)/librnc.a

# The main library
$(LIBDIR)/librnc.a : $(BINDIR)/RNC_NetCDF.o $(BINDIR)/RNC_ReadConfig.o
	ar r $(LIBDIR)/librnc.a $(BINDIR)/RNC_NetCDF.o \
		$(BINDIR)/RNC_ReadConfig.o

$(BINDIR)/RNC_NetCDF.o : $(SRCDIR)/RNC_NetCDF.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RNC_NetCDF.c

$(BINDIR)/RNC_ReadConfig.o : $(SRCDIR)/RNC_ReadConfig.c $(INC)
	$(CC) $(CFLAGS) -o $@ -c $(SRCDIR)/RNC_ReadConfig.c

clean :
	$(RM) $(BINDIR)/*.[doa]
	$(RM) $(LIBDIR)/*.[doa]
