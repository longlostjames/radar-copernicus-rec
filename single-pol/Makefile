# make file for radar-copernicus-rec.c
# Owain Davies 22/5/2003
# Ed Pavelin 18/5/2004

# Explanation of the flags
# -Wall prints lots of extra warning messages
# -O3 does as much optimising as it can (trade off, longer compilation time)
# -ffast-math breaks some rules in the quest for execution speed
# -ggdb make the code debuggable
# -lm includes the maths library
CC     = gcc
CFLAGS = -Wall -O3 -ffast-math
LIBS   = -lnetcdf -lfftw3 -lm -lpthread -lrt

URC_PATH = urc
PATH_RSP = $(URC_PATH)/RSP/lib
PATH_RDQ = $(URC_PATH)/RDQ/lib
PATH_RNC = $(URC_PATH)/RNC/lib
PATH_RSM = $(URC_PATH)/RSM/lib
PATH_REL = $(URC_PATH)/REL/lib

INSTALL_DIR = /usr/local/bin

override CFLAGS += -I$(URC_PATH)/RDQ/include
override CFLAGS += -I$(URC_PATH)/RSP/include
override CFLAGS += -I$(URC_PATH)/RNC/include
override CFLAGS += -I$(URC_PATH)/RSM/include
override CFLAGS += -I$(URC_PATH)/REL/include
override CFLAGS += -I$(URC_PATH)/include

EXE = radar-copernicus-rec

# Universal radar code libraries
URC_LIBS = $(PATH_RDQ)/librdq12.a $(PATH_RSM)/librsm.a \
	   $(PATH_RSP)/librsp.a $(PATH_RNC)/librnc.a $(PATH_REL)/librel.a
LDFLAGS  = -L $(PATH_RDQ) -L $(PATH_RSM) -L $(PATH_RSP) -L $(PATH_RNC) \
	   -L $(PATH_REL) -lrdq12 -lrsm -lrsp -lrnc -lrel

all: copernicus

.PHONY: all clean help copernicus install

help:
	@echo
	@echo "make copernicus  Make Copernicus runtime data acquisition program"
	@echo
	@echo "make install     Install Copernicus runtime data acquisition program"
	@echo
	@echo "make clean       Cleanup build"
	@echo "make help        This help"
	@echo

copernicus : $(EXE)

$(EXE) : radar-copernicus-rec.o median.o $(URC_LIBS)
	$(CC) $(CFLAGS) -o $@ radar-copernicus-rec.o median.o \
		$(LDFLAGS) $(LIBS)

median.o : median.c median.h
	$(CC) $(CFLAGS) -c median.c

radar-copernicus-rec.o : radar-copernicus-rec.c radar-copernicus-rec.h
	$(CC) $(CFLAGS) -c radar-copernicus-rec.c

$(PATH_RDQ)/librdq12.a:
	$(MAKE) -C $(URC_PATH)/RDQ

$(PATH_RSM)/librsm.a:
	$(MAKE) -C $(URC_PATH)/RSM

$(PATH_RSP)/librsp.a:
	$(MAKE) -C $(URC_PATH)/RSP

$(PATH_RNC)/librnc.a:
	$(MAKE) -C $(URC_PATH)/RNC

$(PATH_REL)/librel.a:
	$(MAKE) -C $(URC_PATH)/REL

clean : 
	$(RM) *.[doa] $(EXE)
	$(MAKE) -C $(URC_PATH)/RDQ $@
	$(MAKE) -C $(URC_PATH)/RSM $@
	$(MAKE) -C $(URC_PATH)/RSP $@
	$(MAKE) -C $(URC_PATH)/RNC $@
	$(MAKE) -C $(URC_PATH)/REL $@
	@find . $(URC_PATH) -name "*~" -type f -print -exec rm \{\} \;

install : 
	echo "Installing $(EXE)"
	echo $(INSTALL_DIR)
	install -d $(INSTALL_DIR)
	install -o jif -g radar_dt -m 6550 $(EXE)                       $(INSTALL_DIR)

	echo "Installing radar-acrobat-rec"
	echo $(INSTALL_DIR) 
	cp $(EXE) $(INSTALL_DIR).
	install -o jif -g radar_dt -m 6550 run_copernicus_data_acquisition.pl $(INSTALL_DIR)
	install -o jif -g radar_dt -m 6550 start_radar-copernicus.sh          $(INSTALL_DIR)
	install -o jif -g radar_dt -m 6550 stop_radar-copernicus.sh           $(INSTALL_DIR)
