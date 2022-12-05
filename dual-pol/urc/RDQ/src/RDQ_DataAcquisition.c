/*
   Purpose: 	To manage the data acqusition for the radar using the 12 bit data acquisition cards
		( PCICARD, ISACTRL, ISAADC )

   Created on:  02/02/2003
   Created by:  Owain Davies (OTD)
 */


#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/ioctl.h>

/* the files below are custom include files */
/* there location is defined in the make file */
#include <counter.h>
#include <amccdrv.h>
#include <RDQ.h>

void RDQ_InitialiseISACTRL2(int samples_per_pulse, int clock_divfactor, int delay_clocks) {
	int  timercard_fd = 0;        /* file descriptor for ISACTRL card*/
	int  samplecount;             /* number of undivided clocks before we stop sampling a pulse */

	samplecount = delay_clocks + (samples_per_pulse * clock_divfactor);
	timercard_fd = open("/dev/counter", O_RDONLY);
	if (timercard_fd == -1) printf("Could not open ISACTRL\n");
	reset_counter(timercard_fd);
	set_counter(timercard_fd, COUNTER_COUNT0, delay_clocks - 1, COUNTER_CONTROL_MODE_FIVE);
	set_counter(timercard_fd, COUNTER_COUNT1, clock_divfactor, COUNTER_CONTROL_MODE_TWO);
	set_counter(timercard_fd, COUNTER_COUNT2, samplecount, COUNTER_CONTROL_MODE_FIVE);
	close(timercard_fd);
}

int RDQ_InitialisePCICARD2(void ** dma_buffer , size_t dma_buffer_size) {
	/*
        IN:   dma_banks pointers to the start of the two dma banks
        OUT:  the file descriptor for the pci card
	 */

	int     amcc_fd = 0;     /* file descriptor for the PCICARD */

	/* open the file descriptor to the PCICARD */
	amcc_fd = open("/dev/amcc", O_RDWR);
	if( amcc_fd == -1 ) printf("Could not open PCICARD\n");

	/* allocate the DMA buffer for the PCICARD */
	if( (ioctl(amcc_fd, AMCC_IO_ALLOC_DMA_BUFFER, dma_buffer_size )) == -1) {
		printf("Could not allocate DMA buffer\n");
	}
	*dma_buffer = mmap ( NULL, dma_buffer_size, (PROT_READ | PROT_WRITE),
			MAP_SHARED, amcc_fd, AMCC_MMAP_DMA_BUFFER);
	if(*dma_buffer == MAP_FAILED) printf("Could not mmap DMA buffer\n");
	//dma_banks[ 0 ] = (short int*)  *dma_buffer;
	//dma_banks[ 1 ] = (short int*) (*dma_buffer + (dma_buffer_size/2));
	return (amcc_fd);
}


void RDQ_ClosePCICARD2 (int amcc_fd, void * dma_buffer , size_t dma_buffer_size) {
	/*
        IN:  dma_bank tells us which dma_bank we are going to write to
	 */

	/* unmap DMA buffer */
	if(munmap (dma_buffer, dma_buffer_size) == -1) printf("Could not munmap\n");
	close (amcc_fd);
}



void RDQ_StartAcquisition2 (int amcc_fd, int dma_bank, int tcount) {
	/*
	IN:  dma tells us which dma_bank we are going to write to
	 */

	/* sets the active dma_bank */
	if(ioctl (amcc_fd, AMCC_IO_SET_DMA_BANK, dma_bank)) {
		printf("FATAL ERROR: could not set active DMA bank\n");
	}

	/* start the DMA */
	if(ioctl (amcc_fd, AMCC_IO_START_DMA, tcount)) {
		printf("FATAL ERROR: could not start DMA\n");
	}
}

int RDQ_WaitForAcquisitionToComplete(int amcc_fd) {
	/*
	IN:     amcc_fd
	OUT:
	RETURN: the status value from AMCC_IO_WAIT_FOR_DMA
	 */

	int ret_val = 0;

	ioctl(amcc_fd, AMCC_IO_WAIT_FOR_DMA, &ret_val);
	switch(ret_val) {
	case 1:
		printf("FATAL ERROR: DMA timed out\n");
		break;
	case 2:
		printf("FATAL ERROR: Master abort occured\n");
		break;
	case 3:
		printf("FATAL ERROR: FIFO overflow occured\n");
		break;
	case 4:
		printf("FATAL ERROR: PCICARD did not trigger (HV timeout)\n");
		break;
	}
	return ret_val;
}
