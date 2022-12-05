#ifndef _RDQ_H
#define _RDQ_H

#include <sys/types.h>
#include <stdint.h>
#include <string.h>

#define DMA_BUFFER_SIZE (32*1024*1024*4)

/* Macros to manipulate DMA bank data extraction pointers */
#define GET_POINTER(buffer, num_samples, num_channels, pulse, sample) ((buffer) + (num_channels) * (((num_samples) * (pulse)) + (sample)))
#define INC_POINTER(buffer, num_channels) buffer += (num_channels)

/* Macro to perform de-multiplexing of channels */
#define GET_CHANNEL(buffer, channel) ((buffer)[dmux_table[(channel)]])

void RDQ_InitialiseISACTRL2 (int samples_per_pulse, int clock_divfactor, int delayclocks);
int  RDQ_InitialisePCICARD2 (void ** dma_buffer, size_t dma_buffer_size);
void RDQ_ClosePCICARD2      (int amcc_fd, void * dma_buffer, size_t dma_buffer_size);
void RDQ_StartAcquisition2  (int amcc_fd, int dma_bank, int tcount);
int  RDQ_WaitForAcquisitionToComplete (int amcc_fd);

static inline void
RDQ_InitialiseISACTRL (int pulses_per_daq_cycle __attribute__ ((__unused__)),
		       int samples_per_pulse,
		       int clock_divfactor,
		       int delayclocks)
{
    RDQ_InitialiseISACTRL2 (samples_per_pulse, clock_divfactor, delayclocks);
}

static inline int
RDQ_InitialisePCICARD_New (caddr_t * dma_buffer,
			   size_t    dma_buffer_size)
{
    return RDQ_InitialisePCICARD2 ((void **)dma_buffer, dma_buffer_size);
}

static inline void
RDQ_ClosePCICARD_New (int       amcc_fd,
		      caddr_t * dma_buffer,
		      size_t    dma_buffer_size)
{
    RDQ_ClosePCICARD2 (amcc_fd, *dma_buffer, dma_buffer_size);
}

static inline int
RDQ_InitialisePCICARD (caddr_t *dma_buffer)
{
    return RDQ_InitialisePCICARD_New (dma_buffer, DMA_BUFFER_SIZE);
}

static inline void
RDQ_ClosePCICARD (int       amcc_fd,
		  caddr_t * dma_buffer)
{
    RDQ_ClosePCICARD_New (amcc_fd, dma_buffer, DMA_BUFFER_SIZE);
}

static inline void
RDQ_StartAcquisition  (int         amcc_fd,
		       int         dma_bank,
		       short int * data,
		       int         tcount)
{
    memset (data, -1, tcount);

    RDQ_StartAcquisition2 (amcc_fd, dma_bank, tcount);
}

#endif /* !_RDQ_H */
