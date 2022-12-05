/* this file will contain routines to access and acquire data from the */
/* adlink data acquisition card 2010 */
/* Owain Davies 20070821 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>
#include <unistd.h>
#include <sys/mman.h>
#include <d2kdask.h>
#include <RAT.h>

void
RAT_CheckError (I16    card,
		I16    err,
		char * buffer)
{
    if (err)
    {
	printf ("%s: %d", buffer, err);
	D2K_Release_Card (card);
	exit (1);
    }
}

int RAT_CollectAcquisition (I16 card)
{
    BOOLEAN fStop;
    U32     count;
    U32     startPosition;

    fStop = 0;
    while (!fStop)
    {
        D2K_AI_AsyncCheck (card, &fStop, &count);
    }
    D2K_AI_AsyncClear (card, &startPosition, &count);

    D2K_Release_Card (card);

    return (int)count;
}


I16
RAT_StartAcquisition (int   pulses,
		      int   samples,
		      I16 * data_buffer)
{

    I16 card_number  = 0;
    I16 card;
    I16 err;
    I32 sample_intrval = 8;
    U16 channels    = 4;
    U16 channel[4]  = { 0, 1, 2, 3 };
    U16 range       = AD_U_5_V;
    U16 Id;
    int i;

    /* open access to card */
    card = D2K_Register_Card (DAQ_2010, card_number);
    if (card < 0)
    {
        printf ("Register_Card error: %d\n", card);
        exit (1);
    }

    /* configure the ranges on the analogue inputs */
    for (i = 0; i < channels; i++)
    {
        err = D2K_AI_CH_Config (card, channel[i], range);
	RAT_CheckError (card, err, "D2K_AI_CH_Config error");
    }

    /* select external timebase, but no analogue triggering */
    err = D2K_AIO_Config (card, DAQ2K_ExtTimeBase, 0, 0x80, 0x80);
    RAT_CheckError (card, err, "D2K_AIO_Config error");

    /* select source for ADC conversion and trigger, and retrigger count */
    err = D2K_AI_Config (card,
			 DAQ2K_AI_ADCONVSRC_Int,  /* ADC conversion source */
			 DAQ2K_AI_TRGSRC_ExtD | DAQ2K_AI_TRGMOD_POST |
			 DAQ2K_AI_TrgNegative | DAQ2K_AI_ReTrigEn,
			 0,
			 0,
			 pulses, /* retrigger count */
			 1);
    RAT_CheckError (card, err, "D2K_AI_Config error");

    /* set up the SSI signals */
    err = D2K_SSI_SourceConn (card, SSI_TIME);
    if (err)
    {
	printf ("D2K_SSI_SourceConn error=%d", err);
        D2K_Release_Card (card);
        return -2;
    }

    err = D2K_SSI_SourceConn (card, SSI_CONV);
    if (err)
    {
        printf ("D2K_SSI_SourceConn error=%d", err);
        D2K_Release_Card (card);
        return -2;
    }

    /* switch off double buffer mode */
    err = D2K_AI_AsyncDblBufferMode (card, 0);
    RAT_CheckError (card, err, "D2K_AI_AsyncDblBufferMode error");

    /* set up one of the double buffers */
    err = D2K_AI_ContBufferSetup (card, data_buffer,
				  samples * pulses * channels, &Id);
    RAT_CheckError (card, err, "AI_ContBufferSetup error");

    /* start the acqusiiton off */
    err = D2K_AI_ContReadMultiChannels (card, channels, channel, Id,
					samples * pulses, sample_intrval,
					sample_intrval, ASYNCH_OP);
    RAT_CheckError (card, err, "AI_ContReadChannel error");

    return card;
}
