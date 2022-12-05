/* this program will use the RAT routines */
/* to do some very simple data acquisition */

#include<stdio.h>
#include<d2kdask.h>
#include<RAT.h>

int
main (void)
{

    I16 data_buffer[1000000];
    I16 card;
    int pulses;
    int samples;
    int channels = 4;
    int count;
    int i;

    pulses  = 256;
    samples = 100;

    card  = RAT_StartAcquisition (pulses, samples, data_buffer);
    count = RAT_CollectAcquisition (card);

    printf ("%d\n", pulses);
    printf ("%d\n", samples);
    printf ("%d\n", channels);
    for (i = 0; i < pulses * samples * 4; i += 4)
    {
	printf ("%6d, %6d, %6d, %6d\n",
		data_buffer[i+0],
		data_buffer[i+1],
		data_buffer[i+2],
		data_buffer[i+3]);
    }

    /* printf ("acquisition complete, count: %d\n", count); */

    return 0;
}
