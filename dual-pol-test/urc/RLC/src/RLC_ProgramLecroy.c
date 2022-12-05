/*
 * RLC_SerialMessage.c
 * This code allows the lecroy to be programmed
 * created by : Owain Davies
 * created on : 08/09/2004
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <termios.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/ugpib.h>

#include <RLC.h>

int
RLC_GPIBCommunication (int          pad,
		       const char * msg,
		       char *       reply,
		       bool         ctrl)
{
    Addr4882_t instrument;
    int        gpib_board = 0;

    instrument = MakeAddr (pad, 0);
    DevClear (gpib_board, instrument);
    if (strcmp (msg, ":SYST:ERR?"))
    {
	printf ("%s\n", msg);
    }
    Send (gpib_board, instrument, msg, strlen (msg), NLend);


    if (ctrl)
    {   /* a response is expected */
	Receive (gpib_board, instrument, reply, 255, STOPend);
	if  (ibsta & ERR)
	{
	    printf ("An error has occured while talking to the Lecroy\n");
	    printf ("This could be due to a bus extender being on\n");
	    printf ("The Lecroy LW120 does not like the HP bus extender being on\n");
	    exit (1);
	}
	reply[ibcnt] = '\0';
	//printf ("%s", reply);
    }
    usleep (100000);
    return 1;
}

int
RLC_LecroyErrorCheck (int pad)
{
    char reply[255];

    RLC_GPIBCommunication (pad, ":SYST:ERR?", reply, true);
    if (reply[0] != '0')
    {
	printf ("Lecroy says: %s", reply);
	exit (1);
    }
    return 1;
}

int
RLC_LecroyReset (int pad)
{
    char reply[255];
    int  i;

    RLC_GPIBCommunication (pad, "*RST", reply, false);
    sleep (1);

    RLC_GPIBCommunication (pad, "*CLS", reply, false);
    RLC_LecroyErrorCheck (pad);

    RLC_GPIBCommunication (pad, "*IDN?", reply, true);
    RLC_LecroyErrorCheck (pad);

    for  (i = 1; i <= 2; i++)
    {
	RLC_LecroySelectChannel (pad, i);
	RLC_GPIBCommunication (pad, ":SEQ:DEL:ALL",  reply, false);
	RLC_GPIBCommunication (pad, ":TRAC:DEL:ALL", reply, false);
	RLC_LecroyErrorCheck (pad);
    }
    return 1;
}

int
RLC_LecroySelectChannel (int pad,
			 int channel)
{
    char command[255];
    char reply[255];

    sprintf (command, ":INST:SEL %d", channel);
    RLC_GPIBCommunication (pad, command, reply, false);
    usleep (100000);
    RLC_LecroyErrorCheck (pad);

    RLC_GPIBCommunication (pad, ":INST?", reply, true);
    RLC_LecroyErrorCheck (pad);

    return 1;
}

int
RLC_LecroySequenceOn (int pad)
{
    char reply[255];

    RLC_GPIBCommunication (pad, ":FUNC:MODE SEQ", reply, false);
    RLC_LecroyErrorCheck (pad);

    return 1;
}

int
RLC_LecroyUserOn (int pad)
{
    char reply[255];

    RLC_GPIBCommunication (pad, ":FUNC:MODE USER", reply, false);
    RLC_LecroyErrorCheck (pad);

    return 1;
}

int
RLC_LecroyOutputOn (int pad)
{
    char reply[255];

    RLC_GPIBCommunication (pad, ":OUTPut:STATe 1", reply, false);
    RLC_LecroyErrorCheck (pad);

    return 1;
}


int
RLC_LecroySyncOn (int pad)
{
    char reply[255];

    /* selects the waveform segment for the SYNC pulse */
    RLC_GPIBCommunication (pad, ":TRAC:SEL 1", reply, false);
    RLC_LecroyErrorCheck (pad);

    /* set width */
    RLC_GPIBCommunication (pad, ":OUTPut:SYNC:WIDTh 20", reply, false);
    RLC_LecroyErrorCheck (pad);

    /* set position */
    RLC_GPIBCommunication (pad, ":OUTPut:SYNC:POSition 0", reply, false);
    RLC_LecroyErrorCheck (pad);

    /* switch SYNC on */
    RLC_GPIBCommunication (pad, ":OUTPut:SYNC ON", reply, false);
    RLC_LecroyErrorCheck (pad);

    return 1;
}

int
RLC_LecroyLoadCode (int               pad,
		    int *             lecroy_code,
		    RSP_ParamStruct * param,
		    int               length_lecroy_code,
		    float             vmax,
		    float             vmin)
{
    unsigned char bytes[500];
    int           length_bytes;
    int           i;
    int           j;
    char          buffer[255];
    char          command[255];
    float         range;
    float         offset;
    int           gpib_board = 0;
    char          reply[255];

    length_bytes = length_lecroy_code * 2;

    for (i = 0; i < param->number_of_codes; i++)
    {
	sprintf (command, ":TRAC:DEF %d,%d", i + 1, length_lecroy_code);
	RLC_GPIBCommunication (pad, command, reply, false);
	RLC_LecroyErrorCheck (pad);

	sprintf (command, ":TRAC:SEL %d", i+1);
	RLC_GPIBCommunication (pad, command, reply, false);
	RLC_LecroyErrorCheck (pad);

	for (j = 0; j < length_bytes; j++)
	{
	    bytes[j] = 0;
	}
	for  (j=0; j < length_lecroy_code; j++)
	{
	    //printf ("%d ", lecroy_code[ (i * length_lecroy_code) + j]);
	    bytes[(j << 1) + 0] = (unsigned char) (lecroy_code[ (i * length_lecroy_code) + j] & 0x00ff);
	    bytes[(j << 1) + 1] = (unsigned char) (lecroy_code[ (i * length_lecroy_code) + j] >> 8);
	}
	//printf ("\n");

	sprintf (buffer, "%d", length_bytes);
	sprintf (command, ":TRAC#%zd%d", strlen (buffer), length_bytes);
	//printf ("%s\n", command);
	Send (gpib_board,  MakeAddr (pad, 0), command, strlen (command), NULLend);
	SendDataBytes (gpib_board, bytes, length_bytes, NULLend);
	usleep (300000);

	RLC_LecroyErrorCheck (pad);
    }

    /* define sequence table */
    if (param->number_of_codes > 1)
    {
	/* SINGle advance mode */
	/* after completing the waveform the lecroy idles (does nothing) while waiting for the next trigger */
	RLC_GPIBCommunication (pad, ":SEQuence:ADV SINGle", reply, false);
	RLC_LecroyErrorCheck (pad);

	/* set source for SEQ ADV to be external */
	RLC_GPIBCommunication (pad, ":SEQuence:ADV:SOURce EXTernal", reply, false);
	RLC_LecroyErrorCheck (pad);

	for (i = 1; i <= param->number_of_codes; i++)
	{
	    sprintf (command, ":SEQ:DEF %d,%d,1,0", i, i);
	    RLC_GPIBCommunication (pad, command, reply, false);
	    RLC_LecroyErrorCheck (pad);
	}
    }

    /* set the voltage output of the waveform */
    range  = vmax - vmin;
    offset = vmin + (range / 2);

    sprintf (command, ":VOLT %04.3f", range);
    RLC_GPIBCommunication (pad, command, reply, false);
    RLC_LecroyErrorCheck (pad);

    sprintf (command, ":VOLT:OFFS %04.3f", offset);
    RLC_GPIBCommunication (pad, command, reply, false);
    RLC_LecroyErrorCheck (pad);

    /* set up the trigger for the Lecroy */
    RLC_GPIBCommunication (pad, ":INIT:CONT OFF", reply, false);
    RLC_LecroyErrorCheck (pad);

    RLC_GPIBCommunication (pad, ":TRIGger:SOURce:ADV EXT", reply, false);
    RLC_LecroyErrorCheck (pad);

    /* using the external SCLK */
    RLC_GPIBCommunication (pad, ":FREQ:RAST:SOUR EXT", reply, false);
    RLC_LecroyErrorCheck (pad);

    return 1;
}
