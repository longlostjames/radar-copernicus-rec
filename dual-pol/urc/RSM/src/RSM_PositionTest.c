#include <stdio.h>
#include <unistd.h>

#include <RSM.h>

/*---------------------------------------------------------------------------*
 * Below defines the com port that the serial message arrives on             *
 *---------------------------------------------------------------------------*/
#define SERIALMESSAGE_PORT			"/dev/ttyS1" // Old case

int
main (int   argc,
      char *argv[])
{
    int res;
    RSM_PositionMessageStruct positionmsg;

    res = RSM_InitialisePositionMessage(SERIALMESSAGE_PORT);
    if (res)
    {
	printf("Detected a problem with initialising the position message port\n");
	return 1;
    }

    for (res = 0; res < 100; res++)
    {
	RSM_ReadPositionMessage(&positionmsg);
	printf ("Az: % 8.3f El: % 8.3f %04u/%02u/%02u %02u:%02u:%02u.%02u\n",
	    positionmsg.az,
	    positionmsg.el,
	    positionmsg.year,
	    positionmsg.month,
	    positionmsg.day,
	    positionmsg.hour,
	    positionmsg.min,
	    positionmsg.sec,
	    positionmsg.centi_sec);
	usleep (100000U);
    }

    RSM_ClosePositionMessage();

    return 0;
}
