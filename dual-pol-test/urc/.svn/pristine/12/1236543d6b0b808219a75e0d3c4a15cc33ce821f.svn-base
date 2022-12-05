
/* these functions will talk to the NI ENET gpib controller */
/* Owain Davies 20050422 */

#include <stdio.h>
#include <string.h>
#include <ni488.h>

int
enet_send_gpib_command (int          gpib_board,
			Addr4882_t   instrument,
			const char * buffer)
{
    SendIFC (gpib_board);

    DevClear (gpib_board, instrument);
    Send (gpib_board, instrument, buffer, strlen (buffer), NLend);

    return 0;
}

int
enet_send_receive_gpib_command (int          gpib_board,
				Addr4882_t   instrument,
				const char * buffer,
				char *       reply)
{
    SendIFC (gpib_board);

    DevClear (gpib_board, instrument);
    Send (gpib_board, instrument, buffer, strlen (buffer), NLend);
    //usleep(100000);
    Receive (gpib_board, instrument, reply, 255, STOPend);
    if  (ibsta & ERR)
    {
	printf("An error has occured while talking\n");
	return 1;
    }
    reply[ibcnt] = '\0';
    return 0;
}



