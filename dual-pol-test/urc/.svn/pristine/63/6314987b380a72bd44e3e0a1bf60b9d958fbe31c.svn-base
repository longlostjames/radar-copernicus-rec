/* a simple test of the clinometer */

#include <stdio.h>

#include <REL.h>

#define SERIALMESSAGE_PORT                 "/dev/ttyS1"


int main()
{
	int	temp_int;

	static REL_SerialMessageStruct serialmsg;

	temp_int = REL_InitialiseSerialMessage(SERIALMESSAGE_PORT);
	if ( temp_int != 0) {
        	printf("Detected a problem with initialising the serial port\n");
        	exit(1);
      	}
      	sleep(2);

	while (1) {
		temp_int = REL_ReadSerialMessage(&serialmsg);
		printf("elevation angle : %f (%d)\n", serialmsg.el, temp_int);
		usleep(100000);
	}
	REL_CloseSerialMessage();

	return (0);



}











