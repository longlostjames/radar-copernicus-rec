/*
 * RSM_SerialMessage.c
 * This code is a near copy of some earlier code to read the 20ms
 * serial message
 * created by : Owain Davies
 * created on : 15/04/2004
 */

#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <semaphore.h>
#include <termios.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>

#include <RSM.h>

/* definitions */

/* the length of the serial message without the footer */
#define SERIALMESSAGE_LENGTH 13

/* the length of the serial message footer */
#define FOOTER_LENGTH         3

/* the number of bytes to seach through when looking for the serial message footer */
#define FOOTER_SEARCH_LENGTH 30

/* global variables */

/* file descriptor for the serial port */
static int serialport_fd = -1;

/* 'sm_buffer_semaphore' is a semaphore created during initialisation that is
 * used to control shared access to the receive buffer and its associated
 * valid flag. Access to this buffer and flag is shared by the main thread
 * and the read thread.
 */
static sem_t sm_buffer_semaphore;

/* 'sm_read_thread' is a thread created during initialisation that continuously
 * reads the serial message.
 */
static pthread_t     sm_read_thread;

/* 'sm_buffer' is the receive buffer that is used to store the serial message
 */
static unsigned char sm_buffer[ SERIALMESSAGE_LENGTH ];

/* 'sm_buffer_valid' is a flag that indicates if the receive buffer contains a
 * valid message.
 * NOTE: currently this flag is initialised to TRUE to handle the case where
 * application code is trying to read the latest message before a message
 * has been received by the read thread. It might be wiser to handle this
 * situation by returning a non-fatal error indicating that the first message
 * is currently being processed.
 */
static int           sm_buffer_valid = 1;

/* 'read_raw' read a serial message from the serial port
 * In: none
 * Out: 'valid' is updated to indicate if a message was successfully received
 *  TRUE for success, FALSE for failure
 * 'buffer' is updated with the message
 * Returns: void
 */
static void
read_raw (unsigned char * buffer,
	  int *           valid)
{
    int           rx_count     = 0;
    int           footer_count = 0;
    int           retval;
    int           footer_found = 0;
    unsigned char data;

    /* flush tx and rx buffers */
    tcflush (serialport_fd, TCIOFLUSH);

    /* detect the serial message footer */
    while (rx_count < FOOTER_SEARCH_LENGTH && !footer_found)
    {
	retval = read (serialport_fd, &data, 1);
	if (retval == 0)
	{
	    continue;
	}
	else if (retval == -1)
	{
	    ; /* for some reason we cannot report this error */
	}
	else
	{
	    if (data == 0xFF)
	    {
		footer_count ++;
	    }
	    else
	    {
		footer_count = 0;
	    }
	    if (footer_count == FOOTER_LENGTH)
	    {
		footer_found = 1;
	    }
	}
    }

    if (!footer_found)
    {
	*valid = 0;
    }

    /* read the message */
    while (rx_count < SERIALMESSAGE_LENGTH)
    {
	retval = read (serialport_fd, &data, 1);
	if (retval == 0)
	{
	    continue;
	}
	else if (retval == -1)
	{
	    ; /* for some reason we cannot report this error */
	}
	else
	{
	    buffer[ rx_count ++ ] = data;
	}
    }

    *valid = 1;
}

static void
convert_buffer (unsigned char *             buffer,
		RSM_PositionMessageStruct * msg)
{
    unsigned int temp; // hopefully give us a 32 bit number
    time_t       now;
    struct tm    gmt;

    /* AZ_HI = 0, AZ_MID = 1, AZ_LO = 2 */
    temp    = ((buffer[ 0 ] & 0x03) << 16) + (buffer[ 1 ] << 8) + buffer[ 2 ];
    msg->az = ((double) temp / 480.0) - 90.0;

    /* EL_HI = 3, EL_MID = 4, EL_LO = 5 */
    temp    = ((buffer[ 3 ] & 0x03) << 16) + (buffer[ 4 ] << 8) + buffer[ 5 ];
    msg->el = (double) temp / 480.0;

    /* Special case: if bit 2 is set in the EL_HI byte, then elevation is
     * negative.
     */
    if (buffer[ 3 ] & 0x04)
	msg->el = -msg->el;
    /* YEAR = 6, MONTH = 7, DAY = 8, HOUR = 9, MIN = 10, SEC = 11, CENTI_SEC = 12 */
    msg->day        = buffer[ 8 ];
    msg->month      = buffer[ 7 ];
    //msg->year       = buffer[ 6 ];
    msg->hour       = buffer[ 9 ];
    msg->min        = buffer[ 10 ];
    msg->sec        = buffer[ 11 ];
    msg->centi_sec  = buffer[ 12 ];

    // obtaining year from the local computer clock since this is nearly absolute
    time (&now);
    gmtime_r (&now, &gmt);
    msg->year = gmt.tm_year + 1900;
}

static void *
read_thread (void * arg)
{
    unsigned char rxbuffer [ SERIALMESSAGE_LENGTH ];
    int           valid;

    while (1)
    {
	read_raw (rxbuffer, &valid);
	/* copy over the serial message, but first set the semaphore */
	sem_wait (&sm_buffer_semaphore);
	if (valid)
	{
	    memcpy (sm_buffer, rxbuffer, SERIALMESSAGE_LENGTH);
	}
	sm_buffer_valid = valid;
	sem_post (&sm_buffer_semaphore);
    }
    return NULL;
}



int
RSM_InitialiseSerialMessage (const char * serialport)
{
    struct termios attributes;

    printf ("Opening serial port: %s\n", serialport);
    serialport_fd = open (serialport, O_RDONLY);

    /* open the serial port */
    if (serialport_fd == -1)
    {
	printf ("There has been a problem opening the serial port\n");
	return -1;
    }

    /* obtain, modify and set serial port settings */
    if (tcgetattr (serialport_fd, &attributes) == -1)
    {
	printf ("There has been a problem obtaining the attributes for the serial port\n");
	close (serialport_fd);
	return -1;
    }
    cfmakeraw (&attributes);
    attributes.c_iflag |= PARODD;
    attributes.c_cflag |= PARENB;
    cfsetospeed (&attributes, B115200);
    cfsetispeed (&attributes, B115200);
    if (tcsetattr (serialport_fd, TCSANOW, &attributes) == -1)
    {
	printf ("There has been a problem setting the attributes for the serial port\n");
	close (serialport_fd);
	return -1;
    }

    /* create a semaphore to control shared access to the serial message buffer */
    if (sem_init (&sm_buffer_semaphore, 0, 1) == -1)
    {
	printf ("There has been a problem creating the serial message buffer semaphore\n");
	close (serialport_fd);
	return -1;
    }

    /* create a new thread that continuously reads serial messages */
    if (pthread_create (&sm_read_thread, NULL, read_thread, NULL))
    {
	printf ("There has been a problem with creating the read_thread thread\n");
	close (serialport_fd);
	return -1;
    }
    return 0;
}

int
RSM_ReadSerialMessage (RSM_PositionMessageStruct * serialmessage)
{
    unsigned char buffer [ SERIALMESSAGE_LENGTH ];
    int           valid;

    /* read the serialmessage buffer, but first set the semaphore */
    sem_wait (&sm_buffer_semaphore);
    valid = sm_buffer_valid;
    if (valid)
    {
	memcpy (buffer, sm_buffer, SERIALMESSAGE_LENGTH);
    }
    sem_post (&sm_buffer_semaphore);

    if (!valid)
    {
	return 0;
    }

    /* interpret the data in the buffer */
    convert_buffer (buffer, serialmessage);

    return 1;
}

void
RSM_CloseSerialMessage (void)
{
    pthread_cancel (sm_read_thread);
    sem_destroy (&sm_buffer_semaphore);
    /* close serial port */
    close (serialport_fd);
}
