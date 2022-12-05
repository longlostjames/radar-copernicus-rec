/*
 * REL_SerialMessage.c
 * This code is a near copy of some earlier code to read the 20ms
 * serial message
 * created by : Owain Davies
 * created on : 15/04/2004
 * MOD now does PRO3600 level
 */

#include <stdio.h>
#include <stdbool.h>
#include <pthread.h>
#include <semaphore.h>
#include <termios.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

#include <REL.h>

/* definitions */

/* the length of the serial message without the footer */
#define SERIALMESSAGE_LENGTH 9

/* the length of the serial message footer */
#define FOOTER_LENGTH        1

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
static sem_t         sm_buffer_semaphore;

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
static bool          sm_buffer_valid = true;
static time_t        check_time      = -1;


/* 'read_raw' read a serial message from the serial port
 * In: none
 * Out: 'valid' is updated to indicate if a message was successfully received
 *  TRUE for success, FALSE for failure
 * 'buffer' is updated with the message
 * Returns: void
 */
static void
read_raw (unsigned char * buffer,
	  bool *          valid)
{
    int           rx_count     = 0;
    int           footer_count = 0;
    int           retval;
    bool          footer_found = false;
    unsigned char data = 'a';

    /* flush tx and rx buffers */
    tcflush (serialport_fd, TCIOFLUSH);

    /* detect the serial message footer */
    while (rx_count < FOOTER_SEARCH_LENGTH && !footer_found)
    {
	retval = read(serialport_fd, &data, 1);
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
	    if (data == 0x0A)
	    {
		footer_count ++;
	    }
	    else
	    {
		footer_count = 0;
	    }
	    if (footer_count == FOOTER_LENGTH)
	    {
		footer_found = true;
	    }
	}
    }

    if (!footer_found)
    {
	*valid = false;
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
    *valid = (footer_found) ? true : false;
}

void
convert_buffer (unsigned char *           buffer,
		REL_SerialMessageStruct * msg)
{
    struct timeval tv;
    struct tm      tm;
    char           temp_buffer[8];
    float          temp_double = -999.0;

    if  (!((buffer[0] == '+') || (buffer[0] == '-')))
    {
	temp_double = -999.0;
    }
    else
    {
	memset (temp_buffer, 0, sizeof (temp_buffer));
	strncpy (temp_buffer, (char *)&buffer[1], 6);
	sscanf (temp_buffer, "%f", &temp_double);
	if (buffer[0] == '-')
	{
	    temp_double = temp_double * -1.0;
	}
    }

    msg->el = temp_double;
    msg->az = -999;

    gettimeofday (&tv, NULL);
    gmtime_r (&tv.tv_sec, &tm);
    msg->day        = tm.tm_year + 1900;
    msg->month      = tm.tm_mon  + 1;
    msg->year       = tm.tm_mday;
    msg->hour       = tm.tm_hour;
    msg->min        = tm.tm_min;
    msg->sec        = tm.tm_sec;
    msg->centi_sec  = tv.tv_usec / 10000U;
}

static void *
read_thread (void * arg)
{
    while (1)
    {
	unsigned char rxbuffer [ SERIALMESSAGE_LENGTH ];
	bool          valid;

	valid = false;
	read_raw (rxbuffer, &valid);

	/* copy over the serial message, but first set the semaphore */
	sem_wait (&sm_buffer_semaphore);
	if (valid)
	{
	    memcpy (sm_buffer, rxbuffer, SERIALMESSAGE_LENGTH);
	}
	sm_buffer_valid = valid;
	check_time      = time (NULL);
	sem_post (&sm_buffer_semaphore);
    }
    return NULL;
}

int
REL_InitialiseSerialMessage (const char * serialport)
{
    struct termios attributes;

    printf ("Opening serial port: %s\n", serialport);
    serialport_fd = open (serialport, O_RDONLY);

    /* open the serial port */
    if  (serialport_fd == -1)
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
    cfmakeraw(&attributes);

    // attributes.c_iflag |= PARODD;
    // attributes.c_cflag |= PARENB;
    cfsetospeed(&attributes, B9600);
    cfsetispeed(&attributes, B9600);
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
REL_ReadSerialMessage (REL_SerialMessageStruct * serialmessage)
{
    unsigned char buffer [ SERIALMESSAGE_LENGTH ];
    bool          valid;

    /* read the serialmessage buffer, but first set the semaphore */
    sem_wait(&sm_buffer_semaphore);
    valid = sm_buffer_valid;
    if (!(((check_time - time (NULL)) < 2) && ((check_time - time (NULL)) > -2)))
    {
	valid = false;
    }
    if (valid)
    {
	memcpy (buffer, sm_buffer, SERIALMESSAGE_LENGTH);
    }
    sem_post(&sm_buffer_semaphore);

    if (!valid)
    {
	return 1;
    }

    printf ("valid : %d\n", valid);

    /* interpret the data in the buffer */
    convert_buffer (buffer, serialmessage);

    return 0;
}

void
REL_CloseSerialMessage (void)
{
    pthread_cancel (sm_read_thread);
    sem_destroy (&sm_buffer_semaphore);
    /* close serial port */
    close (serialport_fd);
}
