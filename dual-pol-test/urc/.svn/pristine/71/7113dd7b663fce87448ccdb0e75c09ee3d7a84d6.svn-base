/*
 * RSM_NetworkMessage.c
 * This code is a near copy of some earlier code to read the 20ms
 * network message
 * created by : Owain Davies
 * created on : 15/04/2004
 */

#include <stdio.h>
#include <pthread.h>
#include <semaphore.h>
#include <termios.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <string.h>

#include <RSM.h>

/* definitions */
#define HELLO_PORT 12345
#define HELLO_GROUP "225.0.0.37"
#define UDP_PADDING 12

/* data structure to hold a translated serial message */
typedef struct __attribute__ ((__packed__))
{
    char    header[UDP_PADDING];
    int32_t type;
    int32_t year;
    int32_t month;
    int32_t day;
    int32_t hour;
    int32_t min;
    int32_t sec;
    int32_t millisec;
    double  az;
    double  el;
    char    status;
    char    mode;
    char    vel;
    char    padding;
    double  vel_az;
    double  vel_el;
    double  tg_pos_az;
    double  tg_vel_az;
    double  tg_acc_az;
    int32_t tg_dem_mode_az;
    int32_t tg_mode_az;
    double  err_pos_az;
    double  tg_pos_el;
    double  tg_vel_el;
    double  tg_acc_el;
    int32_t tg_dem_mode_el;
    int32_t tg_mode_el;
    double  err_pos_el;
    double  trq_m1;
    double  trq_m2;
    double  trq_m3;
    double  trq_m4;
    double  trq_m5;
    double  trq_m6;
    int32_t status_d1;
    int32_t status_d2;
    int32_t status_d3;
    int32_t status_d4;
    int32_t status_d5;
    int32_t status_d6;
    int32_t status_az;
    int32_t status_el;
    double  lim_soft_az_ccw;
    double  lim_soft_az_cw;
    double  lim_soft_el_low;
    double  lim_soft_el_high;
    double  vel_m1;
    double  vel_m2;
    double  vel_m3;
    double  vel_m4;
    double  vel_m5;
    double  vel_m6;
} TY_SERIALMSG_MSG;

/* global variables */

/* file descriptor for the network port */
static int	networkport_fd = -1;

/* 'sm_buffer_semaphore' is a semaphore created during initialisation that is
 * used to control shared access to the receive buffer and its associated
 * valid flag. Access to this buffer and flag is shared by the main thread
 * and the read thread.
 */
static sem_t	sm_buffer_semaphore;

/* 'sm_read_thread' is a thread created during initialisation that continuously
 * reads the network message.
 */
static pthread_t        sm_read_thread;

/* 'sm_buffer' is the receive buffer that is used to store the network message
 */
static TY_SERIALMSG_MSG sm_buffer;

/* 'sm_buffer_valid' is a flag that indicates if the receive buffer contains a
 * valid message.
 */
static int             sm_buffer_valid = 0;


static void *read_thread(void *arg)
{
	while (1) {
		int status;
		fd_set ipset;
		struct timeval tv;

		/* Do Select with zero timeout to detect command input */
		FD_ZERO (&ipset);
		FD_SET (networkport_fd, &ipset); /* Serial Port */
		tv.tv_sec  = 0U;
		tv.tv_usec = 100000U;
		status = select (networkport_fd+1, &ipset, NULL, NULL, &tv);
		if (status < 0) break;
		if (status > 0 && FD_ISSET (networkport_fd, &ipset))
		{
			TY_SERIALMSG_MSG rxbuffer;

			status = recv (networkport_fd, &rxbuffer, sizeof (rxbuffer), 0);
			if (status == sizeof (rxbuffer) &&
			    (rxbuffer.type == 0 || rxbuffer.type == 1) &&
			    !strcmp (rxbuffer.header, "POSITION"))
			{
				struct tm tm;
				struct timeval tv;
				struct timeval tvl = { 0U, 20000U };
				struct timeval tva;

				rxbuffer.az -= 90.0; /* Adjust from Chobs coordinates */

				/* Fixup time because PLC time is not good */
				gettimeofday (&tv, NULL);
				timersub (&tv, &tvl, &tva); /* Adjust for 20ms latency */
				gmtime_r (&tva.tv_sec, &tm);

				rxbuffer.year     = tm.tm_year + 1900;
				rxbuffer.month    = tm.tm_mon + 1;
				rxbuffer.day      = tm.tm_mday;
				rxbuffer.hour     = tm.tm_hour;
				rxbuffer.min      = tm.tm_min;
				rxbuffer.sec      = tm.tm_sec;
				rxbuffer.millisec = tva.tv_usec / 1000U;

				/* copy over the network message, but first set the semaphore */
				sem_wait( &sm_buffer_semaphore );
				memcpy (&sm_buffer, &rxbuffer, sizeof (rxbuffer));
				sm_buffer_valid = 1;
				sem_post( &sm_buffer_semaphore );
			}
		}
	}
	return NULL;
}


int RSM_InitialiseNetworkMessage(void)
{
	struct sockaddr_in addr;
	struct ip_mreqn mreq;
	u_int yes=1;

	memset (&sm_buffer, 0, sizeof (sm_buffer));

	printf("Opening network port: %s:%u\n", HELLO_GROUP, HELLO_PORT);
	networkport_fd = socket(AF_INET,SOCK_DGRAM,0);

	/* open the network port */
	if ( networkport_fd == -1 ) {
		printf("There has been a problem opening the network port\n");
		return (-1);
	}

	/* allow multiple sockets to use the same PORT number */
	if (setsockopt(networkport_fd,SOL_SOCKET,SO_REUSEADDR,&yes,sizeof(yes)) < 0)
	{
	    perror("Reusing ADDR failed");
	    return -1;
	}

	/* set up destination address */
	memset(&addr,0,sizeof(addr));
	addr.sin_family=AF_INET;
	addr.sin_addr.s_addr=htonl(INADDR_ANY); /* N.B.: differs from sender */
	addr.sin_port=htons(HELLO_PORT);

	/* bind to receive address */
	if (bind(networkport_fd,(struct sockaddr *) &addr,sizeof(addr)) < 0)
	{
	    perror("bind");
	    return -1;
	}

	/* use setsockopt() to request that the kernel join a multicast group */
	mreq.imr_multiaddr.s_addr=inet_addr(HELLO_GROUP);
	mreq.imr_address.s_addr=htonl(INADDR_ANY);
	mreq.imr_ifindex = 0;
	if (setsockopt(networkport_fd,IPPROTO_IP,IP_ADD_MEMBERSHIP,&mreq,sizeof(mreq))  < 0)
	{
	    perror("setsockopt");
	    return -1;
	}

	/* create a semaphore to control shared access to the network message buffer */
	if ( sem_init(&sm_buffer_semaphore, 0, 1) == -1){
		printf("There has been a problem creating the network message buffer semaphore\n");
		setsockopt(networkport_fd,IPPROTO_IP,IP_DROP_MEMBERSHIP,&mreq,sizeof(mreq));
		close( networkport_fd);
		return (-1);
	}

	/* create a new thread that continuously reads network messages */
	if ( pthread_create( &sm_read_thread, NULL, read_thread, NULL)) {
		printf("There has been a problem with creating the read_thread thread\n");
		setsockopt(networkport_fd,IPPROTO_IP,IP_DROP_MEMBERSHIP,&mreq,sizeof(mreq));
		close( networkport_fd);
		return (-1);
	}
	return (0);
}

int RSM_ReadNetworkMessage ( RSM_PositionMessageStruct *networkmessage)
{
	int valid;

	/* read the networkmessage buffer, but first set the semaphore */
	sem_wait(&sm_buffer_semaphore);
	valid = sm_buffer_valid;
	if (valid)
	{
	    networkmessage->az        = sm_buffer.az;
	    networkmessage->el        = sm_buffer.el;
	    networkmessage->year      = sm_buffer.year;
	    networkmessage->month     = sm_buffer.month;
	    networkmessage->day       = sm_buffer.day;
	    networkmessage->hour      = sm_buffer.hour;
	    networkmessage->min       = sm_buffer.min;
	    networkmessage->sec       = sm_buffer.sec;
	    networkmessage->centi_sec = sm_buffer.millisec / 10;
	}
	else
	{
	    memset (networkmessage, 0, sizeof (RSM_PositionMessageStruct));
	}
	sem_post(&sm_buffer_semaphore);

	return valid ? 1 : 0;
}

void RSM_CloseNetworkMessage (void)
{
	struct ip_mreqn mreq;

	pthread_cancel( sm_read_thread );
	sem_destroy( &sm_buffer_semaphore );

	/* use setsockopt() to request that the kernel leave a multicast group */
	mreq.imr_multiaddr.s_addr=inet_addr(HELLO_GROUP);
	mreq.imr_address.s_addr=htonl(INADDR_ANY);
	mreq.imr_ifindex = 0;
	setsockopt(networkport_fd,IPPROTO_IP,IP_DROP_MEMBERSHIP,&mreq,sizeof(mreq));

	/* close network port */
	close(networkport_fd);
}



