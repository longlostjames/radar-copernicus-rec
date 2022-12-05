/*
   Module for serial port communication with PLC

   Call OpenWDPort() at the start of the program, and
   CloseWDPort() at the end.

   D.N. Ladd 11/02/2010 - Position message not received is now fatal error.
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <termios.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <signal.h>

#include <RSM_SerialPLC.h>

/* ========== MACRO DEFINITIONS ========== */

/* the length of a message WITHOUT the footer */
#define MSG_LENGTH		31
#define MSG_BUFFER	MSG_LENGTH + 1

#define SEND_BUFFER_SIZE 40
#define RECEIVE_BUFFER_SIZE 256

#define READ_TIMEOUT 1000000

#define ACK 0x06
#define NAK 0x15

/* boolean definitions */
#ifndef BOOL
	#define BOOL  unsigned char
	#define TRUE  (1 == 1)
	#define FALSE (!TRUE)
#endif

static int GetPositionMessage(char buffer[]);

/* ========== ERROR REPORTS ========== */

static TY_ERROR ERR_OPEN          = { 1, "Could not open serial port.", TRUE };
static TY_ERROR ERR_GETATTR       = { 2, "Could not get serial port attributes.", TRUE };
static TY_ERROR ERR_SETATTR       = { 3, "Could not set serial port attributes.", TRUE };
static TY_ERROR ERR_SEM_INIT      = { 4, "Could not create semaphore.", TRUE };
static TY_ERROR ERR_THREAD_CREATE = { 5, "Could not create thread.", TRUE };
static TY_ERROR ERR_INVALID	  = { 6, "No position and time message received.", TRUE };

/* ========== STATIC DATA ========== */

/* 'm_read_thread' is a thread created during initialisation that continuously
 * reads WATCHDOG messages.
 */
static pthread_t 	m_read_thread;

/* 'fd_position' and 'fd_control' are the file descriptors for the serial ports
    used by this module */
static int		fd_position = -1;

/* 'm_previous_attributes' is used to preserve the attrbiutes for the serial
 * port used by this module.
 */
static struct termios	m_previous_attributes_position;

/* 'm_rxbuffer' is the receive buffer that is used to store WATCHDOG
 * messages.
 */
static BYTE		m_rxbuffer[ RECEIVE_BUFFER_SIZE ];

/* 'm_rxbuffer_valid' is a flag that indicates if the receive buffer contains a
 * valid message.
 * NOTE: currently this flag is initialised to TRUE to handle the case where
 * application code is trying to read the latest message before a message
 * has been received by the read thread. It might be wiser to handle this
 * situation by returning a non-fatal error indicating that the first message
 * is currently being processed.
 */
static BOOL		m_rxbuffer_valid;

/* 'm_rxbuffer_semaphore' is a semaphore created during initialisation that is
 * used to control shared access to the receive buffer and its associated
 * valid flag. Access to this buffer and flag is shared by the main thread
 * and the read thread.
 */
static sem_t		m_rxbuffer_semaphore;

/* Delay between successive commands to watchdog.  */
/* Value is in microsceonds. 50mS (50000uS) seems to be sufficient */
#define WD_DELAY 50000

/***********************************************************************
Compute the checksum for the string check_string and return it in
check_char
***********************************************************************/
void checksum(char check_string[], char check_char[])
{
    int i, sum, check;

    i = 0;
    sum = 0;
    while (check_string[i])
    {
      sum += (int)check_string[i];
      i++;
    }

    sum %= 256;
    check = sum >> 4;
    if (check < 10) check_char[0] = check + '0';
    else check_char[0] = check - 10 + 'A';

    check = sum % 16;
    if (check < 10) check_char[1] = check + '0';
    else check_char[1] = check - 10 + 'A';

}


static void close_position_port(void)
{
  if (fd_position != -1)
  {
    tcsetattr(fd_position, TCSADRAIN, &m_previous_attributes_position);
    close(fd_position);
    fd_position = -1;
  }
}

static TY_ERROR *open_position_port(char *PLCPositionPort)
{
  struct termios attributes;

//  printf("Opening Position port:  \n");

  if (fd_position != -1)
    close_position_port();
  fd_position = open(PLCPositionPort, O_RDONLY);

  if (fd_position == -1)
    return &ERR_OPEN;

  if (tcgetattr(fd_position, &attributes) == -1)
  {
    close(fd_position);
    return &ERR_GETATTR;
  }

  m_previous_attributes_position = attributes;
  cfmakeraw(&attributes);

  /*  Baud 115200, no parity, 8 bits */
  attributes.c_cflag = CS8 | CLOCAL | CREAD;
  attributes.c_iflag = IGNPAR;
  attributes.c_oflag = 0;
  attributes.c_lflag = 0;

  cfsetospeed(&attributes, B115200);
  cfsetispeed(&attributes, B115200);
  if (tcsetattr(fd_position, TCSANOW, &attributes) == -1)
  {
    close(fd_position);
    return &ERR_SETATTR;
  }

  return NULL;
}


/* 'read_thread' is the thread that continuously reads  PLC position messages.
 * In:		'arg' is a thread argument; this is not used.
 * Out:		None.
 * Returns:	void.
 */
static void *read_thread(void *arg)
{
  int	 i, numbuf;
  BYTE	 rxbuffer[ RECEIVE_BUFFER_SIZE ];
  BOOL 	 valid;

  while (1)
  {
/*  read the serial message into a local rx buffer */
    numbuf = GetPositionMessage(rxbuffer);

    if (numbuf >= MSG_LENGTH )
    {
      valid = TRUE;
      rxbuffer[MSG_LENGTH] = '\0';
//      printf("Message string: %s\n", rxbuffer);
    }
    else
      valid = FALSE;

/*  Copy the message across into the static rx buffer, using the semaphore to prevent access from the main
    thread. */
    sem_wait(&m_rxbuffer_semaphore);

    if (valid)
    {
      for (i = 0; i < MSG_BUFFER; i++)
        m_rxbuffer[ i ] = rxbuffer[ i ];
    }

    m_rxbuffer_valid = valid;

    sem_post(&m_rxbuffer_semaphore);
  }
  return NULL;
}

/* 'convert' interprets the raw bytes received in a serial message, and
 * translates them into real values.
 * In:		'*buffer' contains the raw message bytes.
 * Out:		'*msg' is written to with the converted message values.
 * Returns:	void.
 */
static void convert(BYTE *buffer, TY_SERIALMSG_MSG *msg)
{
  int iaz, iel;

  sscanf(buffer, "%1d%02d%02d%02d%02d%02d%02d%03d%06d%06d%c%c%c",
         &msg->type, &msg->year, &msg->month, &msg->day, &msg->hour,
         &msg->min, &msg->sec, &msg->millisec, &iaz, &iel, &msg->status,
         &msg->mode, &msg->vel);
  msg->az = (double)(iaz / 1000.0);
  msg->el = (double)(iel / 1000.0);

}

/* ========== EXPORTED FUNCTIONS ========== */

/* 'posmsg_init' initialises the position message module. This must be called
 * prior to using any functions exported by this module.
 * In:		None
 * Out:		None.
 * Returns:	TY_ERROR * : an error structure. This will be NULL if
 *		successful.
 */
TY_ERROR *posmsg_init(char *portname)
{
    // TY_ERROR err_val, *error = &err_val;

  open_position_port(portname);

  /* Create a sempahore to control shared access to the rx message buffer. The semaphore is local to this process, and has an initial value of 1. */
  if (sem_init(&m_rxbuffer_semaphore, 0, 1) == -1)
    return &ERR_SEM_INIT;

  /* create a new thread of execution which continuously reads serial messages */
  if (pthread_create(&m_read_thread, NULL, read_thread, NULL))
    return &ERR_THREAD_CREATE;

//  printf("opened\n");
  return NULL;
}

/* 'posmsg_free' frees all resources used by this module.
 * In:		None.
 * Out:		None.
 * Returns:	void.
 */
void posmsg_free(void)
{
  close_position_port();
  pthread_cancel(m_read_thread);
  sem_destroy(&m_rxbuffer_semaphore);
}

/* 'posmsg_readmsg' returns the most recently received message as
 * translated message data.
 * In:		None.
 * Out:		'*msg' is written to with the translated message data.
 * Returns:	TY_ERROR * : an error structure. This will be NULL is
 *		successful.
 */
TY_ERROR *posmsg_readmsg(TY_SERIALMSG_MSG *msg)
{
  int  i;
  BYTE buffer[ MSG_BUFFER+1 ];
  BOOL valid;

  /* read the static rx buffer, using the semaphore to prevent access from the read message thread */
  sem_wait(&m_rxbuffer_semaphore);
  valid = m_rxbuffer_valid;
//  printf("In readmsg routine\n");
  if (valid)
  {
//    printf("Valid buffer:  ");
    for (i = 0; i < MSG_BUFFER; i++)
    {
      buffer[ i ] = m_rxbuffer[ i ];
//      printf("%c", buffer[i]);
    }
    buffer[i] = '\0';
  }

  sem_post(&m_rxbuffer_semaphore);
  if (!valid)
    return &ERR_INVALID;

  /* interpret the raw data */
  convert(buffer, msg);

  return NULL;
}

static int GetPositionMessage(char buffer[])
{
    char c, check_char[3], rcs[3];
    int i, num_buf;
    //BOOL timeout = FALSE;
    int good;

/*  Read in all the bytes of the incoming message until ETX */
    c = 0;

//    printf("Read serial data: ");
    while (c != '$')
    {
	int status = read(fd_position, &c, 1);
	if (status != 1) return 0;
//      printf("%d %c, ", (int)c, c);
    }
//    printf("\nFound STX: ");
//    while (c != ETX && i < RECEIVE_BUFFER_SIZE)
    i = 0;
    while ((c != '\r') && (i < RECEIVE_BUFFER_SIZE))
    {
      int status = read(fd_position, &c, 1);
      if (status != 1) return 0;
      buffer[i] = c;
      i++;
//      printf("%c", c);
    }

    if (i < RECEIVE_BUFFER_SIZE && (c == '\r'))
    {
      num_buf = i - 1;
      buffer[num_buf] = '\0';
      good = 1;
    }
    else
    {
      num_buf = 0;
      good = 0;
    }

//  Check if there is a checksum.
    if (good && (buffer[num_buf - 3] == '*'))
    {
      strncpy(rcs, &buffer[num_buf - 2], 2);
      rcs[2] = '\0';
//    Reduce the end point of the command string, and adjust the length variable.
      num_buf = num_buf - 3;
      buffer[num_buf] = '\0';

//    Check to see if the checksum characters are in the valid range, and
//    if so, then compute the string checksum
      if (isalnum(rcs[0]) && isalnum(rcs[1]))
      {
        checksum(buffer, check_char);
//        printf("Checksum string %s rcs %s %c%c\n", buffer, rcs, check_char[0], check_char[1]);
        if (check_char[0] == rcs[0] && check_char[1] == rcs[1])
          good = 1;
        else
          good = 0;
      }
      else
        good = 0;
    }

    if (good)
      num_buf = num_buf;
    else
      num_buf = 0;

    return(num_buf);
}


#if 0

void sigint_handler(int sig)
{
  printf("Exitting on interrupt\n");
  posmsg_free();
  exit(0);
}

int main(int argc, char *argv[])
{
  int i, numbuf;
  char c;
  char buffer[256];
  TY_SERIALMSG_MSG 	msg;
  TY_ERROR err_val, *err = &err_val;

  if (argc < 2)
  {
    printf("You must enter a serial port. Usage: dispPLC_115200 /dev/ttySx\n\n");
    exit(0);
  }

  if (signal(SIGINT, sigint_handler) == SIG_ERR)
  {
    printf("Error setting up signal interrupt handler\n");
    exit(1);
  }
  if (signal(SIGTERM, sigint_handler) == SIG_ERR)
  {
    printf("Error setting up signal interrupt handler\n");
    exit(1);
  }


  posmsg_init(argv[1]);
  printf("Port opened\n");
  usleep(1000000);

  while (1)
  {

    err = posmsg_readmsg(&msg);
    if (err != NULL)
    {
      printf("Error reading msg\n");
      exit(-1);
    }

    printf("%04d%02d%02d%02d%02d%02d%03d %7.3f %+08.3f %c%c%c",
		msg.year + 2000, msg.month, msg.day,
		msg.hour, msg.min, msg.sec, msg.millisec,
		msg.az, msg.el, msg.mode, msg.status, msg.vel);
//    printf("%s", &buffer[31]);
    printf("\n");
  }
}

#endif
