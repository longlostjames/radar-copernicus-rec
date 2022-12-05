#ifndef __RSM_H__
#define __RSM_H__

#define USE_NETWORK_POSITION_MESSAGE

/* data structure to hold a translated serial message */
typedef struct
{
    double  az;
    double  el;
    int     year;
    int     month;
    int     day;
    int     hour;
    int     min;
    int     sec;
    int     centi_sec;
} RSM_PositionMessageStruct;

extern int RSM_InitialiseSerialMessage (const char * serialport);
extern int RSM_ReadSerialMessage (RSM_PositionMessageStruct * positionmessage);
extern void RSM_CloseSerialMessage (void);

extern int RSM_InitialiseNetworkMessage(void);
extern int RSM_ReadNetworkMessage (RSM_PositionMessageStruct * positionmessage);
extern void RSM_CloseNetworkMessage (void);

#ifdef USE_NETWORK_POSITION_MESSAGE

#define RSM_InitialisePositionMessage(port) RSM_InitialiseNetworkMessage ()
#define RSM_ReadPositionMessage(msgp)       RSM_ReadNetworkMessage (msgp)
#define RSM_ClosePositionMessage()          RSM_CloseNetworkMessage ()

#else /* USE_NETWORK_POSITION_MESSAGE */

#define RSM_InitialisePositionMessage(port) RSM_InitialiseSerialMessage (port)
#define RSM_ReadPositionMessage(msgp)       RSM_ReadSerialMessage (msgp)
#define RSM_ClosePositionMessage()          RSM_CloseSerialMessage ()

#endif /* USE_NETWORK_POSITION_MESSAGE */

#endif /* !__RSM_H__ */
