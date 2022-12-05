#ifndef _REL_H
#define _REL_H

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
} REL_SerialMessageStruct;

extern int  REL_InitialiseSerialMessage (const char * serialport);
extern int  REL_ReadSerialMessage       (REL_SerialMessageStruct * serialmessage);
extern void REL_CloseSerialMessage      (void);


#endif /* !_REL_H */
