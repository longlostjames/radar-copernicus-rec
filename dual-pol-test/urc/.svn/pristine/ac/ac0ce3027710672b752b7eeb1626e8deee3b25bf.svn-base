#ifndef _RLC_H
#define _RLC_H

#include <stdbool.h>
#include <RSP.h>

extern int  RLC_GPIBCommunication    (int pad, const char * msg,
				      char * reply, bool ctrl);
extern int  RLC_LecroyErrorCheck     (int pad);
extern int  RLC_LecroyReset          (int pad);
extern int  RLC_LecroySelectChannel  (int pad, int channel);
extern int  RLC_LecroySequenceOn     (int pad);
extern int  RLC_LecroyUserOn         (int pad);
extern int  RLC_LecroySyncOn         (int pad);
extern int  RLC_LecroyOutputOn       (int pad);
extern int  RLC_LecroyLoadCode       (int pad, int * lecroy_code,
				      RSP_ParamStruct * param,
				      int length_lecroy_code,
				      float vmax, float vmin);

extern int  RLC_InitialiseSerialPort (const char * serialport);
extern void RLC_CloseSerialPort      (int serialport);

#endif /* !_RLC_H */
