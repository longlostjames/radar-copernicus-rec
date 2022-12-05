/*

 */

#ifndef __serialPLC_H
#define __serialPLC_H

/* ========== INCLUDE FILES ========== */

#include "tyerror.h"
#include <time.h>

/* ========== MACRO DEFINITIONS ========== */

#ifndef BYTE
	#define BYTE char
#endif

/* ========== DATA STRUCTURES ========== */

/* data structure to hold a translated serial message */
typedef struct
{
        int     type;
	int	year;
	int	month;
	int	day;
	int	hour;
	int	min;
	int	sec;
	int	millisec;
	double 	az;
	double 	el;
	char    mode;
	char    status;
	char    vel;
} TY_SERIALMSG_MSG;

/* ========== FUNCTION PROTOTYPES ========== */
TY_ERROR *posmsg_readmsg(TY_SERIALMSG_MSG *);
TY_ERROR *posmsg_init(char *);
void posmsg_print(FILE *, TY_SERIALMSG_MSG *);

#endif
