/* File		tyerror.h
 * Project	PCI data acquisition and signal processing
 * Purpose	Header file for Error module
 * History
 * =======
 * V01.00 - 14-Jun-2000 - G.Davies - First version.
 * V01.01 - 27-Jun-2000 - G.Davies - 1. Renamed this module 'tyerror' to make its name more unique.
 *				     2. Removed function protoype for 'print_error'.
 */


#ifndef __TYERROR_H
#define __TYERROR_H


/* ========== MACRO DEFINITIONS ========== */

#ifndef TRUE
  #define TRUE  (1 == 1)
  #define FALSE (!TRUE)
  typedef int   BOOL;
#endif


/* ========== DATA STRUCTURES ========== */

/* Data structure to hold an error report */
typedef struct
{
	int  errno;
	char errmsg[ 256 ];
	BOOL fatal;
} TY_ERROR;


#endif

