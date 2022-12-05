#ifndef _COUNTER_H
#define _COUNTER_H

/**********************************************************************

COUNTER.H

Luke Robertson @ RCRU, Rutherford Appletone Lab, Oxon, UK
lsr62@rcru.rl.ac.uk

This is an INCLUDE FILE FOR LINUX PCI_AMCC DEVICE DRIVER

It includes the device IOCTL definitions (using IOC_MAGIC method to avoid
IOCTL numbering conflicts ), also includes user space ioctl wrapper functions.
In the unlikely event where you need to use one of the IOCTL commands,
you might want to borrow this code to speed up the application.

 *********************************************************************/

/* MODIFICATION RECORD */

/* 20120502 - Updated with changes discovered in:
 * /home/chilbolton_software/universal_radar_code/RDQ/include/counter.h
 * These changes were not released to the original sources used for Galileo
 * testing in 2011.  Discovered when upgrading CAMRa.
 */

#ifdef __KERNEL__

#include <linux/types.h>

#ifndef BYTE
typedef __u8  byte;  /* 8-bit  */
typedef __u16 word;  /* 16-bit */
typedef __u32 dword; /* 32-bit */
#define BYTE
#endif

#else  /* __KERNEL__ */

#include <sys/ioctl.h>
#include <stdint.h>

#ifndef BYTE
typedef uint8_t  byte;  /* 8-bit  */
typedef uint16_t word;  /* 16-bit */
typedef uint32_t dword; /* 32-bit */
#define BYTE
#endif

#endif /* __KERNEL__ */

#define COUNTER_IOC_MAGIC 'd'
#define COUNTER_IOC_RESET                 _IO (COUNTER_IOC_MAGIC, 1)
#define COUNTER_IOC_SET_COUNTER           _IO (COUNTER_IOC_MAGIC, 2)
#define COUNTER_IOC_COUNTER_LATCH         _IO (COUNTER_IOC_MAGIC, 3)
#define COUNTER_IOC_READBACK_COUNTER_MODE _IO (COUNTER_IOC_MAGIC, 4)

#define COUNTER_IOC_MAX 4

/* ISACTRL 82C54 Write Only Addressing */
/* #define COUNT0                  0x300  // COUNTER 0        */
/* #define COUNT1                  0x302  // COUNTER 1        */
/* #define COUNT2                  0x304  // COUNTER 2        */
/* #define CONTROL                 0x306  // CONTROL REGISTER */

/* ISACTRL Read Only Addressing */
/* #define USEC1                   0x308  // USE COUNTER 1 TO DIVIDE CLK      */
/* #define FIFORST                 0x310  // RESET PCICARD ADD-ON FIFOS       */
/* #define BYPASSC1                0x318  // BYPASS COUNTER 1 & USE UNDIV CLK */


/* Define port offsets */
#define COUNTER_COUNT0  0x00
#define COUNTER_COUNT1  0x02
#define COUNTER_COUNT2  0x04
#define COUNTER_CONTROL 0x06

/* For writing the control word...... */
/* counter<<6 for writing mode        */
/* counter<<3 for readback            */

/* ***********FOR SETTING MODES**************** */

/* Read/Write counter modes.... */
#define COUNTER_CONTROL_RW_LSB_ONLY 0x10 /* Read/Write LSB ONLY         */
#define COUNTER_CONTROL_RW_MSB_ONLY 0x20 /* Read/Write MSB ONLY         */
#define COUNTER_CONTROL_RW_LSB_MSB  0x30 /* Read/Write LSB AND THEN MSB */

/* need to clear for MODE_ZERO */
#define COUNTER_CONTROL_MODE_ONE   0x02 /* Select mode 1 */
#define COUNTER_CONTROL_MODE_TWO   0x04 /* Select mode 2 */
#define COUNTER_CONTROL_MODE_THREE 0x06 /* Select mode 3 */
#define COUNTER_CONTROL_MODE_FOUR  0x08 /* Select mode 4 */
#define COUNTER_CONTROL_MODE_FIVE  0x0A /* Select mode 5 */
#define COUNTER_CONTROL_BCD        0xFE /* AND value ( set last bit to 0 ) */

/* ************For READBACK *************************** */
#define COUNTER_CONTROL_READBACK 0xF0 /* Select readback ( and set counter/status latches to 1 -> need to to reset them yourself! */

/* ANDs for readback counters........(set desired bit to 0) */
#define COUNTER_CONTROL_READBACK_LATCHCOUNT 0xDF /* bit 5 should be 0                                    */
#define COUNTER_CONTROL_READBACK_STATUS     0xEF /* bit 4 (0-7) should be 0                              */
#define COUNTER_COUNTER_LATCH_COMMAND_AND   0xC0 /* bits 4-5 should be 0 (0-3 probably should be 0 too ) */

/******************************************************
 * User functions ( wrappers for ioctls )
 ******************************************************/

#ifndef __KERNEL__

static inline int
set_counter (int counter_dev,
	     int counter_offset,
	     int counter_value,
	     int counter_mode)
{
    int test_counter_params[6];
    test_counter_params[0] = counter_offset;
    test_counter_params[1] = counter_value;
    test_counter_params[2] = counter_mode;

    return ioctl (counter_dev, COUNTER_IOC_SET_COUNTER, test_counter_params);
}

static inline int
reset_counter (int counter_dev)
{
    /* reset the counters ??? **CHECK THIS ** */
    return ioctl (counter_dev, COUNTER_IOC_RESET);
}

static inline int
read_counter_latch (int counter_dev,
		    int counter_offset)
{
    int test_counter_params[6], retval;

    test_counter_params[0] = counter_offset;

    retval = ioctl (counter_dev, COUNTER_IOC_COUNTER_LATCH, test_counter_params);
    if (retval < 0)
	return -1;
    return test_counter_params[0];
}

#endif /* !__KERNEL__ */
#endif /* !_COUNTER_H */
