#ifndef __AMCC_H
#define __AMCC_H

/* File		amccdrv.h
 * Project	PCI data acquisition system
 * Pupose	Header file for the AMCC Linux device driver
 * History
 * =======
 * V01.00 - 05-May-2000 - G.Davies - First version.
 *
 * V02.00 - 26-Oct-2000 - G.Davies
 * 1. Changes required for PCI9500. This driver will not work with any previous
 *    PCI cards.
 * 2. Added AMCC_IO_READ_NVRAM.
 * 3. Added AMCC_IO_WRITE_NVRAM.
 * 20101026 CJW Altered AMCC_DEFAULT_WAIT_FOR_DMA_TIMEOUT from 5 to 10
 */

/* ** add legacy no-op define from /usr/include/bits/sigaction.h, not included in kernel src */
# define SA_INTERRUPT 0x20000000 /* Historical no-op.  */

/* Module, device and IRQ line names */
#define AMCC_NAME	"amcc"
#define AMCC_DEV_NAME	"/dev/amcc"
#define AMCC_IRQ        "amccirq"

/* Major device number */
#define AMCC_MAJOR	120

/* Supported IOCTL methods */
#define AMCC_IOCTL_BASE			0x00
#define AMCC_IO_REINIT			_IO  ('z', (AMCC_IOCTL_BASE +  0)        )
#define AMCC_IO_FREE_DMA_BUFFER		_IO  ('z', (AMCC_IOCTL_BASE +  1)        )
#define AMCC_IO_WAIT_FOR_DMA		_IOR ('z', (AMCC_IOCTL_BASE +  2), int * )
#define AMCC_IO_READ_NVRAM		_IOR ('z', (AMCC_IOCTL_BASE + 10), int * )
#define AMCC_IO_SET_DMA_BANK            _IOW ('z', (AMCC_IOCTL_BASE +  3), int   )
#define AMCC_IO_ALLOC_DMA_BUFFER	_IOW ('z', (AMCC_IOCTL_BASE +  4), size_t)
#define AMCC_IO_START_DMA		_IOW ('z', (AMCC_IOCTL_BASE +  5), int   )
#define AMCC_IO_SET_DMA_TIMEOUT		_IOW ('z', (AMCC_IOCTL_BASE +  6), int   )
#define AMCC_IO_SET_TXFIRE_TIMEOUT	_IOW ('z', (AMCC_IOCTL_BASE +  7), int   )
#define AMCC_IO_SET_TEST_DATA		_IOW ('z', (AMCC_IOCTL_BASE +  8), int   )
#define AMCC_IO_FILL_TEST_DATA		_IOW ('z', (AMCC_IOCTL_BASE +  9), int   )
#define AMCC_IO_WRITE_NVRAM		_IOW ('z', (AMCC_IOCTL_BASE + 11), int * )

/* Supported MMAP methods */
#define AMCC_MMAP_DMA_BUFFER		0

/* Define size of NVRAM in bytes */
#define AMCC_NVRAM_SIZE	128

/* Defines for DMA status */
#define DMA_STATUS_FINISHED		0
#define DMA_STATUS_TIMED_OUT		1
#define DMA_STATUS_MASTER_ABORT		2
#define DMA_STATUS_FIFO_OVERFLOW	3
#define DMA_STATUS_NOT_TRIGGERED	4

#ifdef __KERNEL__

/* Boolean definitions */
#ifndef TRUE
#define TRUE  (1 == 1)
#define FALSE (!TRUE)
typedef int   BOOL;
#endif

/* Defines for incoming and outgoing mailbox interrupts */
#define MAILBOX_IN_INT	(WRITE_TC_INT >> 1)
#define MAILBOX_OUT_INT (WRITE_TC_INT >> 2)

/* Define default timeout when waiting for DMA transfer. This is in seconds. */
/*
 * #define AMCC_DEFAULT_WAIT_FOR_DMA_TIMEOUT	5
 * timeouts have been changed to allow radar-acrobat to exit more quickly when
 * it is in the blanking zones was 2 now set to 10 by Owain Davies 20070914
 * this is to allow more scope for collecting data with acrobat vertically pointing
 */
#define AMCC_DEFAULT_WAIT_FOR_DMA_TIMEOUT	10

/* Define default timeout when waiting for a Tx fire. This is in centi-seconds. */
/*
 * #define AMCC_DEFAULT_WAIT_FOR_TXFIRE_TIMEOUT	50
 * timeouts have been changed to allow radar-acrobat to exit more quickly when
 * it is in the blanking zones
 */
#define AMCC_DEFAULT_WAIT_FOR_TXFIRE_TIMEOUT 10

/* Define the timeout when waiting for the PCI card to trigger or report a
 * trigger timeout. If this timeout occurs, then there is a problem with the
 * PCI card or device driver, but we implement the timeout to stop the
 * machine looping forever in kernel mode.
 */
/*
 * #define START_TIMEOUT   100
 * timeouts have been changed to allow radar-acrobat to exit more quickly when
 * it is in the blanking zones. Note this is in the CPLD.
 */
#define START_TIMEOUT 10

/* Macros to set and clear bits in the AMCC registers */
#define AMCC_SET_BITS(val, reg)	 (amcc_write_opreg(val  | amcc_read_opreg(reg), reg))
#define AMCC_CLR_BITS(val, reg)  (amcc_write_opreg(~val & amcc_read_opreg(reg), reg))
#define AMCC_ALL_BITS		 0xFFFFFFFF
#define AMCC_WRITE_NVRAMREG(val) (amcc_write_opreg(((val) << 16) | amcc_read_opreg(AMCC_OP_REG_MCSR), AMCC_OP_REG_MCSR))
#define AMCC_READ_NVRAMREG	 (amcc_read_opreg(AMCC_OP_REG_MCSR) >> 16)

/* Mask associated with enabling/disabling mailbox 4 MSB interrupt in the INTCSR register */
#define AMCC_MAILBOX_BITS	0x1F00

/* Definitions for the NVRAM */
#define AMCC_NVRAM_LOW_ADDRESS	0x8000
#define AMCC_NVRAM_HIGH_ADDRESS	0xA000
#define AMCC_NVRAM_READ		0xE000
#define AMCC_NVRAM_WRITE	0xC000
#define AMCC_NVRAM_BUSY		0x8000

/* Main device driver data structure */
typedef struct
{
	int			open;			 /* has the device been opened ?                         */
	struct pci_dev	        *pdev;			 /* PCI device data structure                            */
	unsigned int	        op_reg_base_address;	 /* base address of operation registers                  */
	size_t		        dma_buffer_size;	 /* size of DMA buffer                                   */
	caddr_t		        physical_dma_bank[ 2 ];	 /* physical addresses of DMA banks                      */
	caddr_t		        virtual_dma_bank [ 2 ];	 /* virtual addresses of DMA banks                       */
	int                     active_dma_bank;         /* bank to perform DMA transfer into                    */
	int			wait_for_dma_timeout;	 /* timeout used when waiting for DMA transfer to finish */
	int			wait_for_txfire_timeout; /* timeout used when waiting for Tx fire                */
	unsigned int		test_data;		 /* test data to fill DMA bank with                      */
	BOOL volatile		dma_finished;
	BOOL volatile		fifo_overflow;
	BOOL volatile		master_abort;
	BOOL volatile		armed;			 /* set when the PCI card has been armed            */
	BOOL volatile		hvtimeout;		 /* set if a hardware timeout occurs waiting for HV */
	BOOL volatile		triggered;
} TY_AMCC_DRV;

/* Device driver waitqueue and spinlock definitions */
static DECLARE_WAIT_QUEUE_HEAD (amcc_drv_wqhead);
static rwlock_t amcc_drv_wqlock = RW_LOCK_UNLOCKED;

/* Function prototypes */
extern int         	init_module(void);
extern void        	cleanup_module(void);

static int  		amcc_ioctl(struct inode *, struct file *, unsigned int, unsigned long);
static int  		amcc_mmap(struct file *, struct vm_area_struct *);
static int  		amcc_open(struct inode *, struct file *);
static int  		amcc_release(struct inode *, struct file *);

static void  		amcc_init_dev_reg(void);
static void		amcc_unregister_device(void);
static void 		amcc_free_dev_reg(void);
static irqreturn_t      amcc_irq_handler (int, void *, struct pt_regs *);

static __u8		amcc_read_nvram(__u16);
static void		amcc_write_nvram(__u16, __u8);

#endif /* __KERNEL__ */
#endif /* !__AMCC_H  */
