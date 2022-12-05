#define VERSION_NUMBER "0.24"

// radar-copernicus-rec.c
// Development Version
// ---------------------------------------------------------------
// Acquisition and signal processing software for 35 GHz radar
//
// Owain Davies / Ed Pavelin - May 2004
// Based on radar-acrobat-rec
//
// ---------------------------------------------------------------
// REVISION HISTORY
// ---------------------------------------------------------------
// 17/05/04: EGP: Created v0.1
// 07/06/04: OTD: converted the program to netCDF
// 08/06/04: OTD: added the ability to output spectra data to netCDF
// 10/06/04: EGP: added support for calibration file
// 05/07/04: EGP: v0.3: Now estimates noise from top range gates
// 19/07/04: EGP: v0.5: Use npsd=nfft
// 27/07/04: EGP: Implemented moments averaging
// 27/07/04: EGP: Added sigma v bar calculation (VEL_HCD)
// 09/09/04: EGP: Started development version for pulse coding
// 18/11/04: OTD: elevation angle from the clinometer added
// 08/12/04: EGP: added sigma-Zbar calculation
// 03/11/06: OTD: added the abiltiy to record SPECTRA_RAPID
// 18/04/07: OTD: improved the precision on time
// 22/05/07: OTD: added moments_averaged to global attributes in netcdf file
// 11/09/07: OTD: i and q, in addition to the spectra (version 0.12)]
// 28/09/07: OTD: debugging the i and q addition to spectra
// 15/02/08: OTD: adding in the dual polarisation capability
// 17/03/08: OTD: streamlining the calibration coefficients
// 15/01/09: OTD: allowing both single pol and co-pol operation
// 23/11/09: JCN: fix bug to synchronise coded pulse order
// 14/09/10: JCN: include new dual-pol parameters phidp and rhohv
// 06/12/13: JCN: committed co and cross channels
// 03/06/15: JCN: fix bug in phidp and rhohv
// 03/08/15: JCN: include interpolated power rhohv
// 14/12/15: JCN: include noise power recording
//
#include <stdio.h>
#include <stdint.h>
#include <fcntl.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>

#include "median.h"
#include <complex.h>
#include <fftw3.h>

//#define PI 3.141592654
#define CLINOMETER_PORT "/dev/ttyS1"


// Master header file for the Universal Radar Code
#include <radar.h>
// Include file for the RSP package
#include <RSP.h>
// Include file for the RDQ package
#include <RDQ.h>
// Include file for the RNC package
#include <RNC.h>
// Include file for the RSM package
#include <RSM.h>
// Include file for the REL package
#include <REL.h>

/* netCDF : netCDF library header file */
#include <netcdf.h>
/* header file */
#include "radar-copernicus-rec.h"
/* below defines the com port that the serial message arrives on */
/* 1 is /dev/ttyS0 or COM1 in MS-DOS langauge */
#define CLINOMETERMESSAGE_PORT                 "/dev/ttyS1"

#define RSP_MOMENTS 5

//-----------------------------
// GLOBAL VARIABLE DEFINITIONS
//-----------------------------

int     dmux_table[8];

/* used for communication between the signal handler */
/* and the main program */
int exit_now = 0;

/* function prototype declaration */
static void sig_handler (int sig);

// Displays a welcome message with version information
static inline void
disp_welcome_message (void)
{
    printf ("\nradar-copernicus-rec: Version %s\n\n", VERSION_NUMBER);
}

/* signal handler */
static void
sig_handler (int sig)
{
    if (exit_now == 0)
    {
        exit_now = 1;
        printf ("***********************************\n");
        printf ("* Received signal %d. Exiting soon *\n", sig);
        printf ("***********************************\n");
    }
}


static inline int
parseargs (int argc, char *argv[], URC_ScanStruct *scan)
{
    //-------------------------
    // Parse command line args
    //-------------------------

    time_t       system_time;
    struct tm    tm;
    const char * operator;

    // Initialise defaults for scan params
    scan->scanType      = SCAN_FIX;
    scan->file_number   = 0;
    scan->scan_number   = 0;
    scan->experiment_id = 0;
    scan->scan_velocity = -9999;
    scan->dwelltime     = -1;

    operator = getenv ("USER");
    if (operator == NULL) operator = getenv ("USERNAME");
    if (operator == NULL) operator = getenv ("LOGNAME");
    if (operator == NULL) operator = "<unknown>";
    strncpy (scan->operator, operator, sizeof (scan->operator) - 1);
    scan->operator[sizeof (scan->operator) - 1] = '\0';

    system_time = time (NULL);
    gmtime_r (&system_time, &tm);
    strftime (scan->date, sizeof (scan->date), "%Y%m%d%H%M%S", &tm);

    return 0;
}


//--------------------------------------------------------------------
// get_config : reads radar config file
static inline void
get_config (char *filename, RSP_ParamStruct *param, URC_ScanStruct *scan, int is_coded)
{
    FILE *file;
    int i, j;
    char codefile[255];
    char dummy[80];
    int tmp_int;

    printf ("Accessing config file: %s\n", filename);

    param->frequency                  = RNC_GetConfigDouble (filename, "radar-frequency");
    param->prf                        = RNC_GetConfigDouble (filename, "prf");
    param->transmit_power             = RNC_GetConfigDouble (filename, "transmit-power");
    param->pulses_per_daq_cycle       = RNC_GetConfigDouble (filename, "pulses");
    param->samples_per_pulse          = RNC_GetConfigDouble (filename, "samples");
    param->ADC_channels               = RNC_GetConfigDouble (filename, "adc-channels");
    param->clock_divfactor            = RNC_GetConfigDouble (filename, "adc-divfactor");
    param->delay_clocks               = RNC_GetConfigDouble (filename, "adc-delayclocks");
    param->pulse_period               = RNC_GetConfigDouble (filename, "chip-length");
    param->pulses_coherently_averaged = RNC_GetConfigDouble (filename, "num-coh-avg");
    param->spectra_averaged           = RNC_GetConfigDouble (filename, "num-spec-avg");
    param->moments_averaged           = RNC_GetConfigDouble (filename, "num-moments-avg");
    param->fft_bins_interpolated      = RNC_GetConfigDouble (filename, "reject-clutter-bins");
    param->clock                      = RNC_GetConfigDouble (filename, "adc-clock");
    param->num_peaks                  = RNC_GetConfigDouble (filename, "num-peaks");
    param->antenna_diameter  = RNC_GetConfigFloat (filename, "antenna_diameter");
    param->beamwidthH  = RNC_GetConfigFloat (filename, "beamwidthH");
    param->beamwidthV  = RNC_GetConfigFloat (filename, "beamwidthV");
    param->height  = RNC_GetConfigFloat (filename, "height");
    param->azimuth_offset = RNC_GetConfigFloat (filename, "azimuth_offset");
    param->dump_spectra    = RNC_GetConfigFloat (filename, "dump_spectra");
    param->dump_spectra_rapid   = RNC_GetConfigFloat (filename, "dump_spectra_rapid");
    param->num_interleave       = RNC_GetConfigFloat (filename, "num-interleave");
    param->num_tx_pol           = RNC_GetConfigFloat (filename, "num-tx-pol");

    scan->min_angle=-9999;
    scan->max_angle=-9999;
    scan->scan_angle=-9999;

    scan->scan_angle  = RNC_GetConfigFloat (filename, "antenna_azimuth");

    RNC_GetConfig (filename, "code-file", codefile, sizeof (codefile));

    strcpy (param->code_name, "NOT YET IMPLEMENTED\0");

    if (is_coded == 1)
    {

	printf ("Reading pulse codes from file: ");
	printf ("%s\n", codefile);

	file=fopen (codefile, "r");
	if (file==NULL)
	{
	    printf ("** ERROR: Unable to open pulse code file!\n");
	    exit (1);
	}
	if (fscanf (file, "%79s %d", dummy, &param->code_length) != 2)
	    exit (3);
	if (fscanf (file, "%79s %d", dummy, &param->number_of_codes) != 2)
	    exit (3);

	for (i=0; i<param->number_of_codes; i++)
	{
	    for (j=0; j<param->code_length; j++)
	    {
		if (fscanf (file, "%d", &tmp_int) != 1)
		    exit (3);
		param->codes[i][j]= (short)tmp_int;
	    }
	}

	fclose (file);
	if (param->num_interleave>1 && param->number_of_codes<2)
	{
	    printf ("** get_config: pulse code file doesn't contain enough pulses!");
	    exit (1);
	}
	else
	{
	    param->number_of_codes = param->number_of_codes / param->num_interleave / param->num_tx_pol;
	}
    }
    else
    {
	// For non-coded pulses
	param->code_length=1;
	param->number_of_codes=1;
    }
}

// Loads calibration information from *.cal file
static inline void
get_cal (RSP_ParamStruct *param, char *calfile)
{
    /* ZED */
    param->ZED_HC_calibration_offset  = RNC_GetConfigDouble (calfile, "ZED_HC_calibration_offset");
    param->ZED_XHC_calibration_offset = RNC_GetConfigDouble (calfile, "ZED_XHC_calibration_offset");
    param->ZED_VC_calibration_offset = RNC_GetConfigDouble (calfile, "ZED_VC_calibration_offset");
    param->ZED_XVC_calibration_offset = RNC_GetConfigDouble (calfile, "ZED_XVC_calibration_offset");
    param->ZED_HCP_calibration_offset       = RNC_GetConfigDouble (calfile, "ZED_HC_calibration_offset");
    param->ZED_XHCP_calibration_offset      = RNC_GetConfigDouble (calfile, "ZED_XHC_calibration_offset");
    param->ZED_VCP_calibration_offset       = RNC_GetConfigDouble (calfile, "ZED_VC_calibration_offset");
    param->ZED_XVCP_calibration_offset      = RNC_GetConfigDouble (calfile, "ZED_XVC_calibration_offset");

    /* ZDR */
    param->ZDR_C_calibration_offset  = RNC_GetConfigDouble (calfile, "ZED_C_calibration_offset");
    param->ZDR_CP_calibration_offset = RNC_GetConfigDouble (calfile, "ZED_CP_calibration_offset");

    /* LDR */
    param->LDR_HC_calibration_offset = RNC_GetConfigDouble (calfile, "LDR_HC_calibration_offset");
    param->LDR_HCP_calibration_offset = RNC_GetConfigDouble (calfile, "LDR_HCP_calibration_offset");
    param->LDR_VC_calibration_offset = RNC_GetConfigDouble (calfile, "LDR_VC_calibration_offset");
    param->LDR_VCP_calibration_offset = RNC_GetConfigDouble (calfile, "LDR_VC_calibration_offset");

    /* range */
    param->range_offset = RNC_GetConfigDouble (calfile, "range_offset");

    printf ("calibration information loaded\n");
}

//---------------------------------------------------------------------
// Routine to find code sync pulse
static inline int
find_sync (uint16_t *Psamps, int ncodes, int num_samples)
{
    int count=4;
    const int thresh=2000;
    int offset=0;

    while ((Psamps[count]<thresh) && (count < ncodes*num_samples))
    {
        count += num_samples;
        offset++;
    }

    if (offset >= ncodes)
    {
        printf ("** WARNING: no sync found, offset : %d\n", offset);
    }
    else
    {
        printf ("Sync found, offset : %d\n", offset);
    }

    return (offset);
}



//--------------------------------------------------------------------
/* make_dmux_table generates the lookup table that is used to extract
 * channels from the DMA buffer
 * IN:  channels  the number of channels
 * OUT: nowt */
static inline void make_dmux_table (int channels)
{
    if (channels == 4)
    {
        dmux_table[0] = 0;
        dmux_table[1] = 2;
        dmux_table[2] = 1;
        dmux_table[3] = 3;
        dmux_table[4] = 3; //0xFFFFFFFF; /* channels 4 to 7 do not exist in the eight channels system */
        dmux_table[5] = 3; //0xFFFFFFFF;
	dmux_table[6] = 3; //0xFFFFFFFF;
        dmux_table[7] = 3; //0xFFFFFFFF;
    }
    if (channels == 8)
    {
        dmux_table[0] = 0;
        dmux_table[1] = 4;
        dmux_table[2] = 1;
        dmux_table[3] = 5;
        dmux_table[4] = 2;
        dmux_table[5] = 6;
        dmux_table[6] = 3;
        dmux_table[7] = 7;
    }
}



/*---------------------------------------------------------------------------*
 * p_interp: Polynomial interpolation                                        *
 *---------------------------------------------------------------------------*/

static inline float
p3_interpf_1 (float y0,
	      float y1,
	      float y2,
	      float y3)
{
    return (y0 - y2) * 0.3125 + y1 * 0.9375 + y3 * 0.0625;
}

static inline float
p3_interpf_2 (float y0,
	      float y1,
	      float y2,
	      float y3)
{
    return (y1 + y2) * 0.5625 - (y0 + y3) * 0.0625;
}

static inline float
p3_interpf_3 (float y0,
	      float y1,
	      float y2,
	      float y3)
{
    return y0 * 0.0625 + y2 * 0.9375 + (y3 - y1) * 0.3125;
}

//-------------------------------------------------------------------
//     calculate correlation coefficient using power time series
/* Calculate the correlation coefficient for H and V samples.
 * POST: 'h' and 'v' are unaltered.
 */
static inline float
corrCoeff (float * h,
	   float * v,
	   int     n)
{
    double num   = 0.0;
    double meanH = 0.0;
    double meanV = 0.0;
    double sumH2 = 0.0;
    double sumV2 = 0.0;
    float  hh;
    float  vv;
    int    i;

    for (i = 0; i < n; i++)
    {
	meanH += h[i];
	meanV += v[i];
    }
    meanH /= n;
    meanV /= n;

    for (i = 0; i < n; i++)
    {
	hh     = h[i] - meanH;
	vv     = v[i] - meanV;
	num   += hh * vv;
	sumH2 += hh * hh;
	sumV2 += vv * vv;
    }

    return num / sqrt (sumH2 * sumV2);
}

//-------------------------------------------------------------------
//     double length of time series and interpolate
/* Double length of series with third degree polynomial interpolation.
 * Input series (x) length n, output series (xx) length 2n-1
 */
static inline void
double_interp (float * x,
	       float * xx,
	       int     n)
{
    int i;
    int j;

    xx[0] = x[0];
    /* Interpolate for the first point between x[0] and x[1]. */
    xx[1] = p3_interpf_1 (x[0], x[1], x[2], x[3]);

    /* In this loop we want to interpolate each 4 points to get
     * the mid-point between the second and third points.
     */
    for (i = 1, j = 2; i < n - 2; i++, j += 2)
    {
	xx[j+0] = x[i];
	xx[j+1] = p3_interpf_2 (x[i-1], x[i], x[i+1], x[i+2]);
    }
    j = n << 1;
    xx[j-4] = x[n-2];
    /* Interpolate for the last point between h[n - 2] and h[n - 1]. */
    xx[j-3] = p3_interpf_3 (x[n-4], x[n-3], x[n-2], x[n-1]);
    xx[j-2] = x[n-1];
}



//-------------------------------------------------------------------
//     interpolate and calculate offset correlation coefficient
/* Calculate the correlation coefficient for H and V samples.
 * Third degree polynomial interpolation.
 * POST: 'h' and 'v' are unaltered.
 */
static inline float
corrCoeffPoly3 (float * h,
		float * v,
		int     n)
{
    double num   = 0.0;
    double sumH2 = 0.0;
    double sumV2 = 0.0;
    float  hh;
    float  vv;
    double meanH = 0.0;
    double meanV = 0.0;

    int i;

    for (i = 0; i < n; i++)
    {
	meanH += h[i];
	meanV += v[i];
    }
    meanH /= n;
    meanV /= n;

    /* Interpolate for the first point between h[0] and h[1]. */
    hh     = p3_interpf_1 (h[0], h[1], h[2], h[3]) - meanH;
    vv     = v[0] - meanV;
    num   += hh * vv;
    sumH2 += hh * hh;
    sumV2 += vv * vv;

    /* In this loop we want to interpolate each 4 points to get
     * the mid-point between the second and third points.
     */
    for (i = 2; i < n-2; i += 2)
    {
	hh     = p3_interpf_2 (h[i-1], h[i], h[i+1], h[i+2]) - meanH;
	vv     = v[i] - meanV;
	num   += hh * vv;
	sumH2 += hh * hh;
	sumV2 += vv * vv;
    }

    for (i = 2; i < n-2; i += 2)
    {
	vv     = p3_interpf_2 (v[i-2], v[i-1], v[i], v[i+1]) - meanV;
	hh     = h[i] - meanH;
	num   += hh * vv;
	sumH2 += hh * hh;
	sumV2 += vv * vv;
    }

    /* And between v[n - 2] and v[n - 1]. */
    vv     = p3_interpf_3 (v[n-4], v[n-3], v[n-2], v[n-1]) - meanV;
    hh     = h[n-1] - meanH;
    num   += hh * vv;
    sumH2 += hh * hh;
    sumV2 += vv * vv;

    return num / sqrt (sumH2 * sumV2);
}

//========================= M A I N   C O D E =======================
//            [ See disp_help () for command-line options ]
//-------------------------------------------------------------------
int
main (int    argc,
      char * argv[])
{
    int       num_pulses;
    int       amcc_fd = 0;        /* file descriptor for the PCICARD */
    caddr_t   dma_buffer = NULL;  /* size of dma buffer */
    uint16_t * dma_banks[2];      /* pointers to the dma banks */
    int       dma_bank = 0;
    int       proc_bank = 1;
    int       tcount;             /* number of bytes to be transferred during the DMA */
    uint16_t * data;
    int       count, sample;
    int       nspectra;
    int       status;
    long      num_data;
    float *   current_PSD;
    register  int  i, j, k;
    int       temp_int = 0;
    float     HH_moments[RSP_MOMENTS];
    float     HV_moments[RSP_MOMENTS];
    float     VV_moments[RSP_MOMENTS];
    float     VH_moments[RSP_MOMENTS];
    float     HHP_moments[RSP_MOMENTS];
    float     HVP_moments[RSP_MOMENTS];
    float     VVP_moments[RSP_MOMENTS];
    float     VHP_moments[RSP_MOMENTS];
    float  HH_noise_level;
    float  HV_noise_level;
    float  VV_noise_level;
    float  VH_noise_level;
    float  HHP_noise_level;
    float  HVP_noise_level;
    float  VVP_noise_level;
    float  VHP_noise_level;
    int noisegate1, noisegate2;
    time_t    system_time;
    time_t  spectra_time = 0;
    time_t  spectra_rapid_time = 0;
    time_t  temp_time_t;
    struct tm tm;
    char      datestring[25];
    float * mean_vsq_HC; // Used in sigma vbar calculation
    float * mean_vsq_HCP;   // Used in sigma vbar calculation
    float * mean_Zsq_HC; // Used in sigma Zbar calculation
    float * mean_Zsq_HCP;   // Used in sigma Zbar calculation
    float * mean_vsq_VC; // Used in sigma vbar calculation
    float * mean_vsq_VCP;   // Used in sigma vbar calculation
    float * mean_Zsq_VC; // Used in sigma Zbar calculation
    float * mean_Zsq_VCP;   // Used in sigma Zbar calculation
    float tempI_HV;
    float tempQ_HV;
    float mod_HH;
    float mod_VV;
    float tempr, tempz;
    float * h, * v, * hh, * vv, * hp, * vp, * hhp, * vvp;
    float * xx, * yy, * xxp, *yyp;

    short int   xcodes[32][64];

    URC_ScanStruct scan;
    RNC_DimensionStruct dimensions;
    PolPSDStruct * PSD;
    IQStruct IQStruct;

    /* netCDF file pointer */
    int ncid;
    int spectra_ncid       = -1;
    int spectra_rapid_ncid = -1;
    int PSD_varid[PSD_varidSize];
    int PSD_rapid_varid[RapidPSD_varidSize];
    int file_stateid = 0;

    int marker = 0;
    int uncoded_marker = 0;
    int raw_marker = 0;

    float norm_coded, norm_uncoded;

    /* signal */
    struct sigaction sig_struct;

    int gate_offset; // Offset of coded gates rel. to uncoded

    int pulse_increment;

    int tempIco_H;
    int tempQco_H;
    int tempIco_V;
    int tempQco_V;
    int tempIcr_H;
    int tempQcr_H;
    int tempIcr_V;
    int tempQcr_V;
    int nm;


    // The following are shortcut pointers to the elements of
    // the obs structure
    //
    // ZED
    float * ZED_HC, * ZED_HCP;
    float * ZED_VC, * ZED_VCP;
    // ZED_X
    float * ZED_XHC, * ZED_XHCP;
    float * ZED_XVC, * ZED_XVCP;
    // SNR
    float * SNR_HC, * SNR_HCP;
    float * SNR_VC, * SNR_VCP;
    // SNR_X
    float * SNR_XHC, * SNR_XHCP;
    float * SNR_XVC, * SNR_XVCP;
    // VEL
    float * VEL_HC, * VEL_HCP;
    float * VEL_VC, * VEL_VCP;
    // VEL D
    float * VEL_HCD, * VEL_HCDP;
    float * VEL_VCD, * VEL_VCDP;
    // ZED D
    float * ZED_HCD, * ZED_HCDP;
    float * ZED_VCD, * ZED_VCDP;
    // SPW
    float * SPW_HC, * SPW_HCP;
    float * SPW_VC, * SPW_VCP;
    // ZDR
    float * ZDR_C, * ZDR_CP;
    // LDR
    float * LDR_HC, * LDR_HCP;
    float * LDR_VC, * LDR_VCP;
    // VEL COS and SIN
    float * VEL_HC_COS, * VEL_HC_SIN, * VEL_HCP_COS, * VEL_HCP_SIN;
    float * VEL_VC_COS, * VEL_VC_SIN, * VEL_VCP_COS, * VEL_VCP_SIN;
    // PHIDP
    float * PHIDP_C, * PHIDP_CP;
    float * real_PHIDP_C, * imag_PHIDP_C, * real_PHIDP_CP, * imag_PHIDP_CP;
    // RHOHV
    float * RHOHV_C, * RHOHV_CP, * RHOHV_P, * RHOHV_PP, * RHOHV_IP, * RHOHV_IPP;
    // zero I and Q
    float * I_UNCH, * Q_UNCH, * I_UNCV, * Q_UNCV;
    float * I_CODH, * Q_CODH, * I_CODV, * Q_CODV;
    float * NPC_H, * NPC_V;
    float * NPC_HP, * NPC_VP;

    float wi;     // Individual weighting value
    float * uncoded_sum_wi; // Sum of weighting values
    float * coded_sum_wi; // Sum of weighting values
    float tot_pulses;

    int pos1, pos2;
    int tot_n_avg, gate;

    int obtain_index;
    int store_index;

    int   sync_pulse_detected;

    long int * I_coded_copolar_H;
    long int * I_coded_copolar_V;
    long int * Q_coded_copolar_H;
    long int * Q_coded_copolar_V;
    long int * I_coded_crosspolar_H;
    long int * I_coded_crosspolar_V;
    long int * Q_coded_crosspolar_H;
    long int * Q_coded_crosspolar_V;
    long int * I_uncoded_copolar_H;
    long int * I_uncoded_copolar_V;
    long int * Q_uncoded_copolar_H;
    long int * Q_uncoded_copolar_V;
    long int * I_uncoded_crosspolar_H;
    long int * I_uncoded_crosspolar_V;
    long int * Q_uncoded_crosspolar_H;
    long int * Q_uncoded_crosspolar_V;
    uint16_t * I_uncorr_coded_copolar_H;
    uint16_t * I_uncorr_coded_copolar_V;
    uint16_t * Q_uncorr_coded_copolar_H;
    uint16_t * Q_uncorr_coded_copolar_V;
    uint16_t * I_uncorr_coded_crosspolar_H;
    uint16_t * I_uncorr_coded_crosspolar_V;
    uint16_t * Q_uncorr_coded_crosspolar_H;
    uint16_t * Q_uncorr_coded_crosspolar_V;
    uint16_t * I_raw_copolar;
    uint16_t * Q_raw_copolar;
    uint16_t * I_raw_crosspolar;
    uint16_t * Q_raw_crosspolar;
    uint16_t * H_not_V;

    int horizontal_first;

    int  collect_spectra_now;
    int collect_spectra_rapid_now;

    RSP_ParamStruct param, paramCoded, paramUncoded;
    RSP_ComplexType * timeseries;
    RSP_ObservablesStruct obs;
    RSP_ObservablesStruct PSD_obs;
    RSP_ObservablesStruct PSD_RAPID_obs;

    RSP_PeakStruct * HH_peaks;
    RSP_PeakStruct * HV_peaks;
    RSP_PeakStruct * VV_peaks;
    RSP_PeakStruct * VH_peaks;
    RSP_PeakStruct * HHP_peaks;
    RSP_PeakStruct * HVP_peaks;
    RSP_PeakStruct * VVP_peaks;
    RSP_PeakStruct * VHP_peaks;

    static REL_SerialMessageStruct clinometermsg;

    fftw_complex * in, * inx;
    fftw_plan p_coded, p_uncoded;

    /* time variables */
    struct timeval          tv;
    struct timezone         tz;


    disp_welcome_message ();


    //---------------------------
    // Set up the signal handler
    //---------------------------
    /* Set up the sig_struct variable */
    sig_struct.sa_handler = sig_handler;
    sigemptyset (&sig_struct.sa_mask);
    sig_struct.sa_flags = 0;
    /* Install signal handler and check for error */
    if (sigaction (SIGINT, &sig_struct, NULL) != 0)
    {
	perror ("Error installing signal handler\n");
	exit (1);
    }



    //------------------------------
    // Parse command line arguments
    //------------------------------
    if (parseargs (argc, argv, &scan) == 1)
	exit (1);


    //------------------------
    // Read radar config file
    //------------------------
    get_config (CONFIG_FILE, &paramCoded,   &scan, 1); // Do it for coded pulses
    get_config (CONFIG_FILE, &paramUncoded, &scan, 0); // Do it for uncoded pulses

    /* match the scan min_angle and max_angle to the measured elevation angle */
    /* initialise the serial port */
    status = REL_InitialiseSerialMessage (CLINOMETERMESSAGE_PORT);
    /* a slight pause to allow things to settle */
    sleep (1);
    if (status != 0)
    {
	printf ("Detected a problem with initialising the serial port\n");
    }
    else
    {
	status = REL_ReadSerialMessage (&clinometermsg);
	if (status == 0)
	{
	    /* we need to apply a 90 deg elevation offset */
	    clinometermsg.el = clinometermsg.el + 90.0;
	    printf ("elevation angle : %f (%d)\n", clinometermsg.el, temp_int);
	    scan.min_angle = clinometermsg.el;
	    scan.max_angle = clinometermsg.el;
	}
    }

    // Read calibration file
    // We will eventaully need separate cal. files for coded and uncoded
    get_cal (&paramCoded,   CAL_FILE);
    get_cal (&paramUncoded, CAL_FILE);

    //------------------------------------
    // Initialise RSP parameter structure
    //------------------------------------
    paramCoded.prf   /= (paramCoded.num_interleave * paramCoded.num_tx_pol);  // Calculate effective PRF
    paramUncoded.prf /= (paramCoded.num_interleave * paramCoded.num_tx_pol);
    RSP_InitialiseParams (&paramCoded);   // This param is used for coded pulses
    RSP_InitialiseParams (&paramUncoded); // This param is used for uncoded pulses

    //// This botches the DAQ to oversample
    //paramCoded.oversample_ratio   = 2;
    //paramUncoded.oversample_ratio = 2;

    printf ("Parameters for Coded mode:\n");
    RSP_DisplayParams (&paramCoded);
    printf ("Parameters for Uncoded mode:\n");
    RSP_DisplayParams (&paramUncoded);

    param = paramUncoded; // This param is used for params that are the same for both modes

    pulse_increment = 2;  // This is because pulses are interleaved

    // Sample extra pulses at end so that we have entire code sequence
    num_pulses = (int)(paramCoded.pulses_per_daq_cycle + (paramCoded.number_of_codes * paramCoded.num_interleave * paramCoded.num_tx_pol));
    tcount     = param.spectra_averaged * num_pulses * param.samples_per_pulse * param.ADC_channels * sizeof (uint16_t);

    // Number of data points to allocate per data stream
    num_data = paramCoded.pulses_per_daq_cycle * paramCoded.samples_per_pulse;

    I_raw_copolar        = malloc (sizeof (uint16_t) * num_pulses * param.samples_per_pulse);
    Q_raw_copolar        = malloc (sizeof (uint16_t) * num_pulses * param.samples_per_pulse);
    I_raw_crosspolar     = malloc (sizeof (uint16_t) * num_pulses * param.samples_per_pulse);
    Q_raw_crosspolar     = malloc (sizeof (uint16_t) * num_pulses * param.samples_per_pulse);
    H_not_V              = malloc (sizeof (uint16_t) * num_pulses * param.samples_per_pulse);


    // Allocate memory for coded and uncoded data streams
    I_uncorr_coded_copolar_H      = malloc (sizeof (uint16_t) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    Q_uncorr_coded_copolar_H      = malloc (sizeof (uint16_t) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    I_uncorr_coded_copolar_V      = malloc (sizeof (uint16_t) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    Q_uncorr_coded_copolar_V      = malloc (sizeof (uint16_t) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    I_coded_copolar_H             = malloc (sizeof (long int) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    Q_coded_copolar_H             = malloc (sizeof (long int) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    I_coded_copolar_V             = malloc (sizeof (long int) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    Q_coded_copolar_V             = malloc (sizeof (long int) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    IQStruct.I_coded_copolar_H    = malloc (sizeof (long int) * paramCoded.samples_per_pulse * paramCoded.nfft * param.spectra_averaged);
    IQStruct.Q_coded_copolar_H    = malloc (sizeof (long int) * paramCoded.samples_per_pulse * paramCoded.nfft * param.spectra_averaged);
    IQStruct.I_coded_copolar_V    = malloc (sizeof (long int) * paramCoded.samples_per_pulse * paramCoded.nfft * param.spectra_averaged);
    IQStruct.Q_coded_copolar_V    = malloc (sizeof (long int) * paramCoded.samples_per_pulse * paramCoded.nfft * param.spectra_averaged);

    I_uncorr_coded_crosspolar_H   = malloc (sizeof (uint16_t) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    Q_uncorr_coded_crosspolar_H   = malloc (sizeof (uint16_t) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    I_uncorr_coded_crosspolar_V   = malloc (sizeof (uint16_t) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    Q_uncorr_coded_crosspolar_V   = malloc (sizeof (uint16_t) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    I_coded_crosspolar_H          = malloc (sizeof (long int) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    Q_coded_crosspolar_H          = malloc (sizeof (long int) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    I_coded_crosspolar_V          = malloc (sizeof (long int) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    Q_coded_crosspolar_V          = malloc (sizeof (long int) * num_data / paramCoded.num_interleave / paramCoded.num_tx_pol);
    IQStruct.I_coded_crosspolar_H = malloc (sizeof (long int) * paramCoded.samples_per_pulse * paramCoded.nfft * param.spectra_averaged);
    IQStruct.Q_coded_crosspolar_H = malloc (sizeof (long int) * paramCoded.samples_per_pulse * paramCoded.nfft * param.spectra_averaged);
    IQStruct.I_coded_crosspolar_V = malloc (sizeof (long int) * paramCoded.samples_per_pulse * paramCoded.nfft * param.spectra_averaged);
    IQStruct.Q_coded_crosspolar_V = malloc (sizeof (long int) * paramCoded.samples_per_pulse * paramCoded.nfft * param.spectra_averaged);

    num_data = paramUncoded.pulses_per_daq_cycle * paramUncoded.samples_per_pulse;

    I_uncoded_copolar_H          = malloc (sizeof (long int) * num_data / paramUncoded.num_interleave / paramUncoded.num_tx_pol);
    Q_uncoded_copolar_H          = malloc (sizeof (long int) * num_data / paramUncoded.num_interleave / paramUncoded.num_tx_pol);
    I_uncoded_copolar_V          = malloc (sizeof (long int) * num_data / paramUncoded.num_interleave / paramUncoded.num_tx_pol);
    Q_uncoded_copolar_V          = malloc (sizeof (long int) * num_data / paramUncoded.num_interleave / paramUncoded.num_tx_pol);
    IQStruct.I_uncoded_copolar_H = malloc (sizeof (uint16_t) * paramUncoded.samples_per_pulse * paramUncoded.nfft * param.spectra_averaged);
    IQStruct.Q_uncoded_copolar_H = malloc (sizeof (uint16_t) * paramUncoded.samples_per_pulse * paramUncoded.nfft * param.spectra_averaged);
    IQStruct.I_uncoded_copolar_V = malloc (sizeof (uint16_t) * paramUncoded.samples_per_pulse * paramUncoded.nfft * param.spectra_averaged);
    IQStruct.Q_uncoded_copolar_V = malloc (sizeof (uint16_t) * paramUncoded.samples_per_pulse * paramUncoded.nfft * param.spectra_averaged);

    I_uncoded_crosspolar_H          = malloc (sizeof (long int) * num_data / paramUncoded.num_interleave / paramUncoded.num_tx_pol);
    Q_uncoded_crosspolar_H          = malloc (sizeof (long int) * num_data / paramUncoded.num_interleave / paramUncoded.num_tx_pol);
    I_uncoded_crosspolar_V          = malloc (sizeof (long int) * num_data / paramUncoded.num_interleave / paramUncoded.num_tx_pol);
    Q_uncoded_crosspolar_V          = malloc (sizeof (long int) * num_data / paramUncoded.num_interleave / paramUncoded.num_tx_pol);
    IQStruct.I_uncoded_crosspolar_H = malloc (sizeof (uint16_t) * paramUncoded.samples_per_pulse * paramUncoded.nfft * param.spectra_averaged);
    IQStruct.Q_uncoded_crosspolar_H = malloc (sizeof (uint16_t) * paramUncoded.samples_per_pulse * paramUncoded.nfft * param.spectra_averaged);
    IQStruct.I_uncoded_crosspolar_V = malloc (sizeof (uint16_t) * paramUncoded.samples_per_pulse * paramUncoded.nfft * param.spectra_averaged);
    IQStruct.Q_uncoded_crosspolar_V = malloc (sizeof (uint16_t) * paramUncoded.samples_per_pulse * paramUncoded.nfft * param.spectra_averaged);

    HH_peaks  = calloc (param.num_peaks, sizeof (RSP_PeakStruct));
    HV_peaks  = calloc (param.num_peaks, sizeof (RSP_PeakStruct));
    VH_peaks  = calloc (param.num_peaks, sizeof (RSP_PeakStruct));
    VV_peaks  = calloc (param.num_peaks, sizeof (RSP_PeakStruct));
    HHP_peaks = calloc (param.num_peaks, sizeof (RSP_PeakStruct));
    HVP_peaks = calloc (param.num_peaks, sizeof (RSP_PeakStruct));
    VHP_peaks = calloc (param.num_peaks, sizeof (RSP_PeakStruct));
    VVP_peaks = calloc (param.num_peaks, sizeof (RSP_PeakStruct));

    in        = fftw_malloc (sizeof (fftw_complex) * param.nfft);
    inx       = fftw_malloc (sizeof (fftw_complex) * param.nfft);
    p_coded   = fftw_plan_dft_1d (paramCoded.nfft,   in, in, FFTW_FORWARD, FFTW_ESTIMATE);
    p_uncoded = fftw_plan_dft_1d (paramUncoded.nfft, in, in, FFTW_FORWARD, FFTW_ESTIMATE);

    timeseries = malloc (param.nfft * sizeof (RSP_ComplexType));

    current_PSD = malloc (param.npsd * sizeof (float));
    PSD         = calloc (param.samples_per_pulse, sizeof (PolPSDStruct));
    for (j = 0; j < param.samples_per_pulse; j++)
    {
	PSD[j].HH  = malloc (param.npsd * sizeof (float)); // H not coded
	PSD[j].HV  = malloc (param.npsd * sizeof (float)); // H not coded
	PSD[j].HHP = malloc (param.npsd * sizeof (float)); // H coded
	PSD[j].HVP = malloc (param.npsd * sizeof (float)); // H coded
	PSD[j].VV  = malloc (param.npsd * sizeof (float)); // V not coded
	PSD[j].VH  = malloc (param.npsd * sizeof (float)); // V not coded
	PSD[j].VVP = malloc (param.npsd * sizeof (float)); // V coded
	PSD[j].VHP = malloc (param.npsd * sizeof (float)); // V coded

    }
    mean_vsq_HCP = malloc (paramCoded.samples_per_pulse   * sizeof (float));
    mean_vsq_HC  = malloc (paramUncoded.samples_per_pulse * sizeof (float));
    mean_Zsq_HCP = malloc (paramCoded.samples_per_pulse   * sizeof (float));
    mean_Zsq_HC  = malloc (paramUncoded.samples_per_pulse * sizeof (float));
    mean_vsq_VCP = malloc (paramCoded.samples_per_pulse   * sizeof (float));
    mean_vsq_VC  = malloc (paramUncoded.samples_per_pulse * sizeof (float));
    mean_Zsq_VCP = malloc (paramCoded.samples_per_pulse   * sizeof (float));
    mean_Zsq_VC  = malloc (paramUncoded.samples_per_pulse * sizeof (float));

    coded_sum_wi   = malloc (paramCoded.samples_per_pulse   * sizeof (float));
    uncoded_sum_wi = malloc (paramUncoded.samples_per_pulse * sizeof (float));



    VEL_HC_COS  = malloc (paramUncoded.samples_per_pulse * sizeof (float));
    VEL_HC_SIN  = malloc (paramUncoded.samples_per_pulse * sizeof (float));
    VEL_VC_COS  = malloc (paramUncoded.samples_per_pulse * sizeof (float));
    VEL_VC_SIN  = malloc (paramUncoded.samples_per_pulse * sizeof (float));
    VEL_HCP_COS = malloc (paramCoded.samples_per_pulse   * sizeof (float));
    VEL_HCP_SIN = malloc (paramCoded.samples_per_pulse   * sizeof (float));
    VEL_VCP_COS = malloc (paramCoded.samples_per_pulse   * sizeof (float));
    VEL_VCP_SIN = malloc (paramCoded.samples_per_pulse   * sizeof (float));

    real_PHIDP_C  = malloc (paramUncoded.samples_per_pulse * sizeof (float));
    imag_PHIDP_C  = malloc (paramUncoded.samples_per_pulse * sizeof (float));
    real_PHIDP_CP = malloc (paramCoded.samples_per_pulse   * sizeof (float));
    imag_PHIDP_CP = malloc (paramCoded.samples_per_pulse   * sizeof (float));

    h   = malloc (paramUncoded.samples_per_pulse * sizeof (float));
    v   = malloc (paramUncoded.samples_per_pulse * sizeof (float));
    hh  = malloc (paramUncoded.samples_per_pulse * sizeof (float) * 2 - 1);
    vv  = malloc (paramUncoded.samples_per_pulse * sizeof (float) * 2 - 1);
    xx  = malloc (paramUncoded.samples_per_pulse * sizeof (float) * 2 - 1);
    yy  = malloc (paramUncoded.samples_per_pulse * sizeof (float) * 2 - 1);
    hp  = malloc (paramCoded.samples_per_pulse   * sizeof (float));
    vp  = malloc (paramCoded.samples_per_pulse   * sizeof (float));
    hhp = malloc (paramCoded.samples_per_pulse   * sizeof (float) * 2 - 1);
    vvp = malloc (paramCoded.samples_per_pulse   * sizeof (float) * 2 - 1);
    xxp = malloc (paramCoded.samples_per_pulse   * sizeof (float) * 2 - 1);
    yyp = malloc (paramCoded.samples_per_pulse   * sizeof (float) * 2 - 1);


    norm_coded   = 1.0 / paramCoded.Wss;
    norm_uncoded = 1.0 / paramUncoded.Wss;

    //----------------------------
    // Initialise RSP Observables
    //----------------------------
    RSP_ObsInit (&obs);
    RSP_ObsInit (&PSD_obs);
    RSP_ObsInit (&PSD_RAPID_obs);
    // Last argument determines whether the parameter will be recorded or not
    // Horizontal
    ZED_HC   = RSP_ObsNew (&obs, "ZED_HC",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZED_HC"));
    ZED_XHC  = RSP_ObsNew (&obs, "ZED_XHC",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZED_XHC"));
    ZED_HCP  = RSP_ObsNew (&obs, "ZED_HCP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZED_HCP"));
    ZED_XHCP = RSP_ObsNew (&obs, "ZED_XHCP", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZED_XHCP"));

    SNR_HC   = RSP_ObsNew (&obs, "SNR_HC",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "SNR_HC"));
    SNR_XHC  = RSP_ObsNew (&obs, "SNR_XHC",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "SNR_XHC"));
    SNR_HCP  = RSP_ObsNew (&obs, "SNR_HCP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "SNR_HCP"));
    SNR_XHCP = RSP_ObsNew (&obs, "SNR_XHCP", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "SNR_XHCP"));

    VEL_HC   = RSP_ObsNew (&obs, "VEL_HC",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "VEL_HC"));
    VEL_HCP  = RSP_ObsNew (&obs, "VEL_HCP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "VEL_HCP"));

    SPW_HC   = RSP_ObsNew (&obs, "SPW_HC",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "SPW_HC"));
    SPW_HCP  = RSP_ObsNew (&obs, "SPW_HCP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "SPW_HCP"));

    VEL_HCD  = RSP_ObsNew (&obs, "VEL_HCD",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "VEL_HCD"));
    ZED_HCD  = RSP_ObsNew (&obs, "ZED_HCD",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZED_HCD"));
    VEL_HCDP = RSP_ObsNew (&obs, "VEL_HCDP", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "VEL_HCDP"));
    ZED_HCDP = RSP_ObsNew (&obs, "ZED_HCDP", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZED_HCDP"));

    // Vertical
    ZED_VC   = RSP_ObsNew (&obs, "ZED_VC",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZED_VC"));
    ZED_XVC  = RSP_ObsNew (&obs, "ZED_XVC",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZED_XVC"));
    ZED_VCP  = RSP_ObsNew (&obs, "ZED_VCP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZED_VCP"));
    ZED_XVCP = RSP_ObsNew (&obs, "ZED_XVCP", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZED_XVCP"));

    SNR_VC   = RSP_ObsNew (&obs, "SNR_VC",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "SNR_VC"));
    SNR_XVC  = RSP_ObsNew (&obs, "SNR_XVC",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "SNR_XVC"));
    SNR_VCP  = RSP_ObsNew (&obs, "SNR_VCP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "SNR_VCP"));
    SNR_XVCP = RSP_ObsNew (&obs, "SNR_XVCP", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "SNR_XVCP"));

    VEL_VC   = RSP_ObsNew (&obs, "VEL_VC",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "VEL_VC"));
    VEL_VCP  = RSP_ObsNew (&obs, "VEL_VCP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "VEL_VCP"));

    SPW_VC   = RSP_ObsNew (&obs, "SPW_VC",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "SPW_VC"));
    SPW_VCP  = RSP_ObsNew (&obs, "SPW_VCP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "SPW_VCP"));

    VEL_VCD  = RSP_ObsNew (&obs, "VEL_VCD",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "VEL_VCD"));
    ZED_VCD  = RSP_ObsNew (&obs, "ZED_VCD",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZED_VCD"));
    VEL_VCDP = RSP_ObsNew (&obs, "VEL_VCDP", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "VEL_VCDP"));
    ZED_VCDP = RSP_ObsNew (&obs, "ZED_VCDP", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZED_VCDP"));

    NPC_H    = RSP_ObsNew (&obs, "NPC_H",    param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "NPC_H"));
    NPC_V    = RSP_ObsNew (&obs, "NPC_V",    param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "NPC_V"));
    NPC_HP   = RSP_ObsNew (&obs, "NPC_HP",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "NPC_HP"));
    NPC_VP   = RSP_ObsNew (&obs, "NPC_VP",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "NPC_VP"));

    // polarisation independent
    ZDR_C   = RSP_ObsNew (&obs, "ZDR_C",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZDR_C"));
    LDR_HC  = RSP_ObsNew (&obs, "LDR_HC",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "LDR_HC"));
    LDR_VC  = RSP_ObsNew (&obs, "LDR_VC",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "LDR_VC"));
    ZDR_CP  = RSP_ObsNew (&obs, "ZDR_CP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "ZDR_CP"));
    LDR_HCP = RSP_ObsNew (&obs, "LDR_HCP", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "LDR_HCP"));
    LDR_VCP = RSP_ObsNew (&obs, "LDR_VCP", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "LDR_VCP"));

    PHIDP_C   = RSP_ObsNew (&obs, "PHIDP_C",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "PHIDP_C"));
    PHIDP_CP  = RSP_ObsNew (&obs, "PHIDP_CP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "PHIDP_CP"));
    RHOHV_C   = RSP_ObsNew (&obs, "RHOHV_C",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "RHOHV_C"));
    RHOHV_CP  = RSP_ObsNew (&obs, "RHOHV_CP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "RHOHV_CP"));
    RHOHV_P   = RSP_ObsNew (&obs, "RHOHV_P",   param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "RHOHV_P"));
    RHOHV_PP  = RSP_ObsNew (&obs, "RHOHV_PP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "RHOHV_PP"));
    RHOHV_IP  = RSP_ObsNew (&obs, "RHOHV_IP",  param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "RHOHV_IP"));
    RHOHV_IPP = RSP_ObsNew (&obs, "RHOHV_IPP", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "RHOHV_IPP"));

    I_UNCH  = RSP_ObsNew (&obs, "I_UNCH", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "I_UNCH"));
    Q_UNCH  = RSP_ObsNew (&obs, "Q_UNCH", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "Q_UNCH"));
    I_UNCV  = RSP_ObsNew (&obs, "I_UNCV", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "I_UNCV"));
    Q_UNCV  = RSP_ObsNew (&obs, "Q_UNCV", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "Q_UNCV"));
    I_CODH  = RSP_ObsNew (&obs, "I_CODH", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "I_CODH"));
    Q_CODH  = RSP_ObsNew (&obs, "Q_CODH", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "Q_CODH"));
    I_CODV  = RSP_ObsNew (&obs, "I_CODV", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "I_CODV"));
    Q_CODV  = RSP_ObsNew (&obs, "Q_CODV", param.samples_per_pulse, (int)RNC_GetConfigDouble (CONFIG_FILE, "Q_CODV"));


    printf ("Recording observables:");
    for (i = 0; i < obs.n_obs; i++)
    {
	if (obs.record_observable[i] == 1)
	{
	    printf (" %s", obs.name[i]);
	}
    }
    printf ("\n");


    //-------------------------------------------
    // Set up the data acquisition
    //-------------------------------------------
    printf ("** Initialising ISACTRL...\n");
    RDQ_InitialiseISACTRL (num_pulses, param.samples_per_pulse,
			   param.clock_divfactor, param.delay_clocks);

    printf ("** Initialising PCICARD...\n");
    amcc_fd = RDQ_InitialisePCICARD_New (&dma_buffer, DMA_BUFFER_SIZE);

    // Initialise pointers to DMA banks
    dma_banks[ 0 ] = (uint16_t *)dma_buffer;
    dma_banks[ 1 ] = (uint16_t *)(dma_buffer + (DMA_BUFFER_SIZE/2));

    make_dmux_table (param.ADC_channels);

    printf ("** Starting acquisition...\n");

    RDQ_StartAcquisition (amcc_fd, dma_bank, (short *)(dma_banks[dma_bank]), tcount);

    /* setup the netCDF file */
    ncid = RNC_OpenNetcdfFile (GetRadarName (COPERNICUS),
			       GetSpectraName (COPERNICUS),
			       scan.date, NULL,
			       GetScanTypeName (scan.scanType),
			       GetSpectraExtension (COPERNICUS), "raw");
    RNC_SetupDimensions (ncid, &param, &dimensions);
    RNC_SetupGlobalAttributes (ncid, COPERNICUS, &scan, &param, argc, argv);
    printf ("netCDF : global attributes have been defined.\n");
    RNC_SetupPulse_Compression_Code (ncid, &param);
    file_stateid = RNC_SetupFile_State (ncid);
    RNC_SetupStaticVariables (ncid, &param);
    printf ("netCDF : static variables have been defined\n");
    RNC_SetupRange (ncid, &param, &dimensions);
    RNC_SetupDynamicVariables (ncid, COPERNICUS, &scan, &param, &dimensions, &obs);
    printf ("netCDF : dynamic variables have been defined\n");
    /* change the mode of netCDF from define to data */
    status = nc_enddef (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    // Set up spectral dump file
    if (param.dump_spectra != 0)
    {
	printf ("setup spectra recording file\n");
	spectra_ncid = RNC_OpenNetcdfFile
	    (GetRadarName (COPERNICUS_SPECTRA),
	     GetSpectraName (COPERNICUS_SPECTRA),
	     scan.date, NULL,
	     GetScanTypeName (scan.scanType),
	     GetSpectraExtension (COPERNICUS_SPECTRA), "raw");
	RNC_SetupDimensions (spectra_ncid, &param, &dimensions);
	RNC_SetupGlobalAttributes (spectra_ncid, COPERNICUS, &scan, &param, argc, argv);
	RNC_SetupPulse_Compression_Code (spectra_ncid, &param);
	RNC_SetupStaticVariables (spectra_ncid, &param);
	RNC_SetupRange (spectra_ncid, &param, &dimensions);
	RNC_SetupDynamicVariables (spectra_ncid, COPERNICUS_SPECTRA, &scan, &param, &dimensions, &PSD_obs);
	RNC_SetupLogPSDVariables (spectra_ncid, COPERNICUS_CODED_SPECTRA, &param, &dimensions, PSD_varid);
	/* change the mode of netCDF from define to data */
	status = nc_enddef (spectra_ncid);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
	time (&temp_time_t);
        spectra_time = param.dump_spectra * (floorl (temp_time_t / param.dump_spectra) + 1);
    }

    if (param.dump_spectra_rapid != 0)
    {
        /* make sure bin_ray_number is 0 */
        PSD_RAPID_obs.bin_ray_number = 0;
	PSD_RAPID_obs.ray_number = 0;

	spectra_rapid_ncid = RNC_OpenNetcdfFile
	    (GetRadarName (COPERNICUS_SPECTRA_RAPID),
	     GetSpectraName (COPERNICUS_SPECTRA_RAPID),
	     scan.date, NULL,
	     GetScanTypeName (scan.scanType),
	     GetSpectraExtension (COPERNICUS_SPECTRA_RAPID), "raw");
        RNC_SetupRapidLogPSDDimensions (spectra_rapid_ncid, COPERNICUS_SPECTRA_RAPID, &param, &dimensions);
        RNC_SetupGlobalAttributes (spectra_rapid_ncid, COPERNICUS, &scan, &param, argc, argv);
	RNC_SetupStaticVariables (spectra_rapid_ncid, &param);
	RNC_SetupRange (spectra_rapid_ncid, &param, &dimensions);
        RNC_SetupDynamicVariables (spectra_rapid_ncid, COPERNICUS_SPECTRA_RAPID, &scan, &param, &dimensions, &PSD_RAPID_obs);
        RNC_SetupLogPSDVariables (spectra_rapid_ncid, COPERNICUS_SPECTRA_RAPID, &param, &dimensions, PSD_rapid_varid);
        /* change the mode of netCDF from define to data */
        status = nc_enddef (spectra_rapid_ncid);
        if (status != NC_NOERR) check_netcdf_handle_error (status);

	time (&temp_time_t);
	spectra_rapid_time = param.dump_spectra_rapid * (floorl (temp_time_t / param.dump_spectra_rapid) + 1);
    }



    //------------------------------------------------------------------
    // "Oversample" the pulse code templates to match the sampling rate
    //------------------------------------------------------------------
    printf ("number_of_codes: %d\n", paramCoded.number_of_codes);
    printf ("num_interleave:  %d\n", paramCoded.num_interleave);
    for (i = 0; i < paramCoded.number_of_codes * paramCoded.num_interleave * paramCoded.num_tx_pol; i++)
    {
	printf ("sequence: %d\n", i);
	for (j = 0; j < paramCoded.code_length; j++)
	    printf ("%i ", paramCoded.codes[i][j]);
	printf ("\n");
	RSP_Oversample (paramCoded.codes[i], xcodes[i], paramCoded.code_length, paramCoded.oversample_ratio);
	for (j = 0; j < paramCoded.code_length * paramCoded.oversample_ratio; j++)
	    printf ("%i ", xcodes[i][j]);
	printf ("\n");
    }

    gate_offset = paramCoded.code_length * paramCoded.oversample_ratio - 1;


    // THIS IS THE START OF THE OUTER RAY LOOP
    while (exit_now == 0)
    {
	printf ("<< PRESS CTRL-C TO EXIT >>\n");

	/* get time of day */
	gettimeofday (&tv, &tz);
	gmtime_r (&tv.tv_sec, &tm);
	printf ("System time: %s\n", asctime (&tm));
	obs.year        = tm.tm_year + 1900;
	obs.month       = tm.tm_mon  + 1;
	obs.day         = tm.tm_mday;
	obs.hour        = tm.tm_hour;
	obs.minute      = tm.tm_min;
	obs.second      = tm.tm_sec;
	obs.centisecond = tv.tv_usec / 10000U;
	sprintf (datestring, "%04d/%02d/%02d %02d:%02d:%02d.%02d",
		 obs.year, obs.month, obs.day,
		 obs.hour, obs.minute, obs.second, obs.centisecond);
	printf ("Date string: %s", datestring);

	obs.azimuth           = scan.scan_angle;
	PSD_obs.azimuth       = obs.azimuth;
	PSD_RAPID_obs.azimuth = obs.azimuth;
	/* read the elevation angle from the clionometer */
	status = REL_ReadSerialMessage (&clinometermsg);
	if (status == 0)
	{
	    /* we need to apply a 90 deg elevation offset */
	    clinometermsg.el = clinometermsg.el + 90.0;
	    printf ("elevation angle : %f (%d)\n", clinometermsg.el, temp_int);
	    obs.elevation = clinometermsg.el;
	}
	else
	{
	    obs.elevation = -999;
	}
	PSD_obs.elevation       = obs.elevation;
	PSD_RAPID_obs.elevation = obs.elevation;

	// Initialise observables to zero
	// (Needed for moments averaging)
	for (i = 0; i < obs.n_obs; i++)
	{
	    printf ("Initialising %s\n", obs.name[i]);
	    for (j = 0; j < obs.n_elements [i]; j++)
	    {
		obs.data[i][j] = 0.0f;
	    }
        }
	for (j = 0; j < param.samples_per_pulse; j++)
        {
	    mean_vsq_HC [j] = 0.0f;
	    mean_vsq_HCP[j] = 0.0f;
	    mean_Zsq_HC [j] = 0.0f;
	    mean_Zsq_HCP[j] = 0.0f;
	    mean_vsq_VC [j] = 0.0f;
	    mean_vsq_VCP[j] = 0.0f;
	    mean_Zsq_VC [j] = 0.0f;
	    mean_Zsq_VCP[j] = 0.0f;

	    uncoded_sum_wi[j] = 0.0f;
	    coded_sum_wi  [j] = 0.0f;


	    VEL_HC_COS [j] = 0.0f;
	    VEL_HC_SIN [j] = 0.0f;
	    VEL_HCP_COS[j] = 0.0f;
	    VEL_HCP_SIN[j] = 0.0f;
	    VEL_VC_COS [j] = 0.0f;
	    VEL_VC_SIN [j] = 0.0f;
	    VEL_VCP_COS[j] = 0.0f;
	    VEL_VCP_SIN[j] = 0.0f;

	    real_PHIDP_C [j] = 0.0f;
	    imag_PHIDP_C [j] = 0.0f;
	    real_PHIDP_CP[j] = 0.0f;
	    imag_PHIDP_CP[j] = 0.0f;
	}

	printf ("Done initialising variables...\n");

	// LOOP THROUGH MOMENTS AVERAGING FROM HERE...
	for (nm = 0; nm < param.moments_averaged; nm++)
	{
	    for (j = 0; j < param.samples_per_pulse; j++)
	    {
		register int bin_no;
		// Initialise spectra to zero
		for (bin_no = 0; bin_no < param.npsd; bin_no++)
		{
		    PSD[j].HH [bin_no] = 0.0f;
		    PSD[j].HV [bin_no] = 0.0f;
		    PSD[j].HHP[bin_no] = 0.0f;
		    PSD[j].HVP[bin_no] = 0.0f;
		    PSD[j].VV [bin_no] = 0.0f;
		    PSD[j].VH [bin_no] = 0.0f;
		    PSD[j].VVP[bin_no] = 0.0f;
		    PSD[j].VHP[bin_no] = 0.0f;
		}
	    }

	    memset (dma_banks[proc_bank], -1, tcount);

	    /* waiting for DMA to complete */
	    status = RDQ_WaitForAcquisitionToComplete (amcc_fd);
	    if (status != 0)
		printf ("There was a problem in WaitForAcquisitionToComplete\n");

	    //----------------------------------------------------------------
	    // Swap around the areas used for storing daq and processing from
	    //----------------------------------------------------------------
	    dma_bank  = 1 - dma_bank;
	    proc_bank = 1 - proc_bank;

	    data = dma_banks[proc_bank];
	    RDQ_StartAcquisition2 (amcc_fd, dma_bank, tcount);

	    //---------------------------------
	    // Loop through spectral averaging
	    //---------------------------------
	    collect_spectra_rapid_now = 0;
	    collect_spectra_now       = 0;
	    system_time = time (NULL);
	    if (param.dump_spectra_rapid != 0)
	    {
                if (spectra_rapid_time <= system_time)
		{
		    collect_spectra_rapid_now = 1;
		}

	    }
	    /* write out spectra to netCDF if required */
	    if (param.dump_spectra != 0)
	    {
                if (spectra_time <= system_time)
		{
		    collect_spectra_now = 1;
		}
	    }

	    for (nspectra = 0; nspectra < param.spectra_averaged; nspectra++)
	    {
		printf ("\nAveraging %d of %d spectra:\n", nspectra + 1, param.spectra_averaged);

		//----------------------------------------------------------------
		// Extract data from DMA memory (still interleaved at this stage)
		//----------------------------------------------------------------
		//
		// find out where the sync pulse is
		marker = 0;
		sync_pulse_detected = -1;
		for (i = 0; i < num_pulses; i++)
		{
		    register int count_reg;
		    for (j = 0; j < param.samples_per_pulse; j++)
		    {
			uint16_t Sync;

			count_reg=marker*param.samples_per_pulse+j;
			I_raw_copolar   [count_reg] = GET_CHANNEL (data, CHAN_Ic);
			Q_raw_copolar   [count_reg] = GET_CHANNEL (data, CHAN_Qc);
			I_raw_crosspolar[count_reg] = GET_CHANNEL (data, CHAN_Ix);
			Q_raw_crosspolar[count_reg] = GET_CHANNEL (data, CHAN_Qx);
			H_not_V         [count_reg] = GET_CHANNEL (data, CHAN_H_not_V);
			Sync                        = GET_CHANNEL (data, CHAN_SYNC);
			// if (j==5)
			// {
			//   printf ("Sync pulse = %d at pulse %i *********\n", Sync, i);
			// }
			// This will discard all data before the first sync pulse
			if (sync_pulse_detected == -1 && Sync > 2000)
			{
			    sync_pulse_detected = i;
			    printf ("Sync pulse found at pulse %i **********\n",
				    sync_pulse_detected + 1);
			}
			INC_POINTER (data, param.ADC_channels);
		    }
		    if (sync_pulse_detected > -1)
		    {
			/* marker will keep a track of the number of pulses we have collected */
			marker = marker + 1;
		    }
		    if (marker == param.pulses_per_daq_cycle)
		    {
			/* if we have all the data we need the exit this loop */
			i = num_pulses;
		    }
		}
		if (sync_pulse_detected == -1)
		{
		    printf ("** ERROR: NO SYNC PULSE FOUND!\n");
		    exit (1);
		}
		// printf ("Pulses transferred: %d, samples transferred: %d\n", marker, count_reg + 1);
		printf ("Pulses transferred: %d\n", marker);

		// NB: Pulse sequence should be as follows:
		// Code A : Code A : Uncoded : Uncoded : Code B : Code B : Uncoded : Uncoded : Code A
		// SYNC ^ :        :         :         :        :        :         :         : SYNC ^ : ...
		// H or V : V or H : H or V  : V or H ......

		/* need to sort out this data */
		pos2 = 0;
		/* first of all let us work out if a H or V went out first */
		// for (j = 0; j < 20; j++)
		// {
		//   printf ("H_not_V = %d \n", H_not_V[j]);
		// }
		if (H_not_V[0] > 2048)
		{
		    /* this means first pulse is horizontal */
		    horizontal_first = 1;
		    printf ("The first pulse is horizontal\n");
		}
		else
		{
		    /* this means first pulse is vertical */
		    horizontal_first = 0;
		    printf ("The first pulse is vertical\n");
		}
		raw_marker     = 0;
		marker         = 0;
		uncoded_marker = 0;
		for (i = 0; i < param.pulses_per_daq_cycle / (paramCoded.num_interleave * paramCoded.num_tx_pol); i++)
		{
		    pos1 = marker;
		    if (horizontal_first == 1)
		    {
			for (j = 0; j < param.samples_per_pulse; j++)
			{
			    I_uncorr_coded_copolar_H   [pos1] = I_raw_copolar   [raw_marker];
			    Q_uncorr_coded_copolar_H   [pos1] = Q_raw_copolar   [raw_marker];
			    I_uncorr_coded_crosspolar_H[pos1] = I_raw_crosspolar[raw_marker];
			    Q_uncorr_coded_crosspolar_H[pos1] = Q_raw_crosspolar[raw_marker];
			    raw_marker ++;
			    pos1 ++;
			}
			if (param.num_tx_pol == 2)
			{
			    pos1 = marker;
			    for (j = 0; j < param.samples_per_pulse; j++)
			    {
				I_uncorr_coded_copolar_V   [pos1] = I_raw_copolar   [raw_marker];
				Q_uncorr_coded_copolar_V   [pos1] = Q_raw_copolar   [raw_marker];
				I_uncorr_coded_crosspolar_V[pos1] = I_raw_crosspolar[raw_marker];
				Q_uncorr_coded_crosspolar_V[pos1] = Q_raw_crosspolar[raw_marker];
				raw_marker ++;
				pos1 ++;
			    }
			}
			pos1 = marker;
			for (j = 0; j < param.samples_per_pulse; j++)
			{
			    I_uncoded_copolar_H   [pos1] = I_raw_copolar   [raw_marker];
			    Q_uncoded_copolar_H   [pos1] = Q_raw_copolar   [raw_marker];
			    I_uncoded_crosspolar_H[pos1] = I_raw_crosspolar[raw_marker];
			    Q_uncoded_crosspolar_H[pos1] = Q_raw_crosspolar[raw_marker];
			    raw_marker ++;
			    pos1 ++;
			}
			if (param.num_tx_pol == 2)
			{
			    pos1 = marker;
			    for (j = 0; j < param.samples_per_pulse; j++)
			    {
				I_uncoded_copolar_V   [pos1] = I_raw_copolar   [raw_marker];
				Q_uncoded_copolar_V   [pos1] = Q_raw_copolar   [raw_marker];
				I_uncoded_crosspolar_V[pos1] = I_raw_crosspolar[raw_marker];
				Q_uncoded_crosspolar_V[pos1] = Q_raw_crosspolar[raw_marker];
				raw_marker ++;
				pos1 ++;
			    }
			}
		    }
		    else
		    {
			for (j = 0; j < param.samples_per_pulse; j++)
			{
			    I_uncorr_coded_copolar_V   [pos1] = I_raw_copolar   [raw_marker];
			    Q_uncorr_coded_copolar_V   [pos1] = Q_raw_copolar   [raw_marker];
			    I_uncorr_coded_crosspolar_V[pos1] = I_raw_crosspolar[raw_marker];
			    Q_uncorr_coded_crosspolar_V[pos1] = Q_raw_crosspolar[raw_marker];
			    raw_marker ++;
			    pos1 ++;
			}
			if (param.num_tx_pol == 2)
			{
			    pos1 = marker;
			    for (j = 0; j < param.samples_per_pulse; j++)
			    {
				I_uncorr_coded_copolar_H   [pos1] = I_raw_copolar   [raw_marker];
				Q_uncorr_coded_copolar_H   [pos1] = Q_raw_copolar   [raw_marker];
				I_uncorr_coded_crosspolar_H[pos1] = I_raw_crosspolar[raw_marker];
				Q_uncorr_coded_crosspolar_H[pos1] = Q_raw_crosspolar[raw_marker];
				raw_marker ++;
				pos1 ++;
			    }
			}
			pos1 = marker;
			for (j = 0; j < param.samples_per_pulse; j++)
			{
			    I_uncoded_copolar_V   [pos1] = I_raw_copolar   [raw_marker];
			    Q_uncoded_copolar_V   [pos1] = Q_raw_copolar   [raw_marker];
			    I_uncoded_crosspolar_V[pos1] = I_raw_crosspolar[raw_marker];
			    Q_uncoded_crosspolar_V[pos1] = Q_raw_crosspolar[raw_marker];
			    raw_marker ++;
			    pos1 ++;
			}
			if (param.num_tx_pol == 2)
			{
			    pos1 = marker;
			    for (j = 0; j < param.samples_per_pulse; j++)
			    {
				I_uncoded_copolar_H   [pos1] = I_raw_copolar   [raw_marker];
				Q_uncoded_copolar_H   [pos1] = Q_raw_copolar   [raw_marker];
				I_uncoded_crosspolar_H[pos1] = I_raw_crosspolar[raw_marker];
				Q_uncoded_crosspolar_H[pos1] = Q_raw_crosspolar[raw_marker];
				raw_marker ++;
				pos1 ++;
			    }
			}
		    }
		    marker += param.samples_per_pulse;
		}

		// DECODE PULSES
		// Loop through coded pulses
		count = paramCoded.code_length * paramCoded.oversample_ratio;
		pos1  = 0;
		k     = 0;
		for (i = 0; i < (param.pulses_per_daq_cycle / (paramCoded.num_interleave * paramCoded.num_tx_pol)); i++ )
		{
		    // diagnostic printf below
		    // printf ("i: %d, pos1: %d, k: %d\n", i, pos1, k);
		    /* correlations */
		    RSP_Correlate (&I_uncorr_coded_copolar_H   [pos1], xcodes[k], param.samples_per_pulse, count, &I_coded_copolar_H   [pos1]);
		    RSP_Correlate (&Q_uncorr_coded_copolar_H   [pos1], xcodes[k], param.samples_per_pulse, count, &Q_coded_copolar_H   [pos1]);
		    RSP_Correlate (&I_uncorr_coded_crosspolar_H[pos1], xcodes[k], param.samples_per_pulse, count, &I_coded_crosspolar_H[pos1]);
		    RSP_Correlate (&Q_uncorr_coded_crosspolar_H[pos1], xcodes[k], param.samples_per_pulse, count, &Q_coded_crosspolar_H[pos1]);
		    RSP_Correlate (&I_uncorr_coded_copolar_V   [pos1], xcodes[k], param.samples_per_pulse, count, &I_coded_copolar_V   [pos1]);
		    RSP_Correlate (&Q_uncorr_coded_copolar_V   [pos1], xcodes[k], param.samples_per_pulse, count, &Q_coded_copolar_V   [pos1]);
		    RSP_Correlate (&I_uncorr_coded_crosspolar_V[pos1], xcodes[k], param.samples_per_pulse, count, &I_coded_crosspolar_V[pos1]);
		    RSP_Correlate (&Q_uncorr_coded_crosspolar_V[pos1], xcodes[k], param.samples_per_pulse, count, &Q_coded_crosspolar_V[pos1]);
		    /* prepare pos1 for the next journey around the for loop */
		    pos1 += param.samples_per_pulse;
		    /* we need to set k to the first code if it equals the number of the codes used */
		    if (k == 4)
		    {
			k = 0;
		    }
		    else
		    {
			/* increment the code */
			k = 4;
		    }
		}

		/* coherently average the Uncoded data */
		/* this offsets the marker to the first sample of the uncoded pulse */
		/* we assume here that the uncoded pulse follows the coded pulses */
		marker = 0;
		pos2   = 0;
		for (i = 0; i < paramUncoded.nfft; i++)
		{
		    for (j = 0; j < paramUncoded.samples_per_pulse; j++)
		    {
			/* set the index to the start of each coherent addition */
			pos1 = marker + j;
			/* look to see if we need to do any coherent averaging */
			/* if we do not then we do not have to do anything */
			if (paramUncoded.pulses_coherently_averaged > 1)
			{
			    /* this is the routine for coherently averaging uncoded data */
			    tempIco_H = 0;
			    tempQco_H = 0;
			    tempIcr_H = 0;
			    tempQcr_H = 0;
			    tempIco_V = 0;
			    tempQco_V = 0;
			    tempIcr_V = 0;
			    tempQcr_V = 0;
			    for (k = 0; k < paramUncoded.pulses_coherently_averaged; k++)
			    {
				// printf ("i: %d, j: %d, k: %d, pos1: %d, pos2: %d\n", i, j, k, pos1, pos2);
				tempIco_H += I_uncoded_copolar_H   [pos1];
				tempQco_H += Q_uncoded_copolar_H   [pos1];
				tempIcr_H += I_uncoded_crosspolar_H[pos1];
				tempQcr_H += Q_uncoded_crosspolar_H[pos1];
				tempIco_V += I_uncoded_copolar_V   [pos1];
				tempQco_V += Q_uncoded_copolar_V   [pos1];
				tempIcr_V += I_uncoded_crosspolar_V[pos1];
				tempQcr_V += Q_uncoded_crosspolar_V[pos1];
				pos1 += paramUncoded.samples_per_pulse;
			    }
			    I_uncoded_copolar_H   [pos2] = tempIco_H / paramUncoded.pulses_coherently_averaged;
			    Q_uncoded_copolar_H   [pos2] = tempQco_H / paramUncoded.pulses_coherently_averaged;
			    I_uncoded_crosspolar_H[pos2] = tempIcr_H / paramUncoded.pulses_coherently_averaged;
			    Q_uncoded_crosspolar_H[pos2] = tempQcr_H / paramUncoded.pulses_coherently_averaged;
			    I_uncoded_copolar_V   [pos2] = tempIco_V / paramUncoded.pulses_coherently_averaged;
			    Q_uncoded_copolar_V   [pos2] = tempQco_V / paramUncoded.pulses_coherently_averaged;
			    I_uncoded_crosspolar_V[pos2] = tempIcr_V / paramUncoded.pulses_coherently_averaged;
			    Q_uncoded_crosspolar_V[pos2] = tempQcr_V / paramUncoded.pulses_coherently_averaged;

			}
			pos2 ++;
		    }
		    marker += paramUncoded.samples_per_pulse * paramUncoded.pulses_coherently_averaged;
		}

		//------------------------------------
		// Add up complementary code sequence
		//------------------------------------
		tot_n_avg = paramCoded.number_of_codes * paramCoded.pulses_coherently_averaged;
		marker    = 0;
		pos2      = 0;
		for (i = 0; i < paramCoded.nfft; i++)
		{
		    for (j = 0; j < paramCoded.samples_per_pulse; j++)
		    {
			pos1      = marker + j;
			tempIco_H = 0;
			tempQco_H = 0;
			tempIcr_H = 0;
			tempQcr_H = 0;
			tempIco_V = 0;
			tempQco_V = 0;
			tempIcr_V = 0;
			tempQcr_V = 0;
			/* set the index to the start of the coherent addition */
			for (k = 0; k < tot_n_avg; k++)
			{
			    // printf ("i: %d, j: %d, k: %d, pos1: %d, pos2: %d\n", i, j, k, pos1, pos2);
			    tempIco_H += I_coded_copolar_H   [pos1];
			    tempQco_H += Q_coded_copolar_H   [pos1];
			    tempIcr_H += I_coded_crosspolar_H[pos1];
			    tempQcr_H += Q_coded_crosspolar_H[pos1];
			    tempIco_V += I_coded_copolar_V   [pos1];
			    tempQco_V += Q_coded_copolar_V   [pos1];
			    tempIcr_V += I_coded_crosspolar_V[pos1];
			    tempQcr_V += Q_coded_crosspolar_V[pos1];
			    pos1 += paramCoded.samples_per_pulse;
			}
			I_coded_copolar_H   [pos2] = tempIco_H / tot_n_avg;
			Q_coded_copolar_H   [pos2] = tempQco_H / tot_n_avg;
			I_coded_crosspolar_H[pos2] = tempIcr_H / tot_n_avg;
			Q_coded_crosspolar_H[pos2] = tempQcr_H / tot_n_avg;
			I_coded_copolar_V   [pos2] = tempIco_V / tot_n_avg;
			Q_coded_copolar_V   [pos2] = tempQco_V / tot_n_avg;
			I_coded_crosspolar_V[pos2] = tempIcr_V / tot_n_avg;
			Q_coded_crosspolar_V[pos2] = tempQcr_V / tot_n_avg;
			pos2 ++;
		    }
		    marker += (paramCoded.samples_per_pulse * paramCoded.number_of_codes * paramCoded.pulses_coherently_averaged);
		}
		/* storing IQ data  */
		if (collect_spectra_now == 1)
		{
		    /* store I and Q for each pulse */
		    /* nspectra defines the spectra number */
		    /* lets do the code first */
		    obtain_index = 0;
		    store_index  = nspectra * paramCoded.samples_per_pulse * paramCoded.nfft;
		    for (i = 0; i < paramCoded.nfft; i++)
		    {
			for (j = 0; j < paramCoded.samples_per_pulse; j++)
			{
			    IQStruct.I_coded_copolar_H   [store_index] = I_coded_copolar_H   [obtain_index];
			    IQStruct.Q_coded_copolar_H   [store_index] = Q_coded_copolar_H   [obtain_index];
			    IQStruct.I_coded_copolar_V   [store_index] = I_coded_copolar_V   [obtain_index];
			    IQStruct.Q_coded_copolar_V   [store_index] = Q_coded_copolar_V   [obtain_index];
			    IQStruct.I_coded_crosspolar_H[store_index] = I_coded_crosspolar_H[obtain_index];
			    IQStruct.Q_coded_crosspolar_H[store_index] = Q_coded_crosspolar_H[obtain_index];
			    IQStruct.I_coded_crosspolar_V[store_index] = I_coded_crosspolar_V[obtain_index];
			    IQStruct.Q_coded_crosspolar_V[store_index] = Q_coded_crosspolar_V[obtain_index];
			    store_index ++;
			    obtain_index ++;
			}
		    }
		    /* and now the uncoded */
		    obtain_index = 0;
		    store_index = nspectra * paramUncoded.samples_per_pulse * paramUncoded.nfft;
		    for (i = 0; i < paramUncoded.nfft; i++)
		    {
			for (j = 0; j < paramUncoded.samples_per_pulse; j++)
			{
			    IQStruct.I_uncoded_copolar_H   [store_index] = I_uncoded_copolar_H   [obtain_index];
			    IQStruct.Q_uncoded_copolar_H   [store_index] = Q_uncoded_copolar_H   [obtain_index];
			    IQStruct.I_uncoded_copolar_V   [store_index] = I_uncoded_copolar_V   [obtain_index];
			    IQStruct.Q_uncoded_copolar_V   [store_index] = Q_uncoded_copolar_V   [obtain_index];
			    IQStruct.I_uncoded_crosspolar_H[store_index] = I_uncoded_crosspolar_H[obtain_index];
			    IQStruct.Q_uncoded_crosspolar_H[store_index] = Q_uncoded_crosspolar_H[obtain_index];
			    IQStruct.I_uncoded_crosspolar_V[store_index] = I_uncoded_crosspolar_V[obtain_index];
			    IQStruct.Q_uncoded_crosspolar_V[store_index] = Q_uncoded_crosspolar_V[obtain_index];
			    store_index ++;
			    obtain_index ++;
			}
		    }
		}

		// Calculate power spectra for each gate
		// Loop through gates
		for (sample = 0; sample < param.samples_per_pulse; sample++)
		{
		    register int ii, index;
		    // 1) CODED H-COPOLAR SPECTRUM (HHP)
		    for (ii = 0; ii < paramCoded.nfft; ii++)
		    {
			index=ii*paramCoded.samples_per_pulse+sample;
			fftw_real_lv (in[ii]) = (float)I_coded_copolar_H[index];
			fftw_imag_lv (in[ii]) = (float)Q_coded_copolar_H[index];
		    }
		    RSP_SubtractOffset_FFTW (in, paramCoded.nfft);
		    RSP_CalcPSD_FFTW (in, paramCoded.nfft, p_coded, paramCoded.window, current_PSD, norm_coded);
		    for (ii = 0; ii < paramCoded.npsd; ii++)
			PSD[sample].HHP[ii] += current_PSD[ii] / paramCoded.spectra_averaged;

		    // 2) CODED H-CROSSPOLAR SPECTRUM (HVP)
		    for (ii = 0; ii < paramCoded.nfft; ii++)
		    {
			index=ii*paramCoded.samples_per_pulse+sample;
			fftw_real_lv (in[ii]) = (float)I_coded_crosspolar_H[index];
			fftw_imag_lv (in[ii]) = (float)Q_coded_crosspolar_H[index];
		    }
		    RSP_SubtractOffset_FFTW (in, paramCoded.nfft);
		    RSP_CalcPSD_FFTW (in, paramCoded.nfft, p_coded, paramCoded.window, current_PSD, norm_coded);
		    for (ii = 0; ii < paramCoded.npsd; ii++)
			PSD[sample].HVP[ii] += current_PSD[ii] / paramCoded.spectra_averaged;

		    // 1) CODED V-COPOLAR SPECTRUM (VVP)
		    for (ii = 0; ii < paramCoded.nfft; ii++)
		    {
			index=ii*paramCoded.samples_per_pulse+sample;
			fftw_real_lv (in[ii]) = (float)I_coded_copolar_V[index];
			fftw_imag_lv (in[ii]) = (float)Q_coded_copolar_V[index];
		    }
		    RSP_SubtractOffset_FFTW (in, paramCoded.nfft);
		    RSP_CalcPSD_FFTW (in, paramCoded.nfft, p_coded, paramCoded.window, current_PSD, norm_coded);
		    for (ii = 0; ii < paramCoded.npsd; ii++)
			PSD[sample].VVP[ii] += current_PSD[ii] / paramCoded.spectra_averaged;

		    // 2) CODED V-CROSSPOLAR SPECTRUM (VHP)
		    for (ii = 0; ii < paramCoded.nfft; ii++)
		    {
			index = ii * paramCoded.samples_per_pulse+sample;
			fftw_real_lv (in[ii]) = (float)I_coded_crosspolar_V[index];
			fftw_imag_lv (in[ii]) = (float)Q_coded_crosspolar_V[index];
		    }
		    RSP_SubtractOffset_FFTW (in, paramCoded.nfft);
		    RSP_CalcPSD_FFTW (in, paramCoded.nfft, p_coded, paramCoded.window, current_PSD, norm_coded);
		    for (ii = 0; ii < paramCoded.npsd; ii++)
			PSD[sample].VHP[ii] += current_PSD[ii] / paramCoded.spectra_averaged;


		    // 3) UNCODED H-COPOLAR SPECTRUM (HH)
		    for (ii = 0; ii < paramUncoded.nfft; ii++)
		    {
			index = ii*paramUncoded.samples_per_pulse + sample;
			//    printf ("I H = %d gate = %d \n", I_uncoded_copolar_H[index], sample);
			//    printf ("Q H = %d gate = %d \n", Q_uncoded_copolar_H[index], sample);

			fftw_real_lv (in[ii]) = (float)I_uncoded_copolar_H[index];
			fftw_imag_lv (in[ii]) = (float)Q_uncoded_copolar_H[index];
		    }
		    RSP_SubtractOffset_FFTW (in, paramUncoded.nfft);
		    RSP_CalcPSD_FFTW (in, paramUncoded.nfft, p_uncoded, paramUncoded.window, current_PSD, norm_uncoded);
		    for (ii = 0; ii < paramUncoded.npsd; ii++)
			PSD[sample].HH[ii] += current_PSD[ii] / paramUncoded.spectra_averaged;

		    // printf ("value = %d gate = %d bin = %d\n", PSD[sample].HH[ii], sample, ii);
		    // for (ii = 0; ii < paramUncoded.npsd; ii++)
		    // {
		    //   printf ("value = %f gate = %i bin = %d\n", current_PSD[ii], sample, ii);
		    // }

		    // 4) UNCODED H-CROSSPOLAR SPECTRUM (HV)
		    for (ii = 0; ii < paramUncoded.nfft; ii++)
		    {
			index = ii * paramUncoded.samples_per_pulse + sample;
			fftw_real_lv (in[ii]) = (float)I_uncoded_crosspolar_H[index];
			fftw_imag_lv (in[ii]) = (float)Q_uncoded_crosspolar_H[index];
		    }
		    RSP_SubtractOffset_FFTW (in, paramUncoded.nfft);
		    RSP_CalcPSD_FFTW (in, paramUncoded.nfft, p_uncoded, paramUncoded.window, current_PSD, norm_uncoded);
		    for (ii = 0; ii < paramUncoded.npsd; ii++)
			PSD[sample].HV[ii] += current_PSD[ii] / paramUncoded.spectra_averaged;

		    // 3) UNCODED V-COPOLAR SPECTRUM (VV)
		    for (ii = 0; ii < paramUncoded.nfft; ii++)
		    {
			index = ii * paramUncoded.samples_per_pulse + sample;
			fftw_real_lv (in[ii]) = (float)I_uncoded_copolar_V[index];
			fftw_imag_lv (in[ii]) = (float)Q_uncoded_copolar_V[index];
		    }
		    RSP_SubtractOffset_FFTW (in, paramUncoded.nfft);
		    RSP_CalcPSD_FFTW (in, paramUncoded.nfft, p_uncoded, paramUncoded.window, current_PSD, norm_uncoded);
		    for (ii = 0; ii < paramUncoded.npsd; ii++)
			PSD[sample].VV[ii] += current_PSD[ii] / paramUncoded.spectra_averaged;

		    // 4) UNCODED V-CROSSPOLAR SPECTRUM (VH)
		    for (ii = 0; ii < paramUncoded.nfft; ii++)
		    {
			index = ii * paramUncoded.samples_per_pulse + sample;
			//    printf ("I V = %d gate = %d \n", I_uncoded_crosspolar_V[index], sample);
			//    printf ("Q V = %d gate = %d \n", Q_uncoded_crosspolar_V[index], sample);
			fftw_real_lv (in[ii]) = (float)I_uncoded_crosspolar_V[index];
			fftw_imag_lv (in[ii]) = (float)Q_uncoded_crosspolar_V[index];
		    }
		    RSP_SubtractOffset_FFTW (in, paramUncoded.nfft);
		    RSP_CalcPSD_FFTW (in, paramUncoded.nfft, p_uncoded, paramUncoded.window, current_PSD, norm_uncoded);
		    for (ii = 0; ii < paramUncoded.npsd; ii++)
			PSD[sample].VH[ii] += current_PSD[ii] / paramUncoded.spectra_averaged;

		    // for (ii = 0; ii < paramUncoded.npsd; ii++)
		    // {
		    //   printf ("value = %f gate = %i bin = %d\n", current_PSD[ii], sample, ii);
		    // }

		    // 5) UNCODED CROSS-POL VARIABLES
		    for (ii = 0; ii < paramUncoded.nfft; ii++)
		    {
			index = ii * paramUncoded.samples_per_pulse + sample;
			fftw_real_lv (in [ii]) = (float)I_uncoded_copolar_H[index];
			fftw_imag_lv (in [ii]) = (float)Q_uncoded_copolar_H[index];
			fftw_real_lv (inx[ii]) = (float)I_uncoded_copolar_V[index];
			fftw_imag_lv (inx[ii]) = (float)Q_uncoded_copolar_V[index];
		    }

		    RSP_SubtractOffset_FFTW (in,  paramUncoded.nfft);
		    RSP_SubtractOffset_FFTW (inx, paramUncoded.nfft);

		    tempI_HV   = 0.0;
		    tempQ_HV   = 0.0;
		    mod_HH     = 0.0;
		    mod_VV     = 0.0;
		    tot_pulses = paramUncoded.nfft * paramUncoded.spectra_averaged * param.moments_averaged;
		    for (ii = 0; ii < paramUncoded.nfft; ii++)
		    {
			I_UNCH[sample] += fftw_real (in [ii]) / tot_pulses;
			Q_UNCH[sample] += fftw_imag (in [ii]) / tot_pulses;
			I_UNCV[sample] += fftw_real (inx[ii]) / tot_pulses;
			Q_UNCV[sample] += fftw_imag (inx[ii]) / tot_pulses;
			tempI_HV += (fftw_real (inx[ii]) * fftw_real (in [ii])) + (fftw_imag (in [ii]) * fftw_imag (inx[ii]));
			tempQ_HV += (fftw_real (in [ii]) * fftw_imag (inx[ii])) - (fftw_real (inx[ii]) * fftw_imag (in [ii]));
			mod_HH   += (fftw_real (in [ii]) * fftw_real (in [ii])) + (fftw_imag (in [ii]) * fftw_imag (in [ii]));
			mod_VV   += (fftw_real (inx[ii]) * fftw_real (inx[ii])) + (fftw_imag (inx[ii]) * fftw_imag (inx[ii]));

			if (ii > 0)
			{
			    if (horizontal_first == 1)
			    {
				tempI_HV += (fftw_real (inx[ii-1]) * fftw_real (in [ii  ])) + (fftw_imag (in [ii  ]) * fftw_imag (inx[ii-1]));
				tempQ_HV += (fftw_real (in [ii  ]) * fftw_imag (inx[ii-1])) - (fftw_real (inx[ii-1]) * fftw_imag (in [ii  ]));
				mod_HH   += (fftw_real (in [ii  ]) * fftw_real (in [ii  ])) + (fftw_imag (in [ii  ]) * fftw_imag (in [ii  ]));
				mod_VV   += (fftw_real (inx[ii-1]) * fftw_real (inx[ii-1])) + (fftw_imag (inx[ii-1]) * fftw_imag (inx[ii-1]));
			    }
			    else
			    {
				tempI_HV += (fftw_real (inx[ii  ]) * fftw_real (in [ii-1])) + (fftw_imag (in [ii-1]) * fftw_imag (inx[ii  ]));
				tempQ_HV += (fftw_real (in [ii-1]) * fftw_imag (inx[ii  ])) - (fftw_real (inx[ii  ]) * fftw_imag (in [ii-1]));
				mod_HH   += (fftw_real (in [ii-1]) * fftw_real (in [ii-1])) + (fftw_imag (in [ii-1]) * fftw_imag (in [ii-1]));
				mod_VV   += (fftw_real (inx[ii  ]) * fftw_real (inx[ii  ])) + (fftw_imag (inx[ii  ]) * fftw_imag (inx[ii  ]));
			    }
			}
		    }

		    real_PHIDP_C[sample] += tempI_HV;
		    imag_PHIDP_C[sample] += tempQ_HV;
		    tempr = sqrt (((tempI_HV * tempI_HV)+ (tempQ_HV * tempQ_HV))/ (mod_HH * mod_VV));
		    tempz = atanh (tempr);
		    RHOHV_C[sample] += tempz / paramUncoded.spectra_averaged / param.moments_averaged;

		    for (ii = 0; ii < paramUncoded.nfft; ii++)
		    {
			h[ii] = fftw_real (in [ii]) * fftw_real (in [ii]) + fftw_imag (in [ii]) * fftw_imag (in [ii]);
			v[ii] = fftw_real (inx[ii]) * fftw_real (inx[ii]) + fftw_imag (inx[ii]) * fftw_imag (inx[ii]);
		    }

		    tempr = corrCoeff (h, v, paramUncoded.nfft);
		    tempz = atanh (tempr);
		    RHOHV_P[sample] += tempz / paramUncoded.spectra_averaged / param.moments_averaged;

		    double_interp (h, hh, paramUncoded.nfft);
		    double_interp (v, vv, paramUncoded.nfft);

		    for (ii = 0; ii < paramUncoded.nfft * 2 - 1; ii++)
		    {
			if (horizontal_first)
			{
			    xx[ii] = hh[ii];
			    yy[ii] = vv[ii];
			}
			else
			{
			    xx[ii] = vv[ii];
			    yy[ii] = hh[ii];
			}
		    }

		    tempr = corrCoeffPoly3 (xx, yy, paramUncoded.nfft);
		    tempz = atanh (tempr);
		    RHOHV_IP[sample] += tempz / paramUncoded.spectra_averaged / param.moments_averaged;

		    // 5) CODED CROSS-POL VARIABLES
		    for (ii = 0; ii < paramCoded.nfft; ii++)
		    {
			index = ii * paramCoded.samples_per_pulse + sample;
			fftw_real_lv (in [ii]) = (float)I_coded_copolar_H[index];
			fftw_imag_lv (in [ii]) = (float)Q_coded_copolar_H[index];
			fftw_real_lv (inx[ii]) = (float)I_coded_copolar_V[index];
			fftw_imag_lv (inx[ii]) = (float)Q_coded_copolar_V[index];
		    }

		    RSP_SubtractOffset_FFTW (in,  paramCoded.nfft);
		    RSP_SubtractOffset_FFTW (inx, paramCoded.nfft);

		    tempI_HV   = 0.0;
		    tempQ_HV   = 0.0;
		    mod_HH     = 0.0;
		    mod_VV     = 0.0;
		    tot_pulses = paramCoded.nfft * paramCoded.spectra_averaged * param.moments_averaged;
		    for (ii = 0; ii < paramCoded.nfft; ii++)
		    {
			I_CODH[sample] += fftw_real (in [ii]) / tot_pulses;
			Q_CODH[sample] += fftw_imag (in [ii]) / tot_pulses;
			I_CODV[sample] += fftw_real (inx[ii]) / tot_pulses;
			Q_CODV[sample] += fftw_imag (inx[ii]) / tot_pulses;
			tempI_HV += (fftw_real (inx[ii]) * fftw_real (in [ii])) + (fftw_imag (in [ii]) * fftw_imag (inx[ii]));
			tempQ_HV += (fftw_real (inx[ii]) * fftw_imag (in [ii])) - (fftw_real (in [ii]) * fftw_imag (inx[ii]));
			mod_HH   += (fftw_real (in [ii]) * fftw_real (in [ii])) + (fftw_imag (in [ii]) * fftw_imag (in [ii]));
			mod_VV   += (fftw_real (inx[ii]) * fftw_real (inx[ii])) + (fftw_imag (inx[ii]) * fftw_imag (inx[ii]));

			// if (ii > 0)
			// {
			//   if (horizontal_first == 1)
			//   {
			//     tempI_HV += (inx[ii-1][0] * in [ii  ][0]) + (in [ii  ][1] * inx[ii-1][1]);
			//     tempQ_HV +=-(in [ii  ][0] * inx[ii-1][1]) + (inx[ii-1][0] * in [ii  ][1]);
			//     mod_HH   += (in [ii  ][0] * in [ii  ][0]) + (in [ii  ][1] * in [ii  ][1]);
			//     mod_VV   += (inx[ii-1][0] * inx[ii-1][0]) + (inx[ii-1][1] * inx[ii-1][1]);
			//   }
			//   else
			//   {
			//     tempI_HV += (inx[ii  ][0] * in [ii-1][0]) + (in [ii-1][1] * inx[ii  ][1]);
			//     tempQ_HV +=-(in [ii-1][0] * inx[ii  ][1]) + (inx[ii  ][0] * in [ii-1][1]);
			//     mod_HH   += (in [ii-1][0] * in [ii-1][0]) + (in [ii-1][1] * in [ii-1][1]);
			//     mod_VV   += (inx[ii  ][0] * inx[ii  ][0]) + (inx[ii  ][1] * inx[ii  ][1]);
			//   }
			// }
		    }

		    real_PHIDP_CP[sample] += tempI_HV;
		    imag_PHIDP_CP[sample] += tempQ_HV;
		    tempr = sqrt (((tempI_HV * tempI_HV) + (tempQ_HV * tempQ_HV)) / (mod_HH * mod_VV));
		    tempz = atanh (tempr);
		    RHOHV_CP[sample] += tempz / paramCoded.spectra_averaged / param.moments_averaged;

		    for (ii = 0; ii < paramCoded.nfft; ii++)
		    {
			hp[ii] = fftw_real (in [ii]) * fftw_real (in [ii]) + fftw_imag (in [ii]) * fftw_imag (in [ii]);
			vp[ii] = fftw_real (inx[ii]) * fftw_real (inx[ii]) + fftw_imag (inx[ii]) * fftw_imag (inx[ii]);
		    }

		    tempr = corrCoeff (hp, vp, paramCoded.nfft);
		    tempz = atanh (tempr);
		    RHOHV_PP[sample] += tempz / paramCoded.spectra_averaged / param.moments_averaged;

		    double_interp (hp, hhp, paramCoded.nfft);
		    double_interp (vp, vvp, paramCoded.nfft);

		    for (ii = 0; ii < paramCoded.nfft * 2 - 1; ii++)
		    {
			if (horizontal_first)
			{
			    xxp[ii] = hhp[ii];
			    yyp[ii] = vvp[ii];
			}
			else
			{
			    xxp[ii] = vvp[ii];
			    yyp[ii] = hhp[ii];
			}
		    }

		    tempr = corrCoeffPoly3 (xxp, yyp, paramCoded.nfft);
		    tempz = atanh (tempr);
		    RHOHV_IPP[sample] += tempz / paramCoded.spectra_averaged / param.moments_averaged;
		}
	    }
	    //---------------------------
	    // END OF SPECTRAL AVERAGING
	    //---------------------------

	    /* update time in spectral information file */
	    /* get time of day */
	    gettimeofday (&tv, &tz);
	    gmtime_r (&tv.tv_sec, &tm);
	    PSD_obs.year        = tm.tm_year + 1900;
	    PSD_obs.month       = tm.tm_mon  + 1;
	    PSD_obs.day         = tm.tm_mday;
	    PSD_obs.hour        = tm.tm_hour;
	    PSD_obs.minute      = tm.tm_min;
	    PSD_obs.second      = tm.tm_sec;
	    PSD_obs.centisecond = obs.centisecond = tv.tv_usec / 10000U;

	    PSD_RAPID_obs.year        = PSD_obs.year;
	    PSD_RAPID_obs.month       = PSD_obs.month;
	    PSD_RAPID_obs.day         = PSD_obs.day;
	    PSD_RAPID_obs.hour        = PSD_obs.hour;
	    PSD_RAPID_obs.minute      = PSD_obs.minute;
	    PSD_RAPID_obs.second      = PSD_obs.second;
	    PSD_RAPID_obs.centisecond = PSD_obs.centisecond;

	    system_time = time (NULL);
	    if (collect_spectra_rapid_now == 1)
	    {
		printf ("Writing Rapid PSD Variables... ***************************\n");
		RNC_WriteRapidLogPSDVariables (spectra_rapid_ncid, COPERNICUS_SPECTRA_RAPID, &param, &PSD_RAPID_obs, PSD, PSD_rapid_varid);
		status = nc_sync (spectra_ncid);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
		spectra_rapid_time = spectra_rapid_time + param.dump_spectra_rapid;
		printf ("Written Rapid PSD Variables\n");
	    }
	    /* write out spectra to netCDF if required */
	    if (collect_spectra_now == 1)
	    {
		/* end of section that was not thought out */
		printf ("Writing LogPSD Variables... ***************************\n");
		RNC_WriteLogPSDVariables (spectra_ncid, COPERNICUS_CODED_SPECTRA, &param, &PSD_obs, PSD, &IQStruct, PSD_varid);
		status = nc_sync (spectra_ncid);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
		spectra_time = spectra_time + param.dump_spectra;
	    }
	    // Calculate noise from upper range gates
	    noisegate1 = param.samples_per_pulse - 50 - (paramCoded.code_length * paramCoded.oversample_ratio);
	    noisegate2 = param.samples_per_pulse -  1 - (paramCoded.code_length * paramCoded.oversample_ratio);
	    count      = (noisegate2 - noisegate1) + 1;
	    HH_noise_level  = 0;
	    HV_noise_level  = 0;
	    VV_noise_level  = 0;
	    VH_noise_level  = 0;
	    HHP_noise_level = 0;
	    HVP_noise_level = 0;
	    VVP_noise_level = 0;
	    VHP_noise_level = 0;
	    for (i = noisegate1; i <= noisegate2; i++)
	    {
		HH_noise_level  += median (PSD[i].HH,  paramUncoded.npsd);
		HV_noise_level  += median (PSD[i].HV,  paramUncoded.npsd);
		VV_noise_level  += median (PSD[i].VV,  paramUncoded.npsd);
		VH_noise_level  += median (PSD[i].VH,  paramUncoded.npsd);
		HHP_noise_level += median (PSD[i].HHP, paramCoded.npsd);
		HVP_noise_level += median (PSD[i].HVP, paramCoded.npsd);
		VVP_noise_level += median (PSD[i].VVP, paramCoded.npsd);
		VHP_noise_level += median (PSD[i].VHP, paramCoded.npsd);
	    }
	    HH_noise_level  /= count;
	    HV_noise_level  /= count;
	    VV_noise_level  /= count;
	    VH_noise_level  /= count;
	    HHP_noise_level /= count;
	    HVP_noise_level /= count;
	    VVP_noise_level /= count;
	    VHP_noise_level /= count;

	    NPC_H [0] += HH_noise_level;
	    NPC_HP[0] += HHP_noise_level;
	    NPC_V [0] += HV_noise_level;
	    NPC_VP[0] += HVP_noise_level;

	    // Now calculate the Doppler parameters
	    printf ("** Calculating Doppler parameters...\n");
	    // Loop through all spectra and get parameters
	    for (i = 0; i < param.samples_per_pulse; i++)
	    {
		float noise_power, tempPower, tempVel, tempZED;

		// interpolate over clutter
		RSP_ClutterInterp (PSD[i].HH, paramUncoded.npsd, paramUncoded.fft_bins_interpolated);
		RSP_ClutterInterp (PSD[i].HV, paramUncoded.npsd, paramUncoded.fft_bins_interpolated);
		RSP_ClutterInterp (PSD[i].VV, paramUncoded.npsd, paramUncoded.fft_bins_interpolated);
		RSP_ClutterInterp (PSD[i].VH, paramUncoded.npsd, paramUncoded.fft_bins_interpolated);
		RSP_ClutterInterp (PSD[i].HHP, paramCoded.npsd, paramCoded.fft_bins_interpolated);
		RSP_ClutterInterp (PSD[i].HVP, paramCoded.npsd, paramCoded.fft_bins_interpolated);
		RSP_ClutterInterp (PSD[i].VVP, paramCoded.npsd, paramCoded.fft_bins_interpolated);
		RSP_ClutterInterp (PSD[i].VHP, paramCoded.npsd, paramCoded.fft_bins_interpolated);


		// Find HH peak
		RSP_FindPeaksMulti_Destructive (PSD[i].HH, paramUncoded.npsd, paramUncoded.num_peaks, HH_noise_level, HH_peaks);
		// Calculate HH moments
		RSP_CalcSpecMom (PSD[i].HH, paramUncoded.npsd, HH_peaks, HH_noise_level, HH_moments, RSP_MOMENTS);
		// Find HV peak
		RSP_FindPeaksMulti_Destructive (PSD[i].HV, paramUncoded.npsd, paramUncoded.num_peaks, HV_noise_level, HV_peaks);
		// Calculate HV moments
		RSP_CalcSpecMom (PSD[i].HV, paramUncoded.npsd, HV_peaks, HV_noise_level, HV_moments, RSP_MOMENTS);
		// Find VV peak
		RSP_FindPeaksMulti_Destructive (PSD[i].VV, paramUncoded.npsd, paramUncoded.num_peaks, VV_noise_level, VV_peaks);
		// Calculate VV moments
		RSP_CalcSpecMom (PSD[i].VV, paramUncoded.npsd, VV_peaks, VV_noise_level, VV_moments, RSP_MOMENTS);
		// Find VH peak
		RSP_FindPeaksMulti_Destructive (PSD[i].VH, paramUncoded.npsd, paramUncoded.num_peaks, VH_noise_level, VH_peaks);
		// Calculate HV moments
		RSP_CalcSpecMom (PSD[i].VH, paramUncoded.npsd, VH_peaks, VH_noise_level, VH_moments, RSP_MOMENTS);


		// Find HHP peak
		RSP_FindPeaksMulti_Destructive (PSD[i].HHP, paramCoded.npsd, paramCoded.num_peaks, HHP_noise_level, HHP_peaks);
		// Calculate HHP moments
		RSP_CalcSpecMom (PSD[i].HHP, paramCoded.npsd, HHP_peaks, HHP_noise_level, HHP_moments, RSP_MOMENTS);
		// Find HVP peak
		RSP_FindPeaksMulti_Destructive (PSD[i].HVP, paramCoded.npsd, paramCoded.num_peaks, HVP_noise_level, HVP_peaks);
		// Calculate HVP moments
		RSP_CalcSpecMom (PSD[i].HVP, paramCoded.npsd, HVP_peaks, HVP_noise_level, HVP_moments, RSP_MOMENTS);
		// Find VVP peak
		RSP_FindPeaksMulti_Destructive (PSD[i].VVP, paramCoded.npsd, paramCoded.num_peaks, VVP_noise_level, VVP_peaks);
		// Calculate VVP moments
		RSP_CalcSpecMom (PSD[i].VVP, paramCoded.npsd, VVP_peaks, VVP_noise_level, VVP_moments, RSP_MOMENTS);
		// Find VHP peak
		RSP_FindPeaksMulti_Destructive (PSD[i].VHP, paramCoded.npsd, paramCoded.num_peaks, VHP_noise_level, VHP_peaks);
		// Calculate VHP moments
		RSP_CalcSpecMom (PSD[i].VHP, paramCoded.npsd, VHP_peaks, VHP_noise_level, VHP_moments, RSP_MOMENTS);


		// ----------------------------
		//  PROCESS UNCODED PARAMETERS
		// ----------------------------

		// Calculate weighting coefficient
		//wi = peaks[0].peakPSD - HH_noise_level;
		//wi = tempPower / noise_power;
		//if (wi > 1.0f) wi = 1.0f;
		wi = 1.0f; // Turn off weighting
		uncoded_sum_wi[i] += wi;

		//****************************************
		// HH and HV or should I say _HC and _XHC
		//****************************************

		gate= i + gate_offset;
		noise_power = RSP_CalcNoisePower (HH_noise_level, HH_peaks, &paramUncoded);
		tempPower   = HH_moments[0] * paramUncoded.frequency_bin_width;
		// printf ("tempPower = %d %5.2f %f\n", i, tempPower, paramUncoded.frequency_bin_width);
		// printf ("noisePower = %d %5.2f\n", i, noise_power);

		// COPOLAR
		tempZED         = 10.0 * log10 (tempPower);
		tempVel         = RSP_BinToVelocity (HH_moments[1], &paramUncoded);
		SNR_HC     [i] += tempPower / noise_power * wi;
		ZED_HC     [i] += tempPower               * wi;
		VEL_HC_COS [i] += cos (tempVel / paramUncoded.folding_velocity * PI) * wi;
		VEL_HC_SIN [i] += sin (tempVel / paramUncoded.folding_velocity * PI) * wi;
		SPW_HC     [i] += HH_moments[2] * paramUncoded.frequency_bin_width / paramUncoded.hz_per_mps * wi;
		mean_vsq_HC[i] += tempVel   * tempVel   * wi;
		mean_Zsq_HC[i] += tempPower * tempPower * wi;

		// CROSSPOLAR
		noise_power = RSP_CalcNoisePower (HV_noise_level, HV_peaks, &paramUncoded);
		tempPower   = HV_moments[0] * paramUncoded.frequency_bin_width;
		SNR_XHC[i] += tempPower / noise_power * wi;
		ZED_XHC[i] += tempPower               * wi;

		//****************************************
		// VV and VH or should I say _VC and _XVC
		//****************************************

		gate        = i + gate_offset;
		noise_power = RSP_CalcNoisePower (VV_noise_level, VV_peaks, &paramUncoded);
		tempPower   = VV_moments[0] * paramUncoded.frequency_bin_width;

		// COPOLAR
		tempZED         = 10.0 * log10 (tempPower);
		tempVel         = RSP_BinToVelocity (VV_moments[1], &paramUncoded);
		SNR_VC     [i] += tempPower / noise_power * wi;
		ZED_VC     [i] += tempPower               * wi;
		VEL_VC_COS [i] += cos (tempVel / paramUncoded.folding_velocity * PI) * wi;
		VEL_VC_SIN [i] += sin (tempVel / paramUncoded.folding_velocity * PI) * wi;
		SPW_VC     [i] += VV_moments[2] * paramUncoded.frequency_bin_width / paramUncoded.hz_per_mps * wi;
		mean_vsq_VC[i] += tempVel   * tempVel   * wi;
		mean_Zsq_VC[i] += tempPower * tempPower * wi;

		// CROSSPOLAR
		noise_power = RSP_CalcNoisePower (VH_noise_level, VH_peaks, &paramUncoded);
		tempPower   = VH_moments[0] * paramUncoded.frequency_bin_width;
		SNR_XVC[i] += tempPower / noise_power * wi;
		ZED_XVC[i] += tempPower               * wi;


		// ----------------------------
		//  PROCESS CODED PARAMETERS
		// ----------------------------
		if (gate < paramCoded.samples_per_pulse)
		{
		    // Calculate weighting coefficient
		    //wi = peaks[0].peakPSD - HHP_noise_level;
		    //wi = tempPower / noise_power;
		    //if (wi > 1.0f) wi = 1.0f;
		    wi = 1.0f; // Turn off weighting
		    coded_sum_wi[gate] += wi;

		    //********************************************
		    // HHP and HVP or should I say _HCP and _XHCP
		    //********************************************

		    noise_power = RSP_CalcNoisePower (HHP_noise_level, HHP_peaks, &paramCoded);
		    tempPower   = HHP_moments[0] * paramCoded.frequency_bin_width;

		    // COPOLAR
		    tempVel             = RSP_BinToVelocity (HHP_moments[1], &paramCoded);
		    tempZED             = 10.0 * log10 (tempPower);
		    SNR_HCP     [gate] += tempPower / noise_power * wi;
		    ZED_HCP     [gate] += tempPower               * wi;
		    VEL_HCP_COS [gate] += cos (tempVel / paramCoded.folding_velocity * PI) * wi;
		    VEL_HCP_SIN [gate] += sin (tempVel / paramCoded.folding_velocity * PI) * wi;
		    SPW_HCP     [gate] += HHP_moments[2] * paramCoded.frequency_bin_width / paramCoded.hz_per_mps * wi;
		    mean_vsq_HCP[gate] += tempVel   * tempVel   * wi;
		    mean_Zsq_HCP[gate] += tempPower * tempPower * wi;

		    // CROSSPOLAR
		    noise_power     = RSP_CalcNoisePower (HVP_noise_level, HVP_peaks, &paramCoded);
		    tempPower       = HVP_moments[0] * paramCoded.frequency_bin_width;
		    SNR_XHCP[gate] += tempPower / noise_power * wi;
		    ZED_XHCP[gate] += tempPower               * wi;

		    //********************************************
		    // VVP and VHP or should I say _VCP and _XVCP
		    //********************************************

		    noise_power = RSP_CalcNoisePower (VVP_noise_level, VVP_peaks, &paramCoded);
		    tempPower   = VVP_moments[0]*paramCoded.frequency_bin_width;

		    // COPOLAR
		    tempVel             = RSP_BinToVelocity (VVP_moments[1], &paramCoded);
		    tempZED             = 10.0 * log10 (tempPower);
		    SNR_VCP     [gate] += tempPower / noise_power * wi;
		    ZED_VCP     [gate] += tempPower               * wi;
		    VEL_VCP_COS [gate] += cos (tempVel / paramCoded.folding_velocity * PI) * wi;
		    VEL_VCP_SIN [gate] += sin (tempVel / paramCoded.folding_velocity * PI) * wi;
		    SPW_VCP     [gate] += VVP_moments[2] * paramCoded.frequency_bin_width / paramCoded.hz_per_mps * wi;
		    mean_vsq_VCP[gate] += tempVel   * tempVel   * wi;
		    mean_Zsq_VCP[gate] += tempPower * tempPower * wi;

		    // CROSSPOLAR
		    noise_power     = RSP_CalcNoisePower (VHP_noise_level, VHP_peaks, &paramCoded);
		    tempPower       = VHP_moments[0] * paramCoded.frequency_bin_width;
		    SNR_XVCP[gate] += tempPower / noise_power * wi;
		    ZED_XVCP[gate] += tempPower               * wi;
		}
	    }

	} // End of moments averaging loop

	for (i = 0; i < param.samples_per_pulse; i++)
	{
	    // COMPLETE THE WEIGHTED AVERAGING WITH DIVISION
	    // Horizontal
	    SNR_HC     [i] /= uncoded_sum_wi[i];
	    ZED_HC     [i] /= uncoded_sum_wi[i];
	    VEL_HC_COS [i] /= uncoded_sum_wi[i];
	    VEL_HC_SIN [i] /= uncoded_sum_wi[i];
	    SPW_HC     [i] /= uncoded_sum_wi[i];
	    SNR_XHC    [i] /= uncoded_sum_wi[i];
	    ZED_XHC    [i] /= uncoded_sum_wi[i];
	    SNR_HCP    [i] /= coded_sum_wi[i];
	    ZED_HCP    [i] /= coded_sum_wi[i];
	    VEL_HCP_COS[i] /= coded_sum_wi[i];
	    VEL_HCP_SIN[i] /= coded_sum_wi[i];
	    SPW_HCP    [i] /= coded_sum_wi[i];
	    SNR_XHCP   [i] /= coded_sum_wi[i];
	    ZED_XHCP   [i] /= coded_sum_wi[i];
	    // Vertical
	    SNR_VC     [i] /= uncoded_sum_wi[i];
	    ZED_VC     [i] /= uncoded_sum_wi[i];
	    VEL_VC_COS [i] /= uncoded_sum_wi[i];
	    VEL_VC_SIN [i] /= uncoded_sum_wi[i];
	    SPW_VC     [i] /= uncoded_sum_wi[i];
	    SNR_XVC    [i] /= uncoded_sum_wi[i];
	    ZED_XVC    [i] /= uncoded_sum_wi[i];
	    SNR_VCP    [i] /= coded_sum_wi[i];
	    ZED_VCP    [i] /= coded_sum_wi[i];
	    VEL_VCP_COS[i] /= coded_sum_wi[i];
	    VEL_VCP_SIN[i] /= coded_sum_wi[i];
	    SPW_VCP    [i] /= coded_sum_wi[i];
	    SNR_XVCP   [i] /= coded_sum_wi[i];
	    ZED_XVCP   [i] /= coded_sum_wi[i];

	    // Horizontal
	    VEL_HC [i] = atan2 (VEL_HC_SIN [i], VEL_HC_COS [i]) / PI * paramUncoded.folding_velocity;
	    VEL_HCP[i] = atan2 (VEL_HCP_SIN[i], VEL_HCP_COS[i]) / PI * paramCoded.folding_velocity;
	    // Vertical
	    VEL_VC [i] = atan2 (VEL_VC_SIN [i], VEL_VC_COS [i]) / PI * paramUncoded.folding_velocity;
	    VEL_VCP[i] = atan2 (VEL_VCP_SIN[i], VEL_VCP_COS[i]) / PI * paramCoded.folding_velocity;

	    // Calculate sigma Z bar (this is the relative linear standard deviation)
	    // Horizontal
	    mean_Zsq_HC [i] /= uncoded_sum_wi[i];
	    mean_Zsq_HCP[i] /= coded_sum_wi[i];
	    ZED_HCD [i] = sqrt (mean_Zsq_HC [i] - (ZED_HC [i] * ZED_HC [i]));
	    ZED_HCDP[i] = sqrt (mean_Zsq_HCP[i] - (ZED_HCP[i] * ZED_HCP[i]));
	    ZED_HCD [i] = 10.0 * log10 (ZED_HCD [i] / ZED_HC [i]);
	    ZED_HCDP[i] = 10.0 * log10 (ZED_HCDP[i] / ZED_HCP[i]);
	    // Vertical
	    mean_Zsq_VC [i] /= uncoded_sum_wi[i];
	    mean_Zsq_VCP[i] /= coded_sum_wi[i];
	    ZED_VCD [i] = sqrt (mean_Zsq_VC [i] - (ZED_VC [i] * ZED_VC [i]));
	    ZED_VCDP[i] = sqrt (mean_Zsq_VCP[i] - (ZED_VCP[i] * ZED_VCP[i]));
	    ZED_VCD [i] = 10.0 * log10 (ZED_VCD [i] / ZED_VC [i]);
	    ZED_VCDP[i] = 10.0 * log10 (ZED_VCDP[i] / ZED_VCP[i]);

	    // Convert SNRs and ZEDs to dB
	    // Horizontal
	    SNR_HC  [i] = 10.0 * log10 (SNR_HC[i]);
	    ZED_HC  [i] = 10.0 * log10 (ZED_HC[i]);
	    SNR_HCP [i] = 10.0 * log10 (SNR_HCP[i]);
	    ZED_HCP [i] = 10.0 * log10 (ZED_HCP[i]);
	    SNR_XHC [i] = 10.0 * log10 (SNR_XHC[i]);
	    ZED_XHC [i] = 10.0 * log10 (ZED_XHC[i]);
	    SNR_XHCP[i] = 10.0 * log10 (SNR_XHCP[i]);
	    ZED_XHCP[i] = 10.0 * log10 (ZED_XHCP[i]);
	    // Vertical
	    SNR_VC  [i] = 10.0 * log10 (SNR_VC[i]);
	    ZED_VC  [i] = 10.0 * log10 (ZED_VC[i]);
	    SNR_VCP [i] = 10.0 * log10 (SNR_VCP[i]);
	    ZED_VCP [i] = 10.0 * log10 (ZED_VCP[i]);
	    SNR_XVC [i] = 10.0 * log10 (SNR_XVC[i]);
	    ZED_XVC [i] = 10.0 * log10 (ZED_XVC[i]);
	    SNR_XVCP[i] = 10.0 * log10 (SNR_XVCP[i]);
	    ZED_XVCP[i] = 10.0 * log10 (ZED_XVCP[i]);

	    // Calculate sigma v bar
	    // Horizontal polarisation
	    mean_vsq_HC [i] /= uncoded_sum_wi[i];
	    mean_vsq_HCP[i] /= coded_sum_wi[i];
	    VEL_HCD [i] = sqrt (mean_vsq_HC [i] - (VEL_HC [i] * VEL_HC [i]));
	    VEL_HCDP[i] = sqrt (mean_vsq_HCP[i] - (VEL_HCP[i] * VEL_HCP[i]));
	    // Vertical polarisation
	    mean_vsq_VC [i] /= uncoded_sum_wi[i];
	    mean_vsq_VCP[i] /= coded_sum_wi[i];
	    VEL_VCD [i] = sqrt (mean_vsq_VC [i] - (VEL_VC [i] * VEL_VC [i]));
	    VEL_VCDP[i] = sqrt (mean_vsq_VCP[i] - (VEL_VCP[i] * VEL_VCP[i]));

	    // Calculate LDR
	    LDR_HC [i] += ZED_XHC [i] - ZED_HC [i] + paramUncoded.LDR_HC_calibration_offset;
	    LDR_HCP[i] += ZED_XHCP[i] - ZED_HCP[i] + paramCoded.LDR_HCP_calibration_offset;
	    LDR_VC [i] += ZED_XVC [i] - ZED_VC [i] + paramUncoded.LDR_VC_calibration_offset;
	    LDR_VCP[i] += ZED_XVCP[i] - ZED_VCP[i] + paramCoded.LDR_VCP_calibration_offset;

	    // Calculate ZDR
	    ZDR_C [i] += ZED_HC [i] - ZED_VC [i] + paramUncoded.ZDR_C_calibration_offset;
	    ZDR_CP[i] += ZED_HCP[i] - ZED_VCP[i] + paramCoded.ZDR_CP_calibration_offset;

	    // Calculate PHIDP
	    PHIDP_C [i] = atan2 (imag_PHIDP_C [i], real_PHIDP_C [i]) * 180.0 / PI;
	    PHIDP_CP[i] = atan2 (imag_PHIDP_CP[i], real_PHIDP_CP[i]) * 180.0 / PI;

	    // Linearize correlation coefficients
	    RHOHV_C  [i] = tanh (RHOHV_C  [i]);
	    RHOHV_CP [i] = tanh (RHOHV_CP [i]);
	    RHOHV_P  [i] = tanh (RHOHV_P  [i]);
	    RHOHV_PP [i] = tanh (RHOHV_PP [i]);
	    RHOHV_IP [i] = tanh (RHOHV_IP [i]);
	    RHOHV_IPP[i] = tanh (RHOHV_IPP[i]);

	    // Do range correction
	    // Horizontal
	    ZED_HC  [i] += 10.0 * log10 (paramUncoded.range[i] * paramUncoded.range[i]) + paramUncoded.ZED_HC_calibration_offset;
	    ZED_XHC [i] += 10.0 * log10 (paramUncoded.range[i] * paramUncoded.range[i]) + paramUncoded.ZED_XHC_calibration_offset;
	    ZED_HCP [i] += 10.0 * log10 (paramCoded.range  [i] * paramCoded.range  [i]) + paramCoded.ZED_HCP_calibration_offset;
	    ZED_XHCP[i] += 10.0 * log10 (paramCoded.range  [i] * paramCoded.range  [i]) + paramCoded.ZED_XHCP_calibration_offset;
	    // Vertical
	    ZED_VC  [i] += 10.0 * log10 (paramUncoded.range[i] * paramUncoded.range[i]) + paramUncoded.ZED_VC_calibration_offset;
	    ZED_XVC [i] += 10.0 * log10 (paramUncoded.range[i] * paramUncoded.range[i]) + paramUncoded.ZED_XVC_calibration_offset;
	    ZED_VCP [i] += 10.0 * log10 (paramCoded.range  [i] * paramCoded.range  [i]) + paramCoded.ZED_VCP_calibration_offset;
	    ZED_XVCP[i] += 10.0 * log10 (paramCoded.range  [i] * paramCoded.range  [i]) + paramCoded.ZED_XVCP_calibration_offset;

	}

	NPC_H [0] /= uncoded_sum_wi[0];
	NPC_V [0] /= uncoded_sum_wi[0];
	NPC_HP[0] /= coded_sum_wi[gate_offset];
	NPC_VP[0] /= coded_sum_wi[gate_offset];

	NPC_H [0] = 10.0 * log10 (NPC_H [0]);
	NPC_V [0] = 10.0 * log10 (NPC_V [0]);
	NPC_HP[0] = 10.0 * log10 (NPC_HP[0]);
	NPC_VP[0] = 10.0 * log10 (NPC_VP[0]);

	// write out variables to netCDF
	printf ("Writing dynamic variables to NetCDF...\n");
	RNC_WriteDynamicVariables (ncid, &param, &obs);
	status = nc_sync (ncid);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/* check to see if we have started a new day */
	/* the + 2L covers us for the time it takes to wrap around this loop */
	system_time = time (NULL) + 2L;
	gmtime_r (&system_time, &tm);
	if (tm.tm_mday != obs.day)
	{
	    exit_now = 1;
	    printf ("we are at the end of the day, exiting now\n");
	}
    }

    // Finish off
    printf ("*** Closing PCICARD...\n");
    RDQ_ClosePCICARD_New (amcc_fd, &dma_buffer, DMA_BUFFER_SIZE);


    /* netCDF : close the netCDF file */
    status = nc_sync (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);
    status = nc_close (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    if (param.dump_spectra  != 0)
    {
	status = nc_sync (spectra_ncid);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
	status = nc_close (spectra_ncid);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
    }

    if (param.dump_spectra_rapid  != 0)
    {
	status = nc_sync (spectra_rapid_ncid);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
	status = nc_close (spectra_rapid_ncid);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
    }

    //---------------------------
    // Unallocate all the memory
    //---------------------------
    RSP_FreeMemory (&paramCoded);  // Free memory allocated by RSP package
    RSP_FreeMemory (&paramUncoded);  // Free memory allocated by RSP package
    RSP_ObsFree (&obs); // Free observables memory
    free (timeseries);
    free (current_PSD);
    for (i = 0; i < param.samples_per_pulse; i++)
    {
	free (PSD[i].HH);
	free (PSD[i].HV);
	free (PSD[i].VV);
	free (PSD[i].VH);
	free (PSD[i].HHP);
	free (PSD[i].HVP);
	free (PSD[i].VVP);
	free (PSD[i].VHP);
    }
    free (PSD);
    free (mean_vsq_HC);
    free (mean_vsq_HCP);
    free (mean_Zsq_HC);
    free (mean_Zsq_HCP);
    free (mean_vsq_VC);
    free (mean_vsq_VCP);
    free (mean_Zsq_VC);
    free (mean_Zsq_VCP);

    free (uncoded_sum_wi);
    free (coded_sum_wi);

    free (VEL_HC_COS);
    free (VEL_HC_SIN);
    free (VEL_HCP_COS);
    free (VEL_HCP_SIN);
    free (VEL_VC_COS);
    free (VEL_VC_SIN);
    free (VEL_VCP_COS);
    free (VEL_VCP_SIN);

    free (I_coded_copolar_H);
    free (I_coded_copolar_V);
    free (Q_coded_copolar_H);
    free (Q_coded_copolar_V);
    free (I_coded_crosspolar_H);
    free (I_coded_crosspolar_V);
    free (Q_coded_crosspolar_H);
    free (Q_coded_crosspolar_V);
    free (I_uncoded_copolar_H);
    free (I_uncoded_copolar_V);
    free (Q_uncoded_copolar_H);
    free (Q_uncoded_copolar_V);
    free (I_uncoded_crosspolar_H);
    free (I_uncoded_crosspolar_V);
    free (Q_uncoded_crosspolar_H);
    free (Q_uncoded_crosspolar_V);
    free (I_uncorr_coded_copolar_H);
    free (I_uncorr_coded_copolar_V);
    free (Q_uncorr_coded_copolar_H);
    free (Q_uncorr_coded_copolar_V);
    free (I_uncorr_coded_crosspolar_H);
    free (I_uncorr_coded_crosspolar_V);
    free (Q_uncorr_coded_crosspolar_H);
    free (Q_uncorr_coded_crosspolar_V);


    free (IQStruct.I_coded_copolar_H);
    free (IQStruct.I_coded_copolar_V);
    free (IQStruct.I_coded_crosspolar_H);
    free (IQStruct.I_coded_crosspolar_V);
    free (IQStruct.Q_coded_copolar_H);
    free (IQStruct.Q_coded_copolar_V);
    free (IQStruct.Q_coded_crosspolar_H);
    free (IQStruct.Q_coded_crosspolar_V);
    free (IQStruct.I_uncoded_copolar_H);
    free (IQStruct.I_uncoded_copolar_V);
    free (IQStruct.I_uncoded_crosspolar_H);
    free (IQStruct.I_uncoded_crosspolar_V);
    free (IQStruct.Q_uncoded_copolar_H);
    free (IQStruct.Q_uncoded_copolar_V);
    free (IQStruct.Q_uncoded_crosspolar_H);
    free (IQStruct.Q_uncoded_crosspolar_V);

    fftw_destroy_plan (p_coded);
    fftw_destroy_plan (p_uncoded);
    fftw_free (in);

    //=========
    // THE END
    //=========
    printf ("All done.\n");
    exit (0);
}
