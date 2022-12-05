#define VERSION_NUMBER "0.21"

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
// 02/12/09: OTD: this is a branch from version 0.12, parameter added for switching off/on iq recording
// 28/11/11: JCN: added skewness and kurtosis
// 10/11/14: JCN: included recording of noise levels
 
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
int	exit_now = 0; 

/* function prototype declaration */
void 	sig_handler(int sig);

// Displays a welcome message with version information
void disp_welcome_message(void)
{
  char buffer[80];

  strcpy(buffer,VERSION_NUMBER);
  printf("\n");
  printf("radar-copernicus-rec: Version %s\n\n",buffer);
}

/* signal handler */
void sig_handler(int sig)
{
   if (exit_now == 0)
   {
        exit_now = 1;
        printf("***********************************\n");
        printf("* Received signal %d. Exiting soon *\n", sig);
        printf("***********************************\n");
   }
}


int parseargs(int argc, char *argv[], URC_ScanStruct *scan)
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
    strftime (scan->date, sizeof (scan->date),"%Y%m%d%H%M%S", &tm);

    return 0;
}


//--------------------------------------------------------------------
// get_config : reads radar config file
void get_config(char *filename,RSP_ParamStruct *param, URC_ScanStruct *scan, int is_coded)
{
    FILE *file;
    int i,j;
    char codefile[255];
    char dummy[80];
    int tmp_int;

    printf("Accessing config file: %s\n", filename);

    param->frequency                  = RNC_GetConfigDouble(filename,"radar-frequency");
    param->prf                        = RNC_GetConfigDouble(filename,"prf");
    param->transmit_power             = RNC_GetConfigDouble(filename,"transmit-power");
    param->pulses_per_daq_cycle       = RNC_GetConfigDouble(filename,"pulses");
    param->samples_per_pulse          = RNC_GetConfigDouble(filename,"samples");
    param->ADC_channels               = RNC_GetConfigDouble(filename,"adc-channels");
    param->clock_divfactor            = RNC_GetConfigDouble(filename,"adc-divfactor");
    param->delay_clocks               = RNC_GetConfigDouble(filename,"adc-delayclocks");
    param->pulse_period               = RNC_GetConfigDouble(filename,"chip-length");
    param->pulses_coherently_averaged = RNC_GetConfigDouble(filename,"num-coh-avg");
    param->spectra_averaged           = RNC_GetConfigDouble(filename,"num-spec-avg");
    param->moments_averaged           = RNC_GetConfigDouble(filename,"num-moments-avg");
    param->fft_bins_interpolated      = RNC_GetConfigDouble(filename,"reject-clutter-bins");
    param->clock                      = RNC_GetConfigDouble(filename,"adc-clock");
    param->num_peaks                  = RNC_GetConfigDouble(filename,"num-peaks");
    param->antenna_diameter 	      = RNC_GetConfigFloat(filename, "antenna_diameter");
    param->beamwidthH		      = RNC_GetConfigFloat(filename, "beamwidthH");
    param->beamwidthV		      = RNC_GetConfigFloat(filename, "beamwidthV");
    param->height		      = RNC_GetConfigFloat(filename, "height");
    param->azimuth_offset	      = RNC_GetConfigFloat(filename, "azimuth_offset");
    param->dump_spectra   	      = RNC_GetConfigFloat(filename, "dump_spectra");
    param->dump_spectra_rapid         = RNC_GetConfigFloat(filename, "dump_spectra_rapid");
    param->num_interleave             = RNC_GetConfigFloat(filename, "num-interleave");
    param->num_tx_pol                 = RNC_GetConfigFloat(filename, "num-tx-pol");
    param->include_iq_in_spectra      = RNC_GetConfigFloat(filename, "include_iq_in_spectra");

    scan->min_angle=-9999;
    scan->max_angle=-9999;
    scan->scan_angle=-9999;
	
    scan->scan_angle		      = RNC_GetConfigFloat(filename, "antenna_azimuth");

    RNC_GetConfig(filename,"code-file",codefile, sizeof (codefile));

    strcpy(param->code_name, "NOT YET IMPLEMENTED\0"); 
    if (strchr(codefile, '\n') != NULL)
      {
        codefile[strlen(codefile)-1] = 0x0;
      }

    if(is_coded == 1)
      {
	printf("Reading pulse codes from file: ");
	printf("*%s*\n",codefile);
	
	file=fopen(codefile,"r");
	if(file==NULL)
	  {
	    printf("** ERROR: Unable to open pulse code file!\n");
	    exit(1);
	  }
	if (fscanf(file,"%79s %d",dummy,&param->code_length) != 2)
	    exit (3);
	if (fscanf(file,"%79s %d",dummy,&param->number_of_codes) != 2)
	    exit (3);
	
	for(i=0; i<param->number_of_codes; i++) 
	  {
	    for(j=0; j<param->code_length; j++) 
	      {
		  if (fscanf(file,"%d",&tmp_int) != 1)
		      exit (3);
		param->codes[i][j]=(short)tmp_int;
	      }
	  }
	
	fclose(file);
	if(param->num_interleave>1 && param->number_of_codes<2)
	  {
	    printf("** get_config: pulse code file doesn't contain enough pulses!");
	    exit(1);
	  } 
	else 
	  {
	    param->number_of_codes = param->number_of_codes / param->num_interleave;
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
void get_cal(RSP_ParamStruct *param, char *calfile)
{
  param->ZED_calibration_offset=RNC_GetConfigDouble(calfile,"z-calibration");
  param->ZDR_calibration_offset=RNC_GetConfigDouble(calfile,"zdr-calibration");
  param->LDR_calibration_offset=RNC_GetConfigDouble(calfile,"ldr-calibration");
  param->ZED_incoherent_calibration_offset=RNC_GetConfigDouble(calfile,"z-incoherent-calibration");
  param->ZED_incoherent_noise=RNC_GetConfigDouble(calfile,"z-incoherent-noise");
  param->range_offset=RNC_GetConfigDouble(calfile,"range-offset");
}


//--------------------------------------------------------------------
/* make_dmux_table generates the lookup table that is used to extract
 * channels from the DMA buffer 
 * IN:  channels  the number of channels
 * OUT: nowt */
void make_dmux_table( int channels)
{
  if (channels == 4)
    {
      dmux_table[0] = 0;
      dmux_table[1] = 2;
      dmux_table[2] = 1;
      dmux_table[3] = 3;
      dmux_table[4] = 0xFFFFFFFF; /* channels 4 to 7 do not exist in the eight channels system */
      dmux_table[5] = 0xFFFFFFFF;
      dmux_table[6] = 0xFFFFFFFF;
      dmux_table[7] = 0xFFFFFFFF;
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


//========================= M A I N   C O D E =======================
//            [ See disp_help() for command-line options ]
//-------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int       num_pulses;
  int       amcc_fd = 0;        /* file descriptor for the PCICARD */
  caddr_t   dma_buffer = NULL;  /* size of dma buffer */
  uint16_t *dma_banks[2];      /* pointers to the dma banks */
  int       dma_bank = 0;
  int       proc_bank = 1;
  int       tcount;             /* number of bytes to be transferred during the DMA */
  uint16_t *data;
  int       count,sample;
  int       nspectra;
  int       status;
  long      num_data;
  float     *current_PSD;
  register  int  i,j,k;
  int       temp_int = 0;
  float     HH_moments[RSP_MOMENTS];
  float     HV_moments[RSP_MOMENTS];
  float     HHP_moments[RSP_MOMENTS];
  float     HVP_moments[RSP_MOMENTS];
  float     HH_noise_level, HV_noise_level;
  float     HHP_noise_level, HVP_noise_level;
  float     TX1_level, TX2_level;
  int       noisegate1, noisegate2;
  int       TX1_gate1, TX1_gate2, TX2_gate1, TX2_gate2;
  time_t    system_time;
  time_t    spectra_time = 0;
  time_t    spectra_rapid_time = 0;
  time_t    temp_time_t;
  
  struct tm tm;
  char   datestring[25];
  float *uncoded_mean_vsq; // Used in sigma vbar calculation
  float *coded_mean_vsq; // Used in sigma vbar calculation
  float *uncoded_mean_Zsq; // Used in sigma Zbar calculation
  float *coded_mean_Zsq; // Used in sigma Zbar calculation
  short int   xcodes[32][64];

  URC_ScanStruct scan;
  RNC_DimensionStruct dimensions;
  PolPSDStruct * PSD;
  IQStruct	IQStruct;

  /* netCDF file pointer */
  int ncid;
  int spectra_ncid       = -1;
  int spectra_rapid_ncid = -1;
  int PSD_varid[PSD_varidSize];
  int PSD_rapid_varid[RapidPSD_varidSize];
  int file_stateid = 0;

  int pulse = 0;
  int marker = 0;
  int counter = 0;

  float norm_coded, norm_uncoded;

  /* signal */
  struct      sigaction sig_struct;

  int gate_offset; // Offset of coded gates rel. to uncoded

  int pulse_increment;

  int tempIco,tempQco,tempIcr,tempQcr;
  int nm;

  // The following are shortcut pointers to the elements of
  // the obs structure
  float *ZED_H;
  float *ZED_HC, *ZED_HCP;
  //float *ZED_VC, *ZED_VCP;
  float *ZED_XHC, *ZED_XHCP;
  float *SNR_HC, *SNR_HCP;
  //float *SNR_VC, *SNR_VCP;
  float *SNR_XHC, *SNR_XHCP;
  float *VEL_HC, *VEL_HCP;
  float *VEL_HCD, *VEL_HCDP;
  float *ZED_HCD, *ZED_HCDP;
  //float *VEL_VC;
  float *SPW_HC, *SPW_HCP;
  float *SKW_HC, *SKW_HCP;
  float *KRT_HC, *KRT_HCP;
  //float *SPW_VC;
  float *ZDR_C, *ZDR_CP;
  float *LDR_C, *LDR_CP;

  float *VEL_HC_COS, *VEL_HC_SIN, *VEL_HCP_COS, *VEL_HCP_SIN;
  float *NPC_H, *NPC_HP, *NPC_V, *NPC_VP;

  float *TX_1H, *TX_2H;

  float wi;              // Individual weighting value
  float *uncoded_sum_wi; // Sum of weighting values
  float *coded_sum_wi;   // Sum of weighting values

  int offset, ipulse;
  int pos1,pos2;
  int tot_n_avg,gtna,gate;

  int	obtain_index;
  int	store_index;
  long int *I_coded_copolar_H;
  long int *Q_coded_copolar_H;
  long int *I_coded_crosspolar_H;
  long int *Q_coded_crosspolar_H;
  uint16_t *I_uncoded_copolar_H;
  uint16_t *Q_uncoded_copolar_H;
  uint16_t *I_uncoded_crosspolar_H;
  uint16_t *Q_uncoded_crosspolar_H;
  uint16_t *log_raw;
  uint16_t *TX1_raw, *TX2_raw;

  int 	collect_spectra_now;
  int	collect_spectra_rapid_now;

  RSP_ParamStruct param, paramCoded, paramUncoded;
  RSP_ComplexType *timeseries;
  RSP_ObservablesStruct obs;
  RSP_ObservablesStruct PSD_obs;
  RSP_ObservablesStruct PSD_RAPID_obs;

  RSP_PeakStruct * HH_peaks;
  RSP_PeakStruct * HV_peaks;
  RSP_PeakStruct * HHP_peaks;
  RSP_PeakStruct * HVP_peaks;

  static REL_SerialMessageStruct clinometermsg;

  fftw_complex *in;
  fftw_plan p_coded,p_uncoded;

  /* time variables */
  struct timeval          tv;
  struct timezone         tz;


  disp_welcome_message();


  //---------------------------
  // Set up the signal handler
  //---------------------------
  /* Set up the sig_struct variable */
  sig_struct.sa_handler = sig_handler;
  sigemptyset( &sig_struct.sa_mask );
  sig_struct.sa_flags = 0;
  /* Install signal handler and check for error */
  if (sigaction(SIGINT, &sig_struct, NULL) != 0)
    {
      perror ("Error installing signal handler\n");
      exit(1);
    }
  

  //------------------------------
  // Parse command line arguments
  //------------------------------
  if( parseargs(argc,argv,&scan) == 1 )
    exit(1);


  //------------------------
  // Read radar config file
  //------------------------
  get_config(CONFIG_FILE,&paramCoded,&scan,1);   // Do it for coded pulses
  get_config(CONFIG_FILE,&paramUncoded,&scan,0); // Do it for uncoded pulses

  /* match the scan min_angle and max_angle to the measured elevation angle */
  /* initialise the serial port */
  status = REL_InitialiseSerialMessage(CLINOMETERMESSAGE_PORT);
  /* a slight pause to allow things to settle */
  sleep(1);
  if ( status != 0) 
    {
      printf("Detected a problem with initialising the serial port\n");
    } 
  else 
    {
      status = REL_ReadSerialMessage(&clinometermsg);
      if (status == 0) 
	{
	  /* we need to apply a 90 deg elevation offset */
	  clinometermsg.el = clinometermsg.el + 90.0;
	  printf("elevation angle : %f (%d)\n", clinometermsg.el, temp_int);
	  scan.min_angle = clinometermsg.el;
	  scan.max_angle = clinometermsg.el;
	}
    }

  // Read calibration file
  // We will eventaully need separate cal. files for coded and uncoded
  get_cal(&paramCoded, CAL_FILE_CODED);
  get_cal(&paramUncoded, CAL_FILE_UNCODED);
  
  //------------------------------------
  // Initialise RSP parameter structure
  //------------------------------------
  paramCoded.prf /= (paramCoded.num_interleave * paramCoded.num_tx_pol);  // Calculate effective PRF
  paramUncoded.prf /= (paramCoded.num_interleave * paramCoded.num_tx_pol);
  RSP_InitialiseParams(&paramCoded);   // This param is used for coded pulses
  RSP_InitialiseParams(&paramUncoded); // This param is used for uncoded pulses

  //// This botches the DAQ to oversample
  //paramCoded.oversample_ratio=2;
  //paramUncoded.oversample_ratio=2;

  printf("Parameters for Coded mode:\n");
  RSP_DisplayParams(&paramCoded);
  printf("Parameters for Uncoded mode:\n");
  RSP_DisplayParams(&paramUncoded);

  param = paramUncoded; // This param is used for params that are the same for both modes

  pulse_increment=2;  // This is because pulses are interleaved
 
  // Sample extra pulses at end so that we have entire code sequence
  num_pulses = (int)(paramCoded.pulses_per_daq_cycle+(paramCoded.number_of_codes-1)*paramCoded.num_interleave);
  tcount = param.spectra_averaged * num_pulses * param.samples_per_pulse * param.ADC_channels * sizeof(uint16_t);

  // Number of data points to allocate per data stream
  num_data = paramCoded.pulses_per_daq_cycle*paramCoded.samples_per_pulse;

  // Allocate memory for coded and uncoded data streams
  I_coded_copolar_H          = malloc(sizeof(long int) * num_data/paramCoded.num_interleave);
  Q_coded_copolar_H          = malloc(sizeof(long int) * num_data/paramCoded.num_interleave);
  IQStruct.I_coded_copolar_H = malloc(sizeof(long int) * paramCoded.samples_per_pulse * paramCoded.nfft * param.spectra_averaged);
  IQStruct.Q_coded_copolar_H = malloc(sizeof(long int) * paramCoded.samples_per_pulse * paramCoded.nfft * param.spectra_averaged);

  I_coded_crosspolar_H          = malloc(sizeof(long int) * num_data/paramCoded.num_interleave);
  Q_coded_crosspolar_H          = malloc(sizeof(long int) * num_data/paramCoded.num_interleave);
  IQStruct.I_coded_crosspolar_H = malloc(sizeof(long int) * paramCoded.samples_per_pulse * paramCoded.nfft * param.spectra_averaged);
  IQStruct.Q_coded_crosspolar_H = malloc(sizeof(long int) * paramCoded.samples_per_pulse * paramCoded.nfft * param.spectra_averaged);

  I_uncoded_copolar_H          = malloc(sizeof(uint16_t) * num_data);
  Q_uncoded_copolar_H          = malloc(sizeof(uint16_t) * num_data);
  IQStruct.I_uncoded_copolar_H = malloc(sizeof(uint16_t) * paramUncoded.samples_per_pulse * paramUncoded.nfft * param.spectra_averaged);
  IQStruct.Q_uncoded_copolar_H = malloc(sizeof(uint16_t) * paramUncoded.samples_per_pulse * paramUncoded.nfft * param.spectra_averaged);

  I_uncoded_crosspolar_H          = malloc(sizeof(uint16_t) * num_data);
  Q_uncoded_crosspolar_H          = malloc(sizeof(uint16_t) * num_data);
  IQStruct.I_uncoded_crosspolar_H = malloc(sizeof(uint16_t) * paramUncoded.samples_per_pulse * paramUncoded.nfft * param.spectra_averaged);
  IQStruct.Q_uncoded_crosspolar_H = malloc(sizeof(uint16_t) * paramUncoded.samples_per_pulse * paramUncoded.nfft * param.spectra_averaged);


  log_raw = malloc(sizeof(uint16_t) * num_data);
  TX1_raw = malloc(sizeof(uint16_t) * num_data);
  TX2_raw = malloc(sizeof(uint16_t) * num_data);

  HH_peaks  = calloc (param.num_peaks, sizeof (RSP_PeakStruct));
  HV_peaks  = calloc (param.num_peaks, sizeof (RSP_PeakStruct));
  HHP_peaks = calloc (param.num_peaks, sizeof (RSP_PeakStruct));
  HVP_peaks = calloc (param.num_peaks, sizeof (RSP_PeakStruct));

  in        = fftw_malloc(sizeof(fftw_complex)*param.nfft);
  p_coded   = fftw_plan_dft_1d(paramCoded.nfft,in,in,FFTW_FORWARD,FFTW_ESTIMATE);
  p_uncoded = fftw_plan_dft_1d(paramUncoded.nfft,in,in,FFTW_FORWARD,FFTW_ESTIMATE);

  timeseries=malloc(param.nfft*sizeof(RSP_ComplexType));

  current_PSD=malloc(param.npsd*sizeof(float));
  PSD         = calloc (param.samples_per_pulse, sizeof (PolPSDStruct));

  for(j=0;j<param.samples_per_pulse;j++) 
    {
      PSD[j].HH=malloc(param.npsd*sizeof(float));   // not coded
      PSD[j].HV=malloc(param.npsd*sizeof(float));   // not coded
      PSD[j].HHP=malloc(param.npsd*sizeof(float));  // coded
      PSD[j].HVP=malloc(param.npsd*sizeof(float));  // coded
    }
  coded_mean_vsq   = malloc(paramCoded.samples_per_pulse*sizeof(float));
  uncoded_mean_vsq = malloc(paramUncoded.samples_per_pulse*sizeof(float));
  coded_mean_Zsq   = malloc(paramCoded.samples_per_pulse*sizeof(float));
  uncoded_mean_Zsq = malloc(paramUncoded.samples_per_pulse*sizeof(float));
  coded_sum_wi     = malloc(paramCoded.samples_per_pulse*sizeof(float));
  uncoded_sum_wi   = malloc(paramUncoded.samples_per_pulse*sizeof(float));
  VEL_HC_COS       = malloc(paramUncoded.samples_per_pulse*sizeof(float));
  VEL_HC_SIN       = malloc(paramUncoded.samples_per_pulse*sizeof(float));
  VEL_HCP_COS      = malloc(paramCoded.samples_per_pulse*sizeof(float));
  VEL_HCP_SIN      = malloc(paramCoded.samples_per_pulse*sizeof(float));

  norm_coded=1/paramCoded.Wss;
  norm_uncoded=1/paramUncoded.Wss;

  //----------------------------
  // Initialise RSP Observables
  //----------------------------
  RSP_ObsInit(&obs);
  RSP_ObsInit(&PSD_obs);
  RSP_ObsInit(&PSD_RAPID_obs);

  // Power monitoring
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"TX_1H");
  TX_1H = RSP_ObsNew(&obs, "TX_1H", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"TX_2H");
  TX_2H = RSP_ObsNew(&obs, "TX_2H", param.samples_per_pulse, temp_int);

  // Last argument determines whether the parameter will be recorded or not
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"ZED_H");
  ZED_H    = RSP_ObsNew(&obs, "ZED_H", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"ZED_HC");
  ZED_HC   = RSP_ObsNew(&obs, "ZED_HC", param.samples_per_pulse, temp_int);
  //temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"ZED_VC");
  //ZED_VC  = RSP_ObsNew(&obs, "ZED_VC", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"ZED_XHC");
  ZED_XHC  = RSP_ObsNew(&obs, "ZED_XHC", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"SNR_HC");
  SNR_HC   = RSP_ObsNew(&obs, "SNR_HC", param.samples_per_pulse, temp_int);
  //temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"SNR_VC");
  //SNR_VC  = RSP_ObsNew(&obs, "SNR_VC", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"SNR_XHC");
  SNR_XHC  = RSP_ObsNew(&obs, "SNR_XHC", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"VEL_HC");
  VEL_HC   = RSP_ObsNew(&obs, "VEL_HC", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"VEL_HCD");
  VEL_HCD  = RSP_ObsNew(&obs, "VEL_HCD", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"ZED_HCD");
  ZED_HCD  = RSP_ObsNew(&obs, "ZED_HCD", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"SPW_HC");
  SPW_HC   = RSP_ObsNew(&obs, "SPW_HC", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"SKW_HC");
  SKW_HC   = RSP_ObsNew(&obs, "SKW_HC", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"KRT_HC");
  KRT_HC   = RSP_ObsNew(&obs, "KRT_HC", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"ZDR_C");
  ZDR_C    = RSP_ObsNew(&obs, "ZDR_C", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"LDR_C");
  LDR_C    = RSP_ObsNew(&obs, "LDR_C", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"ZED_HCP");
  ZED_HCP  = RSP_ObsNew(&obs, "ZED_HCP", param.samples_per_pulse, temp_int);
  //temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"ZED_VCP");
  //ZED_VCP  = RSP_ObsNew(&obs, "ZED_VCP", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"ZED_XHCP");
  ZED_XHCP = RSP_ObsNew(&obs, "ZED_XHCP", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"SNR_HCP");
  SNR_HCP  = RSP_ObsNew(&obs, "SNR_HCP", param.samples_per_pulse, temp_int);
  //temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"SNR_VCP");
  //SNR_VCP  = RSP_ObsNew(&obs, "SNR_VCP", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"SNR_XHCP");
  SNR_XHCP = RSP_ObsNew(&obs, "SNR_XHCP", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"VEL_HCP");
  VEL_HCP  = RSP_ObsNew(&obs, "VEL_HCP", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"VEL_HCDP");
  VEL_HCDP = RSP_ObsNew(&obs, "VEL_HCDP", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"ZED_HCDP");
  ZED_HCDP = RSP_ObsNew(&obs, "ZED_HCDP", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"SPW_HCP");
  SPW_HCP  = RSP_ObsNew(&obs, "SPW_HCP", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"SKW_HCP");
  SKW_HCP  = RSP_ObsNew(&obs, "SKW_HCP", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"KRT_HCP");
  KRT_HCP  = RSP_ObsNew(&obs, "KRT_HCP", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"ZDR_CP");
  ZDR_CP   = RSP_ObsNew(&obs, "ZDR_CP", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"LDR_CP");
  LDR_CP   = RSP_ObsNew(&obs, "LDR_CP", param.samples_per_pulse, temp_int);
  
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"NPC_H");
  NPC_H    = RSP_ObsNew(&obs, "NPC_H", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"NPC_V");
  NPC_V    = RSP_ObsNew(&obs, "NPC_V", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"NPC_HP");
  NPC_HP   = RSP_ObsNew(&obs, "NPC_HP", param.samples_per_pulse, temp_int);
  temp_int = (int)RNC_GetConfigDouble(CONFIG_FILE,"NPC_VP");
  NPC_VP   = RSP_ObsNew(&obs, "NPC_VP", param.samples_per_pulse, temp_int);

  //if(args.quiet==0)
  //  {
  printf("Recording observables:");
  for(i=0; i<obs.n_obs; i++)
    {
      if(obs.record_observable[i]==1)
	{
	  printf(" %s",obs.name[i]);
	}
    }
  printf("\n");
  //  }
  

  //-------------------------------------------
  // Set up the data acquisition 
  //-------------------------------------------
  printf("** Initialising ISACTRL...\n");
  RDQ_InitialiseISACTRL( num_pulses, param.samples_per_pulse,
			 param.clock_divfactor, param.delay_clocks );
  
  printf("** Initialising PCICARD...\n");
  amcc_fd = RDQ_InitialisePCICARD_New( &dma_buffer, DMA_BUFFER_SIZE );
  
  // Initialise pointers to DMA banks
  dma_banks[ 0 ] = (uint16_t *) dma_buffer;
  dma_banks[ 1 ] = (uint16_t *) (dma_buffer + (DMA_BUFFER_SIZE/2));
  
  make_dmux_table(param.ADC_channels);
  
  printf("** Starting acquisition...\n");

  RDQ_StartAcquisition( amcc_fd, dma_bank, (short *)(dma_banks[dma_bank]), tcount);

  /* setup the netCDF file */
  ncid = RNC_OpenNetcdfFile (GetRadarName (COPERNICUS),
			     GetSpectraName (COPERNICUS),
			     scan.date, NULL,
			     GetScanTypeName (scan.scanType),
			     GetSpectraExtension (COPERNICUS), "raw");
  RNC_SetupDimensions( ncid, &param, &dimensions );
  RNC_SetupGlobalAttributes( ncid, COPERNICUS, &scan, &param, argc, argv );
  printf("netCDF : global attributes have been defined.\n");
  RNC_SetupPulse_Compression_Code( ncid, &param );
  file_stateid = RNC_SetupFile_State(ncid);
  RNC_SetupStaticVariables( ncid, &param);
  printf("netCDF : static variables have been defined\n");
  RNC_SetupRange(ncid, &param, &dimensions);
  RNC_SetupDynamicVariables(ncid, COPERNICUS, &scan, &param, &dimensions, &obs); 
  printf("netCDF : dynamic variables have been defined\n");
  /* change the mode of netCDF from define to data */
  status = nc_enddef(ncid);
  if (status != NC_NOERR) check_netcdf_handle_error(status);

  // Set up spectral dump file
  if (param.dump_spectra != 0) 
    {
      printf("setup spectra recording file\n");
      spectra_ncid = RNC_OpenNetcdfFile
	  (GetRadarName (COPERNICUS_SPECTRA),
	   GetSpectraName (COPERNICUS_SPECTRA),
	   scan.date, NULL,
	   GetScanTypeName (scan.scanType),
	   GetSpectraExtension (COPERNICUS_SPECTRA), "raw");
      RNC_SetupDimensions( spectra_ncid, &param, &dimensions );
      RNC_SetupGlobalAttributes( spectra_ncid, COPERNICUS, &scan, &param, argc, argv );
      RNC_SetupPulse_Compression_Code( spectra_ncid, &param );
      RNC_SetupStaticVariables( spectra_ncid, &param);
      RNC_SetupRange(spectra_ncid, &param, &dimensions);
      RNC_SetupDynamicVariables( spectra_ncid, COPERNICUS_SPECTRA, &scan, &param, &dimensions, &PSD_obs);
      RNC_SetupLogPSDVariables( spectra_ncid, COPERNICUS_CODED_SPECTRA, &param, &dimensions, PSD_varid );
      /* change the mode of netCDF from define to data */
      status = nc_enddef(spectra_ncid);
      if (status != NC_NOERR) check_netcdf_handle_error(status);
      time(&temp_time_t);
      spectra_time = param.dump_spectra * (floorl(temp_time_t/param.dump_spectra) + 1);
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
      RNC_SetupRapidLogPSDDimensions( spectra_rapid_ncid, COPERNICUS_SPECTRA_RAPID, &param, &dimensions );
      RNC_SetupGlobalAttributes( spectra_rapid_ncid, COPERNICUS, &scan, &param, argc, argv );
      RNC_SetupStaticVariables( spectra_rapid_ncid, &param);
      RNC_SetupRange(spectra_rapid_ncid, &param, &dimensions);
      RNC_SetupDynamicVariables( spectra_rapid_ncid, COPERNICUS_SPECTRA_RAPID, &scan, &param, &dimensions, &PSD_RAPID_obs);
      RNC_SetupLogPSDVariables( spectra_rapid_ncid, COPERNICUS_SPECTRA_RAPID, &param, &dimensions, PSD_rapid_varid );
      /* change the mode of netCDF from define to data */
      status = nc_enddef(spectra_rapid_ncid);
      if (status != NC_NOERR) check_netcdf_handle_error(status);
      
      time(&temp_time_t);
      spectra_rapid_time = param.dump_spectra_rapid * (floorl(temp_time_t/param.dump_spectra_rapid) + 1);
    }



  //------------------------------------------------------------------
  // "Oversample" the pulse code templates to match the sampling rate
  //------------------------------------------------------------------
  for(i=0; i<paramCoded.number_of_codes*paramCoded.num_interleave; i++)
    {
    for(j=0; j<paramCoded.code_length; j++) printf("%i ",paramCoded.codes[i][j]);
    printf("\n");
    RSP_Oversample(paramCoded.codes[i],xcodes[i],paramCoded.code_length,paramCoded.oversample_ratio);
    for(j=0; j<paramCoded.code_length*paramCoded.oversample_ratio; j++) printf("%i ",xcodes[i][j]);
    printf("\n");
    }
  

  gate_offset = paramCoded.code_length * paramCoded.oversample_ratio - 1; 
  

  // THIS IS THE START OF THE OUTER RAY LOOP
  while( exit_now == 0)
    {
      printf("<< PRESS CTRL-C TO EXIT >>\n");
      
      /* get time of day */
      gettimeofday( &tv, &tz);
      gmtime_r (&tv.tv_sec, &tm);
      printf("System time: %s\n", asctime (&tm));
      obs.year = tm.tm_year+1900;
      obs.month= tm.tm_mon+1;
      obs.day  = tm.tm_mday;
      obs.hour = tm.tm_hour;
      obs.minute=tm.tm_min;
      obs.second=tm.tm_sec;
      obs.centisecond=(int)tv.tv_usec/10000;
      sprintf(datestring,"%04d/%02d/%02d %02d:%02d:%02d.%02d",obs.year,obs.month,obs.day,obs.hour,obs.minute,obs.second, obs.centisecond);
      printf("Date string: %s",datestring);
      
      obs.azimuth = scan.scan_angle;
      PSD_obs.azimuth = obs.azimuth;
      PSD_RAPID_obs.azimuth = obs.azimuth;
      /* read the elevation angle from the clionometer */
      status = REL_ReadSerialMessage(&clinometermsg);

      if (status == 0) 
	{
	  /* we need to apply a 90 deg elevation offset */
	  clinometermsg.el = clinometermsg.el + 90.0;
	  printf("elevation angle : %f (%d)\n", clinometermsg.el, temp_int);
	  obs.elevation = clinometermsg.el;
	} 
      else 
	{
	  obs.elevation = -999;
	}

      PSD_obs.elevation = obs.elevation;
      PSD_RAPID_obs.elevation = obs.elevation;
      
      
      
      // Initialise observables to zero
      // (Needed for moments averaging)
      for(i=0; i<obs.n_obs; i++)
	{
	  printf("Initialising %s\n",obs.name[i]);
	  for(j=0; j<obs.n_elements [i]; j++)
	    {
	      obs.data[i][j]=0.0f;
	    }
	}
      for(j=0; j<param.samples_per_pulse; j++)
        {
	  uncoded_mean_vsq[j] = 0.0f;
	  coded_mean_vsq[j]   = 0.0f;
	  uncoded_mean_Zsq[j] = 0.0f;
	  coded_mean_Zsq[j]   = 0.0f;
	  uncoded_sum_wi[j]   = 0.0f;
	  coded_sum_wi[j]     = 0.0f;
          VEL_HC_COS[j]       = 0.0f;
          VEL_HC_SIN[j]       = 0.0f;
          VEL_HCP_COS[j]      = 0.0f;
          VEL_HCP_SIN[j]      = 0.0f;
	}
      
      printf("Done initialising variables...\n");
      
      int TX_counter = 0;

      // LOOP THROUGH MOMENTS AVERAGING FROM HERE...
      for(nm=0; nm<param.moments_averaged; nm++) 
	{
	  
	  for(j=0; j<param.samples_per_pulse; j++)
	    {
	      register int bin_no;
	      // Initialise spectra to zero
	      for(bin_no=0; bin_no<param.npsd; bin_no++)
		{
		  PSD[j].HH[bin_no]  = 0.0f;
		  PSD[j].HV[bin_no]  = 0.0f;
		  PSD[j].HHP[bin_no] = 0.0f;
		  PSD[j].HVP[bin_no] = 0.0f;
		}
	    }
	  memset (dma_banks[proc_bank], -1, tcount);

	  /* waiting for DMA to complete */
	  status = RDQ_WaitForAcquisitionToComplete( amcc_fd );
	  if (status != 0) printf("There was a problem in WaitForAcquisitionToComplete\n");
	  
	  //----------------------------------------------------------------
	  // Swap around the areas used for storing daq and processing from
	  //----------------------------------------------------------------
	  dma_bank  = 1 - dma_bank;
	  proc_bank = 1 - proc_bank;
	  
	  data = dma_banks[proc_bank];
	  RDQ_StartAcquisition2( amcc_fd, dma_bank, tcount);
	  
	  //---------------------------------
	  // Loop through spectral averaging
	  //---------------------------------
	  collect_spectra_rapid_now = 0;
	  collect_spectra_now = 0;
	  system_time = time(NULL);
	  if (param.dump_spectra_rapid != 0)
	    {
	      if ( spectra_rapid_time <= system_time ) 
		{
		  collect_spectra_rapid_now = 1;
		}
	      
	    }
	  /* write out spectra to netCDF if required */
	  if (param.dump_spectra != 0)
	    {
	      if ( spectra_time <= system_time ) 
		{
		  collect_spectra_now = 1;
		}
	    }
	  
	  for(nspectra=0;nspectra<param.spectra_averaged;nspectra++)
	    {
	      printf("\nAveraging %d of %d spectra:\n",nspectra+1,param.spectra_averaged);
	      
	      ipulse = 0;
	      offset = -1;
	      
	      //----------------------------------------------------------------
	      // Extract data from DMA memory (still interleaved at this stage)
	      //----------------------------------------------------------------
	      for (i = 0; i < num_pulses; i++) 
		{
		  register int count_reg;
		  for (j = 0; j < param.samples_per_pulse; j++) 
		    {
			uint16_t Sync;

		      count_reg=ipulse*param.samples_per_pulse+j;
		      I_uncoded_copolar_H[count_reg]    = GET_CHANNEL( data, CHAN_Ic );
		      Q_uncoded_copolar_H[count_reg]    = GET_CHANNEL( data, CHAN_Qc );
		      I_uncoded_crosspolar_H[count_reg] = GET_CHANNEL( data, CHAN_Ix );
		      Q_uncoded_crosspolar_H[count_reg] = GET_CHANNEL( data, CHAN_Qx );
		      
		      log_raw[count_reg]                = GET_CHANNEL( data, CHAN_INC);
		      TX1_raw[count_reg]                = GET_CHANNEL( data, CHAN_T1);
		      TX2_raw[count_reg]                = GET_CHANNEL( data, CHAN_T2);
		      
		      Sync = GET_CHANNEL( data, CHAN_SYNC );
		      //	printf("Sync pulse = %d at pulse %i **********\n",Sync,i);
		      
		      // Only increment the pulse counter if a sync has been found
		      // This will discard all data before the first sync pulse
		      if(offset == -1 && Sync > 2000)
			{
			  offset = i;
			  //	      printf("Sync pulse found at pulse %i **********\n",offset+1);
			}
		      
		      INC_POINTER ( data, param.ADC_channels);
		    }
		  if(offset>-1) ipulse++;
		} 
	      
	      if(offset==-1)
		{
		  printf("** ERROR: NO SYNC PULSE FOUND!\n");
		  exit(1);
		}
	      
	      
	      // NB: Pulse sequence should be as follows:
	      // Code A : Uncoded : Code B : Uncoded : Code A : ...
	      // SYNC ^ :         :        :         : SYNC ^ : ...
	      
	      printf("perform RSP_Correlate\n");
	      // DECODE PULSES AND AVERAGE THE LOG CHANNEL
	      // Loop through coded pulses
	      
	      count=paramCoded.code_length*paramCoded.oversample_ratio;
	      for(i=0; i<paramCoded.pulses_per_daq_cycle; i=i+paramCoded.num_interleave)
		{
		  pos1=param.samples_per_pulse*i;
		  pos2=param.samples_per_pulse*(i/paramCoded.num_interleave);
		  k=i%(paramCoded.number_of_codes*paramCoded.num_interleave);
		  RSP_Correlate(&I_uncoded_copolar_H[pos1],xcodes[k],param.samples_per_pulse,count,&I_coded_copolar_H[pos2]);
		  RSP_Correlate(&Q_uncoded_copolar_H[pos1],xcodes[k],param.samples_per_pulse,count,&Q_coded_copolar_H[pos2]);
		  RSP_Correlate(&I_uncoded_crosspolar_H[pos1],xcodes[k],param.samples_per_pulse,count,&I_coded_crosspolar_H[pos2]);
		  RSP_Correlate(&Q_uncoded_crosspolar_H[pos1],xcodes[k],param.samples_per_pulse,count,&Q_coded_crosspolar_H[pos2]);
		  
		  // Average the incoherent channel (uncoded only)
		  for(j=0; j<paramUncoded.samples_per_pulse; j++)
		    ZED_H[j] += (float)log_raw[i+1] / ((float)paramUncoded.pulses_per_daq_cycle*paramUncoded.spectra_averaged*paramUncoded.moments_averaged/2);
		  
		}
	      
	      printf("RSP_Correlate complete\n");
	      
	      // Up-sample to 60m gates here
	      
	      //------------------------------------
	      // Add up complementary code sequence
	      //------------------------------------
	      tot_n_avg=paramCoded.number_of_codes*paramCoded.pulses_coherently_averaged;
	      for(i=0; i<paramCoded.nfft; i++)
		{
		  register int pos2_reg;
		  gtna=i*tot_n_avg;  // Work this out first to speed things up later
		  count=paramCoded.samples_per_pulse*i;
		  for(j=0; j<paramCoded.samples_per_pulse; j++)
		    {
		      pos2_reg=count+j;
		      tempIco=0; tempQcr=0; tempIcr=0; tempQco=0;
		      for(k=0; k<tot_n_avg; k++)
			{
			  register int pos1_reg;
			  pos1_reg=paramCoded.samples_per_pulse*(gtna+k)+j;
			  tempIco += I_coded_copolar_H[pos1_reg];
			  tempQco += Q_coded_copolar_H[pos1_reg];
			  tempIcr += I_coded_crosspolar_H[pos1_reg];
			  tempQcr += Q_coded_crosspolar_H[pos1_reg];
			}
		      I_coded_copolar_H[pos2_reg]    = tempIco/tot_n_avg;
		      Q_coded_copolar_H[pos2_reg]    = tempQco/tot_n_avg;
		      I_coded_crosspolar_H[pos2_reg] = tempIcr/tot_n_avg;
		      Q_coded_crosspolar_H[pos2_reg] = tempQcr/tot_n_avg;
		    }
		}
	      
	      if (collect_spectra_now == 1) 
		{
		  printf("storing IQs\n");
		  /* store I and Q for each pulse */
		  /* nspectra defines the spectra number */
		  /* lets do the code first */
		  for (i = 0; i < paramCoded.nfft; i++) 
		    {
		      obtain_index = i * paramCoded.samples_per_pulse;
		      store_index = (i * paramCoded.samples_per_pulse) + (nspectra * paramCoded.samples_per_pulse * paramCoded.nfft);
		      for (j=0; j<paramCoded.samples_per_pulse; j++) 
			{
			  IQStruct.I_coded_copolar_H[store_index]    = I_coded_copolar_H[obtain_index];
			  IQStruct.Q_coded_copolar_H[store_index]    = Q_coded_copolar_H[obtain_index];
			  IQStruct.I_coded_crosspolar_H[store_index] = I_coded_crosspolar_H[obtain_index];
			  IQStruct.Q_coded_crosspolar_H[store_index] = Q_coded_crosspolar_H[obtain_index];		
			  store_index = store_index + 1;
			  obtain_index = obtain_index + 1;
			}
		    }
		  for (i = 0; i < paramUncoded.nfft; i++) 
		    {
		      obtain_index = ( i * paramCoded.num_interleave + 1) * paramUncoded.samples_per_pulse;
		      store_index = (i * paramUncoded.samples_per_pulse) + (nspectra * paramUncoded.samples_per_pulse * paramUncoded.nfft);
		      for (j=0; j<paramUncoded.samples_per_pulse; j++) 
			{	
			  IQStruct.I_uncoded_copolar_H[store_index]    = I_uncoded_copolar_H[obtain_index];
			  IQStruct.Q_uncoded_copolar_H[store_index]    = Q_uncoded_copolar_H[obtain_index];
			  IQStruct.I_uncoded_crosspolar_H[store_index] = I_uncoded_crosspolar_H[obtain_index];
			  IQStruct.Q_uncoded_crosspolar_H[store_index] = Q_uncoded_crosspolar_H[obtain_index];
			  store_index  = store_index + 1;
			  obtain_index = obtain_index + 1;
			}
		    }
		  printf("completed storing IQs\n");
		}
	      
	      
	      
	      // Calculate power spectra for each gate
	      // Loop through gates
	      for (sample=0; sample<param.samples_per_pulse; sample++)
		{
		  register int ii,index;
		  // 1) CODED H-COPOLAR SPECTRUM (HHP)
		  for(ii=0;ii<paramCoded.nfft;ii++)
		    {
		      index=ii*paramCoded.samples_per_pulse+sample;
		      fftw_real_lv (in[ii])=(float)I_coded_copolar_H[index];
		      fftw_imag_lv (in[ii])=(float)Q_coded_copolar_H[index];
		    }
		  RSP_SubtractOffset_FFTW(in,paramCoded.nfft);
		  RSP_CalcPSD_FFTW(in, paramCoded.nfft, p_coded, paramCoded.window, current_PSD, norm_coded);         
		  for(ii=0; ii<paramCoded.npsd; ii++) PSD[sample].HHP[ii]+=current_PSD[ii]/paramCoded.spectra_averaged;
		  
		  // 2) CODED H-CROSSPOLAR SPECTRUM (HVP)
		  for(ii=0;ii<paramCoded.nfft;ii++)
		    {
		      index=ii*paramCoded.samples_per_pulse+sample;
		      fftw_real_lv (in[ii])=(float)I_coded_crosspolar_H[index];
		      fftw_imag_lv (in[ii])=(float)Q_coded_crosspolar_H[index];
		    }
		  RSP_SubtractOffset_FFTW(in,paramCoded.nfft);
		  RSP_CalcPSD_FFTW(in, paramCoded.nfft, p_coded, paramCoded.window, current_PSD, norm_coded);         
		  for(ii=0; ii<paramCoded.npsd; ii++) PSD[sample].HVP[ii]+=current_PSD[ii]/paramCoded.spectra_averaged;
		  
		  // 3) UNCODED H-COPOLAR SPECTRUM (HH)
		  for(ii=0;ii<paramUncoded.nfft;ii++)
		    {
		      index=(ii*paramCoded.num_interleave+1)*paramUncoded.samples_per_pulse+sample;
		      fftw_real_lv (in[ii])=(float)I_uncoded_copolar_H[index];
		      fftw_imag_lv (in[ii])=(float)Q_uncoded_copolar_H[index];
		    }
		  RSP_SubtractOffset_FFTW(in,paramUncoded.nfft);
		  RSP_CalcPSD_FFTW(in, paramUncoded.nfft, p_uncoded, paramUncoded.window, current_PSD, norm_uncoded);         
		  for(ii=0; ii<paramUncoded.npsd; ii++) PSD[sample].HH[ii]+=current_PSD[ii]/paramUncoded.spectra_averaged;
		  
		  // 4) UNCODED H-CROSSPOLAR SPECTRUM (HV)
		  for(ii=0;ii<paramUncoded.nfft;ii++)
		    {
		      index=(ii*paramCoded.num_interleave+1)*paramUncoded.samples_per_pulse+sample;
		      fftw_real_lv (in[ii])=(float)I_uncoded_crosspolar_H[index];
		      fftw_imag_lv (in[ii])=(float)Q_uncoded_crosspolar_H[index];
		    }
		  RSP_SubtractOffset_FFTW(in,paramUncoded.nfft);
		  RSP_CalcPSD_FFTW(in, paramUncoded.nfft, p_uncoded, paramUncoded.window, current_PSD, norm_uncoded);         
		  for(ii=0; ii<paramUncoded.npsd; ii++) PSD[sample].HV[ii]+=current_PSD[ii]/paramUncoded.spectra_averaged;
		  
		  
		}
	    }
	  //---------------------------
	  // END OF SPECTRAL AVERAGING
	  //---------------------------
	  


	  // TX

	  TX1_level = 0;
	  TX1_gate1 = 15;
	  TX1_gate2 = 20;
	  int navg1 = TX1_gate2-TX1_gate1;
	  TX2_level = 0;
	  TX2_gate1 = 13;
	  TX2_gate2 = 18;
	  int navg2 = TX2_gate2-TX2_gate1;
	  
	  TX_counter += 1;

          int np = paramCoded.pulses_per_daq_cycle/paramCoded.num_interleave;

	  printf("np = %5d\n", np);

	  for(j=TX1_gate1; j<TX1_gate2; j++) 
	    {
	      for(i=0; i<paramCoded.pulses_per_daq_cycle; i=i+paramCoded.num_interleave)
		{
		  pos1=param.samples_per_pulse*i+j;
		  TX1_level += (float)TX1_raw[pos1];
		}
	    }
	  TX1_level /= (float)(navg1*np);

	  for(j=TX2_gate1; j<TX2_gate2; j++) 
	    {
	      for(i=0; i<paramCoded.pulses_per_daq_cycle; i=i+paramCoded.num_interleave)
		{
		  pos1=param.samples_per_pulse*i+j;
		  TX2_level += (float)TX2_raw[pos1];
		}
	    }
	  TX2_level /= (float)(navg2*np);

	  //for(j=TX2_gate1; j<TX2_gate2; j++)
	  //  {
	  //    register int ii,index;
	  //    for(ii=0;ii<paramCoded.nfft;ii++)
	  //	{
	  //  index=ii*paramCoded.samples_per_pulse+j;
	  //  TX2_level += (float)TX2_raw[index]/((float)navg2*paramCoded.nfft);
	  //}
	  //}
	  
	  
	  /* update time in spectral information file */	
	  /* get time of day */
	  gettimeofday( &tv, &tz);
	  gmtime_r (&tv.tv_sec, &tm);
	  PSD_obs.year 		= tm.tm_year+1900;
	  PSD_obs.month		= tm.tm_mon+1;
	  PSD_obs.day  		= tm.tm_mday;
	  PSD_obs.hour 		= tm.tm_hour;
	  PSD_obs.minute		= tm.tm_min;
	  PSD_obs.second		= tm.tm_sec;
	  PSD_obs.centisecond	= obs.centisecond=(int)tv.tv_usec/10000;
	  
	  PSD_RAPID_obs.year        = PSD_obs.year;
	  PSD_RAPID_obs.month       = PSD_obs.month;
	  PSD_RAPID_obs.day         = PSD_obs.day;
	  PSD_RAPID_obs.hour        = PSD_obs.hour;
	  PSD_RAPID_obs.minute      = PSD_obs.minute;
	  PSD_RAPID_obs.second      = PSD_obs.second;
	  PSD_RAPID_obs.centisecond = PSD_obs.centisecond;
	  
	  system_time = time(NULL);	
	  if ( collect_spectra_rapid_now == 1 ) 
	    {
	      printf("Writing Rapid PSD Variables... ***************************\n");
	      RNC_WriteRapidLogPSDVariables(spectra_rapid_ncid, COPERNICUS_SPECTRA_RAPID, &param, &PSD_RAPID_obs, PSD, PSD_rapid_varid);
	      status = nc_sync(spectra_ncid);
	      if (status != NC_NOERR) check_netcdf_handle_error(status);
	      spectra_rapid_time =  spectra_rapid_time + param.dump_spectra_rapid;
	      printf("Written Rapid PSD Variables\n");
	    }
	  /* write out spectra to netCDF if required */
	  if ( collect_spectra_now == 1) 
	    {
	      /* reorder the elements in the uncoded (raw) IQ arrays */
	      /* I think what is below is just a carry over from something that was not thought out */
	      /* I think that this can be deleted */
	      counter = 0;
	      for (pulse = 0; pulse < paramUncoded.nfft; pulse++ ) 
		{
		  marker = (pulse * paramCoded.num_interleave + 1) * paramUncoded.samples_per_pulse;
		  for (sample=0; sample<param.samples_per_pulse; sample++) 
		    {
		      I_uncoded_copolar_H[counter]    = I_uncoded_copolar_H[marker + sample];
		      I_uncoded_crosspolar_H[counter] = I_uncoded_crosspolar_H[marker + sample];
		      Q_uncoded_copolar_H[counter]    = Q_uncoded_copolar_H[marker + sample];
		      Q_uncoded_crosspolar_H[counter] = Q_uncoded_crosspolar_H[marker + sample];
		      counter = counter + 1;
		    }
		}
	      /* end of section that was not thought out */
	      printf("Writing PSD Variables... ***************************\n");
	      RNC_WriteLogPSDVariables(spectra_ncid, COPERNICUS_CODED_SPECTRA, &param, &PSD_obs, PSD, &IQStruct, PSD_varid);
	      status = nc_sync(spectra_ncid);
	      if (status != NC_NOERR) check_netcdf_handle_error(status);
	      spectra_time =  spectra_time + param.dump_spectra;
	    }
	  
	  // Calculate noise from upper range gates
	  noisegate1=param.samples_per_pulse-50-(paramCoded.code_length*paramCoded.oversample_ratio);
	  noisegate2=param.samples_per_pulse-1-(paramCoded.code_length*paramCoded.oversample_ratio);
	  count=(noisegate2-noisegate1)+1;
	  HH_noise_level  = 0;
	  HV_noise_level  = 0;
	  HHP_noise_level = 0;
	  HVP_noise_level = 0;
	  for(i=noisegate1; i<=noisegate2; i++)
	    {
	      HH_noise_level  += median(PSD[i].HH,paramUncoded.npsd);
	      HV_noise_level  += median(PSD[i].HV,paramUncoded.npsd);
	      HHP_noise_level += median(PSD[i].HHP,paramCoded.npsd);
	      HVP_noise_level += median(PSD[i].HVP,paramCoded.npsd);
	    }
	  HH_noise_level  /= count;
	  HV_noise_level  /= count;
	  HHP_noise_level /= count;
	  HVP_noise_level /= count;
	  
	  // Now calculate the Doppler parameters
	  printf("** Calculating Doppler parameters...\n");
	  NPC_H[0]  += HH_noise_level; 
	  NPC_HP[0] += HHP_noise_level;
	  NPC_V[0]  += HV_noise_level;
	  NPC_VP[0] += HVP_noise_level;
	  //      printf("Noise = %5.2f\n",NPC_H[0]);
	  
	  // TX

          printf("TX_counter = %5d\n",TX_counter);

	  printf("TX1_level = %5.2f\n",TX1_level);
	  TX_1H[0] += TX1_level/((float)param.moments_averaged);
	  printf("TX2_level = %5.2f\n",TX2_level);
	  TX_2H[0] += TX2_level/((float)param.moments_averaged);
	  
	  printf("TX_1H = %5.2f\n",TX_1H[0]);
	  printf("TX_2H = %5.2f\n",TX_2H[0]);
	  // Loop through all spectra and get parameters
	  for(i=0; i<param.samples_per_pulse; i++)
	    {
	      float noise_power, tempPower, tempVel, tempZED;
	      
	      // interpolate over clutter
	      RSP_ClutterInterp(PSD[i].HH,paramUncoded.npsd,paramUncoded.fft_bins_interpolated);
	      RSP_ClutterInterp(PSD[i].HV,paramUncoded.npsd,paramUncoded.fft_bins_interpolated);
	      RSP_ClutterInterp(PSD[i].HHP,paramCoded.npsd,paramCoded.fft_bins_interpolated);
	      RSP_ClutterInterp(PSD[i].HVP,paramCoded.npsd,paramCoded.fft_bins_interpolated);
	      
	      
	      // Find HH peak
	      RSP_FindPeaksMulti_Destructive(PSD[i].HH,paramUncoded.npsd,paramUncoded.num_peaks,HH_noise_level,HH_peaks);
	      // Calculate HH moments
	      RSP_CalcSpecMom(PSD[i].HH,paramUncoded.npsd,HH_peaks,HH_noise_level,HH_moments, RSP_MOMENTS);
	      // Find HV peak
	      RSP_FindPeaksMulti_Destructive(PSD[i].HV,paramUncoded.npsd,paramUncoded.num_peaks,HV_noise_level,HV_peaks);
	      // Calculate HV moments
	      RSP_CalcSpecMom(PSD[i].HV,paramUncoded.npsd,HV_peaks,HV_noise_level,HV_moments, RSP_MOMENTS);
	      // Find HHP peak
	      RSP_FindPeaksMulti_Destructive(PSD[i].HHP,paramCoded.npsd,paramCoded.num_peaks,HHP_noise_level,HHP_peaks);
	      // Calculate HHP moments
	      RSP_CalcSpecMom(PSD[i].HHP,paramCoded.npsd,HHP_peaks,HHP_noise_level,HHP_moments, RSP_MOMENTS);
	      // Find HVP peak
	      RSP_FindPeaksMulti_Destructive(PSD[i].HVP,paramCoded.npsd,paramCoded.num_peaks,HVP_noise_level,HVP_peaks);
	      // Calculate HVP moments
	      RSP_CalcSpecMom(PSD[i].HVP,paramCoded.npsd,HVP_peaks,HVP_noise_level,HVP_moments, RSP_MOMENTS);
	      
	      
	      
	      // ----------------------------
	      //  PROCESS UNCODED PARAMETERS
	      // ----------------------------
	      gate=i+gate_offset;
	      noise_power = RSP_CalcNoisePower(HH_noise_level,HH_peaks,&paramUncoded);
	      tempPower   = HH_moments[0]*paramUncoded.frequency_bin_width;
	      // printf("tempPower = %d %5.2f\n",i,tempPower);
	      
	      // Calculate weighting coefficient
	      //wi = (peaks[0].peakPSD - HH_noise_level);
	      wi=1; // Turn off weighting
	      //wi = tempPower/noise_power;
	      //if(wi>1) wi=1;
	      uncoded_sum_wi[i] += wi;
	      
	      // COPOLAR
	      SNR_HC[i] += tempPower/noise_power * wi;
	      ZED_HC[i] += tempPower * wi;
	      tempZED    = 10*log10(tempPower);
	      tempVel    = RSP_BinToVelocity(HH_moments[1],&paramUncoded);
	      //VEL_HC[i] += tempVel * wi;
	      VEL_HC_COS[i]       += cos(tempVel/paramUncoded.folding_velocity*PI) * wi;
	      VEL_HC_SIN[i]       += sin(tempVel/paramUncoded.folding_velocity*PI) * wi;
	      SPW_HC[i]           += HH_moments[2]*paramUncoded.frequency_bin_width/paramUncoded.hz_per_mps * wi;
	      SKW_HC[i]           += HH_moments[3] * wi;
	      KRT_HC[i]           += HH_moments[4] * wi;
	      uncoded_mean_vsq[i] += (tempVel*tempVel) * wi;
	      uncoded_mean_Zsq[i] += (tempPower*tempPower) * wi;
	      
	      // CROSSPOLAR
	      noise_power = RSP_CalcNoisePower(HV_noise_level,HV_peaks,&paramUncoded);	      
	      tempPower   = HV_moments[0]*paramUncoded.frequency_bin_width;
	      SNR_XHC[i] += tempPower/noise_power * wi;
	      ZED_XHC[i] += tempPower * wi;	  
	      
	      
	      // ----------------------------
	      //  PROCESS CODED PARAMETERS
	      // ----------------------------
	      if(gate<paramCoded.samples_per_pulse) 
		{
		  noise_power = RSP_CalcNoisePower(HHP_noise_level,HHP_peaks,&paramCoded);
		  tempPower   = HHP_moments[0]*paramCoded.frequency_bin_width;
		  
		  // Calculate weighting coefficient
		  //wi = (peaks[0].peakPSD - HHP_noise_level);
		  wi=1; // Turn off weighting
		  //wi = tempPower/noise_power;
		  //if(wi>1) wi=1;
		  coded_sum_wi[gate] += wi;
		  
		  // COPOLAR
		  SNR_HCP[gate] += tempPower/noise_power * wi;
		  ZED_HCP[gate] += tempPower * wi;
		  tempVel        = RSP_BinToVelocity(HHP_moments[1],&paramCoded);
		  tempZED        = 10*log10(tempPower);
		  //VEL_HCP[gate] += tempVel * wi;
		  VEL_HCP_COS[gate]    += cos(tempVel/paramCoded.folding_velocity*PI) * wi;
		  VEL_HCP_SIN[gate]    += sin(tempVel/paramCoded.folding_velocity*PI) * wi;
		  SPW_HCP[gate]        += HHP_moments[2]*paramCoded.frequency_bin_width/paramCoded.hz_per_mps * wi;
		  SKW_HCP[gate]        += HHP_moments[3] * wi;
		  KRT_HCP[gate]        += HHP_moments[4] * wi;
		  coded_mean_vsq[gate] += (tempVel*tempVel) * wi;
		  coded_mean_Zsq[gate] += (tempPower*tempPower) * wi;
		  
		  // CROSSPOLAR
		  noise_power     = RSP_CalcNoisePower(HVP_noise_level,HVP_peaks,&paramCoded);
		  tempPower       = HVP_moments[0]*paramCoded.frequency_bin_width;
		  SNR_XHCP[gate] += tempPower/noise_power * wi;
		  ZED_XHCP[gate] += tempPower * wi;
		}
	    }
	  
	} // End of moments averaging loop
      
      for(i=0; i<param.samples_per_pulse; i++) 
	{
	  // COMPLETE THE WEIGHTED AVERAGING WITH DIVISION
	  SNR_HC[i]      /= uncoded_sum_wi[i];
	  ZED_HC[i]      /= uncoded_sum_wi[i];
	  //VEL_HC[i] /= uncoded_sum_wi[i];
	  VEL_HC_COS[i]  /= uncoded_sum_wi[i];
	  VEL_HC_SIN[i]  /= uncoded_sum_wi[i];
	  SPW_HC[i]      /= uncoded_sum_wi[i];
	  SKW_HC[i]      /= uncoded_sum_wi[i];
	  KRT_HC[i]      /= uncoded_sum_wi[i];
	  SNR_XHC[i]     /= uncoded_sum_wi[i];
	  ZED_XHC[i]     /= uncoded_sum_wi[i];
	  SNR_HCP[i]     /= coded_sum_wi[i];
	  ZED_HCP[i]     /= coded_sum_wi[i];
	  //VEL_HCP[i] /= coded_sum_wi[i];
	  VEL_HCP_COS[i] /= coded_sum_wi[i];
	  VEL_HCP_SIN[i] /= coded_sum_wi[i];
	  SPW_HCP[i]     /= coded_sum_wi[i];
	  SKW_HCP[i]     /= coded_sum_wi[i];
	  KRT_HCP[i]     /= coded_sum_wi[i];
	  SNR_XHCP[i]    /= coded_sum_wi[i];
	  ZED_XHCP[i]    /= coded_sum_wi[i];
	  
	  VEL_HC[i]  = atan2(VEL_HC_SIN[i],VEL_HC_COS[i])/PI*paramUncoded.folding_velocity;
	  VEL_HCP[i] = atan2(VEL_HCP_SIN[i],VEL_HCP_COS[i])/PI*paramCoded.folding_velocity;
	  
	  // Calculate sigma Z bar (this is the relative linear standard deviation)
	  uncoded_mean_Zsq[i] /= uncoded_sum_wi[i];
	  ZED_HCD[i]           = sqrt( uncoded_mean_Zsq[i] - (ZED_HC[i]*ZED_HC[i]) );
	  coded_mean_Zsq[i]   /= coded_sum_wi[i];
	  ZED_HCDP[i]          = sqrt( coded_mean_Zsq[i] - (ZED_HCP[i]*ZED_HCP[i]) );
	  ZED_HCD[i]           = 10*log10(ZED_HCD[i]/ZED_HC[i]);
	  ZED_HCDP[i]          = 10*log10(ZED_HCDP[i]/ZED_HCP[i]);
	  
	  // Convert SNRs and ZEDs to dB
	  SNR_HC[i]   = 10*log10(SNR_HC[i]);
	  ZED_HC[i]   = 10*log10(ZED_HC[i]);
	  SNR_HCP[i]  = 10*log10(SNR_HCP[i]);
	  ZED_HCP[i]  = 10*log10(ZED_HCP[i]);
	  SNR_XHC[i]  = 10*log10(SNR_XHC[i]);
	  ZED_XHC[i]  = 10*log10(ZED_XHC[i]);
	  SNR_XHCP[i] = 10*log10(SNR_XHCP[i]);
	  ZED_XHCP[i] = 10*log10(ZED_XHCP[i]);
	  
	  // Calculate sigma v bar
	  uncoded_mean_vsq[i] /= uncoded_sum_wi[i];
	  VEL_HCD[i]           = sqrt( uncoded_mean_vsq[i] - (VEL_HC[i]*VEL_HC[i]) );
	  coded_mean_vsq[i]   /= coded_sum_wi[i];
	  VEL_HCDP[i]          = sqrt( coded_mean_vsq[i] - (VEL_HCP[i]*VEL_HCP[i]) );
	  
	  // Calculate LDR
	  LDR_C[i]  += ZED_XHC[i]-ZED_HC[i]+paramUncoded.LDR_calibration_offset;
	  LDR_CP[i] += ZED_XHCP[i]-ZED_HCP[i]+paramCoded.LDR_calibration_offset;
	  
	  // Do range correction
	  ZED_H[i]    += 10*log10(paramUncoded.range[i]*paramUncoded.range[i])+paramUncoded.ZED_incoherent_calibration_offset;
	  ZED_HC[i]   += 10*log10(paramUncoded.range[i]*paramUncoded.range[i])+paramUncoded.ZED_calibration_offset;
	  ZED_XHC[i]  += 10*log10(paramUncoded.range[i]*paramUncoded.range[i])+paramUncoded.ZED_calibration_offset;
	  ZED_HCP[i]  += 10*log10(paramCoded.range[i]*paramCoded.range[i])+paramCoded.ZED_calibration_offset;
	  ZED_XHCP[i] += 10*log10(paramCoded.range[i]*paramCoded.range[i])+paramCoded.ZED_calibration_offset;
	  
	}
      
      NPC_H[0]  /= uncoded_sum_wi[0];
      NPC_V[0]  /= uncoded_sum_wi[0];
      NPC_HP[0] /= coded_sum_wi[gate_offset];
      NPC_VP[0] /= coded_sum_wi[gate_offset];
      
      NPC_H[0]  = 10*log10(NPC_H[0]);
      NPC_V[0]  = 10*log10(NPC_V[0]);
      NPC_HP[0] = 10*log10(NPC_HP[0]);
      NPC_VP[0] = 10*log10(NPC_VP[0]);
      
      // TX
      TX_2H[0] *= 3000.0/4096.0;
      TX_1H[0] *= 3000.0/4096.0;

      // write out variables to netCDF 
      printf("Writing dynamic variables to NetCDF...\n");
      RNC_WriteDynamicVariables(ncid, &param, &obs);
      status = nc_sync(ncid);
      if (status != NC_NOERR) check_netcdf_handle_error(status);

      /* check to see if we have started a new day */
      /* the + 2L covers us for the time it takes to wrap around this loop */
      system_time = time(NULL) + 2L;
      gmtime_r (&system_time, &tm);
      if (tm.tm_mday != obs.day) 
	{
	  exit_now = 1;
	  printf("we are at the end of the day, exiting now\n");
	}
  
    }
  
  // Finish off
  printf("*** Closing PCICARD...\n");
  RDQ_ClosePCICARD_New( amcc_fd, &dma_buffer, DMA_BUFFER_SIZE );
  
    
  /* netCDF : close the netCDF file */
  status = nc_sync(ncid);
  if (status != NC_NOERR) check_netcdf_handle_error(status);
  status = nc_close(ncid);
  if (status != NC_NOERR) check_netcdf_handle_error(status);

  if ( param.dump_spectra  != 0) 
    {
      status = nc_sync(spectra_ncid);
      if (status != NC_NOERR) check_netcdf_handle_error(status);
      status = nc_close(spectra_ncid);
      if (status != NC_NOERR) check_netcdf_handle_error(status);
    }
  
  if ( param.dump_spectra_rapid  != 0)
    {
      status = nc_sync(spectra_rapid_ncid);
      if (status != NC_NOERR) check_netcdf_handle_error(status);
      status = nc_close(spectra_rapid_ncid);
      if (status != NC_NOERR) check_netcdf_handle_error(status);
    }
  
  //---------------------------
  // Unallocate all the memory
  //---------------------------
  RSP_FreeMemory(&paramCoded);  // Free memory allocated by RSP package
  RSP_FreeMemory(&paramUncoded);  // Free memory allocated by RSP package
  RSP_ObsFree(&obs); // Free observables memory
  free(timeseries);
  free(current_PSD);
  for(i=0; i<param.samples_per_pulse; i++)
    {
      free(PSD[i].HH);
      free(PSD[i].HV);
      free(PSD[i].VV);
      free(PSD[i].VH);
    }
  free (PSD);
  free(uncoded_mean_vsq);
  free(coded_mean_vsq);
  free(uncoded_mean_Zsq);
  free(coded_mean_Zsq);
  free(uncoded_sum_wi);
  free(coded_sum_wi);

  free(VEL_HC_COS);
  free(VEL_HC_SIN);
  free(VEL_HCP_COS);
  free(VEL_HCP_SIN);

  free(IQStruct.I_coded_copolar_H);
  free(IQStruct.Q_coded_copolar_H);
  free(IQStruct.I_coded_crosspolar_H);
  free(IQStruct.Q_coded_crosspolar_H);
  free(IQStruct.I_uncoded_copolar_H);
  free(IQStruct.Q_uncoded_copolar_H);
  free(IQStruct.I_uncoded_crosspolar_H);
  free(IQStruct.Q_uncoded_crosspolar_H);

  fftw_destroy_plan(p_coded);
  fftw_destroy_plan(p_uncoded);
  fftw_free(in);

  //=========
  // THE END
  //=========
  printf("All done.\n");
  exit(0);
}
