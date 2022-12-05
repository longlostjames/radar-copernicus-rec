#ifndef __RSP_H
#define __RSP_H

// RSP.h
// General header file for the Radar Signal Processing Package.
// Contains definitions required by more than one module.

#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
//#define FFTW_NO_Complex
#include <fftw3.h>

#if !defined(FFTW_NO_Complex) && defined(_Complex_I) && defined(complex) && defined(I)
#define fftw_real(x)        creal (x)
#define fftw_imag(x)        cimag (x)
#define fftw_real_lv(x)     __extension__ __real__ (x)
#define fftw_imag_lv(x)     __extension__ __imag__ (x)
#define fftw_cassign(x,r,i) (x) = (r) + _Complex_I * (i)
#define fftw_is_complex 1
#else
#define fftw_real(x)        (x)[0]
#define fftw_imag(x)        (x)[1]
#define fftw_real_lv(x)     (x)[0]
#define fftw_imag_lv(x)     (x)[1]
#define fftw_cassign(x,r,i) (x)[0] = (r); (x)[1] = (i)
#undef  fftw_is_complex
#endif

// Define useful constants
#define SPEED_LIGHT 299792458
#ifdef M_PI
#define PI      M_PI
#define TWOPI   (2.0 * M_PI)
#define RAD2DEG (180.0 / M_PI)
#else
#define PI      3.141592654
#define TWOPI   6.283185307
#define RAD2DEG 57.29577951
#endif

// Blackman Window
#define BLACKMAN_WINDOW(j,N) (0.42 - 0.5 * cos (TWOPI * (j)/(N)) + 0.08 * cos (2.0 * TWOPI * (j)/(N)))

#define MAX_OBSERVABLES 100
#define MAX_NAME_LENGTH 15
#define RSP_MAX_MOMENTS 5
#define RSP_MIN_MOMENTS 3

// The RSP_PeakStruct type definition defines the structure used to store
// information about spectral peaks.
typedef struct
{
    int   peakBin;  // Bin no. of spectral peak
    int   leftBin;  // Bin no. of left edge of peak
    int   rightBin; // Bin no. of right edge of peak
    float peakPSD;  // Peak linear PSD value
} RSP_PeakStruct;

// The RSP_ParamStruct type definition defines the structure used to store
// basic operating parameters of the radar.
// Parameters marked + must be set by the main program; others are calculated by
// the RSP_InitialiseParams module.
typedef struct
{
    float   frequency;                                    // + Radar frequency (GHz)
    float   prf;                                          // + Pulse repetition frequency (Hz)
    float   prt;                                          //   Pulse repetition time (s)
    float   pulse_period;                                 // + Duration of TX pulse or subpulse (nanosec)
    int     pulses_per_daq_cycle;                         // + Number of pulses collected by ADC per ray
    int     samples_per_pulse;                            // + Number of samples collected per pulse
    int     clock_divfactor;                              // + The sample clock divide factor
    int     delay_clocks;                                 // + The number of ADC delay clock cycles
    int     ADC_channels;                                 // + Number of ADC channels to use
    int     pulses_coherently_averaged;                   // + Number of pulses coherently averaged
    int     spectra_averaged;                             // + Number of spectra averaged
    int     moments_averaged;                             // + Number of moments to average (if required)
    int     nfft;                                         //   Number of points per fft
    int     npsd;                                         //   Number of bins in PSD
    int     code_length;                                  // + Number of bits per pulse code (1 = no coding)
    int     number_of_codes;                              // + Number of complementary codes (1 = no complementary codes)
    float   frequency_bin_width;                          //   Width of FFT frequency bins (Hz)
    float   hz_per_mps;                                   //   Doppler frequency per m/s (Hz)
    float * frequency_axis;                               //   PSD frequency axis (Hz)
    float * velocity_axis;                                //   PSD velocity axis (m/s)
    float   folding_velocity;                             //   Nyquist folding velocity (m/s)
    float   folding_frequency;                            //   Nyquist folding frequency (Hz)
    float   range_gate_width;                             //   Width of range gates (m)
    float   range_resolution;                             //   Range resolution (m)
    float * range;                                        //   Distance to centre of each range gate (m)
    float   range_offset;                                 // + Distance to centre of first range gate (m)
    float   dwell_time;                                   //   Total length of dwell, inc. all averaging (s)
    float   daq_time;                                     //   Time for one DAQ cycle (s)
    float   clock;                                        // + Clock input to DAQ cards (Hz)
    float   sample_frequency;                             //   Sample frequency (Hz) after clock division
    float   sample_period;                                //   Sample period (s) after clock division
    float * window;                                       //   FFT window of length nfft
    int     long_pulse_mode;                              //   Flag to indicate long pulse mode
    int     alternate_modes;                              //   Flag to indicate if will alternate between mode0 and mode1
    int     mode0;                                        //   Wivern mode0
    int     mode1;                                        //   Wivern mode1
    int     nrays_mode0;                                  //   Number of rays in Wivern mode0
    int     nrays_mode1;                                  //   Number of rays in Wivern mode1
    float   Wss;                                          //   Window Squared and Summed (for fft)
    int     fft_bins_interpolated;                        // + Number of bins either side of 0 to interpolate over
    int     oversample_ratio;                             //   Ratio of sampling period to pulse (chip) duration
    int     include_iq_in_spectra;                        //   Flag to indicate if IQ data is stored within the spectra file
    short   codes [32][32];                               // + Pulse code data
    float   ZED_calibration_offset;                       // + Reflectivity calibration offset
    float   ZDR_calibration_offset;                       // + Differential reflectivity calibration offset
    float   LDR_calibration_offset;                       // + Linear depolarisation ratio calibration offset
    float   ZED_incoherent_calibration_offset;            // + Incoherent reflectivity calibration offset
    float   ZED_crosspolar_incoherent_calibration_offset; // + Crosspolar incoherent reflectivity calibration offset
    float   ZED_incoherent_noise;                         // + Incoherent channel mean noise level
    float   pulse_offset;                                 // + Delay between HV pulse pairs (microsec)
    float   phidp_offset;                                 // + Phidp offset (degrees)

    /* the calibration offsets are now going to be much more belt and braces */
    /* I have not deleted the previous method of defining calibration offsets, so as to provide a means of backward */
    /* compatability */
    /* ZED */
    float    ZED_HC_calibration_offset;    // ZED_HC
    float    ZED_XHC_calibration_offset;   // ZED_XHC
    float    ZED_VC_calibration_offset;    // ZED_VC
    float    ZED_XVC_calibration_offset;   // ZED_XVC
    float    ZED_HCP_calibration_offset;   // ZED_HCP
    float    ZED_XHCP_calibration_offset;  // ZED_XHCP
    float    ZED_VCP_calibration_offset;   // ZED_VCP
    float    ZED_XVCP_calibration_offset;  // ZED_XVCP
    /* ZDR */
    float    ZDR_C_calibration_offset;     // ZDR_C
    float    ZDR_CP_calibration_offset;    // ZDR_CP
    /* LDR */
    float    LDR_HC_calibration_offset;    // LDR_HC
    float    LDR_HCP_calibration_offset;   // LDR_HCP
    float    LDR_VC_calibration_offset;    // LDR_VC
    float    LDR_VCP_calibration_offset;   // LDR_VCP

    float *  range_correction;             //   Range correction factor
    float    azimuth_offset;               // + Beam azimuth offset relative to prime focus beam
    char     code_name [255];              //   The name of the code
    float    transmit_power;               // + TX power (Watts)
    int      num_peaks;                    // + Number of Doppler peaks to process
    float    antenna_diameter;             // + The diameter of the antenna (m)
    float    beamwidthH;                   // + The horizontal beamwidth of the antenna
    float    beamwidthV;                   // + The vertical beamwidth of the antenna
    float    height;                       // + The height of the antenna above OSGB mean sea level
    int      dump_spectra;                 // + allows for the dump of spectra every dump_spectra seconds(if 0 then it does not dump spectra)
    int      dump_spectra_rapid;           // + allows for the rapid dump of spectra in a similar manner to dump_spectra
    int      dump_spectra_every;           // + allows radar-acrobat to dump spectra every so many rays
    int      num_interleave;               // + number of pulse modes interleaved
    int      num_tx_pol;                   // + number of TX polarisations used
    float    mod_pulse_length;             //   the length of the mod pulse (for copernicus only);
    int      real_time_spectra_display;    // denotes if real time spectra is to be displayed
} RSP_ParamStruct;


// The RSP_ObservablesStruct structure is used to store the measured parameters
// for each ray. The pointers and tags should be defined using the
// RSP_ObsNew function.
typedef struct
{
    float   azimuth;
    int     azimuthid;
    float   elevation;
    int     elevationid;
    int     bin_numberid;
    int     ray_numberid;
    /* FIXME: Make a time structure [or use struct timespec] */
    int     year;
    int     month;
    int     day;
    int     hour;
    int     minute;
    int     second;
    int     centisecond;
    /* */
    int     dish_year;
    int     dish_month;
    int     dish_day;
    int     dish_hour;
    int     dish_minute;
    int     dish_second;
    int     dish_centisecond;
    /* */
    int     ts;
    int     tsid;
    int     dish_tsid;
    int     ray_number;
    int     PSD_ray_number;
    int     bin_ray_number;
    int     n_obs;
    /* FIXME: Make a structure */
    int     n_elements        [MAX_OBSERVABLES];
    int     varid             [MAX_OBSERVABLES];
    char    name              [MAX_OBSERVABLES][MAX_NAME_LENGTH];
    float * data              [MAX_OBSERVABLES];
    int     record_observable [MAX_OBSERVABLES];
} RSP_ObservablesStruct;


// Define a complex number type (tried C99 complex.h method, but wouldn't compile)
typedef struct
{
    float real; // Real component of complex number
    float imag; // Imaginary component of complex number
} RSP_ComplexType;

extern void    RSP_InitialiseParams (RSP_ParamStruct * param);
extern void    RSP_FreeMemory (RSP_ParamStruct * param);

extern void    RSP_FindPeaks (float * psd, int nBins, int nPeaks, float noiseLevel, RSP_PeakStruct * peaks);
extern int     RSP_FindPeaksMulti (const float * psd, int nBins, int nPeaks, float noiseLevel, RSP_PeakStruct * peaks);
extern void    RSP_FindPeaksMulti_Destructive (float * psd, int nBins, int nPeaks, float noiseLevel, RSP_PeakStruct * peaks);
extern void    RSP_FindEdges (const float * psd, int nBins, float noiseLevel, RSP_PeakStruct * peak);

extern int     RSP_CalcSpecMom (const float * psd, int nBins, const RSP_PeakStruct * peak, float noiseLevel, float * moments, size_t num_moments);
extern float   RSP_BinToVelocity (float bin, const RSP_ParamStruct * param);
extern float   RSP_CalcNoisePower (float noiseLevel, const RSP_PeakStruct * peak, const RSP_ParamStruct * param);

extern void    RSP_CalcPSD (RSP_ComplexType * IQ, int nfft, const float * window, float * psd, float norm);
extern void    RSP_CalcPSD_FFTW (fftw_complex * in, int nfft, const fftw_plan p, const float * window, float * psd, float norm);
extern void    RSP_SubtractOffset (RSP_ComplexType * IQ, int nfft);
extern void    RSP_SubtractOffset_FFTW (fftw_complex * IQ, int nfft);
extern void    RSP_FFT (float * data, unsigned long nn, int isign);
extern void    RSP_FFT2PowerSpec (const float * data, float * PSD, int nfft, float norm);
extern void    RSP_FFT2PowerSpec_FFTW (fftw_complex * data, float * PSD, int nfft, float norm);

extern void    RSP_Correlate (const uint16_t * data, const short * code, int samples, int bits, long int * corr);
extern void    RSP_Oversample (const short * code, short * newcode, int numel, int n);

extern void    RSP_ClutterInterp (float * PSD, int nBins, int nInterp);

extern void    RSP_CalcPhase (const RSP_ComplexType * IQ, float * phi, float * sdphi, int nfft);

extern void    RSP_ObsInit (RSP_ObservablesStruct * obs);
extern float * RSP_ObsNew (RSP_ObservablesStruct * obs, const char * name, int n_elements, int record_observable);
extern void    RSP_ObsFree (RSP_ObservablesStruct * obs);
extern float * RSP_ObsGet (const RSP_ObservablesStruct * obs, const char * name);
extern int     RSP_ObsIndex (const RSP_ObservablesStruct * obs, const char * name);

extern void    RSP_DisplayParams (const RSP_ParamStruct * param);

#endif /* !__RSP_H */
