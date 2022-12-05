#ifndef __RADAR_H
#define __RADAR_H

#include <stdint.h>

#define RADAR_DATA_PATH   "/radar/" /* "/mnt/focus_radar_data/" */

#define MAX_OPERATOR_NAME 32  /* MTF: Updated from 4 to 32 characters */

/*
 * Define radar ID numbers
 * (inherited from Andrew's old CF format)
 * CJW renamed HIGHLANDER to VIPER - 20080716
 */
enum RADAR_en
{
    CAMRA                       = 1,
    RABELAIS,                  /* 2  */
    GALILEO,                   /* 3  */
    BISTATIC,                  /* 4  */
    VIPER,                     /* 5  */
    SINGAPORE,                 /* 6  */
    COPERNICUS,                /* 7  */
    ACROBAT,                   /* 8  */
    TEST,                      /* 9  */ /* Galileo */
    COPERNICUS_SPECTRA          = 71,
    COPERNICUS_CODED_SPECTRA,  /* 72 */
    ACROBAT_CODED_SPECTRA,     /* 73 */
    GALILEO_SPECTRA,           /* 74 */
    GALILEO_SPECTRA_RAPID,     /* 75 */
    COPERNICUS_SPECTRA_RAPID,  /* 76 */
    TEST_SPECTRA                = 91,   /* Galileo */
    TEST_SPECTRA_RAPID         /* 92 */ /* Galileo */
};

/*
 * Define scan types
 * CJW added scan type SCAN_MAN - 20101020
 */
enum ScanType_en
{
    SCAN_PPI   = 0,   /* PPI Scan */
    SCAN_RHI, /* 1 */ /* RHI Scan */
    SCAN_FIX, /* 2 */ /* Fixed dwell */
    SCAN_CSP, /* 3 */ /* CSP scan */
    SCAN_CAL, /* 4 */ /* Calibration scan */
    SCAN_r5,  /* 5 */
    SCAN_MAN, /* 6 */ /* Manual scan (also for target tracking) */
    SCAN_TRK, /* 7 */ /* Gandalf tracking scan */
    SCAN_r8,  /* 8 */
    SCAN_SGL  /* 9 */ /* Single ray */
};

/*
 * URC_ScanStruct : Holds information about the scan
 * (mostly from the information record)
 * MTF: 'date' field size increased from 14 to 15 to fix buffer overrun
 */
typedef struct
{
    enum ScanType_en scanType;
    float            min_angle;
    float            max_angle;
    float            scan_angle;
    float            scan_velocity;
    float            dwelltime;
    int              file_number; /* Obsolete */
    int              scan_number;
    int              experiment_id;
    char             operator [MAX_OPERATOR_NAME];
    char             date     [15];
    int              min_gate;
    int              max_gate;
} URC_ScanStruct;

/*
 * File name creation helpers:
 *  <radar-data-path>/<radar-name>/<spectra-name>/<date-ccyymmdd>/
 *  <radar-name><host-ext>_<time-ccyymmddhhmmss>_<scan-type>
 *  <specta-ext>-raw.<ext>
 *
 * <host-ext>: "-felixx", "-bolt" or ""
 */

static inline const char *
GetRadarName (enum RADAR_en radar)
{
    switch (radar)
    {
    case CAMRA:
	return "radar-camra";

    case ACROBAT_CODED_SPECTRA:
    case ACROBAT:
	return "radar-acrobat";

    case COPERNICUS_SPECTRA:
    case COPERNICUS_CODED_SPECTRA:
    case COPERNICUS_SPECTRA_RAPID:
    case COPERNICUS:
	return "radar-copernicus";

    case GALILEO_SPECTRA:
    case GALILEO_SPECTRA_RAPID:
    case GALILEO:
	return "radar-galileo";

    case TEST_SPECTRA:
    case TEST_SPECTRA_RAPID:
    case TEST:
	return "radar-test";

    default:
	break;
    }
    return "radar-new";
}

static inline const char *
GetSpectraName (enum RADAR_en radar)
{
    switch (radar)
    {
    case COPERNICUS_SPECTRA_RAPID:
    case GALILEO_SPECTRA_RAPID:
    case TEST_SPECTRA_RAPID:
	return "spectra-rapid";

    case ACROBAT_CODED_SPECTRA:
    case COPERNICUS_CODED_SPECTRA:
    case COPERNICUS_SPECTRA:
    case GALILEO_SPECTRA:
    case TEST_SPECTRA:
	return "spectra";

    case ACROBAT:
    case COPERNICUS:
    case GALILEO:
    case CAMRA:
    case TEST:
    default:
	break;
    }
    return "raw";
}

static inline const char *
GetScanTypeName (enum ScanType_en scan_type)
{
    switch (scan_type)
    {
    case SCAN_PPI: return "ppi";
    case SCAN_RHI: return "rhi";
    case SCAN_CSP: return "s-p";
    case SCAN_FIX: return "fix";
    case SCAN_SGL: return "ray";
    case SCAN_CAL: return "cal";
    case SCAN_MAN: return "man";
    case SCAN_TRK: return "trk";
    default:       break;
    }
    return "tst";
}

static inline const char *
GetSpectraExtension (enum RADAR_en radar)
{
    switch (radar)
    {
    case ACROBAT_CODED_SPECTRA :
    case COPERNICUS_CODED_SPECTRA :
    case COPERNICUS_SPECTRA :
    case GALILEO_SPECTRA :
    case TEST_SPECTRA :
	return "-fft";

    case COPERNICUS_SPECTRA_RAPID :
    case GALILEO_SPECTRA_RAPID :
    case TEST_SPECTRA_RAPID :
	return "-fft-rapid";

    default:
	break;
    }
    return "";
}

// This structure holds pointers to the various PSDs
typedef struct
{
    float *HH;  // H copolar
    float *HV;  // H crosspolar
    float *VV;  // V copolar
    float *VH;  // V crosspolar
    float *HHP; // H copolar (pulse coded)
    float *HVP; // H crosspolar (pulse coded)
    float *VVP; // V copolar (pulse coded)
    float *VHP; // V crosspolar (pulse coded)
} PolPSDStruct;

/* this will hold the IQ data for a spectra */
typedef struct
{
    long int * I_coded_copolar_H;
    long int * Q_coded_copolar_H;
    long int * I_coded_copolar_V;
    long int * Q_coded_copolar_V;
    long int * I_coded_crosspolar_H;
    long int * Q_coded_crosspolar_H;
    long int * I_coded_crosspolar_V;
    long int * Q_coded_crosspolar_V;
    uint16_t * I_uncoded_copolar_H;
    uint16_t * Q_uncoded_copolar_H;
    uint16_t * I_uncoded_copolar_V;
    uint16_t * Q_uncoded_copolar_V;
    uint16_t * I_uncoded_crosspolar_H;
    uint16_t * Q_uncoded_crosspolar_H;
    uint16_t * I_uncoded_crosspolar_V;
    uint16_t * Q_uncoded_crosspolar_V;
} IQStruct;


// This structure is used to address the various pulses in the
// transmit sequence
typedef struct
{
    int HP;  // H pulse coded
    int VP;  // V pulse coded
    int H;   // H uncoded
    int V;   // V uncoded
} URC_SampleStruct;

#endif
