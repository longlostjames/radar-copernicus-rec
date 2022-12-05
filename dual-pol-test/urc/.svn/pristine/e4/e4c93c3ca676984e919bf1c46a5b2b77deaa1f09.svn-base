#ifndef _RNC_H
#define _RNC_H

#include <RSP.h>
#include <radar.h>

/*
 * This will stop any screws up with the PSD_varid.
 * This only works if you allocate PSD_varid to be largenough.
 * If re-ordered, can improve things for Galileo a bit -
 * not worth the efort yet.
 */
enum PSD_varidx_en
{
    PSD_HH    = 0,
    PSD_HV,  /* 1  */
    PSD_VV,  /* 2  */
    PSD_VH,  /* 3  */
    /*             */
    PSD_HHP, /* 4  */
    PSD_HVP, /* 5  */
    PSD_VVP, /* 6  */
    PSD_VHP, /* 7  */
    /*             */
    IPF_HH,  /* 8  */
    QPF_HH,  /* 9  */
    IPF_HV,  /* 10 */
    QPF_HV,  /* 11 */
    IPF_VV,  /* 12 */
    QPF_VV,  /* 13 */
    IPF_VH,  /* 14 */
    QPF_VH,  /* 15 */
    /*             */
    IPF_HHP, /* 16 */
    QPF_HHP, /* 17 */
    IPF_HVP, /* 18 */
    QPF_HVP, /* 19 */
    IPF_VVP, /* 20 */
    QPF_VVP, /* 21 */
    IPF_VHP, /* 22 */
    QPF_VHP, /* 23 */
    PSD_varidSize /* Must always be last */
};

#define RapidPSD_varidSize   2

typedef struct
{
    int time_dim;
    int range_dim;
    int unaveraged_range_dim;
    int peak_dim;
    int time_range_dim;
    int fft_bin_dim;
    int coded_fft_bin_dim;
    int spectra_number_dim;
} RNC_DimensionStruct;

extern int  RNC_OpenNetcdfFile (const char * radar_name,
				const char * sectra_name,
				const char * date,
				const char * host_ext,
				const char * scan_name,
				const char * sectra_ext,
				const char * rectype);
extern void   RNC_SetupGlobalAttributes (int ncid, int radar,
					 const URC_ScanStruct * scan,
					 const RSP_ParamStruct * param,
					 int argc, char * const argv[]);
extern void   RNC_SetupPulse_Compression_Code (int ncid,
					       const RSP_ParamStruct * param);
extern int    RNC_SetupFile_State (int ncid);
extern void   RNC_SetupStaticVariables (int ncid,
					const RSP_ParamStruct * param);
extern void   RNC_SetupDimensions (int ncid, const RSP_ParamStruct * param,
				   RNC_DimensionStruct * dimensions);
extern void   RNC_SetupRapidLogPSDDimensions (int ncid, int radar,
					      const RSP_ParamStruct * param,
					      RNC_DimensionStruct * dimensions);
extern void   RNC_SetupRange (int ncid, const RSP_ParamStruct * param,
			      const RNC_DimensionStruct * dimensions);
extern void   RNC_SetupDynamicVariables (int ncid, int radar,
					 const URC_ScanStruct * scan,
					 const RSP_ParamStruct * param,
					 const RNC_DimensionStruct * dimensions,
					 RSP_ObservablesStruct * obs);
extern void   RNC_SetupPSDVariables (int ncid, int radar,
				     const RSP_ParamStruct * param,
				     const RNC_DimensionStruct * dimensions,
				     int * PSD_varid);
extern void   RNC_SetupLogPSDVariables (int ncid, int radar,
					const RSP_ParamStruct * param,
					const RNC_DimensionStruct * dimensions,
					int * PSD_varid);
extern void   RNC_WriteDynamicVariables (int ncid,
					 const RSP_ParamStruct * param,
					 RSP_ObservablesStruct * obs);
extern void   RNC_WritePSDVariables (int ncid, int radar,
				     const RSP_ParamStruct * param,
				     RSP_ObservablesStruct * obs,
				     const PolPSDStruct * PSD,
				     const int * PSD_varid);
extern void   RNC_WriteLogPSDVariables (int ncid, int radar,
					const RSP_ParamStruct * param,
					RSP_ObservablesStruct * obs,
					const PolPSDStruct * PSD,
					const IQStruct * IQStruct,
					const int * PSD_varid);
extern void   RNC_WriteRapidLogPSDVariables (int ncid, int radar,
					     const RSP_ParamStruct * param,
					     RSP_ObservablesStruct * obs,
					     const PolPSDStruct * PSD,
					     const int * PSD_varid);
extern void   check_netcdf_handle_error (int status);
extern double RNC_GetConfigDouble (const char * filename, const char * keyword);
extern float  RNC_GetConfigFloat  (const char * filename, const char * keyword);
extern int    RNC_GetConfig       (const char * filename, const char * keyword,
				   char * value, size_t value_size);

#endif /* !_RNC_H */
