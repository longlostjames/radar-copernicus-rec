/*
  Purpose:  To oversee the generate of netCDF files

  Created on:  12/02/2003
  Created by:  Owain Davies (OTD)
  Modified on: 29/08/2007
  Modified by: Chris Walden (CJW)
  Modified on: 13/09/2007
  Modified by: Owain Davies
  added in the ability to dump out IQ data prior to fft
  Modified on: 29/09/2010
  Modified by: John Nicol (JCN)
  included dual-pol. parameters; phi dp and rho hv
  Modified on 04/10/2010 (JCN)
  included uncoded rapid spectra dump
  Modified on 28/11/2011 (JCN)
  included skewness and kurtosis
  Modified on 08/04/2014 (JCN)
  fixed bug in Galileo rapid spectra dump
  Modified on 13/11/2014 (JCN)
  included noise power counts
  Modified on 02/06/2015 (JCN)
  fixed bug in VV spectra
  Modified on 03/08/20115 (JCN)
  include interpolated power correlation
*/

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <alloca.h>

#include <RNC.h>
#include <RSP.h>

/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
int
RNC_OpenNetcdfFile (const char * radar_name,
		    const char * sectra_name,
		    const char * date,
		    const char * host_ext,
		    const char * scan_name,
		    const char * spectra_ext,
		    const char * rectype)
{
    /*--------------------------------------------------------------------------*
     * IN:                                                                      *
     * OUT: ncid : the netCDF file id                                           *
     * RETURN:                                                                  *
     *--------------------------------------------------------------------------*/
    char * netcdf_pathfile;
    char * pt;
    int	   status;
    int    ncid;
    mode_t mask = 002;
    size_t size;

    umask (mask);

    size = strlen (RADAR_DATA_PATH) + 21;
    size += strlen (radar_name) << 1;
    size += strlen (sectra_name);
    if (host_ext != NULL)
	size += strlen (host_ext);
    size += strlen (date);
    size += strlen (scan_name);
    if (spectra_ext != NULL)
	size += strlen (spectra_ext);
    netcdf_pathfile = alloca (size);

    pt    = stpcpy (netcdf_pathfile, RADAR_DATA_PATH);
    pt    = stpcpy (pt, radar_name);
    *pt++ = '/';
    pt    = stpcpy (pt, sectra_name);
    *pt++ = '/';
    stpncpy (pt, date, 8);
    pt[8] = '\0';
    pt   += strlen (pt);

    if (!mkdir (netcdf_pathfile, (S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)))
    {
	printf ("directory created : %s\n", netcdf_pathfile);
    }
    else
    {
	if (errno != EEXIST)
	{
	    printf ("mkdir ERROR, errno : %d (%m), netcdf_pathfile : %s\n",
		    errno,
		    netcdf_pathfile);
	}
    }

    *pt++ = '/';
    pt    = stpcpy (pt, radar_name);
    if (host_ext != NULL)
	pt = stpcpy (pt, host_ext);
    *pt++ = '_';
    pt    = stpcpy (pt, date);
    *pt++ = '_';
    pt    = stpcpy (pt, scan_name);
    if (spectra_ext != NULL)
	pt = stpcpy (pt, spectra_ext);
    *pt++ = '-';
    pt    = stpcpy (pt, rectype);
    strcpy (pt, ".nc");

    printf ("netCDF creating : %s\n", netcdf_pathfile);

    /* NC_64BIT_OFFSET not available on all systems */
    status = nc_create (netcdf_pathfile,
			NC_NOCLOBBER | NC_SHARE,
			&ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    return ncid;
}

/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
void
RNC_SetupGlobalAttributes (int                     ncid,
			   int                     radar,
			   const URC_ScanStruct *  scan,
			   const RSP_ParamStruct * param,
			   int                     argc,
			   char * const            argv [])
{

    char         buffer[1024];
    int          temp;
    float        temp_float;
    int          status;
    int          n;
    time_t       now;
    struct tm    tm;
    char       * pt;
    const char * scan_type;
    const char * variable;

    printf ("netCDF : setup global attributes\n");


    /*--------------------------------------------------------------------------*
     * radar                                                                    *
     *--------------------------------------------------------------------------*/
    switch (radar)
    {
    case ACROBAT_CODED_SPECTRA:
    case ACROBAT:
	variable = "ACROBAT";
	break;

    case COPERNICUS_SPECTRA:
    case COPERNICUS_SPECTRA_RAPID:
    case COPERNICUS :
	variable = "COPERNICUS";
	break;

    case GALILEO_SPECTRA:
    case GALILEO_SPECTRA_RAPID:
    case GALILEO :
	variable = "GALILEO";
	break;

    case TEST_SPECTRA:
    case TEST_SPECTRA_RAPID:
    case TEST :
    default:
	variable = "TEST-GALILEO";
	break;

    case CAMRA :
	variable = "CAMRA";
	break;
    }

    status = nc_put_att_text (ncid, NC_GLOBAL, "radar",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * source                                                                   *
     *--------------------------------------------------------------------------*/
    switch (radar)
    {
    case ACROBAT_CODED_SPECTRA:
    case ACROBAT :
	variable = "Advanced Clear-air Radar for Observing the Boundary layer And Troposphere (ACROBAT)";
	break;

    case COPERNICUS_SPECTRA:
    case COPERNICUS_SPECTRA_RAPID:
    case COPERNICUS :
	variable = "35 GHz radar (COPERNICUS)";
	break;

    case GALILEO_SPECTRA:
    case GALILEO_SPECTRA_RAPID:
    case GALILEO :
	variable = "94 GHz radar (GALILEO)";
	break;

    case TEST_SPECTRA:
    case TEST_SPECTRA_RAPID:
    case TEST:
    default:
	variable = "94 GHz radar (TEST-GALILEO)";
	break;

    case CAMRA :
	variable = "3 GHz Advanced Meteorological Radar (CAMRa)";
	break;
    }
    status = nc_put_att_text (ncid, NC_GLOBAL, "source",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * history                                                                  *
     *--------------------------------------------------------------------------*/
    time (&now);
    gmtime_r (&now, &tm);
    n = strftime (buffer, sizeof (buffer), "%a %b %d %T %Y", &tm);
    if (n < 0) n = 0;
    pt = stpcpy (buffer + n, " - ");
    for (n = 0; n < argc; n++)
    {
	pt = stpcpy (pt, argv[n]);
	pt = stpcpy (pt, " ");
    }
    status = nc_put_att_text (ncid, NC_GLOBAL, "history",
			      strlen (buffer), buffer);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * file number                                                              *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_int (ncid, NC_GLOBAL, "file_number",
			     NC_INT, 1, &scan->file_number);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * scan number                                                              *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_int (ncid, NC_GLOBAL, "scan_number",
			     NC_INT, 1, &scan->scan_number);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * scan type                                                                *
     *--------------------------------------------------------------------------*/
    switch (scan->scanType)
    {
    case SCAN_PPI: scan_type = "PPI";         break;
    case SCAN_RHI: scan_type = "RHI";         break;
    case SCAN_CSP: scan_type = "S-P";         break;
    case SCAN_FIX: scan_type = "Fixed";       break;
    case SCAN_SGL: scan_type = "Fixed";       break;
    case SCAN_CAL: scan_type = "Calibration"; break;
    case SCAN_MAN: scan_type = "Manual";      break;
    case SCAN_TRK: scan_type = "Tracking";    break;
    default:       scan_type = "\?";          break;
    }
    status = nc_put_att_text (ncid, NC_GLOBAL, "scantype",
			      strlen (scan_type), scan_type);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * experiment id                                                            *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_int (ncid, NC_GLOBAL, "experiment_id",
			     NC_INT, 1, &scan->experiment_id);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * operator                                                                 *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_text (ncid, NC_GLOBAL, "operator",
			      strlen (scan->operator), scan->operator);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * scan velocity                                                            *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_float (ncid, NC_GLOBAL, "scan_velocity",
			       NC_FLOAT, 1, &scan->scan_velocity);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * min range                                                                *
     *--------------------------------------------------------------------------*/
    temp_float = param->range[0];
    status = nc_put_att_float (ncid, NC_GLOBAL, "min_range",
			       NC_FLOAT, 1, &temp_float);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * max range                                                                *
     *--------------------------------------------------------------------------*/
    temp_float = param->range[param->samples_per_pulse - 1];
    status = nc_put_att_float (ncid, NC_GLOBAL, "max_range",
			       NC_FLOAT, 1, &temp_float);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * min angle                                                                *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_float (ncid, NC_GLOBAL, "min_angle",
			       NC_FLOAT, 1, &scan->min_angle);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * max angle                                                                *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_float (ncid, NC_GLOBAL, "max_angle",
			       NC_FLOAT, 1, &scan->max_angle);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * scan angle                                                               *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_float (ncid, NC_GLOBAL, "scan_angle",
			       NC_FLOAT, 1, &scan->scan_angle);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * scan_datetime                                                            *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_text (ncid, NC_GLOBAL, "scan_datetime",
			      strlen (scan->date), scan->date);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * ADC_bits_per_sample                                                      *
     *--------------------------------------------------------------------------*/
    temp = 12;
    status = nc_put_att_int (ncid, NC_GLOBAL, "ADC_bits_per_sample",
			     NC_INT, 1, &temp);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * samples per pulse                                                        *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_int (ncid, NC_GLOBAL, "samples_per_pulse",
			     NC_INT, 1, &param->samples_per_pulse);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * pulses per daq cycle                                                     *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_int (ncid, NC_GLOBAL, "pulses_per_daq_cycle",
			     NC_INT, 1, &param->pulses_per_daq_cycle);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * ADC channels                                                             *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_int (ncid, NC_GLOBAL, "ADC_channels",
			     NC_INT, 1, &param->ADC_channels);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * delay clocks                                                             *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_int (ncid, NC_GLOBAL, "delay_clocks",
			     NC_INT, 1, &param->delay_clocks);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * pulses per ray                                                           *
     *--------------------------------------------------------------------------*/
    if (radar == CAMRA)
    {
	temp = param->pulses_per_daq_cycle;
    }
    else
    {
	temp = param->spectra_averaged * param->pulses_per_daq_cycle;
    }
    status = nc_put_att_int (ncid, NC_GLOBAL, "pulses_per_ray",
			     NC_INT, 1, &temp);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * pulse compression                                                        *
     *--------------------------------------------------------------------------*/
    temp   = (radar == ACROBAT) ? 1 : 0;
    status = nc_put_att_int (ncid, NC_GLOBAL, "pulse_compression",
			     NC_INT, 1, &temp);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    if (radar != CAMRA)
    {
	/* pulses_coherently_averaged */
	temp = param->pulses_coherently_averaged * param->number_of_codes;
	status = nc_put_att_int (ncid, NC_GLOBAL, "pulses_coherently_averaged",
				 NC_INT, 1, &temp);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/* moments_averaged */
	status = nc_put_att_int (ncid, NC_GLOBAL, "moments_averaged",
				 NC_INT, 1, &param->moments_averaged);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/* spectra_averaged */
	status = nc_put_att_int (ncid, NC_GLOBAL, "spectra_averaged",
				 NC_INT, 1, &param->spectra_averaged);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/* fft_bins_interpolated */
	status = nc_put_att_int (ncid, NC_GLOBAL, "fft_bins_interpolated",
				 NC_INT, 1, &param->fft_bins_interpolated);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/* dump spectra */
	status = nc_put_att_int (ncid, NC_GLOBAL, "dump_spectra",
				 NC_INT, 1, &param->dump_spectra);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/* dump spectra rapid */
	status = nc_put_att_int (ncid, NC_GLOBAL, "dump_spectra_rapid",
				 NC_INT, 1, &param->dump_spectra_rapid);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
    }
    else
    {
#if 0
	/*--------------------------------------------------------------------------*
	 * extra_attenuation                                                        *
	 *--------------------------------------------------------------------------*/
	status = nc_put_att_float (ncid, NC_GLOBAL, "extra_attenuation",
				   NC_FLOAT, 1, &param->extra_attenuation);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
#endif

	/*--------------------------------------------------------------------------*
	 * radar constant                                                           *
	 *--------------------------------------------------------------------------*/
	temp_float = 64.7f;
	status = nc_put_att_float (ncid, NC_GLOBAL, "radar_constant",
				   NC_FLOAT, 1, &temp_float);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/*--------------------------------------------------------------------------*
	 * receiver gain                                                            *
	 *--------------------------------------------------------------------------*/
	temp_float = 45.5f;
	status = nc_put_att_float (ncid, NC_GLOBAL, "receiver_gain",
				   NC_FLOAT, 1, &temp_float);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/*--------------------------------------------------------------------------*
	 * cable losses                                                             *
	 *--------------------------------------------------------------------------*/
	temp_float = 4.8f;
	status = nc_put_att_float (ncid, NC_GLOBAL, "cable_losses",
				   NC_FLOAT, 1, &temp_float);
    }

    /*--------------------------------------------------------------------------*
     * year                                                                     *
     * extract year from date field in scan structure                           *
     *--------------------------------------------------------------------------*/
    strncpy (buffer, scan->date, 4);
    buffer[4] = '\0';
    temp = atoi (buffer);
    status = nc_put_att_int (ncid, NC_GLOBAL, "year", NC_INT, 1, &temp);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * month                                                                    *
     * extract month from date field in scan structure                          *
     *--------------------------------------------------------------------------*/
    strncpy (buffer, &scan->date[4], 2);
    buffer[2] = '\0';
    temp = atoi (buffer);
    status = nc_put_att_int (ncid, NC_GLOBAL, "month", NC_INT, 1, &temp);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*-------------------------------------------------------------------------*
     * day                                                                     *
     * extract day from date field in scan structure                           *
     *-------------------------------------------------------------------------*/
    strncpy (buffer, &scan->date[6], 2);
    buffer[2] = '\0';
    temp = atoi (buffer);
    status = nc_put_att_int (ncid, NC_GLOBAL, "day", NC_INT, 1, &temp);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*-------------------------------------------------------------------------*
     * British_National_Grid_Reference                                         *
     * using the standard location for Chilbolton                              *
     *-------------------------------------------------------------------------*/
    variable = "SU394386";
    status = nc_put_att_text (ncid, NC_GLOBAL, "British_National_Grid_Reference",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);
}

/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
void
RNC_SetupPulse_Compression_Code (int                     ncid,
				 const RSP_ParamStruct * param)
{

    int variable_shape[3];
    size_t variable_count[] = { 0, 0, 0 };
    size_t variable_start[] = { 0, 0, 0 };

    int code_channel_dim, code_dim, code_length_dim;
    int pccid;
    int  status;
    int  n, q;
    float temp;
    const char * variable;

    /* dimension definitions */
    status = nc_def_dim (ncid, "code_channel",
			 2, &code_channel_dim);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    status = nc_def_dim (ncid, "code",
			 param->number_of_codes, &code_dim);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    status = nc_def_dim (ncid, "code_length",
			 param->code_length, &code_length_dim);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* variable definition */
    variable_shape[0] = code_channel_dim;
    variable_shape[1] = code_dim;
    variable_shape[2] = code_length_dim;
    status = nc_def_var (ncid, "pulse_compression_code",
			 NC_FLOAT, 3, variable_shape, &pccid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "pulse compression code";
    status = nc_put_att_text (ncid, pccid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    status = nc_put_att_text (ncid, pccid, "code_name",
			      strlen (param->code_name), param->code_name);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    status = nc_put_att_float (ncid, pccid, "code_bit_period",
			       NC_FLOAT, 1, &param->pulse_period);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* variable data */

    /* go from define to data */
    status = nc_enddef (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* the code will be translated in to the netCDF variable code by code */
    /* this is to avoid any confusion */
    variable_count[0] = 1;
    variable_count[1] = 1;
    variable_count[2] = param->code_length;
    for (n = 0; n < param->number_of_codes; n++)
    {
	status = nc_put_vara_short (ncid, pccid,
				    variable_start, variable_count,
				    &param->codes[n][0]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	variable_start[1]++;
    }

    temp = 0;
    variable_start[0] = 1;
    for (n = 0; n < param->number_of_codes; n++)
    {
	variable_start[1] = n;
	for (q = 0; q < param->code_length; q++)
	{
	    variable_start[2] = q;
	    status = nc_put_var1_float (ncid, pccid,
					variable_start, &temp);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);
	}
    }

    /* go from data to define */
    status = nc_redef (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);
}

int
RNC_SetupFile_State (int ncid)
{
    int   file_stateid;
    int          variable_shape = 1;
    int   status;
    int          temp;
    const char * variable;

    /*--------------------------------------------------------------------------*
     * define variable file_state                                               *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "file_state",
			 NC_INT, 0, &variable_shape, &file_stateid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "current_state_of_file";
    status = nc_put_att_text (ncid, file_stateid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "current state of file";
    status = nc_put_att_text (ncid, file_stateid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * go from define to data                                                   *
     *--------------------------------------------------------------------------*/
    status = nc_enddef (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * set to 0 : scan in progress                                              *
     *--------------------------------------------------------------------------*/
    temp = 0;
    status = nc_put_var_int (ncid, file_stateid, &temp);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * go from data to define                                                   *
     *--------------------------------------------------------------------------*/
    status = nc_redef (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    return (file_stateid);
}

/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
void
RNC_SetupStaticVariables (int                     ncid,
			  const RSP_ParamStruct * param)
{
    int   variable_shape = 1;
    int   clockid, transmit_powerid, pulse_periodid, antenna_diameterid;
    int   beamwidthVid, beamwidthHid, prfid, frequencyid;
    int   heightid, longitudeid, latitudeid;
    int   status;
    float temp_float;
    const char * variable;

    /*--------------------------------------------------------------------------*
     * define variables                                                         *
     *--------------------------------------------------------------------------*/

    /*--------------------------------------------------------------------------*
     * latitude                                                                 *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "latitude",
			 NC_FLOAT, 0, &variable_shape, &latitudeid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "latitude";
    status = nc_put_att_text (ncid, latitudeid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "latitude of the antenna";
    status = nc_put_att_text (ncid, latitudeid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "%.4f";
    status = nc_put_att_text (ncid, latitudeid, "C_format",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "degree_north";
    status = nc_put_att_text (ncid, latitudeid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * longitude                                                                *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "longitude",
			 NC_FLOAT, 0, &variable_shape, &longitudeid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "longitude";
    status = nc_put_att_text (ncid, longitudeid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "longitude of the antenna";
    status = nc_put_att_text (ncid, longitudeid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "%.4f";
    status = nc_put_att_text (ncid, longitudeid, "C_format",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "degree_east";
    status = nc_put_att_text (ncid, longitudeid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * height                                                                   *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "height",
			 NC_FLOAT, 0, &variable_shape, &heightid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "height";
    status = nc_put_att_text (ncid, heightid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "height of the elevation axis above mean sea level (Ordnance Survey Great Britain)";
    status = nc_put_att_text (ncid, heightid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "%.0f";
    status = nc_put_att_text (ncid, heightid, "C_format",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "m";
    status = nc_put_att_text (ncid, heightid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * frequency                                                                *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "frequency",
			 NC_FLOAT, 0, &variable_shape, &frequencyid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "frequency";
    status = nc_put_att_text (ncid, frequencyid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "frequency of transmitted radiation";
    status = nc_put_att_text (ncid, frequencyid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "%.4f";
    status = nc_put_att_text (ncid, frequencyid, "C_format",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "GHz";
    status = nc_put_att_text (ncid, frequencyid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * prf                                                                      *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "prf",
			 NC_FLOAT, 0, &variable_shape, &prfid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "pulse_repetition_frequency";
    status = nc_put_att_text (ncid, prfid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "pulse repetition frequency";
    status = nc_put_att_text (ncid, prfid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "Hz";
    status = nc_put_att_text (ncid, prfid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * beamwidthH                                                               *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "beamwidthH",
			 NC_FLOAT, 0, &variable_shape, &beamwidthHid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "horizontal_angular_beamwidth";
    status = nc_put_att_text (ncid, beamwidthHid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "horizontal angular beamwidth";
    status = nc_put_att_text (ncid, beamwidthHid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "degree";
    status = nc_put_att_text (ncid, beamwidthHid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * beamwidthV                                                               *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "beamwidthV",
			 NC_FLOAT, 0, &variable_shape, &beamwidthVid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "vertical angular beamwidth";
    status = nc_put_att_text (ncid, beamwidthVid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "vertical angular beamwidth";
    status = nc_put_att_text (ncid, beamwidthVid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "degree";
    status = nc_put_att_text (ncid, beamwidthVid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * antenna diameter                                                         *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "antenna_diameter",
			 NC_FLOAT, 0, &variable_shape, &antenna_diameterid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "antenna_diameter";
    status = nc_put_att_text (ncid, antenna_diameterid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "antenna diameter";
    status = nc_put_att_text (ncid, antenna_diameterid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "m";
    status = nc_put_att_text (ncid, antenna_diameterid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* special attribute: antennae_separation */
    temp_float = 0;
    status = nc_put_att_float (ncid, antenna_diameterid, "antennae_separation",
			       NC_FLOAT, 1, &temp_float);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * pulse period                                                             *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "pulse_period",
			 NC_FLOAT, 0, &variable_shape, &pulse_periodid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "pulse_period";
    status = nc_put_att_text (ncid, pulse_periodid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "pulse period";
    status = nc_put_att_text (ncid, pulse_periodid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "us";
    status = nc_put_att_text (ncid, pulse_periodid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * transmit power                                                           *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "transmit_power",
			 NC_FLOAT, 0, &variable_shape, &transmit_powerid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "peak_transmitted_power";
    status = nc_put_att_text (ncid, transmit_powerid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "peak transmitted power";
    status = nc_put_att_text (ncid, transmit_powerid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "W";
    status = nc_put_att_text (ncid, transmit_powerid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * clock                                                                    *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "clock",
			 NC_FLOAT, 0, &variable_shape, &clockid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "clock";
    status = nc_put_att_text (ncid, clockid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "clock input to ISACTRL";
    status = nc_put_att_text (ncid, clockid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "Hz";
    status = nc_put_att_text (ncid, clockid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * special attribute: clock_divfactor                                       *
     *--------------------------------------------------------------------------*/
    status = nc_put_att_int (ncid, clockid, "clock_divfactor",
			     NC_INT, 1, &param->clock_divfactor);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * variable data                                                            *
     *--------------------------------------------------------------------------*/

    /*--------------------------------------------------------------------------*
     * go from define to data                                                   *
     *--------------------------------------------------------------------------*/
    status = nc_enddef (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* latitude */
    temp_float = 51.1445;
    status = nc_put_var_float (ncid, latitudeid, &temp_float);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* longitude */
    temp_float = 358.5630;
    status = nc_put_var_float (ncid, longitudeid, &temp_float);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* height */
    status = nc_put_var_float (ncid, heightid, &param->height);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* frequency */
    status = nc_put_var_float (ncid, frequencyid, &param->frequency);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* prf */
    status = nc_put_var_float (ncid, prfid, &param->prf);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* beamwidthH */
    status = nc_put_var_float (ncid, beamwidthHid, &param->beamwidthH);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* beamwidthV */
    status = nc_put_var_float (ncid, beamwidthVid, &param->beamwidthV);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* antenna_diameter */
    status = nc_put_var_float (ncid, antenna_diameterid, &param->antenna_diameter);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* pulse period */
    status = nc_put_var_float (ncid, pulse_periodid, &param->pulse_period);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* transmit power */
    status = nc_put_var_float (ncid, transmit_powerid, &param->transmit_power);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* clock */
    status = nc_put_var_float (ncid, clockid, &param->clock);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * go from data to define                                                   *
     *--------------------------------------------------------------------------*/
    status = nc_redef (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);
}

/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
void
RNC_SetupDimensions (int                     ncid,
		     const RSP_ParamStruct * param,
		     RNC_DimensionStruct *   dimensions)
{
    int status;

    /*--------------------------------------------------------------------------*
     * define time dimension                                                    *
     *--------------------------------------------------------------------------*/
    status = nc_def_dim (ncid, "time", NC_UNLIMITED, &dimensions->time_dim);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * define range dimension                                                   *
     *--------------------------------------------------------------------------*/
    status = nc_def_dim (ncid, "range",
			 param->samples_per_pulse,
			 &dimensions->range_dim);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* define fft_bin_dim */
    status = nc_def_dim (ncid, "fft_bin_dim",
			 param->npsd, &dimensions->fft_bin_dim);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* define coded_fft_bin_dim */
    /* this line is going to come back and haunt me */
    status = nc_def_dim (ncid, "coded_fft_bin_dim",
			 param->npsd/param->num_interleave,
			 &dimensions->coded_fft_bin_dim);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* define spectra_number_dim */
    status = nc_def_dim (ncid, "spectra_number_dim",
			 param->spectra_averaged,
			 &dimensions->spectra_number_dim);
    if (status != NC_NOERR) check_netcdf_handle_error (status);
}

/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
void
RNC_SetupRapidLogPSDDimensions (int                     ncid,
				int                     radar,
				const RSP_ParamStruct * param,
				RNC_DimensionStruct *   dimensions)
{
    int     status;

    /*--------------------------------------------------------------------------*
     * define time dimension                                                    *
     *--------------------------------------------------------------------------*/
    status = nc_def_dim (ncid, "time", NC_UNLIMITED, &dimensions->time_dim);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * define range dimension                                                   *
     *--------------------------------------------------------------------------*/
    status = nc_def_dim (ncid, "range",
			 param->samples_per_pulse, &dimensions->range_dim);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* define fft_bin_dim */
    if (radar == COPERNICUS_SPECTRA_RAPID)
    {
	/*
	 * status = nc_def_dim (ncid, "fft_bin_dim",
	 *                      param->npsd / 2, &dimensions->fft_bin_dim);
	 */
	status = nc_def_dim (ncid, "fft_bin_dim",
			     param->npsd, &dimensions->fft_bin_dim);
    }
    else
    {
	status = nc_def_dim (ncid, "fft_bin_dim",
			     param->npsd, &dimensions->fft_bin_dim);
    }
    if (status != NC_NOERR) check_netcdf_handle_error (status);
}


/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
void
RNC_SetupRange (int                         ncid,
		const RSP_ParamStruct *     param,
		const RNC_DimensionStruct * dimensions)
{

    int          variable_shape = 1;
    size_t       variable_count;
    size_t       variable_start;
    int          rangeid;
    int          status;
    const char * variable;
    int          temp_int;
    float  temp_float;

    /*--------------------------------------------------------------------------*
     * variable definition                                                      *
     *--------------------------------------------------------------------------*/
    variable_shape = dimensions->range_dim;
    status = nc_def_var (ncid, "range",
			 NC_FLOAT, 1, &variable_shape, &rangeid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "distance from the antenna to the middle of each range gate";
    status = nc_put_att_text (ncid, rangeid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "m";
    status = nc_put_att_text (ncid, rangeid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    status = nc_put_att_float (ncid, rangeid, "range_offset",
			       NC_FLOAT, 1, &param->range_offset);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    temp_int = 1;
    status = nc_put_att_int (ncid, rangeid, "gates_range_averaged",
			     NC_INT, 1, &temp_int);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    temp_float = 66;
    status = nc_put_att_float (ncid, rangeid, "gate_width",
			       NC_FLOAT, 1, &temp_float);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * go from define to data                                                   *
     *--------------------------------------------------------------------------*/
    status = nc_enddef (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /* writing out the range */
    variable_count = param->samples_per_pulse;
    variable_start = 0;
    status = nc_put_vara_float (ncid, rangeid,
				&variable_start, &variable_count,
				param->range);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * go from data to define                                                   *
     *--------------------------------------------------------------------------*/
    status = nc_redef (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

}

/*****************************************************************************/
/*!
 *
 * \fn int ObservableTemplate (int ncid, const char * name,
 *                             const int  * variable_shape,
 *                             const char * std_name,
 *                             const char * long_name,
 *                             const char * units)
 *
 * \breif Helper function for setting up observables.
 *
 * \param [in] ncid            NetCDF file handle.
 * \param [in] name            Observable/variable name.
 * \param [in] variable_shape  Shape array for variable (all have 2 dims)
 * \param [in] std_name        Chilbolton standard name.
 * \param [in] long_name       Long name.
 * \param [in] units           Units string.
 *
 * \return NetCFG veriable id.
 *
 * This function assumes that the variable_shape has been correctly set up
 * for the first two dimensions: time,range or time,unaveraged_range as
 * required.
 *
 *****************************************************************************/
static int
ObservableTemplate (int          ncid,
		    const char * name,
		    const int  * variable_shape,
		    const char * std_name,
		    const char * long_name,
		    const char * units)
{
    int varid;
    int status;

    status = nc_def_var (ncid, name, NC_FLOAT, 2, variable_shape, &varid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    if (std_name != NULL)
    {
	status = nc_put_att_text (ncid, varid, "chilbolton_standard_name",
				  strlen (std_name), std_name);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
    }
    if (long_name != NULL)
    {
	status = nc_put_att_text (ncid, varid, "long_name",
				  strlen (long_name), long_name);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
    }
    if (units != NULL)
    {
	status = nc_put_att_text (ncid, varid, "units",
				  strlen (units), units);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
    }
    return varid;
}


static int
ObservableTemplateScalar (int          ncid,
			  const char * name,
			  const int  * variable_shape,
			  const char * std_name,
			  const char * long_name,
			  const char * units)
{
    int varid;
    int status;

    status = nc_def_var (ncid, name, NC_FLOAT, 1, variable_shape, &varid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    if (std_name != NULL)
    {
	status = nc_put_att_text (ncid, varid, "chilbolton_standard_name",
				  strlen (std_name), std_name);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
    }
    if (long_name != NULL)
    {
	status = nc_put_att_text (ncid, varid, "long_name",
				  strlen (long_name), long_name);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
    }
    if (units != NULL)
    {
	status = nc_put_att_text (ncid, varid, "units",
				  strlen (units), units);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
    }
    return varid;
}

/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
void
RNC_SetupDynamicVariables (int                         ncid,
			   int                         radar,
			   const URC_ScanStruct *      scan,
			   const RSP_ParamStruct *     param,
			   const RNC_DimensionStruct * dimensions,
			   RSP_ObservablesStruct *     obs)
{
    int          variable_shape[2];
    int          status;
    int          n;
    char         buffer [80];
    const char * variable;

    float temp_float = -999.0f;

    variable_shape[0] = dimensions->time_dim;
    variable_shape[1] = dimensions->range_dim;

    /*--------------------------------------------------------------------------*
     * time definition                                                          *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "time",
			 NC_FLOAT, 1, variable_shape, &obs->tsid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "time";
    status = nc_put_att_text (ncid, obs->tsid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "time";
    status = nc_put_att_text (ncid, obs->tsid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "%.2f";
    status = nc_put_att_text (ncid, obs->tsid, "C_format",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    sprintf (buffer, "seconds since %c%c%c%c-%c%c-%c%c 00:00:00 +00:00",
	     scan->date[0], scan->date[1], scan->date[2], scan->date[3],
	     scan->date[4], scan->date[5],
	     scan->date[6], scan->date[7]);
    status = nc_put_att_text (ncid, obs->tsid, "units",
			      strlen (buffer), buffer);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * dish_time definition                                                     *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "dish_time",
			 NC_FLOAT, 1, variable_shape, &obs->dish_tsid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "dish_time";
    status = nc_put_att_text (ncid, obs->dish_tsid, "chilbolton_standard_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "dish_time";
    status = nc_put_att_text (ncid, obs->dish_tsid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "%.2f";
    status = nc_put_att_text (ncid, obs->tsid, "C_format",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    sprintf (buffer, "seconds since %04d-%02d-%02d 00:00:00 +00:00",
	     obs->dish_year,
	     obs->dish_month,
	     obs->dish_day);
    status = nc_put_att_text (ncid, obs->dish_tsid, "units",
			      strlen (buffer), buffer);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * elevation definition                                                     *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "elevation",
			 NC_FLOAT, 1, variable_shape, &obs->elevationid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "elevation angle above the horizon at the start of the beamwidth";
    status = nc_put_att_text (ncid, obs->elevationid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "%.3f";
    status = nc_put_att_text (ncid, obs->elevationid, "C_format",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "degree";
    status = nc_put_att_text (ncid, obs->elevationid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);


    /*--------------------------------------------------------------------------*
     * azimuth definition                                                       *
     *--------------------------------------------------------------------------*/
    status = nc_def_var (ncid, "azimuth",
			 NC_FLOAT, 1, variable_shape, &obs->azimuthid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "azimuth angle clockwise from the grid north at the start of the beamwidth";
    status = nc_put_att_text (ncid, obs->azimuthid, "long_name",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "%.3f";
    status = nc_put_att_text (ncid, obs->azimuthid, "C_format",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    variable = "degree";
    status = nc_put_att_text (ncid, obs->azimuthid, "units",
			      strlen (variable), variable);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    status = nc_put_att_float (ncid, obs->azimuthid, "azimuth_offset",
			       NC_FLOAT, 1, &param->azimuth_offset);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    switch (radar)
    {
    case COPERNICUS_SPECTRA_RAPID :
    case GALILEO_SPECTRA_RAPID :
    case TEST_SPECTRA_RAPID :

	/* bin_number */
	status = nc_def_var (ncid, "bin_number",
			     NC_INT, 1, variable_shape, &obs->bin_numberid);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	variable = "bin_number";
	status = nc_put_att_text (ncid, obs->bin_numberid, "long_name",
				  strlen (variable), variable);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	variable = "1";
	status = nc_put_att_text (ncid, obs->bin_numberid, "units",
				  strlen (variable), variable);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/* ray_number */
	status = nc_def_var (ncid, "ray_number",
			     NC_INT, 1, variable_shape, &obs->ray_numberid);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	variable = "ray_number";
	status = nc_put_att_text (ncid, obs->ray_numberid, "long_name",
				  strlen (variable), variable);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	variable = "1";
	status = nc_put_att_text (ncid, obs->ray_numberid, "units",
				  strlen (variable), variable);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
	break;

    case ACROBAT :
    case COPERNICUS :
    case GALILEO :
    case CAMRA :
    case TEST :
	printf ("Number of observations %d \n", obs->n_obs);
	for (n = 0; n < obs->n_obs; n++)
	{
	    if (!obs->record_observable[n])
		continue;

	    if (!strcmp (obs->name[n], "SNR_HC"))
	    {
		/*--------------------------------------------------------------*
		 * SNR_HC definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "signal_to_noise_ratio_at_horizontal_polarisation_from_spectral_processing",
		    "signal to noise ratio at horizontal polarisation from spectral processing",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "SNR_HCP"))
	    {
		/*--------------------------------------------------------------*
		 * SNR_HCP definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "signal_to_noise_ratio_at_horizontal_polarisation_from_spectral_processing_using_pulse_compression",
		    "signal to noise ratio at horizontal polarisation from spectral processing using pulse compression",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "SNR_VC"))
	    {
		/*--------------------------------------------------------------*
		 * SNR_VC definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "signal_to_noise_ratio_at_vertical_polarisation_from_spectral_processing",
		    "signal to noise ratio at vertical polarisation from spectral processing",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "SNR_VCP"))
	    {
		/*--------------------------------------------------------------*
		 * SNR_VCP definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "signal_to_noise_ratio_at_vertical_polarisation_from_spectral_processing_using_pulse_compression",
		    "signal to noise ratio at vertical polarisation from spectral processing using pulse compression",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "SNR_XHC"))
	    {
		/*--------------------------------------------------------------*
		 * SNR_XHC definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "crosspolar_signal_to_noise_ratio_at_horizontal_polarisation_from_spectral_processing",
		    "crosspolar signal to noise ratio at horizontal polarisation from spectral processing",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "SNR_XHCP"))
	    {
		/*--------------------------------------------------------------*
		 * SNR_XHCP definition                                          *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "crosspolar_signal_to_noise_ratio_at_horizontal_polarisation_from_spectral_processing_using_pulse_compression",
		    "crosspolar signal to noise ratio at horizontal polarisation from spectral processing using pulse compression",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "SNR_XVC"))
	    {
		/*--------------------------------------------------------------*
		 * SNR_XVC definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "crosspolar_signal_to_noise_ratio_at_vertical_polarisation_from_spectral_processing",
		    "crosspolar signal to noise ratio at vertical polarisation from spectral processing",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "SNR_XVCP"))
	    {
		/*--------------------------------------------------------------*
		 * SNR_XVCP definition                                          *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "crosspolar_signal_to_noise_ratio_at_vertical_polarisation_from_spectral_processing_using_pulse_compression",
		    "crosspolar signal to noise ratio at vertical polarisation from spectral processing using pulse compression",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "ZED_H"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_H definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radar_reflectivity_factor_at_horizontal_polarisation",
		    "radar reflectivity factor at horizontal polarisation",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_incoherent_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_V"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_V definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radar_reflectivity_factor_at_vertical_polarisation",
		    "radar reflectivity factor at vertical polarisation",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_incoherent_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_XV"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_XV definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "crosspolar_radar_reflectivity_factor_at_vertical_polarisation",
		    "crosspolar radar reflectivity factor at vertical polarisation",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_crosspolar_incoherent_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_XH"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_XH definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "crosspolar_radar_reflectivity_factor_at_horizontal_polarisation",
		    "crosspolar radar reflectivity factor at horizontal polarisation",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_crosspolar_incoherent_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "CXC"))
	    {
		/*--------------------------------------------------------------*
		 * CXC definition                                               *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "copolar_cross_correlation_coefficient",
		    "copolar cross correlation coefficient",
		    "");
	    }
	    else if (!strcmp (obs->name[n], "ZED_HP"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_HP definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radar_reflectivity_factor_at_horizontal_polarisation_using_pulse_compression",
		    "radar reflectivity factor at horizontal polarisation using pulse compression",
		    "dBZ");
	    }
	    else if (!strcmp (obs->name[n], "ZED_HC"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_HC definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radar_reflectivity_factor_at_horizontal_polarisation_from_spectral_processing",
		    "radar reflectivity factor at horizontal polarisation from spectral processing",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_calibration_offset); // ZED_HC_calibration_offset
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_HCD"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_HCD definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "fractional_standard_deviation_of_radar_reflectivity_factor_at_horizontal_polarisation_from_spectral_processing",
		    "fractional standard deviation of radar reflectivity factor at horizontal polarisation from spectral processing",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_VCD"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_VCD definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "fractional_standard_deviation_of_radar_reflectivity_factor_at_vertical_polarisation_from_spectral_processing",
		    "fractional standard deviation of radar reflectivity factor at vertical polarisation from spectral processing",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_HCP"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_HCP definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radar_reflectivity_factor_at_horizontal_polarisation_from_spectral_processing_using_pulse_compression",
		    "radar reflectivity factor at horizontal polarisation from spectral processing using pulse compression",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_calibration_offset);  // ZED_HPC_calibration_offset
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_HCDP"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_HCDP definition                                          *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "fractional_standard_deviation_of_radar_reflectivity_factor"
		    "_at_horizontal_polarisation_from_spectral_processing_using_pulse_compression",
		    "fractional standard deviation of radar reflectivity factor"
		    " at horizontal polarisation from spectral processing using pulse compression",
		    "dB");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_VCDP"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_VCDP definition                                          *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "fractional_standard_deviation_of_radar_reflectivity_factor"
		    "_at_vertical_polarisation_from_spectral_processing_using_pulse_compression",
		    "fractional standard deviation of radar reflectivity factor"
		    " at vertical polarisation from spectral processing using pulse compression",
		    "dB");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_XHC"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_XHC definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "crosspolar_radar_reflectivity_factor_at_horizontal_polarisation_from_spectral_processing",
		    "crosspolar radar reflectivity factor at horizontal polarisation from spectral processing",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_calibration_offset); // ZED_XHC_calibration_offset
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_XHCP"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_XHCP definition                                          *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "crosspolar_radar_reflectivity_factor_at_horizontal_polarisation_from_spectral_processing_using_pulse_compression",
		    "crosspolar radar reflectivity factor at horizontal polarisation from spectral processing using pulse compression",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_calibration_offset); // ZED_XHCP_calibration_offset
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_VC"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_VC definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radar_reflectivity_factor_at_vertical_polarisation_from_spectral_processing",
		    "radar reflectivity factor at vertical polarisation from spectral processing",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset"
					   , NC_FLOAT, 1,
					   &param->ZED_calibration_offset); // ZED_VC_calibration_offset
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_VCP"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_VCP definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radar_reflectivity_factor_at_vertical_polarisation_from_spectral_processing_using_pulse_compression",
		    "radar reflectivity factor at vertical polarisation from spectral processing using pulse compression",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_calibration_offset); // ZED_VCP_calibration_offset
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_XVC"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_XVC definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "crosspolar_radar_reflectivity_factor_at_vertical_polarisation_from_spectral_processing",
		    "crosspolar radar reflectivity factor at vertical polarisation from spectral processing",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_calibration_offset); // ZED_XVC_calibration_offset
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZED_XVCP"))
	    {
		/*--------------------------------------------------------------*
		 * ZED_XVCP definition                                          *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "crosspolar_radar_reflectivity_factor_at_vertical_polarisation_from_spectral_processing_using_pulse_compression",
		    "crosspolar radar reflectivity factor at vertical polarisation from spectral processing using pulse compression",
		    "dBZ");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZED_calibration_offset); // ZED_XVCP_calibration_offset
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "POW_H"))
	    {
		/*--------------------------------------------------------------*
		 * POW_H definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "copolar_incoherent_power_at_horizontal_polarisation",
		    "copolar incoherent power at horizontal polarisation",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "POW_HX"))
	    {
		/*--------------------------------------------------------------*
		 * POW_HX definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "crosspolar_incoherent_power_at_horizontal_polarisation",
		    "crosspolar incoherent power at horizontal polarisation",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "POW_V"))
	    {
		/*--------------------------------------------------------------*
		 * POW_V definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "copolar_incoherent_power_at_vertical_polarisation",
		    "copolar incoherent power at vertical polarisation",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "POW_VX"))
	    {
		/*--------------------------------------------------------------*
		 * POW_VX definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "crosspolar_incoherent_power_at_vertical_polarisation",
		    "crosspolar incoherent power at vertical polarisation",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "SPW_H"))
	    {
		/*--------------------------------------------------------------*
		 * SPW_H definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "doppler_spectral_width_at_horizontal_polarisation",
		    "doppler spectral width at horizontal polarisation",
		    "m s-1");
	    }
	    else if (!strcmp (obs->name[n], "SPW_V"))
	    {
		/*--------------------------------------------------------------*
		 * SPW_V definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "doppler_spectral_width_at_vertical_polarisation",
		    "doppler spectral width at vertical polarisation",
		    "m s-1");
	    }
	    else if (!strcmp (obs->name[n], "SPW_HV"))
	    {
		/*--------------------------------------------------------------*
		 * SPW_HV definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "doppler_spectral_width_at_horizontal_and_vertical_polarisation",
		    "doppler spectral width at horizontal and vertical polarisation",
		    "m s-1");
	    }
	    else if (!strcmp (obs->name[n], "SPW_HC"))
	    {
		/*--------------------------------------------------------------*
		 * SPW_HC definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "doppler_spectral_width_at_horizontal_polarisation_from_spectral_processing",
		    "doppler spectral width at horizontal polarisation from spectral processing",
		    "m s-1");
	    }
	    else if (!strcmp (obs->name[n], "SPW_VC"))
	    {
		/*--------------------------------------------------------------*
		 * SPW_VC definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "doppler_spectral_width_at_vertical_polarisation_from_spectral_processing",
		    "doppler spectral width at vertical polarisation from spectral processing",
		    "m s-1");
	    }
	    else if (!strcmp (obs->name[n], "SPW_HCP"))
	    {
		/*--------------------------------------------------------------*
		 * SPW_HCP definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "doppler_spectral_width_at_horizontal_polarisation_from_spectral_processing_using_pulse_compression",
		    "doppler spectral width at horizontal polarisation from spectral processing using pulse compression",
		    "m s-1");
	    }
	    else if (!strcmp (obs->name[n], "SPW_VCP"))
	    {
		/*--------------------------------------------------------------*
		 * SPW_VCP definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "doppler_spectral_width_at_vertical_polarisation_from_spectral_processing_using_pulse_compression",
		    "doppler spectral width at vertical polarisation from spectral processing using pulse compression",
		    "m s-1");
	    }
	    else if (!strcmp (obs->name[n], "VEL_H"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_H definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radial_velocity_of_scatterers_away_from_instrument_at_horizontal_polarisation",
		    "radial velocity of scatterers away from instrument at horizontal polarisation",
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "VEL_V"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_V definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radial_velocity_of_scatterers_away_from_instrument_at_vertical_polarisation",
		    "radial velocity of scatterers away from instrument at vertical polarisation",
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "VEL_HV"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_HV definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radial_velocity_of_scatterers_away_from_instrument_at_horizontal_and_vertical_polarisation",
		    "radial velocity of scatterers away from instrument at horizontal and vertical polarisation",
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "VEL_HC"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_HC definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radial_velocity_of_scatterers_away_from_instrument_at_horizontal_polarisation_from_spectral_processing",
		    "radial velocity of scatterers away from instrument at horizontal polarisation from spectral processing",
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "VEL_VC"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_VC definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radial_velocity_of_scatterers_away_from_instrument_at_vertical_polarisation_from_spectral_processing",
		    "radial velocity of scatterers away from instrument at vertical polarisation from spectral processing",
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "VEL_HCP"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_HCP definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radial_velocity_of_scatterers_away_from_instrument"
		    "_at_horizontal_polarisation_from_spectral_processing_using_pulse_compression",
		    "radial velocity of scatterers away from instrument"
		    " at horizontal polarisation from spectral processing using pulse compression",
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "VEL_VCP"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_VCP definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radial_velocity_of_scatterers_away_from_instrument"
		    "_at_vertical_polarisation_from_spectral_processing_using_pulse_compression",
		    "radial velocity of scatterers away from instrument"
		    " at vertical polarisation from spectral processing using pulse compression",
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "VEL_VD"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_VD definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radial_velocity_of_scatterers_away_from_instrument_using_variable_delay_HV_pulse_pair_processing",
		    "radial velocity of scatterers away from instrument using variable delay pulse pair processing",
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "VEL_FD"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_FD definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "radial_velocity_of_scatterers_away_from_instrument_using_fixed_delay_HV_pulse_pair_processing",
		    "radial velocity of scatterers away from instrument using fixed delay HV pulse pair processing",
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "VEL_HCD"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_HCD definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "standard_deviation_of_radial_velocity_of_scatterers_away"
		    "_from_instrument_at_horizontal_polarisation_from_spectral_processing",
		    "standard deviation of radial velocity of scatterers away"
		    " from instrument at horizontal polarisation from spectral processing",
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "VEL_VCD"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_VCD definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "standard_deviation_of_radial_velocity_of_scatterers_away"
		    "_from_instrument_at_vertical_polarisation_from_spectral_processing",
		    "standard deviation of radial velocity of scatterers away"
		    " from instrument at vertical polarisation from spectral processing",
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }

	    else if (!strcmp (obs->name[n], "VEL_HCDP"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_HCDP definition                                          *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "standard_deviation_of_radial_velocity_of_scatterers_away_from"
		    "_instrument_at_horizontal_polarisation_from_spectral_processing_using_pulse_compression",
		    "standard deviation of radial velocity of scatterers away from"
		    " instrument at horizontal polarisation from spectral processing using pulse compression",
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "VEL_VCDP"))
	    {
		/*--------------------------------------------------------------*
		 * VEL_VCDP definition                                          *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "standard_deviation_of_radial_velocity_of_scatterers_away_from"
		    "_instrument_at_vertical_polarisation_from_spectral_processing_using_pulse_compression",
		    "standard deviation of radial velocity of scatterers away from"
		    " instrument at vertical polarisation from spectral processing using pulse compression"      ,
		    "m s-1");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "folding_velocity",
					   NC_FLOAT, 1,
					   &param->folding_velocity);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "PHI_H"))
	    {
		/*--------------------------------------------------------------*
		 * PHI_H definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "absolute_phase_at_horizontal_polarisation",
		    "absolute phase at horizontal polarisation",
		    "degree");
	    }
	    else if (!strcmp (obs->name[n], "PHI_V"))
	    {
		/*--------------------------------------------------------------*
		 * PHI_V definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "absolute_phase_at_vertical_polarisation",
		    "absolute phase at vertical polarisation",
		    "degree");
	    }
	    else if (!strcmp (obs->name[n], "PHI_HV"))
	    {
		/*--------------------------------------------------------------*
		 * PHI_HV definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "absolute_phase_at_horizontal_and_vertical_polarisation",
		    "absolute phase at horizontal and vertical polarisation",
		    "degree");
	    }
	    else if (!strcmp (obs->name[n], "PHI_HD"))
	    {
		/*--------------------------------------------------------------*
		 * PHI_HD definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "standard_deviation_of_absolute_phase_at_horizontal_polarisation",
		    "standard_deviation_of_absolute phase at horizontal polarisation",
		    "degree");
	    }
	    else if (!strcmp (obs->name[n], "PHI_VD"))
	    {
		/*--------------------------------------------------------------*
		 * PHI_VD definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "standard_deviation_of_absolute_phase_at_vertical_polarisation",
		    "standard_deviation_of_absolute phase at vertical polarisation",
		    "degree");
	    }
	    else if (!strcmp (obs->name[n], "PHI_HVD"))
	    {
		/*--------------------------------------------------------------*
		 * PHI_HVD definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "standard_deviation_of_absolute_phase_at_horizontal_and_vertical_polarisation",
		    "standard_deviation_of_absolute phase at horizontal and vertical polarisation",
		    "degree");
	    }
	    else if (!strcmp (obs->name[n], "PDP"))
	    {
		/*--------------------------------------------------------------*
		 * PDP definition                                               *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "differential_phase_shift",
		    "differential phase shift",
		    "degree");
	    }
	    else if (!strcmp (obs->name[n], "ZDR"))
	    {
		/*--------------------------------------------------------------*
		 * ZDR definition                                               *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "differential_reflectivity",
		    "differential reflectivity",
		    "dB");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZDR_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZDR_C"))
	    {
		/*--------------------------------------------------------------*
		 * ZDR_C definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "differential_reflectivity_from_spectral_processing",
		    "differential reflectivity from spectral processing",
		    "dB");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZDR_C_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "ZDR_CP"))
	    {
		/*--------------------------------------------------------------*
		 * ZDR_CP definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "differential_reflectivity_from_spectral_processing_using_pulse_compression",
		    "differential reflectivity from spectral processing using pulse compression",
		    "dB");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->ZDR_CP_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "LDR"))
	    {
		/*--------------------------------------------------------------*
		 * LDR definition                                               *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "linear_depolarisation_ratio",
		    "linear depolarisation ratio",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "LDR_C"))
	    {
		/*--------------------------------------------------------------*
		 * LDR_C definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "linear_depolarisation_ratio_from_spectral_processing",
		    "linear depolarisation ratio from spectral processing",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "LDR_HC"))
	    {
		/*--------------------------------------------------------------*
		 * LDR_HC definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "linear_depolarisation_ratio_at_horizontal_polarisation_from_spectral_processing",
		    "linear depolarisation ratio at horizontal polarisation from spectral processing",
		    "dB");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->LDR_HC_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "LDR_VC"))
	    {
		/*--------------------------------------------------------------*
		 * LDR_VC definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "linear_depolarisation_ratio_at_vertical_polarisation_from_spectral_processing",
		    "linear depolarisation ratio at vertical polarisation from spectral processing",
		    "dB");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->LDR_VC_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "LDR_CP"))
	    {
		/*--------------------------------------------------------------*
		 * LDR_CP definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "linear_depolarisation_ratio_from_spectral_processing_using_pulse_compression",
		    "linear depolarisation ratio from spectral processing using pulse compression",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "LDR_HCP"))
	    {
		/*--------------------------------------------------------------*
		 * LDR_HCP definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "linear_depolarisation_ratio_at_horizontal_polarisation_from_spectral_processing_using_pulse_compression",
		    "linear depolarisation ratio at horizontal polarisation from spectral processing using pulse compression",
		    "dB");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->LDR_HCP_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "LDR_VCP"))
	    {
		/*--------------------------------------------------------------*
		 * LDR_VCP definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "linear_depolarisation_ratio_at_vertical_polarisation_from_spectral_processing_using_pulse_compression",
		    "linear depolarisation ratio at vertical polarisation from spectral processing using pulse compression",
		    "dB");

		status = nc_put_att_float (ncid, obs->varid[n],
					   "applied_calibration_offset",
					   NC_FLOAT, 1,
					   &param->LDR_VCP_calibration_offset);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	    else if (!strcmp (obs->name[n], "PHIDP_C"))
	    {
		/*--------------------------------------------------------------*
		 * PHIDP_C definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "differential_phase_shift",
		    "differential phase shift",
		    "degree");
	    }
	    else if (!strcmp (obs->name[n], "PHIDP_CP"))
	    {
		/*--------------------------------------------------------------*
		 * PHIDP_CP definition                                          *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "differential_phase_shift_using_pulse_compression",
		    "differential phase shift using pulse compression",
		    "degree");
	    }
	    else if (!strcmp (obs->name[n], "PHIDP_FD"))
	    {
		/*--------------------------------------------------------------*
		 * PHIDP_FD definition                                          *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "differential_phase_shift (fixed delay)",
		    "differential phase shift (fixed delay)",
		    "degree");
	    }
	    else if (!strcmp (obs->name[n], "PHIDP_VD"))
	    {
		/*--------------------------------------------------------------*
		 * PHIDP_VD definition                                          *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "differential_phase_shift (variable delay)",
		    "differential phase shift (variable delay)",
		    "degree");
	    }
	    else if (!strcmp (obs->name[n], "RHOHV_C"))
	    {
		/*--------------------------------------------------------------*
		 * RHOHV_C definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "copolar_cross_correlation_coefficient",
		    "copolar cross correlation coefficient",
		    "");
	    }
	    else if (!strcmp (obs->name[n], "RHOHV_CP"))
	    {
		/*--------------------------------------------------------------*
		 * RHOHV_CP definition                                          *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "copolar_cross_correlation_coefficient_using_pulse_compression",
		    "copolar cross correlation coefficient using pulse compression",
		    "");
	    }
	    else if (!strcmp (obs->name[n], "RHO_VD"))
	    {
		/*--------------------------------------------------------------*
		 * RHO_VD definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "cross_correlation_coefficient_using_variable_delay_HV_pulse_pair_processing",
		    "cross correlation coefficient using variable delay pulse pair processing",
		    "");
	    }
	    else if (!strcmp (obs->name[n], "RHO_FD"))
	    {
		/*--------------------------------------------------------------*
		 * RHO_FD definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "cross_correlation_coefficient_using_fixed_delay_HV_pulse_pair_processing",
		    "cross correlation coefficient using fixed delay pulse pair processing",
		    "");
	    }
	    else if (!strcmp (obs->name[n], "RHO_VDS"))
	    {
		/*--------------------------------------------------------------*
		 * RHO_VDS definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "cross_correlation_coefficient_using_variable_delay_HV_pulse_pair_processing",
		    "cross correlation coefficient using variable delay pulse pair processing (HV+VH)",
		    "");
	    }
	    else if (!strcmp (obs->name[n], "RHO_FDS"))
	    {
		/*--------------------------------------------------------------*
		 * RHO_FDS definition                                           *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "cross_correlation_coefficient_using_fixed_delay_HV_pulse_pair_processing",
		    "cross correlation coefficient using fixed delay pulse pair processing (HV+VH)",
		    "");
	    }
	    else if (!strcmp (obs->name[n], "NPC_H"))
	    {
		/*--------------------------------------------------------------*
		 * NPC_H definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "copolar_noise_power_counts",
		    "copolar noise power ADC counts",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "NPC_V"))
	    {
		/*--------------------------------------------------------------*
		 * NPC_V definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    "cross_polar_noise_power_counts",
		    "crosspolar noise power ADC counts",
		    "dB");
	    }
	    else if (!strcmp (obs->name[n], "I_UNCH"))
	    {
		/*--------------------------------------------------------------*
		 * I_UNCH definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    NULL,
		    NULL,
		    "");
	    }
	    else if (!strcmp (obs->name[n], "Q_UNCH"))
	    {
		/*--------------------------------------------------------------*
		 * Q_UNCH definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    NULL,
		    NULL,
		    "");
	    }
	    else if (!strcmp (obs->name[n], "I_UNCV"))
	    {
		/*--------------------------------------------------------------*
		 * I_UNCV definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    NULL,
		    NULL,
		    "");
	    }
	    else if (!strcmp (obs->name[n], "Q_UNCV"))
	    {
		/*--------------------------------------------------------------*
		 * Q_UNCV definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    NULL,
		    NULL,
		    "");
	    }
	    else if (!strcmp (obs->name[n], "I_CODH"))
	    {
		/*--------------------------------------------------------------*
		 * I_CODH definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    NULL,
		    NULL,
		    "");
	    }
	    else if (!strcmp (obs->name[n], "Q_CODH"))
	    {
		/*--------------------------------------------------------------*
		 * Q_CODH definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    NULL,
		    NULL,
		    "");
	    }
	    else if (!strcmp (obs->name[n], "I_CODV"))
	    {
		/*--------------------------------------------------------------*
		 * I_CODV definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    NULL,
		    NULL,
		    "");
	    }
	    else if (!strcmp (obs->name[n], "Q_CODV"))
	    {
		/*--------------------------------------------------------------*
		 * Q_CODV definition                                            *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n],
		    variable_shape,
		    NULL,
		    NULL,
		    "");
	    }
	    else if (!strcmp (obs->name[n], "TX_1A"))
	    {
		/*--------------------------------------------------------------*
		 * TX_1A definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplateScalar (
		    ncid, obs->name[n],
		    variable_shape,
		    "tx_power_eika_A",
		    "tx power eika A",
		    "mv");
	    }
	    else if (!strcmp (obs->name[n], "TX_1B"))
	    {
		/*--------------------------------------------------------------*
		 * TX_1B definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplateScalar (
		    ncid, obs->name[n],
		    variable_shape,
		    "tx_power_eika_B",
		    "tx power eika B",
		    "mv");
	    }
	    else if (!strcmp (obs->name[n], "TX_2A"))
	    {
		/*--------------------------------------------------------------*
		 * TX_2A definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplateScalar (
		    ncid, obs->name[n],
		    variable_shape,
		    "tx_power_after_radome_A",
		    "tx power after radome A",
		    "mv");
	    }
	    else if (!strcmp (obs->name[n], "TX_2B"))
	    {
		/*--------------------------------------------------------------*
		 * TX_2B definition                                             *
		 *--------------------------------------------------------------*/
		obs->varid[n] = ObservableTemplateScalar (
		    ncid, obs->name[n],
		    variable_shape,
		    "tx_power_after_radome_B",
		    "tx power after radome B",
		    "mv");
	    }
	    else
	    {
		printf ("Unkown observable %s. Using default template\n",
			obs->name[n]);

		obs->varid[n] = ObservableTemplate (
		    ncid, obs->name[n], variable_shape,
		    obs->name[n], obs->name[n], "");
	    }

	    /*------------------------------------------------------------------*
	     * apply missing value                                              *
	     *------------------------------------------------------------------*/
	    printf ("applying missing value to: %s\n", obs->name[n]);

	    status = nc_put_att_float (ncid, obs->varid[n], "missing_value",
				       NC_FLOAT, 1, &temp_float);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    status = nc_put_att_float (ncid, obs->varid[n], "_FillValue",
				       NC_FLOAT, 1, &temp_float);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);
	}
	break;
    }
}

/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
void
RNC_SetupLogPSDVariables (int                         ncid,
			  int                         radar,
			  const RSP_ParamStruct *     param,
			  const RNC_DimensionStruct * dimensions,
			  int *                       PSD_varid)
{
    int     variable_shape[4];
    int     status;

    variable_shape[0] = dimensions->time_dim;
    variable_shape[1] = dimensions->range_dim;
    variable_shape[2] = dimensions->fft_bin_dim;

    switch (radar)
    {
    case COPERNICUS_CODED_SPECTRA :
	status = nc_def_var (ncid, "PSD_HH",
			     NC_SHORT, 3, variable_shape,
			     &PSD_varid[PSD_HH]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	status = nc_def_var (ncid, "PSD_HV",
			     NC_SHORT, 3, variable_shape,
			     &PSD_varid[PSD_HV]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	if (param->num_tx_pol == 2)
	{
	    status = nc_def_var (ncid, "PSD_VV",
				 NC_SHORT, 3, variable_shape,
				 &PSD_varid[PSD_VV]);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    status = nc_def_var (ncid, "PSD_VH",
				 NC_SHORT, 3, variable_shape,
				 &PSD_varid[PSD_VH]);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);
	}

	status = nc_def_var (ncid, "PSD_HHP",
			     NC_SHORT, 3, variable_shape,
			     &PSD_varid[PSD_HHP]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	status = nc_def_var (ncid, "PSD_HVP",
			     NC_SHORT, 3, variable_shape,
			     &PSD_varid[PSD_HVP]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	if (param->num_tx_pol == 2)
	{
	    status = nc_def_var (ncid, "PSD_VVP",
				 NC_SHORT, 3, variable_shape,
				 &PSD_varid[PSD_VVP]);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    status = nc_def_var (ncid, "PSD_VHP",
				 NC_SHORT, 3, variable_shape,
				 &PSD_varid[PSD_VHP]);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);
	}

	if (param->include_iq_in_spectra == 1)
	{
	    variable_shape[1] = dimensions->spectra_number_dim;
	    variable_shape[2] = dimensions->fft_bin_dim;
	    variable_shape[3] = dimensions->range_dim;
	    status = nc_def_var (ncid, "IPF_HH",
				 NC_SHORT, 4, variable_shape,
				 &PSD_varid[IPF_HH]);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    status = nc_def_var (ncid, "IPF_HV",
				 NC_SHORT, 4, variable_shape,
				 &PSD_varid[IPF_HV]);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    variable_shape[2] = dimensions->coded_fft_bin_dim;
	    status = nc_def_var (ncid, "IPF_HHP",
				 NC_LONG, 4, variable_shape,
				 &PSD_varid[IPF_HHP]);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    status = nc_def_var (ncid, "IPF_HVP",
				 NC_LONG, 4, variable_shape,
				 &PSD_varid[IPF_HVP]);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    variable_shape[2] = dimensions->fft_bin_dim;
	    status = nc_def_var (ncid, "QPF_HH",
				 NC_SHORT, 4, variable_shape,
				 &PSD_varid[QPF_HH]);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    status = nc_def_var (ncid, "QPF_HV",
				 NC_SHORT, 4, variable_shape,
				 &PSD_varid[QPF_HV]);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    variable_shape[2] = dimensions->coded_fft_bin_dim;
	    status = nc_def_var (ncid, "QPF_HHP",
				 NC_LONG, 4, variable_shape,
				 &PSD_varid[QPF_HHP]);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    status = nc_def_var (ncid, "QPF_HVP",
				 NC_LONG, 4, variable_shape,
				 &PSD_varid[QPF_HVP]);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    if (param->num_tx_pol == 2)
	    {
		variable_shape[2] = dimensions->fft_bin_dim;
		status = nc_def_var (ncid, "IPF_VV",
				     NC_SHORT, 4, variable_shape,
				     &PSD_varid[IPF_VV]);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		status = nc_def_var (ncid, "IPF_VH",
				     NC_SHORT, 4, variable_shape,
				     &PSD_varid[IPF_VH]);
		variable_shape[2] = dimensions->coded_fft_bin_dim;
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		status = nc_def_var (ncid, "IPF_VVP",
				     NC_LONG, 4, variable_shape,
				     &PSD_varid[IPF_VVP]);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		status = nc_def_var (ncid, "IPF_VHP",
				     NC_LONG, 4, variable_shape,
				     &PSD_varid[IPF_VHP]);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		variable_shape[2] = dimensions->fft_bin_dim;
		status = nc_def_var (ncid, "QPF_VV",
				     NC_SHORT, 4, variable_shape,
				     &PSD_varid[QPF_VV]);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		status = nc_def_var (ncid, "QPF_VH",
				     NC_SHORT, 4, variable_shape,
				     &PSD_varid[QPF_VH]);
		variable_shape[2] = dimensions->coded_fft_bin_dim;
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		status = nc_def_var (ncid, "QPF_VVP",
				     NC_LONG, 4, variable_shape,
				     &PSD_varid[QPF_VVP]);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		status = nc_def_var (ncid, "QPF_VHP",
				     NC_LONG, 4, variable_shape,
				     &PSD_varid[QPF_VHP]);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	}
	break;

    case COPERNICUS_SPECTRA_RAPID :
	variable_shape[0] = dimensions->time_dim;
	variable_shape[1] = dimensions->fft_bin_dim;
	status = nc_def_var (ncid, "PSD_RAPID_HHP",
			     NC_SHORT, 2, variable_shape,
			     &PSD_varid[0]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	status = nc_def_var (ncid, "PSD_RAPID_HH",
			     NC_SHORT, 2, variable_shape,
			     &PSD_varid[1]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/* See header file for max allocation size for PSD_varid */

	break;

    case ACROBAT_CODED_SPECTRA :
	status = nc_def_var (ncid, "PSD_HH",
			     NC_SHORT, 3, variable_shape,
			     &PSD_varid[0]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
	break;

    case GALILEO_SPECTRA :
    case TEST_SPECTRA :
	status = nc_def_var (ncid, "PSD_HH",
			     NC_SHORT, 3, variable_shape,
			     &PSD_varid[PSD_HH]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	status = nc_def_var (ncid, "PSD_HV",
			     NC_SHORT, 3, variable_shape,
			     &PSD_varid[PSD_HV]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	status = nc_def_var (ncid, "PSD_VV",
			     NC_SHORT, 3, variable_shape,
			     &PSD_varid[PSD_VV]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	status = nc_def_var (ncid, "PSD_VH",
			     NC_SHORT, 3, variable_shape,
			     &PSD_varid[PSD_VH]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	variable_shape[1] = dimensions->spectra_number_dim;
	variable_shape[2] = dimensions->fft_bin_dim;
	variable_shape[3] = dimensions->range_dim;
	status = nc_def_var (ncid, "IPF_HH",
			     NC_SHORT, 4, variable_shape,
			     &PSD_varid[IPF_HH]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	status = nc_def_var (ncid, "QPF_HH",
			     NC_SHORT, 4, variable_shape,
			     &PSD_varid[QPF_HH]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	status = nc_def_var (ncid, "IPF_HV",
			     NC_SHORT, 4, variable_shape,
			     &PSD_varid[IPF_HV]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	status = nc_def_var (ncid, "QPF_HV",
			     NC_SHORT, 4, variable_shape,
			     &PSD_varid[QPF_HV]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
	break;

    case GALILEO_SPECTRA_RAPID :
    case TEST_SPECTRA_RAPID :
	variable_shape[0] = dimensions->time_dim;
	variable_shape[1] = dimensions->fft_bin_dim;
	status = nc_def_var (ncid, "PSD_RAPID_HH",
			     NC_SHORT, 2, variable_shape,
			     &PSD_varid[0]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/* See header file for max allocation size for PSD_varid */

	break;
    }
}

/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
void
RNC_WriteDynamicVariables (int                     ncid,
			   const RSP_ParamStruct * param,
			   RSP_ObservablesStruct * obs)
{
    size_t     variable_count[2];
    size_t     variable_start[2];
    int        status;
    int        n;
    float      temp_float;

    /*--------------------------------------------------------------------------*
     * write time                                                               *
     *--------------------------------------------------------------------------*/
    variable_start[0] = obs->ray_number;
    temp_float = (((int)obs->hour * 3600) + ((int)obs->minute * 60) +
		  obs->second + ((float)obs->centisecond / 100.0));
    status = nc_put_var1_float (ncid, obs->tsid, variable_start, &temp_float);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * write dish_time                                                          *
     *--------------------------------------------------------------------------*/
    temp_float = (((int)obs->dish_hour * 3600) +
		  ((int)obs->dish_minute * 60) + obs->dish_second +
		  ((float)obs->dish_centisecond / 100.0));
    status = nc_put_var1_float (ncid, obs->dish_tsid,
				variable_start, &temp_float);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * write elevation                                                          *
     *--------------------------------------------------------------------------*/
    status = nc_put_var1_float (ncid, obs->elevationid,
				variable_start, &obs->elevation);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * write azimuth                                                            *
     *--------------------------------------------------------------------------*/
    temp_float = obs->azimuth + param->azimuth_offset;
    status = nc_put_var1_float (ncid, obs->azimuthid,
				variable_start, &temp_float);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * write radar observables                                                  *
     *--------------------------------------------------------------------------*/
    variable_start[0] = obs->ray_number ;
    variable_start[1] = 0;
    variable_count[0] = 1;
    variable_count[1] = param->samples_per_pulse;

    for (n = 0; n < obs->n_obs; n++)
    {
	if (!obs->record_observable[n])
	    continue;

	status = nc_put_vara_float (ncid, obs->varid[n], variable_start,
				    variable_count, obs->data[n]);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
    }

    obs->ray_number++;
}

/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
void
RNC_WriteLogPSDVariables (int                     ncid,
			  int                     radar,
			  const RSP_ParamStruct * param,
			  RSP_ObservablesStruct * obs,
			  const PolPSDStruct *    PSD,
			  const IQStruct *        IQStruct,
			  const int *             PSD_varid)
{
    size_t     variable_count[4];
    size_t     variable_start[4];
    int     status;
    int     n,j;
    float   temp_float;
    short int *log_psd;


    variable_start[0] = obs->PSD_ray_number;
    variable_start[1] = 0;
    variable_start[2] = 0;
    variable_start[3] = 0;

    variable_count[0] = 1;
    variable_count[1] = 1;
    variable_count[2] = param->npsd;
    variable_count[3] = 0;

    /* create a short integer array to allow the results of log10 to be stored */
    log_psd = (short int *)calloc (param->npsd, sizeof (short int));
    if (log_psd == NULL)
    {
	printf ("memory request for log_psd failed\n");
	exit (1);
    }

    /*--------------------------------------------------------------------------*
     * write time                                                               *
     *--------------------------------------------------------------------------*/
    temp_float = (((int)obs->hour * 3600) + ((int)obs->minute * 60) +
		  obs->second + ((float)obs->centisecond / 100.0));
    status = nc_put_var1_float (ncid, obs->tsid, variable_start, &temp_float);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

#if 0
    /*--------------------------------------------------------------------------*
     * write dish_time                                                          *
     *--------------------------------------------------------------------------*/
    temp_float = (((int)obs->dish_hour * 3600) +
		  ((int)obs->dish_minute * 60) + obs->dish_second +
		  ((float)obs->dish_centisecond / 100.0));
    status = nc_put_var1_float (ncid, obs->dish_tsid,
				variable_start, &temp_float);
    if (status != NC_NOERR) check_netcdf_handle_error (status);
#endif

    /*--------------------------------------------------------------------------*
     * write elevation                                                          *
     *--------------------------------------------------------------------------*/
    status = nc_put_var1_float (ncid, obs->elevationid,
				variable_start, &obs->elevation);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * write azimuth                                                            *
     *--------------------------------------------------------------------------*/
    temp_float = obs->azimuth + param->azimuth_offset;
    status = nc_put_var1_float (ncid, obs->azimuthid,
				variable_start, &temp_float);
    if (status != NC_NOERR) check_netcdf_handle_error (status);

    /*--------------------------------------------------------------------------*
     * write radar observables                                                  *
     *--------------------------------------------------------------------------*/
    switch (radar)
    {
    case COPERNICUS_CODED_SPECTRA :
	if (param->include_iq_in_spectra == 1)
	{
	    /* write out IQ data first for all gates in one go */
	    variable_count[1] = param->spectra_averaged;
	    variable_count[2] = param->npsd;
	    variable_count[3] = param->samples_per_pulse;

	    printf ("ipf_hh\n");
	    /* write out IPF_HH */
	    status = nc_put_vara_short (ncid, PSD_varid[IPF_HH],
					variable_start, variable_count,
					(short *)IQStruct->I_uncoded_copolar_H);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    printf ("ipf_hv\n");
	    /* write out IPF_HV */
	    status = nc_put_vara_short (ncid, PSD_varid[IPF_HV],
					variable_start, variable_count,
					(short *)IQStruct->I_uncoded_crosspolar_H);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    variable_count[2] = param->npsd/param->num_interleave;

	    printf ("ipf_hhp\n");
	    /* write out IPF_HHP */
	    status = nc_put_vara_long (ncid, PSD_varid[IPF_HHP],
				       variable_start, variable_count,
				       IQStruct->I_coded_copolar_H);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    printf ("ipf_hvp\n");
	    /* write out IPF_HVP */
	    status = nc_put_vara_long (ncid, PSD_varid[IPF_HVP],
				       variable_start, variable_count,
				       IQStruct->I_coded_crosspolar_H);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    variable_count[2] = param->npsd;
	    printf ("qpf_hh\n");
	    /* write out QPF_HH */
	    status = nc_put_vara_short (ncid, PSD_varid[QPF_HH],
					variable_start, variable_count,
					(short *)IQStruct->Q_uncoded_copolar_H);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    printf ("qpf_hv\n");
	    /* write out QPF_HV */
	    status = nc_put_vara_short (ncid, PSD_varid[QPF_HV],
					variable_start, variable_count,
					(short *)IQStruct->Q_uncoded_crosspolar_H);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    variable_count[2] = param->npsd/param->num_interleave;

	    printf ("qpf_hhp\n");
	    /* write out QPF_HHP */
	    status = nc_put_vara_long (ncid, PSD_varid[QPF_HHP],
				       variable_start, variable_count,
				       IQStruct->Q_coded_copolar_H);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    printf ("ipf_hvp\n");
	    /* write out QPF_HVP */
	    status = nc_put_vara_long (ncid, PSD_varid[QPF_HVP],
				       variable_start, variable_count,
				       IQStruct->Q_coded_crosspolar_H);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);

	    printf ("we are done with i and q dump\n");
	    if (param->num_tx_pol == 2)
	    {
		variable_count[2] = param->npsd;
		/* write out IPF_VV */
		status = nc_put_vara_short (ncid, PSD_varid[IPF_VV],
					    variable_start, variable_count,
					    (short *)IQStruct->I_uncoded_copolar_V);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		/* write out IPF_VH */
		status = nc_put_vara_short (ncid, PSD_varid[IPF_VH],
					    variable_start, variable_count,
					    (short *)IQStruct->I_uncoded_crosspolar_V);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		variable_count[2] = param->npsd/param->num_interleave;

		/* write out IPF_VVP */
		status = nc_put_vara_long (ncid, PSD_varid[IPF_VVP],
					   variable_start, variable_count,
					   IQStruct->I_coded_copolar_V);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		/* write out IPF_VHP */
		status = nc_put_vara_long (ncid, PSD_varid[IPF_VHP],
					   variable_start, variable_count,
					   IQStruct->I_coded_crosspolar_V);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		variable_count[2] = param->npsd;

		/* write out QPF_VV */
		status = nc_put_vara_short (ncid, PSD_varid[QPF_VV],
					    variable_start, variable_count,
					    (short *)IQStruct->Q_uncoded_copolar_V);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		/* write out QPF_VH */
		status = nc_put_vara_short (ncid, PSD_varid[QPF_VH],
					    variable_start, variable_count,
					    (short *)IQStruct->Q_uncoded_crosspolar_V);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		variable_count[2] = param->npsd/param->num_interleave;

		/* write out QPF_VVP */
		status = nc_put_vara_long (ncid, PSD_varid[QPF_VVP],
					   variable_start, variable_count,
					   IQStruct->Q_coded_copolar_V);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		/* write out QPF_VHP */
		status = nc_put_vara_long (ncid, PSD_varid[QPF_VHP],
					   variable_start, variable_count,
					   IQStruct->Q_coded_crosspolar_V);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	}

	/* one gate at a time */
	variable_count[1] = 1;
	variable_count[2] = param->npsd;

	/* PSD_HH */
	for (n = 0; n < param->samples_per_pulse; n++)
	{
	    /* calculate the log10 of the PSD */
	    for (j = 0; j < param->npsd; j++)
	    {
		log_psd[j] = (short int) 1000 * log10 (PSD[n].HH[j]);
	    }
	    variable_start[1] = n;
	    status = nc_put_vara_short (ncid, PSD_varid[PSD_HH],
					variable_start, variable_count,
					log_psd);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);
	}

	/* PSD_HV */
	for (n = 0; n < param->samples_per_pulse; n++)
	{
	    /* calculate the log10 of the PSD */
	    for (j = 0; j < param->npsd; j++)
	    {
		log_psd[j] = (short int) 1000 * log10 (PSD[n].HV[j]);
	    }

	    variable_start[1] = n;
	    status = nc_put_vara_short (ncid, PSD_varid[PSD_HV],
					variable_start, variable_count,
					log_psd);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);
	}

	if (param->num_tx_pol == 2)
	{
	    /* PSD_VV */
	    for (n = 0; n < param->samples_per_pulse; n++)
	    {
		/* calculate the log10 of the PSD */
		for (j = 0; j < param->npsd; j++)
		{
		    log_psd[j] = (short int) 1000 * log10 (PSD[n].VV[j]);
		}
		variable_start[1] = n;
		status = nc_put_vara_short (ncid, PSD_varid[PSD_VV],
					    variable_start, variable_count,
					    log_psd);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }

	    /* PSD_VH */
	    for (n = 0; n < param->samples_per_pulse; n++)
	    {
		/* calculate the log10 of the PSD */
		for (j = 0; j < param->npsd; j++)
		{
		    log_psd[j] = (short int) 1000 * log10 (PSD[n].VH[j]);
		}
		variable_start[1] = n;
		status = nc_put_vara_short (ncid, PSD_varid[PSD_VH],
					    variable_start, variable_count,
					    log_psd);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	}

	/* PSD_HHP */
	for (n = 0; n < param->samples_per_pulse; n++)
	{
	    /* calculate the log10 of the PSD */
	    for (j = 0; j < param->npsd; j++)
	    {
		log_psd[j] = (short int) 1000 * log10 (PSD[n].HHP[j]);
	    }
	    variable_start[1] = n;
	    status = nc_put_vara_short (ncid, PSD_varid[PSD_HHP],
					variable_start, variable_count,
					log_psd);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);
	}

	/* PSD_HVP */
	for (n = 0; n < param->samples_per_pulse; n++)
	{
	    /* calculate the log10 of the PSD */
	    for (j = 0; j < param->npsd; j++)
	    {
		log_psd[j] = (short int) 1000 * log10 (PSD[n].HVP[j]);
	    }
	    variable_start[1] = n;
	    status = nc_put_vara_short (ncid, PSD_varid[PSD_HVP],
					variable_start, variable_count,
					log_psd);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);
	}

	if (param->num_tx_pol == 2)
	{
	    /* PSD_VVP */
	    for (n = 0; n < param->samples_per_pulse; n++)
	    {
		/* calculate the log10 of the PSD */
		for (j = 0; j < param->npsd; j++)
		{
		    log_psd[j] = (short int) 1000 * log10 (PSD[n].VVP[j]);
		}
		variable_start[1] = n;
		status = nc_put_vara_short (ncid, PSD_varid[PSD_VVP],
					    variable_start, variable_count,
					    log_psd);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }

	    /* PSD_VHP */
	    for (n = 0; n < param->samples_per_pulse; n++)
	    {
		/* calculate the log10 of the PSD */
		for (j = 0; j < param->npsd; j++)
		{
		    log_psd[j] = (short int) 1000 * log10 (PSD[n].VHP[j]);
		}
		variable_start[1] = n;
		status = nc_put_vara_short (ncid, PSD_varid[PSD_VHP],
					    variable_start, variable_count,
					    log_psd);
		if (status != NC_NOERR) check_netcdf_handle_error (status);
	    }
	}
	break;

    case ACROBAT_CODED_SPECTRA :
	for (n = 0; n < param->samples_per_pulse; n++)
	{
	    /* calculate the log10 of the PSD */
	    for (j = 0; j < param->npsd; j++)
	    {
		log_psd[j] = (short int) 1000 * log10 (PSD[n].HH[j]);
	    }
	    variable_start[1] = n;
	    status = nc_put_vara_short (ncid, PSD_varid[0],
					variable_start, variable_count,
					log_psd);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);
	}
	break;

    case GALILEO_SPECTRA :
    case TEST_SPECTRA :
	/* write out IQ data first for all gates in one go */
	variable_start[0] = obs->PSD_ray_number;
	variable_start[1] = 0;
	variable_start[2] = 0;
	variable_start[3] = 0;
	variable_count[0] = 1;
	variable_count[1] = param->spectra_averaged;
	variable_count[2] = param->npsd;
	variable_count[3] = param->samples_per_pulse;

	/* write out IPF_HH */
	status = nc_put_vara_short (ncid, PSD_varid[IPF_HH],
				    variable_start, variable_count,
				    (short *)IQStruct->I_uncoded_copolar_H);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/* write out QPF_HH */
	status = nc_put_vara_short (ncid, PSD_varid[QPF_HH],
				    variable_start, variable_count,
				    (short *)IQStruct->Q_uncoded_copolar_H);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/* one gate at a time */
	variable_count[0] = 1;
	variable_count[1] = 1;
	variable_count[2] = param->npsd;

	/* PSD_HH */
	for (n = 0; n < param->samples_per_pulse; n++)
	{
	    /* calculate the log10 of the PSD */
	    for (j = 0; j < param->npsd; j++)
	    {
		log_psd[j] = (short int) 1000 * log10 (PSD[n].HH[j]);
	    }
	    variable_start[1] = n;
	    status = nc_put_vara_short (ncid, PSD_varid[PSD_HH],
					variable_start, variable_count,
					log_psd);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);
	}

	status = nc_sync (ncid);
	if (status != NC_NOERR) check_netcdf_handle_error (status);

	/* PSD_HV */
	for (n = 0; n < param->samples_per_pulse; n++)
	{
	    /* calculate the log10 of the PSD */
	    for (j = 0; j < param->npsd; j++)
	    {
		log_psd[j] = (short int) 1000 * log10 (PSD[n].HV[j]);
	    }
	    variable_start[1] = n;
	    status = nc_put_vara_short (ncid, PSD_varid[PSD_HV],
					variable_start, variable_count,
					log_psd);
	    if (status != NC_NOERR) check_netcdf_handle_error (status);
	}
	status = nc_sync (ncid);
	if (status != NC_NOERR) check_netcdf_handle_error (status);
	break;
    }

    free (log_psd);

    obs->PSD_ray_number++;

    status = nc_sync (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);
}

/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
void
RNC_WriteRapidLogPSDVariables (int                     ncid,
			       int                     radar,
			       const RSP_ParamStruct * param,
			       RSP_ObservablesStruct * obs,
			       const PolPSDStruct *    PSD,
			       const int *             PSD_varid)
{
    size_t     variable_count[3];
    size_t     variable_start[3];
    int     status;
    int     n,j;
    float   timestamp;
    float  azimuth;
    short int       *log_psd;
    short int *log_psd_coded;
    int temp_int = 0;
    int save_data = 0;
    int saved_ray = 0;
    int start_bin = 0;
    int noise_bin = 0;

    variable_count[0] = 1;
    variable_count[1] = param->npsd;

    /* create a short integer array to allow the results of log10 to be stored */
    log_psd = (short int *)calloc (param->npsd, sizeof (short int));
    if (log_psd == NULL)
    {
	printf ("memory request for log_psd failed\n");
	exit (1);
    }

    /* create a short integer array to allow the results of log10 to be stored */
    log_psd_coded = (short int *)calloc (param->npsd, sizeof (short int));
    if (log_psd_coded == NULL)
    {
	printf ("memory request for log_psd_coded failed\n");
	exit (1);
    }

    /* calculate time */
    timestamp = (((int)obs->hour * 3600) + ((int)obs->minute * 60) +
		 obs->second + ((float)obs->centisecond / 100.0));

    /* obtain azimuth */
    azimuth = obs->azimuth + param->azimuth_offset;

    /* start point of copernicus fft */
    if (radar == GALILEO_SPECTRA_RAPID || radar == TEST_SPECTRA_RAPID)
    {
	temp_int = 0;
	variable_count[1] = param->npsd;
	start_bin = 6;
	noise_bin = 0;
    }
    else if (radar == COPERNICUS_SPECTRA_RAPID)
    {
	temp_int = 0;
	/* since I am dumping out the spectra for a range bin */
	/* can I make the */
	variable_count[1] = param->npsd;
	start_bin = 26;
	noise_bin = 469;
    }

    switch (radar)
    {
    case COPERNICUS_SPECTRA_RAPID :
    case GALILEO_SPECTRA_RAPID :
    case TEST_SPECTRA_RAPID :
	for (n = 0; n < param->samples_per_pulse; n++)
	{
	    save_data = 0;
	    if (radar == GALILEO_SPECTRA_RAPID || radar == TEST_SPECTRA_RAPID)
	    {
		for (j = 0; j < param->npsd; j++)
		{
		    log_psd[j] = (short int) 1000 * log10 (PSD[n].HH[j]);
		    if (log_psd[j] > 5000)
		    {
			if (n >= start_bin)
			{
			    /* this is an attempt to jump over the central bin */
			    if (param->fft_bins_interpolated > 0)
			    {
				if ((j < (param->npsd/2 - (param->fft_bins_interpolated - 1))) ||
				    (j > (param->npsd/2 + (param->fft_bins_interpolated - 1))))
				{
				    save_data = 1;
				}
			    }
			    else
			    {
				save_data = 1;
			    }
			}
		    }
		}
		if ((n == noise_bin) || (n == noise_bin + 1) ||
		    (n == noise_bin + 2) || (n == noise_bin + 3))
		{
		    save_data = 1;
		}
	    }
	    else if (radar == COPERNICUS_SPECTRA_RAPID)
	    {
		for (j = 0; j < param->npsd; j++)
		{
		    log_psd_coded[j] = (short int) 1000 * log10 (PSD[n].HHP[j]);
		    log_psd[j] = (short int) 1000 * log10 (PSD[n].HH[j]);
		    if (log_psd_coded[j] > 4700)
		    {
			if (n >= start_bin)
			{
			    save_data = 1;
			}
		    }
		}
		if ((n == noise_bin + 0) || (n == noise_bin + 1) ||
		    (n == noise_bin + 2) || (n == noise_bin + 3))
		{
		    save_data = 1;
		}
	    }
	    if (save_data == 1)
	    {
		variable_start[0] = obs->bin_ray_number;
		variable_start[1] = 0;

		/* write out variables */
		/* time */
		status = nc_put_var1_float (ncid, obs->tsid,
					    variable_start, &timestamp);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		/* elevation */
		status = nc_put_var1_float (ncid, obs->elevationid,
					    variable_start, &obs->elevation);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		/* azimuth */
		status = nc_put_var1_float (ncid, obs->azimuthid,
					    variable_start, &azimuth);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		/* ray_number */
		status = nc_put_var1_int (ncid, obs->ray_numberid,
					  variable_start, &obs->ray_number);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		/* bin_number */
		status = nc_put_var1_int (ncid, obs->bin_numberid,
					  variable_start, &n);
		if (status != NC_NOERR) check_netcdf_handle_error (status);

		/* write out fft for this ray */
		if (radar == GALILEO_SPECTRA_RAPID || radar == TEST_SPECTRA_RAPID)
		{
		    status = nc_put_vara_short (ncid, PSD_varid[0],
						variable_start, variable_count,
						&log_psd[temp_int]);
		    if (status != NC_NOERR) check_netcdf_handle_error (status);
		}
		else if (radar == COPERNICUS_SPECTRA_RAPID)
		{
		    status = nc_put_vara_short (ncid, PSD_varid[0],
						variable_start, variable_count,
						&log_psd_coded[temp_int]);
		    if (status != NC_NOERR) check_netcdf_handle_error (status);

		    status = nc_put_vara_short (ncid, PSD_varid[1],
						variable_start, variable_count,
						&log_psd[temp_int]);
		    if (status != NC_NOERR) check_netcdf_handle_error (status);
		}
		obs->bin_ray_number = obs->bin_ray_number + 1;
		saved_ray = 1;
	    }
	}
	break;
    }

    /* if we have stored some PSD of a ray let us inc. the ray_number */
    if (saved_ray == 1)
    {
	obs->ray_number++;
    }

    free (log_psd);
    free (log_psd_coded);

    status = nc_sync (ncid);
    if (status != NC_NOERR) check_netcdf_handle_error (status);
}

/*****************************************************************************
 *                                                                           *
 *****************************************************************************/
void
check_netcdf_handle_error (int status)
{
    if (status != NC_NOERR)
    {
	fprintf (stderr, "netCDF ERROR : %s\n", nc_strerror (status));
	exit (1);
    }
}
