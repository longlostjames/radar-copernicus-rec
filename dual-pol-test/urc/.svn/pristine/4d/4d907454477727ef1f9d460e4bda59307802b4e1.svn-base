// RSP_DisplayParams.c
// -------------------
// Part of the Chilbolton Radar Signal Processing Package
//
// Purpose: To display the contents of the param structure
//
// Created on: 19/02/04
// Created by: E. Pavelin (EGP)
// -------------------------------------------------------

#include <stdio.h>

#include <RSP.h>

void
RSP_DisplayParams (const RSP_ParamStruct * param)
{
    printf ("\n");
    printf ("=============================================\n");
    printf ("RADAR PARAMETERS\n");
    printf ("---------------------------------------------\n");
    printf ("Radar Frequency            : %f GHz\n",  param->frequency);
    printf ("Pulse repetition frequency : %f Hz\n",   param->prf);
    printf ("Transmit power             : %f W\n",    param->transmit_power);
    printf ("Inter-pulse period         : %f msec\n", param->prt * 1000.0);
    printf ("Pulse/subpulse length      : %f usec\n", param->pulse_period);
    printf ("Range resolution           : %f m\n",    param->range_resolution);
    printf ("Maximum unambiguous range  : %f km\n",   param->prt * SPEED_LIGHT / 2000.0);

    if (!param->long_pulse_mode)
    {
	printf ("Long pulse mode            : Off\n");
	if (param->alternate_modes)
	{
	    printf ("Alternating modes           : On\n");
	    printf ("Mode0                       : %d\n", param->mode0);
	    printf ("Mode1                       : %d\n", param->mode1);
	    printf ("Rays in Mode0               : %d\n", param->nrays_mode0);
	    printf ("Rays in Mode1               : %d\n", param->nrays_mode1);
	}
	else
	{
	    printf ("Alternating modes           : Off\n");
	    printf ("Mode                        : %d\n", param->mode0);
	}
    }
    else
    {
	printf ("Long pulse mode            : On\n");
    }

    printf ("Pulse coding               : ");
    if (param->code_length > 1)
    {
	printf ("%d bits\n", param->code_length);
	printf ("Complementary codes        : ");
	if (param->number_of_codes > 1)
	{
	    printf ("%d\n", param->number_of_codes);
	}
	else
	{
	    printf ("OFF\n");
	}
    }
    else
    {
	printf ("OFF\n");
    }
    printf ("\n");
    printf ("DATA ACQUISITION PARAMETERS\n");
    printf ("---------------------------------------------\n");
    printf ("Pulses per DAQ cycle       : %d\n",     param->pulses_per_daq_cycle);
    printf ("DAQ cycle time             : %f sec\n", param->daq_time);
    printf ("Number of range gates      : %d\n",     param->samples_per_pulse);
    printf ("Range gate width           : %f m\n",   param->range_gate_width);
    printf ("Maximum sampled range      : %f km\n",  param->range[param->samples_per_pulse - 1] / 1000.0);
    printf ("DAQ clock input            : %f MHz\n", param->clock/1.0e6);
    printf ("DAQ divfactor              : %d\n",     param->clock_divfactor);
    printf ("DAQ sample frequency       : %f MHz\n", param->sample_frequency / 1.0e6);
    printf ("Pulse width oversampling   : %d x\n",   param->oversample_ratio);
    printf ("DAQ delay clock cycles     : %d\n",     param->delay_clocks);
    printf ("DAQ channels in use        : %d\n",     param->ADC_channels);
    printf ("\n");
    printf ("SIGNAL PROCESSING PARAMETERS\n");
    printf ("---------------------------------------------\n");
    printf ("Pulses coherently averaged : %d\n",      param->pulses_coherently_averaged);
    printf ("Spectra averaged           : %d\n",      param->spectra_averaged);
    printf ("Total dwell time           : %f sec\n",  param->dwell_time);
    printf ("Number of points in FFT    : %d\n",      param->nfft);
    printf ("Number of bins in PSD      : %d\n",      param->npsd);
    printf ("FFT frequency bin width    : %f Hz\n",   param->frequency_bin_width);
    printf ("Doppler shift per m/s      : %f Hz\n",   param->hz_per_mps);
    printf ("Folding velocity           : %f m/s\n",  param->folding_velocity);
    printf ("Folding frequency          : %f Hz\n",   param->folding_frequency);
    printf ("Clutter interpolation      : %d bins\n", param->fft_bins_interpolated);
    printf ("Number of peaks processed  : %d\n",      param->num_peaks);
    printf ("\n");
    printf ("CALIBRATION PARAMETERS\n");
    printf ("---------------------------------------------\n");
    printf ("Range offset added         : %f m\n",       param->range_offset);
    printf ("Beam azimuth offset        : %f degrees\n", param->azimuth_offset);
    printf ("Z calibration in .cal file : %f dB\n",      param->ZED_calibration_offset);
    printf ("=============================================\n");
    printf ("\n");
}
