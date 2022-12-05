// RSP_Initialise.c
// --------------
// Part of the Chilbolton Radar Signal Processing Package
//
// Purpose: To initialise operating parameters and allocate
//          memory required by all modules.
//
// Created on: 17/11/03
// Created by: E. Pavelin (EGP)
// -------------------------------------------------------

#include <stdlib.h>
#include <math.h>

#include <RSP.h>

void
RSP_InitialiseParams (RSP_ParamStruct * param)
{
    int i;
    int tot_avg;

    tot_avg = param->pulses_coherently_averaged * param->number_of_codes;

    // Calculate radar operating parameters
    // See RSP.h for explanations
    param->prt              = 1.0 / param->prf;
    param->sample_frequency = param->clock / param->clock_divfactor;
    param->sample_period    = 1.0 / param->sample_frequency;
    param->daq_time         = param->pulses_per_daq_cycle * param->prt;
    param->dwell_time       = param->daq_time * param->spectra_averaged;
    param->range_gate_width = param->sample_period * SPEED_LIGHT / 2.0;
    param->range_resolution = param->pulse_period * 1e-9 * SPEED_LIGHT / 2.0;

    param->nfft                = param->pulses_per_daq_cycle / tot_avg / param->num_interleave / param->num_tx_pol;
    param->npsd                = param->nfft;
    param->frequency_bin_width = 1.0 / ((float)param->nfft * param->prt * (float)tot_avg);
    param->hz_per_mps          = 2.0 * param->frequency * 1e9 / SPEED_LIGHT;
    param->folding_frequency   = param->prf / tot_avg / 2.0;
    param->folding_velocity    = param->folding_frequency / param->hz_per_mps;
    param->oversample_ratio    = (int) (param->pulse_period * 1e-9 * param->sample_frequency);

    // Allocate memory blocks. Allocate as single blocks where size is the same.
    param->range            = malloc ((param->samples_per_pulse << 1) * sizeof (float));
    param->range_correction = param->range + param->samples_per_pulse;
    param->frequency_axis   = malloc ((param->npsd  << 1) * sizeof (float));
    param->velocity_axis    = param->frequency_axis + param->npsd;
    param->window           = malloc (param->nfft * sizeof (float));


    // Fill array variables
    for (i = 0; i < param->npsd; i++)
    {
	param->frequency_axis[i] = param->frequency_bin_width * (i - (param->nfft / 2 - 1));
	param->velocity_axis[i] = param->frequency_axis[i] / param->hz_per_mps;
    }


    for (i = 0; i < param->samples_per_pulse; i++)
    {
	// Original 	param->range[i] = i * param->range_gate_width + param->range_offset;
	// New	 	param->range[i] = (i + param->delay_clocks/param->clock_divfactor) * param->range_gate_width + param->range_offset;

  	param->range            [i] = (i + param->delay_clocks / param->clock_divfactor) * param->range_gate_width + param->range_offset;
  	param->range_correction [i]= (param->range [i] / 1000.0) * (param->range[i] / 1000.0);
    }

    param->Wss = 0;
    for (i = 0; i < param->nfft; i++)
    {
	param->window[i] = BLACKMAN_WINDOW (i, param->nfft);
	param->Wss      += param->window[i] * param->window[i];
    }
    param->Wss /= param->nfft;
}
