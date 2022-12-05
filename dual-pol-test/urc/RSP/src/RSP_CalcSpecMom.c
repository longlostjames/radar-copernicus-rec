// RSP_CalcSpecMom.c
// -----------------
// Part of the Chilbolton Radar Signal Processing Package
//
// Purpose: To estimate spectral moments from a radar power
//          sepctrum.
//
// Created on: 05/11/03
// Created by: E. Pavelin (EGP)
// --------------------------------------------------------
//
// Revision History
// --------------------------------------------------------
// When     : Who : What
// --------------------------------------------------------
// 05-11-03 : EGP : Created
// 03-03-04 : EGP : Algorithm streamlined
// 28-11-11 : JCN : Added 3rd and 4th order moments

#include <stdio.h>
#include <math.h>

#include <RSP.h>

int
RSP_CalcSpecMom (const float *          psd,
		 int                    nBins,
		 const RSP_PeakStruct * peak,
		 float                  noiseLevel,
		 float *                moments,
		 size_t                 num_moments)
{

    register int k, j;
    int b2;

    if (num_moments < 3)
	return -1;

    if (peak->leftBin > peak->rightBin) // Folded spectrum
    {
	b2 = peak->rightBin + nBins;
    }
    else
    {
	b2 = peak->rightBin;
    }

    moments[0] = 0.0f;
    moments[1] = 0.0f;
    moments[2] = 0.0f;

    for (k = peak->leftBin; k <= b2; k++)
    {
	float pxx;

	j   = k % nBins;
	pxx = psd[j] - noiseLevel;
	moments[0] += pxx;
	moments[1] += k * pxx;
	moments[2] += k * k * pxx;
    }

    moments[1] = moments[1] / moments[0];
    moments[2] = sqrt (moments[2] / moments[0] - moments[1] * moments[1]);

    // Deal with case of folding from -ve to +ve frequency
    if (peak->leftBin > peak->rightBin && peak->peakBin < peak->leftBin)
    {
	moments[1] -= nBins;
    }

    if (num_moments < 5)
	return 0;

    moments[3] = 0.0f;
    moments[4] = 0.0f;

    for (k = peak->leftBin; k <= b2; k++)
    {
	float pxx;

        j   = k % nBins;
        pxx = psd[j] - noiseLevel;
        moments[3] += pxx * (k - moments[1]) * (k - moments[1]) * (k - moments[1]);
        moments[4] += pxx * (k - moments[1]) * (k - moments[1]) * (k - moments[1]) * (k - moments[1]);
    }

    moments[3] = (moments[3] / moments[0]);
    moments[3] = moments[3] / (moments[2] * moments[2] * moments[2]);
    moments[4] = (moments[4] / moments[0]);
    moments[4] = moments[4] / (moments[2] * moments[2] * moments[2] * moments[2]);
    moments[4] = moments[4] - 3.0;

    return 0;
}


//----------------------------------------------------------
// Converts bin number (floating point) to Doppler velocity
//----------------------------------------------------------
float
RSP_BinToVelocity (float                   bin,
		   const RSP_ParamStruct * param)
{
    float mid;

    mid = param->nfft / 2 - 1;
    return (param->frequency_bin_width * (bin - mid) / param->hz_per_mps);
}


// Calculate total noise power under a peak
float
RSP_CalcNoisePower (float                   noiseLevel,
		    const RSP_PeakStruct *  peak,
		    const RSP_ParamStruct * param)
{
    float noise_power;

    noise_power = (peak->rightBin - peak->leftBin + 1);
    if (peak->leftBin > peak->rightBin) // If folded
    {
	noise_power = noise_power + param->npsd;
    }
    noise_power = noise_power * param->frequency_bin_width * noiseLevel;

    return noise_power;
}
