// RSP_CalcPhase.c
// --------------
// Part of the Chilbolton Radar Signal Processing Package
//
// Purpose: To calculate phase parameters from I and Q time
//          series.
//
// Created on: 12/02/04
// Created by: E. Pavelin (EGP)
// --------------------------------------------------------

#include <stdio.h>
#include <math.h>

#include <RSP.h>

// Calculates phase and std. dev. of phase in degrees

void
RSP_CalcPhase (const RSP_ComplexType * IQ,
	       float *                 phi,
	       float *                 sdphi,
	       int                     nfft)
{
    int j;
    float sum_I = 0.0f, sum_Q = 0.0f;
    float phase, diff;
    float sumsdphi;
    float phirad;

    // Calculate mean phase
    for (j = 0; j < nfft; j++)
    {
	sum_I += IQ[j].real;
	sum_Q += IQ[j].imag;
    }

    phirad = atan2 (sum_I, sum_Q);

    // CALC STD DEVIATION
    sumsdphi = 0.0;
    for (j = 0; j < nfft; j++)
    {
	phase = atan2 (IQ[j].real, IQ[j].imag);
	diff  = phase - phirad;
	if      (diff < -PI) diff += (2.0 * PI);
	else if (diff >  PI) diff -= (2.0 * PI);
	sumsdphi += (diff * diff);
    }

    *sdphi = sqrt (sumsdphi / (double)nfft) * RAD2DEG;
    *phi   = phirad * RAD2DEG;
}
