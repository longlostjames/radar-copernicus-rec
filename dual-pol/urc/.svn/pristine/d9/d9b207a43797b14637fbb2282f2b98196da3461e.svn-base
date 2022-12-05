// RSP_ClutterInterp.c
// -----------------
// Part of the Chilbolton Radar Signal Processing Package
//
// Purpose: Interpolates over clutter near zero Doppler
//
// Created on: 21/11/03
// Created by: E. Pavelin (EGP)
// --------------------------------------------------------
//
// Revision History
// --------------------------------------------------------
// When     : Who : What
// --------------------------------------------------------
// 21-11-03 : EGP : Created

#include <RSP.h>

void
RSP_ClutterInterp (float * PSD,
		   int     nBins,
		   int     nInterp)
{
    // Interpolates nInterp bins either size of zero Doppler.
    // nInterp == 0 : Don't interpolate at all
    // nInterp == 1 : Just interpolate over central bin
    // nInterp == 2 : Interpolate over central bin +/- 1 bin  (3 bins total)
    // nInterp == 3 : Interpolate over central bin +/- 2 bins (5 bins total)
    // etc...

    int i;
    float dpdb;
    int ib1, ib2;

    if (nInterp < 1)
	return;

    ib1 = (nBins - 1) / 2 - nInterp;
    ib2 = (nBins - 1) / 2 + nInterp;

    // Calculate gradient
    dpdb = (PSD[ib2] - PSD[ib1]) / (float)(ib2 - ib1);
    for (i = ib1 + 1; i < ib2; i++)
    {
	PSD[i] = PSD[ib1] + dpdb * (i - ib1);
    }
}
