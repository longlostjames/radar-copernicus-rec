// RSP_Correlate.c
// -----------------
// Part of the Chilbolton Radar Signal Processing Package
//
// Purpose: Performs integer correlation to decode binary
//          phase codes.
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

#include <stdint.h>

#include <RSP.h>

void
RSP_Correlate (const uint16_t * data,
	       const short *    code,
	       int              samples,
	       int              bits,
	       long int *       corr)
{
    int nmax;
    int j , i;

    nmax = samples - bits + 1;

    for (j = 0; j < nmax; j++)
	corr [j] = 0;

    for (j = 0; j < bits; j++)
    {
        if (code[j] == -1)
	{
            for (i = 0; i < nmax; i++)
	    {
                corr[i] -= (int)data[i+j];
            }
        }
	else
	{
            for (i = 0; i < nmax; i++)
	    {
                corr[i] += (int)data[i+j];
	    }
        }
    }
}

void
RSP_Oversample (const short * code,
		short *       newcode,
		int           numel,
		int           n)
/* Oversamples array by n times */
{
    short int i,j,ni;

    for (i = 0; i < numel; i++)
    {
        ni = n * i;
        for (j = 0; j < n; j++)
	{
            newcode [ni+j] = code [i];
        }
    }
}
