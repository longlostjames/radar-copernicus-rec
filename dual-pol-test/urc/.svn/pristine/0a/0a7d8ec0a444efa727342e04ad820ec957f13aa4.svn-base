// RSP_CalcPSD.c
// --------------
// Part of the Chilbolton Radar Signal Processing Package
//
// Purpose: To calculate the power spectral density (PSD) from
//          a complex video time series (I and Q).
//
// NB: RSP_CalcPSD acts destructively on the input data IQ.
//
// Created on: 03/11/03
// Created by: E. Pavelin (EGP)
// --------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
//#define FFTW_NO_Complex
#include <fftw3.h>

#include <RSP.h>

void
RSP_CalcPSD (RSP_ComplexType * IQ,
	     int               nfft,
	     const float *     window,
	     float *           psd,
	     float             norm)
{
    register int i;

    // Multiply by window
    for (i = 0; i < nfft; i++)
    {
	IQ[i].real = IQ[i].real * window[i];
	IQ[i].imag = IQ[i].imag * window[i];
    }

    // Do FFT
    RSP_FFT ((float *)IQ, nfft, 1);

    // Calculate PSD
    RSP_FFT2PowerSpec ((float *)IQ, psd, nfft, norm);
}

void
RSP_CalcPSD_FFTW (fftw_complex *  in,
		  int             nfft,
		  const fftw_plan p,
		  const float *   window,
		  float *         psd,
		  float           norm)
{
    register int i;

    // Multiply by window
    for (i = 0; i < nfft; i++)
    {
#ifdef fftw_is_complex
	in[i] *= window[i];
#else
	fftw_real (in[i]) *= window[i];
	fftw_imag (in[i]) *= window[i];
#endif
    }

    // Do FFT
    fftw_execute (p);

    // Calculate PSD
    RSP_FFT2PowerSpec_FFTW (in, psd, nfft, norm);
}

void
RSP_SubtractOffset (RSP_ComplexType * IQ,
		    int               nfft)
{
    register int i;
    float Imean = 0.0f, Qmean = 0.0f;

    // Calculate DC offsets
    for (i = 0; i < nfft; i++)
    {
	Imean += IQ[i].real;
	Qmean += IQ[i].imag;
    }
    Imean /= nfft;
    Qmean /= nfft;

    // Subtract offsets
    for (i = 0; i < nfft; i++)
    {
	IQ[i].real = (IQ[i].real - Imean);
	IQ[i].imag = (IQ[i].imag - Qmean);
    }
}

void
RSP_SubtractOffset_FFTW (fftw_complex * IQ,
			 int            nfft)
{
    register int i;
    float Imean = 0.0f, Qmean = 0.0f;

    // Calculate DC offsets
    for (i = 0; i < nfft; i++)
    {
	Imean += fftw_real (IQ[i]);
	Qmean += fftw_imag (IQ[i]);
    }
    Imean /= nfft;
    Qmean /= nfft;

    // Subtract offsets
    for (i = 0; i < nfft; i++)
    {
	fftw_real_lv (IQ[i]) -= Imean;
	fftw_imag_lv (IQ[i]) -= Qmean;
    }
}

//-------------------------------------------------------------------
/* FFT routine adapted from Numerical Recipes in C */

#define SWAP(a,b) tempr = (a); (a) = (b); (b) = tempr
void
RSP_FFT (float *       data_in,
	 unsigned long nn,
	 int           isign)
{
    register unsigned long i, m;
    unsigned long n, mmax, j, istep;
    double wtemp, wr, wpr, wpi, wi, theta;
    float tempr, tempi;
    float * data;

    data = data_in - 1; // Shift to zero-based array

    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2)
    {
        if (j > i)
	{
            SWAP (data[j], data[i]);
            SWAP (data[j+1], data[i+1]);
        }
        m = n >> 1;
        while (m >= 2 && j > m)
	{
            j  -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax = 2;
    while (n > mmax)
    {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin (0.5 * theta);
        wpr   = -2.0 * wtemp * wtemp;
        wpi   = sin (theta);
        wr    = 1.0;
        wi    = 0.0;
        for (m = 1; m < mmax; m += 2)
	{
            for (i = m; i <= n; i += istep)
	    {
                j          = i + mmax;
                tempr      = wr * data[j]   - wi * data[j+1];
                tempi      = wr * data[j+1] + wi * data[j];
                data[j]    = data[i]   - tempr;
                data[j+1]  = data[i+1] - tempi;
                data[i]   += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}
#undef SWAP

//----------------------------------------------------------------------
void
RSP_FFT2PowerSpec (const float * data,
		   float *       PSD,
		   int           nfft,
		   float         norm)
/* Calculates power spectrum from complex FFT output */
{
    int i, j, nn, i2;
    float * PP;

    PP = malloc (sizeof (float) * nfft);

    for (i = 0; i < nfft; i++)
    {
        i2 = i * 2;
        PP[i] = data[i2] * data[i2] + data[i2+1] * data[i2+1];
    }

    j  = 0;
    nn = nfft >> 1;
    for (i = nn + 1; i < nfft; i++)
    {
        PSD[j] = PP[i] * norm;
        j++;
    }
    for (i = 0; i <= nn; i++)
    {
        PSD[j] = PP[i] * norm;
        j++;
    }
    free (PP);
}

//----------------------------------------------------------------------
void
RSP_FFT2PowerSpec_FFTW (fftw_complex * data,
			float *        PSD,
			int            nfft,
			float          norm)
/* Calculates power spectrum from complex FFT output */
{
    int i, j, nn;
    float * PP;

    PP = malloc (sizeof (float) * nfft);

    for (i = 0; i < nfft; i++)
    {
        PP[i] = (fftw_real (data[i]) * fftw_real (data[i]) +
		 fftw_imag (data[i]) * fftw_imag (data[i]));
    }

    j  = 0;
    nn = nfft >> 1;
    for (i = nn + 1; i < nfft; i++)
    {
        PSD[j] = PP[i] * norm;
        j++;
    }
    for (i = 0; i <= nn; i++)
    {
        PSD[j] = PP[i] * norm;
        j++;
    }
    free (PP);
}
