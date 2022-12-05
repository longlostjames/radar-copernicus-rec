// RSP_Correlate.c
// -----------------
// Part of the Chilbolton Radar Signal Processing Package
//
// Purpose: Manages RSP_Observables structures
//
// Created on: 11/02/04
// Created by: E. Pavelin (EGP)
// --------------------------------------------------------
//
// Revision History
// --------------------------------------------------------
// When     : Who : What
// --------------------------------------------------------
// 11-02-04 : EGP : Created
// 18-02-04 : EGP : RSP_ObsNew now returns the new pointer

#include <string.h>
#include <stdlib.h>

#include <RSP.h>

// Initialise observables structure
void
RSP_ObsInit (RSP_ObservablesStruct * obs)
{
    obs->n_obs          = 0;
    obs->ray_number     = 0;
    obs->PSD_ray_number = 0;
}


// Initialise a new observable
float *
RSP_ObsNew (RSP_ObservablesStruct * obs,
	    const char *            name,
	    int                     n_elements,
	    int                     record_observable)
{
    float * data;

    data = malloc (sizeof (float) * n_elements);
    if (data != NULL)
    {
	strcpy (obs->name [obs->n_obs], name);
	obs->data              [obs->n_obs] = data;
	obs->n_elements        [obs->n_obs] = n_elements;
	obs->record_observable [obs->n_obs] = record_observable;
	obs->n_obs++;
    }
    return data;
}

// Free all the memory
void
RSP_ObsFree (RSP_ObservablesStruct * obs)
{
    int i;

    for (i = 0; i < obs->n_obs; i++)
    {
	free (obs->data [i]);
    }
    obs->n_obs = 0;
}


// Return pointer corresponding to a given name
float *
RSP_ObsGet (const RSP_ObservablesStruct * obs,
	    const char *                  name)
{
    int i;

    for (i = 0; i < obs->n_obs; i++)
    {
	if (!strcmp (obs->name [i], name))
	    return obs->data [i];
    }

    return (NULL);
}

// Return index corresponding to a given name
int
RSP_ObsIndex (const RSP_ObservablesStruct * obs,
	      const char *                  name)
{
    int i;

    for (i = 0; i < obs->n_obs; i++)
    {
	if (!strcmp (obs->name [i], name))
	    return i;
    }
    return -1;
}
