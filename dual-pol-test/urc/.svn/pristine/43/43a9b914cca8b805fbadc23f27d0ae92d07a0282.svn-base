// RSP_FreeMemory.c
// --------------
// Part of the Chilbolton Radar Signal Processing Package
//
// Purpose: To free the memory blocks previously allocated
//          by the RSP_InitialiseParams module.
//
// Created on: 19/01/04
// Created by: E. Pavelin (EGP)
// -------------------------------------------------------

#include <stdlib.h>

#include <RSP.h>

void
RSP_FreeMemory (RSP_ParamStruct * param)
{
    // Free memory blocks
    free (param->range);          // also frees parm->range_correction
    free (param->frequency_axis); // also frees parm->velocity_axis
    free (param->window);
}
