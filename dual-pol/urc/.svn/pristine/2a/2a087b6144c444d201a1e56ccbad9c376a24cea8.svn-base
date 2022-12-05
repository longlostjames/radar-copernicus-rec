#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <RNC.h>

double
RNC_GetConfigDouble(const char *filename,
		    const char *keyword)
{
    /*-----------------------------------------------------------------------*
     * Reads a double from a space-delimited config file                     *
     * Reads the double preceded by the given keyword in the given filename  *
     *-----------------------------------------------------------------------*/
    char   valuestr [80];
    double value;
    char * end;

    if (RNC_GetConfig (filename, keyword, valuestr, sizeof (valuestr)))
    {
	printf ("** getconf_double: Keyword %s not found.\n", keyword);
	return 0.0;
    }

    value = strtod (valuestr, &end);
    if (end == valuestr)
    {
	printf ("** getconf_double: Bad value for Keyword %s.\n", keyword);
	return 0.0;
    }

    return value;
}

float
RNC_GetConfigFloat(const char *filename,
		   const char *keyword)
{
    /*-----------------------------------------------------------------------*
     * Reads a float from a space-delimited config file                      *
     * Reads the float preceded by the given keyword in the given filename   *
     *-----------------------------------------------------------------------*/
    char   valuestr [80];
    float  value;
    char * end;

    if (RNC_GetConfig (filename, keyword, valuestr, sizeof (valuestr)))
    {
	printf ("** getconf_float: Keyword %s not found.\n", keyword);
	return 0.0;
    }
 
    value = strtod (valuestr, &end);
    if (end == valuestr)
    {
	printf ("** getconf_float: Bad value for Keyword %s.\n", keyword);
	return 0.0;
    }

    return value;
}


int
RNC_GetConfig(const char *filename,
	 const char * keyword,
	 char *       value,
	 size_t       size)
{
    /*-----------------------------------------------------------------------*
     * Reads a string from a space-delimited config file                     *
     * Reads the string preceded by the given keyword in the given filename  *
     *-----------------------------------------------------------------------*/
    char   buffer[255];
    char * p;
    FILE * fptr;
    int    status = 1;

    fptr = fopen (filename, "r");
    if (fptr == NULL)
    {
	printf ("** getconf: Unable to open %s: %m\n", filename);
	return 1;
    }

    while (feof (fptr) != EOF)
    {
	if (fgets (buffer, 255, fptr) == NULL)
	{
	    if (feof (fptr) != EOF)
		printf ("** getconf: Error reading file.\n");
	    return 1;
	}

	if (*buffer == '#')
	    continue;

	p = strstr (buffer, keyword);
	if (p == buffer)
	{
	    p += strlen (keyword);

	    /* Ensure that the delimiter has been reached */
	    if (*p != ' ' && *p != '\t' && *p != '\n' && *p != '\r')
		continue;

	    /* Skip delimeters to data/end of line */
	    while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') p++;

	    strncpy (value, p, size);
	    value [size - 1] = '\0';
	    status = 0;
	    break;
	}
    }

    fclose (fptr);
    return status;
}
