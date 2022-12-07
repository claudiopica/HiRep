/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file error.c
 * @brief Error handling functions
 */

#include <stdlib.h>
#include <stdio.h>
#include "logger.h"
#include "setup.h"

/**
 * @brief Print message to error file defined on startup.
 *
 * @param test              Condition on whether an error should be raised.
 *                          0 for no error and continue
 *                          1 for error, stop and print error message
 * @param no                Exit Code
 *                          Value smaller than zero exits immediately with code 0.
 *                          Value larger or equal then zero exits with code given
 *                          after finalizing.
 * @param name              Function name, where the error was raised
 * @param text              Error message text
 */
void error(int test, int no, const char *name, const char *text)
{
	if(test != 0)
	{
		lprintf("ERROR", 0, "%s:\n%s\n", name, text);
		lprintf("ERROR", 0, "Exiting program...\n");
		if(no < 0)
		{
			exit(0);
		}
		else
		{
			finalize_process();
			exit(no);
		}
	}
}
