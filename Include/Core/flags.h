/***************************************************************************\
* Copyright (c) 2024, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file flags.h
 * @brief Global flags that indicate for example that a certain field 
 *        has been updated.
 */

#ifndef FLAGS_H
#define FLAGS_H

GLB_VAR(int, stale_clover, = 0);

#ifdef WITH_EXPCLOVER
GLB_VAR(int, stale_expclover, = 0);
#endif

#endif