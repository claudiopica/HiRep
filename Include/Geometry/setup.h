/***************************************************************************\
 * Copyright (c) 2008-2014, Vincent Drach                                   *
 * All rights reserved.                                                     *
 \**************************************************************************/

/// Headerfile for:
/// - process_init.c

/**
 * @file setup.h
 * @brief Setup and finalize to run at the beginning and the end of every 
 *        program
 */

#ifndef SETUP_H
#define SETUP_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Read input filename from command line
 */
char *get_input_filename();

/**
 * @brief Get output filename, default out_0 in the local directory
 */
char *get_output_filename();

/**
 * @brief Get filename of file to print errors to, default is err_0 in
 *        the current directory.
 */
char *get_error_filename();

/**
 * @brief Setup the process at the beginning of each run
 *
 * @param argc			Command line input from main
 * @param argv			Command line input from main
 */
int setup_process(int *argc, char ***argv);

/**
 * @brief Finalize process at the end of each run
 */
void finalize_process(void);

/**
 * @brief Initialize gauge fields
 */
void setup_gauge_fields();

#ifdef __cplusplus
}
#endif
#endif
