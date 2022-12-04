/***************************************************************************\
 * Copyright (c) 2008-2014, Vincent Drach                                   *
 * All rights reserved.                                                     *
 \**************************************************************************/

/**
 * @file setup.h
 * @brief Setup and finalize to run at the beginning and the end of every 
 *        program
 */ 

#ifndef SETUP_H
#define SETUP_H

#ifdef WITH_GPU
    #include "gpu.h"
#endif

/**
 * @brief Read input filename from command line
 */
char* get_input_filename();

/**
 * @brief Get output filename, default out_0 in the local directory
 */
char* get_output_filename();

/**
 * @brief Get filename of file to print errors to, default is err_0 in
 *        the current directory.
 */
char* get_error_filename();

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
int finalize_process(void);

/**
 * @brief Initialize gauge fields
 */
void setup_gauge_fields();

#ifdef WITH_GPU
    /**
     * @brief Call this in an init function to setup available graphics card for
     *        use. This also logs information on available software and hardware.
     * 
     * @param input_gpu             A struct containing information on the current active
     *                              GPU
     */
    void init_gpu(input_gpu gpu_var);
#endif


#endif
