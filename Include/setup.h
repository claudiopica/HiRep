/***************************************************************************\
 * Copyright (c) 2008-2014, Vincent Drach                                   *
 * All rights reserved.                                                     *
 \**************************************************************************/

#ifndef SETUP_H
#define SETUP_H


char* get_input_filename();
char* get_output_filename();
char* get_error_filename();

int setup_process(int *argc, char ***argv);
int finalize_process(void);
void setup_gauge_fields();

#endif
