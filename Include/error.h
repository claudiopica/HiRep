/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File error.h
* 
* Error handling functions
*
*******************************************************************************/

#ifndef ERROR_H
#define ERROR_H

void error(int test,int no,char *name,char *text);
#define null_error() error(0, 0, 0, 0)

#endif
