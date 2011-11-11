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

#ifdef __cplusplus
extern "C" {
#endif

void error(int test,int no,char *name,char *text);

#ifdef __cplusplus
}
#endif

#endif
