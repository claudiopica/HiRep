/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef SAFE_MOD_H
#define SAFE_MOD_H

#include <stdlib.h>

static int safe_mod(int x,int y)
{
   if (x>=0)
      return(x%y);
   else
      return((y-(abs(x)%y))%y);
}

#endif
