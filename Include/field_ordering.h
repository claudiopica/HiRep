/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef FIELD_ORDERING_H
#define FIELD_ORDERING_H

/* NB: it is assumed in the code that different directions are contiguous in memory */
#define coord_to_index(ix,mu) (((ix)<<2)|(mu))
#define index_to_coord(i,ix,mu) (mu)=((i)&3);(ix)=((i)>>2)

#endif /* FIELD_ORDERING_H */
