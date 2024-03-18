/***************************************************************************\
 * Copyright (c) 2008-2024, Claudio Pica, Sofie Martins                     *
 * All rights reserved.                                                     *
 \**************************************************************************/

/**
 * @file geometry_type.h
 * @brief This file contains a geometry type for the new geometry
 */

#ifndef GEOMETRY_TYPE_H
#define GEOMETRY_TYPE_H

#ifdef WITH_NEW_GEOMETRY
typedef enum { EVEN = 1, ODD = 2, GLOBAL = 3 } gd_type;
#endif

#endif