/**
 * @file core_utils.h
 * @brief Macros to declare global variables 
 */

#ifndef CORE_UTILS_H
#define CORE_UTILS_H

#ifdef DEF_HIREP_GBL_VAR
#define GLB_VAR(type, name, ...) type name __VA_ARGS__
#else
#define GLB_VAR(type, name, ...) extern type name
#endif

#endif
