#ifndef CHECK_OPTIONS_H
#define CHECK_OPTIONS_H

#include "suN_types.h"

/* Boundary conditions */

#ifdef BC_XYZ_TWISTED

# ifndef REPR_ADJOINT
#   error Twisted boundary conditions can be used only with the adjoint representation!!!
# endif

# undef BC_X_PERIODIC
# undef BC_X_ANTIPERIODIC
# undef BC_X_THETA
# define BC_X_PERIODIC
# define BC_X_ALREADY

# undef BC_Y_PERIODIC
# undef BC_Y_ANTIPERIODIC
# undef BC_Y_THETA
# define BC_Y_PERIODIC
# define BC_Y_ALREADY

# undef BC_Z_PERIODIC
# undef BC_Z_ANTIPERIODIC
# undef BC_Z_THETA
# define BC_Z_PERIODIC
# define BC_Z_ALREADY

# define PLAQ_WEIGHTS

# if defined(BASIC_SF) || defined(ROTATED_SF)
#   error (BC_XYZ_TWISTED) Twisted BCs cannot be used with Schroedinger functional!!!
# endif

#endif



#ifdef BASIC_SF

# undef BC_T_PERIODIC
# undef BC_T_ANTIPERIODIC
# undef BC_T_OPEN
# undef BC_T_THETA
# define BC_T_ALREADY

# define PLAQ_WEIGHTS

# ifdef ROTATED_SF
#   error (BASIC_SF) BASIC_SF and ROTATED_SF cannot be used at the same time!!!
# endif

#endif



#ifdef ROTATED_SF

# undef BC_T_PERIODIC
# undef BC_T_ANTIPERIODIC
# undef BC_T_OPEN
# undef BC_T_THETA
# define BC_T_ALREADY

# undef BC_X_PERIODIC
# undef BC_X_ANTIPERIODIC
# undef BC_X_THETA
# define BC_X_THETA

# undef BC_Y_PERIODIC
# undef BC_Y_ANTIPERIODIC
# undef BC_Y_THETA
# define BC_Y_THETA

# undef BC_Z_PERIODIC
# undef BC_Z_ANTIPERIODIC
# undef BC_Z_THETA
# define BC_Z_THETA

# define PLAQ_WEIGHTS

# ifdef BASIC_SF
#   error (ROTATED_SF) BASIC_SF and ROTATED_SF cannot be used at the same time!!!
# endif

#endif

#if defined(HALFBG_SF) && !( (NG == 2) && ( defined(BASIC_SF) || defined(ROTATED_SF) ) )
#   error (HALFBG_SF) can be defined only if NG=2 and or BASIC_SF or ROTATED_SF is used!!!
# endif


#ifdef BC_T_ANTIPERIODIC

# ifdef BC_T_ALREADY
#   error (BC_T_ANTIPERIODIC) BC_T already defined!!!
# endif
# define BC_T_ALREADY

#endif



#ifdef BC_T_OPEN

# ifdef BC_T_ALREADY
#   error (BC_T_OPEN) BC_T already defined!!!
# endif
# define BC_T_ALREADY

# undef PLAQ_WEIGHTS
# define PLAQ_WEIGHTS

#endif



#ifdef BC_T_THETA

# ifdef BC_T_ALREADY
#   error (BC_T_OPEN) BC_T already defined!!!
# endif
# define BC_T_ALREADY
# define FERMION_THETA

#endif



#if !defined(BC_T_ALREADY) && !defined(BC_T_PERIODIC)
# define BC_T_PERIODIC
# define BC_T_ALREADY
#endif



#ifdef BC_X_ANTIPERIODIC

# ifdef BC_X_ALREADY
#   error (BC_X_ANTIPERIODIC) BC_X already defined!!!
# endif
# define BC_X_ALREADY

#endif



#ifdef BC_X_THETA

# ifdef BC_X_ALREADY
#   error (BC_X_OPEN) BC_X already defined!!!
# endif
# define BC_X_ALREADY
# define FERMION_THETA

#endif



#if !defined(BC_X_ALREADY) && !defined(BC_X_PERIODIC)
# define BC_X_PERIODIC
# define BC_X_ALREADY
#endif



#ifdef BC_Y_ANTIPERIODIC

# ifdef BC_Y_ALREADY
#   error (BC_Y_ANTIPERIODIC) BC_Y already defined!!!
# endif
# define BC_Y_ALREADY

#endif



#ifdef BC_Y_THETA

# ifdef BC_Y_ALREADY
#   error (BC_Y_OPEN) BC_Y already defined!!!
# endif
# define BC_Y_ALREADY
# define FERMION_THETA

#endif



#if !defined(BC_Y_ALREADY) && !defined(BC_Y_PERIODIC)
# define BC_Y_PERIODIC
# define BC_Y_ALREADY
#endif



#ifdef BC_Z_ANTIPERIODIC

# ifdef BC_Z_ALREADY
#   error (BC_Z_ANTIPERIODIC) BC_Z already defined!!!
# endif
# define BC_Z_ALREADY

#endif



#ifdef BC_Z_THETA

# ifdef BC_Z_ALREADY
#   error (BC_Z_OPEN) BC_Z already defined!!!
# endif
# define BC_Z_ALREADY
# define FERMION_THETA

#endif



#if !defined(BC_Z_ALREADY) && !defined(BC_Z_PERIODIC)
# define BC_Z_PERIODIC
# define BC_Z_ALREADY
#endif


#endif /* CHECK_OPTIONS_H */

