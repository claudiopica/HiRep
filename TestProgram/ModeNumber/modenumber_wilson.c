/*******************************************************************************
*
* Helper for check_modenumber_1 (not a test on its own: excluded in build.pl).
*
* Compiles the production ModeNumber/modenumber.c a second time, with the
* Wilson operator g5Dphi_sq in place of g5Cphi_sq and all public symbols
* renamed with a _wilson suffix. In a clover/expclover build this gives the
* Wilson reference from the very same code path.
*
*******************************************************************************/

/* library headers first, with the real names (they are include-guarded, so the
   #include in modenumber.c below is then a no-op) */
#include "libhr.h"

#define init_modenumber init_modenumber_wilson
#define free_modenumber free_modenumber_wilson
#define ModeNumber ModeNumber_wilson
#define g5Cphi_sq g5Dphi_sq

#include "../../ModeNumber/modenumber.c"
