#ifndef MRE_H
#define MRE_H

#include "spinor_field.h"
#include "Inverters/linear_solvers.h"

#ifdef __cplusplus
	extern "C" {
#endif

/* functions and structures for the MRE algorithm */
typedef struct mre_par {
  spinor_field *s[2];
  int num[2];
  int max;
  int init;
} mre_par;

void mre_guess(mre_par *, int, spinor_field *, spinor_operator, spinor_field *);
void mre_store(mre_par *, int, spinor_field *);
void mre_init(mre_par *, int, double);


#ifdef __cplusplus
	}
#endif
#endif //MRE_H
