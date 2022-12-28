#ifndef READ_ACTION_H
#define READ_ACTION_H
#include "Update/integrators.h"

#ifdef __cplusplus
	extern "C" {
#endif

void read_action(char *filename, integrator_par **ipp);

#ifdef __cplusplus
	}
#endif
#endif //READ_ACTION_H
