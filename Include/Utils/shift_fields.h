#ifndef SHIFT_FIELDS_H
#define SHIFT_FIELDS_H

#include "spinor_field.h"

#ifdef __cplusplus
	extern "C" {
#endif

/*Global shift for fields, the routine accepts also NULL entries in which case it does nothing*/
void shift_fields(int *shift, spinor_field *sin, suNg_field *uin, spinor_field *sout, suNg_field *uout);

#ifdef __cplusplus
	}
#endif
#endif //SHIFT_FIELDS_H
