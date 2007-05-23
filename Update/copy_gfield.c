#include "global.h"
#include "update.h"
#include <string.h>

/* g1=g2 */
void suNg_field_copy(suNg *g1, suNg *g2)
{
	memcpy(g1,g2,4*VOLUME*sizeof(*g1));
}

