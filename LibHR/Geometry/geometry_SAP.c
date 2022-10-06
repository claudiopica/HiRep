#include "geometry.h"
#include "global.h"
#include "error.h"
#include "logger.h"
#include <string.h>

static geometry_descriptor empty_gd = {0};

static void copy_gd_no_comm(geometry_descriptor *dptr, geometry_descriptor *sptr)
{
	*dptr = empty_gd;
	dptr->inner_master_pieces = sptr->inner_master_pieces;
	dptr->local_master_pieces = sptr->local_master_pieces;
	dptr->total_spinor_master_pieces = sptr->total_spinor_master_pieces;
	dptr->total_gauge_master_pieces = sptr->total_gauge_master_pieces;
	dptr->master_start = sptr->master_start;
	dptr->master_end = sptr->master_end;
	dptr->master_shift = sptr->master_shift;
	dptr->gsize_spinor = sptr->gsize_spinor;
	dptr->gsize_gauge = sptr->gsize_gauge;
}

void init_geometry_SAP()
{
	int parity = (COORD[0]+COORD[1]+COORD[2]+COORD[3])&1;

	// Parity 0 is red
	// Parity 1 is black
	if(parity)
	{
		copy_gd_no_comm(&glat_black, &glattice);
		glat_red = empty_gd;
		copy_gd_no_comm(&glat_even_black, &glat_even);
		glat_even_red = empty_gd;
		copy_gd_no_comm(&glat_odd_black, &glat_odd);
		glat_odd_red = empty_gd;
	}
	else
	{
		copy_gd_no_comm(&glat_red, &glattice);
		glat_black = empty_gd;
		copy_gd_no_comm(&glat_even_red, &glat_even);
		glat_even_black = empty_gd;
		copy_gd_no_comm(&glat_odd_red, &glat_odd);
		glat_odd_black = empty_gd;
	}
}

void empty_buffers(spinor_field *s)
{
	int shift, size;
	for(int i = 0; i < s->type->nbuffers_spinor; i++)
	{
		shift = s->type->rbuf_start[i] - s->type->master_shift;
		size = sizeof(suNf_spinor) * s->type->rbuf_len[i];
		memset(s->ptr+shift, 0, size);
	}
}
