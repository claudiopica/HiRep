#include "geometry.h" 
#include "global.h" 
#include "error.h"
#include "logger.h"
#include <string.h>

static geometry_descriptor empty_gd={0};

static void copy_gd_no_comm(geometry_descriptor * to, geometry_descriptor * from){
*to=empty_gd;
to->inner_master_pieces = from->inner_master_pieces;
to->local_master_pieces = from->local_master_pieces;
to->total_spinor_master_pieces = from->total_spinor_master_pieces;
to->total_gauge_master_pieces = from->total_gauge_master_pieces;
to->master_start = from->master_start;
to->master_end = from->master_end;
to->master_shift = from->master_shift;
to->gsize_spinor = from->gsize_spinor;
to->gsize_gauge = from->gsize_gauge;
}


void init_geometry_SAP(){
int parity = (COORD[0]+COORD[1]+COORD[2]+COORD[3])&1;
/* Parity 0 is black */
	if(parity){
		copy_gd_no_comm(&glat_black,&glattice);
		glat_red = empty_gd;
		copy_gd_no_comm(&glat_even_black,&glat_even);
		glat_even_red = empty_gd;
		copy_gd_no_comm(&glat_odd_black,&glat_odd);
		glat_odd_red = empty_gd;
	}else{
		copy_gd_no_comm(&glat_red,&glattice);
		glat_black = empty_gd;
		copy_gd_no_comm(&glat_even_red,&glat_even);
		glat_even_black = empty_gd;
		copy_gd_no_comm(&glat_odd_red,&glat_odd);
		glat_odd_black = empty_gd;
	}	
}

void empty_buffers(spinor_field *s){
	for(int i=0;i<s->type->nbuffers_spinor;++i){
		memset(s->ptr + (s->type->rbuf_start[i]-s->type->master_shift),0,sizeof(suNf_spinor)*s->type->rbuf_len[i]);
	}
}

