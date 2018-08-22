#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>

#include "geometry.h" 
#include "global.h" 
#include "error.h"
#include "logger.h"



#define true 1
#define even 0
#define false 0
#define odd 1
#define no_eo 2

#define TDIR 0
#define XDIR 1
#define YDIR 2
#define ZDIR 3

#define REPORTLVL 200000

#define local_index(nt,nx,ny,nz)  (((nt)+2*T_BORDER+T)%(2*T_BORDER+T)+	\
				   (T+2*T_BORDER)*(((nx)+2*X_BORDER+X)%(2*X_BORDER+X))+ \
				   (T+2*T_BORDER)*(X+2*X_BORDER)*(((ny)+2*Y_BORDER+Y)%(2*Y_BORDER+Y))+ \
				   (T+2*T_BORDER)*(X+2*X_BORDER)*(Y+2*Y_BORDER)*(((nz)+2*Z_BORDER+Z)%(2*Z_BORDER+Z)))


typedef struct
{
  int b_start[4];
  int b_width[4];
  int id_mask;
  int checksum[8];
  int index_start;
  int index_end;
  int size;
  int eo_type;
  int id_zone;
  int level;
  int id_proc; /*id of the processor to which i should send or from whom i would like to recive data*/
} border_id;

static int OVERSIZE_VOLUME;
static int TOTAL_VOLUME;
static int N_BORDER;

static int ** dir_mask;
static int ** inv_mask;

static int index_position=0;
static int index_border=-1;
static int index_zone=0;
static int first_l2_border=-1;

static int *map_id2overlexi=NULL;
static int *memory=NULL;
static int *map_overlexi2id=NULL;
static border_id *border=NULL;
static char dc[4]={'T','X','Y','Z'};
static int inner[2]={-1,-1};



static void walk_on_lattice(int id_mask,int eotype,int level,int id_zone, int* bl_start,int* incr,int* bl_width);
static int init_border(int id_mask, int eotype,int level, int id_zone, int* bl_start,int* incr,int* bl_width,int * match_length);
static void close_border();
static void set_border_pointer(int actualn, int match_p);


static void set_block(int id_mask, int* bl_s,int* incr, int* bl_w, int pts, int wts, int pxs, int wxs,int pys, int wys, int pzs, int wzs){
  bl_s[dir_mask[id_mask][TDIR]]=pts;
  bl_s[dir_mask[id_mask][XDIR]]=pxs;
  bl_s[dir_mask[id_mask][YDIR]]=pys;
  bl_s[dir_mask[id_mask][ZDIR]]=pzs;

  incr[dir_mask[id_mask][TDIR]]=(pts+wts<=T+2*T_BORDER)?+1:-1;
  incr[dir_mask[id_mask][XDIR]]=(pxs+wxs<=X+2*X_BORDER)?+1:-1;
  incr[dir_mask[id_mask][YDIR]]=(pys+wys<=Y+2*Y_BORDER)?+1:-1;
  incr[dir_mask[id_mask][ZDIR]]=(pzs+wzs<=Z+2*Z_BORDER)?+1:-1;

  bl_w[dir_mask[id_mask][TDIR]]=(pts+wts<=T+2*T_BORDER)?wts:-wts;
  bl_w[dir_mask[id_mask][XDIR]]=(pxs+wxs<=X+2*X_BORDER)?wxs:-wxs;
  bl_w[dir_mask[id_mask][YDIR]]=(pys+wys<=Y+2*Y_BORDER)?wys:-wys;
  bl_w[dir_mask[id_mask][ZDIR]]=(pzs+wzs<=Z+2*Z_BORDER)?wzs:-wzs;

}


static int zone_border(int ix){
  int i;
  int output=false;
  for(i=1;i<index_border;i++)
    if(border[i-1].id_zone!=border[i].id_zone)
      if(ix==border[i].index_start) output=true;
  return output;
}

static unsigned int block_cond(int b1,int b2,int x){
  if (x>=b1 && x<b2) return true;
  if (x<=b1 && x>b2) return true;
  return false;
} 


static void geometry_mpi_init()
{
  int i,BOR_CUBE,BOR_SQUARE,L3_BORDER,L2_BORDER;
  
  /* X_BORDER=(NP_X>1)?BORDERSIZE:0; */
  /* Y_BORDER=(NP_Y>1)?BORDERSIZE:0; */
  /* Z_BORDER=(NP_Z>1)?BORDERSIZE:0; */
  /* T_BORDER=(NP_T>1)?BORDERSIZE:0; */

 
  
  /* error(T-2*T_BORDER<0,1,"geometry_mpi_init [geometry_mpi.c]","Too large T Border in the geometry"); */
  /* error(X-2*X_BORDER<0,1,"geometry_mpi_init [geometry_mpi.c]","Too large X Border in the geometry"); */
  /* error(Y-2*Y_BORDER<0,1,"geometry_mpi_init [geometry_mpi.c]","Too large Y Border in the geometry"); */
  /* error(Z-2*Z_BORDER<0,1,"geometry_mpi_init [geometry_mpi.c]","Too large Z Border in the geometry"); */
  
  /* X_EXT=X+2*X_BORDER; */
  /* Y_EXT=Y+2*Y_BORDER; */
  /* Z_EXT=Z+2*Z_BORDER; */
  /* T_EXT=T+2*T_BORDER; */
  
  
  OVERSIZE_VOLUME=((X+2*X_BORDER)*(Y+2*Y_BORDER)*(Z+2*Z_BORDER)*(T+2*T_BORDER));
  
  BOR_CUBE=2*(X*Y*Z*T_BORDER+
	      X*Y*Z_BORDER*T+
	      X*Y_BORDER*Z*T+
	      X_BORDER*Y*Z*T);

  BOR_SQUARE=4*(X*Y*Z_BORDER*T_BORDER+
		X*Y_BORDER*Z*T_BORDER+
		X*Y_BORDER*Z_BORDER*T+
		X_BORDER*Y*Z*T_BORDER+
		X_BORDER*Y*Z_BORDER*T+
		X_BORDER*Y_BORDER*Z*T);

   
  TOTAL_VOLUME=(X-2*X_BORDER)*(Y-2*Y_BORDER)*(Z-2*Z_BORDER)*(T-2*T_BORDER)+
    2*BOR_CUBE+2*BOR_SQUARE;
  
  L3_BORDER = 2*( (X_BORDER!=0 ? 1 : 0) + (Y_BORDER!=0 ? 1 : 0)  + (Z_BORDER!=0 ? 1 : 0)  + (T_BORDER!=0 ? 1 : 0) );
  
  L2_BORDER = L3_BORDER*(L3_BORDER-2)/2;

  N_BORDER = 2*L3_BORDER+2*L2_BORDER; /* this 2 is due to Even/Odd */

  lprintf("GEOMETRY",200,"NBORDER=%d L3BORDER=%d\n",N_BORDER,L3_BORDER);
  

  /*   printf("volume di un bordo %d\n",BOR_CUBE/4+BOR_SQUARE/4);   */
  /*   printf("volume totale %d\n",TOTAL_VOLUME); */

  map_id2overlexi=malloc(TOTAL_VOLUME*sizeof(int));
  error((map_id2overlexi==NULL),0,"geometry_mpi_init [geometry_mpi.c]",
	"Cannot allocate memory for map_id2overlexi");
  
  map_overlexi2id=malloc(OVERSIZE_VOLUME*sizeof(int));
  error((map_overlexi2id==NULL),0,"geometry_mpi_init [geometry_mpi.c]",
	"Cannot allocate memory for map_overlexi2id");
 
  memory=malloc(TOTAL_VOLUME*sizeof(int));
  error((memory==NULL),0,"geometry_mpi_init [geometry_mpi.c]",
	"Cannot allocate memory for memory");

  border=malloc((2+2*N_BORDER)*sizeof(border_id));
  error((border==NULL),1,"geometry_mpi_init [geometry_mpi.c]",
	"Cannot allocate memory for border");

  if(N_BORDER>0){
    glattice.rbuf_len=malloc((N_BORDER)*sizeof(int));
    error((glattice.rbuf_len==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glattice.sbuf_len=malloc((N_BORDER)*sizeof(int));
    error((glattice.sbuf_len==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glattice.sbuf_to_proc=malloc((N_BORDER)*sizeof(int));
    error((glattice.sbuf_to_proc==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glattice.sbuf_start=malloc((N_BORDER)*sizeof(int));
    error((glattice.sbuf_start==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glattice.rbuf_from_proc=malloc((N_BORDER)*sizeof(int));
    error((glattice.rbuf_from_proc==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glattice.rbuf_start=malloc((N_BORDER)*sizeof(int));
    error((glattice.rbuf_start==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
  }
  glattice.nbuffers_gauge = N_BORDER;

  /*Setting glat_even & glat_odd values*/

  if(L3_BORDER>0){
    glat_even.rbuf_len=malloc((L3_BORDER)*sizeof(int));
    error((glat_even.rbuf_len==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glat_even.sbuf_len=malloc((L3_BORDER)*sizeof(int));
    error((glat_even.sbuf_len==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glat_even.sbuf_to_proc=malloc((L3_BORDER)*sizeof(int));
    error((glat_even.sbuf_to_proc==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glat_even.sbuf_start=malloc((L3_BORDER)*sizeof(int));
    error((glat_even.sbuf_start==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glat_even.rbuf_from_proc=malloc((L3_BORDER)*sizeof(int));
    error((glat_even.rbuf_from_proc==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glat_even.rbuf_start=malloc((L3_BORDER)*sizeof(int));
    error((glat_even.rbuf_start==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glat_odd.rbuf_len=malloc((L3_BORDER)*sizeof(int));
    error((glat_odd.rbuf_len==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glat_odd.sbuf_len=malloc((L3_BORDER)*sizeof(int));
    error((glat_odd.sbuf_len==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glat_odd.sbuf_to_proc=malloc((L3_BORDER)*sizeof(int));
    error((glat_odd.sbuf_to_proc==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glat_odd.sbuf_start=malloc((L3_BORDER)*sizeof(int));
    error((glat_odd.sbuf_start==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glat_odd.rbuf_from_proc=malloc((L3_BORDER)*sizeof(int));
    error((glat_odd.rbuf_from_proc==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
    
    glat_odd.rbuf_start=malloc((L3_BORDER)*sizeof(int));
    error((glat_odd.rbuf_start==NULL),1,"fix_geometry_descriptor [geometry_mpi.c]",
	  "Cannot allocate memory");
  }
  glat_odd.nbuffers_spinor = L3_BORDER;

  glat_even.nbuffers_spinor = L3_BORDER;

  glat_odd.nbuffers_gauge = -1;

  glat_even.nbuffers_gauge = -1;

  glattice.nbuffers_spinor = 2*L3_BORDER;

  glattice.inner_master_pieces=0;
  glat_even.inner_master_pieces=0;  
  glat_odd.inner_master_pieces=0;

  glat_odd.gsize_gauge=-1;
  glat_even.gsize_gauge=-1;

  for(i=0;i<OVERSIZE_VOLUME;i++)
    {
      map_overlexi2id[i]=-1;
    }

  for(i=0;i<TOTAL_VOLUME;i++)
    {
      map_id2overlexi[i]=-1;
      memory[i]=-1;
    }

  dir_mask=malloc(2*sizeof(int*));
  inv_mask=malloc(2*sizeof(int*));
  for(i=0;i<2;i++){
    dir_mask[i]=malloc(4*sizeof(int));
    error((dir_mask[i]==NULL),1,"geometry_mpi_init [geometry_mpi.c]",
	  "Cannot allocate memory");
    inv_mask[i]=malloc(4*sizeof(int));
    error((inv_mask[i]==NULL),1,"geometry_mpi_init [geometry_mpi.c]",
	  "Cannot allocate memory");
  }

  i=0;
  if(NP_X==1) {
    dir_mask[0][XDIR]=i;
    i++;}
  
  if(NP_Y==1) {
    dir_mask[0][YDIR]=i;
    i++;}
  
  if(NP_Z==1) {
    dir_mask[0][ZDIR]=i;
    i++;}

  if(NP_T==1) {
    dir_mask[0][TDIR]=i;
    i++;}
  
  i=3;
  if(NP_T>1) {
    dir_mask[0][TDIR]=i;
    i--;}

  if(NP_Z>1) {
    dir_mask[0][ZDIR]=i;
    i--;}

  if(NP_Y>1) {
    dir_mask[0][YDIR]=i;
    i--;}

  if(NP_X>1) {
    dir_mask[0][XDIR]=i;
    i--;}


  for(i=0;i<4;i++)
    dir_mask[1][i]=(dir_mask[0][i]+2)%4;

  for(i=0;i<4;i++){
    inv_mask[0][dir_mask[0][i]]=i;
    inv_mask[1][dir_mask[1][i]]=i;
  }

  for(i=0;i<4;i++)
    lprintf("GEOMETRY" ,REPORTLVL,"dirmask[%c] =  %d \tinvmask[%d] =  %c\n",dc[i], dir_mask[0][i],i,dc[inv_mask[0][i]]);

  for(i=0;i<4;i++)
    lprintf("GEOMETRY" ,REPORTLVL,"dirmask[%c] =  %d \tinvmask[%d] =  %c\n",dc[i], dir_mask[1][i],i,dc[inv_mask[1][i]]);



}

static void geometry_mpi_finalize(){
  if(map_id2overlexi!=NULL) free(map_id2overlexi);
  if(map_overlexi2id!=NULL) free(map_overlexi2id);
  if(border!=NULL) free(border);   
  if(memory!=NULL) free(memory);   

  if(dir_mask!=NULL){
    int i;
    for(i=0;i<2;i++){
      free(dir_mask[i]);
      free(inv_mask[i]);
    }
    free(dir_mask);    
    free(inv_mask);    
  }

  lprintf("GEOMETRY" ,0,"Gauge field: size %d nbuffer %d \n",glattice.gsize_gauge,glattice.nbuffers_gauge);
  lprintf("GEOMETRY" ,0,"Spinor field (EO): size %d nbuffer %d \n",glattice.gsize_spinor,glattice.nbuffers_spinor);
  lprintf("GEOMETRY" ,0,"Even Spinor field: size %d nbuffer %d \n",glat_even.gsize_spinor,glat_even.nbuffers_spinor);
  lprintf("GEOMETRY" ,0,"Odd Spinor field: size %d nbuffer %d \n",glat_odd.gsize_spinor,glat_odd.nbuffers_spinor);





}

static void set_inner(int eotype){
  int bl_start[4];
  int bl_width[4];
  int incr[4];
  lprintf("INNER",REPORTLVL,"START\n");
  set_block(0,bl_start,incr,bl_width,
	    2*T_BORDER,T-2*T_BORDER,
	    2*X_BORDER,X-2*X_BORDER,
	    2*Y_BORDER,Y-2*Y_BORDER,
	    2*Z_BORDER,Z-2*Z_BORDER);

  /* if(inner[0]==-1) inner[0]=index_border+1; */
  /* else inner[1]=index_border+1; */

  inner[eotype]=index_border+1;

  walk_on_lattice(0,eotype,4,index_zone,bl_start,incr,bl_width);

  if( border[index_border].size !=0){
    if(eotype==even){
      glattice.inner_master_pieces++;
      glat_even.inner_master_pieces++;
    }  
    if(eotype==odd){
      glattice.inner_master_pieces++;
      glat_odd.inner_master_pieces++;
    }
  } 
  index_zone++;
  lprintf("INNER",REPORTLVL,"END info: start%d end %d \n",border[index_border].index_start,border[index_border].index_end);
}


static void walk_on_lattice(int id_mask,int eotype,int level,int id_zone, int* bl_start,int* incr,int* bl_width)
{
  int x0,x1,x2,x3;
  int x[4];
  int match_length,match_point,steps=0,last_point;
  lprintf("GEOMETRY",REPORTLVL,"\n\n\nwalk_on_lattice The START (%d,%d,%d,%d) \n",bl_start[dir_mask[id_mask][0]],bl_start[dir_mask[id_mask][1]],bl_start[dir_mask[id_mask][2]],bl_start[dir_mask[id_mask][3]]);
  lprintf("GEOMETRY",REPORTLVL,"walk_on_lattice The WIDTH (%d,%d,%d,%d) \n",bl_width[dir_mask[id_mask][0]],bl_width[dir_mask[id_mask][1]],bl_width[dir_mask[id_mask][2]],bl_width[dir_mask[id_mask][3]]);
  lprintf("GEOMETRY",REPORTLVL,"walk_on_lattice The incr (%d,%d,%d,%d) \n",incr[dir_mask[id_mask][0]],incr[dir_mask[id_mask][1]],incr[dir_mask[id_mask][2]],incr[dir_mask[id_mask][3]]);
  lprintf("GEOMETRY",REPORTLVL,"walk_on_lattice The EO (%d) \n",eotype);
  if(init_border(id_mask,eotype,level,id_zone,bl_start,incr,bl_width,&match_length)){
    for (x3=bl_start[3];block_cond(bl_start[3],bl_start[3]+bl_width[3],x3);x3+=incr[3]) 
      for (x2=bl_start[2];block_cond(bl_start[2],bl_start[2]+bl_width[2],x2);x2+=incr[2]) 
	for (x1=bl_start[1];block_cond(bl_start[1],bl_start[1]+bl_width[1],x1);x1+=incr[1]) 
	  for (x0=bl_start[0];block_cond(bl_start[0],bl_start[0]+bl_width[0],x0);x0+=incr[0]){
	    x[inv_mask[id_mask][0]]=x0;
	    x[inv_mask[id_mask][1]]=x1;
	    x[inv_mask[id_mask][2]]=x2;
	    x[inv_mask[id_mask][3]]=x3;
	    if(eotype==no_eo || eotype==(x0+x1+x2+x3+T_BORDER+X_BORDER+Y_BORDER+Z_BORDER+PSIGN)%2){
	      lprintf("GEOMETRY",REPORTLVL,"walk_on_lattice T[%d] X[%d] Y[%d] Z[%d] \n",x[0],x[1],x[2],x[3]); 
	      if(match_length>steps) match_point=true;
	      else match_point=false;
	      steps++;
	      set_border_pointer(local_index(x[0],x[1],x[2],x[3]),match_point);

	    }
	  }


    last_point=border[index_border].index_start+steps;
    close_border(last_point);
 
  }
  
}

static void border_checksum(int* checksum,int id_mask, int* bl_start,int* bl_width )
{
  checksum[0] = bl_start[dir_mask[id_mask][0]];
  checksum[1] = bl_start[dir_mask[id_mask][1]];
  checksum[2] = bl_start[dir_mask[id_mask][2]];
  checksum[3] = bl_start[dir_mask[id_mask][3]];
  checksum[4] = bl_width[dir_mask[id_mask][0]];
  checksum[5] = bl_width[dir_mask[id_mask][1]];
  checksum[6] = bl_width[dir_mask[id_mask][2]];
  checksum[7] = bl_width[dir_mask[id_mask][3]];
  if( checksum[4]<0) {
    checksum[4]=-checksum[4];
    checksum[0]-=(checksum[4]-T_BORDER);
  }
  if( checksum[5]<0) {
    checksum[5]=-checksum[5];
    checksum[1]-=(checksum[5]-X_BORDER);
  }
  if( checksum[6]<0) {
    checksum[6]=-checksum[6];
    checksum[2]-=(checksum[6]-Y_BORDER);
  }
  if( checksum[7]<0) {
    checksum[7]=-checksum[7];
    checksum[3]-=(checksum[7]-Z_BORDER);
  }
}

static void close_border(int i ){
  border[index_border].index_end = i;
}


static int init_border(int id_mask, int eotype,int level,int id_zone,  int* bl_start,int *incr,int* bl_width,int* match_length)
{
  int i,retval,check_index,match_step,match_interval,length_match=0,id,found_start,x1,x2,x3,x0,start=0,index_start;
  int x[4];
  int checksum[8];
  int bl_size,foundstart;
  retval=true;
  bl_size= (bl_width[0]*bl_width[1]*bl_width[2]*bl_width[3]>0)?bl_width[0]*bl_width[1]*bl_width[2]*bl_width[3]:-bl_width[0]*bl_width[1]*bl_width[2]*bl_width[3];
  border_checksum(checksum,id_mask,bl_start,bl_width );



  for(i=0;i<index_border+1;i++){
    if( border[i].checksum[0] == checksum[0] &&
	border[i].checksum[1] == checksum[1] &&
	border[i].checksum[2] == checksum[2] &&
	border[i].checksum[3] == checksum[3] &&
	border[i].checksum[4] == checksum[4] &&
	border[i].checksum[5] == checksum[5] &&
	border[i].checksum[6] == checksum[6] &&
	border[i].checksum[7] == checksum[7] &&
	border[i].eo_type == eotype )  retval = false;
  }
  
  if( bl_size == 0 && inner[0] != index_border &&  inner[1] != index_border ) retval = false;
  
  if(!retval){
    lprintf("GEOMETRY",REPORTLVL,"init_border Border already evaluated or empty \n");
    return(retval);
  }

  lprintf("GEOMETRY",REPORTLVL,"init_border Setting a new border \n");
  x[inv_mask[id_mask][0]]=bl_start[0];
  x[inv_mask[id_mask][1]]=bl_start[1];
  x[inv_mask[id_mask][2]]=bl_start[2];
  x[inv_mask[id_mask][3]]=bl_start[3];
  if(eotype == no_eo) 
    start=map_overlexi2id[local_index(x[0],x[1],x[2],x[3])];
  else {
    bl_size=0;
    foundstart=false;
    for (x3=bl_start[3];block_cond(bl_start[3],bl_start[3]+bl_width[3],x3);x3+=incr[3])
      for (x2=bl_start[2];block_cond(bl_start[2],bl_start[2]+bl_width[2],x2);x2+=incr[2])
	for (x1=bl_start[1];block_cond(bl_start[1],bl_start[1]+bl_width[1],x1);x1+=incr[1])
	  for (x0=bl_start[0];block_cond(bl_start[0],bl_start[0]+bl_width[0],x0);x0+=incr[0]){
	    x[inv_mask[id_mask][0]]=x0;
	    x[inv_mask[id_mask][1]]=x1;
	    x[inv_mask[id_mask][2]]=x2;
	    x[inv_mask[id_mask][3]]=x3;
	    if((x0+x1+x2+x3+T_BORDER+X_BORDER+Y_BORDER+Z_BORDER+PSIGN)%2==eotype){
	      if(!foundstart){
		start=map_overlexi2id[local_index(x[0],x[1],x[2],x[3])];
		foundstart=true;
	      }
	      bl_size++;
	    }
	  }
  }


  
  /* if can append */
  found_start=false;
  int index_start_out=0,length_match_out=0;
  for(i=index_position-1; i>=0;i--){
    if(memory[i]==start){
      index_start=i;
      found_start=true;
      length_match=(index_position-i>bl_size)?bl_size:index_position-i;
      lprintf("GEOMETRY",REPORTLVL,"init_border Found a candidate start:%d for a total of %d steps\n",i,length_match);
      match_interval=true;
      match_step=true;
      check_index=index_start;
      for (x3=bl_start[3];block_cond(bl_start[3],bl_start[3]+bl_width[3],x3);x3+=incr[3])
	for (x2=bl_start[2];block_cond(bl_start[2],bl_start[2]+bl_width[2],x2);x2+=incr[2])
	  for (x1=bl_start[1];block_cond(bl_start[1],bl_start[1]+bl_width[1],x1);x1+=incr[1])
	    for (x0=bl_start[0];block_cond(bl_start[0],bl_start[0]+bl_width[0],x0);x0+=incr[0]){
	      x[inv_mask[id_mask][0]]=x0;
	      x[inv_mask[id_mask][1]]=x1;
	      x[inv_mask[id_mask][2]]=x2;
	      x[inv_mask[id_mask][3]]=x3;
	      if((x0+x1+x2+x3+T_BORDER+X_BORDER+Y_BORDER+Z_BORDER+PSIGN)%2==eotype || eotype ==no_eo ) {
		if( map_overlexi2id[local_index(x[0],x[1],x[2],x[3])]!= memory[check_index]) match_step=false;
		check_index++;
		if(check_index-index_start==length_match && !match_step) match_interval=false;
	      }
	    }
      
      if(match_interval && length_match_out<length_match){
	length_match_out=length_match;
	index_start_out=index_start;
      }
      else 
	lprintf("GEOMETRY",REPORTLVL,"init_border Not a good match\n");
    }
  }
  
  
  
  if(!found_start ||length_match_out==0 ){
    lprintf("GEOMETRY",REPORTLVL,"init_border Unable to find a match\n");
    index_start=index_position;
    *match_length=0;
  } else {
    *match_length=length_match_out;
    length_match=length_match_out;
    index_start=index_start_out;
    lprintf("GEOMETRY",REPORTLVL,"init_border Good match -> start %d of %d steps \n",index_start,length_match);
  }


  index_border++;

  
  if( bl_size == 0 ) index_start=index_position;

  border[index_border].b_start[0] = bl_start[dir_mask[id_mask][0]];
  border[index_border].b_start[1] = bl_start[dir_mask[id_mask][1]];
  border[index_border].b_start[2] = bl_start[dir_mask[id_mask][2]];
  border[index_border].b_start[3] = bl_start[dir_mask[id_mask][3]];
  border[index_border].b_width[0] = bl_width[dir_mask[id_mask][0]];
  border[index_border].b_width[1] = bl_width[dir_mask[id_mask][1]];
  border[index_border].b_width[2] = bl_width[dir_mask[id_mask][2]];
  border[index_border].b_width[3] = bl_width[dir_mask[id_mask][3]];
  border[index_border].checksum[0] = checksum[0];
  border[index_border].checksum[1] = checksum[1];
  border[index_border].checksum[2] = checksum[2];
  border[index_border].checksum[3] = checksum[3];
  border[index_border].checksum[4] = checksum[4];
  border[index_border].checksum[5] = checksum[5];
  border[index_border].checksum[6] = checksum[6];
  border[index_border].checksum[7] = checksum[7];
  border[index_border].id_mask = id_mask;
  border[index_border].id_zone = id_zone;
  border[index_border].eo_type = eotype;
  border[index_border].size = bl_size;
  border[index_border].level = level;
  border[index_border].index_start = index_start ;
  
  id=CID;

  if(bl_width[dir_mask[id_mask][TDIR]] == T_BORDER )
    {
      if(bl_start[dir_mask[id_mask][TDIR]] == T_BORDER || bl_start[dir_mask[id_mask][TDIR]] == 0 ) id = proc_dn(id,0);
      else if(bl_start[dir_mask[id_mask][TDIR]] == T || bl_start[dir_mask[id_mask][TDIR]] == T + T_BORDER) id = proc_up(id,0);
    }
  
  if(bl_width[dir_mask[id_mask][XDIR]] == X_BORDER )
    {
      if(bl_start[dir_mask[id_mask][XDIR]] == X_BORDER || bl_start[dir_mask[id_mask][XDIR]] == 0 ) id = proc_dn(id,1);
      else if(bl_start[dir_mask[id_mask][XDIR]] == X || bl_start[dir_mask[id_mask][XDIR]] == X + X_BORDER ) id = proc_up(id,1);
    }
  
  if(bl_width[dir_mask[id_mask][YDIR]] == Y_BORDER )
    {
      if(bl_start[dir_mask[id_mask][YDIR]] == Y_BORDER || bl_start[dir_mask[id_mask][YDIR]] == 0 ) id = proc_dn(id,2);
      else if(bl_start[dir_mask[id_mask][YDIR]] == Y || bl_start[dir_mask[id_mask][YDIR]] == Y + Y_BORDER ) id = proc_up(id,2);
    }
  
  if(bl_width[dir_mask[id_mask][ZDIR]] == Z_BORDER )
    {
      if(bl_start[dir_mask[id_mask][ZDIR]] == Z_BORDER || bl_start[dir_mask[id_mask][ZDIR]] == 0 ) id = proc_dn(id,3);
      else if(bl_start[dir_mask[id_mask][ZDIR]] == Z || bl_start[dir_mask[id_mask][ZDIR]] == Z + Z_BORDER ) id = proc_up(id,3);
    }
      
  border[index_border].id_proc = id;
  
  lprintf("GEOMETRY",REPORTLVL,"id related %d\n",id); 

  return retval;
}

static void set_border_pointer(int actualn,int match_point)
{
  if( map_overlexi2id[actualn]==-1 ) {
    map_overlexi2id[actualn]=index_position;
    map_id2overlexi[index_position]=actualn;
    memory[index_position]=index_position;
    lprintf("GEOMETRY",REPORTLVL,"set_border_pointer  memory[%d] =  %d \n",memory[index_position],index_position); 
    index_position++;
  }
  else    
    if(!match_point){
      map_id2overlexi[index_position]=actualn;
      memory[index_position]=map_overlexi2id[actualn];
      index_position++;
    }
  
}





static void set_border_l3(int eotype){
  int pdir=0,fdir,sdir;
  int bl_start[4];
  int bl_width[4];
  int tmp_start[4];
  int tmp_width[4];
  int incr[4];
  
  fdir = 3;
  sdir = 2;
  
  if(NP_X>1) pdir++;
  if(NP_Y>1) pdir++;
  if(NP_Z>1) pdir++;
  if(NP_T>1) pdir++;


  int ind[4][2]={{T_BORDER,T},
  		 {X_BORDER,X},
  		 {Y_BORDER,Y},
  		 {Z_BORDER,Z}};

#define reset_tmp  tmp_start[0]=T_BORDER;  tmp_start[1]=X_BORDER; tmp_start[2]=Y_BORDER; tmp_start[3]=Z_BORDER; tmp_width[0]=T; tmp_width[1]=X; tmp_width[2]=Y; tmp_width[3]=Z;
 
  /* int tmp_ind=index_position; */

  /* tmp_ind=index_position; */
  
  lprintf("GEOMETRY",REPORTLVL,"set_border_l3 After setting only the inner piece we have %d boders\n",index_border);
  /*(Level=3)*/
  if(pdir>0){
 
    reset_tmp;
    lprintf("GEOMETRY",REPORTLVL,"set_border_l3 first dir %c          second dir %c\n",dc[inv_mask[0][fdir]],dc[inv_mask[0][sdir]]);

    tmp_width[inv_mask[0][fdir]]=ind[inv_mask[0][fdir]][0];
    lprintf("GEOMETRY",REPORTLVL,"set_border_l3 cambio dir %d\n",inv_mask[0][fdir]);

    set_block(0,bl_start,incr,bl_width, 
    	      tmp_start[0],tmp_width[0],
    	      tmp_start[1],tmp_width[1],
    	      tmp_start[2],tmp_width[2],
    	      tmp_start[3],tmp_width[3]);

    lprintf("GEOMETRY",REPORTLVL,"set_border_l3 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[0][0]],bl_start[dir_mask[0][1]],bl_start[dir_mask[0][2]],bl_start[dir_mask[0][3]]);
    lprintf("GEOMETRY",REPORTLVL,"set_border_l3 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[0][0]],bl_width[dir_mask[0][1]],bl_width[dir_mask[0][2]],bl_width[dir_mask[0][3]]);

    walk_on_lattice(0,eotype,3,index_zone,bl_start,incr,bl_width);
    
    if(pdir>1){
      
      reset_tmp;
      
      tmp_width[inv_mask[0][sdir]]=ind[inv_mask[0][sdir]][0];
      tmp_start[inv_mask[0][sdir]]=ind[inv_mask[0][sdir]][1];
      
      set_block(0,bl_start,incr,bl_width, 
		tmp_start[0],tmp_width[0],
		tmp_start[1],tmp_width[1],
		tmp_start[2],tmp_width[2],
		tmp_start[3],tmp_width[3]);
      
      lprintf("GEOMETRY",REPORTLVL,"set_border_l3 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[0][0]],bl_start[dir_mask[0][1]],bl_start[dir_mask[0][2]],bl_start[dir_mask[0][3]]);
      lprintf("GEOMETRY",REPORTLVL,"set_border_l3 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[0][0]],bl_width[dir_mask[0][1]],bl_width[dir_mask[0][2]],bl_width[dir_mask[0][3]]);
      
      walk_on_lattice(0,eotype,3,index_zone,bl_start,incr,bl_width);
    }

    
    reset_tmp;

    tmp_width[inv_mask[0][fdir]]=ind[inv_mask[0][fdir]][0];
    tmp_start[inv_mask[0][sdir]]=ind[inv_mask[0][sdir]][1];
    if(pdir==1) tmp_start[inv_mask[0][sdir]]--;
    tmp_start[inv_mask[0][fdir]]=ind[inv_mask[0][fdir]][1];
    
    set_block(0,bl_start,incr,bl_width, 
    	      tmp_start[0],tmp_width[0],
    	      tmp_start[1],tmp_width[1],
    	      tmp_start[2],tmp_width[2],
    	      tmp_start[3],tmp_width[3]);

    lprintf("GEOMETRY",REPORTLVL,"set_border_l3 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[0][0]],bl_start[dir_mask[0][1]],bl_start[dir_mask[0][2]],bl_start[dir_mask[0][3]]);
    lprintf("GEOMETRY",REPORTLVL,"set_border_l3 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[0][0]],bl_width[dir_mask[0][1]],bl_width[dir_mask[0][2]],bl_width[dir_mask[0][3]]);

    walk_on_lattice(0,eotype,3,index_zone,bl_start,incr,bl_width);


    if(pdir>1){
      
      reset_tmp;
      
      tmp_width[inv_mask[0][sdir]]=ind[inv_mask[0][sdir]][0];
      tmp_start[inv_mask[0][fdir]]=ind[inv_mask[0][fdir]][1];
      
      set_block(0,bl_start,incr,bl_width, 
		tmp_start[0],tmp_width[0],
		tmp_start[1],tmp_width[1],
		tmp_start[2],tmp_width[2],
		tmp_start[3],tmp_width[3]);
      
      lprintf("GEOMETRY",REPORTLVL,"set_border_l3 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[0][0]],bl_start[dir_mask[0][1]],bl_start[dir_mask[0][2]],bl_start[dir_mask[0][3]]);
      lprintf("GEOMETRY",REPORTLVL,"set_border_l3 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[0][0]],bl_width[dir_mask[0][1]],bl_width[dir_mask[0][2]],bl_width[dir_mask[0][3]]);
      
      walk_on_lattice(0,eotype,3,index_zone,bl_start,incr,bl_width);
    }

    if(pdir>2){
   
      reset_tmp;
      
      tmp_width[inv_mask[1][fdir]]=ind[inv_mask[1][fdir]][0];
      lprintf("GEOMETRY",REPORTLVL,"set_border_l3 cambio dir %d\n",inv_mask[1][fdir]);

      set_block(1,bl_start,incr,bl_width, 
		tmp_start[0],tmp_width[0],
		tmp_start[1],tmp_width[1],
		tmp_start[2],tmp_width[2],
		tmp_start[3],tmp_width[3]);

      lprintf("GEOMETRY",REPORTLVL,"set_border_l3 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[1][0]],bl_start[dir_mask[1][1]],bl_start[dir_mask[1][2]],bl_start[dir_mask[1][3]]);
      lprintf("GEOMETRY",REPORTLVL,"set_border_l3 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[1][0]],bl_width[dir_mask[1][1]],bl_width[dir_mask[1][2]],bl_width[dir_mask[1][3]]);

      walk_on_lattice(1,eotype,3,index_zone,bl_start,incr,bl_width);
      lprintf("GEOMETRY",REPORTLVL,"set_border_l3 \n");

      if(pdir>3){
	
	reset_tmp;
	
	tmp_width[inv_mask[1][sdir]]=ind[inv_mask[1][sdir]][0];
	tmp_start[inv_mask[1][sdir]]=ind[inv_mask[1][sdir]][1];
      
	set_block(1,bl_start,incr,bl_width, 
		  tmp_start[0],tmp_width[0],
		  tmp_start[1],tmp_width[1],
		  tmp_start[2],tmp_width[2],
		  tmp_start[3],tmp_width[3]);
	
	lprintf("GEOMETRY",REPORTLVL,"set_border_l3 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[1][0]],bl_start[dir_mask[1][1]],bl_start[dir_mask[1][2]],bl_start[dir_mask[1][3]]);
	lprintf("GEOMETRY",REPORTLVL,"set_border_l3 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[1][0]],bl_width[dir_mask[1][1]],bl_width[dir_mask[1][2]],bl_width[dir_mask[1][3]]);
	
	walk_on_lattice(1,eotype,3,index_zone,bl_start,incr,bl_width);
	lprintf("GEOMETRY",REPORTLVL,"set_border_l3 \n");
      }
      
      
      reset_tmp;
      
      tmp_width[inv_mask[1][fdir]]=ind[inv_mask[1][fdir]][0];
      tmp_start[inv_mask[1][sdir]]=ind[inv_mask[1][sdir]][1];
      if(pdir==3) tmp_start[inv_mask[1][sdir]]--;
      tmp_start[inv_mask[1][fdir]]=ind[inv_mask[1][fdir]][1];
        
      set_block(1,bl_start,incr,bl_width, 
		tmp_start[0],tmp_width[0],
		tmp_start[1],tmp_width[1],
		tmp_start[2],tmp_width[2],
		tmp_start[3],tmp_width[3]);

      lprintf("GEOMETRY",REPORTLVL,"set_border_l3 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[1][0]],bl_start[dir_mask[1][1]],bl_start[dir_mask[1][2]],bl_start[dir_mask[1][3]]);
      lprintf("GEOMETRY",REPORTLVL,"set_border_l3 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[1][0]],bl_width[dir_mask[1][1]],bl_width[dir_mask[1][2]],bl_width[dir_mask[1][3]]);

      walk_on_lattice(1,eotype,3,index_zone,bl_start,incr,bl_width);
      lprintf("GEOMETRY",REPORTLVL,"set_border_l3 \n");

      if(pdir>3){
	
	reset_tmp;
	
	tmp_width[inv_mask[1][sdir]]=ind[inv_mask[1][sdir]][0];
	tmp_start[inv_mask[1][fdir]]=ind[inv_mask[1][fdir]][1];
	
	set_block(1,bl_start,incr,bl_width, 
		  tmp_start[0],tmp_width[0],
		  tmp_start[1],tmp_width[1],
		  tmp_start[2],tmp_width[2],
		  tmp_start[3],tmp_width[3]);
	
	lprintf("GEOMETRY",REPORTLVL,"set_border_l3 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[1][0]],bl_start[dir_mask[1][1]],bl_start[dir_mask[1][2]],bl_start[dir_mask[1][3]]);
	lprintf("GEOMETRY",REPORTLVL,"set_border_l3 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[1][0]],bl_width[dir_mask[1][1]],bl_width[dir_mask[1][2]],bl_width[dir_mask[1][3]]);
	
	walk_on_lattice(1,eotype,3,index_zone,bl_start,incr,bl_width);
	lprintf("GEOMETRY",REPORTLVL,"set_border_l3 \n");
      }
    }
  }
  index_zone++;
}

static void set_border_l2(int eotype){
  int pdir=0,fdir,sdir;
  int bl_start[4];
  int bl_width[4];
  int tmp_start[4];
  int tmp_width[4];
  int incr[4];

  if(first_l2_border==-1) first_l2_border=index_border;



  fdir = 3;
  sdir = 2;

  if(NP_X>1) pdir++;
  if(NP_Y>1) pdir++;
  if(NP_Z>1) pdir++;
  if(NP_T>1) pdir++;


  int ind[4][2]={{T_BORDER,T},
  		 {X_BORDER,X},
  		 {Y_BORDER,Y},
  		 {Z_BORDER,Z}};


  if(pdir==0) return;
  /* (Level=2) */
  lprintf("GEOMETRY",REPORTLVL,"set_border_l2 After level 3 we have %d boders\n",index_border);
  if(pdir>1){    
    
    reset_tmp;
    
    tmp_width[inv_mask[0][fdir]]=ind[inv_mask[0][fdir]][0];
    tmp_width[inv_mask[0][sdir]]=ind[inv_mask[0][sdir]][0];
    
    set_block(0,bl_start,incr,bl_width, 
    	      tmp_start[0],tmp_width[0],
    	      tmp_start[1],tmp_width[1],
    	      tmp_start[2],tmp_width[2],
    	      tmp_start[3],tmp_width[3]);
    
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[0][0]],bl_start[dir_mask[0][1]],bl_start[dir_mask[0][2]],bl_start[dir_mask[0][3]]);
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[0][0]],bl_width[dir_mask[0][1]],bl_width[dir_mask[0][2]],bl_width[dir_mask[0][3]]);

    walk_on_lattice(0,eotype,2,index_zone,bl_start,incr,bl_width);

    reset_tmp;
    
    tmp_start[inv_mask[0][fdir]]=ind[inv_mask[0][fdir]][1];
    tmp_width[inv_mask[0][fdir]]=ind[inv_mask[0][fdir]][0];
    tmp_width[inv_mask[0][sdir]]=ind[inv_mask[0][sdir]][0];
    
    set_block(0,bl_start,incr,bl_width, 
    	      tmp_start[0],tmp_width[0],
    	      tmp_start[1],tmp_width[1],
    	      tmp_start[2],tmp_width[2],
    	      tmp_start[3],tmp_width[3]);
    
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[0][0]],bl_start[dir_mask[0][1]],bl_start[dir_mask[0][2]],bl_start[dir_mask[0][3]]);
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[0][0]],bl_width[dir_mask[0][1]],bl_width[dir_mask[0][2]],bl_width[dir_mask[0][3]]);

    walk_on_lattice(0,eotype,2,index_zone,bl_start,incr,bl_width);

    reset_tmp;
    
    tmp_start[inv_mask[0][fdir]]=ind[inv_mask[0][fdir]][1];
    tmp_start[inv_mask[0][sdir]]=ind[inv_mask[0][sdir]][1];
    tmp_width[inv_mask[0][fdir]]=ind[inv_mask[0][fdir]][0];
    tmp_width[inv_mask[0][sdir]]=ind[inv_mask[0][sdir]][0];
    
    set_block(0,bl_start,incr,bl_width, 
    	      tmp_start[0],tmp_width[0],
    	      tmp_start[1],tmp_width[1],
    	      tmp_start[2],tmp_width[2],
    	      tmp_start[3],tmp_width[3]);
    
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[0][0]],bl_start[dir_mask[0][1]],bl_start[dir_mask[0][2]],bl_start[dir_mask[0][3]]);
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[0][0]],bl_width[dir_mask[0][1]],bl_width[dir_mask[0][2]],bl_width[dir_mask[0][3]]);

    walk_on_lattice(0,eotype,2,index_zone,bl_start,incr,bl_width);
    
    reset_tmp;
  
    tmp_start[inv_mask[0][sdir]]=ind[inv_mask[0][sdir]][1];
    tmp_width[inv_mask[0][fdir]]=ind[inv_mask[0][fdir]][0];
    tmp_width[inv_mask[0][sdir]]=ind[inv_mask[0][sdir]][0];
    
    set_block(0,bl_start,incr,bl_width, 
    	      tmp_start[0],tmp_width[0],
    	      tmp_start[1],tmp_width[1],
    	      tmp_start[2],tmp_width[2],
    	      tmp_start[3],tmp_width[3]);
    
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[0][0]],bl_start[dir_mask[0][1]],bl_start[dir_mask[0][2]],bl_start[dir_mask[0][3]]);
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[0][0]],bl_width[dir_mask[0][1]],bl_width[dir_mask[0][2]],bl_width[dir_mask[0][3]]);

    walk_on_lattice(0,eotype,2,index_zone,bl_start,incr,bl_width);

  }

  if(pdir>3){    
    
    reset_tmp;
    
    tmp_width[inv_mask[1][fdir]]=ind[inv_mask[1][fdir]][0];
    tmp_width[inv_mask[1][sdir]]=ind[inv_mask[1][sdir]][0];
    
    set_block(1,bl_start,incr,bl_width, 
    	      tmp_start[0],tmp_width[0],
    	      tmp_start[1],tmp_width[1],
    	      tmp_start[2],tmp_width[2],
    	      tmp_start[3],tmp_width[3]);
    
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[1][0]],bl_start[dir_mask[1][1]],bl_start[dir_mask[1][2]],bl_start[dir_mask[1][3]]);
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[1][0]],bl_width[dir_mask[1][1]],bl_width[dir_mask[1][2]],bl_width[dir_mask[1][3]]);

    walk_on_lattice(1,eotype,2,index_zone,bl_start,incr,bl_width);

    reset_tmp;
    
    tmp_start[inv_mask[1][fdir]]=ind[inv_mask[1][fdir]][1];
    tmp_width[inv_mask[1][fdir]]=ind[inv_mask[1][fdir]][0];
    tmp_width[inv_mask[1][sdir]]=ind[inv_mask[1][sdir]][0];
    
    set_block(1,bl_start,incr,bl_width, 
    	      tmp_start[0],tmp_width[0],
    	      tmp_start[1],tmp_width[1],
    	      tmp_start[2],tmp_width[2],
    	      tmp_start[3],tmp_width[3]);
    
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[1][0]],bl_start[dir_mask[1][1]],bl_start[dir_mask[1][2]],bl_start[dir_mask[1][3]]);
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[1][0]],bl_width[dir_mask[1][1]],bl_width[dir_mask[1][2]],bl_width[dir_mask[1][3]]);

    walk_on_lattice(1,eotype,2,index_zone,bl_start,incr,bl_width);

    reset_tmp;
    
    tmp_start[inv_mask[1][fdir]]=ind[inv_mask[1][fdir]][1];
    tmp_start[inv_mask[1][sdir]]=ind[inv_mask[1][sdir]][1];
    tmp_width[inv_mask[1][fdir]]=ind[inv_mask[1][fdir]][0];
    tmp_width[inv_mask[1][sdir]]=ind[inv_mask[1][sdir]][0];
    
    set_block(1,bl_start,incr,bl_width, 
    	      tmp_start[0],tmp_width[0],
    	      tmp_start[1],tmp_width[1],
    	      tmp_start[2],tmp_width[2],
    	      tmp_start[3],tmp_width[3]);
    
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[1][0]],bl_start[dir_mask[1][1]],bl_start[dir_mask[1][2]],bl_start[dir_mask[1][3]]);
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[1][0]],bl_width[dir_mask[1][1]],bl_width[dir_mask[1][2]],bl_width[dir_mask[1][3]]);

    walk_on_lattice(1,eotype,2,index_zone,bl_start,incr,bl_width);
    
    reset_tmp;
  
    tmp_start[inv_mask[1][sdir]]=ind[inv_mask[1][sdir]][1];
    tmp_width[inv_mask[1][fdir]]=ind[inv_mask[1][fdir]][0];
    tmp_width[inv_mask[1][sdir]]=ind[inv_mask[1][sdir]][0];
    
    set_block(1,bl_start,incr,bl_width, 
    	      tmp_start[0],tmp_width[0],
    	      tmp_start[1],tmp_width[1],
    	      tmp_start[2],tmp_width[2],
    	      tmp_start[3],tmp_width[3]);
    
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[1][0]],bl_start[dir_mask[1][1]],bl_start[dir_mask[1][2]],bl_start[dir_mask[1][3]]);
    lprintf("GEOMETRY",REPORTLVL,"set_border_l2 width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[1][0]],bl_width[dir_mask[1][1]],bl_width[dir_mask[1][2]],bl_width[dir_mask[1][3]]);
       
    walk_on_lattice(1,eotype,2,index_zone,bl_start,incr,bl_width);

  } 
  int id0,id1,id2,id3;
  int wd0,wd1,wd2,wd3;
  lprintf("GEOMETRY",REPORTLVL,"set_border_l2 After the intersecting level 2 we have %d boders\n",index_border);


  for(id0=0;id0<2;id0++){
    tmp_start[inv_mask[0][0]]=ind[inv_mask[0][0]][id0];
    for(id1=0;id1<2;id1++){
      tmp_start[inv_mask[0][1]]=ind[inv_mask[0][1]][id1];
      for(id2=0;id2<2;id2++){
	tmp_start[inv_mask[0][2]]=ind[inv_mask[0][2]][id2];
	for(id3=0;id3<2;id3++){
	  tmp_start[inv_mask[0][3]]=ind[inv_mask[0][3]][id3];
	  
	  for(wd0=0;wd0<2;wd0++){
	    tmp_width[inv_mask[0][0]]=ind[inv_mask[0][0]][wd0];
	    for(wd1=0;wd1<2;wd1++){
	      tmp_width[inv_mask[0][1]]=ind[inv_mask[0][1]][wd1];
	      for(wd2=0;wd2<2;wd2++){
		tmp_width[inv_mask[0][2]]=ind[inv_mask[0][2]][wd2];
		for(wd3=0;wd3<2;wd3++){
		  tmp_width[inv_mask[0][3]]=ind[inv_mask[0][3]][wd3];
		  if(wd0+wd1+wd2+wd3==2){
		    
		    set_block(0,bl_start,incr,bl_width, 
			      tmp_start[0],tmp_width[0],
			      tmp_start[1],tmp_width[1],
			      tmp_start[2],tmp_width[2],
			      tmp_start[3],tmp_width[3]);
		    
		    lprintf("GEOMETRY",REPORTLVL,"set_border_l2_auto start T[%d] X[%d] Y[%d] Z[%d] \n",bl_start[dir_mask[0][0]],bl_start[dir_mask[0][1]],bl_start[dir_mask[0][2]],bl_start[dir_mask[0][3]]);
		    lprintf("GEOMETRY",REPORTLVL,"set_border_l2_auto width T[%d] X[%d] Y[%d] Z[%d] \n",bl_width[dir_mask[0][0]],bl_width[dir_mask[0][1]],bl_width[dir_mask[0][2]],bl_width[dir_mask[0][3]]);
		    
		    walk_on_lattice(0,eotype,2,index_zone,bl_start,incr,bl_width);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  index_zone++;
 
  lprintf("GEOMETRY",REPORTLVL,"set_border In the end we have %d borders\n",index_border+1);

}


inline int abs(int n){
  return (n>0)?n:-n;
}




static void buffer_size(int i,int *tmp_width, int* tmp_start){
  tmp_width[0]= abs(border[i].b_width[0]);
  tmp_width[1]= abs(border[i].b_width[1]);
  tmp_width[2]= abs(border[i].b_width[2]);
  tmp_width[3]= abs(border[i].b_width[3]);
  lprintf("GEOMETRY",REPORTLVL,"fix_buffer width[%d]=%d\n",0,tmp_width[0]);
  lprintf("GEOMETRY",REPORTLVL,"fix_buffer width[%d]=%d\n",1,tmp_width[1]);
  lprintf("GEOMETRY",REPORTLVL,"fix_buffer width[%d]=%d\n",2,tmp_width[2]);
  lprintf("GEOMETRY",REPORTLVL,"fix_buffer width[%d]=%d\n",3,tmp_width[2]);
  
  
  if(tmp_width[TDIR] == T_BORDER )
    {
      if( border[i].b_start[TDIR] == T_BORDER )
	tmp_start[TDIR]= T+T_BORDER;
      else if( border[i].b_start[TDIR] == T)
	tmp_start[TDIR]= 0;
    }
  else
    tmp_start[TDIR] = border[i].b_start[TDIR];
  
  if(tmp_width[XDIR] == X_BORDER )
    {
      if( border[i].b_start[XDIR] == X_BORDER )
	tmp_start[XDIR]= X+X_BORDER;
      else if( border[i].b_start[XDIR] == X)
	tmp_start[XDIR]= 0;
    }
  else
    tmp_start[XDIR] = border[i].b_start[XDIR];
  
  if(tmp_width[YDIR] == Y_BORDER )
    {
      if( border[i].b_start[YDIR] == Y_BORDER )
	tmp_start[YDIR]= Y+Y_BORDER;
      else if( border[i].b_start[YDIR] == Y)
	tmp_start[YDIR]= 0;
    }
  else
    tmp_start[YDIR] = border[i].b_start[YDIR];
  
  if(tmp_width[ZDIR] == Z_BORDER )
    {
      if( border[i].b_start[ZDIR] == Z_BORDER )
	tmp_start[ZDIR]= Z+Z_BORDER;
      else if( border[i].b_start[ZDIR] == Z)
	tmp_start[ZDIR]= 0;
    }
  else
    tmp_start[ZDIR] = border[i].b_start[ZDIR];
}
  
  
static void fix_buffer(int zone_id)
{	
  lprintf("GEOMETRY",REPORTLVL,"fix_buffer start on zone %d\n",zone_id);

  int i,done_border=index_border+1,id_mask;
  int tmp_start[4],bl_start[4],bl_width[4];
  int tmp_width[4];
  int incr[4];

  for(i=0;i<done_border;i++){
    if(border[i].id_zone==zone_id){
      if( inner[0]==i || inner[1]==i ) continue;
      id_mask=border[i].id_mask;
      
      buffer_size(i,tmp_width,tmp_start);
      
      set_block(id_mask,bl_start,incr,bl_width,
		tmp_start[0],tmp_width[0],
		tmp_start[1],tmp_width[1],
		tmp_start[2],tmp_width[2],
		tmp_start[3],tmp_width[3]);


      /* walk_on_lattice(id_mask,buffer_sign,border[i].level,-zone_id,bl_start,incr,bl_width); */
      walk_on_lattice(id_mask,border[i].eo_type,border[i].level,-zone_id,bl_start,incr,bl_width); 

    }
  }
}

static int correnspondig_buffer(int * bf , int i){
  int tmp_start[4],bl_start[4],bl_width[4];
  int tmp_width[4],done_border=index_border+1,checksum[8],buffer_sign;
  int id_mask=border[i].id_mask,found;
  int incr[4],j;

  buffer_size(i,tmp_width,tmp_start);
  
  set_block(id_mask,bl_start,incr,bl_width,
	    tmp_start[0],tmp_width[0],
	    tmp_start[1],tmp_width[1],
	    tmp_start[2],tmp_width[2],
	    tmp_start[3],tmp_width[3]);
  
  border_checksum(checksum,id_mask,bl_start,bl_width);

  
  if((border[i].checksum[0]+border[i].checksum[1]+border[i].checksum[2]+border[i].checksum[3])%2==(checksum[0]+checksum[1]+checksum[2]+checksum[3])%2)
    buffer_sign= border[i].eo_type;
  else
    buffer_sign= (border[i].eo_type+1)%2;



  //int myfound[10];

  found=0;
  for(j=0;j<done_border;j++){
    if(border[j].id_zone<0 &&
       border[j].checksum[0] == checksum[0] &&
       border[j].checksum[1] == checksum[1] &&
       border[j].checksum[2] == checksum[2] &&
       border[j].checksum[3] == checksum[3] &&
       border[j].checksum[4] == checksum[4] &&
       border[j].checksum[5] == checksum[5] &&
       border[j].checksum[6] == checksum[6] &&
       border[j].checksum[7] == checksum[7]){
	found++;

	//myfound[found]=j;
	if(border[j].eo_type == buffer_sign )
	  bf[0]=j;
	else
	  bf[1]=j;
      }
  }

  error(found!=2,1,"correnspondig_buffer [geometry_mpi_eo.c]",
	"found more than one buffer corrensponding to one border");

  return 0;
}


 


static void  fix_next_neightbours()
{
  int x0,x1,x2,x3,ix;
  for (x3=0;x3<Z+2*Z_BORDER;x3++)
    for (x2=0;x2<Y+2*Y_BORDER;x2++)
      for (x1=0;x1<X+2*X_BORDER;x1++)
  	for (x0=0;x0<T+2*T_BORDER;x0++)
	  {
	    
	    ix = map_overlexi2id[local_index(x0,x1,x2,x3)];
	    
	    ipt(x0-T_BORDER,x1-X_BORDER,x2-Y_BORDER,x3-Z_BORDER)=ix ;

	    if(ix != -1)
	      {
	    	iup(ix,0)=map_overlexi2id[local_index(x0+1,x1,x2,x3)];
	    	idn(ix,0)=map_overlexi2id[local_index(x0-1,x1,x2,x3)];
	    	iup(ix,1)=map_overlexi2id[local_index(x0,x1+1,x2,x3)];
	    	idn(ix,1)=map_overlexi2id[local_index(x0,x1-1,x2,x3)];
	    	iup(ix,2)=map_overlexi2id[local_index(x0,x1,x2+1,x3)];
	    	idn(ix,2)=map_overlexi2id[local_index(x0,x1,x2-1,x3)];
	    	iup(ix,3)=map_overlexi2id[local_index(x0,x1,x2,x3+1)];
	    	idn(ix,3)=map_overlexi2id[local_index(x0,x1,x2,x3-1)];
	      }
   
	   	    
	  }
  
}





static void set_memory_order()
{
  int i0,i1;

  int master_piece_number=0;
  int master_piece_number_e=0;
  int master_piece_number_o=0;
  int *master_piece_start_list=malloc(index_position*sizeof(int));
  int *master_piece_length_list=malloc(index_position*sizeof(int));
  int *master_piece_start_list_e=malloc(index_position*sizeof(int));
  int *master_piece_length_list_e=malloc(index_position*sizeof(int));
  int *master_piece_start_list_o=malloc(index_position*sizeof(int));
  int *master_piece_length_list_o=malloc(index_position*sizeof(int));

  int local_master_piece_number_e=0;
  int local_master_piece_number_o=0;


  int test_master=1;
  int test_master_e=1;
  int test_master_o=1;

  int * copy_list_from=malloc(index_position*sizeof(int));
  int * copy_list_to=malloc(index_position*sizeof(int));
  int * copy_list_len=malloc(index_position*sizeof(int));
  
  int * copy_list_from_e=malloc(index_position*sizeof(int));
  int * copy_list_to_e=malloc(index_position*sizeof(int));
  int * copy_list_len_e=malloc(index_position*sizeof(int));

  int * copy_list_from_o=malloc(index_position*sizeof(int));
  int * copy_list_to_o=malloc(index_position*sizeof(int));
  int * copy_list_len_o=malloc(index_position*sizeof(int));
  
  
  int copy_piece_number=0;
  int copy_piece_number_e=0;
  int copy_piece_number_o=0;

  int local_copy_piece_number_e=0;
  int local_copy_piece_number_o=0;

  int index_start=0;
  int first_entrance=true;
  int first_odd=true;

  int shift;

  glattice.gsize_gauge=border[index_border].index_end;

  glat_even.total_gauge_master_pieces=-1;
  glat_odd.total_gauge_master_pieces=-1;
  glat_even.ncopies_gauge=-1;
  glat_odd.ncopies_gauge=-1;
 


  /*Border & Buffer*/
  for(i0=0;i0<=index_border;i0++){
    /* lprintf("GEOMETRY",REPORTLVL,"Border :%d ----------------------------\n",i0); */

    if(border[i0].level>2 && border[i0].size>0) {

      glattice.gsize_spinor=border[i0].index_end;

      if(border[i0].eo_type==even ){   
	/* lprintf("GEOMETRY",REPORTLVL,"Even L4L3 ----------------------------\n",i0); */
	glat_even.gsize_spinor=border[i0].index_end;


     	for (i1=index_start;i1<border[i0].index_end;i1++){

	  if(memory[i1]==i1) {
	    /*master even*/
	    
	    if(test_master_e==0 && ! zone_border(i1) ) {
	      master_piece_length_list_e[master_piece_number_e-1]++;
	      /* lprintf("GEOMETRY",REPORTLVL,"\t\t | \n"); */
	    }
	    else {
	      master_piece_start_list_e[master_piece_number_e]=i1;
	      master_piece_length_list_e[master_piece_number_e]=1;
	      master_piece_number_e++;
	      /* lprintf("GEOMETRY",REPORTLVL,"\t\t---\n"); */
	    }
	    
	    test_master_e=0;
	    
	  }
	  else {
	    /*copy even*/
	    
	    if(memory[i1]==memory[i1-1]+1 && test_master_e==1) {
	      copy_list_len_e[copy_piece_number_e-1]++;
	      /* lprintf("GEOMETRY",REPORTLVL,"\t\t\t\t * \n"); */
	    }
	    else {
	      copy_list_from_e[copy_piece_number_e]=memory[i1];
	      copy_list_to_e[copy_piece_number_e]=i1;
	      copy_list_len_e[copy_piece_number_e]=1;
	      copy_piece_number_e++;
	      /* lprintf("GEOMETRY",REPORTLVL,"\t\t\t\t~~~\n"); */
	    }
	    test_master_e=1;
	    
	  }
	  /* lprintf("GEOMETRY",REPORTLVL,"Considero %d\n",i1); */
	}
      }

      if(border[i0].eo_type==odd){   
	/* lprintf("GEOMETRY",REPORTLVL,"Odd L4L3 ----------------------------\n"); */
	for (i1=index_start;i1<border[i0].index_end;i1++){
	  if(memory[i1]==i1) {
	    /*master*/
	    
	    if(test_master_o==0 && ! zone_border(i1) ) {
	      master_piece_length_list_o[master_piece_number_o-1]++;
	    }
	    else {
	      master_piece_start_list_o[master_piece_number_o]=i1;
	      master_piece_length_list_o[master_piece_number_o]=1;
	      master_piece_number_o++;
	    }
	    
	    test_master_o=0;
	    
	  }
	  else {
	    /*copy*/
	    
	    if(memory[i1]==memory[i1-1]+1 && test_master_o==1) {
	      copy_list_len_o[copy_piece_number_o-1]++;
	    }
	    else {
	      copy_list_from_o[copy_piece_number_o]=memory[i1];
	      copy_list_to_o[copy_piece_number_o]=i1;
	      copy_list_len_o[copy_piece_number_o]=1;
	      copy_piece_number_o++;
	    }
	    test_master_o=1;
	    
	  }
	}
   	/* lprintf("GEOMETRY",REPORTLVL,"Out  ----------------------------\n",i0);  */
      }     
    }

    if(border[i0].level>2 && border[i0].eo_type==odd && border[i0].id_zone >=0 ) {
      local_copy_piece_number_o=copy_piece_number_o;      
      local_master_piece_number_o=master_piece_number_o;
    } 

    if(border[i0].level>2 && border[i0].eo_type==even && border[i0].id_zone >=0 ) {
      local_copy_piece_number_e=copy_piece_number_e;      
      local_master_piece_number_e=master_piece_number_e;
    } 

    if(border[i0].eo_type==odd && first_odd && border[i0].id_zone>=0) {
      glat_odd.master_shift=border[i0].index_start;
      first_odd=false;
    }
  
    /* lprintf("GEOMETRY",REPORTLVL,"L2 ----------------------------\n",i0); */
    if((i0==first_l2_border || i0==index_border ) && first_entrance){
      glat_odd.gsize_spinor=glattice.gsize_spinor-glat_even.gsize_spinor;
      

      glat_even.master_shift=0;	
      glattice.master_shift=0;
      
      glat_odd.copy_shift=copy_piece_number_e;
      glat_even.copy_shift=0;	
      glattice.copy_shift=0;
      
      glat_odd.local_master_pieces=local_master_piece_number_o;
      glat_even.local_master_pieces=local_master_piece_number_e;
      glattice.local_master_pieces=local_master_piece_number_o+local_master_piece_number_e;

      int i_inner=0;
      int i_inner_e=0;
      int i_inner_o=0;

      if(border[inner[even]].size!=0){
	master_piece_start_list[i_inner]=master_piece_start_list_e[0];
	master_piece_length_list[i_inner]=master_piece_length_list_e[0];
	i_inner++; i_inner_e++;
      }

      if(border[inner[odd]].size!=0){
	master_piece_start_list[i_inner]=master_piece_start_list_o[0];
	master_piece_length_list[i_inner]=master_piece_length_list_o[0];
	i_inner++; i_inner_o++;
      }
      

      
      shift=i_inner_e;
      /* lprintf("GEOMETRY",REPORTLVL,"1 shift=%d\n",shift); */
      for(i1=shift;i1<local_master_piece_number_e;i1++){
	master_piece_start_list[i_inner]=master_piece_start_list_e[i_inner_e];
	master_piece_length_list[i_inner]=master_piece_length_list_e[i_inner_e];
	i_inner_e++;	
	i_inner++;
      }
      shift=i_inner_o;
      /* lprintf("GEOMETRY",REPORTLVL,"2 shift=%d\n",shift); */
      for(i1=shift;i1<local_master_piece_number_o;i1++){
	master_piece_start_list[i_inner]=master_piece_start_list_o[i_inner_o];
	master_piece_length_list[i_inner]=master_piece_length_list_o[i_inner_o];
	i_inner_o++;	
	i_inner++;
      }
      /* lprintf("GEOMETRY",REPORTLVL,"3 shift=%d\n",shift); */
      for(i1=local_master_piece_number_e;i1<master_piece_number_e;i1++){
	master_piece_start_list[i_inner]=master_piece_start_list_e[i1];
	master_piece_length_list[i_inner]=master_piece_length_list_e[i1];
	i_inner++;
      }
      shift=master_piece_number_e;
      /* lprintf("GEOMETRY",REPORTLVL,"4 shift=%d\n",shift); */
      for(i1=local_master_piece_number_o;i1<master_piece_number_o;i1++){
	master_piece_start_list[i_inner]=master_piece_start_list_o[i1];
	master_piece_length_list[i_inner]=master_piece_length_list_o[i1];
	i_inner++;
      }


      for(i1=0;i1<local_copy_piece_number_e;i1++){
	copy_list_from[i1]=copy_list_from_e[i1];
	copy_list_to[i1]=copy_list_to_e[i1];
	copy_list_len[i1]=copy_list_len_e[i1];
	
      }
      shift=local_copy_piece_number_e;
      for(i1=0;i1<local_copy_piece_number_o;i1++){
	copy_list_from[i1+shift]=copy_list_from_o[i1];
	copy_list_to[i1+shift]=copy_list_to_o[i1];
	copy_list_len[i1+shift]=copy_list_len_o[i1];
      }
      shift=local_copy_piece_number_o;
      for(i1=local_copy_piece_number_e;i1<copy_piece_number_e;i1++){
	copy_list_from[i1+shift]=copy_list_from_e[i1];
	copy_list_to[i1+shift]=copy_list_to_e[i1];
	copy_list_len[i1+shift]=copy_list_len_e[i1];
      }
      shift=copy_piece_number_e;
      for(i1=local_copy_piece_number_o;i1<copy_piece_number_o;i1++){
	copy_list_from[i1+shift]=copy_list_from_o[i1];
	copy_list_to[i1+shift]=copy_list_to_o[i1];
	copy_list_len[i1+shift]=copy_list_len_o[i1];
      }
      
      master_piece_number=master_piece_number_e+master_piece_number_o;
      copy_piece_number=copy_piece_number_e+copy_piece_number_o;
      
      glattice.total_spinor_master_pieces=master_piece_number;
      glattice.ncopies_spinor=copy_piece_number;
      
      lprintf("GEOMETRY",REPORTLVL,"glattice MPN%d CPN%d \n", glattice.total_spinor_master_pieces,glattice.ncopies_spinor);
      first_entrance=false;
    }

    if(border[i0].level==2){ 
      
      

      for (i1=index_start;i1<border[i0].index_end;i1++){

	if(memory[i1]==i1) {
	  /*master*/
	  
	  if(test_master==0 && ! zone_border(i1) ) {
	    master_piece_length_list[master_piece_number-1]++;
	  }
	  else {
	    master_piece_start_list[master_piece_number]=i1;
	    master_piece_length_list[master_piece_number]=1;
	    master_piece_number++;
	  }
	  
	  test_master=0;
	  
	}
	else {
	  /*copy*/
	  
	  if(memory[i1]==memory[i1-1]+1 && test_master==1) {
	    copy_list_len[copy_piece_number-1]++;
	  }
	  else {
	    copy_list_from[copy_piece_number]=memory[i1];
	    copy_list_to[copy_piece_number]=i1;
	    copy_list_len[copy_piece_number]=1;
	    copy_piece_number++;
	  }
	  test_master=1;
	  
	}
	
      }
    }

    if(border[i0].index_end>index_start) index_start=border[i0].index_end;   
  }

  if(master_piece_number!=0)
    {
      glattice.master_start=malloc(master_piece_number*sizeof(unsigned int));
      error((glattice.master_start==NULL),1,"set_memory_order [geometry_mpi.c]",
  	    "Cannot allocate memory");
      
      glattice.master_end=malloc(master_piece_number*sizeof(unsigned int));
      error((glattice.master_end==NULL),1,"set_memory_order [geometry_mpi.c]",
  	    "Cannot allocate memory");
      
      for (i1=0;i1<master_piece_number;i1++)
  	{
  	  glattice.master_start[i1]=master_piece_start_list[i1];
  	  glattice.master_end[i1]=master_piece_start_list[i1]+master_piece_length_list[i1]-1;
	  lprintf("[GEOMETRY]",REPORTLVL,"START END i:%d start:%d end:%d\n",i1,glattice.master_start[i1],glattice.master_end[i1]);
 	}

      glattice.total_gauge_master_pieces=master_piece_number;
      
    }

  if(master_piece_number_e!=0)
    {
      glat_even.master_start=malloc(master_piece_number_e*sizeof(unsigned int));
      error((glat_even.master_start==NULL),1,"set_memory_order [geometry_mpi.c]",
  	    "Cannot allocate memory");
      
      glat_even.master_end=malloc(master_piece_number_e*sizeof(unsigned int));
      error((glat_even.master_end==NULL),1,"set_memory_order [geometry_mpi.c]",
  	    "Cannot allocate memory");
      
      for (i1=0;i1<master_piece_number_e;i1++)
  	{
  	  glat_even.master_start[i1]=master_piece_start_list_e[i1];
  	  glat_even.master_end[i1]=master_piece_start_list_e[i1]+master_piece_length_list_e[i1]-1;
  	}
      
      glat_even.total_spinor_master_pieces=master_piece_number_e;
      
    }
  
  if(master_piece_number_o!=0)
    {
      glat_odd.master_start=malloc(master_piece_number_e*sizeof(unsigned int));
      error((glat_odd.master_start==NULL),1,"set_memory_order [geometry_mpi.c]",
  	    "Cannot allocate memory");

      glat_odd.master_end=malloc(master_piece_number_o*sizeof(unsigned int));
      error((glat_odd.master_end==NULL),1,"set_memory_order [geometry_mpi.c]",
  	    "Cannot allocate memory");

      for (i1=0;i1<master_piece_number_o;i1++)
  	{
  	  glat_odd.master_start[i1]=master_piece_start_list_o[i1];
  	  glat_odd.master_end[i1]=master_piece_start_list_o[i1]+master_piece_length_list_o[i1]-1;
  	}

      glat_odd.total_spinor_master_pieces=master_piece_number_o;
 
    }

  if(copy_piece_number!=0)
    {
      glattice.copy_from=malloc(copy_piece_number*sizeof(int));
      error((glattice.copy_from==NULL),1,"set_memory_order [geometry_mpi.c]",
	    "Cannot allocate memory");
    
      glattice.copy_to=malloc(copy_piece_number*sizeof(int));
      error((glattice.copy_to==NULL),1,"set_memory_order [geometry_mpi.c]",
	    "Cannot allocate memory");
    
      glattice.copy_len=malloc(copy_piece_number*sizeof(int));
      error((glattice.copy_len==NULL),1,"set_memory_order [geometry_mpi.c]",
	    "Cannot allocate memory");
    
  
      for (i1=0;i1<copy_piece_number;i1++)
	{
	  glattice.copy_from[i1]=copy_list_from[i1];
	  glattice.copy_to[i1]=copy_list_to[i1];
	  glattice.copy_len[i1]=copy_list_len[i1];
	  /* lprintf("GEOMETRY",REPORTLVL,"Copy From:%d Len:%d To:%d\n",copy_list_from[i1],copy_list_len[i1],copy_list_to[i1]); */
	  /* if(i1==copy_piece_number_o+copy_piece_number_e-1 ) lprintf("GEOMETRY",REPORTLVL,"----------------------------\n"); */
	}
      glattice.ncopies_gauge=copy_piece_number;
    }

  if(copy_piece_number_e!=0)
    {
      glat_even.copy_from=malloc(copy_piece_number_e*sizeof(int));
      error((glat_even.copy_from==NULL),1,"set_memory_order [geometry_mpi.c]",
	    "Cannot allocate memory");
    
      glat_even.copy_to=malloc(copy_piece_number_e*sizeof(int));
      error((glat_even.copy_to==NULL),1,"set_memory_order [geometry_mpi.c]",
	    "Cannot allocate memory");
    
      glat_even.copy_len=malloc(copy_piece_number_e*sizeof(int));
      error((glat_even.copy_len==NULL),1,"set_memory_order [geometry_mpi.c]",
	    "Cannot allocate memory");
    
    
      for (i1=0;i1<copy_piece_number_e;i1++)
	{
	  glat_even.copy_from[i1]=copy_list_from_e[i1];
	  glat_even.copy_to[i1]=copy_list_to_e[i1];
	  glat_even.copy_len[i1]=copy_list_len_e[i1];
	}
      glat_even.ncopies_spinor=copy_piece_number_e;
    }

  if(copy_piece_number_o!=0)
    {
      glat_odd.copy_from=malloc(copy_piece_number_o*sizeof(int));
      error((glat_odd.copy_from==NULL),1,"set_memory_order [geometry_mpi.c]",
	    "Cannot allocate memory");
    
      glat_odd.copy_to=malloc(copy_piece_number_o*sizeof(int));
      error((glat_odd.copy_to==NULL),1,"set_memory_order [geometry_mpi.c]",
	    "Cannot allocate memory");
    
      glat_odd.copy_len=malloc(copy_piece_number_o*sizeof(int));
      error((glat_odd.copy_len==NULL),1,"set_memory_order [geometry_mpi.c]",
	    "Cannot allocate memory");
    
    
      for (i1=0;i1<copy_piece_number_o;i1++)
	{
	  glat_odd.copy_from[i1]=copy_list_from_o[i1];
	  glat_odd.copy_to[i1]=copy_list_to_o[i1];
	  glat_odd.copy_len[i1]=copy_list_len_o[i1];
	}
      glat_odd.ncopies_spinor=copy_piece_number_o;
    }
   
  free(copy_list_from);
  free(copy_list_to);
  free(copy_list_len);
  free(copy_list_from_e);
  free(copy_list_to_e);
  free(copy_list_len_e);
  free(copy_list_from_o);
  free(copy_list_to_o);
  free(copy_list_len_o);
  free(master_piece_start_list);
  free(master_piece_length_list);
  free(master_piece_start_list_e);
  free(master_piece_length_list_e);
  free(master_piece_start_list_o);
  free(master_piece_length_list_o);
 
  int fbuf=0;
  int sfbuf_e=-1;
  int sfbuf_o=-1;
  int done_border=index_border+1,j,i,j1 ;
  int bf[2] = {0,0};

  for(i=0;i<done_border;i++){
    if(border[i].id_zone>=0 &&  inner[0]!=i && inner[1]!=i && border[i].size > 0){
      correnspondig_buffer(bf,i);
      j=bf[0];
      j1=bf[1];

      (glattice.sbuf_to_proc)[fbuf] = border[i].id_proc;
      (glattice.sbuf_start)[fbuf] = border[i].index_start;
      (glattice.rbuf_from_proc)[fbuf] = border[j].id_proc;
      (glattice.sbuf_len)[fbuf] = border[i].size;
      if(border[i].eo_type== border[j].eo_type){
	(glattice.rbuf_start)[fbuf] = border[j].index_start;
	(glattice.rbuf_len)[fbuf] = border[j].size;
      }
      else 
	{
	  (glattice.rbuf_start)[fbuf] = border[j1].index_start;
	  (glattice.rbuf_len)[fbuf] = border[j1].size;
	}
      

      fbuf++;
      
      
      if(border[i].level==3){
      	if(border[i].eo_type==even){
      	  sfbuf_e++;
      	  (glat_even.sbuf_to_proc)[sfbuf_e] = border[i].id_proc;
      	  (glat_even.sbuf_start)[sfbuf_e] = border[i].index_start;
      	  (glat_even.sbuf_len)[sfbuf_e] =  border[i].size;
      	  lprintf("GEOMETRY",REPORTLVL,"glattice_eo (%d) E send to %d start %d(%d)",sfbuf_e, (glat_even.sbuf_to_proc)[sfbuf_e],(glat_even.sbuf_start)[sfbuf_e],(glat_even.sbuf_len)[sfbuf_e]);
	  
	  (glat_even.rbuf_from_proc)[sfbuf_e] = border[j].id_proc;
      	  if(border[j].eo_type==even){
      	    (glat_even.rbuf_start)[sfbuf_e] = border[j].index_start;
      	    (glat_even.rbuf_len)[sfbuf_e] = border[j].size;
	  }
	  else{
      	    (glat_even.rbuf_start)[sfbuf_e] = border[j1].index_start;
	    (glat_even.rbuf_len)[sfbuf_e] = border[j1].size;
	  }
	  
	  lprintf("GEOMETRY",REPORTLVL," -> E receive from %d start %d(%d) \n", (glat_even.rbuf_from_proc)[sfbuf_e],(glat_even.rbuf_start)[sfbuf_e],(glat_even.rbuf_len)[sfbuf_e]);
      	}
	
      	if(border[i].eo_type==odd){
      	  sfbuf_o++;
      	  (glat_odd.sbuf_to_proc)[sfbuf_o] = border[i].id_proc;
      	  (glat_odd.sbuf_start)[sfbuf_o] = border[i].index_start;
      	  (glat_odd.sbuf_len)[sfbuf_o] = border[i].size;
      	  lprintf("GEOMETRY",REPORTLVL,"glattice_eo (%d) O send to %d start %d(%d)",sfbuf_o, (glat_odd.sbuf_to_proc)[sfbuf_o],(glat_odd.sbuf_start)[sfbuf_o],(glat_odd.sbuf_len)[sfbuf_o]);
	  
	  (glat_odd.rbuf_from_proc)[sfbuf_o] = border[j].id_proc;	  
      	  if(border[j].eo_type==odd){
      	    (glat_odd.rbuf_start)[sfbuf_o] = border[j].index_start;
	    (glat_odd.rbuf_len)[sfbuf_o] = border[j].size;
	  }
	  else{
      	    (glat_odd.rbuf_start)[sfbuf_o] = border[j1].index_start;
	    (glat_odd.rbuf_len)[sfbuf_o] = border[j1].size;
	  }
	  lprintf("GEOMETRY",REPORTLVL," -> O receive from %d start %d(%d) \n", (glat_odd.rbuf_from_proc)[sfbuf_o],(glat_odd.rbuf_start)[sfbuf_o],(glat_odd.rbuf_len)[sfbuf_o]);
	  
	}
      }
      
      
      
      lprintf("GEOMETRY",REPORTLVL,"\n");
      


      
    }
  }
}


void geometry_mpi_eo(void){

  geometry_mpi_init();

  set_inner(even);

  set_border_l3(even);
  
  fix_buffer(index_zone-1);
  
  set_inner(odd);

  set_border_l3(odd);

  fix_buffer(index_zone-1);

  set_border_l2(even);

  set_border_l2(odd); 
  
  fix_buffer(index_zone-2);
  
  fix_buffer(index_zone-1);
  
  //  logger_setlevel(0,11);

  set_memory_order();
  
  geometry_mem_alloc();
  
  fix_next_neightbours();
  
  geometry_mpi_finalize(); 

      
}












