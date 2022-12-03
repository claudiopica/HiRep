#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include "geometry.h"
#include "global.h"
#include "error.h"
#include "logger.h"
#include "utils.h"

// #define _INCLUDE_UP_UP_L2

char *const LOGTAG = "GEOMETRY DEFINE";

#define _PRINT_BYTE "%c%c%c%c%c%c%c%c"
#define _BINARY(byte)  \
  (byte & 0x80 ? '1' : '0'), \
  (byte & 0x40 ? '1' : '0'), \
  (byte & 0x20 ? '1' : '0'), \
  (byte & 0x10 ? '1' : '0'), \
  (byte & 0x08 ? '1' : '0'), \
  (byte & 0x04 ? '1' : '0'), \
  (byte & 0x02 ? '1' : '0'), \
  (byte & 0x01 ? '1' : '0') 


static inline int lexi(int b0, int b1, int b2, int b3, int x0, int x1, int x2, int x3) {
    return x0 + x1*b0 + x2*b1*b0 + x3*b2*b1*b0;
}

//find a suitale block size for the lenght L
//input b: requested block size
static int set_block_size(int L, int b) {
#ifndef NDEBUG
    error(b>0,1, "GEOMETRY: "__FILE__,"Incorrect inner blocking dimensions used");
#endif
    if (L%b==0) return b;
    //keep decreasing b to find a good value
    while(b>2) {
        b--;
        if (L%b==0) return b;
    }
    return 1;
}

// b[0-3] = block size in each direction
// X[0-3] = total size of the box
// x[0-3] = coordinate of site in the box
// sign = the parity of orgin site of the box
// return an index from 0 to (X0*X1*X2*X3)-1
static int index_blocked(
    int b0, int b1, int b2, int b3, 
    int X0, int X1, int X2, int X3, 
    int x0, int x1, int x2, int x3,
    int sign
    ) 
{
    int xb0,xb1,xb2,xb3; //coordinates inside the block
    int xn0,xn1,xn2,xn3; //coordinates of the block

    //Check that block sizes divive box size
    b0=set_block_size(X0,b0);
    b1=set_block_size(X1,b1);
    b2=set_block_size(X2,b2);
    b3=set_block_size(X3,b3);

    //compute coordinates inside block and block coordinates
    xb0=x0%b0; xn0=x0/b0;
    xb1=x1%b1; xn1=x1/b1;
    xb2=x2%b2; xn2=x2/b2;
    xb3=x3%b3; xn3=x3/b3;   
    // lprintf("INDEX",1,"BLK size used=[%d,%d,%d,%d]\n",b0,b1,b2,b3);
    // lprintf("INDEX",1,"coord =[%d,%d,%d,%d]\n",x0,x1,x2,x3);
    // lprintf("INDEX",1,"inner =[%d,%d,%d,%d]\n",xb0,xb1,xb2,xb3);
    // lprintf("INDEX",1,"blkco =[%d,%d,%d,%d]\n",xn0,xn1,xn2,xn3);

    int parity = (x0+x1+x2+x3+sign)%2; //beware of global parity
    int blk_vol = b0*b1*b2*b3;

    //this counts within a block
    //goes from 0 to blk_vol-1
    int inner_index = lexi(b0, b1, b2, b3, xb0, xb1, xb2, xb3);

    //this counts the blocks
    //goes from 0 to num_blocks-1
    int block_index = lexi(X0/b0, X1/b1, X2/b2, X3/b3, xn0, xn1, xn2, xn3);

    int idx = inner_index + block_index*blk_vol;

    // lprintf("INDEX",1,"sign=%d idx=%d, inner idx=%d, blk idx=%d, blk vol=%d\n",parity, idx, inner_index, block_index, blk_vol);

    //correct for parity
    //all even sites go in the first half
    //NB: box_volume can be odd
    idx = idx / 2; 
    const int box_vol= X0*X1*X2*X3;
    // lprintf("TEST",1,"BOX=vol=%d ODD start=%d\n",box_vol,box_vol/2 + (box_vol&1));
    if (parity) idx += box_vol/2 + (box_vol&(1^sign)); //beware for overflows!!!

    return idx;
}

static void gd_free_mem(geometry_descriptor *gd) {
    if(gd->rbuf_len) {
        free(gd->rbuf_len);
        gd->rbuf_len         = NULL;
        gd->sbuf_len         = NULL;
        gd->sbuf_to_proc     = NULL;
        gd->sbuf_start       = NULL;
        gd->rbuf_from_proc   = NULL;
        gd->rbuf_start       = NULL;
    }    
}

// if N>0 this allocates memory for the send/rec index buffers
// otherwise sets them to NULL
static void gd_alloc_mem(geometry_descriptor *gd, int N) {
    int *buf = NULL;
    if (N>0) {
        buf = malloc(6 * N * sizeof(int));
        error((buf == NULL), 1, __func__ , "Cannot allocate memory");
    }

    gd->rbuf_len         = buf; if (buf) buf += N;
    gd->sbuf_len         = buf; if (buf) buf += N;
    gd->sbuf_to_proc     = buf; if (buf) buf += N;
    gd->sbuf_start       = buf; if (buf) buf += N;
    gd->rbuf_from_proc   = buf; if (buf) buf += N;
    gd->rbuf_start       = buf;

}

// in this geometry we don't need to copy pieces of geometry
// this functions sets the number of copies to zero e NULLs the arrays for copies 
static void gd_set_copy(geometry_descriptor *gd) {
    gd->ncopies_gauge = 0;
    gd->ncopies_spinor = 0;
    gd->copy_from = NULL;
    gd->copy_to = NULL;
    gd->copy_len = NULL;
    gd->copy_shift = 0;
}

// TODO: this should be in geometry.h and geometry descriptor should contain it
// enum to define geometry type
enum gd_type {GLOBAL, EVEN, ODD};

// in this geometry we always have only 2 (global) or 1 (even/odd) master piece(s)
static void gd_set_master_pieces(geometry_descriptor *gd, enum gd_type gd_t) {
    //this contains the start and end of the master pieces
    // 0 -> global geometry
    // 1 -> even geometry
    // 2 -> odd geometry
    static int master_start[3], master_end[3];

    const int even_volume = VOLUME/2 + ((PSIGN==0)? (VOLUME&1) : 0);
    master_start[0] = 0; master_end[0] = VOLUME-1; //TODO: remove the -1
    master_start[1] = 0; master_end[1] = even_volume-1; //TODO: remove the -1
    master_start[2] = master_end[1]+1; master_end[2] = master_end[0]; //TODO: remove the -1

    gd->inner_master_pieces = 1;
    gd->local_master_pieces = 1;
    gd->total_gauge_master_pieces = 1;
    gd->total_spinor_master_pieces = 1;

    if (gd_t==GLOBAL) {
        gd->master_start = &master_start[0];
        gd->master_end = &master_end[0];
        gd->master_shift = 0;
    } else if (gd_t==EVEN) {
        gd->master_start = &master_start[1];
        gd->master_end = &master_end[1];
        gd->master_shift = 0;
    } else if (gd_t==ODD) {
        gd->master_start = &master_start[2];
        gd->master_end = &master_end[2];
        gd->master_shift = even_volume;
    } else {
        error(1,1, "GEOMETRY: "__FILE__,"Incorrect geometry type used in gd_set_master_pieces");
    }
}

static void gd_set_size(geometry_descriptor *gd, enum gd_type gd_t){
    const int halo_volume = T_EXT*X_EXT*Y_EXT*Z_EXT-VOLUME;
    const int l3_border = (X * Y * Z * T_BORDER +
                     X * Y * Z_BORDER * T +
                     X * Y_BORDER * Z * T +
                     X_BORDER * Y * Z * T);
    //NB: the total number of sites in the local lattice can be odd
    //e.g. if T=X=Y=Z=3 (or any odd number)
    //we can't assume that the even and odd part of the lattice have the same number of points!
    if (gd_t==GLOBAL) {
        gd->gsize_gauge = VOLUME + 2*halo_volume; //master points + 2x halo for buffers TODO: does it makes sense if non global?
        gd->gsize_spinor = VOLUME + 2*l3_border; // master points + 2x L3(=dimension 3) borders for buffers
    } else if (gd_t==EVEN) {
        int even_volume = VOLUME/2 + ((PSIGN==0)? (VOLUME&1) : 0); //add one if VOLUME IS ODD and parity of local lattice is even
        gd->gsize_gauge = even_volume + halo_volume; //master points + 2x halo for buffers TODO: does it makes sense if non global?
        gd->gsize_spinor = even_volume + l3_border; // master points + 2x L3(=dimension 3) borders for buffers
    } else if (gd_t==ODD) {
        int odd_volume = VOLUME/2 + ((PSIGN==0)? 0 : (VOLUME&1) ); //add one if VOLUME IS ODD and parity of local lattice is even
        gd->gsize_gauge = odd_volume + halo_volume; //master points + 2x halo for buffers TODO: does it makes sense if non global?
        gd->gsize_spinor = odd_volume + l3_border; // master points + 2x L3(=dimension 3) borders for buffers
    } else {
        error(1,1, "GEOMETRY: "__FILE__,"Incorrect geometry type used in gd_set_size");
    }
    // lprintf("TEST",1,"HALO VOL=%d L3 VOL=%d\n",halo_volume,l3_border);
    // lprintf("TEST",1,"GSIZE_GAUGE=%d GSIZE_SPINOR=%d\n",gd->gsize_gauge,gd->gsize_spinor);
}

//allocate memory for global indexes
// ipt, iup, idn, imask
static void index_alloc() {
    const int ext_volume = T_EXT*X_EXT*Y_EXT*Z_EXT;
    const int N = (1+2*4)*ext_volume; // 1 for ipt + 4 dir x 2 for iup and idn
    ipt = malloc(N * sizeof(*ipt));
    error((ipt == NULL), 1, __func__ , "Cannot allocate memory for ipt, iup, idn");
    iup = ipt + ext_volume;
    idn = iup + (4*ext_volume);

    imask = malloc(ext_volume * sizeof(*imask));
    error((imask == NULL), 1, __func__ , "Cannot allocate memory for imask");
}

static void index_free() {
    if (ipt) {
        free(ipt);
        ipt = NULL;
        iup = NULL;
        idn = NULL;
    }
    if (imask) {
        free(imask);
        imask = NULL;
    }
}

static void geometry_index_init() {
    //number of parallel directions
    const int NPAR = (NP_T>1 ? 1 : 0) + (NP_X>1? 1 : 0) + (NP_Y>1 ? 1 : 0) + (NP_Z>1 ? 1 : 0);
    // Dimension 3 borders. These are for both gauge fields and spinors
    const int L3_BORDER = 2 * NPAR;
    // Dimension 2 borders. These are for the gauge fields
#ifdef _INCLUDE_UP_UP_L2
    const int L2_BORDER = 4*NPAR*(NPAR-1)/2; //TODO: change from 4->3. the last one is not needed actually
#else
    const int L2_BORDER = 3*NPAR*(NPAR-1)/2; //TODO: change from 4->3. the last one is not needed actually
#endif

    lprintf(LOGTAG,200,"L3 BORDERS=%d  L2 BORDERS=%d\n", L3_BORDER, L2_BORDER);

    //set up global lattice    
    glattice.nbuffers_gauge = L3_BORDER+L2_BORDER;
    glattice.nbuffers_spinor = L3_BORDER;
    gd_alloc_mem(&glattice, L3_BORDER+L2_BORDER);
    gd_set_copy(&glattice);
    gd_set_master_pieces(&glattice, GLOBAL);
    gd_set_size(&glattice, GLOBAL);

    //set up even lattice
    glat_even.nbuffers_gauge = L3_BORDER+L2_BORDER; //TODO: why is this not allowed? it would make no difference for this geometry
    glat_even.nbuffers_spinor = L3_BORDER;
    gd_alloc_mem(&glat_even, L3_BORDER+L2_BORDER); //TODO: do we want only L3 here?
    gd_set_copy(&glat_even);
    gd_set_master_pieces(&glat_even, EVEN);
    gd_set_size(&glat_even, EVEN);
    
    //set up odd lattice
    glat_odd.nbuffers_gauge = L3_BORDER+L2_BORDER;
    glat_odd.nbuffers_spinor = L3_BORDER;
    gd_alloc_mem(&glat_odd, L3_BORDER+L2_BORDER); //TODO: do we want only L3 here?
    gd_set_copy(&glat_odd);
    gd_set_master_pieces(&glat_odd, ODD);
    gd_set_size(&glat_odd, ODD);

    //alloc memory for global indexes
    index_alloc();

}

static void geometry_index_free() {
    gd_free_mem(&glattice);
    gd_free_mem(&glat_even);
    gd_free_mem(&glat_odd);
    index_free();
}

enum MaskState {
    T_UP_MASK = (1u << 0),
    T_DN_MASK = (1u << 1),
    X_UP_MASK = (1u << 2),
    X_DN_MASK = (1u << 3),
    Y_UP_MASK = (1u << 4),
    Y_DN_MASK = (1u << 5),
    Z_UP_MASK = (1u << 6),
    Z_DN_MASK = (1u << 7),
    FULL_MASK = (1u << 8)-1
};

char invertMask(char mask) {
    return mask ^ FULL_MASK;
}

enum box_type {
    INNER = 4,
    L3 = 3,
    L2 = 2
};

//  ----h
//  ....|  NB: the h[4] is not in the box,   
//  ....|  i.e. coordinates needs to be l[]<= p[] <h[] 
//  l...|
typedef struct _box_t {
    int l[4]; // the lower left corner, the box base point (e.g. in the extended lattice)
    int h[4]; // the upper right corner
    int base_index;
    char mask; //tells if the box is a border, e.g. if T_UP_MASK is set the box is in top T border of the extended lattice
    enum box_type type; // tell the type of the box, just a convenience for testing
    struct _box_t *sendBox; //if this is a border corresponding to a Recv buffer, this is the box to copy data from, i.e. corresponding to the Send buffer
    struct _box_t *next; // link to next box. NULL if last
} box_t;

// This is to keep a list of boxes to fill out the field buffers
static box_t geometryBoxes;

box_t* duplicateBox(box_t *B) {
    box_t *b=malloc(sizeof(box_t));
    error((b == NULL), 1, __func__ , "Cannot allocate memory for box");
    *b=*B; //make copy
    return b;
}

void insertBox(box_t *B, box_t *A){
    A->next = B->next;
    B->next = A;
}

//split a box in two with an hyperplane at level ortogonal to dir
//returns two boxes: the original one is changed, plus a second one
//the hyperplane must cut the box 
//the larger box in the split direction is returned in the orginal argument
//dir = 0, 1, 2, 3
// Returns: 
// 0 => no split
// 1 => the box was split, a new box was appended to the list
int splitBox(box_t *B, int level, int dir) {
    // l . . (level) . . . h
    // e.g. l=0 (level=1) 2 3 ... h=T_EXT-1
    // => dl = 1
    // => dh = T_EXT-1-1
    //the hyperplane level becomes the h of one box and the l of the other box
    const int dl=level-(B->l[dir]);
    const int dh=(B->h[dir])-level;

    //hyperplane outsite the box => do nothing
    if (dl<1 || dh<1) return 0;
    
    //make copy of original box
    box_t *B2 = duplicateBox(B);

    //the larger box becomes B
    if (dl<dh) {
        B2->h[dir]=level;
        B->l[dir]=level;
    } else {
        B->h[dir]=level;
        B2->l[dir]=level;
    }

    //append new box to list
    insertBox(B, B2);

    return 1;
}

//split a list of Boxes
//assumes the list is NULL terminated!
void splitBoxList(box_t *B, int level, int dir) {
    while(B) {
        int isSplit = splitBox(B, level, dir);
        B=B->next;
        //skip next if it was just splitted
        if (isSplit) B=B->next; 
    }
}

//returns how many new boxes were created
int splitBorder(box_t *B, int border_size, int dir) {
    int i = splitBox(B, B->l[dir]+border_size, dir);
    i += splitBox(B, B->h[dir]-border_size, dir);
    return i;
}

void splitBorderList(box_t *B, int border_size, int dir) {
    while(B) {
        int newBoxes = splitBorder(B, border_size, dir);
        B=B->next;
        while(newBoxes--) B=B->next;    
    }
}

//Free memory for boxex in the list, 
//but not the first one which is passed
void freeBoxList(box_t *B) {
    while(B->next) {
        box_t *next = B->next->next;
        if (B->next->sendBox) { free(B->next->sendBox); B->next->sendBox=NULL; } 
        free(B->next);
        B->next = next;
    }
}

void freeGeometryBoxes(){
    freeBoxList(&geometryBoxes);
}

void printBox(box_t *B, int level) {
    lprintf(LOGTAG, level, "[%2d,%2d,%2d,%2d]..[%2d,%2d,%2d,%2d] sizes=[%2d,%2d,%2d,%2d] baseIndex=%9d",
        B->l[0],B->l[1],B->l[2],B->l[3],
        B->h[0],B->h[1],B->h[2],B->h[3],
        B->h[0]-B->l[0],B->h[1]-B->l[1],B->h[2]-B->l[2],B->h[3]-B->l[3],
        B->base_index
    );
    if (B->sendBox) {
        box_t *S=B->sendBox;
        lprintf(LOGTAG, level, " Send Box [%2d,%2d,%2d,%2d]..[%2d,%2d,%2d,%2d] sizes=[%2d,%2d,%2d,%2d]",
        S->l[0],S->l[1],S->l[2],S->l[3],
        S->h[0],S->h[1],S->h[2],S->h[3],
        S->h[0]-S->l[0],S->h[1]-S->l[1],S->h[2]-S->l[2],S->h[3]-S->l[3]
        );
    }
    lprintf(LOGTAG, 1, "\n");
}

void printBoxList(box_t *B, int level){
    int i=0;
    while(B) {
        i++;
        printBox(B, level);
        B=B->next;
    }
    lprintf(LOGTAG,level,"Number of boxes = %d\n", i);
}

int boxVolume(box_t *B){
    return (B->h[0]-B->l[0])*(B->h[1]-B->l[1])*(B->h[2]-B->l[2])*(B->h[3]-B->l[3]);

}
void setBoxListBaseIndex(box_t *B) {
    int i=B->base_index;
    while(B) {
        B->base_index = i;
        i += boxVolume(B);
        B=B->next;
    }
}

//calculate if the base point of the box is even (0) or odd (1)
//the box coordinates refers to the extended lattice
int boxParity(box_t *B) {
    return (B->l[0]+B->l[1]+B->l[2]+B->l[3]+PSIGN+T_BORDER+X_BORDER+Y_BORDER+Z_BORDER)%2; //remember global parity!
}
int safe_ipt_ext(int x0_ext, int x1_ext, int x2_ext, int x3_ext) {
    const int ix = ipt_ext(safe_mod(x0_ext,T_EXT),safe_mod(x1_ext,X_EXT),safe_mod(x2_ext,Y_EXT),safe_mod(x3_ext,Z_EXT)); 
#ifndef NDEBUG
    if (ix<0 || ix>(T_EXT*X_EXT*Y_EXT*Z_EXT-1)) { 
        lprintf("GEOMETRY",1,"Error3 ix=%d site [%d, %d, %d, %d]\n", ix,safe_mod(x0_ext,T_EXT),safe_mod(x1_ext,X_EXT),safe_mod(x2_ext,Y_EXT),safe_mod(x3_ext,Z_EXT));
    }
#endif
    return ix;
}

int nearProc(char mask) {
    int mycoord[4] = {COORD[0], COORD[1], COORD[2], COORD[3]};
    if (mask & T_UP_MASK) { mycoord[0]++; }
    if (mask & T_DN_MASK) { mycoord[0]--; }
    if (mask & X_UP_MASK) { mycoord[1]++; }
    if (mask & X_DN_MASK) { mycoord[1]--; }
    if (mask & Y_UP_MASK) { mycoord[2]++; }
    if (mask & Y_DN_MASK) { mycoord[2]--; }
    if (mask & Z_UP_MASK) { mycoord[3]++; }
    if (mask & Z_DN_MASK) { mycoord[3]--; }
    return proc_id(mycoord);
}

box_t *makeLocalBox(){
    box_t B = { .l={T_BORDER,X_BORDER,Y_BORDER,Z_BORDER}, 
                .h={T+T_BORDER,X+X_BORDER,Y+Y_BORDER,Z+Z_BORDER}, 
                .base_index=0, 
                .mask=0, 
                .type=INNER,
                .sendBox=NULL,
                .next=NULL
              };
    return duplicateBox(&B);
}

//this function create a new box starting from a box B
//and a mask of directions
//the box returned corrsponds to a border of the box in the 
//directions indicated by the mask
box_t *makeBorderBox(char mask) {
        box_t *b=makeLocalBox(); 
        b->mask=mask; b->type = 4;
        if (mask & T_UP_MASK) { b->l[0] = T+T_BORDER; b->h[0] = T_EXT;    b->type--;} else     
        if (mask & T_DN_MASK) { b->l[0] = 0;          b->h[0] = T_BORDER; b->type--;}     
        if (mask & X_UP_MASK) { b->l[1] = X+X_BORDER; b->h[1] = X_EXT;    b->type--;} else
        if (mask & X_DN_MASK) { b->l[1] = 0;          b->h[1] = X_BORDER; b->type--;}     
        if (mask & Y_UP_MASK) { b->l[2] = Y+Y_BORDER; b->h[2] = Y_EXT;    b->type--;} else
        if (mask & Y_DN_MASK) { b->l[2] = 0;          b->h[2] = Y_BORDER; b->type--;}     
        if (mask & Z_UP_MASK) { b->l[3] = Z+Z_BORDER; b->h[3] = Z_EXT;    b->type--;} else
        if (mask & Z_DN_MASK) { b->l[3] = 0;          b->h[3] = Z_BORDER; b->type--;} 
        //TODO: we don't check explicitly that only one of UP and DN bits per directions are set
        //make sendBox
        box_t *s=makeLocalBox(); 
        s->mask=FULL_MASK; s->type = 4; //we don't need these
        if (mask & T_UP_MASK) { s->h[0] = 2*T_BORDER; } else     
        if (mask & T_DN_MASK) { s->l[0] = T; }
        if (mask & X_UP_MASK) { s->h[1] = 2*X_BORDER; } else
        if (mask & X_DN_MASK) { s->l[1] = X; }     
        if (mask & Y_UP_MASK) { s->h[2] = 2*Y_BORDER; } else
        if (mask & Y_DN_MASK) { s->l[2] = Y; }     
        if (mask & Z_UP_MASK) { s->h[3] = 2*Z_BORDER; } else
        if (mask & Z_DN_MASK) { s->l[3] = Z; } 

        //link sendBox in b
        b->sendBox=s;
        return b;
}

// fill out ipt, iup, idn, imask
static void enumerate_lattice() {

    // const int size[4] = { T, X, Y, Z };
    // const int borders[4] = { T_BORDER, X_BORDER, Y_BORDER, Z_BORDER };
    // const int ext_size[4] = { T_EXT, X_EXT, Y_EXT, Z_EXT };
    const int ext_vol = T_EXT*X_EXT*Y_EXT*Z_EXT;
    const char up_masks[4] = { T_UP_MASK, X_UP_MASK, Y_UP_MASK, Z_UP_MASK };
    const char dn_masks[4] = { T_DN_MASK, X_DN_MASK, Y_DN_MASK, Z_DN_MASK };

    int n_parallel[4] = { -1 }; //the list of parallel directions
    int par_dirs=0;
    if (NP_T>1) { n_parallel[par_dirs] = 0; par_dirs++; }
    if (NP_X>1) { n_parallel[par_dirs] = 1; par_dirs++; }
    if (NP_Y>1) { n_parallel[par_dirs] = 2; par_dirs++; }
    if (NP_Z>1) { n_parallel[par_dirs] = 3; par_dirs++; }
    lprintf(LOGTAG,10,"Number of parallel directions: %d\n", par_dirs);

    //define local lattice box
    box_t B = *makeLocalBox();
    box_t *L=&B; //this is our list of boxes. first box is the local lattice B

    //add L3 borders
    for (int i=0; i<par_dirs; i++){
        int d1=n_parallel[i];
        box_t *b=NULL;
        char mask=0;
        //up dir
        mask=up_masks[d1];
        b=makeBorderBox(mask);
        insertBox(L,b); L=b; //insert new box in the list. Move tail of the list to new element
        //dn dir
        mask=dn_masks[d1];
        b=makeBorderBox(mask);
        insertBox(L,b); L=b; 
    }

    //add L2 borders
    for (int i=0; i<par_dirs; i++){
        for (int j=i+1; j<par_dirs; j++){
            int d1=n_parallel[i];
            int d2=n_parallel[j];
            box_t *b=NULL;
            char mask=0;
            //up dn dir
            mask=up_masks[d1] | dn_masks[d2];
            b=makeBorderBox(mask);
            insertBox(L,b); L=b; //insert new box in the list. Move tail of the list to new element
            //dn up dir
            mask=dn_masks[d1] | up_masks[d2];
            b=makeBorderBox(mask);
            insertBox(L,b); L=b; //insert new box in the list. Move tail of the list to new element
            //dn dn dir
            mask=dn_masks[d1] | dn_masks[d2];
            b=makeBorderBox(mask);
            insertBox(L,b); L=b; //insert new box in the list. Move tail of the list to new element
#ifdef _INCLUDE_UP_UP_L2
            //CP: WE DON'T NEED THIS PIECE. KEEP IT INLY FOR TESTING
            //up up dir
            mask=up_masks[d1] | up_masks[d2];
            b=makeBorderBox(mask);
            insertBox(L,b); L=b; //insert new box in the list. Move tail of the list to new element
#endif
        }
    }

#if 0
    //CP:OLD code. Delete.
    //define extended lattice box
    box_t B = {.l={0,0,0,0}, .h={T_EXT,X_EXT,Y_EXT,Z_EXT}, .base_index=0, .next=NULL};
    
    if (NP_T>1) { splitBorderList(&B, T_BORDER, 0); }
    if (NP_X>1) { splitBorderList(&B, X_BORDER, 1); }
    if (NP_Y>1) { splitBorderList(&B, Y_BORDER, 2); }
    if (NP_Z>1) { splitBorderList(&B, Z_BORDER, 3); }
#endif

    setBoxListBaseIndex(&B); //don't need to do this before hand, but for now we do it here for convenience
    printBoxList(&B, 10);

    const int INVALID_POINT = -1;
    //now the list of boxes is not going to cover the whole
    //extended lattice (e.g. L1 or L0 blocks or even some L2 might not there)
    //so we prefill the ipt, iup, idn array with INVALID_POINT
    // we prefill imask = 0 
    for (int ix=0, x0_ext=0;x0_ext<T_EXT;x0_ext++){
    for (int x1_ext=0;x1_ext<X_EXT;x1_ext++){
    for (int x2_ext=0;x2_ext<Y_EXT;x2_ext++){
    for (int x3_ext=0;x3_ext<Z_EXT;x3_ext++){
        //prefill tables with INVALID_POINT
        ipt_ext(x0_ext,x1_ext,x2_ext,x3_ext) = INVALID_POINT;
#ifndef NDEBUG
        if (ix<0 || ix>(ext_vol-1)) lprintf(LOGTAG,1,"ERROR ix=%d @ (%d, %d, %d, %d)\n", ix, x0_ext, x1_ext, x2_ext, x3_ext);
#endif
        imask(ix)=0;
        iup(ix,0) = INVALID_POINT;
        iup(ix,1) = INVALID_POINT;
        iup(ix,2) = INVALID_POINT;
        iup(ix,3) = INVALID_POINT;
        idn(ix,0) = INVALID_POINT;
        idn(ix,1) = INVALID_POINT;
        idn(ix,2) = INVALID_POINT;
        idn(ix,3) = INVALID_POINT;
        ix++;
    }}}}

    //enumerate points in each box
    //fill out ipt_ext
    box_t *b = &B;
    while(b) {
        const int sign = boxParity(b);
        for (int x0_ext=b->l[0];x0_ext<b->h[0];x0_ext++){
        for (int x1_ext=b->l[1];x1_ext<b->h[1];x1_ext++){
        for (int x2_ext=b->l[2];x2_ext<b->h[2];x2_ext++){
        for (int x3_ext=b->l[3];x3_ext<b->h[3];x3_ext++){
            int ix = b->base_index; 
            const int idx = index_blocked(BLK_T,BLK_X,BLK_Y,BLK_Z,
                                b->h[0]-b->l[0],b->h[1]-b->l[1],b->h[2]-b->l[2],b->h[3]-b->l[3],
                                x0_ext-b->l[0],x1_ext-b->l[1],x2_ext-b->l[2],x3_ext-b->l[3], sign);
            ix += idx;
#ifndef NDEBUG
            if (idx<0 || !(idx<boxVolume(b))) lprintf(LOGTAG,1,"ERROR in ipt: Volume=%d idx=%d  @ (%d, %d, %d, %d)\n", boxVolume(b), idx, x0_ext, x1_ext, x2_ext, x3_ext);
            if (ix<0 || ix>(ext_vol-1)) lprintf(LOGTAG,1,"ERROR2 ix=%d @ (%d, %d, %d, %d)\n", ix, x0_ext, x1_ext, x2_ext, x3_ext);
#endif
            ipt_ext(x0_ext,x1_ext,x2_ext,x3_ext) = ix;
        }}}}

        b=b->next;
    }

    //fill out iup, idn, imask from ipt_ext
    //note that some ix can be INVALID_POINT
    //moving from a valid point in a direction indicated by the mask
    //is always a valid point on the local lattice
    for (int x0_ext=0;x0_ext<T_EXT;x0_ext++){
    for (int x1_ext=0;x1_ext<X_EXT;x1_ext++){
    for (int x2_ext=0;x2_ext<Y_EXT;x2_ext++){
    for (int x3_ext=0;x3_ext<Z_EXT;x3_ext++){
        const int ix = ipt_ext(x0_ext,x1_ext,x2_ext,x3_ext);
        if (ix == INVALID_POINT) continue;
#ifndef NDEBUG
        if (ix<0 || ix>(ext_vol-1)) {
            lprintf("GEOMETRY",1,"Error ix=%d site [%d, %d, %d, %d]\n", ix,x0_ext,x1_ext,x2_ext,x3_ext);
        }
#endif

        char xmask=0; int iy=0;
        iy=safe_ipt_ext(x0_ext+1,x1_ext,x2_ext,x3_ext); iup(ix,0)=iy; if(iy>=0 && iy<VOLUME) xmask |= T_UP_MASK; 
        iy=safe_ipt_ext(x0_ext-1,x1_ext,x2_ext,x3_ext); idn(ix,0)=iy; if(iy>=0 && iy<VOLUME) xmask |= T_DN_MASK; 
        iy=safe_ipt_ext(x0_ext,x1_ext+1,x2_ext,x3_ext); iup(ix,1)=iy; if(iy>=0 && iy<VOLUME) xmask |= X_UP_MASK; 
        iy=safe_ipt_ext(x0_ext,x1_ext-1,x2_ext,x3_ext); idn(ix,1)=iy; if(iy>=0 && iy<VOLUME) xmask |= X_DN_MASK; 
        iy=safe_ipt_ext(x0_ext,x1_ext,x2_ext+1,x3_ext); iup(ix,2)=iy; if(iy>=0 && iy<VOLUME) xmask |= Y_UP_MASK; 
        iy=safe_ipt_ext(x0_ext,x1_ext,x2_ext-1,x3_ext); idn(ix,2)=iy; if(iy>=0 && iy<VOLUME) xmask |= Y_DN_MASK; 
        iy=safe_ipt_ext(x0_ext,x1_ext,x2_ext,x3_ext+1); iup(ix,3)=iy; if(iy>=0 && iy<VOLUME) xmask |= Z_UP_MASK; 
        iy=safe_ipt_ext(x0_ext,x1_ext,x2_ext,x3_ext-1); idn(ix,3)=iy; if(iy>=0 && iy<VOLUME) xmask |= Z_DN_MASK; 

        imask(ix)=xmask;


    }}}}

    //set rbuf_len, sbuf_len, rbuf_from_proc, rbuf_start, sbuf_to_proc, sbuf_start;
    int sbuf_base=ext_vol; //send buffers are at the end of the extended lattice. If we move to a global send buffer this should be set to 0.
    int ib=0; //index of buffer in geometry descriptor
    L=B.next; //first border buffer. L3 are first, then L2
    while(L) {
        //each box contains an even and odd part
        //these are two separate buffers for spinors
        const int vol = boxVolume(L);
        const int sign = boxParity(L);
        const int even_vol = vol/2 + ((sign==0)? (vol&1) : 0); //add one if parity is even and vol is odd
        const int odd_vol = vol - even_vol;
        const int recv_proc=nearProc(L->mask);
        const int send_proc=nearProc(invertMask(L->mask));
        lprintf(LOGTAG,20,"ID=%d Base=%d vol=%d Level=%d mask=" _PRINT_BYTE " Recv proc=%d Send proc=%d\n", ib,L->base_index, vol, (int)(L->type),_BINARY(L->mask), recv_proc, send_proc);
        
        glattice.rbuf_len[ib]=vol;
        glattice.sbuf_len[ib]=vol;
        glattice.rbuf_start[ib]=L->base_index;
        glattice.sbuf_start[ib]=sbuf_base;
        glattice.rbuf_from_proc[ib]=recv_proc;
        glattice.sbuf_to_proc[ib]=send_proc;

        glat_even.rbuf_len[ib]=even_vol;
        glat_even.sbuf_len[ib]=even_vol;
        glat_even.rbuf_start[ib]=L->base_index;
        glat_even.sbuf_start[ib]=sbuf_base;
        glat_even.rbuf_from_proc[ib]=recv_proc;
        glat_even.sbuf_to_proc[ib]=send_proc;

        glat_odd.rbuf_len[ib]=odd_vol;
        glat_odd.sbuf_len[ib]=odd_vol;
        glat_odd.rbuf_start[ib]=L->base_index+even_vol;
        glat_odd.sbuf_start[ib]=sbuf_base+even_vol;
        glat_odd.rbuf_from_proc[ib]=recv_proc;
        glat_odd.sbuf_to_proc[ib]=send_proc;

        sbuf_base += vol;
        L=L->next; ib++;
    }

    //save geometry boxes
    geometryBoxes = B;
    atexit(&freeGeometryBoxes);

}

//this function fill out MPI send buffers
//it copies memory corresponding to a box of sites on the lattice 
//to send buffer memory area which matches the geometry of a dst box
// lattice : memory of gauge or spinor field
// gd_t : global, even or odd
// src, dst : source and destination geometry boxes. Must have same size.
// bytes_per_site : size in bytes of the local object to copy. Could be 4 gauge fields or a (half)spinor
// sendbuf : memory destination to be send over MPI to a matching recv buffer in the extended lattice (dst)
void syncBoxToBuffer(enum gd_type gd_t, int bytes_per_site, box_t *src, box_t *dst, void *lattice, void *sendbuf) {
    const int srcparity = boxParity(src);
    const int dstparity = boxParity(dst);
    const int isParityDifferent = (srcparity==dstparity)?0:1;
    const int vol = boxVolume(dst); //src is assumed to be equal
    const int dst_even_vol = vol/2 + (vol&(1^dstparity)); //if the number of sites in volume is odd and parity is even, then add 1 to vol/2
    const int base_dix=dst->base_index;

#ifndef NDEBUG
    //Check that src and dst have the same geometry!!
    error(((src->h[0]-src->l[0]) != (dst->h[0]-dst->l[0])), 1, __func__ , "Send/Recv buffer dimension 0 do not match!");
    error(((src->h[1]-src->l[1]) != (dst->h[1]-dst->l[1])), 1, __func__ , "Send/Recv buffer dimension 1 do not match!");
    error(((src->h[2]-src->l[2]) != (dst->h[2]-dst->l[2])), 1, __func__ , "Send/Recv buffer dimension 2 do not match!");
    error(((src->h[3]-src->l[3]) != (dst->h[3]-dst->l[3])), 1, __func__ , "Send/Recv buffer dimension 3 do not match!");

    const int src_even_vol = vol/2 + (vol&(1^srcparity)); //if the number of sites in volume is odd and parity is even, then add 1 to vol/2
    lprintf("SYNC",200,"VOLUME/PARITY: SRC=%d[even=%d]/%d  DST=%d[even=%d]/%d\n", vol, src_even_vol, srcparity, vol, dst_even_vol, dstparity);
#endif 

    for (int x0=0;x0<(src->h[0]-src->l[0]);x0++){
    for (int x1=0;x1<(src->h[1]-src->l[1]);x1++){
    for (int x2=0;x2<(src->h[2]-src->l[2]);x2++){
    for (int x3=0;x3<(src->h[3]-src->l[3]);x3++){
        //compute parity of site to skip even or odd sites when needed
        int parity = (x0+x1+x2+x3+srcparity)%2; //0->even ; 1->odd
        if (gd_t==EVEN && parity==1) continue;
        if (gd_t==ODD && parity==0) continue;
        int six=ipt_ext(x0+src->l[0],x1+src->l[1],x2+src->l[2],x3+src->l[3]);
        int dix=ipt_ext(x0+dst->l[0],x1+dst->l[1],x2+dst->l[2],x3+dst->l[3])-base_dix;
        // if srcparity!=dstparity, the destination parity needs to be changed
        // i.e. the dstparity must always match the srcparity
        // this is because on the receiving neighbor processor the parity
        // will be the one of srcparity
        // we here force that six and dix have same parity
        // if srcparity != dstparity we must correct it
#if 0
        if (srcparity!=dstparity) {
            //NB: boxes can have an odd number of points e.g. 
            // srcbox | vol/2+1 | vol/2 | srcparity even
            // dstbox | vol/2 | vol/2+1 |
            // srcbox | vol/2 | vol/2+1 | srcparity odd
            // dstbox | vol/2+1 | vol/2 |
            if(parity==0) dix -= dst_even_vol;//subtract even volume of srcbox 
            if(parity==1) dix += src_even_vol;//add even volume of srcbox. 
            //Note: since the volumes are equal and parity different: src_even_vol=vol-dst_even_vol
        }
#else
        //this avoid a branch otherwise it's equivalent to the code above
        dix += isParityDifferent*((parity*vol)-dst_even_vol);
#endif
        // if ( !(parity) && !(dix<vol/2)) dix -= vol/2; //must be even
        // if (  (parity) &&  (dix<vol/2)) dix += vol/2; //must be odd
        char *srcbuf = ((char*)lattice)+six*bytes_per_site;
        char *dstbuf = ((char*)sendbuf)+dix*bytes_per_site;
#ifndef NDEBUG
        lprintf("SYNC",200,"Writing [%d] bytes: %d -> %d [%d,%d,%d,%d] sizeof(suNg)=%d\n", bytes_per_site, six, dix, x0, x1, x2, x3,sizeof(suNg));
        if (dix<0 || !(dix<vol)) lprintf("SYNC",1,"ERROR in SYNC: dix=%d vol=%d\n", dix, vol);
        if (six<0 || !(six<VOLUME)) lprintf("SYNC",1,"ERROR in SYNC: six=%d\n", six);
        if (parity!=((dix<src_even_vol)?0:1)) lprintf("SYNC",1,"ERROR in SYNC: parity mismatch: six=%d[parity=%d] dix=%d[parity=%d]\n", six, parity, dix, ((dix<src_even_vol)?0:1));
#endif

        memcpy(dstbuf, srcbuf, bytes_per_site);

    }}}}
}

void sync_field(geometry_descriptor *gd, int bytes_per_site, int is_spinor_like, void *latticebuf) {
    enum gd_type gd_t = GLOBAL;
    //TODO: the type should really be in the descriptor itself
    // we shouldn't compare pointers...
    if (gd == &glat_even) gd_t = EVEN;
    if (gd == &glat_odd) gd_t = ODD;

    void *const sendbuf_base=latticebuf; //TODO: we could move the sendbuffers outside the fields. for now the memory is together with lattice data.
    const int n_buffers = is_spinor_like ? gd->nbuffers_spinor : gd->nbuffers_gauge;
#ifndef NDEBUG
    lprintf(LOGTAG,200,"Geometry size: %d \n", is_spinor_like? gd->gsize_spinor : gd->gsize_gauge );
    lprintf(LOGTAG,200,"SPINOR_LIKE: %d  NBUFFERS: %d\n", is_spinor_like, n_buffers );
#endif
    int i=0;
    box_t *L=geometryBoxes.next; //first border buffer. L3 are first, then L2
    while(L && i<n_buffers) {
#ifndef NDEBUG
        lprintf(LOGTAG,200,"Copying TO RECV box:\n");
        printBox(L,200);
        lprintf(LOGTAG,200,"Copying to index %d -> %d\n", gd->sbuf_start[i], gd->sbuf_start[i]+boxVolume(L));
#endif
        syncBoxToBuffer(gd_t, bytes_per_site, L->sendBox, L, latticebuf, ((char*)sendbuf_base)+bytes_per_site*gd->sbuf_start[i]);
        L=L->next; i++;
    }

}

// this function fills out a number of global variables:
// glattice, glat_even, glat_odd: geometry descriptors
// ipt, iup, idn: look up tables to move around the lattice
// imask: bitfield to encode some properties of the point 
void define_geometry() {
    lprintf(LOGTAG,1,"Define Lattice Geometry...\n");
    geometry_index_init();
    enumerate_lattice();
    atexit(&geometry_index_free);
    lprintf(LOGTAG,1,"Define Lattice Geometry... Done.\n");
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//             TESTING FUNCTIONS
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

void resetArray(char *array, int len) {
    for(int i=0;i<len;i++) array[i]=0;
}


static char *LOGTESTTAG="GEOMETRY TEST";

int checkArray(char *array, int len, int local_len) {
    int errors=0;
    for(int i=0;i<len;i++) {
        if (i<local_len && array[i]!=1) {
            errors++;
            lprintf(LOGTESTTAG,1,"ERROR: Local site %d used %d!=1 times\n", i, array[i]);
        }
        if (!(i<local_len) && array[i]!=0) {
            errors++;
            lprintf(LOGTESTTAG,1,"ERROR: Wrong local site used %d!=1 times\n", i, array[i]);
        }
    }
    return errors;
}

//check if a point is inside the local lattice
int isLocal(int x0_ext, int x1_ext, int x2_ext, int x3_ext) {
    x0_ext=safe_mod(x0_ext,T_EXT);
    x1_ext=safe_mod(x1_ext,X_EXT);
    x2_ext=safe_mod(x2_ext,Y_EXT);
    x3_ext=safe_mod(x3_ext,Z_EXT);
    // | 0 .. | T_BORDER ... T+T_BORDER-1 | T+T_BORDER .. T_EXT-1 |
    return (x0_ext>(T_BORDER-1)) && (x0_ext<(T+T_BORDER)) &&
           (x1_ext>(X_BORDER-1)) && (x1_ext<(X+X_BORDER)) &&
           (x2_ext>(Y_BORDER-1)) && (x2_ext<(Y+Y_BORDER)) &&
           (x3_ext>(Z_BORDER-1)) && (x3_ext<(Z+Z_BORDER));
}

//we test that moving int he direction active in the mask
//we end up inside the LOCAL lattice (no halos allowed)
//NB: the starting point can be in the halo
int checkMask(char mask, int x0_ext, int x1_ext, int x2_ext, int x3_ext ) {
    int errors=0;
    if ((mask & T_UP_MASK) && !isLocal(x0_ext+1, x1_ext, x2_ext, x3_ext)) {
        errors++;
        lprintf(LOGTESTTAG,1,"ERROR: Mask T_UP at site ext[%d,%d,%d,%d]\n", x0_ext, x1_ext, x2_ext, x3_ext);

    }
    if ((mask & T_DN_MASK) && !isLocal(x0_ext-1, x1_ext, x2_ext, x3_ext)) {
        errors++;
        lprintf(LOGTESTTAG,1,"ERROR: Mask T_DN at site ext[%d,%d,%d,%d]\n", x0_ext, x1_ext, x2_ext, x3_ext);

    }
    if ((mask & X_UP_MASK) && !isLocal(x0_ext, x1_ext+1, x2_ext, x3_ext)) {
        errors++;
        lprintf(LOGTESTTAG,1,"ERROR: Mask X_UP at site ext[%d,%d,%d,%d]\n", x0_ext, x1_ext, x2_ext, x3_ext);

    }
    if ((mask & X_DN_MASK) && !isLocal(x0_ext, x1_ext-1, x2_ext, x3_ext)) {
        errors++;
        lprintf(LOGTESTTAG,1,"ERROR: Mask X_DN at site ext[%d,%d,%d,%d]\n", x0_ext, x1_ext, x2_ext, x3_ext);

    }
    if ((mask & Y_UP_MASK) && !isLocal(x0_ext, x1_ext, x2_ext+1, x3_ext)) {
        errors++;
        lprintf(LOGTESTTAG,1,"ERROR: Mask Y_UP at site ext[%d,%d,%d,%d]\n", x0_ext, x1_ext, x2_ext, x3_ext);

    }
    if ((mask & Y_DN_MASK) && !isLocal(x0_ext, x1_ext, x2_ext-1, x3_ext)) {
        errors++;
        lprintf(LOGTESTTAG,1,"ERROR: Mask Y_DN at site ext[%d,%d,%d,%d]\n", x0_ext, x1_ext, x2_ext, x3_ext);

    }
    if ((mask & Z_UP_MASK) && !isLocal(x0_ext, x1_ext, x2_ext, x3_ext+1)) {
        errors++;
        lprintf(LOGTESTTAG,1,"ERROR: Mask Z_UP at site ext[%d,%d,%d,%d]\n", x0_ext, x1_ext, x2_ext, x3_ext);

    }
    if ((mask & Z_DN_MASK) && !isLocal(x0_ext, x1_ext, x2_ext, x3_ext-1)) {
        errors++;
        lprintf(LOGTESTTAG,1,"ERROR: Mask Z_DN at site ext[%d,%d,%d,%d]\n", x0_ext, x1_ext, x2_ext, x3_ext);

    }

    return errors;

}

int test_define_geometry() {

    char *idxcheck;
    const int ext_vol=T_EXT*X_EXT*Y_EXT*Z_EXT;
    idxcheck=malloc(ext_vol*sizeof(char));

    int total_errors=0;
    
    //check that ipt returns values within range
    //check that values for even and odd sites are correct
    int boxnumber=0;
    box_t *b = &geometryBoxes;
    while(b) {
        lprintf(LOGTESTTAG,1,"Testing Box #%d\n",boxnumber++);
        printBox(b,1);

        const int sign = boxParity(b);
        const int base_idx = b->base_index;
        const int vol = boxVolume(b);
        const int even_vol = vol/2 + (vol&(1^sign));
        lprintf(LOGTESTTAG,1,"Box: sign=%d vol=%d even_vol=%d\n",sign, vol, even_vol);

        //reset array to check local idx in geometry boxes
        resetArray(idxcheck,ext_vol);
        int errors=0;

        for (int x0_ext=b->l[0];x0_ext<b->h[0];x0_ext++){
        for (int x1_ext=b->l[1];x1_ext<b->h[1];x1_ext++){
        for (int x2_ext=b->l[2];x2_ext<b->h[2];x2_ext++){
        for (int x3_ext=b->l[3];x3_ext<b->h[3];x3_ext++){
            const int site_parity=((x0_ext-b->l[0])+(x1_ext-b->l[1])+(x2_ext-b->l[2])+(x3_ext-b->l[3])+sign)%2;
            const int ix = ipt_ext(x0_ext,x1_ext,x2_ext,x3_ext);
            const char mask = imask(ix);
            const int idx = ix-base_idx; //local idx 
            int site_errors=0;

            //add ix to idxcheck
            idxcheck[idx]++;

            //check that ix has correct range
            if (ix<0 || !(ix<ext_vol)) { 
                site_errors++; 
                lprintf(LOGTESTTAG,1,"Global index out of range: ix=%d ext_vol=%d\n", ix, ext_vol);
            }

            //check that idx has correct range
            if (idx<0 || !(idx<vol)) { 
                site_errors++; 
                lprintf(LOGTESTTAG,1,"Local index out of range: idx=%d vol=%d\n", idx, vol);
            }

            //check that idx has correct parity
            if (site_parity==0 && !(idx<even_vol)) {
                site_errors++; 
                lprintf(LOGTESTTAG,1,"Local index has incorrect parity: sign=%d idx=%d vol=%d even_vol=%d\n", site_parity, idx, vol, even_vol);
            }
            if (site_parity==1 &&  (idx<even_vol)) {
                site_errors++; 
                lprintf(LOGTESTTAG,1,"Local index has incorrect parity: sign=%d idx=%d vol=%d even_vol=%d\n", site_parity, idx, vol, even_vol);
            }

            //test mask
            site_errors += checkMask(mask, x0_ext, x1_ext, x2_ext, x3_ext);

            if (site_errors) {lprintf(LOGTESTTAG,1,"ERROR: %d errors  @ (%d, %d, %d, %d)\n", site_errors, x0_ext, x1_ext, x2_ext, x3_ext);}

            if (idx==4) {lprintf(LOGTESTTAG,1,"TEST: idx=%d (%d, %d, %d, %d)\n", idx, x0_ext, x1_ext, x2_ext, x3_ext);}


            errors += site_errors;
        }}}}

        errors+=checkArray(idxcheck,ext_vol,vol);
        lprintf(LOGTESTTAG,1,"Box Errors: %d errors\n", errors);

        total_errors += errors;

        b=b->next;
    }
    

    free(idxcheck);

    lprintf(LOGTESTTAG,1,"Geometry test finished with  %d errors\n", total_errors);

    return total_errors;
}
