#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include "geometry.h"
#include "global.h"
#include "error.h"
#include "logger.h"
#include "utils.h"
#include "memory.h"
#include "new_geometry.h"

#ifdef __cplusplus
    extern "C" {
#endif

inline int safe_mod(int x, int y)
{
  if (x >= 0)
    return (x % y);
  else
    return ((y - (abs(x) % y)) % y);
}

// This is to keep a list of boxes to fill out the field buffers
static box_t *geometryBoxes;
// TODO: make this into an array instead of a list?

// this will include the 4th L2 corner in the geometry
// it's not needed, was here for testng purposed
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

// ██ ███    ██ ██████  ███████ ██   ██ ██ ███    ██  ██████  
// ██ ████   ██ ██   ██ ██       ██ ██  ██ ████   ██ ██       
// ██ ██ ██  ██ ██   ██ █████     ███   ██ ██ ██  ██ ██   ███ 
// ██ ██  ██ ██ ██   ██ ██       ██ ██  ██ ██  ██ ██ ██    ██ 
// ██ ██   ████ ██████  ███████ ██   ██ ██ ██   ████  ██████  

///////////////////////////////////////////////////////////////////
// Basic functions to index points in a 4D box
///////////////////////////////////////////////////////////////////

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
// return an index from 0 to (X0*X1*X2*X3)/2-1
static int index_blocked(
    int b0, int b1, int b2, int b3, 
    int X0, int X1, int X2, int X3, 
    int x0, int x1, int x2, int x3
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
    idx = idx / 2; 

    return idx;
}


// ██████   ██████  ██   ██ 
// ██   ██ ██    ██  ██ ██  
// ██████  ██    ██   ███   
// ██   ██ ██    ██  ██ ██  
// ██████   ██████  ██   ██ 

///////////////////////////////////////////////////////////////////
//  Functions for handling BOXes
///////////////////////////////////////////////////////////////////

char invertMask(char mask) {
    return mask ^ FULL_MASK;
}

static const char*bt_names[]={"L0","L1","L2","L3","INNER","SENDBUF"};

static coord4 *rb_icoord=NULL; //global lookup table for recv buffers (L3+L2 borders)
static coord4 *sb_icoord=NULL; //global lookup table for send buffers

box_t* duplicateBox(box_t *B) {
    box_t *b=(box_t*)malloc(sizeof(box_t));
    error((b == NULL), 1, __func__ , "Cannot allocate memory for box");
    *b=*B; //make copy
    return b;
}

void insertBox(box_t *B, box_t *A){
    A->next = B->next;
    B->next = A;
}

//Free memory for boxes in the list
void freeBoxList(box_t *B) {
    while(B) {
        box_t *b=B->next;
        if (B->sendBox) { free(B->sendBox); }
        free(B); 
        B=b;
    }
}

static void freeGeometryBoxes(){
    freeBoxList(geometryBoxes);
    geometryBoxes=NULL;
}

int boxVolume(box_t *B){
    return (B->h[0]-B->l[0])*(B->h[1]-B->l[1])*(B->h[2]-B->l[2])*(B->h[3]-B->l[3]);

}
int boxEvenVolume(box_t *B){
    const int box_vol=boxVolume(B);
    int even_vol=box_vol/2 + (box_vol&(1^(B->parity)));
    return even_vol;
}
int boxOddVolume(box_t *B){
    const int odd_vol=boxVolume(B)-boxEvenVolume(B);
    return odd_vol;
}

//calculate if the base point of the box is even (0) or odd (1)
//the box coordinates refers to the extended lattice
int boxParity(box_t *B) {
    return B->parity;
}
void setBoxParity(box_t *B) {
    int p = (B->l[0]+B->l[1]+B->l[2]+B->l[3]+PSIGN+T_BORDER+X_BORDER+Y_BORDER+Z_BORDER)%2; //remember global parity!
    B->parity=p;
}

void printBox(box_t *B, int level) {
    static int i=1;
    lprintf(LOGTAG, level, "%s %c [%2d,%2d,%2d,%2d]..[%2d,%2d,%2d,%2d] = [%2d,%2d,%2d,%2d] Idx=(%d+%d / %d+%d)",
        bt_names[B->type], 
        boxParity(B) ? '-' : '+',
        B->l[0],B->l[1],B->l[2],B->l[3],
        B->h[0],B->h[1],B->h[2],B->h[3],
        B->h[0]-B->l[0],B->h[1]-B->l[1],B->h[2]-B->l[2],B->h[3]-B->l[3],
        B->base_index,
        boxEvenVolume(B),
        B->base_index_odd,
        boxOddVolume(B)
    );
    if (B->sendBox) {
        box_t *S=B->sendBox;
        lprintf(LOGTAG, level, "\t");
        i--; //disable print \n
        printBox(S, level);
        i++; //re-enable print \n
    }
    if(i) lprintf(LOGTAG, level, "\n");
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

//Compute the base index and base_index_odd of each box in the list
//The order is: EVEN INNER, EVEN L3, ODD INNER, ODD L3, EVEN L2, ODD L2 
//assume the first box in the list is the INNER volume
void setBoxListBaseIndex(box_t *B) {
    error((B->type != INNER), 1, __func__ , "Box list deos not start with INNER box");
    int i=B->base_index;
    int sendbox_i=0; //base index of send buffers
    box_t *saveBox=B; //save start of list
    // count INNER and L3 EVEN part
    while(B && (B->type != L2)) {
        B->base_index = i;
        const int vol = boxEvenVolume(B);
        i += vol;
        //set sendbox
        if(B->type==L3) {
            if (B->sendBox) B->sendBox->base_index=sendbox_i;
            sendbox_i += boxEvenVolume(B->sendBox);
        }
        B=B->next;
    }
    //count INNER and L3 ODD part
    B=saveBox; //restart from the beginning of the list
    while(B && (B->type != L2)) {
        B->base_index_odd = i;
        const int vol = boxOddVolume(B);
        i += vol;
        //set sendbox
        if(B->type==L3) {
            if (B->sendBox) B->sendBox->base_index_odd=sendbox_i;
            sendbox_i += boxOddVolume(B->sendBox);
        }
        B=B->next;
    }
    saveBox=B; //save start of L2 list 
    //count L2 EVEN part
    B=saveBox; //restart from the beginning of the list
    while(B) {
        B->base_index = i;
        const int vol = boxEvenVolume(B);
        i += vol;
        if (B->sendBox) B->sendBox->base_index=sendbox_i;
        sendbox_i += boxEvenVolume(B->sendBox);
        B=B->next;
    }
    //count L2 ODD part
    B=saveBox; //restart from the beginning of the list
    while(B) {
        B->base_index_odd = i;
        const int vol = boxOddVolume(B);
        i += vol;
        if (B->sendBox) B->sendBox->base_index_odd=sendbox_i;
        sendbox_i += boxOddVolume(B->sendBox);
        B=B->next;
    }

}

//check if a point at index idx is inside the box
//returns 1 (true) or 0 (false)
int isInsideBox(box_t *B, int ix) {
    return (ix>=B->base_index && ix<(B->base_index+boxEvenVolume(B))) || 
            (ix>=B->base_index_odd && ix<(B->base_index_odd+boxOddVolume(B)));
}

int safe_ipt_ext(int x0_ext, int x1_ext, int x2_ext, int x3_ext) {
    const int ix = ipt_ext(safe_mod(x0_ext,T_EXT),safe_mod(x1_ext,X_EXT),safe_mod(x2_ext,Y_EXT),safe_mod(x3_ext,Z_EXT)); 
#ifndef NDEBUG
    if (ix<0 || ix>(T_EXT*X_EXT*Y_EXT*Z_EXT-1)) { 
        lprintf(LOGTAG,1,"Error3 ix=%d site [%d, %d, %d, %d]\n", ix,safe_mod(x0_ext,T_EXT),safe_mod(x1_ext,X_EXT),safe_mod(x2_ext,Y_EXT),safe_mod(x3_ext,Z_EXT));
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

static box_t *makeLocalBox(){
    box_t B = { .l={T_BORDER,X_BORDER,Y_BORDER,Z_BORDER}, 
                .h={T+T_BORDER,X+X_BORDER,Y+Y_BORDER,Z+Z_BORDER}, 
                .base_index=0, 
                .base_index_odd=0, 
                .mask=0, 
                .type=INNER,
                .sendBox=NULL,
                .ipt_ext=NULL,
                .icoord=NULL,
                .next=NULL
              };
    setBoxParity(&B);
    return duplicateBox(&B);
}

//this function create a new box starting from a box B
//and a mask of directions
//the box returned corrsponds to a border of the box in the 
//directions indicated by the mask
static box_t *makeBorderBox(char mask) {
        //check that this is actually a border box otherwise return NULL
        if (mask==0) return NULL;

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
        setBoxParity(b);

        //make sendBox
        box_t *s=makeLocalBox(); 
        s->mask=0; s->type = SENDBUF;
        if (mask & T_UP_MASK) { s->h[0] = 2*T_BORDER; } else     
        if (mask & T_DN_MASK) { s->l[0] = T; }
        if (mask & X_UP_MASK) { s->h[1] = 2*X_BORDER; } else
        if (mask & X_DN_MASK) { s->l[1] = X; }     
        if (mask & Y_UP_MASK) { s->h[2] = 2*Y_BORDER; } else
        if (mask & Y_DN_MASK) { s->l[2] = Y; }     
        if (mask & Z_UP_MASK) { s->h[3] = 2*Z_BORDER; } else
        if (mask & Z_DN_MASK) { s->l[3] = Z; } 
        setBoxParity(s);

        //link sendBox in b
        b->sendBox=s;

        return b;
}

// ███████ ███    ██ ██    ██ ███    ███ ███████ ██████   █████  ████████ ███████ 
// ██      ████   ██ ██    ██ ████  ████ ██      ██   ██ ██   ██    ██    ██      
// █████   ██ ██  ██ ██    ██ ██ ████ ██ █████   ██████  ███████    ██    █████   
// ██      ██  ██ ██ ██    ██ ██  ██  ██ ██      ██   ██ ██   ██    ██    ██      
// ███████ ██   ████  ██████  ██      ██ ███████ ██   ██ ██   ██    ██    ███████ 
////////////////////////////////////////////////////////////
// Functions to assign to each 4D lattice point an id
// This is the core part of the geometry definition
// Defines a global lookup table called ipt_ext
// based on ipt_ext, other look up tables are defined:
// iup, idn, imask
////////////////////////////////////////////////////////////

// Lattice memory layout is as follows:
// |<- 0 index EVEN   
// | EVEN LOCAL VOL | EVEN L3 #1 | ... | EVEN L3 #n |
// | ODD  LOCAL VOL | ODD  L3 #1 | ... | ODD  L3 #n |
// and only for the gauge (glattice)
// | EVEN L2 #1 | ... | EVEN L2 #m |
// | ODD  L2 #1 | ... | ODD  L2 #m |
//
// Send buffers are outside the fields
// they are organized similarly:
// | SBUF EVEN L3 #1 | ... | SBUF EVEN L3 #n | 
// | SBUF ODD  L3 #1 | ... | SBUF ODD  L3 #n | 
// | SBUF EVEN L2 #1 | ... | SBUF EVEN L2 #m |
// | SBUF ODD  L2 #1 | ... | SBUF ODD  L2 #m |
//
// Send buffers are shared among all fields of same size
// they are always allocated as glattice 
// i.e. they contain all L3 and L2 buffers
// TODO: change this: make sendbuf_alloc take the full size of the sendbuffer, not the size per site

//allocate memory for global indexes
// ipt, iup, idn, imask
static void index_alloc() {
    const int ext_volume = T_EXT*X_EXT*Y_EXT*Z_EXT;
    const int N = (1+2*4)*ext_volume; // 1 for ipt + 4 dir x 2 for iup and idn
    ipt = (int*)malloc(N * sizeof(*ipt));
    error((ipt == NULL), 1, __func__ , "Cannot allocate memory for ipt, iup, idn");
    iup = ipt + ext_volume;
    idn = iup + (4*ext_volume);

    imask = (char*)malloc(ext_volume * sizeof(*imask));
    error((imask == NULL), 1, __func__ , "Cannot allocate memory for imask");

    //allocate memory for send/recv buffers icoord
    int sb_volume=0;
    box_t *L=geometryBoxes->next; //first border buffer. L3 are first, then L2
    while(L) { sb_volume += boxVolume(L); L=L->next; }
    const int rb_volume = sb_volume+boxVolume(geometryBoxes); //TODO: can we not count the INNER volume?
    rb_icoord=(coord4*)malloc((rb_volume+sb_volume)*sizeof(coord4));
    error((rb_icoord == NULL), 1, __func__ , "Cannot allocate memory for send/recv icoord");
    sb_icoord=rb_icoord+rb_volume;
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
    if (rb_icoord) {
        free(rb_icoord);
        rb_icoord=NULL;
        sb_icoord=NULL;
    }
}

// fill out ipt, iup, idn, imask
static void enumerate_lattice() {

    // const int size[4] = { T, X, Y, Z };
    // const int borders[4] = { T_BORDER, X_BORDER, Y_BORDER, Z_BORDER };
    // const int ext_size[4] = { T_EXT, X_EXT, Y_EXT, Z_EXT };
    // const int ext_vol = T_EXT*X_EXT*Y_EXT*Z_EXT;
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
    geometryBoxes = makeLocalBox();
    atexit(&freeGeometryBoxes); //free list of geometry boxes at exit

    box_t *L=geometryBoxes; //this is our list of boxes. first box is the local lattice B
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

    setBoxListBaseIndex(geometryBoxes);
    printBoxList(geometryBoxes, 10);

    //allocate memory for global indexes
    //requires boxes to be defined for sendbuffer icoord
    index_alloc();
    atexit(&index_free);

    const int INVALID_POINT = -1;
    //now the list of boxes is not going to cover the whole
    //extended lattice (e.g. L1 or L0 blocks or even some L2 might not there)
    //so we prefill the ipt, iup, idn array with INVALID_POINT
    // we prefill imask = 0 
#if 0
    //plain implementation
    for (int ix=0, x0_ext=0;x0_ext<T_EXT;x0_ext++){
    for (int x1_ext=0;x1_ext<X_EXT;x1_ext++){
    for (int x2_ext=0;x2_ext<Y_EXT;x2_ext++){
    for (int x3_ext=0;x3_ext<Z_EXT;x3_ext++){
        //prefill tables with INVALID_POINT
        ipt_ext(x0_ext,x1_ext,x2_ext,x3_ext) = INVALID_POINT;
#ifndef NDEBUG
        const int ext_vol = T_EXT*X_EXT*Y_EXT*Z_EXT;
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
#else
    //use memset instead
    //NB: memset only works with 0 and -1 (assuming 2's complement numbers)
    memset(ipt, -1, sizeof(int)*T_EXT*X_EXT*Y_EXT*Z_EXT*5); //assume ipt, iup and idn are contiguous in memory
    memset(imask, 0, T_EXT*X_EXT*Y_EXT*Z_EXT); 
#endif

    //enumerate points in each box
    //fill out ipt_ext
    box_t *b = geometryBoxes;
    while(b) {
        const int sign = boxParity(b);
        const int sb_sign = b->sendBox ? boxParity(b->sendBox) : -1;
        for (int x0=0;x0<(b->h[0]-b->l[0]);x0++){
        for (int x1=0;x1<(b->h[1]-b->l[1]);x1++){
        for (int x2=0;x2<(b->h[2]-b->l[2]);x2++){
        for (int x3=0;x3<(b->h[3]-b->l[3]);x3++){
            // const int ix = index_blocked(BLK_T,BLK_X,BLK_Y,BLK_Z,
            //                              b->h[0]-b->l[0],b->h[1]-b->l[1],b->h[2]-b->l[2],b->h[3]-b->l[3],
            //                              x0_ext-b->l[0],x1_ext-b->l[1],x2_ext-b->l[2],x3_ext-b->l[3], 
            //                              sign,
            //                              b->base_index,
            //                              b->base_index_odd
            //                             );

            const int idx = index_blocked(BLK_T,BLK_X,BLK_Y,BLK_Z,
                                         b->h[0]-b->l[0],b->h[1]-b->l[1],b->h[2]-b->l[2],b->h[3]-b->l[3],
                                         x0,x1,x2,x3
                                        );
            const int parity = (x0+x1+x2+x3+sign)%2;
            const int ix = idx + (parity ? b->base_index_odd : b->base_index );
#ifndef NDEBUG
            // if (idx<0 || !(idx<boxVolume(b))) lprintf(LOGTAG,1,"ERROR in ipt: Volume=%d idx=%d  @ (%d, %d, %d, %d)\n", boxVolume(b), idx, x0_ext, x1_ext, x2_ext, x3_ext);
            if (ix<0 || ix>(T_EXT*X_EXT*Y_EXT*Z_EXT-1)) lprintf(LOGTAG,1,"ERROR2 ix=%d @ (%d, %d, %d, %d)\n", ix, x0_ext, x1_ext, x2_ext, x3_ext);
#endif
            ipt_ext(x0+b->l[0],x1+b->l[1],x2+b->l[2],x3+b->l[3]) = ix;

            //fill out rb_icoord
            const coord4 rc = { .x={x0+b->l[0],x1+b->l[1],x2+b->l[2],x3+b->l[3]} } ;
            rb_icoord[ix] = rc;
            b->icoord = rb_icoord;           

            if(b->sendBox) {
                //fill out sb_icoord
                box_t *SB=b->sendBox;
                const int sb_parity = (x0+x1+x2+x3+sb_sign)%2;
                const int sb_ix = idx + (sb_parity ? SB->base_index_odd : SB->base_index );
                const coord4 sc = { .x={x0+SB->l[0],x1+SB->l[1],x2+SB->l[2],x3+SB->l[3]} } ;
                sb_icoord[sb_ix] = sc;           
                SB->icoord = sb_icoord;           
            }
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
        if (ix<0 || ix>(T_EXT*X_EXT*Y_EXT*Z_EXT-1)) {
            lprintf("GEOMETRY",1,"Error ix=%d site [%d, %d, %d, %d]\n", ix,x0_ext,x1_ext,x2_ext,x3_ext);
        }
#endif

        char xmask=0; int iy=0;
        iy=safe_ipt_ext(x0_ext+1,x1_ext,x2_ext,x3_ext); iup(ix,0)=iy; if(isInsideBox(geometryBoxes, iy)) xmask |= T_UP_MASK; 
        iy=safe_ipt_ext(x0_ext-1,x1_ext,x2_ext,x3_ext); idn(ix,0)=iy; if(isInsideBox(geometryBoxes, iy)) xmask |= T_DN_MASK; 
        iy=safe_ipt_ext(x0_ext,x1_ext+1,x2_ext,x3_ext); iup(ix,1)=iy; if(isInsideBox(geometryBoxes, iy)) xmask |= X_UP_MASK; 
        iy=safe_ipt_ext(x0_ext,x1_ext-1,x2_ext,x3_ext); idn(ix,1)=iy; if(isInsideBox(geometryBoxes, iy)) xmask |= X_DN_MASK; 
        iy=safe_ipt_ext(x0_ext,x1_ext,x2_ext+1,x3_ext); iup(ix,2)=iy; if(isInsideBox(geometryBoxes, iy)) xmask |= Y_UP_MASK; 
        iy=safe_ipt_ext(x0_ext,x1_ext,x2_ext-1,x3_ext); idn(ix,2)=iy; if(isInsideBox(geometryBoxes, iy)) xmask |= Y_DN_MASK; 
        iy=safe_ipt_ext(x0_ext,x1_ext,x2_ext,x3_ext+1); iup(ix,3)=iy; if(isInsideBox(geometryBoxes, iy)) xmask |= Z_UP_MASK; 
        iy=safe_ipt_ext(x0_ext,x1_ext,x2_ext,x3_ext-1); idn(ix,3)=iy; if(isInsideBox(geometryBoxes, iy)) xmask |= Z_DN_MASK; 

        imask(ix)=xmask;

    }}}}


}


//  ██████  ███████  ██████  ███    ███ ███████ ████████ ██████  ██    ██ 
// ██       ██      ██    ██ ████  ████ ██         ██    ██   ██  ██  ██  
// ██   ███ █████   ██    ██ ██ ████ ██ █████      ██    ██████    ████   
// ██    ██ ██      ██    ██ ██  ██  ██ ██         ██    ██   ██    ██    
//  ██████  ███████  ██████  ██      ██ ███████    ██    ██   ██    ██    

// ██████  ███████ ███████  ██████ ██████  ██ ██████  ████████  ██████  ██████  
// ██   ██ ██      ██      ██      ██   ██ ██ ██   ██    ██    ██    ██ ██   ██ 
// ██   ██ █████   ███████ ██      ██████  ██ ██████     ██    ██    ██ ██████  
// ██   ██ ██           ██ ██      ██   ██ ██ ██         ██    ██    ██ ██   ██ 
// ██████  ███████ ███████  ██████ ██   ██ ██ ██         ██     ██████  ██   ██ 

///////////////////////////////////////////////////////////////////
// Fills out global variables for lattice geometry description:
// glattice, glat_even and glat_odd
///////////////////////////////////////////////////////////////////

static void gd_free_mem(geometry_descriptor *gd) {
    if(gd->rbuf_len) {
        free(gd->rbuf_len);
        gd->rbuf_len       = NULL;
        gd->sbuf_len       = NULL;
        gd->sbuf_to_proc   = NULL;
        gd->sbuf_start     = NULL;
        gd->rbuf_from_proc = NULL;
        gd->rbuf_start     = NULL;
    }    
    if(gd->master_start) {
        free(gd->master_start);
        gd->master_start   = NULL;
        gd->master_end     = NULL;
    }
}

// if N>0 this allocates memory for the send/rec index buffers and 
// master pieces indices. The number of master pieces is NBuf+NInner
// otherwise sets them to NULL
static void gd_alloc_mem(geometry_descriptor *gd, int NBuf, int NInner) {
    int *buf = NULL;
    if (NBuf>0) {
        buf = (int*)malloc((6*NBuf) * sizeof(int));
        error((buf == NULL), 1, __func__ , "Cannot allocate memory");
    }

    gd->rbuf_len         = buf; if (buf) buf += NBuf;
    gd->sbuf_len         = buf; if (buf) buf += NBuf;
    gd->sbuf_to_proc     = buf; if (buf) buf += NBuf;
    gd->sbuf_start       = buf; if (buf) buf += NBuf;
    gd->rbuf_from_proc   = buf; if (buf) buf += NBuf;
    gd->rbuf_start       = buf;

    buf = NULL;
    if (NInner>0) {
        buf = (int*)malloc((2*(NInner+NBuf)) * sizeof(int));
        error((buf == NULL), 1, __func__ , "Cannot allocate memory");
    }
    gd->master_start     = buf; if (buf) buf += NBuf+NInner; //+1 because of the local lattice box
    gd->master_end       = buf; //same size as master_start

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
// this is a simple bitmask with GLOBAL = EVEN | ODD
enum gd_type {
    EVEN   = 1, 
    ODD    = 2,
    GLOBAL = 3
};

static void gd_free() {
    gd_free_mem(&glattice);
    gd_free_mem(&glat_even);
    gd_free_mem(&glat_odd);
}

static void gd_init() {
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

    const int volume = boxVolume(geometryBoxes);
    const int even_volume = boxEvenVolume(geometryBoxes); //VOLUME/2 + ((PSIGN==0)? (VOLUME&1) : 0);
    const int odd_volume = volume-even_volume;

    lprintf(LOGTAG,50,"L3 BORDERS=%d  L2 BORDERS=%d\n", L3_BORDER, L2_BORDER);

    //set up global lattice    
    gd_alloc_mem(&glattice, 2*(L3_BORDER+L2_BORDER), 2);
    gd_set_copy(&glattice);
    glattice.nbuffers_gauge = 2*(L3_BORDER+L2_BORDER);
    glattice.nbuffers_spinor = 2*L3_BORDER;
    glattice.inner_master_pieces = 2;
    glattice.local_master_pieces = 2;
    glattice.master_start[0]=geometryBoxes->base_index;
    glattice.master_end[0]=glattice.master_start[0]+even_volume-1;    
    glattice.master_start[1]=geometryBoxes->base_index_odd;
    glattice.master_end[1]=glattice.master_start[1]+odd_volume-1;
    glattice.master_shift=0;
    glattice.total_gauge_master_pieces=2+2*(L3_BORDER+L2_BORDER);
    glattice.total_spinor_master_pieces=2+2*(L3_BORDER);
    glattice.gsize_gauge=volume;  //this will be increased in gd_set_boxIndices()
    glattice.gsize_spinor=volume;  //this will be increased in gd_set_boxIndices()

    //set up even lattice
    gd_alloc_mem(&glat_even, L3_BORDER, 1);
    gd_set_copy(&glat_even);
    glat_even.nbuffers_gauge = -1; //TODO: can we do better?
    glat_even.nbuffers_spinor = L3_BORDER;
    glat_even.inner_master_pieces = 1;
    glat_even.local_master_pieces = 1;
    glat_even.master_start[0]=geometryBoxes->base_index;
    glat_even.master_end[0]=glat_even.master_start[0]+even_volume-1;
    glat_even.master_shift=0;
    glat_even.total_gauge_master_pieces=-1; //TODO: can we do better?
    glat_even.total_spinor_master_pieces=1+L3_BORDER;
    glat_even.gsize_gauge=-1;  //TODO: can we do better?
    glat_even.gsize_spinor=even_volume;  //this will be increased in gd_set_boxIndices()
    
    //set up odd lattice
    gd_alloc_mem(&glat_odd, L3_BORDER, 1);
    gd_set_copy(&glat_odd);
    glat_odd.nbuffers_gauge = -1;
    glat_odd.nbuffers_spinor = L3_BORDER;
    glat_odd.inner_master_pieces = 1;
    glat_odd.local_master_pieces = 1;
    glat_odd.master_start[0]=geometryBoxes->base_index_odd;
    glat_odd.master_end[0]=glat_odd.master_start[0]+odd_volume-1;
    glat_odd.master_shift=geometryBoxes->base_index_odd;
    glat_odd.total_gauge_master_pieces=-1; //TODO: can we do better?
    glat_odd.total_spinor_master_pieces=1+L3_BORDER;
    glat_odd.gsize_gauge=-1;  //TODO: can we do better?
    glat_odd.gsize_spinor=odd_volume;  //this will be increased in gd_set_boxIndices()

    //free memory at exit
    atexit(&gd_free);

}

//sets: rbuf_len, sbuf_len, rbuf_from_proc, rbuf_start, sbuf_to_proc, sbuf_start,
// master_start, master_end
void gd_set_boxIndices() {
    int ib=0; //index of buffer in geometry descriptor
    box_t *L=geometryBoxes->next; //first border buffer. L3 are first, then L2
    while(L) {
        //each box contains an even and odd part
        const int vol = boxVolume(L);
        //const int sign = boxParity(L);
        const int even_vol = boxEvenVolume(L);
        const int odd_vol = vol - even_vol;
        const int send_even_vol = boxEvenVolume(L->sendBox);
        const int send_odd_vol = boxOddVolume(L->sendBox);
        const int recv_proc=nearProc(L->mask);
        const int send_proc=nearProc(invertMask(L->mask));
        lprintf(LOGTAG,50,"ID=%d Level=%d mask=" _PRINT_BYTE " Recv proc=%d Send proc=%d even_vol=%d odd_vol=%d send_even_vol=%d send_odd_vol=%d\n", ib, (int)(L->type),_BINARY(L->mask), recv_proc, send_proc, even_vol, odd_vol, send_even_vol, send_odd_vol);
        
        //There are two send/recv buffers per box: EVEN and ODD
        //EVEN part
        const int ib_even = 2*ib;
        glattice.rbuf_len[ib_even]=even_vol;
        glattice.sbuf_len[ib_even]=send_even_vol;
        glattice.rbuf_start[ib_even]=L->base_index;
        glattice.sbuf_start[ib_even]=L->sendBox->base_index;
        glattice.rbuf_from_proc[ib_even]=recv_proc;
        glattice.sbuf_to_proc[ib_even]=send_proc;
        glattice.master_start[ib_even+2]=L->base_index; //+2 because index 0 and 1 is the local master piece
        glattice.master_end[ib_even+2]=L->base_index+even_vol-1;
        //ODD part
        // const int ib_odd=ib+glattice.nbuffers_gauge/2; //nbuffers_gauge is always even
        const int ib_odd = 2*ib+1;
        glattice.rbuf_len[ib_odd]=odd_vol;
        glattice.sbuf_len[ib_odd]=send_odd_vol;
        glattice.rbuf_start[ib_odd]=L->base_index_odd;
        glattice.sbuf_start[ib_odd]=L->sendBox->base_index_odd;
        glattice.rbuf_from_proc[ib_odd]=recv_proc;
        glattice.sbuf_to_proc[ib_odd]=send_proc;
        glattice.master_start[ib_odd+2]=L->base_index_odd; //+2 because index 0 and 1 is the local master piece
        glattice.master_end[ib_odd+2]=L->base_index_odd+odd_vol-1;
        //add to size 
        glattice.gsize_gauge+=vol; //send buffers are not included
        if(L->type==L3) { glattice.gsize_spinor+=vol; }

        //EVEN / ODD geometries contain only L3 borders
        if (L->type==L3) {
            glat_even.rbuf_len[ib]=even_vol;
            glat_even.sbuf_len[ib]=send_even_vol;
            glat_even.rbuf_start[ib]=L->base_index;
            glat_even.sbuf_start[ib]=L->sendBox->base_index;
            glat_even.rbuf_from_proc[ib]=recv_proc;
            glat_even.sbuf_to_proc[ib]=send_proc;
            glat_even.master_start[ib+1]=L->base_index; //+1 because index 0 is the local master piece
            glat_even.master_end[ib+1]=L->base_index+even_vol-1;
            glat_even.gsize_spinor+=even_vol;

            glat_odd.rbuf_len[ib]=odd_vol;
            glat_odd.sbuf_len[ib]=send_odd_vol;
            glat_odd.rbuf_start[ib]=L->base_index_odd;
            glat_odd.sbuf_start[ib]=L->sendBox->base_index_odd;
            glat_odd.rbuf_from_proc[ib]=recv_proc;
            glat_odd.sbuf_to_proc[ib]=send_proc;
            glat_odd.master_start[ib+1]=L->base_index_odd; //+1 because index 0 is the local master piece
            glat_odd.master_end[ib+1]=L->base_index_odd+odd_vol-1;
            glat_odd.gsize_spinor+=odd_vol; 
        }

        L=L->next; ib++;
    }

}


// ███████ ███████ ███    ██ ██████      ██████  ██    ██ ███████ ███████ ███████ ██████  ███████ 
// ██      ██      ████   ██ ██   ██     ██   ██ ██    ██ ██      ██      ██      ██   ██ ██      
// ███████ █████   ██ ██  ██ ██   ██     ██████  ██    ██ █████   █████   █████   ██████  ███████ 
//      ██ ██      ██  ██ ██ ██   ██     ██   ██ ██    ██ ██      ██      ██      ██   ██      ██ 
// ███████ ███████ ██   ████ ██████      ██████   ██████  ██      ██      ███████ ██   ██ ███████ 

/////////////////////////////////////////////////////
// functions for sendbuffers
// send buffers are always allocated with a glattice geometry
// they are reused for all instances of a given field
/////////////////////////////////////////////////////

//we keep a list of send buffers created
typedef struct _SB_t {
    size_t bytes_per_site;
    void *buf;
} SB_t;

#define _MAX_SENDBUF 8
static SB_t send_buffers[_MAX_SENDBUF] = { 0 };
typedef int sb_handle;

#ifndef NDEBUG
static int isValidHandle(sb_handle h) {
    return (h>=0) && (h<_MAX_SENDBUF);
}
#endif

static void sendbuf_free() {
    for (int i=0; i<_MAX_SENDBUF; i++) {
        void *buf=send_buffers[i].buf;
        if (buf) afree(buf); //this was allocated with amalloc
        const SB_t empty_sb = { 0 };
        send_buffers[i] = empty_sb;
    }
}

static void sendbuf_init() {
    atexit(&sendbuf_free);
}

static SB_t* getSB(sb_handle h){
#ifndef NDEBUG
    if (!isValidHandle(h)) 
        error(1, 1, __func__ , "Attempt to access invalid sendbuf handle");
#endif
    return &send_buffers[h];
}

//find send buf of given size
//if not already existing returns first free slot
//if no free slots left, exit with error
static sb_handle findSB(size_t bytes_per_site){
    for (int i=0; i<_MAX_SENDBUF; i++) {
        if ((send_buffers[i].buf == NULL) || (send_buffers[i].bytes_per_site==bytes_per_site)) {
            //this is the first empty slot
            return i;
        }       
    }
    //not enough slots for sendbuffers
    error(1, 1, __func__ , "Not enough slots for sendbuffers: increase _MAX_SENDBUF.");
    return -1; //invalid handle 

}

//allocate a new sendbuffer on the glattice geometry
static sb_handle allocSB(size_t bytes_per_site) {
    //find slot for send buffer
    sb_handle sb_h = findSB(bytes_per_site);
    SB_t *new_sb=getSB(sb_h);
    //if this is an already filled slot, just return it
    if (new_sb->buf) return sb_h;

    //otherwise allocate memory
    //gsize_gauge is the local volume + rsend buffers. 
    //We subtract the local volume to obtain rsend buffer size = send buffer size
    const int nsites = glattice.gsize_gauge - boxVolume(geometryBoxes);
    lprintf(LOGTAG, 10, "Allocating SEND buf: %.1fkB [bytes per site:%zu sites=%d]\n", ((float)(bytes_per_site*nsites)/1024.), bytes_per_site, nsites);
    void *buf = amalloc(nsites*bytes_per_site, ALIGN); 
    error(buf == NULL, 1, __func__ , "Could not allocate memory space for sendbuffer"); 
    // fill our new structure
    new_sb->bytes_per_site = bytes_per_site;
    new_sb->buf=buf;

    return sb_h;
}

// function to allocate a sendbuffer
// returns memory buffer
// TODO: do we need EVEN/ODD here?
void *sendbuf_alloc(size_t bytes_per_site) {
    const sb_handle sb_h = allocSB(bytes_per_site);
    return getSB(sb_h)->buf;
}

//this prints out the state of the send buffer list
void sendbuf_report() {
    for (int i=0; i<_MAX_SENDBUF; i++) {
        SB_t * const sb=getSB(i);
        if (sb->buf==NULL) break;
        lprintf(LOGTAG, 1, "%d size=%zu\n", i, sb->bytes_per_site);
    }
}

#ifndef NDEBUG
static int isSiteParityCorrect(int idx, int parity, box_t *B);  //this is defined below in the testing section
#endif

//this function fill out MPI send buffers
//it copies memory corresponding to a box of sites on the lattice 
//to send buffer memory area which matches the geometry of a dst box
// lattice : memory of gauge or spinor field
// gd_t : global, even or odd
// src : source geometry box.
// bytes_per_site : size in bytes of the local object to copy. Could be 4 gauge fields or a (half)spinor
// sendbuf : memory destination to be send over MPI to a matching recv buffer in the extended lattice (dst)
static void syncBoxToBuffer(enum gd_type gd_t, int bytes_per_site, box_t *src, void *lattice, void *sendbuf) {
#ifndef NDEBUG
    lprintf("SYNC",1,"VOLUME/PARITY: SRC=%d[even=%d]/%d\n", boxVolume(src), boxEvenVolume(src), boxParity(src));
#endif 
    //EVEN part
    if (gd_t & EVEN) {
        const int vol = boxEvenVolume(src);
        // lprintf("SYNC",1," Copy EVEN\n");
        // printBox(src,1);
        for (int dix=src->base_index; dix<(src->base_index+vol);dix++){
            // lprintf("SYNC",1,"ptr=%p %p\n",src->icoord,sb_icoord);
            coord4 c = src->icoord[dix];
            // lprintf("SYNC",1,"c=[%d,%d,%d,%d]\n",c.x[0],c.x[1],c.x[2],c.x[3]);
            int six=ipt_ext(c.x[0],c.x[1],c.x[2],c.x[3]);
            char *srcbuf = ((char*)lattice)+six*bytes_per_site;
            char *dstbuf = ((char*)sendbuf)+dix*bytes_per_site;
            memcpy(dstbuf, srcbuf, bytes_per_site);
        }
    }
    //ODD part
    if (gd_t & ODD) {
        const int vol = boxOddVolume(src);
        // lprintf("SYNC",1," Copy ODD\n");
        for (int dix=src->base_index_odd; dix<(src->base_index_odd+vol);dix++){
            coord4 c = src->icoord[dix];
            int six=ipt_ext(c.x[0],c.x[1],c.x[2],c.x[3]);
            char *srcbuf = ((char*)lattice)+six*bytes_per_site;
            char *dstbuf = ((char*)sendbuf)+dix*bytes_per_site;
            memcpy(dstbuf, srcbuf, bytes_per_site);
        }
    }
}

void sync_field(geometry_descriptor *gd, int bytes_per_site, int is_spinor_like, void *latticebuf, void *sb_ptr) {
    enum gd_type gd_t = GLOBAL;
    //TODO: the type should really be in the descriptor itself
    // we shouldn't compare pointers...
    if (gd == &glat_even) gd_t = EVEN;
    if (gd == &glat_odd) gd_t = ODD;

    latticebuf = ((char*) latticebuf)-(bytes_per_site*gd->master_shift); //shift buffer by master_shift
    char *const sendbuf_base=(char*)sb_ptr;
    int n_buffers = is_spinor_like ? gd->nbuffers_spinor : gd->nbuffers_gauge;
    if (gd_t == GLOBAL) n_buffers /= 2; //we fill even and odd buffers at the same time
#ifndef NDEBUG
    lprintf(LOGTAG,50,"Geometry size: %d \n", is_spinor_like? gd->gsize_spinor : gd->gsize_gauge );
    lprintf(LOGTAG,50,"SPINOR_LIKE: %d  NBUFFERS: %d\n", is_spinor_like, n_buffers );
#endif
    int i=0;
    box_t *L=geometryBoxes->next; //first border buffer. L3 are first, then L2
    while(L && i<n_buffers) {
#ifndef NDEBUG
        lprintf(LOGTAG,50,"Copying to RECV box:\n");
        printBox(L,50);
        if (gd_t == GLOBAL) {
            lprintf(LOGTAG,50,"EVEN copying to memory index %d -> %d [len=%d]\n", gd->sbuf_start[2*i], gd->sbuf_start[2*i]+gd->sbuf_len[2*i], gd->sbuf_len[2*i]);
            lprintf(LOGTAG,50,"ODD copying to memory index %d -> %d [len=%d]\n", gd->sbuf_start[2*i+1], gd->sbuf_start[2*i+1]+gd->sbuf_len[2*i+1], gd->sbuf_len[2*i+1]);
        } else {
            lprintf(LOGTAG,50,"%s copying to memory index %d -> %d [len=%d]\n", (gd_t==EVEN)?"EVEN":"ODD",gd->sbuf_start[i], gd->sbuf_start[i]+gd->sbuf_len[i], gd->sbuf_len[i]);
        }
#endif
        // syncBoxToBuffer_old(gd_t, bytes_per_site, L->sendBox, L, latticebuf, sendbuf_base);
        syncBoxToBuffer(gd_t, bytes_per_site, L->sendBox, latticebuf, sendbuf_base);

        L=L->next; i++;
    }
}


#ifdef __cplusplus
    }
#endif

#ifdef WITH_GPU
#define _DECLARE_KERNEL(_name, _type) \
    __global__ void box_to_buffer_kernel_##_name(void *dst, _type *lattice, int base_index, int stride, coord4* icoord, int bytes_per_site) \
    { \
        _type* src_uncast; \
        int dix = threadIdx.x*BLOCK_SIZE + blockIdx.x + base_index; \
        coord4 c = icoord[dix]; \
        int six = ipt_ext_gpu(c.x[0], c.x[1], c.x[2], c.x[3]); \
        read_gpu_##_type(stride, (*src_uncast), lattice, six, 0);\
        char *srcbuf = (char*)src_uncast;\
        char *dstbuf = ((char*)dst) + dix*bytes_per_site; \
        /*Make this work for CUDA-aware MPI by adding a sendbuf that is allocated on the device! (SAM)*/ \
        cudaMemcpy(dstbuf, srcbuf, bytes_per_site, cudaMemcpyDeviceToHost); \
    }

#define _DECLARE_SYNC_BOX(_name, _type, _size) \
    static void sync_box_to_buffer_gpu_##_name(geometry_descriptor *gd, \
                                    box_t *src, \
                                    _type *lattice, \
                                    void *sendbuf) \
    { \
        const int bytes_per_site = (_size)*sizeof(*lattice); \
        _PIECE_FOR(gd, ixp) \
        { \
            int vol = 0; \
            int stride = 0; \
            if ((ixp%2)==0) vol = boxEvenVolume(src); \
            else vol = boxOddVolume(src); \
            coord4 *c = src->icoord; \
            coord4 *d_c; \
            cudaMalloc((void**)&d_c, sizeof(c)); \
            cudaMemcpy(d_c, c, sizeof(c), cudaMemcpyHostToDevice); \
            int grid = (vol -1)/BLOCK_SIZE +1; \
            stride = gd->master_start[ixp] - gd->master_end[ixp] + 1; \
            box_to_buffer_kernel_##_name<<<grid, BLOCK_SIZE>>>(sendbuf, lattice, src->base_index, stride, d_c, bytes_per_site); \
        } \
    }

#define _DECLARE_SYNC_FIELD(_name, _type, _geom) \
    void sync_field_gpu_##_name(geometry_descriptor *gd, \
                        _type *lattice,  \
                        void *sendbuf) \
    { \
        if (geometryBoxes==NULL) printf("geometryBoxes are not initialized.\n");\
        box_t *L = geometryBoxes->next; \
        int i = 0; \
        while (L && i < gd->nbuffers_##_geom) \
        { \
            sync_box_to_buffer_gpu_##_name(gd, L->sendBox, lattice, sendbuf); \
            L=L->next; i++; \
        } \
    }

#define _DECLARE_NEW_GEOM(_name, _type, _geom, _size) \
    _DECLARE_KERNEL(_name, _type) \
    _DECLARE_SYNC_BOX(_name, _type, _size) \
    _DECLARE_SYNC_FIELD(_name, _type, _geom)

_DECLARE_NEW_GEOM(gfield_f, suNf, gauge, 4);
_DECLARE_NEW_GEOM(spinor_field_f, suNf_spinor, spinor, 1);

#undef _DECLARE_NEW_GEOM
#undef _DECLARE_SYNC_BOX
#undef _DECLARE_SYNC_FIELD

#endif

#ifdef __cplusplus
    extern "C" {
#endif

// ███    ███      █████      ██     ███    ██ 
// ████  ████     ██   ██     ██     ████   ██ 
// ██ ████ ██     ███████     ██     ██ ██  ██ 
// ██  ██  ██     ██   ██     ██     ██  ██ ██ 
// ██      ██     ██   ██     ██     ██   ████ 

// this function fills out a number of global variables:
// ipt, iup, idn: look up tables to move around the lattice
// imask: bitfield to encode some properties of the point 
// glattice, glat_even, glat_odd: geometry descriptors
void define_geometry() {
    lprintf(LOGTAG,1,"Define Lattice Geometry...\n");
    enumerate_lattice();
    gd_init();
    gd_set_boxIndices();
    sendbuf_init();
    lprintf(LOGTAG,1,"Define Lattice Geometry... Done.\n"); 
}


// ████████ ███████ ███████ ████████ ██ ███    ██  ██████  
//    ██    ██      ██         ██    ██ ████   ██ ██       
//    ██    █████   ███████    ██    ██ ██ ██  ██ ██   ███ 
//    ██    ██           ██    ██    ██ ██  ██ ██ ██    ██ 
//    ██    ███████ ███████    ██    ██ ██   ████  ██████  

#ifndef NDEBUG

//checks if a point with given idx and parity has correct ranges
//in the box.
static int isSiteParityCorrect(int idx, int parity, box_t *B) {
    if (parity==0) {
        //even site
        return (idx>=B->base_index && idx<(B->base_index+boxEvenVolume(B)));
    }
    //odd site
    return (idx>=B->base_index_odd && idx<(B->base_index_odd+boxOddVolume(B)));
}
#endif

static void resetArray(char *array, int len) {
    for(int i=0;i<len;i++) array[i]=0;
}

static char *LOGTESTTAG="GEOMETRY TEST";

static int checkArray(char *array, int len, int local_len) {
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
static int isLocal(int x0_ext, int x1_ext, int x2_ext, int x3_ext) {
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
static int checkMask(char mask, int x0_ext, int x1_ext, int x2_ext, int x3_ext ) {
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

static int checkBoxNumbers(int boxnumber) {
        //number of parallel directions
    const int NPAR = (NP_T>1 ? 1 : 0) + (NP_X>1? 1 : 0) + (NP_Y>1 ? 1 : 0) + (NP_Z>1 ? 1 : 0);
    // Dimension 3 borders. These are for both gauge fields and spinors
    const int L3_BORDER = 2 * NPAR;
    // Dimension 2 borders. These are for the gauge fields
#ifdef _INCLUDE_UP_UP_L2
    const int L2_BORDER = 4*NPAR*(NPAR-1)/2; 
#else
    const int L2_BORDER = 3*NPAR*(NPAR-1)/2; 
#endif

    int correct_boxnumber=1+L3_BORDER+L2_BORDER;

    int errors = 0;

    if (boxnumber!=correct_boxnumber) {
        errors++;
        lprintf(LOGTESTTAG,1,"ERROR: Number of boxes incorrect. Expected %d but %d boxes found\n", correct_boxnumber, boxnumber);

    }

    return errors;


}

int test_define_geometry() {

    char *idxcheck;
    const int ext_vol=T_EXT*X_EXT*Y_EXT*Z_EXT;
    idxcheck=(char*)malloc(ext_vol*sizeof(char));

    int total_errors=0;
    
    //check that ipt returns values within range
    //check that values for even and odd sites are correct
    int boxnumber=0;
    box_t *b = geometryBoxes;
    while(b) {
        lprintf(LOGTESTTAG,1,"Testing Box #%d\n",++boxnumber);
        printBox(b,1);

        const int sign = boxParity(b);
        const int base_idx = b->base_index;
        const int base_idx_odd = b->base_index_odd;
        const int vol = boxVolume(b);
        const int even_vol = vol/2 + (vol&(1^sign));
        const int odd_vol = vol-even_vol;
        lprintf(LOGTESTTAG,1,"Box: sign=%d vol=%d even_vol=%d odd_vol=%d\n",sign, vol, even_vol, odd_vol);

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
            const int shift = (site_parity) ? base_idx_odd-even_vol : base_idx;
            const int idx = ix-shift; //local idx 
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

            if (site_errors) {lprintf(LOGTESTTAG,1,"ERROR: %d errors  ix=%d idx=%d shift=%d @ (%d, %d, %d, %d)\n", site_errors, ix, idx, shift, x0_ext, x1_ext, x2_ext, x3_ext);}

            errors += site_errors;
        }}}}

        errors+=checkArray(idxcheck,ext_vol,vol);
        lprintf(LOGTESTTAG,1,"Box Errors: %d errors\n", errors);

        total_errors += errors;

        b=b->next;
    }
    
    total_errors+=checkBoxNumbers(boxnumber);

    free(idxcheck);

    lprintf(LOGTESTTAG,1,"Geometry test finished with  %d errors\n", total_errors);

    return total_errors;
}

#ifdef __cplusplus
    }
#endif
