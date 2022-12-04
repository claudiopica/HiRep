/**
 * @file new_geometry.h
 * @brief Functions needed for the new geometry implementation that
 *        will replace the current geometry in the future
 */

 #ifndef NEW_GEOMETRY_H
 #define NEW_GEOMETRY_H

static inline int lexi(int b0, int b1, int b2, int b3, int x0, int x1, int x2, int x3);
static int index_blocked(int b0, int b1, int b2, int b3, 
                         int X0, int X1, int X2, int X3, 
                         int x0, int x1, int x2, int x3,
                         int sign);
static void gd_free_mem(geometry_descriptor *gd);
static void gd_alloc_mem(geometry_descriptor *gd, int N);
static void gd_set_copy(geometry_descriptor *gd);
static void index_alloc();
static void index_free();
static void geometry_index_init();
static void geometry_index_free();
char invertMask(char mask);
void freeGeometryBoxes();
int safe_ipt_ext(int x0_ext, int x1_ext, int x2_ext, int x3_ext);
int nearProc(char mask);
static void enumerate_lattice();
void sync_field(geometry_descriptor *gd, int bytes_per_site, int is_spinor_like, void *latticebuf);
void define_geometry();
void resetArray(char *array, int len);
int checkArray(char *array, int len, int local_len);
int isLocal(int x0_ext, int x1_ext, int x2_ext, int x3_ext);
int checkMask(char mask, int x0_ext, int x1_ext, int x2_ext, int x3_ext);

 #endif