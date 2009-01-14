#ifndef QDP_CONFIG_A
#define QDP_CONFIG_A

#define QDP_Precision 'F'

#define QDP_Nc NG

#include <stdio.h>

typedef struct{
  double beta;
  int X;
  int Y;
  int Z;
  int T;
  int Nc;
}params_struct;


#include <qdp.h>
#include <qop_qdp.h>

#define NDIM 4

#define printf0 if (QDP_this_node == 0) printf

QLA_Real plaquette(QDP_ColorMatrix *U[]);
void translate(QLA_ColorMatrix *m, int coords[]);
int read_gauge(QDP_ColorMatrix *U[], char *name);
void chroma_to_HiRep(QDP_ColorMatrix *link[]);
void read_tag(char* from, char* to, char * tag);

params_struct * params;


/*  extern int write_gauge(char *name, QDP_ColorMatrix *U[]); */



#endif
