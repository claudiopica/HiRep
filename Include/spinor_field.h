typedef struct {
	suNf_spinor* ptr;
} spinor_field;

typedef struct {
	suNf_spinor_flt* ptr;
} spinor_field_flt;

#define _SPINOR_ADDR(s) ((s)->ptr)

/* x e' l'indice di sito, i e' l'indice di componente */
#define _SPINOR_AT_SITE(s,x) (((s)->ptr)+x)
#define _SPINOR_AT(s,i) (((s)->ptr)+i)
#define _SITE2SC(x) (x)
#define _SC2SITE(i) (i)

#define EVENSITES 0
#define ODDSITES 1
#define ALLSITES 2

#define FOR_LOCAL_SC(x,i) for((i)=0,(x)=(i);(i)<VOLUME;(i)++,(x)=(i))
#define FOR_EVEN_LOCAL_SC(x,i) for((i)=0,(x)=(i);(i)<VOLUME/2;(i)++,(x)=(i))
#define FOR_ODD_LOCAL_SC(x,i) for((i)=VOLUME/2,(x)=(i);(i)<VOLUME;(i)++,(x)=(i))

#define FOR_SOME_SC(block,x,i)                                  \
	for( (i) = ( ((block)==ODDSITES) ? VOLUME/2 : 0 ), (x)=(i); \
	     (i) < ( ((block)==EVENSITES) ? VOLUME/2 : VOLUME );    \
	     (i)++, (x)=(i)                                         \
	   )

#define FOR_BOUNDARY_SC
#define FOR_BUFFER_SC


