TOPDIR = ../../

MAIN := chroma_converter

default: chroma_converter

MKDIR=$(TOPDIR)/Make
include $(MKDIR)/MkRules

QDIR=/opt/usqcd/single/
QOPTOP=${QDIR}/qop
QDPTOP=${QDIR}/qdp
QLATOP=${QDIR}/qla
QIOTOP=${QDIR}/qio
QMPTOP=${QDIR}/qmp

MODDIR := Geometry
MODS := geometry_init geometry_mpi_eo communications
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := archive logger read_input
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := suN_utils time boundary_conditions
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Observables
MODS := avr_plaquette
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss ranlxs ranlxd
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

main.o: main.c qdp-config.h
	echo -n Compiling [`pwd`/main.o]... ;\
	$(CC) -DNG=${NG} $(CFLAGS) -I${TOPDIR}/Include -DQDP_PROFILE -I. -I${QOPTOP}/include -I${QDPTOP}/include -I${QLATOP}/include -I${QIOTOP}/include -I${QMPTOP}/include -I${QLATOP}/include -I${QIOTOP}/include -I${QMPTOP}/include -c -o main.o main.c &&\
	echo " done." || echo " failed."

xml_utils.o: xml_utils.c qdp-config.h
	echo -n Compiling [`pwd`/xml_utils.o]... ;\
	$(CC) -DNG=${NG} $(CFLAGS) -c -o xml_utils.o xml_utils.c &&\
	echo " done." || echo " failed."

plaquette.o: plaquette.c qdp-config.h
	echo -n Compiling [`pwd`/plaquette.o]... ;\
	$(CC) -DNG=${NG} $(CFLAGS) -DQDP_PROFILE -I. -I${QOPTOP}/include -I${QDPTOP}/include -I${QLATOP}/include -I${QIOTOP}/include -I${QMPTOP}/include -I${QLATOP}/include -I${QIOTOP}/include -I${QMPTOP}/include -c -o plaquette.o plaquette.c &&\
	echo " done." || echo " failed."

read-gauge.o: read-gauge.c qdp-config.h
	echo -n Compiling [`pwd`/read-gauge.o]... ;\
	$(CC) -DNG=${NG} $(CFLAGS) -DQDP_PROFILE -I. -I${QOPTOP}/include -I${QDPTOP}/include -I${QLATOP}/include -I${QIOTOP}/include -I${QMPTOP}/include -I${QLATOP}/include -I${QIOTOP}/include -I${QMPTOP}/include -c -o read-gauge.o read-gauge.c &&\
	echo " done." || echo " failed."


translate.o: translate.c qdp-config.h
	echo -n Compiling [`pwd`/translate.o]... ;\
	$(CC) -DNG=${NG} $(CFLAGS) -I${TOPDIR}/Include -I${TOPDIR}/include -DQDP_PROFILE -I. -I${QOPTOP}/include -I${QDPTOP}/include -I${QLATOP}/include -I${QIOTOP}/include -I${QMPTOP}/include -I${QLATOP}/include -I${QIOTOP}/include -I${QMPTOP}/include -c -o translate.o translate.c &&\
	echo " done." || echo " failed."


chroma_converter: $(OBJS) $(MOD_OBJS) main.o plaquette.o read-gauge.o translate.o xml_utils.o
	echo -n Linking [$@]... ;\
	$(CC) -DNG=${NG} $(CFLAGS) $(LDFLAGS) -DQDP_PROFILE -I. -I${TOPDIR}/Include -I${QDPTOP}/include -I${QIOTOP}/include -I${QLATOP}/include -I${QMPTOP}/include -I${QOPTOP}/include -o $@ $(OBJS) $(MOD_OBJS) main.o xml_utils.o plaquette.o read-gauge.o translate.o -L${QOPTOP}/lib -lqopqdp  -L${QDPTOP}/lib -L${QLATOP}/lib -L${QIOTOP}/lib -L${QMPTOP}/lib -lqdp_f${NG} -lqdp_d${NG} -lqdp_df${NG} -lqla_f${NG} -lqla_d${NG} -lqla_q${NG} -lqla_df${NG} -lqla_dq${NG} -lqdp_f -lqdp_d -lqdp_df -lqla_f -lqla_d -lqla_q -lqla_df -lqla_dq -lqdp_common -lqdp_int -lqla_int -lqla_random -lqla_cmath -lqla_c99 -lqio -llime -lqmp -lm -L${QLATOP}/lib -lqla_f${NG} -lqla_d${NG} -lqla_q${NG} -lqla_df${NG} -lqla_dq${NG} -lqla_f -lqla_d -lqla_q -lqla_df -lqla_dq -lqla_int -lqla_random -lqla_cmath -lm -lqla_c99 -L${QIOTOP}/lib -lqio -llime -L${QMPTOP}/lib -lqmp -lm &&\
        echo " done." || echo " failed."


clean::
	rm -rf *.o *~ 
