#REPRESENTATION = REPR_FUNDAMENTAL
#REPRESENTATION = REPR_ADJOINT
#REPRESENTATION = REPR_ANTISYMMETRIC
#REPRESENTATION = REPR_SYMMETRIC

TOPDIR = ../../..
MKDIR = $(TOPDIR)/Make

#SRCS = 

EXEC = writeREPR

$(EXEC): $(MKDIR)/MkFlags main.cc adjoint.h antisymmetric.h complex.h fundamental.h list.h matrix.h polynomial.h print.h representation.h sparse.h style.h sun.h symmetric.h
	echo -n "Compiling [$(EXEC)]..." ;\
	$(CXX) main.cc -o ${EXEC} ${MACRO} &&\
	echo " done." || ( echo " failed." && false )

include $(MKDIR)/MkRules

MACRO = -D_${REPR}_ -D_${GAUGE_GROUP}_ -D_PICA_STYLE_ -DNDEBUG

clean::
	rm -f writeREPR
