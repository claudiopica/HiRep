#REPRESENTATION = REPR_FUNDAMENTAL
#REPRESENTATION = REPR_ADJOINT
#REPRESENTATION = REPR_ANTISYMMETRIC
#REPRESENTATION = REPR_SYMMETRIC

TOPDIR = ../../..

SRCS = 

MKDIR = $(TOPDIR)/Make
include $(MKDIR)/MkRules

MACRO = -D_${REPRESENTATION}_ -D_PICA_STYLE_ -DNDEBUG
EXEC = writeREPR

$(EXEC): $(MKDIR)/MkFlags main.cc adjoint.h antisymmetric.h complex.h fundamental.h list.h matrix.h polynomial.h print.h representation.h sparse.h style.h sun.h symmetric.h
	$(CXX) -Wall main.cc -o ${EXEC} ${MACRO}

clean::
	rm -f writeREPR
