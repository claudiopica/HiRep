CPP=clang++
CFLAGS= -O3 -Wall -Wextra
LDFLAGS= -L/opt/local/lib -lefence 
INC = -I ./.

# Location of Gnu Science Library
GSL_LIB=/opt/local


all: gevp test_mat  test_eig 

gevp: gevp.cc read_arg.o utils.o matrix_algebra.o fit.o type.h Makefile
	${CPP} fit.o read_arg.o utils.o matrix_algebra.o ${CFLAGS} ${LDFLAGS} $< -o $@

read_arg.o: read_arg.cc read_arg.h type.h utils.h Makefile
	${CPP} ${CFLAGS} -c read_arg.cc

utils.o: utils.cc type.h utils.h Makefile
	${CPP} ${CFLAGS} -c utils.cc

fit.o: fit.cc fit.h type.h Makefile
	${CPP} ${CFLAGS} -c fit.cc

matrix_algebra.o: matrix_algebra.cc matrix_algebra.h  Makefile
	${CPP} ${CFLAGS} -c matrix_algebra.cc

test_mat: test_mat.cc matrix_algebra.o Makefile
	${CPP} matrix_algebra.o  ${CFLAGS} $< -o $@

test_eig: test_eig.cc matrix_algebra.o Makefile
	${CPP} -o $@  -I ${GSL_LIB}/include/   test_eig.cc matrix_algebra.o utils.o  -L ${GSL_LIB}/lib  -lm -lgsl -lgslcblas -std=c++11

clean:
	@rm -f gevp *~ *.o test_eig test_mat
