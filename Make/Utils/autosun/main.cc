#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace std;


#include "style.h"

#include "list.h"
#include "complex.h"
#include "sparse.h"
#include "polynomial.h"
#include "matrix.h"
#include "sun.h"

#ifdef _REPR_FUNDAMENTAL_
#include "fundamental.h"
#elif _REPR_ADJOINT_
#include "adjoint.h"
#elif _REPR_ANTISYMMETRIC_
#include "antisymmetric.h"
#elif _REPR_SYMMETRIC_
#include "symmetric.h"
#endif

#include "representation.h"
#include "print.h"

int main(int argc, char* argv[]) {
	if(argc != 3) {
		cerr << "Usage: "<<argv[0]<<" N tmplfile\n";
		cerr << "N => number of colors\n";
		cerr << "tmpldir => template file\n";
		return 1;
	}
	
	int N = atoi(argv[1]);
	if(N < 2) {
		cerr << argv[0] <<": ERROR: number of color must be >1 !\n";
		return 2;
	}
	
	group::init(N);
	representation::init();

#ifndef NDEBUG
	cerr << argv[0] << ": Opening file..... ";
#endif
//#ifdef _PICA_STYLE_
	//ostringstream h_out_name;
	//h_out_name << "SU" << group::N << "_" << representation::name << ".h";
	/*
	string h_out_name = "suN_repr_func.h";
	ofstream h_out(h_out_name.c_str());
	if(!h_out) {
		cerr << "ERROR (" << h_out_name.c_str() << ")\n";
		return 0;
	}
	*/
	string h_tmpl_name = string(argv[2]);// + string(_TMPL_FILE_) + ".h.tmpl";
	ifstream h_tmpl(h_tmpl_name.c_str());
	if(!h_tmpl) {
		cerr << argv[0]<<": ERROR: cannot open template file (" << h_tmpl_name << ")!\n";
		return 3;
	}
	
#ifndef NDEBUG
	cerr << "OK\n";
#endif
	
	printfile(cout, h_tmpl);

	/*
	printfile(h_out, h_tmpl);
	
	h_out.close();
	h_tmpl.close();

	ofstream h_out2("suN_exp.c");
	if(!h_out2) {
		cerr << "ERROR (" << "suN_exp.c" << ")\n";
		return 0;
	}

	string h_tmpl_name2 = string(argv[2]) + "suN_exp.c.tmpl";
	ifstream h_tmpl2(h_tmpl_name2.c_str());
	if(!h_tmpl2) {
		cerr << "ERROR (" << h_tmpl_name << ")\n";
		return 0;
	}
	
	cerr << "OK\n";
	
	printfile(h_out2, h_tmpl2);
	
	h_out2.close();
	h_tmpl2.close();
	*/
//#endif
	
	return 0;
}
