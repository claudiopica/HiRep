namespace representation
{
	int DIM;
	const int PHI_FLAVORS = 4;
	typedef complex TYPE;

	smatrix* iT;
	string name;
	FLOATING iTnorm;

	static smatrix* e;

	void init();
};


void representation::init()
{
#ifndef NDEBUG
	cerr << "Initializing ANTISYMMETRIC representation..... ";
#endif

	int N = group::N;
	
	if(N == 0)
	{
		cerr << "Initialization of group needed.";
		exit(1);
	}

	int A, B, C;
	smatrix tmp(N), tmp1(N);
	
	name = "ANTISYMMETRIC";
	DIM = N*(N-1)/2;
	iT = new smatrix[group::DIM];
	e = new smatrix[DIM];
	
	C = 0;
	for(A = 1; A < N; A++)
		for(B = 0; B < A; B++)
		{
			e[C].size = N;
			e[C].set(A,B, complex(sqrt(.5),0.0));
			e[C].set(B,A, complex(-sqrt(.5),0.0));
			C++;
		}
	
	for(C = 0; C < group::DIM; C++)
	{
		iT[C].size = DIM;
		for(A = 0; A < DIM; A++)
		{
			tmp.mult(e[A], group::T[C]);
			for(B = 0; B < DIM; B++)
			{
				tmp1.mult(tmp, e[B]);
				complex ctmp;
				trace(ctmp, tmp1);
				ctmp *= complex(0.0,-2.0);
				iT[C].set(A,B, ctmp);
			}
		}
	}
	
	iTnorm = (N-2.)*group::Tnorm;

#ifndef NDEBUG
	cerr << "OK\n";
#endif
}


string group_represent(const char* vname, const char* uname)
{
	string RET;
	cmatrix U(group::N,uname);
	pmatrix trU(group::N);
	pmatrix rU(representation::DIM);
	pmatrix *Ue;
    
    Ue = new pmatrix[representation::DIM];
	
	trU = U;
	trU.transpose();

	for(int A = 0; A < representation::DIM; A++)
	{
		pmatrix e(representation::e[A]);
		pmatrix mtmp(group::N);
		mtmp.mult(e, trU);
		Ue[A].mult(U, mtmp);
	}
	
	for(int B = 0; B < representation::DIM; B++)
		for(int A = 0; A < representation::DIM; A++)
		{
			polynomial v;
			pmatrix mtmp(group::N);
			pmatrix e(representation::e[A]);
			e.adjoint();
			mtmp.mult(e, Ue[B]);
			trace(v, mtmp);
			rU.set(A,B, v);
		}

	RET += rU.assignment("=", vname);
	
    delete[] Ue;
    
	return RET;
}

string debug_group_represent(const char* vname, const char* uname)
{
	string RET = string("copy(") + vname + "," + uname + ");\n\
	int A, C, a, b, c, d;\n\
	A = 0;\n\
	for(a = 1; a < GROUP::N; a++) for(b = 0; b < a; b++) {\n\
		C = 0;\n\
		for(c = 1; c < GROUP::N; c++) for(d = 0; d < c; d++) {\n\
			RET(A,C) = U(a,c)*U(b,d)-U(a,d)*U(b,c);\n\
			C++;\n\
		}\n\
		A++;\n\
	}\n";
	return RET;
}
