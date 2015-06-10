namespace representation
{
	int DIM;
	const int PHI_FLAVORS = 2;
	typedef FLOATING TYPE;

	smatrix* iT;
	string name;
	FLOATING iTnorm;

	static smatrix* e;

	void init();
};


void representation::init()
{
#ifndef NDEBUG
	cerr << "Initializing ADJOINT representation..... ";
#endif
	int N = group::N;
	
	if(N == 0)
	{
		cerr << "Initialization of group needed.";
		exit(1);
	}

	int A, B, C;
	smatrix tmp(N), tmp1(N);

	name = "ADJOINT";
#ifdef _GAUGE_SON
        DIM = N*(N-1)/2;
#else
	DIM = N*N-1;
#endif
	iT = new smatrix[group::DIM];
	e = new smatrix[DIM];
	
	for(C = 0; C < DIM; C++)
	{
		e[C] = group::T[C];
		e[C].scale(1./sqrt(group::Tnorm));
	}
	
	for(C = 0; C < group::DIM; C++)
	{
		iT[C].size = DIM;
		for(A = 0; A < DIM; A++)
		{
			tmp.mult(e[A], group::T[C]);
			tmp.minus();
			tmp.add_mult(group::T[C], e[A]);
			
			for(B = 0; B < DIM; B++)
			{
				complex ctmp;
				tmp1.mult(tmp, e[B]);
				trace(ctmp, tmp1);
				ctmp.re = -ctmp.im;
				ctmp.im = 0.0;
				iT[C].set(B,A, ctmp);
			}
		}
	}
	
	iTnorm = 2.*group::Tnorm*N;

#ifndef NDEBUG
	cerr << "OK\n";
#endif
}


string group_represent(const char* vname, const char* uname)
{
	string RET;
	cmatrix U(group::N,uname);
	pmatrix adjU(group::N);
	pmatrix rU(representation::DIM);
	pmatrix *Ue;
    
        Ue = new pmatrix[representation::DIM];
	
	adjU = U;
	adjU.adjoint();

	for(int A = 0; A < representation::DIM; A++)
	{
		pmatrix e(representation::e[A]);
		pmatrix mtmp(group::N);
		mtmp.mult(e, adjU);
		Ue[A].mult(U, mtmp);
	}
	
	for(int B = 0; B < representation::DIM; B++)
		for(int A = 0; A < representation::DIM; A++)
		{
			polynomial v;
			pmatrix mtmp(group::N);
			pmatrix e(representation::e[A]);
			mtmp.mult(e, Ue[B]);
			trace(v, mtmp);
			v.real();
			rU.set(A,B, v);
		}
	
	RET += rU.assignment("=", vname);
	
    delete[] Ue;
    
	return RET;
}

string debug_group_represent(const char* vname, const char* uname)
{
	string RET = string("copy(") + vname + "," + uname + ");\n\
	int A, B;\n\
	GROUP::matrix e[DIM];\n\
	GROUP::matrix tmp[2];\n\
	REAL enorm2 = .5;\n\
	A = 0;\n\
	for(int a = 0; a < GROUP::N; a++) for(int b = 0; b < GROUP::N; b++)\n\
		if(a > b)\n\
		{\n\
			setzero(e[A]);\n\
			e[A](a,b) = COMPLEX(.5,0.);\n\
			e[A](b,a) = COMPLEX(.5,0.);\n\
			A++;\n\
		}\n\
		else if(a < b)\n\
		{\n\
			setzero(e[A]);\n\
			e[A](a,b) = COMPLEX(0.,.5);\n\
			e[A](b,a) = COMPLEX(0.,-.5);\n\
			A++;\n\
		}\n\
		else if(a == b && a != 0)\n\
		{\n\
			setzero(e[A]);\n\
			REAL n = sqrt(.5/(REAL)(a*(a+1)));\n\
			for(int k = 0; k < a; k++)\n\
				e[A](k,k) = n;\n\
			e[A](a,a) = -a*n;\n\
			A++;\n\
		}\n\
	for(A = 0; A < DIM; A++)\n\
	{\n\
		_mult(tmp[0], U, e[A]);\n\
		_mult_2adj(tmp[1], tmp[0], U);\n\
		for(B = 0; B < DIM; B++)\n\
			RET(B,A) = _herm_real(e[B],tmp[1])/enorm2;\n\
	}\n";
	return RET;
}

