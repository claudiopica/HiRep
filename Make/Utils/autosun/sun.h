string infinitesimal_evolution(const char* vname, const char* hname, const char* uname, const char* dtname);
string fundamental_algebra_represent(const char* mname, const char* hname);
string fundamental_algebra_project(const char* hname, const char* mname);


namespace group
{
	int N;
	int DIM;

	smatrix* T;
	string name;
	FLOATING Tnorm;
	
	void init(int n);
};


void group::init(int n)
{
#ifdef NDEBUG 
#ifdef _GAUGE_SON_
	cerr << " Initializing group SO(" << n << ")..... ";
#else
	cerr << " Initializing group SU(" << n << ")..... ";
#endif 
#endif 

	int A;
	int a, b;
	
	N = n;

#ifdef _GAUGE_SON_
	DIM = N*(N-1)/2;
	T = new smatrix[group::DIM];
	
	A = 0;
	for(a = 0; a < N; a++) for(b = a+1; b < N; b++){
	    T[A].size = N;
	    T[A].set(a,b, complex(.0,1.));
	    T[A].set(b,a, complex(.0,-1.));
	    A++;
	  }

#else
	DIM = N*N-1;
	T = new smatrix[group::DIM];
	
	A = 0;
	for(a = 0; a < N; a++) for(b = 0; b < N; b++)
		if(a > b)
		{
			T[A].size = N;
			T[A].set(a,b, complex(1.,.0));
			T[A].set(b,a, complex(1.,.0));
			A++;
		}
		else if(a < b)
		{
			T[A].size = N;
			T[A].set(a,b, complex(.0,1.));
			T[A].set(b,a, complex(.0,-1.));
			A++;
		}
		else if(a == b && a != 0)
		{
			T[A].size = N;
			for(int k = 0; k < a; k++)
				T[A].set(k,k, complex(sqrt(2./(a*(a+1.))),.0));
			T[A].set(a,a, complex(-a*sqrt(2./(a*(a+1.))),.0));
			A++;
		}
#endif
		
	Tnorm = 2.0;

        //my changes below
        for (A=0;A<group::DIM;A++)
           T[A].scale(1.0/(Tnorm));
        Tnorm=0.5;
        
#ifndef NDEBUG 
	cerr << "OK\n";
#endif
}


string infinitesimal_evolution(const char* vname, const char* hname, const char* uname, const char* dtname)
{
	string RET;
	
	rvector H(group::DIM,hname);
	cmatrix U(group::N,uname);
	pmatrix M(group::N);
	pmatrix V(group::N);
	rvariable dt(dtname);
	
	dt.scale(complex(0.0,1.0));
	H.scale(dt);
	
	for(int A = 0; A < group::DIM; A++)
	{
		pmatrix T(group::T[A]);
		T.scale(H.get(A));
		M.add(T);
	}
	
	V.mult(M, U);
	
	RET = V.assignment("+=", vname);
	
	return RET;
}


string ExpX(const char* dtname,  const char* hname, const char* uname)
{
	ostringstream RET;
	
	rvector H(group::DIM,hname);
	pmatrix M(group::N);
	rvariable dt(dtname);
	
	dt.scale(complex(0.0,1.0));
	H.scale(dt);
	
	for(int A = 0; A < group::DIM; A++)
	{
		pmatrix T(group::T[A]);
		T.scale(H.get(A));
		M.add(T);
	}
	
	RET << 
	"\tdouble y[3];\n" << 
	"\tdouble s[" << group::N*(group::N-1)/2 << "][4];\n";

#ifdef _GAUGE_SON_
	RET <<
	  "\tsuNgc ut, *u;\n\n"
	  "\tfor (int i=0; i<NG*NG; ++i) { ut.c[i].re=r->c[i]; ut.c[i].im=0.; }\n"
	  "\tu=&ut;\n\n";
#endif
	
	int k = 0;
	for(int j = 1; j < group::N; j++)
		for(int i = 0; i < j; i++)
		{
			polynomial tmp;
			pconstant ntmp(1.0/group::N);
			tmp = M.get(j,j);
			tmp.minus();
			tmp += M.get(i,i);
			tmp *= ntmp;
			RET <<
			"\ty[0] = " << M.get(i,j).str_imag() << ";\n" <<
			"\ty[1] = " << M.get(i,j).str_real() << ";\n" <<
			"\ty[2] = " << tmp.str_imag() << ";\n" <<
			"\tYtoU(s[" << k << "],y);\n";
			for(int p = 0; p < group::N; p++)
				RET << "\tsu2_rotate(s[" << k << "],&("
					<< uname << mindex(i,p,group::N) << "),&("
					<< uname << mindex(j,p,group::N) << "));\n";
			k++;
		}

	k = group::N*(group::N-1)/2 - 1;
	for(int j = group::N-1; j >= 1; j--)
		for(int i = j-1; i >= 0; i--)
		{
			for(int p = 0; p < group::N; p++)
				RET << "\tsu2_rotate(s[" << k << "],&("
					<< uname << mindex(i,p,group::N) << "),&("
					<< uname << mindex(j,p,group::N) << "));\n";
			k--;
		}
#ifdef _GAUGE_SON_
    RET<<"\n\tfor (int i=0; i<NG*NG; ++i) { r->c[i]=ut.c[i].re; }\n";
#endif
	
	return RET.str();
}


string fundamental_algebra_represent(const char* mname, const char* hname)
{
	string RET;
	rvector H(group::DIM,hname);
	pmatrix M(group::N);
	pconstant I(complex(0.0,1.0));
	
	for(int A = 0; A < group::DIM; A++)
	{
		pmatrix iT(group::T[A]);
		iT.scale(H.get(A));
		iT.scale(I);
		M.add(iT);
	}
	
	RET = M.assignment("=", mname);

	return RET;
}


string fundamental_algebra_project(const char* hname, const char* mname)
{
	string RET;
	pvector H(group::DIM);
	pmatrix *M;
//	pmatrix adjM(group::N);
	pconstant I(complex(0.0,1.0));

#ifdef _GAUGE_SON_
	M = new rmatrix(group::N,mname);
#else
	M = new cmatrix(group::N,mname);
#endif

//	adjM = *M;
//	adjM.adjoint();
//	M->sub(adjM);
	
	for(int A = 0; A < group::DIM; A++)
	{
		pmatrix iT(group::T[A]);
		iT.scale(I);
		polynomial iTM;
		herm(iTM,iT,*M);
		iTM.real();
		iTM.scale(1.0/group::Tnorm);
		H.set(A, iTM);
	}
	
	RET = H.assignment("=", hname);

	delete M;
	
	return RET;
}

