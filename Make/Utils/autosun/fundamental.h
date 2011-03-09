namespace representation
{
	int DIM;
	const int PHI_FLAVORS = 4;
	typedef complex TYPE;

	smatrix* iT;
	string name;
	FLOATING iTnorm;

	void init();
};


void representation::init()
{
#ifndef NDEBUG
	cerr << "Initializing FUNDAMENTAL representation..... ";
#endif

	int A;
	//int a, b;
	
	int N = group::N;
	
	if(N == 0)
	{
		cerr << "Initialization of group needed.";
		exit(1);
	}
	
	name = "FUNDAMENTAL";
	DIM = N;
	iT = new smatrix[group::DIM];
	
	for(A = 0; A < group::DIM; A++)
	{
		iT[A] = group::T[A];
		iT[A].scale(complex(0.0,1.0));
	}
	
	iTnorm = group::Tnorm;

#ifndef NDEBUG
	cerr << "OK\n";
#endif
}

string group_represent(const char* vname, const char* uname)
{
	string RET;
	cmatrix U(group::N,uname);

	RET += U.assignment("=", vname);
	
	return RET;
}

string debug_group_represent(const char* vname, const char* uname)
{
	string RET;
	RET = string("copy(") + vname + "," + uname + ");\n";
	return RET;
}
