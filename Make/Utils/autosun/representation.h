string group_represent(const char* vname, const char* uname);
string algebra_represent(const char* mname, const char* hname);
string algebra_project(const char* hname, const char* mname);


string algebra_represent(const char* mname, const char* hname)
{
	string RET;
	rvector H(group::DIM,hname);
	pmatrix M(representation::DIM);
	
	for(int A = 0; A < group::DIM; A++)
	{
		pmatrix iT(representation::iT[A]);
		iT.scale(H.get(A));
		M.add(iT);
	}
	
	RET = M.assignment("=", mname);

	return RET;
}


string algebra_project(const char* hname, const char* mname)
{
	string RET;
	pvector H(group::DIM);
	pmatrix *M;
	//pmatrix adjM(representation::DIM);

	if(sizeof(representation::TYPE)==sizeof(FLOATING))
		M = new rmatrix(representation::DIM,mname);
	else
#ifdef _GAUGE_SON_
	        M = new rmatrix(representation::DIM,mname);
#else
		M = new cmatrix(representation::DIM,mname);
#endif

	//adjM = *M;
	//adjM.adjoint();
	//M->sub(adjM);
	
	for(int A = 0; A < group::DIM; A++)
	{
		pmatrix iT(representation::iT[A]);
		polynomial iTM;
		herm(iTM,iT,*M);
		iTM.real();
		iTM.scale(1.0/representation::iTnorm);
		H.set(A, iTM);
	}
	
	RET = H.assignment("=", hname);

	delete M;
	
	return RET;
}
