#define _PARAMETER_CHECK(str) \
	if(str == string::npos) \
	{ \
		cerr << "ERROR. Bad parameters: $" << tag << "$.\n"; \
		return; \
	}



void printfile(ostream& out, istream& tmpl)
{
	while(!tmpl.eof())
	{
		string txt, tag;
		getline(tmpl, txt, '$');
		getline(tmpl, tag, '$');
		
		out << txt;
		if(tag == "N")
			out << group::N;
		else if(tag == "GROUP::DIM")
			out << group::DIM;
		else if(tag == "GROUP::NORM2")
			out << ftos(group::Tnorm);
		else if(tag.find("GROUP::ALGEBRA_REPRESENT") == 0)
		{
#ifndef NDEBUG
			cerr << "Writing fundamental_algebra_represent..... ";
#endif

			size_t start, end;
			
			start = tag.find(":", 8);
			_PARAMETER_CHECK(start)
			end = tag.find(":",start+1);
			_PARAMETER_CHECK(end)
			string mname = tag.substr(start+1, end-start-1);

			start = end;
			string hname = tag.substr(start+1);
			
			out << fundamental_algebra_represent(mname.c_str(), hname.c_str());
#ifndef NDEBUG
			cerr << "OK\n";
#endif
		}
		else if(tag.find("GROUP::ALGEBRA_PROJECT") == 0)
		{
#ifndef NDEBUG
			cerr << "Writing fundamental_algebra_project..... ";
#endif

			size_t start, end;
			
			start = tag.find(":", 8);
			_PARAMETER_CHECK(start)
			end = tag.find(":",start+1);
			_PARAMETER_CHECK(end)
			string hname = tag.substr(start+1, end-start-1);

			start = end;
			string mname = tag.substr(start+1);
			
			out << fundamental_algebra_project(hname.c_str(), mname.c_str());
#ifndef NDEBUG
			cerr << "OK\n";
#endif
		}
		else if(tag.find("GROUP::INFINITESIMAL") == 0)
		{
#ifndef NDEBUG
			cerr << "Writing infinitesimal_evolution..... ";
#endif

			size_t start, end;
			
			start = tag.find(":", 8);
			_PARAMETER_CHECK(start)
			end = tag.find(":",start+1);
			_PARAMETER_CHECK(end)
			string vname = tag.substr(start+1, end-start-1);

			start = end;
			end = tag.find(":",start+1);
			_PARAMETER_CHECK(end)
			string hname = tag.substr(start+1, end-start-1);

			start = end;
			end = tag.find(":",start+1);
			_PARAMETER_CHECK(end)
			string uname = tag.substr(start+1, end-start-1);

			start = end;
			string dtname = tag.substr(start+1);
			
			out << infinitesimal_evolution(vname.c_str(), hname.c_str(), uname.c_str(), dtname.c_str());
#ifndef NDEBUG
			cerr << "OK\n";
#endif
		}
		else if(tag.find("GROUP::EXPX") == 0)
		{
#ifndef NDEBUG
			cerr << "Writing ExpX..... ";
#endif

			size_t start, end;
			
			start = tag.find(":", 8);
			_PARAMETER_CHECK(start)
			end = tag.find(":",start+1);
			_PARAMETER_CHECK(end)
			string dtname = tag.substr(start+1, end-start-1);

			start = end;
			end = tag.find(":",start+1);
			_PARAMETER_CHECK(end)
			string hname = tag.substr(start+1, end-start-1);

			start = end;
			string uname = tag.substr(start+1);
			
			out << ExpX(dtname.c_str(), hname.c_str(), uname.c_str());
#ifndef NDEBUG
			cerr << "OK\n";
#endif
		}
		else if(tag == "REPR::DIM")
			out << representation::DIM;
		else if(tag == "REPR::NORM2")
			out << ftos(representation::iTnorm);
		else if(tag == "REPR::NAME")
			out << representation::name;
		else if(tag == "REPR::PHI_FLAVORS")
			out << representation::PHI_FLAVORS;
		else if(tag == "REPR::TYPE")
		{
			if(sizeof(representation::TYPE) == sizeof(float))
				out << "float";
			else if(sizeof(representation::TYPE) == sizeof(double))
				out << "double";
			else if(sizeof(representation::TYPE) == sizeof(complex))
				out << "COMPLEX";
		}
		else if(tag.find("REPR::GROUP_REPRESENT") == 0)
		{
#ifndef NDEBUG
			cerr << "Writing group_represent..... ";
#endif

			size_t start, end;
			
			start = tag.find(":", 8);
			_PARAMETER_CHECK(start)
			end = tag.find(":",start+1);
			_PARAMETER_CHECK(end)
			string vname = tag.substr(start+1, end-start-1);

			start = end;
			string uname = tag.substr(start+1);
			
			out << group_represent(vname.c_str(), uname.c_str());
#ifndef NDEBUG
			cerr << "OK\n";
#endif
		}
		else if(tag.find("REPR::DEBUG_GROUP_REPRESENT") == 0)
		{
#ifndef NDEBUG
			cerr << "Writing debug_group_represent..... ";
#endif

			size_t start, end;
			
			start = tag.find(":", 8);
			_PARAMETER_CHECK(start)
			end = tag.find(":",start+1);
			_PARAMETER_CHECK(end)
			string vname = tag.substr(start+1, end-start-1);

			start = end;
			string uname = tag.substr(start+1);
			
			out << debug_group_represent(vname.c_str(), uname.c_str());
#ifndef NDEBUG
			cerr << "OK\n";
#endif
		}
		else if(tag.find("REPR::ALGEBRA_REPRESENT") == 0)
		{
#ifndef NDEBUG
			cerr << "Writing algebra_represent..... ";
#endif
			size_t start, end;
			
			start = tag.find(":", 8);
			_PARAMETER_CHECK(start)
			end = tag.find(":",start+1);
			_PARAMETER_CHECK(end)
			string mname = tag.substr(start+1, end-start-1);

			start = end;
			string hname = tag.substr(start+1);
			
			out << algebra_represent(mname.c_str(), hname.c_str());
#ifndef NDEBUG
			cerr << "OK\n";
#endif
		}
		else if(tag.find("REPR::ALGEBRA_PROJECT") == 0)
		{
#ifndef NDEBUG
			cerr << "Writing algebra_project..... ";
#endif

			size_t start, end;
			
			start = tag.find(":", 8);
			_PARAMETER_CHECK(start)
			end = tag.find(":",start+1);
			_PARAMETER_CHECK(end)
			string hname = tag.substr(start+1, end-start-1);

			start = end;
			string mname = tag.substr(start+1);
			
			out << algebra_project(hname.c_str(), mname.c_str());
#ifndef NDEBUG
			cerr << "OK\n";
#endif
		}
		else if(tag.find("REPR::GAUSSIAN_RG") == 0)
		{
#ifndef NDEBUG
			cerr << "Writing gaussian_vector..... ";
#endif

			if(sizeof(representation::TYPE) == sizeof(complex))
				out << "COMPLEX(gaussian_rg(.5), gaussian_rg(.5))";
			else
				out << "gaussian_rg(.5)";

#ifndef NDEBUG
			cerr << "OK\n";
#endif
		}
		else if(tag != "")
		{
			cerr << "ERROR. Bad tag: $" << tag << "$.\n";
			return;
		}
		out.flush();
	}
}


#undef _PARAMETER_CHECK
