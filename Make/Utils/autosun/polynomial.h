#define _ZERO_ (1e-10)

bool equal(FLOATING a,FLOATING b)
{
	return fabs(a-b) < _ZERO_;
}


class rmonomial : public orderedlist<string,int>
{
public:
	int& getzero() const
	{
		static int zero = 0;
		return zero;
	}
	
	using orderedlist<string,int>::operator=;
	
	friend bool operator<(const rmonomial& a, const rmonomial& b) { return a.getstring() < b.getstring(); }
	friend bool operator<=(const rmonomial& a, const rmonomial& b) { return a.getstring() <= b.getstring(); }
	friend bool operator==(const rmonomial& a, const rmonomial& b) { return a.getstring() == b.getstring(); }

	virtual string getstring() const
	{
		string RET;
		for(KEYTYPE i = 0; i < length; i++)
			for(int j = 0; j < data[i]->value; j++)
			{
				if(i+j != 0)
					RET += "*";
				RET += data[i]->index;
			}
		return RET;
	}

	void mult(const rmonomial& a, const rmonomial& b)
	{
		clear();
		for(KEYTYPE i = 0; i < a.length; i++)
			add(a[i].index, a[i].value);
		for(KEYTYPE i = 0; i < b.length; i++)
			add(b[i].index, b[i].value);
	}
};


class polynomial : public orderedlist<rmonomial,complex>
{
public:
	polynomial() : orderedlist<rmonomial,complex>() {}
	polynomial(const rmonomial& index, const complex& value) : orderedlist<rmonomial,complex>(index, value) {}
	
	complex& getzero() const
	{
		static complex zero = complex(0.0,0.0);
		return zero;
	}
	
	using orderedlist<rmonomial,complex>::operator=;

	polynomial& operator+=(const polynomial& b)
	{
		for(KEYTYPE i = 0; i < b.length; i++)
			add(b.data[i]->index, b.data[i]->value);
		return *this;
	}
	polynomial& operator*=(const polynomial& b)
	{
		polynomial tmp;
		tmp.mult(*this, b);
		*this = tmp;
		return *this;
	}

	bool isreal() const
	{
		bool RET = true;
		for(KEYTYPE i = 0; i < length; i++)
			RET = RET && data[i]->value.im == 0.0;
		return RET;
	}
	
	string str_real() const
	{
		string RET;
		FLOATING a = 0.0;
		int counter = 0;
		
		for(KEYTYPE i = 0; i < length; i++)
		{
			if(equal(data[i]->value.re,0.0)) continue;
			if(!equal(a,0.0) && !equal(a,data[i]->value.re) && !equal(a,-data[i]->value.re)) a = 1.0;
			else a = data[i]->value.re;
			counter++;
		}
		if(counter <= 1) a = 1.0;

		if(equal(a,-1.0)) RET = "-(";
		else if(!equal(a,1.0)) RET = ftos(a) + "*(";
		
		for(KEYTYPE i = 0; i < length; i++)
			if(!equal(data[i]->value.re,0.0))
			{
				FLOATING ct = data[i]->value.re/a;
				string coeff = ftos(ct);
				string rmon = data[i]->index.getstring();
				if(equal(ct,1.) && rmon != "") RET += "+";
				else if(equal(ct,-1) && rmon != "") RET += "-";
				else if (rmon != "") RET += coeff + "*";
				else RET += coeff;
				RET += rmon;
			}
		
		if(!equal(a,1.0)) RET += ")";
		
		if(RET == "") RET = "0.0";

		return RET;
	}
	string str_imag() const
	{
		string RET;
		FLOATING a = 0.0;
		int counter = 0;

		for(KEYTYPE i = 0; i < length; i++)
		{
			if(equal(data[i]->value.im,0.0)) continue;
			if(!equal(a,0.0) && !equal(a,data[i]->value.im) && !equal(a,-data[i]->value.im)) a = 1.0;
			else a = data[i]->value.im;
			counter++;
		}
		if(counter <= 1) a = 1.0;

		if(equal(a,-1.0)) RET = "-(";
		else if(!equal(a,1.0)) RET = ftos(a) + "*(";

		for(KEYTYPE i = 0; i < length; i++)
			if(!equal(data[i]->value.im,0.0))
			{
				FLOATING ct = data[i]->value.im/a;
				string coeff = ftos(ct);
				string rmon = data[i]->index.getstring();
				if(equal(ct,1.) && rmon != "") RET += "+";
				else if(equal(ct,-1) && rmon != "") RET += "-";
				else if (rmon != "") RET += coeff + "*";
				else RET += coeff;
				RET += rmon;
			}
		
		if(!equal(a,1.0)) RET += ")";
		
		if(RET == "") RET = "0.0";

		return RET;
	}

	virtual string getstring() const
	{
		return "(" + str_real() + "," + str_imag() + ")";
	}
	
	void minus()
	{
		for(int i = 0; i < length; i++)
			data[i]->value = -data[i]->value;
	}
	void conjugate()
	{
		for(int i = 0; i < length; i++)
			data[i]->value = conj(data[i]->value);
	}
	void real()
	{
		for(int i = 0; i < length; i++)
			data[i]->value.im = 0.0;
	}
	void mult(const polynomial& a, const polynomial& b)
	{
		clear();
		for(KEYTYPE i = 0; i < a.length; i++)
			for(KEYTYPE j = 0; j < b.length; j++)
			{
				rmonomial index;
				index.mult(a[i].index, b[j].index);
				add(index, a[i].value*b[j].value);
			}
	}
	void add_mult(const polynomial& a, const polynomial& b)
	{
		for(KEYTYPE i = 0; i < a.length; i++)
			for(KEYTYPE j = 0; j < b.length; j++)
			{
				rmonomial index;
				index.mult(a[i].index, b[j].index);
				add(index, a[i].value*b[j].value);
			}
	}
};


class rvariable : public polynomial
{
public:
	rvariable(const char* name) : polynomial()
	{
		rmonomial index;
		index.add(string(name), 1);
		orderedlist<rmonomial,complex>::add(index,complex(1.0,0.0));
	}

	using polynomial::operator=;
};


class cvariable : public polynomial
{
public:
	cvariable(const char* name) : polynomial()
	{
		rmonomial index1;
		index1.add(string(name)+".re", 1);
		orderedlist<rmonomial,complex>::add(index1,complex(1.0,0.0));
		rmonomial index2;
		index2.add(string(name)+".im", 1);
		orderedlist<rmonomial,complex>::add(index2,complex(0.0,1.0));
	}

	using polynomial::operator=;
};


class pconstant : public polynomial
{
public:
	pconstant(const double z) : polynomial()
	{
		rmonomial index;
		orderedlist<rmonomial,complex>::add(index,complex(z,0.0));
	}
	pconstant(const complex& z) : polynomial()
	{
		rmonomial index;
		orderedlist<rmonomial,complex>::add(index,z);
	}

	using polynomial::operator=;
};
