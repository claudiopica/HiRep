#define FLOATING double


class complex
{
public:
	FLOATING re, im;

	complex() { re = 0.0; im = 0.0; }
	complex(const FLOATING& r) { re = r; im = 0.0; }
	complex(const FLOATING& r, const FLOATING& i) { re = r; im = i; }
	complex(const complex& z) { re = z.re; im = z.im; }

	FLOATING real() { return re; }
	FLOATING imag() { return im; }

	void clear()
		{ re = 0.0; im = 0.0; }
	void minus()
		{ re = -re; im = -im; }
	void conjugate()
		{ im = -im; }
	void add(const complex& a)
		{ *this += a; }
	void mult(const complex& a, const complex& b)
		{ *this = a*b; }
	void add_mult(const complex& a, const complex& b)
		{ *this += a*b; }

	complex operator=(const complex& a)
		{ re = a.re; im = a.im; return *this; }
	complex operator=(const FLOATING& a)
		{ re = a; im = 0.0; return *this; }

	complex operator+=(const complex& a)
		{ re += a.re; im += a.im; return *this; }
	complex operator+=(const FLOATING& a)
		{ re += a; return *this; }

	complex operator-=(const complex& a)
		{ re -= a.re; im -= a.im; return *this; }
	complex operator-=(const FLOATING& a)
		{ re -= a; return *this; }

	complex operator*=(const complex& a)
		{ FLOATING tmp = re*a.re-im*a.im; im = re*a.im+im*a.re; re = tmp; return *this; }
	complex operator*=(const FLOATING& a)
		{ re *= a; im *= a; return *this; }

	friend complex operator-(const complex& a)
		{ return complex(-a.re,-a.im); }

	friend bool operator==(const complex& a, const complex& b)
		{ return a.re==b.re && a.im==b.im; }
	friend bool operator==(const FLOATING& a, const complex& b)
		{ return a==b.re && 0.==b.im; }
	friend bool operator==(const complex& a, const FLOATING& b)
		{ return a.re==b && a.im==0.; }

	friend bool operator!=(const complex& a, const complex& b)
		{ return a.re!=b.re || a.im!=b.im; }
	friend bool operator!=(const FLOATING& a, const complex& b)
		{ return a!=b.re || 0.!=b.im; }
	friend bool operator!=(const complex& a, const FLOATING& b)
		{ return a.re!=b && a.im!=0.; }

	friend complex operator+(const complex& a, const complex& b)
		{ return complex(a.re+b.re,a.im+b.im); }
	friend complex operator+(const FLOATING& a, const complex& b)
		{ return complex(a+b.re,b.im); }
	friend complex operator+(const complex& a, const FLOATING& b)
		{ return complex(a.re+b,a.im); }

	friend complex operator-(const complex& a, const complex& b)
		{ return complex(a.re-b.re,a.im-b.im); }
	friend complex operator-(const FLOATING& a, const complex& b)
		{ return complex(a-b.re,-b.im); }
	friend complex operator-(const complex& a, const FLOATING& b)
		{ return complex(a.re-b,a.im); }

	friend complex operator*(const complex& a, const complex& b)
		{ return complex(a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re); }
	friend complex operator*(const FLOATING& a, const complex& b)
		{ return complex(a*b.re,a*b.im); }
	friend complex operator*(const complex& a, const FLOATING& b)
		{ return complex(a.re*b,a.im*b); }
 
	friend complex operator/(const complex& a, const complex& b)
	{
		FLOATING den = b.re*b.re+b.im*b.im;
		return complex(
			(a.re*b.re+a.im*b.im)/den,
			(-a.re*b.im+a.im*b.re)/den
		);
	}
	friend complex operator/(const FLOATING& a, const complex& b)
	{
		FLOATING den = b.re*b.re+b.im*b.im;
		return complex(
			a*b.re/den,
			-a*b.im/den
		);
	}
	friend complex operator/(const complex& a, const FLOATING& b)
		{ return complex(a.re/b,a.im/b); }
   
	friend ostream& operator<<(ostream& os, const complex& z)
	{
		os << "(" << z.re << "," << z.im << ")";
		return os;
	}
	
	friend complex conj(const complex& z) { return complex(z.re,-z.im); }
	friend FLOATING abs(const complex& z) { return sqrt(z.re*z.re+z.im*z.im); }
	friend FLOATING arg(const complex& z) { return atan2(z.im, z.re); }
};


string ftos(FLOATING x)
{
	static char tmp[100];
	if(x >= 0.0) sprintf(tmp, _PNUMBER_, x);
	else sprintf(tmp, _NUMBER_, x);
	return string(tmp);
}
