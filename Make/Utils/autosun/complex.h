#define FLOATING double

class complex_t {
public:
    FLOATING re, im;

    complex_t() {
        re = 0.0;
        im = 0.0;
    }
    complex_t(const FLOATING &r) {
        re = r;
        im = 0.0;
    }
    complex_t(const FLOATING &r, const FLOATING &i) {
        re = r;
        im = i;
    }
    complex_t(const complex_t &z) {
        re = z.re;
        im = z.im;
    }

    FLOATING real() {
        return re;
    }
    FLOATING imag() {
        return im;
    }

    void clear() {
        re = 0.0;
        im = 0.0;
    }
    void minus() {
        re = -re;
        im = -im;
    }
    void conjugate() {
        im = -im;
    }
    void add(const complex_t &a) {
        *this += a;
    }
    void mult(const complex_t &a, const complex_t &b) {
        *this = a * b;
    }
    void add_mult(const complex_t &a, const complex_t &b) {
        *this += a * b;
    }

    complex_t operator=(const complex_t &a) {
        re = a.re;
        im = a.im;
        return *this;
    }
    complex_t operator=(const FLOATING &a) {
        re = a;
        im = 0.0;
        return *this;
    }

    complex_t operator+=(const complex_t &a) {
        re += a.re;
        im += a.im;
        return *this;
    }
    complex_t operator+=(const FLOATING &a) {
        re += a;
        return *this;
    }

    complex_t operator-=(const complex_t &a) {
        re -= a.re;
        im -= a.im;
        return *this;
    }
    complex_t operator-=(const FLOATING &a) {
        re -= a;
        return *this;
    }

    complex_t operator*=(const complex_t &a) {
        FLOATING tmp = re * a.re - im * a.im;
        im = re * a.im + im * a.re;
        re = tmp;
        return *this;
    }
    complex_t operator*=(const FLOATING &a) {
        re *= a;
        im *= a;
        return *this;
    }

    friend complex_t operator-(const complex_t &a) {
        return complex_t(-a.re, -a.im);
    }

    friend bool operator==(const complex_t &a, const complex_t &b) {
        return a.re == b.re && a.im == b.im;
    }
    friend bool operator==(const FLOATING &a, const complex_t &b) {
        return a == b.re && 0. == b.im;
    }
    friend bool operator==(const complex_t &a, const FLOATING &b) {
        return a.re == b && a.im == 0.;
    }

    friend bool operator!=(const complex_t &a, const complex_t &b) {
        return a.re != b.re || a.im != b.im;
    }
    friend bool operator!=(const FLOATING &a, const complex_t &b) {
        return a != b.re || 0. != b.im;
    }
    friend bool operator!=(const complex_t &a, const FLOATING &b) {
        return a.re != b && a.im != 0.;
    }

    friend complex_t operator+(const complex_t &a, const complex_t &b) {
        return complex_t(a.re + b.re, a.im + b.im);
    }
    friend complex_t operator+(const FLOATING &a, const complex_t &b) {
        return complex_t(a + b.re, b.im);
    }
    friend complex_t operator+(const complex_t &a, const FLOATING &b) {
        return complex_t(a.re + b, a.im);
    }

    friend complex_t operator-(const complex_t &a, const complex_t &b) {
        return complex_t(a.re - b.re, a.im - b.im);
    }
    friend complex_t operator-(const FLOATING &a, const complex_t &b) {
        return complex_t(a - b.re, -b.im);
    }
    friend complex_t operator-(const complex_t &a, const FLOATING &b) {
        return complex_t(a.re - b, a.im);
    }

    friend complex_t operator*(const complex_t &a, const complex_t &b) {
        return complex_t(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
    }
    friend complex_t operator*(const FLOATING &a, const complex_t &b) {
        return complex_t(a * b.re, a * b.im);
    }
    friend complex_t operator*(const complex_t &a, const FLOATING &b) {
        return complex_t(a.re * b, a.im * b);
    }

    friend complex_t operator/(const complex_t &a, const complex_t &b) {
        FLOATING den = b.re * b.re + b.im * b.im;
        return complex_t((a.re * b.re + a.im * b.im) / den, (-a.re * b.im + a.im * b.re) / den);
    }
    friend complex_t operator/(const FLOATING &a, const complex_t &b) {
        FLOATING den = b.re * b.re + b.im * b.im;
        return complex_t(a * b.re / den, -a * b.im / den);
    }
    friend complex_t operator/(const complex_t &a, const FLOATING &b) {
        return complex_t(a.re / b, a.im / b);
    }

    friend ostream &operator<<(ostream &os, const complex_t &z) {
        os << "(" << z.re << "," << z.im << ")";
        return os;
    }

    friend complex_t conj(const complex_t &z) {
        return complex_t(z.re, -z.im);
    }
    friend FLOATING abs(const complex_t &z) {
        return sqrt(z.re * z.re + z.im * z.im);
    }
    friend FLOATING arg(const complex_t &z) {
        return atan2(z.im, z.re);
    }
};

string ftos(FLOATING x) {
    static char tmp[100];
    if (x >= 0.0) {
        snprintf(tmp, 100, _PNUMBER_, x);
    } else {
        snprintf(tmp, 100, _NUMBER_, x);
    }
    return string(tmp);
}
