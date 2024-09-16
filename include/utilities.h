#ifndef TAU_OO_UTILITIES_H_
#define TAU_OO_UTILITIES_H_

#include <complex>
#include <vector>

typedef std::vector<double> vd_t;
typedef std::vector<std::complex<double> > vcd_t;

vd_t real(const vcd_t &v) {
    vd_t res(v.size());
    for (size_t i = 0; i < v.size(); i++) {
        res[i] = real(v[i]);
    }
    return res;
}
vd_t imag(const vcd_t &v) {
    vd_t res(v.size());
    for (size_t i = 0; i < v.size(); i++) {
        res[i] = imag(v[i]);
    }
    return res;
}
vcd_t conjugate(const vcd_t &v) {
    vcd_t res(v.size());
    for (size_t i = 0; i < res.size(); i++) {
        res[i] = conj(v[i]);
    }
    return res;
}

vcd_t tocomplex(const vd_t &v) {
    vcd_t res(v.size());
    for (size_t i = 0; i < res.size(); i++) {
        res[i] = v[i];
    }
    return res;
}

template <typename T>
std::complex<T> operator*(const std::vector<std::complex<T> > &left, const std::vector<T> &right) {
    assert(left.size() == right.size());
    assert(left.size() == 4);
    std::complex<T> res = 0;
    T sign[4] = {1, -1, -1, -1};
    for (size_t i = 0; i < left.size(); i++) {
        res += sign[i] * left[i] * right[i];
    }
    return res;
}

template <typename T>
std::complex<T> operator*(const std::vector<T> &left, const std::vector<std::complex<T> > &right) {
    assert(left.size() == right.size());
    assert(left.size() == 4);
    std::complex<T> res = 0;
    T sign[4] = {1, -1, -1, -1};
    for (size_t i = 0; i < left.size(); i++) {
        res += sign[i] * left[i] * right[i];
    }
    return res;
}

template <typename T>
T operator*(const std::vector<T> &left, const std::vector<T> &right) {
    assert(left.size() == right.size());
    assert(left.size() == 4);
    T sign[4] = {1, -1, -1, -1};
    T res = 0;
    for (size_t i = 0; i < left.size(); i++) {
        res += sign[i] * left[i] * right[i];
    }
    return res;
}

template <typename T>
std::vector<T> operator+(const std::vector<T> &left, const std::vector<T> &right) {
    assert(left.size() == right.size());
    std::vector<T> res = left;
    for (size_t i = 0; i < res.size(); i++) {
        res[i] += right[i];
    }
    return res;
}

template <typename T>
std::vector<T> operator-(const std::vector<T> &left, const std::vector<T> &right) {
    assert(left.size() == right.size());
    std::vector<T> res = left;
    for (size_t i = 0; i < res.size(); i++) {
        res[i] -= right[i];
    }
    return res;
}

template <typename T>
std::vector<T> operator*(const T &scale, const std::vector<T> &right) {
    std::vector<T> res = right;
    for (size_t i = 0; i < res.size(); i++) {
        res[i] *= scale;
    }
    return res;
}

template <typename T>
std::vector<T> operator*(const std::vector<T> &left, const T &scale) {
    std::vector<T> res = left;
    for (size_t i = 0; i < res.size(); i++) {
        res[i] *= scale;
    }
    return res;
}

template <typename T>
std::vector<std::complex<T> > operator*(const T &scale, const std::vector<std::complex<T> > &right) {
    std::vector<std::complex<T> > res = right;
    for (size_t i = 0; i < res.size(); i++) {
        res[i] *= scale;
    }
    return res;
}

template <typename T>
std::vector<std::complex<T> > operator*(const std::vector<std::complex<T> > &left, const T &scale) {
    std::vector<std::complex<T> > res = left;
    for (size_t i = 0; i < res.size(); i++) {
        res[i] *= scale;
    }
    return res;
}

template <typename T>
std::vector<std::complex<T> > operator*(const std::complex<T> &scale, const std::vector<T> &right) {
    std::vector<std::complex<T> > res = tocomplex(right);
    for (size_t i = 0; i < res.size(); i++) {
        res[i] *= scale;
    }
    return res;
}

template <typename T>
std::vector<std::complex<T> > operator*(const std::vector<T> &left, const std::complex<T> &scale) {
    std::vector<std::complex<T> > res = tocomplex(left);
    for (size_t i = 0; i < res.size(); i++) {
        res[i] *= scale;
    }
    return res;
}

template <typename T>
std::vector<T> operator/(const std::vector<T> &left, const T &fac) {
    std::vector<T> res = left;
    for (size_t i = 0; i < res.size(); i++) {
        res[i] /= fac;
    }
    return res;
}

template <typename T>
std::vector<T> epsilon(const std::vector<T> &a, const std::vector<T> &b, const std::vector<T> &c) {
    // * epsilon^{mu,nu,rho,sigma}a_{nu}b_{rho}c_{sigma}
    assert(a.size() == 4);
    assert(b.size() == 4);
    assert(c.size() == 4);
    std::vector<T> res(4);

    // *  0 1 2 3 - 0 1 3 2 + 0 2 3 1 - 0 2 1 3 + 0 3 1 2 - 0 3 2 1;
    res[0] = a[1] * b[2] * c[3] - a[1] * b[3] * c[2] + a[2] * b[3] * c[1] - a[2] * b[1] * c[3] + a[3] * b[1] * c[2] -
             a[3] * b[2] * c[1];

    // * 1 0 3 2 - 1 0 2 3 + 1 2 0 3 - 1 2 3 0 + 1 3 2 0 - 1 3 0 2;
    res[1] = -1.0 * (a[0] * b[3] * c[2] - a[0] * b[2] * c[3] + a[2] * b[0] * c[3] - a[2] * b[3] * c[0] +
                     a[3] * b[2] * c[0] - a[3] * b[0] * c[2]);

    // * 2 0 1 3 - 2 0 3 1 + 2 1 3 0 - 2 1 0 3 + 2 3 1 0 - 2 3 0 1;
    res[2] = -1.0 * (a[0] * b[1] * c[3] - a[0] * b[3] * c[1] + a[1] * b[3] * c[0] - a[1] * b[0] * c[3] +
                     a[3] * b[1] * c[0] - a[3] * b[0] * c[1]);

    // * 3 0 2 1 - 3 0 1 2 + 3 1 0 2 - 3 1 2 0 + 3 2 1 0 - 3 2 0 1;
    res[3] = -1.0 * (a[0] * b[2] * c[1] - a[0] * b[1] * c[2] + a[1] * b[0] * c[2] - a[1] * b[2] * c[0] +
                     a[2] * b[1] * c[0] - a[2] * b[0] * c[1]);

    return res;
}

template <typename T>
std::vector<std::complex<T> > epsilon(const std::vector<std::complex<T> > &a, const std::vector<std::complex<T> > &b,
                                      const std::vector<T> &c) {
    // * epsilon^{mu,nu,rho,sigma}a_{nu}b_{rho}c_{sigma}
    assert(a.size() == 4);
    assert(b.size() == 4);
    assert(c.size() == 4);
    std::vector<std::complex<T> > res(4);

    // *  0 1 2 3 - 0 1 3 2 + 0 2 3 1 - 0 2 1 3 + 0 3 1 2 - 0 3 2 1;
    res[0] = a[1] * b[2] * c[3] - a[1] * b[3] * c[2] + a[2] * b[3] * c[1] - a[2] * b[1] * c[3] + a[3] * b[1] * c[2] -
             a[3] * b[2] * c[1];

    // * 1 0 3 2 - 1 0 2 3 + 1 2 0 3 - 1 2 3 0 + 1 3 2 0 - 1 3 0 2;
    res[1] = -1.0 * (a[0] * b[3] * c[2] - a[0] * b[2] * c[3] + a[2] * b[0] * c[3] - a[2] * b[3] * c[0] +
                     a[3] * b[2] * c[0] - a[3] * b[0] * c[2]);

    // * 2 0 1 3 - 2 0 3 1 + 2 1 3 0 - 2 1 0 3 + 2 3 1 0 - 2 3 0 1;
    res[2] = -1.0 * (a[0] * b[1] * c[3] - a[0] * b[3] * c[1] + a[1] * b[3] * c[0] - a[1] * b[0] * c[3] +
                     a[3] * b[1] * c[0] - a[3] * b[0] * c[1]);

    // * 3 0 2 1 - 3 0 1 2 + 3 1 0 2 - 3 1 2 0 + 3 2 1 0 - 3 2 0 1;
    res[3] = -1.0 * (a[0] * b[2] * c[1] - a[0] * b[1] * c[2] + a[1] * b[0] * c[2] - a[1] * b[2] * c[0] +
                     a[2] * b[1] * c[0] - a[2] * b[0] * c[1]);

    return res;
}

template <typename T>
T epsilon(const std::vector<T> &a, const std::vector<T> &b, const std::vector<T> &c, const std::vector<T> &d) {
    // * epsilon^{mu,nu,rho,sigma}a_{mu}b_{nu}c_{rho}d_{sigma}
    assert(a.size() == 4);
    assert(b.size() == 4);
    assert(c.size() == 4);
    assert(d.size() == 4);
    return a[0] * b[1] * c[2] * d[3] + a[0] * b[2] * c[3] * d[1] + a[0] * b[3] * c[1] * d[2] +
           a[1] * b[0] * c[3] * d[2] + a[1] * b[2] * c[0] * d[3] + a[1] * b[3] * c[2] * d[0] +
           a[2] * b[0] * c[1] * d[3] + a[2] * b[1] * c[3] * d[0] + a[2] * b[3] * c[0] * d[1] +
           a[3] * b[0] * c[2] * d[1] + a[3] * b[2] * c[1] * d[0] + a[3] * b[1] * c[0] * d[2] -
           a[0] * b[1] * c[3] * d[2] - a[0] * b[2] * c[1] * d[3] - a[0] * b[3] * c[2] * d[1] -
           a[1] * b[0] * c[2] * d[3] - a[1] * b[2] * c[3] * d[0] - a[1] * b[3] * c[0] * d[2] -
           a[2] * b[0] * c[3] * d[1] - a[2] * b[1] * c[0] * d[3] - a[2] * b[3] * c[1] * d[0] -
           a[3] * b[0] * c[1] * d[2] - a[3] * b[2] * c[0] * d[1] - a[3] * b[1] * c[2] * d[0];
    // return a * epsilon(b, c, d);
}

#endif  // TAU_OO_UTILITIES_H_
