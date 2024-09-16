#ifndef TAU_OO_UTILITIES_H_
#define TAU_OO_UTILITIES_H_

#include <complex>
#include <vector>

#include "Math/Vector4D.h"

typedef std::complex<double> cd_t;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> dlv_t;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<cd_t>> clv_t;

cd_t operator*(const int scale, const cd_t &v) { return double(scale) * v; }
cd_t operator*(const cd_t &v, const int scale) { return double(scale) * v; }

clv_t ToComplex(const dlv_t &v) {
    cd_t x = v.x();
    cd_t y = v.y();
    cd_t z = v.z();
    cd_t t = v.t();
    return clv_t(x, y, z, t);
}

dlv_t real(const clv_t &v) {
    double x = real(v.x());
    double y = real(v.y());
    double z = real(v.z());
    double t = real(v.t());
    return dlv_t(x, y, z, t);
}
dlv_t imag(const clv_t &v) {
    double x = imag(v.x());
    double y = imag(v.y());
    double z = imag(v.z());
    double t = imag(v.t());
    return dlv_t(x, y, z, t);
}

clv_t TrimImaginary(const clv_t &v) {
    // * Trim Imaginary part to zero
    // * Usually used when we know the imaginary should be zero
    cd_t x = real(v.x());
    cd_t y = real(v.y());
    cd_t z = real(v.z());
    cd_t t = real(v.t());
    return clv_t(x, y, z, t);
}

clv_t conjugate(const clv_t &v) {
    cd_t x = conj(v.x());
    cd_t y = conj(v.y());
    cd_t z = conj(v.z());
    cd_t t = conj(v.t());
    return clv_t(x, y, z, t);
}

clv_t epsilon(const clv_t &b, const clv_t &c, const clv_t &d) {
    // * g_{mu tau}epsilon^{tau nu rho sigma}b_nu c_rho d_sigma
    // * epsilon^{0123}=1
    cd_t t = b.x() * c.y() * d.z();  // * 0 1 2 3 +
    t -= b.x() * c.z() * d.y();      // * 0 1 3 2 -
    t -= b.y() * c.x() * d.z();      // * 0 2 1 3 -
    t -= b.y() * c.z() * d.x();      // * 0 2 3 1 +
    t += b.z() * c.x() * d.y();      // * 0 3 1 2 +
    t -= b.z() * c.y() * d.x();      // * 0 3 2 1 -

    cd_t x = -1.0 * b.t() * c.y() * d.z();  // * 1 0 2 3 -
    x += b.t() * c.z() * d.y();             // * 1 0 3 2 +
    x += b.y() * c.t() * d.z();             // * 1 2 0 3 +
    x -= b.y() * c.z() * d.t();             // * 1 2 3 0 -
    x -= b.z() * c.t() * d.y();             // * 1 3 0 2 -
    x += b.z() * c.y() * d.t();             // * 1 3 2 0 +

    cd_t y = b.t() * c.x() * d.z();  // * 2 0 1 3 +
    y -= b.t() * c.z() * d.x();      // * 2 0 3 1 -
    y -= b.x() * c.t() * d.z();      // * 2 1 0 3 -
    y += b.x() * c.z() * d.t();      // * 2 1 3 0 +
    y += b.z() * c.t() * d.x();      // * 2 3 0 1 +
    y -= b.z() * c.x() * d.t();      // * 2 3 1 0 -

    cd_t z = -b.t() * c.x() * d.y();  // * 3 0 1 2 -
    z += b.t() * c.y() * d.x();       // * 3 0 2 1 +
    z += b.x() * c.t() * d.y();       // * 3 1 0 2 +
    z -= b.x() * c.y() * d.t();       // * 3 1 2 0 -
    z -= b.y() * c.t() * d.x();       // * 3 2 0 1 -
    z += b.y() * c.x() * d.t();       // * 3 2 1 0 +

    return clv_t(t, -x, -y, -z);  // * Extra g_{mu nu}
}
cd_t epsilon(const clv_t &a, const clv_t &b, const clv_t &c, const clv_t &d) {
    // * epsilon^{mu nu rho sigma}a_mu b_nu c_rho d_sigma
    // * epsilon^{0123}=1
    // * 0 for t, 1 for x, 2 for y, 3 for z
    cd_t res = a.t() * b.x() * c.y() * d.z();  // * 0 1 2 3 +
    res -= a.t() * b.x() * c.z() * d.y();      // * 0 1 3 2 -
    res -= a.t() * b.y() * c.x() * d.z();      // * 0 2 1 3 -
    res += a.t() * b.y() * c.z() * d.x();      // * 0 2 3 1 +
    res += a.t() * b.z() * c.x() * d.y();      // * 0 3 1 2 +
    res -= a.t() * b.z() * c.y() * d.x();      // * 0 3 2 1 -
    res -= a.x() * b.t() * c.y() * d.z();      // * 1 0 2 3 -
    res += a.x() * b.t() * c.z() * d.y();      // * 1 0 3 2 +
    res += a.x() * b.y() * c.t() * d.z();      // * 1 2 0 3 +
    res -= a.x() * b.y() * c.z() * d.t();      // * 1 2 3 0 -
    res -= a.x() * b.z() * c.t() * d.y();      // * 1 3 0 2 -
    res += a.x() * b.z() * c.y() * d.t();      // * 1 3 2 0 +
    res += a.y() * b.t() * c.x() * d.z();      // * 2 0 1 3 +
    res -= a.y() * b.t() * c.z() * d.x();      // * 2 0 3 1 -
    res -= a.y() * b.x() * c.t() * d.z();      // * 2 1 0 3 -
    res += a.y() * b.x() * c.z() * d.t();      // * 2 1 3 0 +
    res += a.y() * b.z() * c.t() * d.x();      // * 2 3 0 1 +
    res -= a.y() * b.z() * c.x() * d.t();      // * 2 3 1 0 -
    res -= a.z() * b.t() * c.x() * d.y();      // * 3 0 1 2 -
    res += a.z() * b.t() * c.y() * d.x();      // * 3 0 2 1 +
    res += a.z() * b.x() * c.t() * d.y();      // * 3 1 0 2 +
    res -= a.z() * b.x() * c.y() * d.t();      // * 3 1 2 0 -
    res -= a.z() * b.y() * c.t() * d.x();      // * 3 2 0 1 -
    res += a.z() * b.y() * c.x() * d.t();      // * 3 2 1 0 +
    return res;
}

typedef std::vector<double> vd_t;
typedef std::vector<std::complex<double>> vcd_t;

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
std::complex<T> operator*(const std::vector<std::complex<T>> &left, const std::vector<T> &right) {
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
std::complex<T> operator*(const std::vector<T> &left, const std::vector<std::complex<T>> &right) {
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
std::vector<std::complex<T>> operator*(const T &scale, const std::vector<std::complex<T>> &right) {
    std::vector<std::complex<T>> res = right;
    for (size_t i = 0; i < res.size(); i++) {
        res[i] *= scale;
    }
    return res;
}

template <typename T>
std::vector<std::complex<T>> operator*(const std::vector<std::complex<T>> &left, const T &scale) {
    std::vector<std::complex<T>> res = left;
    for (size_t i = 0; i < res.size(); i++) {
        res[i] *= scale;
    }
    return res;
}

template <typename T>
std::vector<std::complex<T>> operator*(const std::complex<T> &scale, const std::vector<T> &right) {
    std::vector<std::complex<T>> res = tocomplex(right);
    for (size_t i = 0; i < res.size(); i++) {
        res[i] *= scale;
    }
    return res;
}

template <typename T>
std::vector<std::complex<T>> operator*(const std::vector<T> &left, const std::complex<T> &scale) {
    std::vector<std::complex<T>> res = tocomplex(left);
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
std::vector<std::complex<T>> epsilon(const std::vector<std::complex<T>> &a, const std::vector<std::complex<T>> &b,
                                     const std::vector<T> &c) {
    // * epsilon^{mu,nu,rho,sigma}a_{nu}b_{rho}c_{sigma}
    assert(a.size() == 4);
    assert(b.size() == 4);
    assert(c.size() == 4);
    std::vector<std::complex<T>> res(4);

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
