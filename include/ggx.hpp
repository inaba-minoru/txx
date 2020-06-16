#ifndef GGX_H
#define GGX_H

#include <assert.h>
#include <vecmath.h>

// i 入射
// o 出射
// n 宏观法向量
// m 微观法向量

// ok
double F(const Vector3f& i, const Vector3f& m, const double& eta_t,
        const double& eta_i) {
    double c = fabs(Vector3f::dot(i, m));
    double g2 = eta_t * eta_t / (eta_i * eta_i) - 1 + c * c;
    if (g2 < 0) return 1;  // if g is imaginary, total reflection
    double g = sqrtf(g2);
    double ret = pow(g - c, 2) *
                (1 + pow(c * (g + c) - 1, 2) / pow(c * (g - c) + 1, 2)) /
                (2 * pow(g + c, 2));
    assert(ret >= 0 && ret <= 1);
    return ret;
}

inline double ksip(const double& x) { return x > 0 ? 1 : 0; }
inline double sign(const double& x) { return x >= 0 ? 1 : -1; }

double D(const Vector3f& m, const Vector3f& n, const double& alpha_g2) {
    double mn = Vector3f::dot(m, n);
    double ret =
        alpha_g2 * ksip(mn) / (M_PI * pow(mn * mn * (alpha_g2 - 1) + 1, 2));
    return ret;
}

double G1(const Vector3f& v, const Vector3f& m, const Vector3f& n,
         const double& alpha_g2) {
    // double vm = Vector3f::dot(v, m);
    // double vn = Vector3f::dot(v, n);
    // double ret =
    //     ksip(vm * vn) * 2 / (1 + sqrt(1 + alpha_g2 * (1 / (vn * vn) - 1)));
    // // if (!(ret >= 0 && ret <= 1)) printf("%f %f %f\n", vm, vn, ret);
    // return ret;

    double vm = Vector3f::dot(v, m);
    double vn = Vector3f::dot(v, n);
    double ret =
        ksip(vm * vn) * std::min(1., 2 * fabs(Vector3f::dot(m, n) * vn / vm));
    return ret;
}

double G(const Vector3f& i, const Vector3f& o, const Vector3f& m,
        const Vector3f& n, const double& alpha_g2) {
    // double ret = G1(i, m, n, alpha_g2) * G1(o, m, n, alpha_g2);
    // // assert(ret >= 0 && ret <= 1);
    // // if (!(ret >= 0 && ret <= 1)) printf("%f\n", ret);
    // return ret;

    double ret;
    if (Vector3f::dot(o, n) > 0)
        ret = std::min(G1(i, m, n, alpha_g2), G1(o, m, n, alpha_g2));
    else
        ret = std::max(G1(i, m, n, alpha_g2) + G1(o, m, n, alpha_g2) - 1, 0.);
    return ret;
}

double fr(const Vector3f& i, const Vector3f& o, const Vector3f& n,
         const double& eta_o, const double& eta_i, const double& alpha_g2,
         const double& in, const double& on) {
    if (i == -o) return 1;

    Vector3f hr = (i + o).normalized();
    if (Vector3f::dot(i, hr) < 0) hr = -hr;

    return F(i, hr, eta_o, eta_i) * G(i, o, hr, n, alpha_g2) *
           D(hr, n, alpha_g2) / (4 * in * on);
}

double ft(const Vector3f& i, const Vector3f& o, const Vector3f& n,
         const double& eta_o, const double& eta_i, const double& alpha_g2,
         const double& in, const double& on) {
    // Vector3f ht = -(eta_i * i + eta_o * o).normalized();
    Vector3f ht = (eta_i * i + eta_o * o).normalized();
    if (Vector3f::dot(i, ht) < 0) ht = -ht;
    assert(Vector3f::dot(i, ht) >= 0);
    // printf("%f\n", Vector3f::dot(i, ht));

    double iht = Vector3f::dot(i, ht);
    double oht = Vector3f::dot(o, ht);
    return iht * oht * eta_o * eta_o * (1 - F(i, ht, eta_o, eta_i)) *
           G(i, o, ht, n, alpha_g2) * D(ht, n, alpha_g2) /
           (in * on * pow(eta_i * iht + eta_o * oht, 2));
}

double fs(const Vector3f& i, const Vector3f& o, const Vector3f& n,
         const double& eta_o, const double& eta_i, const double& alpha_g2) {
    double in = Vector3f::dot(i, n);
    double on = Vector3f::dot(o, n);
    return fr(i, o, n, eta_o, eta_i, alpha_g2, in, on) +
           ft(i, o, n, eta_o, eta_i, alpha_g2, in, on);
}

void sample(unsigned short* Xi, const Vector3f& i, const Vector3f& n,
            Vector3f& o, double& weight, double& new_eta, const double& eta_t,
            const double& eta_i, const double& alpha_g, const double& alpha_g2) {
    double epsilon = erand48(Xi);
    double theta_m = atanf(alpha_g * sqrt(epsilon) / sqrt(1 - epsilon));
    double phi_m = 2 * M_PI * erand48(Xi);

    Vector3f alpha =
        (fabs(n.x()) < 0.7f ? Vector3f::cross(n, Vector3f(1, 0, 0))
                             : Vector3f::cross(n, Vector3f(0, 1, 0)))
            .normalized();
    Vector3f beta = Vector3f::cross(n, alpha).normalized();

    Vector3f m = (cos(theta_m) * n +
                  sin(theta_m) * (cos(phi_m) * alpha + sin(phi_m) * beta))
                     .normalized();
    double mn = Vector3f::dot(m, n);
    Vector3f m_ = 2 * mn * n - m;
    double im = Vector3f::dot(i, m);
    double im_ = Vector3f::dot(i, m_);

    if (erand48(Xi) * (im + im_) < im_) {
        m = m_;
        im = im_;
        mn = Vector3f::dot(m, n);
    }

    assert(im >= 0);

    // n可以是任意的法向量，得到的o, weight都正确
    double in = Vector3f::dot(i, n);

    if (erand48(Xi) < F(i, m, eta_t, eta_i)) {
        // o = (2 * fabs(im) * m - i).normalized();
        o = (2 * im * m - i).normalized();
        new_eta = eta_i;
    } else {
        double eta = eta_i / eta_t;
        // o = (eta * im - sign(in) * sqrt(1 + eta * (im * im - 1))) * m - eta *
        // i;
        // 原先少了一个eta?
        o = ((eta * im - sign(in) * sqrt(1 + eta * eta * (im * im - 1))) * m -
             eta * i)
                .normalized();
        new_eta = eta_t;
    }
    // weight = fabs(im) * G(i, o, m, n, alpha_g2) / (fabs(in) * fabs(mn));
    weight = G(i, o, m, n, alpha_g2) / G1(i, m, n, alpha_g2);
}

#endif  // GGX_H