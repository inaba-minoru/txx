#ifndef GGX_H
#define GGX_H

#include <assert.h>
#include <vecmath.h>

// i 入射
// o 出射
// n 宏观法向量
// m 微观法向量

inline float sign(const float& x) { return x >= 0 ? 1 : -1; }

float D(const Vector3f& m, const Vector3f& n, const float& alpha_g2) {
    float mn = Vector3f::dot(m, n);
    if (mn > 0) {
        return alpha_g2 / (M_PI * pow(mn * mn * (alpha_g2 - 1) + 1, 2));
    }
    return 0;
}

float G1(const Vector3f& v, const Vector3f& m, const Vector3f& n,
         const float& alpha_g2) {
    float vn = std::max(-0.999999f, std::min(0.999999f, Vector3f::dot(v, n)));
    float vn2 = vn * vn;
    if (vn * Vector3f::dot(v, m) > 0) {
        float ret = 2 / (1 + sqrt(1 + alpha_g2 * vn2 / (1 - vn2)));
        assert(ret >= 0 && ret <= 1);
        return ret;
    }
    return 0;
}

float G(const Vector3f& i, const Vector3f& o, const Vector3f& m,
        const Vector3f& n, const float& alpha_g2) {
    float ret = G1(i, m, n, alpha_g2) * G1(o, m, n, alpha_g2);
    assert(ret >= 0 && ret <= 1);
    return ret;
}

float F(const Vector3f& i, const Vector3f& m, const float& eta_t,
        const float& eta_i) {
    float c = Vector3f::dot(i, m);

    // assert(c >= 0);

    float g2 = eta_t * eta_t / (eta_i * eta_i) - 1 + c * c;
    if (g2 < 0) return 1;

    float g = sqrt(g2);
    float ret = pow(g - c, 2) *
                (1 + pow(c * (g + c) - 1, 2) / pow(c * (g - c) + 1, 2)) /
                (2 * pow(g + c, 2));

    // assert(ret >= 0 && ret <= 1);
    if (!(ret >= 0 && ret <= 1)) {
#pragma omp critical
        {
            // printf("F %f %f %f\n", eta_t / eta_i, c, ret);
            // i.print();
            // m.print();
        }
    }

    return ret;
}

float fr(const Vector3f& i, const Vector3f& o, const Vector3f& n,
         const float& eta_o, const float& eta_i, const float& alpha_g2,
         const float& in, const float& on) {
    if (i == -o) return 1;

    Vector3f hr = (i + o).normalized();
    if (Vector3f::dot(i, hr) < 0) hr = -hr;

    return F(i, hr, eta_o, eta_i) * G(i, o, hr, n, alpha_g2) *
           D(hr, n, alpha_g2) / (4 * in * on);
}

float ft(const Vector3f& i, const Vector3f& o, const Vector3f& n,
         const float& eta_o, const float& eta_i, const float& alpha_g2,
         const float& in, const float& on) {
    // Vector3f ht = -(eta_i * i + eta_o * o).normalized();
    Vector3f ht = (eta_i * i + eta_o * o).normalized();
    if (Vector3f::dot(i, ht) < 0) ht = -ht;
    assert(Vector3f::dot(i, ht) >= 0);
    // printf("%f\n", Vector3f::dot(i, ht));

    float iht = Vector3f::dot(i, ht);
    float oht = Vector3f::dot(o, ht);
    return iht * oht * eta_o * eta_o * (1 - F(i, ht, eta_o, eta_i)) *
           G(i, o, ht, n, alpha_g2) * D(ht, n, alpha_g2) /
           (in * on * pow(eta_i * iht + eta_o * oht, 2));
}

float fs(const Vector3f& i, const Vector3f& o, const Vector3f& n,
         const float& eta_o, const float& eta_i, const float& alpha_g2) {
    float in = Vector3f::dot(i, n);
    float on = Vector3f::dot(o, n);
    return fr(i, o, n, eta_o, eta_i, alpha_g2, in, on) +
           ft(i, o, n, eta_o, eta_i, alpha_g2, in, on);
}

void sample(unsigned short* Xi, const Vector3f& i, const Vector3f& n,
            Vector3f& o, float& weight, float& new_eta, const float& eta_t,
            const float& eta_i, const float& alpha_g, const float& alpha_g2) {
    float epsilon = erand48(Xi);
    float theta_m = atanf(alpha_g * sqrt(epsilon) / sqrt(1 - epsilon));
    float phi_m = 2 * M_PI * erand48(Xi);

    Vector3f alpha =
        (fabsf(n.x()) < 0.7f ? Vector3f::cross(n, Vector3f(1, 0, 0))
                             : Vector3f::cross(n, Vector3f(0, 1, 0)))
            .normalized();
    Vector3f beta = Vector3f::cross(n, alpha).normalized();

    Vector3f m = (cosf(theta_m) * n +
                  sinf(theta_m) * (cosf(phi_m) * alpha + sinf(phi_m) * beta))
                     .normalized();

    // if (!(m.x() >= -1 && m.x() <= 1)) {
    //     printf("m %f %f\n", theta_m, phi_m);
    //     m.print();
    // }

    // n可以是任意的法向量，得到的o, weight都正确
    float im = Vector3f::dot(i, m);
    if (im < 0) m = -m, im = -im;
    float in = Vector3f::dot(i, n);
    float mn = Vector3f::dot(m, n);
    // assert(im >= 0);
    // printf("%f\n", F(i, m, eta_t, eta_i));

    //     if (n.y() < -0.9) {
    // #pragma omp critical
    //         {
    //             // printf("%f\n", F(i, m, eta_t, eta_i));
    //             puts("haha");
    //             i.print();
    //             m.print();
    //         }
    //     }
    //     if (n.z() > 0.9) {
    // #pragma omp critical
    //         {
    //             // printf("%f\n", F(i, m, eta_t, eta_i));
    //             puts("haha");
    //             n.print();
    //             i.print();
    //             m.print();
    //         }
    //     }

    float on = Vector3f::dot(o, n);
    float k;

    if (erand48(Xi) < F(i, m, eta_t, eta_i)) {
        o = 2 * im * m - i;
        new_eta = eta_i;

        k = fr(o, i, n, eta_t, eta_i, alpha_g2, in, on);
    } else {
        float eta = eta_i / eta_t;
        // o = (eta * im - sign(in) * sqrt(1 + eta * (im * im - 1))) * m - eta *
        // i;
        // 原先少了一个eta?
        o = (eta * im - sign(in) * sqrt(1 + eta * eta * (im * im - 1))) * m -
            eta * i;
        new_eta = eta_t;

        k = ft(o, i, n, eta_t, eta_i, alpha_g2, in, on);
    }
    weight = im * G(i, o, m, n, alpha_g2) / (in * mn);

    weight = std::min(weight, 1.f);
}

#endif  // GGX_H