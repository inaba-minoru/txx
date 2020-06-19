#ifndef MATERIAL_H
#define MATERIAL_H

#include <vecmath.h>

#include <cassert>
#include <iostream>
#include <texture.hpp>

#include "hit.hpp"
#include "ray.hpp"

// TODO: Implement Shade function that computes Phong introduced in class.
class Material {
   public:
    double alpha_g, alpha_g2;

    explicit Material(const Vector3f &d_color, const Vector3f &s_color,
                      const double &s, const double &eta, const double &alpha_g,
                      const Vector3f &emission)
        : diffuseColor(d_color),
          specularColor(s_color),
          shininess(s),
          eta(eta),
          alpha_g(alpha_g),
          emission(emission) {
        alpha_g2 = alpha_g * alpha_g;
    }

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const { return diffuseColor; }

    Vector3f Shade(const Ray &ray, const Hit &hit, const Vector3f &dirToLight,
                   const Vector3f &lightColor) {
        // Vector3f shaded = Vector3f::ZERO;

        Vector3f kdre =
            diffuseColor *
            std::max(double(0), Vector3f::dot(dirToLight.normalized(),
                                              hit.getNormal().normalized()));
        Vector3f ksre =
            specularColor *
            std::pow(std::max(double(0),
                              Vector3f::dot(
                                  -ray.getDirection().normalized(),
                                  2 *
                                          Vector3f::dot(
                                              dirToLight.normalized(),
                                              hit.getNormal().normalized()) *
                                          hit.getNormal().normalized() -
                                      dirToLight.normalized())),
                     shininess);

        return lightColor * (kdre + ksre);

        // return lightColor *
        //        (diffuseColor *
        //             std::max(double(0),
        //                      Vector3f::dot(dirToLight, hit.getNormal())) +
        //         specularColor *
        //             pow(std::max(
        //                     double(0),
        //                     Vector3f::dot(
        //                         -ray.getDirection(),
        //                         2 * Vector3f::dot(dirToLight,
        //                         hit.getNormal()) *
        //                                 hit.getNormal() -
        //                             dirToLight)),
        //                 shininess));

        // return shaded;
    }

    Texture texture;
    void loadTexture(const char *filename) { texture.load(filename); }

   public:
    // double kd = 0.0;  // 漫反射占比
    // double ks = 1.0;  // 镜面反射占比
    // double kr = 0;    // 透射占比
    // double alpha = 0.01;
    // double alpha2 = 0.0001;
    // double ior = 30;

    Vector3f diffuseColor;
    Vector3f specularColor;
    double shininess;

   public:
    Material(const Vector3f &color, const Vector3f &emission, const double &eta)
        : color(color), emission(emission), eta(eta) {}

    Vector3f color, emission;
    double eta;
    virtual void sample(unsigned short *Xi, const Vector3f &i,
                        const Vector3f &n, const double &eta_i,
                        const double &eta_o, Vector3f &o, double &w,
                        double &eta, bool &accept_light) = 0;
    virtual double f(const Vector3f &i, const Vector3f &n, const Vector3f &o,
                     const double &eta_i, const double &eta_o) = 0;
};

class SpecularMaterial : public Material {
   public:
    SpecularMaterial(const Vector3f &color, const Vector3f &emission = 0,
                     const double &eta = 0)
        : Material(color, emission, eta) {}
    ~SpecularMaterial() {}

    virtual void sample(unsigned short *Xi, const Vector3f &i,
                        const Vector3f &n, const double &eta_i,
                        const double &eta_o, Vector3f &o, double &w,
                        double &eta, bool &accept_light) {
        double in = Vector3f::dot(i, n);
        assert(in >= 0);
        o = 2 * in * n - i;
        w = 1;
        accept_light = 1;
    }

    virtual double f(const Vector3f &i, const Vector3f &n, const Vector3f &o,
                     const double &eta_i, const double &eta_o) {
        return 0;
    }
};

class DiffuseMaterial : public Material {
   public:
    DiffuseMaterial(const Vector3f &color, const Vector3f &emission = 0,
                    const double &eta = 0)
        : Material(color, emission, eta) {}
    ~DiffuseMaterial() {}

    virtual void sample(unsigned short *Xi, const Vector3f &i,
                        const Vector3f &n, const double &eta_i,
                        const double &eta_o, Vector3f &o, double &w,
                        double &eta, bool &accept_light) {
        // double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s =
        // sqrt(r2);
        double theta = asin(erand48(Xi));
        // double theta = erand48(Xi) * M_PI / 2;
        double phi = erand48(Xi) * 2 * M_PI;
        Vector3f alpha =
            (fabs(n.x()) < 0.707 ? Vector3f::cross(n, Vector3f(1, 0, 0))
                                 : Vector3f::cross(n, Vector3f(0, 1, 0)))
                .normalized();
        Vector3f beta = Vector3f::cross(alpha, n);
        o = cos(theta) * n + sin(theta) * (cos(phi) * alpha + sin(phi) * beta);
        // w = 2 * cos(theta);
        // o = (alpha * cos(r1) * r2s + beta * sin(r1) * r2s + n * sqrt(1 -
        // r2))
        //         .normalized();
        w = 1;
        accept_light = 0;
    }

    virtual double f(const Vector3f &i, const Vector3f &n, const Vector3f &o,
                     const double &eta_i, const double &eta_o) {
        return Vector3f::dot(o, n) < 0 ? 0 : Vector3f::dot(i, n) / M_PI;
    }
};

class DielectricMaterial : public Material {
   public:
    DielectricMaterial(const Vector3f &color, const Vector3f &emission = 0,
                       const double &eta = 0)
        : Material(color, emission, eta) {}
    ~DielectricMaterial() {}

    double FrDielectric(const double &cosThetaI, const double &etaI,
                        const double &etaT) {
        double sinThetaI =
            std::sqrt(std::max((double)0, 1 - cosThetaI * cosThetaI));
        double sinThetaT = etaI / etaT * sinThetaI;
        if (sinThetaT >= 1) return 1;

        double cosThetaT =
            std::sqrt(std::max((double)0, 1 - sinThetaT * sinThetaT));

        double Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                       ((etaT * cosThetaI) + (etaI * cosThetaT));
        double Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                       ((etaI * cosThetaI) + (etaT * cosThetaT));
        return (Rparl * Rparl + Rperp * Rperp) / 2;
    }

    virtual void sample(unsigned short *Xi, const Vector3f &i,
                        const Vector3f &n, const double &eta_i,
                        const double &eta_o, Vector3f &o, double &w,
                        double &eta, bool &accept_light) {
        double in = Vector3f::dot(i, n);
        double Fr = FrDielectric(in, eta_i, eta_o);
        if (erand48(Xi) < Fr) {
            o = 2 * in * n - i;
        } else {
            double eta_ido = eta_i / eta_o;
            o = ((eta_ido * in - sqrt(1 + eta_ido * eta_ido * (in * in - 1))) *
                     n -
                 eta_ido * i)
                    .normalized();
            eta = eta_o;
        }
        w = 1;
        accept_light = 1;
    }

    virtual double f(const Vector3f &i, const Vector3f &n, const Vector3f &o,
                     const double &eta_i, const double &eta_o) {
        return 0;
    }
};

class GGXMaterial : public Material {
   public:
    double alpha_g, alpha_g2;

    GGXMaterial(const Vector3f &color, const Vector3f &emission,
                const double &eta, const double &alpha_g)
        : Material(color, emission, eta) {
        alpha_g2 = alpha_g * alpha_g;
    }

    virtual void sample(unsigned short *Xi, const Vector3f &i,
                        const Vector3f &n, const double &eta_i,
                        const double &eta_o, Vector3f &o, double &w,
                        double &eta, bool &accept_light) {
        sample(Xi, i, n, o, w, eta, eta_o, eta_i, alpha_g, alpha_g2);
        accept_light = 0;
    }
    virtual double f(const Vector3f &i, const Vector3f &n, const Vector3f &o,
                     const double &eta_i, const double &eta_o) {
        return fs(i, o, n, eta_o, eta_i, alpha_g2);
    }

    double F(const Vector3f &i, const Vector3f &m, const double &eta_t,
             const double &eta_i) {
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

    inline double ksip(const double &x) { return x > 0 ? 1 : 0; }
    inline double sign(const double &x) { return x >= 0 ? 1 : -1; }

    double D(const Vector3f &m, const Vector3f &n, const double &alpha_g2) {
        double mn = Vector3f::dot(m, n);
        double ret =
            alpha_g2 * ksip(mn) / (M_PI * pow(mn * mn * (alpha_g2 - 1) + 1, 2));
        return ret;
    }

    double G1(const Vector3f &v, const Vector3f &m, const Vector3f &n,
              const double &alpha_g2) {
        // double vm = Vector3f::dot(v, m);
        // double vn = Vector3f::dot(v, n);
        // double ret =
        //     ksip(vm * vn) * 2 / (1 + sqrt(1 + alpha_g2 * (1 / (vn * vn) -
        //     1)));
        // // if (!(ret >= 0 && ret <= 1)) printf("%f %f %f\n", vm, vn, ret);
        // return ret;

        double vm = Vector3f::dot(v, m);
        double vn = Vector3f::dot(v, n);
        double ret = ksip(vm * vn) *
                     std::min(1., 2 * fabs(Vector3f::dot(m, n) * vn / vm));
        return ret;
    }

    double G(const Vector3f &i, const Vector3f &o, const Vector3f &m,
             const Vector3f &n, const double &alpha_g2) {
        // double ret = G1(i, m, n, alpha_g2) * G1(o, m, n, alpha_g2);
        // // assert(ret >= 0 && ret <= 1);
        // // if (!(ret >= 0 && ret <= 1)) printf("%f\n", ret);
        // return ret;

        double ret;
        if (Vector3f::dot(o, n) > 0)
            ret = std::min(G1(i, m, n, alpha_g2), G1(o, m, n, alpha_g2));
        else
            ret =
                std::max(G1(i, m, n, alpha_g2) + G1(o, m, n, alpha_g2) - 1, 0.);
        return ret;
    }

    double fr(const Vector3f &i, const Vector3f &o, const Vector3f &n,
              const double &eta_o, const double &eta_i, const double &alpha_g2,
              const double &in, const double &on) {
        return 0;

        Vector3f hr = (i + o).normalized();
        // if (Vector3f::dot(i, hr) < 0) hr = -hr;

        return F(i, hr, eta_o, eta_i) * G(i, o, hr, n, alpha_g2) *
               D(hr, n, alpha_g2) / (4 * in * on);
    }

    double ft(const Vector3f &i, const Vector3f &o, const Vector3f &n,
              const double &eta_o, const double &eta_i, const double &alpha_g2,
              const double &in, const double &on) {
        return 0;

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

    double fs(const Vector3f &i, const Vector3f &o, const Vector3f &n,
              const double &eta_o, const double &eta_i,
              const double &alpha_g2) {
        double in = Vector3f::dot(i, n);
        double on = Vector3f::dot(o, n);
        return fr(i, o, n, eta_o, eta_i, alpha_g2, in, on) +
               ft(i, o, n, eta_o, eta_i, alpha_g2, in, on);
    }

    void sample(unsigned short *Xi, const Vector3f &i, const Vector3f &n,
                Vector3f &o, double &weight, double &new_eta,
                const double &eta_t, const double &eta_i, const double &alpha_g,
                const double &alpha_g2) {
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
            // o = (eta * im - sign(in) * sqrt(1 + eta * (im * im - 1))) * m -
            // eta * i; 原先少了一个eta?
            o = ((eta * im - sign(in) * sqrt(1 + eta * eta * (im * im - 1))) *
                     m -
                 eta * i)
                    .normalized();
            new_eta = eta_t;
        }
        // weight = fabs(im) * G(i, o, m, n, alpha_g2) / (fabs(in) * fabs(mn));
        weight = G(i, o, m, n, alpha_g2) / G1(i, m, n, alpha_g2);
    }
};

#endif  // MATERIAL_H
