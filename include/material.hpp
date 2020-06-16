#ifndef MATERIAL_H
#define MATERIAL_H

#include <vecmath.h>

#include <cassert>
#include <iostream>

#include "hit.hpp"
#include "ray.hpp"

// TODO: Implement Shade function that computes Phong introduced in class.
class Material {
   public:
    double eta, alpha_g, alpha_g2;
    Vector3f emission;

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
        // printf("%f\n", alpha_g);
        // diffuseColor.print();
        emission.print();
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

    // Vector3f CookTorrance(const Ray &ray, const Hit &hit, const Vector3f &l,
    //                       const Vector3f &light_color) {
    //     Vector3f n = hit.getNormal();
    //     Vector3f v = -ray.getDirection();
    //     Vector3f h = (l + v).normalized();

    //     double ln = Vector3f::dot(l, n);
    //     double s = 0.0;
    //     if (ln > 0.0) {
    //         double hn = Vector3f::dot(h, n);
    //         double vn = Vector3f::dot(v, n);
    //         double hv = Vector3f::dot(h, v);

    //         // G项
    //         double hn2 = 2.0 * hn;
    //         double g1 = (hn2 * ln) / hv;
    //         double g2 = (hn2 * vn) / hv;
    //         double g = std::min((double)1.0, std::min(g1, g2));

    //         // D项：beckmann distribution function
    //         double m2 = roughness * roughness;
    //         double r1 = 1.0 / (4.0 * m2 * pow(hn, 4.0));
    //         double r2 = (hn * hn - 1.0) / (m2 * hn * hn);
    //         double d = r1 * exp(r2);

    //         // F项
    //         double f = pow(1.0 - hv, 5.0);
    //         f *= (1.0 - fresnel);
    //         f += fresnel;

    //         s = (f * g * d) / (vn * ln * M_PI);

    //         return light_color * ln *
    //                (k * diffuseColor * 2 + s * (1 - k) * specularColor * 2);
    //     }
    //     return Vector3f::ZERO;
    // }

    Vector3f getfd() { return diffuseColor / M_PI; }
    Vector3f getfs(const Vector3f &v, const Vector3f &n, const Vector3f &l) {
        Vector3f h = (v + l).normalized();

        double ln = Vector3f::dot(l, n);
        double vn = Vector3f::dot(v, n);
        double hn = Vector3f::dot(h, n);
        double hv = Vector3f::dot(h, v);

        // D
        double temp = hn * hn * (alpha2 - 1) + 1;
        double D = alpha2 / (M_PI * temp * temp);
        // D = std::min(D, (double)1);
        //

        // G
        double k = (alpha2 + 2 * alpha + 1) / 8;
        double G1 = vn / (vn * (1 - k) + k);
        double G2 = ln / (ln * (1 - k) + k);
        double G = G1 * G2;
        //

        // F
        double F0 = pow(ior - 1, 2) / pow(ior + 1, 2);
        double F = F0 + (1 - F0) * pow(1 - hv, 5);
        //

        return D * F * G / (4 * ln * vn) * diffuseColor;
        // return Vector3f::ZERO;
    }
    Vector3f GGX(const Vector3f &v, const Vector3f &n, const Vector3f &l,
                 const Vector3f &light_color) {
        Vector3f h = (v + l).normalized();

        double ln = Vector3f::dot(l, n);
        if (ln >= 0) {
            double vn = Vector3f::dot(v, n);
            double hn = std::min(Vector3f::dot(h, n), .9999);
            double hv = Vector3f::dot(h, v);

            // D

            double temp = hn * hn * (alpha2 - 1) + 1;
            double D = alpha2 / (M_PI * temp * temp);
            //

            // G
            double k = (alpha2 + 2 * alpha + 1) / 8;
            double G1 = vn / (vn * (1 - k) + k);
            double G2 = ln / (ln * (1 - k) + k);
            double G = G1 * G2;
            //

            // F
            double F0 = pow(ior - 1, 2) / pow(ior + 1, 2);
            double F = F0 + (1 - F0) * pow(1 - hv, 5);
            //

            Vector3f id = kd * diffuseColor / M_PI;
            Vector3f is = ks * D * F * G / (4 * ln * vn) * diffuseColor;

            return light_color * (id + is);
        }

        return Vector3f::ZERO;
    }

   public:
    // double k = 0.2;
    // double roughness = 0.1;
    // double fresnel = 0.1;
    double kd = 0.0;  // 漫反射占比
    double ks = 1.0;  // 镜面反射占比
    double kr = 0;    // 透射占比
    double alpha = 0.01;
    double alpha2 = 0.0001;
    double ior = 30;

    Vector3f diffuseColor;
    Vector3f specularColor;
    double shininess;
};

#endif  // MATERIAL_H
