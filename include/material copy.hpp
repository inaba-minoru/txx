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
    explicit Material(const Vector3f &d_color,
                      const Vector3f &s_color = Vector3f::ZERO, float s = 0)
        : diffuseColor(d_color), specularColor(s_color), shininess(s) {}

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const { return diffuseColor; }

    Vector3f Shade(const Ray &ray, const Hit &hit, const Vector3f &dirToLight,
                   const Vector3f &lightColor) {
        // Vector3f shaded = Vector3f::ZERO;

        Vector3f kdre =
            diffuseColor *
            std::max(float(0), Vector3f::dot(dirToLight.normalized(),
                                             hit.getNormal().normalized()));
        Vector3f ksre =
            specularColor *
            std::pow(std::max(float(0),
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
        //             std::max(float(0),
        //                      Vector3f::dot(dirToLight, hit.getNormal())) +
        //         specularColor *
        //             pow(std::max(
        //                     float(0),
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

    //     float ln = Vector3f::dot(l, n);
    //     float s = 0.0;
    //     if (ln > 0.0) {
    //         float hn = Vector3f::dot(h, n);
    //         float vn = Vector3f::dot(v, n);
    //         float hv = Vector3f::dot(h, v);

    //         // G项
    //         float hn2 = 2.0 * hn;
    //         float g1 = (hn2 * ln) / hv;
    //         float g2 = (hn2 * vn) / hv;
    //         float g = std::min((float)1.0, std::min(g1, g2));

    //         // D项：beckmann distribution function
    //         float m2 = roughness * roughness;
    //         float r1 = 1.0 / (4.0 * m2 * pow(hn, 4.0));
    //         float r2 = (hn * hn - 1.0) / (m2 * hn * hn);
    //         float d = r1 * exp(r2);

    //         // F项
    //         float f = pow(1.0 - hv, 5.0);
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

        float ln = Vector3f::dot(l, n);
        float vn = Vector3f::dot(v, n);
        float hn = Vector3f::dot(h, n);
        float hv = Vector3f::dot(h, v);

        // D
        float temp = hn * hn * (alpha2 - 1) + 1;
        float D = alpha2 / (M_PI * temp * temp);
        // D = std::min(D, (float)1);
        //

        // G
        float k = (alpha2 + 2 * alpha + 1) / 8;
        float G1 = vn / (vn * (1 - k) + k);
        float G2 = ln / (ln * (1 - k) + k);
        float G = G1 * G2;
        //

        // F
        float F0 = pow(ior - 1, 2) / pow(ior + 1, 2);
        float F = F0 + (1 - F0) * pow(1 - hv, 5);
        //

        return D * F * G / (4 * ln * vn) * diffuseColor;
        // return Vector3f::ZERO;
    }
    Vector3f GGX(const Vector3f &v, const Vector3f &n, const Vector3f &l,
                 const Vector3f &light_color) {
        Vector3f h = (v + l).normalized();

        float ln = Vector3f::dot(l, n);
        if (ln >= 0) {
            float vn = Vector3f::dot(v, n);
            float hn = Vector3f::dot(h, n);
            float hv = Vector3f::dot(h, v);

            // D
            float temp = hn * hn * (alpha2 - 1) + 1;
            float D = alpha2 / (M_PI * temp * temp);
            //

            // G
            float k = (alpha2 + 2 * alpha + 1) / 8;
            float G1 = vn / (vn * (1 - k) + k);
            float G2 = ln / (ln * (1 - k) + k);
            float G = G1 * G2;
            //

            // F
            float F0 = pow(ior - 1, 2) / pow(ior + 1, 2);
            float F = F0 + (1 - F0) * pow(1 - hv, 5);
            //

            Vector3f id = kd * diffuseColor / M_PI;
            Vector3f is = ks * D * F * G / (4 * ln * vn) * diffuseColor;

            return light_color * (id + is);
        }

        return Vector3f::ZERO;
    }

   public:
    // float k = 0.2;
    // float roughness = 0.1;
    // float fresnel = 0.1;
    float kd = 0.6;  // 漫反射占比
    float ks = 0.4;  // 镜面反射占比
    float kr = 0;    // 透射占比
    float alpha = 0.01;
    float alpha2 = 0.0001;
    float ior = 30;

    Vector3f diffuseColor;
    Vector3f specularColor;
    float shininess;
};

#endif  // MATERIAL_H
