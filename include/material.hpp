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
};

#endif  // MATERIAL_H
