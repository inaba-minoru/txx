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
        Vector3f shaded = Vector3f::ZERO;

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

        return lightColor *
               (diffuseColor *
                    std::max(float(0),
                             Vector3f::dot(dirToLight, hit.getNormal())) +
                specularColor *
                    pow(std::max(
                            float(0),
                            Vector3f::dot(
                                -ray.getDirection(),
                                2 * Vector3f::dot(dirToLight, hit.getNormal()) *
                                        hit.getNormal() -
                                    dirToLight)),
                        shininess));

        return shaded;
    }

   protected:
    Vector3f diffuseColor;
    Vector3f specularColor;
    float shininess;
};

#endif  // MATERIAL_H
