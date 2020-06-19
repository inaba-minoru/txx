#ifndef OBJECT3D_H
#define OBJECT3D_H

#include "hit.hpp"
#include "material.hpp"
#include "ray.hpp"

// Base class for all 3d entities.
class Object3D {
   public:
    Object3D() : material(nullptr) {}

    virtual ~Object3D() = default;

    explicit Object3D(Material *material) { this->material = material; }

    // Intersect Ray with this object. If hit, store information in hit
    // structure.
    virtual bool intersect(const Ray &r, Hit &h, double tmin) = 0;

    double area = 0;
    virtual void sampleLight(unsigned short *Xi, Vector3f &p, Vector3f &light) {
        //   n = Vector3f(1);
        light = Vector3f::ZERO;
    }

    Material *material;

    static Vector3f transformDirection(const Matrix4f &mat,
                                       const Vector3f &dir) {
        return (mat * Vector4f(dir, 0)).xyz();
    }
    static Vector3f transformPoint(const Matrix4f &mat, const Vector3f &point) {
        return (mat * Vector4f(point, 1)).xyz();
    }

    virtual void transform(const Matrix4f &mat) {}

   protected:
};

#endif
