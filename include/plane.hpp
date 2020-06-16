#ifndef PLANE_H
#define PLANE_H

#include <vecmath.h>

#include <cmath>

#include "object3d.hpp"

// TODO: Implement Plane representing an infinite plane ------ OK
// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions
// 平面
class Plane : public Object3D {
   private:
    Vector3f normal;  // 法线
    double d;          // ax + by + cz = d

   public:
    Plane() = delete;
    explicit Plane(const Vector3f &n, const double &d)
        : normal(n), d(d), Object3D() {}
    explicit Plane(const Vector3f &n, const double &d, Material *m)
        : normal(n), d(d), Object3D(m) {}

    ~Plane() override = default;

    bool intersect(const Ray &r, Hit &h, double tmin) override {
        if (Vector3f::dot(normal, r.getDirection()) == 0) {  // 光线与平面平行
            return false;
        }

        double tt = (d - Vector3f::dot(normal, r.getOrigin())) /
                   Vector3f::dot(normal, r.getDirection());

        if (tt >= tmin && tt < h.getT()) {
            h.set(tt, material,
                  normal *
                      (Vector3f::dot(normal, r.getDirection()) > 0 ? -1 : 1));

            return true;
        }

        return false;
    }
};

#endif  // PLANE_H
