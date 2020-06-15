#ifndef SPHERE_H
#define SPHERE_H

#include <vecmath.h>

#include <cmath>

#include "object3d.hpp"

// TODO: Implement functions and add more fields as necessary ------ OK
// 球
class Sphere : public Object3D {
   public:
    Vector3f centre;  // 球心
    float radius;     // 半径

    Sphere() = delete;
    explicit Sphere(const Vector3f &c, const float &r)
        : centre(c), radius(r), Object3D() {}
    explicit Sphere(const Vector3f &c, const float &r, Material *m)
        : centre(c), radius(r), Object3D(m) {}

    ~Sphere() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        Vector3f l = centre - r.getOrigin();  // 光源到球心
        float t =
            Vector3f::dot(l, r.getDirection().normalized());  // 光源到垂线距离

        float d2 = l.squaredLength() - t * t;  // 球心到光线距离平方
        if (d2 >= radius * radius) {           // 未击中球，舍弃
            return false;
        }

        float t_ = sqrt(radius * radius - d2);  // 垂足到球面
        float tt;
        if (t - t_ >= tmin)
            tt = t - t_;
        else
            t = t + t_;
        // if (l.squaredLength() < radius * radius) {  // 若光源在球内
        //     tt = t + t_;
        // } else {  // 光源在球外
        //     tt = t - t_;
        // }

        tt /= r.getDirection().length();

        if (tt >= tmin && tt < h.getT()) {  // tt合法并且比之前的小
            h.set(tt, material,
                  (r.pointAtParameter(tt) - centre).normalized() *
                      (l.squaredLength() < radius * radius ? -1 : 1));
            return true;
        }

        return false;
    }
};

#endif
