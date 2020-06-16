#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <vecmath.h>

#include <cmath>
#include <iostream>

#include "object3d.hpp"
using namespace std;

// TODO: implement this class and add more fields as necessary, ------ OK
// 三角
class Triangle : public Object3D {
   public:
    Vector3f vertices[3];  // 三顶点
    Vector3f normal;       // 法线

    Triangle() = delete;
    explicit Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c)
        : vertices({a, b, c}), Object3D() {
        normal = Vector3f::cross(b - a, c - a).normalized();
    }
    explicit Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c,
                      Material *m)
        : vertices({a, b, c}), Object3D(m) {
        normal = Vector3f::cross(b - a, c - a).normalized();
    }

    ~Triangle() override = default;

    bool intersect(const Ray &r, Hit &h, double tmin) override {
        Vector3f E1 = vertices[0] - vertices[1], E2 = vertices[0] - vertices[2],
                 S = vertices[0] - r.getOrigin();
        Vector3f tby =
            Vector3f(Matrix3f(S, E1, E2).determinant(),
                     Matrix3f(r.getDirection(), S, E2).determinant(),
                     Matrix3f(r.getDirection(), E1, S).determinant()) /
            Matrix3f(r.getDirection(), E1, E2).determinant();

        if (tby.x() >= tmin && tby.x() < h.getT() && 0 <= tby.y() &&
            0 <= tby.z() && tby.y() + tby.z() <= 1) {
            h.set(tby.x(), material,
                  normal *
                      (Vector3f::dot(normal, r.getDirection()) > 0 ? -1 : 1));
            return true;
        }

        return false;
    }
};

#endif  // TRIANGLE_H
