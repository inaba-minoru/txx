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
    Vector3f normals[3];   // 法线

    Triangle() = delete;
    explicit Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c)
        : vertices({a, b, c}), Object3D() {
        normals[0] = normals[1] = normals[2] =
            Vector3f::cross(b - a, c - a).normalized();
    }
    explicit Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c,
                      Material *m)
        : vertices({a, b, c}), Object3D(m) {
        normals[0] = normals[1] = normals[2] =
            Vector3f::cross(b - a, c - a).normalized();
    }
    explicit Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c,
                      const Vector3f &n0, const Vector3f &n1,
                      const Vector3f &n2, Material *m)
        : vertices({a, b, c}), normals({n0, n1, n2}), Object3D(m) {}

    ~Triangle() override = default;

    bool intersect(const Ray &r, Hit &h, double tmin) override {
        Vector3f E1 = vertices[0] - vertices[1], E2 = vertices[0] - vertices[2],
                 S = vertices[0] - r.getOrigin();
        Vector3f tby =
            Vector3f(Matrix3f(S, E1, E2).determinant(),
                     Matrix3f(r.getDirection(), S, E2).determinant(),
                     Matrix3f(r.getDirection(), E1, S).determinant()) /
            Matrix3f(r.getDirection(), E1, E2).determinant();

        // tby.print();

        if (tby.x() >= tmin && tby.x() < h.getT() && 0 <= tby.y() &&
            0 <= tby.z() && tby.y() + tby.z() <= 1) {
            Vector3f norm = ((1 - tby.y() - tby.z()) * normals[0] +
                             tby.y() * normals[1] + tby.z() * normals[2])
                                .normalized();
            // Vector3f norm = Vector3f::cross(vertices[0] - vertices[1],
            //                                 vertices[0] - vertices[2])
            //                     .normalized();
            // printf("%lf\n", Vector3f::dot(norm, norm_));

            h.set(tby.x(), material,
                  norm * (Vector3f::dot(norm, r.getDirection()) > 0 ? -1 : 1));

            if (hasTex)
                h.setTexCoord(((1 - tby.y() - tby.z()) * texCoords[0] +
                               tby.y() * texCoords[1] +
                               tby.z() * texCoords[2]));

            return true;
        }

        return false;
    }

    bool hasTex = 0;
    Vector2f texCoords[3];  // 贴图
    void setTexCoords(const Vector2f &texCoord0, const Vector2f &texCoord1,
                      const Vector2f &texCoord2) {
        this->texCoords[0] = texCoord0;
        this->texCoords[1] = texCoord1;
        this->texCoords[2] = texCoord2;
        hasTex = 1;
    }
};

#endif  // TRIANGLE_H
