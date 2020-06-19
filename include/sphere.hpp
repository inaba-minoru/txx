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
    double radius;    // 半径

    Sphere() = delete;
    explicit Sphere(const Vector3f &c, const double &r)
        : centre(c), radius(r), Object3D() {
        area = 4 * M_PI * radius * radius;
    }
    explicit Sphere(const Vector3f &c, const double &r, Material *m)
        : centre(c), radius(r), Object3D(m) {
        area = 4 * M_PI * radius * radius;
    }

    ~Sphere() override = default;

    bool intersect(const Ray &r, Hit &h, double tmin) override {
        // Vector3f op = centre - r.getOrigin();
        //   Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        // double t, eps = 1e-4, b = Vector3f::dot(op, r.getDirection()),
        //           det = b * b - op.squaredLength() + radius * radius;
        // if (det < 0)
        //     return 0;
        // else {
        //     det = sqrt(det);
        //     if (b - det >= tmin)
        //         t = b - det;
        //     else
        //         t = b + det;
        //     if (t >= tmin && t < h.getT()) {  // tt合法并且比之前的小
        //         Vector3f n = (r.pointAtParameter(t) - centre).normalized();
        //         if (Vector3f::dot(n, -r.getDirection()) < 0) n = -n;

        //         h.set(t, material, n);
        //         return true;
        //     }
        //     return false;
        // }
        // return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);

        Vector3f l = centre - r.getOrigin();  // 光源到球心
        double t =
            Vector3f::dot(l, r.getDirection().normalized());  // 光源到垂线距离

        double d2 = l.squaredLength() - t * t;  // 球心到光线距离平方
        if (d2 >= radius * radius) {            // 未击中球，舍弃
            return false;
        }

        double t_ = sqrt(radius * radius - d2);  // 垂足到球面
        double tt;
        if (t - t_ >= tmin)
            tt = t - t_;
        else
            tt = t + t_;
        // if (l.squaredLength() < radius * radius) {  // 若光源在球内
        //     tt = t + t_;
        // } else {  // 光源在球外
        //     tt = t - t_;
        // }

        // // #pragma omp critical
        // //         {
        // //             puts("yeahhhh");
        // //             printf("%f %f\n", tt, h.getT());
        // //             r.getDirection().print();
        // //             centre.print();
        // //         }

        tt /= r.getDirection().length();

        if (tt >= tmin && tt < h.getT()) {  // tt合法并且比之前的小
            Vector3f n = (r.pointAtParameter(tt) - centre).normalized();
            if (Vector3f::dot(n, -r.getDirection()) < 0) n = -n;

            h.set(tt, material, n);
            return true;
        }

        return false;
    }

    void sampleLight(unsigned short *Xi, Vector3f &p, Vector3f &light) {
        double theta = erand48(Xi) * M_PI;
        double phi = erand48(Xi) * 2 * M_PI;
        double temp = sin(theta);
        Vector3f n = Vector3f(cos(phi) * temp, sin(phi) * temp, cos(theta));
        p = centre + radius * n;
        light = material->emission;
    }
};

#endif
