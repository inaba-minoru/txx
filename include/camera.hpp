#ifndef CAMERA_H
#define CAMERA_H

// #include <double.h>
#include <vecmath.h>

#include <cmath>

#include "ray.hpp"

class Camera {
   public:
    Camera(const Vector3f &center, const Vector3f &direction,
           const Vector3f &up, int imgW, int imgH) {
        this->center = center;
        this->direction = direction.normalized();
        // this->horizontal = Vector3f::cross(this->direction, up).normalized();
        this->horizontal = Vector3f::cross(this->direction, up).normalized();
        this->up =
            Vector3f::cross(this->horizontal, this->direction).normalized();
        this->width = imgW;
        this->height = imgH;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(unsigned short *Xi, const Vector2f &point) = 0;
    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

   protected:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
};

// TODO: Implement Perspective camera ------ OK
// You can add new functions or variables whenever needed.
class PerspectiveCamera : public Camera {
   private:
    double f;

   public:
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
                      const Vector3f &up, int imgW, int imgH, double angle)
        : Camera(center, direction, up, imgW, imgH) {
        // angle is in radian.
        f = imgH / tan(angle / 2) / 2;
    }

    Ray generateRay(unsigned short *Xi, const Vector2f &point) override {
        Vector3f dc = Matrix3f(horizontal, up, direction) *
                      Vector3f(point.x() - width / 2, point.y() - height / 2, f)
                          .normalized();
        // return Ray(center + dc * 140, dc);
        return Ray(center, dc);
    }
};

class RealisticCamera : public Camera {
   public:
    RealisticCamera(const Vector3f &center, const Vector3f &direction,
                    const Vector3f &up, int imgW, int imgH, const double &d,
                    const double &zi, const double &zo, const double &angle)
        : d(d), zi(zi), zo(zo), Camera(center, direction, up, imgW, imgH) {
        // angle = 60 / 180. * M_PI;

        double f = imgH / tan(angle / 2) / 2;
        // double temp = sqrt(imgH * imgH / 4 + imgW * imgW / 4);
        // double new_angle = atan(temp / f);
        double new_angle = atan(imgW / 2 / f);
        double final_angle = std::max(angle / 2, new_angle);

        zs = d + zi;
        // C = tan(angle / 2) * zs;
        C = tan(final_angle) * zs;
        this->R = zi / d * C;
        printf("%lf %e %e\n", d, C, this->R);
    }

    double d, zi, zs, zo, R, C;
    Ray generateRay(unsigned short *Xi, const Vector2f &point) override {
        double fuck = std::max(height, width);
        Vector3f temp((width - 2 * point.x()) / fuck * C,
                      (height - 2 * point.y()) / fuck * C, 0);
        temp = (Vector3f(0, 0, zs) - temp) / zs * (zo + zs);
        // printf("%lf\n",
        //        acos(Vector3f::dot(temp.normalized(), Vector3f(0, 0, 1))) /
        //        M_PI * 180);

        double theta = erand48(Xi) * 2 * M_PI, r = sqrt(erand48(Xi)) * R;
        Vector3f orig = Vector3f(cos(theta) * r, sin(theta) * r, zs);
        temp = (temp - orig).normalized();

        orig = Matrix3f(horizontal, up, direction) * orig;
        Vector3f dc = Matrix3f(horizontal, up, direction) * temp;
// return Ray(center + dc * 140, dc);
#pragma omp critical
        {
            //             puts("here");
            //             center.print();
            //             orig.print();
            // dc.print();
        }
        // printf("%lf\n",
        //        acos(Vector3f::dot(dc, Vector3f(0, 0, -1))) / M_PI * 180);
        return Ray(center + orig - zs * direction, dc);
    }
};

#endif  // CAMERA_H
